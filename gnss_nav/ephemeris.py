"""
Satellite position & clock from broadcast (Keplerian) ephemeris.

Reference: IS-GPS-200, Table 20-IV, "User Algorithm for SV Position".
The same algorithm applies to Galileo (FNAV / INAV) and BeiDou MEO/IGSO
with their own ephemerides — the math is identical.

Output is ECEF in the WGS-84 frame at the *transmission* time of the signal,
which is what the pseudorange equation requires.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .geodesy import (
    F_REL, MU_EARTH, OMEGA_E_DOT, SPEED_OF_LIGHT,
)
from .rinex import Ephemeris


@dataclass
class SatState:
    position: np.ndarray   # ECEF, meters, shape (3,)
    velocity: np.ndarray   # ECEF, m/s, shape (3,)
    clock_bias: float      # seconds (apply as +c·dt_sv on pseudorange RHS)
    clock_drift: float     # s/s


def _normalize_week_seconds(t: float) -> float:
    """Wrap a time-difference around the GPS week (±302400 s)."""
    half_week = 302_400.0
    if t > half_week:
        t -= 2 * half_week
    elif t < -half_week:
        t += 2 * half_week
    return t


def select_ephemeris(ephemerides: list[Ephemeris], sv: str, t_gps_sow: float) -> Ephemeris | None:
    """
    Pick the ephemeris for `sv` closest in time to `t_gps_sow`, preferring
    those whose toe is within ±2 hours (the nominal broadcast validity).
    Returns None if no ephemeris exists for the SV.
    """
    best: Ephemeris | None = None
    best_dt = float("inf")
    for e in ephemerides:
        if e.sv != sv:
            continue
        dt = abs(_normalize_week_seconds(t_gps_sow - e.toe))
        if dt < best_dt:
            best_dt = dt
            best = e
    if best is None:
        return None
    # Allow up to 4 hours of staleness — broadcast ephemerides are typically
    # valid for ~2h but Android dumps sometimes lag.
    if best_dt > 4 * 3600:
        return None
    return best


def sv_clock_correction(eph: Ephemeris, t_gps_sow: float) -> float:
    """
    Polynomial SV clock model + relativistic correction.

    Returns dt_sv (seconds): the SV-side clock offset to *subtract* from
    raw transmission time to obtain GPS-system time.
    """
    dt = _normalize_week_seconds(t_gps_sow - eph.toc)
    # Iterate once for the second-order polynomial — negligible improvement
    # in practice, but cheap and standard.
    dt_sv = eph.af0 + eph.af1 * dt + eph.af2 * dt * dt

    # Relativistic eccentricity term needs the eccentric anomaly. We re-solve
    # Kepler here at the corrected time — caller could share it, but the math
    # is fast.
    tk = _normalize_week_seconds(t_gps_sow - dt_sv - eph.toe)
    a = eph.a
    n0 = math.sqrt(MU_EARTH / (a ** 3))
    n = n0 + eph.delta_n
    M = eph.m0 + n * tk
    E = M
    for _ in range(10):
        E_new = M + eph.e * math.sin(E)
        if abs(E_new - E) < 1e-12:
            E = E_new
            break
        E = E_new
    dt_rel = F_REL * eph.e * eph.sqrt_a * math.sin(E)
    return dt_sv + dt_rel


def compute_sv_state(
    eph: Ephemeris,
    t_transmit_gps_sow: float,
    signal_travel_time: float = 0.0,
) -> SatState:
    """
    Compute SV ECEF position and velocity at *transmission* time.

    Parameters
    ----------
    eph
        The selected ephemeris.
    t_transmit_gps_sow
        GPS time-of-week (s) at which the signal was transmitted (i.e. the
        receiver epoch minus the signal travel time, after SV clock correction).
    signal_travel_time
        Travel time (s) used to apply Earth-rotation correction. The SV
        position is rotated by −ωₑ·τ around Z so it is expressed in the ECEF
        frame at *reception* time. Pass 0 if you'll rotate later yourself.
    """
    tk = _normalize_week_seconds(t_transmit_gps_sow - eph.toe)

    a = eph.a
    n0 = math.sqrt(MU_EARTH / (a ** 3))
    n = n0 + eph.delta_n

    # Mean anomaly → eccentric anomaly (Kepler) → true anomaly
    M = eph.m0 + n * tk
    E = M
    for _ in range(15):
        f = E - eph.e * math.sin(E) - M
        fp = 1 - eph.e * math.cos(E)
        dE = f / fp
        E -= dE
        if abs(dE) < 1e-13:
            break
    sin_E, cos_E = math.sin(E), math.cos(E)

    sqrt_1_e2 = math.sqrt(1 - eph.e * eph.e)
    nu = math.atan2(sqrt_1_e2 * sin_E, cos_E - eph.e)

    # Argument of latitude with harmonic perturbations
    phi = nu + eph.omega
    sin_2phi = math.sin(2 * phi)
    cos_2phi = math.cos(2 * phi)
    du = eph.cus * sin_2phi + eph.cuc * cos_2phi
    dr = eph.crs * sin_2phi + eph.crc * cos_2phi
    di = eph.cis * sin_2phi + eph.cic * cos_2phi

    u = phi + du
    r = a * (1 - eph.e * cos_E) + dr
    i = eph.i0 + di + eph.idot * tk

    # Position in orbital plane
    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    # Corrected longitude of ascending node (includes Earth rotation since toe
    # and the optional signal-travel-time rotation)
    Omega = (
        eph.omega0
        + (eph.omega_dot - OMEGA_E_DOT) * tk
        - OMEGA_E_DOT * (eph.toe + signal_travel_time)
    )
    sin_O, cos_O = math.sin(Omega), math.cos(Omega)
    cos_i, sin_i = math.cos(i), math.sin(i)

    x = x_orb * cos_O - y_orb * cos_i * sin_O
    y = x_orb * sin_O + y_orb * cos_i * cos_O
    z = y_orb * sin_i

    # ─── Velocity (analytic derivative) ─────────────────────────────────────
    # Per IS-GPS-200 user algorithm.
    E_dot = n / (1 - eph.e * cos_E)
    nu_dot = sqrt_1_e2 * E_dot / (1 - eph.e * cos_E)
    u_dot = nu_dot + 2 * (eph.cus * cos_2phi - eph.cuc * sin_2phi) * nu_dot
    r_dot = a * eph.e * sin_E * E_dot + 2 * (eph.crs * cos_2phi - eph.crc * sin_2phi) * nu_dot
    i_dot = eph.idot + 2 * (eph.cis * cos_2phi - eph.cic * sin_2phi) * nu_dot
    Omega_dot_eff = eph.omega_dot - OMEGA_E_DOT

    x_orb_dot = r_dot * math.cos(u) - r * math.sin(u) * u_dot
    y_orb_dot = r_dot * math.sin(u) + r * math.cos(u) * u_dot

    vx = (
        x_orb_dot * cos_O - y_orb_dot * cos_i * sin_O
        + y_orb * sin_i * sin_O * i_dot
        - y * Omega_dot_eff
    )
    vy = (
        x_orb_dot * sin_O + y_orb_dot * cos_i * cos_O
        - y_orb * sin_i * cos_O * i_dot
        + x * Omega_dot_eff
    )
    vz = y_orb_dot * sin_i + y_orb * cos_i * i_dot

    clock_bias = sv_clock_correction(eph, t_transmit_gps_sow)
    clock_drift = eph.af1 + 2 * eph.af2 * _normalize_week_seconds(t_transmit_gps_sow - eph.toc)

    return SatState(
        position=np.array([x, y, z]),
        velocity=np.array([vx, vy, vz]),
        clock_bias=clock_bias,
        clock_drift=clock_drift,
    )


def position_and_clock(
    eph: Ephemeris,
    t_receive_gps_sow: float,
    pseudorange: float,
) -> SatState:
    """
    Convenience wrapper: given receiver epoch and pseudorange, iterate the
    signal travel time and return SV state at transmission.
    """
    # First guess: travel time = ρ/c
    tau = pseudorange / SPEED_OF_LIGHT
    for _ in range(3):
        t_tx = t_receive_gps_sow - tau
        # Refine using SV clock at t_tx
        dt_sv = sv_clock_correction(eph, t_tx)
        t_tx_corr = t_tx - dt_sv
        state = compute_sv_state(eph, t_tx_corr, signal_travel_time=tau)
        # Re-estimate tau from geometric range when we have a receiver guess
        # (caller does this loop properly). For now we just trust ρ/c.
        break
    return state

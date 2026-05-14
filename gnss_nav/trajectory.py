"""
Trajectory orchestrator: runs the SPP pipeline over every observation epoch
and produces a list of fixes.

For each epoch:
  1. For each SV with a usable pseudorange:
        a. Pick an ephemeris (closest valid toe).
        b. Iterate signal travel time, computing SV ECEF position + clock.
        c. Apply Earth rotation between transmit and receive.
        d. Apply atmosphere corrections (iono + tropo) using the previous
           epoch's position as the user location (or RINEX header approx
           on the first epoch).
        e. Apply an elevation mask once we have a position estimate.
  2. Solve weighted LSQ for (X, Y, Z, c·δt).
  3. Solve Doppler LSQ for velocity if Doppler is available;
     else fall back to finite-differencing positions.
  4. Convert ECEF → LLA, GPS time → UTC, append to output.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from datetime import datetime

import numpy as np

from .corrections import klobuchar_delay, saastamoinen_delay
from .ephemeris import compute_sv_state, select_ephemeris, sv_clock_correction
from .geodesy import (
    OMEGA_E_DOT, SPEED_OF_LIGHT,
    ecef_to_lla, elevation_azimuth,
)
from .rinex import Epoch, RinexData
from .solver import PositionSolution, SatObservation, solve_position, solve_velocity

log = logging.getLogger(__name__)


@dataclass
class TrajectoryFix:
    time_utc: datetime
    lat_deg: float
    lon_deg: float
    alt_m: float
    x_ecef: float
    y_ecef: float
    z_ecef: float
    vx_ecef: float
    vy_ecef: float
    vz_ecef: float
    speed_mps: float
    num_sats: int
    gdop: float
    clock_bias_m: float


@dataclass
class Settings:
    elevation_mask_deg: float = 10.0
    max_iter: int = 10
    systems: tuple[str, ...] = ("G", "E")   # GPS + Galileo by default
    use_iono: bool = True
    use_tropo: bool = True
    rate_hz: float = 1.0


def _earth_rotate(p: np.ndarray, tau: float) -> np.ndarray:
    """Rotate ECEF point by -ωₑ·τ around Z (SV position into reception frame)."""
    angle = OMEGA_E_DOT * tau
    c, s = math.cos(angle), math.sin(angle)
    return np.array([c * p[0] + s * p[1], -s * p[0] + c * p[1], p[2]])


def _build_sat_observations(
    epoch: Epoch,
    data: RinexData,
    settings: Settings,
    user_guess_ecef: np.ndarray | None,
) -> list[SatObservation]:
    sats: list[SatObservation] = []

    for ob in epoch.observations:
        sys = ob.sv[0]
        if sys not in settings.systems:
            continue
        if ob.pseudorange is None:
            continue

        eph = select_ephemeris(data.ephemerides, ob.sv, epoch.gps_tow)
        if eph is None:
            continue

        # ── Iterate signal travel time ──────────────────────────────────────
        tau = ob.pseudorange / SPEED_OF_LIGHT
        sv_pos = np.zeros(3)
        sv_vel = np.zeros(3)
        dt_sv = 0.0
        for _ in range(3):
            t_tx = epoch.gps_tow - tau
            dt_sv = sv_clock_correction(eph, t_tx)
            state = compute_sv_state(eph, t_tx - dt_sv, signal_travel_time=tau)
            sv_pos = state.position
            sv_vel = state.velocity
            if user_guess_ecef is not None:
                new_tau = float(np.linalg.norm(sv_pos - user_guess_ecef)) / SPEED_OF_LIGHT
            else:
                new_tau = ob.pseudorange / SPEED_OF_LIGHT
            if abs(new_tau - tau) < 1e-9:
                tau = new_tau
                break
            tau = new_tau

        # ── Elevation mask ─────────────────────────────────────────────────
        elev = math.pi / 2  # assume zenith if we have no user guess
        if user_guess_ecef is not None:
            try:
                elev, _ = elevation_azimuth(sv_pos, user_guess_ecef)
            except Exception:
                elev = math.pi / 2
            if math.degrees(elev) < settings.elevation_mask_deg:
                continue

        # ── Atmospheric corrections ───────────────────────────────────────
        iono = 0.0
        tropo = 0.0
        if user_guess_ecef is not None:
            if settings.use_iono and data.iono is not None:
                try:
                    iono = klobuchar_delay(
                        data.iono, user_guess_ecef, sv_pos, epoch.gps_tow,
                        frequency_hz=ob.frequency_hz or 1_575_420_000.0,
                    )
                except Exception:
                    iono = 0.0
            if settings.use_tropo:
                try:
                    tropo = saastamoinen_delay(user_guess_ecef, sv_pos)
                except Exception:
                    tropo = 0.0

        # ρ_corr = ρ + c·dt_sv − I − T
        rho_corr = ob.pseudorange + SPEED_OF_LIGHT * dt_sv - iono - tropo

        # Weighting: sin²(elev) is the textbook choice
        weight = max(math.sin(elev) ** 2, 0.05)

        sats.append(SatObservation(
            sv=ob.sv,
            position=sv_pos,
            velocity=sv_vel,
            pseudorange_corrected=rho_corr,
            doppler_hz=ob.doppler,
            frequency_hz=ob.frequency_hz,
            weight=weight,
        ))

    return sats


def run_trajectory(data: RinexData, settings: Settings | None = None) -> list[TrajectoryFix]:
    settings = settings or Settings()
    if not data.obs_epochs:
        log.warning("No observation epochs in RINEX data.")
        return []
    if not data.ephemerides:
        log.warning("No ephemerides in RINEX data — cannot position.")
        return []

    # Initial position guess: header APPROX POSITION XYZ, else Earth's centre.
    guess = (
        np.array(data.approx_position_ecef)
        if data.approx_position_ecef
        else None
    )

    # Decimate to the requested rate.
    target_dt = 1.0 / settings.rate_hz
    fixes: list[TrajectoryFix] = []
    last_kept_t: float | None = None
    prev_position: np.ndarray | None = None
    prev_time: datetime | None = None

    for epoch in data.obs_epochs:
        t_sec = epoch.time_utc.timestamp()
        if last_kept_t is not None and t_sec - last_kept_t < target_dt - 1e-3:
            continue

        # First pass with the current best guess: get SV positions/elevations.
        sats = _build_sat_observations(epoch, data, settings, guess)
        if len(sats) < 4:
            continue

        sol = solve_position(sats, initial_guess=guess, max_iter=settings.max_iter)
        if sol is None or not sol.converged:
            continue

        # Re-run once with the new position so iono/tropo/elev use the right
        # user location. Single re-run is plenty for SPP.
        sats = _build_sat_observations(epoch, data, settings, sol.position)
        if len(sats) >= 4:
            sol2 = solve_position(sats, initial_guess=sol.position, max_iter=settings.max_iter)
            if sol2 is not None and sol2.converged:
                sol = sol2

        # Velocity
        vel = np.zeros(3)
        v_from_doppler = solve_velocity(sats, sol.position)
        if v_from_doppler is not None:
            vel, _drift = v_from_doppler
        elif prev_position is not None and prev_time is not None:
            dt = (epoch.time_utc - prev_time).total_seconds()
            if dt > 0:
                vel = (sol.position - prev_position) / dt

        lla = ecef_to_lla(*sol.position)
        speed = float(np.linalg.norm(vel))
        fixes.append(TrajectoryFix(
            time_utc=epoch.time_utc,
            lat_deg=lla.lat_deg, lon_deg=lla.lon_deg, alt_m=lla.alt_m,
            x_ecef=float(sol.position[0]),
            y_ecef=float(sol.position[1]),
            z_ecef=float(sol.position[2]),
            vx_ecef=float(vel[0]), vy_ecef=float(vel[1]), vz_ecef=float(vel[2]),
            speed_mps=speed,
            num_sats=sol.num_sats,
            gdop=sol.gdop,
            clock_bias_m=sol.clock_bias_m,
        ))

        guess = sol.position
        prev_position = sol.position
        prev_time = epoch.time_utc
        last_kept_t = t_sec

    return fixes

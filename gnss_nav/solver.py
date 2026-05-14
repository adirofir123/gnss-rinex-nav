"""
Single-point positioning (SPP) least-squares solver.

Given pseudoranges ρᵢ and SV positions sᵢ (already corrected for SV clock and
Earth rotation), we solve

    ρᵢ = ‖sᵢ − r‖ + c·δt + Iᵢ + Tᵢ + εᵢ

for the receiver position r ∈ ℝ³ and clock bias δt (in meters: cb = c·δt).
This is the classic 4-unknown linearized least squares (Gauss-Newton).

Velocity is solved separately from Doppler measurements:

    -λ·Dᵢ = (vₛᵢ − vᵣ) · êᵢ + c·δṫ + ε

where êᵢ is the line-of-sight unit vector from receiver to SV.

References:
- Misra & Enge, "Global Positioning System: Signals, Measurements, Performance"
- https://mason.gmu.edu/~treid5/Math447/GPSEquations/
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .geodesy import SPEED_OF_LIGHT


@dataclass
class SatObservation:
    """A single SV's contribution to one epoch's least-squares solve."""
    sv: str
    position: np.ndarray       # ECEF at transmission, m
    velocity: np.ndarray       # ECEF, m/s
    pseudorange_corrected: float   # ρ − atmosphere + c·dt_sv − c·dt_relativistic
    doppler_hz: float | None = None
    frequency_hz: float | None = None
    weight: float = 1.0


@dataclass
class PositionSolution:
    position: np.ndarray       # ECEF, m
    clock_bias_m: float        # c·δt_rx (meters)
    velocity: np.ndarray       # ECEF, m/s; zeros if not estimated
    clock_drift_mps: float = 0.0
    residuals: np.ndarray | None = None
    num_sats: int = 0
    gdop: float = float("nan")
    converged: bool = False


def solve_position(
    sats: list[SatObservation],
    initial_guess: np.ndarray | None = None,
    max_iter: int = 10,
    convergence_m: float = 1e-3,
) -> PositionSolution | None:
    """
    Weighted least-squares solve for (X, Y, Z, c·δt). Needs ≥4 sats.

    Iterates Gauss-Newton from `initial_guess` (defaults to Earth's centre,
    which converges fine in practice — but a better guess is faster).
    """
    if len(sats) < 4:
        return None

    if initial_guess is None:
        x = np.array([0.0, 0.0, 0.0, 0.0])
    else:
        x = np.array([initial_guess[0], initial_guess[1], initial_guess[2], 0.0])

    W = np.diag([s.weight for s in sats])
    sat_pos = np.array([s.position for s in sats])           # (n, 3)
    rho     = np.array([s.pseudorange_corrected for s in sats])  # (n,)

    converged = False
    last_dx = np.inf
    for _ in range(max_iter):
        r = x[:3]
        diff = sat_pos - r
        ranges = np.linalg.norm(diff, axis=1)
        unit_los = diff / ranges[:, None]            # (n, 3); SV - rx, normalized

        predicted = ranges + x[3]
        residuals = rho - predicted

        # Geometry matrix H: ∂ρ/∂(x,y,z,cb) = (-êᵢᵀ, +1)
        H = np.hstack((-unit_los, np.ones((len(sats), 1))))

        # Weighted normal equations
        HtW = H.T @ W
        N = HtW @ H
        try:
            dx = np.linalg.solve(N, HtW @ residuals)
        except np.linalg.LinAlgError:
            return None

        x = x + dx
        last_dx = np.linalg.norm(dx[:3])
        if last_dx < convergence_m:
            converged = True
            break

    # GDOP from unweighted normal matrix
    try:
        H_geom = np.hstack((-unit_los, np.ones((len(sats), 1))))
        Q = np.linalg.inv(H_geom.T @ H_geom)
        gdop = math.sqrt(np.trace(Q))
    except np.linalg.LinAlgError:
        gdop = float("nan")

    return PositionSolution(
        position=x[:3],
        clock_bias_m=float(x[3]),
        velocity=np.zeros(3),
        residuals=residuals,
        num_sats=len(sats),
        gdop=gdop,
        converged=converged,
    )


def solve_velocity(
    sats: list[SatObservation],
    rx_position: np.ndarray,
) -> tuple[np.ndarray, float] | None:
    """
    Estimate receiver velocity (m/s, ECEF) and clock drift (m/s) from Doppler.

    Range-rate model:
        ρ̇ᵢ = -λᵢ · Dᵢ = (vₛᵢ - vᵣ) · êᵢ - c·δṫₛᵢ + c·δṫᵣ + ε
    We treat the SV clock drift as small (broadcast af1 is ~1e-11 s/s — a few
    mm/s — usually folded in elsewhere or ignored at SPP precision).
    """
    usable = [s for s in sats if s.doppler_hz is not None and s.frequency_hz]
    if len(usable) < 4:
        return None

    A = np.zeros((len(usable), 4))
    b = np.zeros(len(usable))
    W = np.diag([s.weight for s in usable])

    for k, s in enumerate(usable):
        diff = s.position - rx_position
        rng = np.linalg.norm(diff)
        e = diff / rng                       # line of sight, rx → sv

        # Pseudorange-rate from Doppler: ρ̇ = -λ·D  (positive Doppler ⇒
        # SV approaching ⇒ range decreasing). λ = c/f.
        lam = SPEED_OF_LIGHT / s.frequency_hz
        rho_dot = -lam * s.doppler_hz

        # rho_dot - vₛ·ê = -vᵣ·ê + c·δṫ
        b[k] = rho_dot - float(np.dot(s.velocity, e))
        A[k, :3] = -e
        A[k, 3] = 1.0

    try:
        AtW = A.T @ W
        sol = np.linalg.solve(AtW @ A, AtW @ b)
    except np.linalg.LinAlgError:
        return None

    return sol[:3], float(sol[3])

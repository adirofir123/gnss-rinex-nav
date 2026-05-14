"""
Quick end-to-end sanity check.

We build a known receiver position, propagate 6 SVs from a hand-rolled
ephemeris, generate fake pseudoranges, feed them through the solver, and
check we recover the receiver within a few cm.

This bypasses the RINEX parser (we tested its types via imports) and
focuses on the math being self-consistent.
"""

import math
import numpy as np
from datetime import datetime, timezone

from gnss_nav.ephemeris import compute_sv_state, Ephemeris
from gnss_nav.solver import SatObservation, solve_position, solve_velocity
from gnss_nav.geodesy import SPEED_OF_LIGHT, lla_to_ecef, ecef_to_lla


# True receiver: Jerusalem, on the ground
rx_true = np.array(lla_to_ecef(31.7683, 35.2137, 754.0))
clock_bias_true = 12345.6  # m  (about 41 µs of receiver clock offset)

# 6 fake "satellites" placed at random-ish unit vectors around the sky,
# at GPS orbital altitude (~20,200 km).
rng = np.random.default_rng(42)
sv_unit = rng.normal(size=(6, 3))
sv_unit /= np.linalg.norm(sv_unit, axis=1)[:, None]
# Push above horizon: positive z-component, then bias toward up
sv_unit[:, 2] = np.abs(sv_unit[:, 2]) + 0.3
sv_unit /= np.linalg.norm(sv_unit, axis=1)[:, None]
sv_positions = rx_true + sv_unit * 20_200_000.0

# Build fake observations
sats = []
for i, sv_p in enumerate(sv_positions):
    geometric = float(np.linalg.norm(sv_p - rx_true))
    rho = geometric + clock_bias_true   # no atmosphere, no noise
    sats.append(SatObservation(
        sv=f"G{i+1:02d}",
        position=sv_p,
        velocity=np.zeros(3),
        pseudorange_corrected=rho,
    ))

sol = solve_position(sats, max_iter=15)
print("Converged:", sol.converged)
print("Num sats:", sol.num_sats)
print("Recovered position error (m):", np.linalg.norm(sol.position - rx_true))
print("Recovered clock bias error (m):", abs(sol.clock_bias_m - clock_bias_true))
print("GDOP:", sol.gdop)

# Velocity test: receiver moving east at 10 m/s
rx_vel_true = np.array([0.0, 10.0, 0.0])
clock_drift_true = 0.5  # m/s

for s in sats:
    e = (s.position - rx_true) / np.linalg.norm(s.position - rx_true)
    rho_dot = float(np.dot(np.zeros(3) - rx_vel_true, e)) + clock_drift_true
    # Convert range-rate to Doppler at L1 (rho_dot = -λ·D ⇒ D = -rho_dot/λ)
    lam = SPEED_OF_LIGHT / 1_575_420_000.0
    s.doppler_hz = -rho_dot / lam
    s.frequency_hz = 1_575_420_000.0

vel_result = solve_velocity(sats, rx_true)
v_est, drift_est = vel_result
print("\nVelocity test:")
print("Estimated velocity (m/s):", v_est)
print("Velocity error (m/s):", np.linalg.norm(v_est - rx_vel_true))
print("Clock drift estimate (m/s):", drift_est, "  true:", clock_drift_true)

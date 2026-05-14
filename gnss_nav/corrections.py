"""
Atmospheric range corrections.

* Klobuchar — single-frequency broadcast ionospheric model. Coefficients come
  from the navigation file header. Models delay on the L1 frequency; we scale
  for other bands using the standard f²-ratio.

* Saastamoinen — dry+wet tropospheric delay with a simple mapping function.
  Adequate for SPP accuracy; doesn't need surface meteorology.

Both return *delay in meters* (positive = path is longer than vacuum geometric
range, i.e. you subtract this from the pseudorange before solving).
"""

from __future__ import annotations

import math

from .geodesy import SPEED_OF_LIGHT, ecef_to_lla, elevation_azimuth
import numpy as np

from .rinex import IonoParams


# ─── Klobuchar ─────────────────────────────────────────────────────────────

def klobuchar_delay(
    iono: IonoParams,
    rx_ecef: np.ndarray,
    sv_ecef: np.ndarray,
    gps_tow: float,
    frequency_hz: float = 1_575_420_000.0,
) -> float:
    """
    Klobuchar single-frequency ionospheric correction.

    Returns delay in meters on the given frequency (defaults to L1).
    """
    rx_lla = ecef_to_lla(*rx_ecef)
    elev, azim = elevation_azimuth(sv_ecef, rx_ecef)
    # Convert to semicircles
    phi_u = rx_lla.lat_deg / 180.0
    lam_u = rx_lla.lon_deg / 180.0
    E = elev / math.pi
    A = azim
    # Earth-centred angle (semicircles)
    psi = 0.0137 / (E + 0.11) - 0.022
    # Sub-ionospheric lat
    phi_i = phi_u + psi * math.cos(A)
    phi_i = max(min(phi_i, 0.416), -0.416)
    # Sub-ionospheric lon
    lam_i = lam_u + psi * math.sin(A) / math.cos(phi_i * math.pi)
    # Geomagnetic lat
    phi_m = phi_i + 0.064 * math.cos((lam_i - 1.617) * math.pi)
    # Local time at sub-iono point
    t = 43200.0 * lam_i + gps_tow
    t = t % 86400.0
    if t < 0:
        t += 86400.0
    # Amplitude / period (clamp negatives to 0 / 72000)
    AMP = sum(iono.alpha[n] * (phi_m ** n) for n in range(4))
    if AMP < 0:
        AMP = 0.0
    PER = sum(iono.beta[n] * (phi_m ** n) for n in range(4))
    if PER < 72000.0:
        PER = 72000.0
    # Phase
    x = 2 * math.pi * (t - 50400.0) / PER
    # Slant factor
    F = 1.0 + 16.0 * (0.53 - E) ** 3

    if abs(x) < 1.57:
        T_iono = F * (5e-9 + AMP * (1 - x * x / 2 + x ** 4 / 24))
    else:
        T_iono = F * 5e-9

    delay_l1 = T_iono * SPEED_OF_LIGHT
    # Scale to actual frequency: delay ∝ 1/f²
    f_l1 = 1_575_420_000.0
    return delay_l1 * (f_l1 / frequency_hz) ** 2


# ─── Saastamoinen ──────────────────────────────────────────────────────────

def saastamoinen_delay(rx_ecef: np.ndarray, sv_ecef: np.ndarray) -> float:
    """
    Saastamoinen tropospheric delay (m), with default atmosphere
    (P=1013.25 mbar, T=15°C, RH=50%).
    """
    rx_lla = ecef_to_lla(*rx_ecef)
    elev, _ = elevation_azimuth(sv_ecef, rx_ecef)
    if elev < math.radians(3):
        return 0.0  # geometry too poor; let elevation mask handle it
    h = max(rx_lla.alt_m, 0.0)
    # Default atmosphere, simple altitude scaling
    P = 1013.25 * (1 - 2.2557e-5 * h) ** 5.2568
    T = 15.0 - 6.5e-3 * h + 273.15
    e = 6.108 * 0.5 * math.exp((17.15 * (T - 273.15) - 4684.0) / (T - 38.45))

    z = math.pi / 2 - elev
    delay = (0.002277 / math.cos(z)) * (P + (1255.0 / T + 0.05) * e - math.tan(z) ** 2)
    return delay

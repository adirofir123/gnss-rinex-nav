"""
Constants and coordinate / time conversions for GNSS positioning.

All angles are radians, all positions ECEF meters, all times seconds
unless explicitly noted otherwise.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone

import numpy as np

# ─── Physical / WGS-84 constants ────────────────────────────────────────────

SPEED_OF_LIGHT = 299_792_458.0          # m/s
OMEGA_E_DOT    = 7.2921151467e-5        # WGS-84 Earth rotation rate, rad/s
MU_EARTH       = 3.986005e14            # WGS-84 gravitational parameter, m^3/s^2
F_REL          = -4.442807633e-10       # Relativistic clock correction factor, s/sqrt(m)

WGS84_A   = 6_378_137.0                  # Semi-major axis, m
WGS84_F   = 1.0 / 298.257223563          # Flattening
WGS84_E2  = WGS84_F * (2.0 - WGS84_F)    # First eccentricity squared
WGS84_B   = WGS84_A * (1.0 - WGS84_F)    # Semi-minor axis

# GPS time origin: 1980-01-06 00:00:00 UTC
GPS_EPOCH_UTC = datetime(1980, 1, 6, 0, 0, 0, tzinfo=timezone.utc)

# Leap-seconds between GPS time and UTC (GPS - UTC). Update when IERS issues
# a new leap second. As of 2017-01-01 the offset is 18 s and no leap seconds
# have been added since.
GPS_UTC_LEAP_SECONDS = 18


# ─── Coordinate transforms ──────────────────────────────────────────────────

@dataclass
class LLA:
    """Geodetic coordinates (WGS-84)."""
    lat_deg: float
    lon_deg: float
    alt_m: float


def ecef_to_lla(x: float, y: float, z: float) -> LLA:
    """
    Convert ECEF (m) to WGS-84 geodetic latitude/longitude/altitude.

    Uses the closed-form Bowring iteration (a few iterations converge to
    sub-millimeter precision in altitude).
    """
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    if p < 1e-9:
        # Polar singularity — latitude is ±90°, altitude is |z| - b
        lat = math.copysign(math.pi / 2, z)
        alt = abs(z) - WGS84_B
        return LLA(math.degrees(lat), math.degrees(lon), alt)

    # Initial guess
    lat = math.atan2(z, p * (1 - WGS84_E2))
    for _ in range(6):
        sin_lat = math.sin(lat)
        N = WGS84_A / math.sqrt(1 - WGS84_E2 * sin_lat * sin_lat)
        alt = p / math.cos(lat) - N
        lat = math.atan2(z, p * (1 - WGS84_E2 * N / (N + alt)))

    sin_lat = math.sin(lat)
    N = WGS84_A / math.sqrt(1 - WGS84_E2 * sin_lat * sin_lat)
    alt = p / math.cos(lat) - N
    return LLA(math.degrees(lat), math.degrees(lon), alt)


def lla_to_ecef(lat_deg: float, lon_deg: float, alt_m: float) -> tuple[float, float, float]:
    """Convert WGS-84 LLA (deg, deg, m) → ECEF (m)."""
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    N = WGS84_A / math.sqrt(1 - WGS84_E2 * sin_lat * sin_lat)
    x = (N + alt_m) * cos_lat * math.cos(lon)
    y = (N + alt_m) * cos_lat * math.sin(lon)
    z = (N * (1 - WGS84_E2) + alt_m) * sin_lat
    return x, y, z


def ecef_to_enu_rotation(lat_deg: float, lon_deg: float) -> np.ndarray:
    """
    Rotation matrix that maps an ECEF vector (relative to a reference point)
    to a local East/North/Up frame anchored at (lat, lon).

    Used for elevation/azimuth of satellites and for plotting.
    """
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    sl, cl = math.sin(lat), math.cos(lat)
    so, co = math.sin(lon), math.cos(lon)
    return np.array([
        [-so,        co,       0.0],
        [-sl * co,  -sl * so,  cl ],
        [ cl * co,   cl * so,  sl ],
    ])


def elevation_azimuth(sv_ecef: np.ndarray, rx_ecef: np.ndarray) -> tuple[float, float]:
    """
    Elevation and azimuth (radians) of an SV as seen from a receiver in ECEF.
    """
    rx_lla = ecef_to_lla(*rx_ecef)
    R = ecef_to_enu_rotation(rx_lla.lat_deg, rx_lla.lon_deg)
    enu = R @ (sv_ecef - rx_ecef)
    horiz = math.hypot(enu[0], enu[1])
    elev = math.atan2(enu[2], horiz) if horiz > 1e-9 else math.copysign(math.pi / 2, enu[2])
    azim = math.atan2(enu[0], enu[1])
    if azim < 0:
        azim += 2 * math.pi
    return elev, azim


# ─── Time helpers ───────────────────────────────────────────────────────────

def gps_time_from_utc(dt_utc: datetime) -> tuple[int, float]:
    """
    Convert a UTC datetime to (GPS week, GPS seconds-of-week).

    Seconds-of-week is in GPS time (i.e. includes the GPS-UTC leap offset).
    """
    if dt_utc.tzinfo is None:
        dt_utc = dt_utc.replace(tzinfo=timezone.utc)
    delta = (dt_utc - GPS_EPOCH_UTC).total_seconds() + GPS_UTC_LEAP_SECONDS
    week = int(delta // (7 * 86400))
    tow = delta - week * 7 * 86400
    return week, tow


def week_tow_from_calendar(
    year: int, month: int, day: int, hour: int, minute: int, second: float,
    sys_to_gps_offset_s: float = 0.0,
) -> tuple[int, float]:
    """
    Convert a *system-time* calendar timestamp to (GPS week, GPS sec-of-week).

    RINEX 3/4 navigation files record each ephemeris epoch (Toc) in the
    *originating constellation's* system time, not UTC. This function takes
    the raw calendar fields and an offset that maps that system to GPS time:

        sys_to_gps_offset_s = (system_time − GPS_time), in seconds.

    Examples:
        GPS / Galileo / QZSS / NavIC: offset 0   (system_time == GPS_time)
        BeiDou:                       offset -14 (BDT runs 14 s behind GPS)
        GLONASS:                      offset = -GPS_UTC_LEAP_SECONDS  (UTC)

    The RINEX observation file's epoch tags are likewise in the time system
    declared by TIME OF FIRST OBS — typically GPS, in which case the offset
    is zero.

    For GPS time inputs (the common case), this is equivalent to
    `gps_time_from_utc(dt) − GPS_UTC_LEAP_SECONDS` — i.e. we do NOT add the
    leap-second correction, because the input is already in GPS time.
    """
    whole_sec = int(second)
    micro = int(round((second - whole_sec) * 1_000_000))
    if micro >= 1_000_000:
        micro = 0
        whole_sec += 1
    dt = datetime(year, month, day, hour, minute, whole_sec, micro, tzinfo=timezone.utc)
    # Naive seconds since the GPS epoch (treating calendar fields as if they
    # were UTC, since datetime arithmetic is what it is) — then subtract the
    # system-to-GPS offset to express in GPS time.
    delta = (dt - GPS_EPOCH_UTC).total_seconds() - sys_to_gps_offset_s
    week = int(delta // (7 * 86400))
    tow = delta - week * 7 * 86400
    return week, tow


def utc_from_gps_time(week: int, tow: float) -> datetime:
    """Convert GPS (week, seconds-of-week) to a UTC datetime."""
    total_seconds = week * 7 * 86400 + tow - GPS_UTC_LEAP_SECONDS
    return GPS_EPOCH_UTC + timedelta(seconds=total_seconds)

"""Sanity tests for the deterministic pieces. Run with `python -m pytest`."""

import math

import numpy as np

from gnss_nav.geodesy import (
    ecef_to_lla, lla_to_ecef, gps_time_from_utc, utc_from_gps_time,
)
from datetime import datetime, timezone


def test_ecef_lla_roundtrip_jerusalem():
    # Jerusalem-ish
    lat, lon, alt = 31.7683, 35.2137, 754.0
    x, y, z = lla_to_ecef(lat, lon, alt)
    back = ecef_to_lla(x, y, z)
    assert abs(back.lat_deg - lat) < 1e-7
    assert abs(back.lon_deg - lon) < 1e-7
    assert abs(back.alt_m - alt) < 1e-3


def test_ecef_lla_roundtrip_equator():
    lat, lon, alt = 0.0, 0.0, 0.0
    x, y, z = lla_to_ecef(lat, lon, alt)
    back = ecef_to_lla(x, y, z)
    assert abs(back.lat_deg) < 1e-9
    assert abs(back.lon_deg) < 1e-9
    assert abs(back.alt_m) < 1e-3


def test_gps_time_roundtrip():
    dt = datetime(2024, 6, 15, 12, 34, 56, tzinfo=timezone.utc)
    week, tow = gps_time_from_utc(dt)
    back = utc_from_gps_time(week, tow)
    assert back == dt

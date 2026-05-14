"""
RINEX 4.0 parser (observation + navigation).

Design notes
------------
* RINEX 4.0 navigation files use the new record framing where each navigation
  message starts with a `> EPH G PRN  LNAV` (or similar) line. We parse that
  framing, but **also** accept legacy RINEX 3.0x nav files (where the message
  header line is just the SV id + epoch) — Android-derived files in the wild
  vary, and we want the pipeline to work on both.
* Observation files follow the standard RINEX 3.0x epoch layout, which RINEX
  4.0 retains unchanged. The header lists observation types per constellation;
  we record the column index of each obs code we care about (pseudorange,
  Doppler).
* This parser is intentionally minimal: it understands GPS (G), Galileo (E),
  GLONASS (R), BeiDou (C) ephemerides in the standard Keplerian form. GLONASS
  uses a different model (position/velocity/accel state vector) and we parse
  it but flag it — the solver currently only propagates Keplerian SVs.

The parser keeps everything in memory; RINEX files are typically <100 MB.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np

from .geodesy import GPS_UTC_LEAP_SECONDS, week_tow_from_calendar


# ─── Data classes ───────────────────────────────────────────────────────────

@dataclass
class Observation:
    """A single observation epoch for one SV on one constellation."""
    sv: str                  # e.g. "G05", "E11"
    pseudorange: float | None = None   # meters
    doppler: float | None = None       # Hz (signed; positive = approaching)
    carrier_phase: float | None = None # cycles
    snr: float | None = None           # dB-Hz
    frequency_hz: float | None = None  # carrier frequency for Doppler→velocity


@dataclass
class Epoch:
    """All observations at one receiver epoch."""
    time_utc: datetime
    gps_week: int
    gps_tow: float                 # seconds of week
    observations: list[Observation] = field(default_factory=list)


@dataclass
class Ephemeris:
    """
    Broadcast ephemeris (Keplerian) for one SV, one issue of data.

    Fields follow RINEX 3/4 LNAV / FNAV layout for GPS/Galileo/BeiDou.
    Times are GPS-time seconds-of-week.
    """
    sv: str
    toc_week: int
    toc: float                # clock reference time, sec of week
    af0: float                # SV clock bias, s
    af1: float                # SV clock drift, s/s
    af2: float                # SV clock drift rate, s/s^2

    # Orbit
    iode: float
    crs: float; delta_n: float; m0: float
    cuc: float; e: float; cus: float; sqrt_a: float
    toe: float                # orbit reference time, sec of week
    cic: float; omega0: float; cis: float
    i0: float; crc: float; omega: float; omega_dot: float
    idot: float

    # Misc (kept for completeness)
    codes_l2: float = 0.0
    gps_week: int = 0
    l2_pflag: float = 0.0
    sv_accuracy: float = 0.0
    sv_health: float = 0.0
    tgd: float = 0.0
    iodc: float = 0.0
    fit_interval: float = 0.0

    @property
    def a(self) -> float:
        return self.sqrt_a * self.sqrt_a


@dataclass
class IonoParams:
    """Klobuchar ionospheric model coefficients (8 floats)."""
    alpha: tuple[float, float, float, float]
    beta:  tuple[float, float, float, float]


@dataclass
class RinexData:
    obs_epochs: list[Epoch] = field(default_factory=list)
    ephemerides: list[Ephemeris] = field(default_factory=list)
    iono: IonoParams | None = None
    approx_position_ecef: tuple[float, float, float] | None = None
    obs_types: dict[str, list[str]] = field(default_factory=dict)  # per system
    leap_seconds: int | None = None   # from LEAP SECONDS header (RINEX 4: mandatory)


# ─── Helpers ────────────────────────────────────────────────────────────────

# RINEX uses Fortran-style 'D' for exponent in nav files.
_NUM_RE = re.compile(r"[-+]?\d*\.?\d+(?:[DEde][-+]?\d+)?")

def _f(s: str) -> float:
    """Parse a RINEX float (handles 'D' exponent and blank → 0.0)."""
    s = s.strip().replace("D", "E").replace("d", "e")
    if not s:
        return 0.0
    return float(s)


def _parse_obs_header(lines: Iterable[str]) -> tuple[dict[str, list[str]], tuple[float, float, float] | None, IonoParams | None]:
    """
    Walk the observation-file header and extract:
      - obs_types per constellation,
      - APPROX POSITION XYZ if present.
    Stops at END OF HEADER.

    Note: ionospheric coefficients live in the *navigation* file (RINEX 3+),
    so we don't search for them here.
    """
    obs_types: dict[str, list[str]] = {}
    approx_pos: tuple[float, float, float] | None = None
    pending_sys: str | None = None
    pending_count: int = 0
    pending_codes: list[str] = []

    for line in lines:
        label = line[60:].rstrip()
        # Some files have slight column misalignment; the canonical labels
        # are unique enough to match by substring as a fallback.
        if "END OF HEADER" in label:
            break

        if "APPROX POSITION XYZ" in label:
            try:
                x = float(line[0:14]); y = float(line[14:28]); z = float(line[28:42])
                approx_pos = (x, y, z)
            except ValueError:
                pass

        elif "SYS / # / OBS TYPES" in label or "OBS TYPES" in label:
            # Either a new system header line, or a continuation
            if pending_sys is None or (line[0] != " " and pending_count == 0):
                pending_sys = line[0]
                try:
                    pending_count = int(line[3:6])
                except ValueError:
                    pending_count = 0
                pending_codes = []
                # codes start at column 7, 4 chars each (3-char code + space)
                codes_part = line[7:60]
            else:
                codes_part = line[7:60]

            for i in range(0, len(codes_part), 4):
                code = codes_part[i:i + 4].strip()
                if code:
                    pending_codes.append(code)
                if len(pending_codes) == pending_count:
                    break

            if len(pending_codes) >= pending_count and pending_count > 0:
                obs_types[pending_sys] = pending_codes[:pending_count]
                pending_sys = None
                pending_count = 0

    return obs_types, approx_pos, None


def _obs_epoch_time(line: str) -> tuple[datetime, int, float] | None:
    """
    Parse a RINEX 3+ epoch line of form:
        > YYYY MM DD HH MM SS.sssssss  FLAG  NSAT  ...

    The RINEX 4 observation header declares `TIME OF FIRST OBS ... GPS`
    (per Table A2 / §4.1) which means the calendar timestamp on every
    epoch line is already in GPS time, not UTC. We therefore compute the
    GPS week / sec-of-week DIRECTLY from these fields without re-adding
    the GPS-UTC leap-second offset.

    For UTC display in our outputs we subtract GPS_UTC_LEAP_SECONDS at
    the writer level.

    Returns (UTC datetime, GPS week, GPS-sec-of-week).
    """
    if not line.startswith(">"):
        return None
    parts = line[2:].split()
    if len(parts) < 6:
        return None
    try:
        year   = int(parts[0])
        month  = int(parts[1])
        day    = int(parts[2])
        hour   = int(parts[3])
        minute = int(parts[4])
        sec    = float(parts[5])
    except ValueError:
        return None
    # GPS-time → GPS week / sec-of-week directly (no leap offset added).
    week, tow = week_tow_from_calendar(year, month, day, hour, minute, sec,
                                        sys_to_gps_offset_s=0.0)
    # Convert the GPS-time calendar fields to a UTC datetime for output:
    # UTC = GPS_time − leap_seconds.
    from datetime import timedelta
    whole = int(sec)
    micro = int(round((sec - whole) * 1_000_000))
    if micro >= 1_000_000:
        micro = 0
        whole += 1
    dt_gps = datetime(year, month, day, hour, minute, whole, micro, tzinfo=timezone.utc)
    dt_utc = dt_gps - timedelta(seconds=GPS_UTC_LEAP_SECONDS)
    return dt_utc, week, tow


# Frequencies (Hz) for the common Android-recorded signals. Used to convert
# Doppler (Hz) into range-rate (m/s) when needed.
_DEFAULT_FREQ = {
    "G_L1": 1_575_420_000.0,  # GPS L1
    "G_L5": 1_176_450_000.0,
    "E_L1": 1_575_420_000.0,  # Galileo E1
    "E_L5": 1_176_450_000.0,  # Galileo E5a
    "C_L1": 1_575_420_000.0,  # BeiDou B1C (close enough for typical Android)
    "R_L1": 1_602_000_000.0,  # GLONASS L1 nominal (channel-dependent in reality)
}

def _freq_for(obs_code: str, sys: str) -> float | None:
    """Best-effort carrier frequency lookup from a 3-char RINEX obs code."""
    if len(obs_code) < 2:
        return None
    band = obs_code[1]   # 1, 2, 5, 6, 7, 8
    key = f"{sys}_L{band}"
    return _DEFAULT_FREQ.get(key) or _DEFAULT_FREQ.get(f"G_L{band}")


# ─── Observation file parser ────────────────────────────────────────────────

def _parse_obs_file(path: Path) -> tuple[list[Epoch], dict[str, list[str]], tuple[float, float, float] | None]:
    with path.open("r", errors="replace") as fh:
        lines = fh.readlines()

    # Split header / body
    header_end = None
    for i, line in enumerate(lines):
        if "END OF HEADER" in line[60:]:
            header_end = i
            break
    if header_end is None:
        raise ValueError(f"No END OF HEADER in {path}")

    obs_types, approx_pos, _ = _parse_obs_header(lines[: header_end + 1])
    epochs: list[Epoch] = []
    i = header_end + 1
    n = len(lines)

    while i < n:
        line = lines[i]
        if not line.startswith(">"):
            i += 1
            continue
        time_info = _obs_epoch_time(line)
        if time_info is None:
            i += 1
            continue
        dt, week, tow = time_info
        try:
            nsat = int(line[32:35])
        except ValueError:
            i += 1
            continue
        epoch = Epoch(time_utc=dt, gps_week=week, gps_tow=tow)
        i += 1

        for _ in range(nsat):
            if i >= n:
                break
            obs_line = lines[i].rstrip("\n")
            i += 1
            if len(obs_line) < 3:
                continue
            sv = obs_line[0:3]
            sys = sv[0]
            codes = obs_types.get(sys, [])
            # Each obs is 16 chars: F14.3 + 2 flag chars.
            values: list[float | None] = []
            for k in range(len(codes)):
                start = 3 + k * 16
                end = start + 14
                chunk = obs_line[start:end]
                if not chunk.strip():
                    values.append(None)
                else:
                    try:
                        values.append(float(chunk))
                    except ValueError:
                        values.append(None)

            # Pick the first pseudorange (C..) and Doppler (D..) and SNR (S..).
            ob = Observation(sv=sv)
            for code, val in zip(codes, values):
                if val is None or not code:
                    continue
                kind = code[0]
                if kind == "C" and ob.pseudorange is None:
                    ob.pseudorange = val
                    ob.frequency_hz = _freq_for(code, sys)
                elif kind == "D" and ob.doppler is None:
                    ob.doppler = val
                    if ob.frequency_hz is None:
                        ob.frequency_hz = _freq_for(code, sys)
                elif kind == "L" and ob.carrier_phase is None:
                    ob.carrier_phase = val
                elif kind == "S" and ob.snr is None:
                    ob.snr = val

            if ob.pseudorange is not None:
                epoch.observations.append(ob)

        if epoch.observations:
            epochs.append(epoch)

    return epochs, obs_types, approx_pos


# ─── Navigation file parser ─────────────────────────────────────────────────

def _parse_nav_records(lines: list[str]) -> tuple[list[Ephemeris], IonoParams | None, int | None]:
    """
    Parse a RINEX 3 *or* 4 navigation file.

    Strategy
    --------
    RINEX 4 introduced a per-record header line of the form

        > EPH G01 LNAV
        > ION G    LNAV
        > STO E08 IFNV
        > EOP J01 CNVX

    where the second token is the record type (EPH/STO/EOP/ION), the third
    is the constellation + optional PRN, and the fourth is the message type.
    The number of data lines that follow depends on (record_type, message_type)
    — see the appendix tables in the spec. Constants in `_EPH_LINES` below
    encode the line counts we need.

    Legacy RINEX 3 navigation files have no `>` header — each record starts
    directly with the SV/clock line. We treat that as if it had been
    preceded by `> EPH <sys> <legacy_msg_type>` so the same code path
    handles both. Klobuchar iono coefficients in v3 live in the file header
    as `GPSA` / `GPSB` lines; in v4 they live in ION records in the body.
    Both are extracted.

    Only the "Keplerian-orbit" EPH messages (LNAV/INAV/FNAV/D1/D2) are
    propagated by the solver; CNAV/CNV1/CNV2/CNV3 are parsed but skipped
    (their orbital model uses semi-major-axis-rate which we don't yet
    implement). GLONASS FDMA is parsed but not propagated for the same
    reason — it uses a position/velocity/acceleration state vector, not
    Keplerian elements.
    """
    # (record_type, msg_type) → number of data lines that follow the > header
    # (i.e. the SV/clock line PLUS the broadcast-orbit lines).
    _EPH_LINES = {
        # GPS / QZSS / NavIC LNAV: clock + 7 orbit lines = 8 total
        ("EPH", "LNAV"): 8,
        # GPS / QZSS CNAV: clock + 8 orbit lines = 9 total
        ("EPH", "CNAV"): 9,
        # GPS / QZSS CNV2: clock + 9 orbit lines = 10 total
        ("EPH", "CNV2"): 10,
        # Galileo I/FNAV: clock + 7 orbit lines = 8 total
        # (last line has only t_tm — we still consume the full line)
        ("EPH", "INAV"): 8,
        ("EPH", "FNAV"): 8,
        # GLONASS FDMA: clock + 4 orbit lines = 5 total
        ("EPH", "FDMA"): 5,
        # BeiDou D1/D2: clock + 7 orbit lines = 8 total
        ("EPH", "D1"):   8,
        ("EPH", "D2"):   8,
        # BeiDou CNV1: clock + 9 orbit lines = 10 total
        ("EPH", "CNV1"): 10,
        # BeiDou CNV3: clock + 9 orbit lines = 10 total
        ("EPH", "CNV3"): 10,
        # SBAS: clock + 3 orbit lines = 4 total
        ("EPH", "SBAS"): 4,

        # ION records: 1 epoch+α0-α2 + 1 (α3,β0-β2) + 1 (β3,region) = 3 total
        ("ION", "LNAV"): 3,
        ("ION", "D1D2"): 3,
        ("ION", "CNVX"): 3,
        # NeQuick-G (Galileo): epoch+ai0-ai2 + disturbance flags = 2 total
        ("ION", "IFNV"): 2,

        # STO records: 1 line of corr/sbas/utc IDs + 1 STO message line = 2 total
        ("STO", "LNAV"): 2, ("STO", "FDMA"): 2, ("STO", "IFNV"): 2,
        ("STO", "D1D2"): 2, ("STO", "SBAS"): 2, ("STO", "CNVX"): 2,

        # EOP records: 1 epoch + 1 EOP-line-1 + 1 EOP-line-2 = 3 total
        ("EOP", "LNAV"): 3, ("EOP", "CNVX"): 3,
    }
    # Message types whose ephemeris we can actually propagate with our
    # Keplerian solver. Everything else is parsed and discarded.
    _SUPPORTED_KEPLER_MSG = {"LNAV", "INAV", "FNAV", "D1", "D2"}

    # ─── Header parsing (legacy iono + leap seconds) ────────────────────────
    iono_alpha: list[float] | None = None
    iono_beta:  list[float] | None = None
    leap_seconds: int | None = None

    i = 0
    while i < len(lines):
        line = lines[i]
        label = line[60:].rstrip()
        if "END OF HEADER" in label:
            i += 1
            break
        # Legacy RINEX 3 iono header lines. RINEX 4 has removed these, but
        # we still parse them for backward compatibility.
        if "IONOSPHERIC CORR" in label:
            tag = line[0:4].strip()
            vals = [
                _f(line[5:17]),
                _f(line[17:29]),
                _f(line[29:41]),
                _f(line[41:53]),
            ]
            if tag in ("GPSA", "GAL"):
                iono_alpha = vals
            elif tag == "GPSB":
                iono_beta = vals
        elif "LEAP SECONDS" in label:
            # I6,I6,I6,I6,A3 — current leap seconds is the first field.
            try:
                leap_seconds = int(line[0:6])
            except ValueError:
                pass
        i += 1

    ephemerides: list[Ephemeris] = []

    # ─── Body parsing ───────────────────────────────────────────────────────
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1
            continue

        # Determine record_type + msg_type for the upcoming record.
        if line.startswith(">"):
            parts = line[1:].split()
            if len(parts) < 3:
                i += 1
                continue
            record_type = parts[0].upper()           # EPH / STO / EOP / ION
            src_field   = parts[1]                    # e.g. "G01" or "G"
            msg_type    = parts[2].upper()            # LNAV / INAV / ...
            i += 1   # advance past the header line; the data block follows
        else:
            # Legacy v3 record: no > line. Infer record_type=EPH and a
            # reasonable msg_type from the constellation.
            if line[0] not in "GERCJIS":
                i += 1
                continue
            record_type = "EPH"
            sys_char = line[0]
            msg_type = {
                "G": "LNAV", "R": "FDMA", "E": "INAV",
                "C": "D1", "J": "LNAV", "I": "LNAV", "S": "SBAS",
            }[sys_char]
            src_field = ""   # PRN is in the data line itself

        n_data_lines = _EPH_LINES.get((record_type, msg_type))
        if n_data_lines is None:
            # Unknown record — skip the > line we already consumed and try
            # to resync by hunting for the next > or recognizable SV line.
            continue

        # Make sure we have enough lines left.
        if i + n_data_lines > len(lines):
            break

        block = [lines[i + k].rstrip("\n") for k in range(n_data_lines)]
        i += n_data_lines

        # ─── ION (Klobuchar) ────────────────────────────────────────────
        if record_type == "ION" and msg_type in ("LNAV", "D1D2", "CNVX"):
            try:
                # Line 0: 4X,I4,5(1X,I2.2) epoch + α0-α2 (3 floats)
                # Line 1: 4X,4D19.12  (α3, β0, β1, β2)
                # Line 2: 4X,2D19.12  (β3, region)
                ln0 = block[0]
                a0 = _f(ln0[23:42]); a1 = _f(ln0[42:61]); a2 = _f(ln0[61:80])
                ln1 = block[1]
                a3 = _f(ln1[4:23]); b0 = _f(ln1[23:42])
                b1 = _f(ln1[42:61]); b2 = _f(ln1[61:80])
                ln2 = block[2]
                b3 = _f(ln2[4:23])
                # Only keep the first Klobuchar set we encounter (typically
                # GPS LNAV — the most widely tracked).
                if iono_alpha is None:
                    iono_alpha = [a0, a1, a2, a3]
                    iono_beta  = [b0, b1, b2, b3]
            except (IndexError, ValueError):
                pass
            continue

        # ─── EPH (only Keplerian flavors propagated) ───────────────────
        if record_type != "EPH" or msg_type not in _SUPPORTED_KEPLER_MSG:
            continue   # STO/EOP and non-Keplerian EPH are skipped

        # Parse the SV/clock line and the broadcast-orbit lines.
        clock_line = block[0]
        if len(clock_line) < 23 or clock_line[0] not in "GERCJIS":
            continue
        try:
            sv = clock_line[0:3]
            year   = int(clock_line[4:8])
            month  = int(clock_line[9:11])
            day    = int(clock_line[12:14])
            hour   = int(clock_line[15:17])
            minute = int(clock_line[18:20])
            sec    = int(clock_line[21:23])
            af0    = _f(clock_line[23:42])
            af1    = _f(clock_line[42:61])
            af2    = _f(clock_line[61:80])
        except (ValueError, IndexError):
            continue

        # The calendar Toc in the nav record is in the *constellation's* own
        # system time (per RINEX 4 §4.1 and §5.4). To express it in GPS-time
        # (which all our solver math uses), we need the appropriate offset:
        #   sys_to_gps_offset = (system_time − GPS_time), seconds.
        # GPS/GAL/QZS/IRN are aligned with GPS time → offset 0.
        # BDS runs 14 s behind GPS time (BDT epoch is 2006-01-01 UTC, vs
        #   GPS's 1980-01-06 UTC, which is 14 leap-second steps later).
        # GLONASS records in UTC.
        sys_char = sv[0]
        if sys_char == "C":              # BeiDou
            sys_to_gps_offset = -14.0
        elif sys_char == "R":            # GLONASS
            sys_to_gps_offset = float(-GPS_UTC_LEAP_SECONDS)
        else:                            # G / E / J / I — same as GPS time
            sys_to_gps_offset = 0.0

        toc_week, toc = week_tow_from_calendar(
            year, month, day, hour, minute, float(sec),
            sys_to_gps_offset_s=sys_to_gps_offset,
        )

        # Collect all numeric fields from the orbit lines (4 per line,
        # column-aligned per RINEX 3+ spec: 4X,4D19.12).
        try:
            fields: list[float] = []
            for ln in block[1:]:
                col_chunks = [ln[4:23], ln[23:42], ln[42:61], ln[61:80]]
                # Try column-based parsing first; if it fails, fall back to
                # whitespace split (some non-conformant files lack proper
                # leading indent).
                ok = True
                for c in col_chunks:
                    s = c.strip()
                    if s and not _NUM_RE.fullmatch(s.replace("D", "E").replace("d", "e")):
                        ok = False
                        break
                if ok:
                    for c in col_chunks:
                        fields.append(_f(c) if c.strip() else 0.0)
                else:
                    parts = ln.replace("D", "E").replace("d", "e").split()
                    while len(parts) < 4:
                        parts.append("0.0")
                    for p in parts[:4]:
                        try:
                            fields.append(float(p))
                        except ValueError:
                            fields.append(0.0)
        except (IndexError, ValueError):
            continue

        # All Keplerian flavors share the same first 17 fields (the orbit
        # parameters proper). Layout, with field index:
        #   0  IODE       1  Crs     2  Δn      3  M0
        #   4  Cuc        5  e       6  Cus     7  √a
        #   8  toe        9  Cic    10  Ω0     11  Cis
        #  12  i0        13  Crc    14  ω      15  Ω̇
        #  16  IDOT      17  codes_l2 / data_src / SatType   18  GPS_week
        #  19  L2_pflag / spare    20  SV_acc  21  SV_health
        #  22  TGD       23  IODC
        try:
            eph = Ephemeris(
                sv=sv, toc_week=toc_week, toc=toc,
                af0=af0, af1=af1, af2=af2,
                iode=fields[0], crs=fields[1], delta_n=fields[2], m0=fields[3],
                cuc=fields[4], e=fields[5], cus=fields[6], sqrt_a=fields[7],
                toe=fields[8], cic=fields[9], omega0=fields[10], cis=fields[11],
                i0=fields[12], crc=fields[13], omega=fields[14], omega_dot=fields[15],
                idot=fields[16],
                codes_l2=fields[17] if len(fields) > 17 else 0.0,
                gps_week=int(fields[18]) if len(fields) > 18 else toc_week,
                l2_pflag=fields[19] if len(fields) > 19 else 0.0,
                sv_accuracy=fields[20] if len(fields) > 20 else 0.0,
                sv_health=fields[21] if len(fields) > 21 else 0.0,
                tgd=fields[22] if len(fields) > 22 else 0.0,
                iodc=fields[23] if len(fields) > 23 else 0.0,
            )
            ephemerides.append(eph)
        except IndexError:
            pass

    iono = None
    if iono_alpha and iono_beta:
        iono = IonoParams(tuple(iono_alpha), tuple(iono_beta))
    return ephemerides, iono, leap_seconds


# ─── Public API ─────────────────────────────────────────────────────────────

def parse_rinex(path: str | Path) -> RinexData:
    """
    Parse a RINEX 4.0 file. Accepts:
      * a single combined file containing both obs + nav records, OR
      * a path whose sibling has the matching nav file (same stem, .nav / .rnx).

    Returns a RinexData populated with whatever was found.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    data = RinexData()

    # Detect file type from the header line "RINEX VERSION / TYPE"
    with path.open("r", errors="replace") as fh:
        head = [next(fh) for _ in range(5)]
    file_type = ""
    for line in head:
        if "RINEX VERSION / TYPE" in line[60:]:
            file_type = line[20:21].strip().upper()  # O, N, M (merged)
            break

    if file_type == "O":
        epochs, obs_types, approx = _parse_obs_file(path)
        data.obs_epochs = epochs
        data.obs_types = obs_types
        data.approx_position_ecef = approx
        # Look for a sibling navigation file
        for ext in (".nav", ".rnx", ".n", ".N", ".P"):
            candidate = path.with_suffix(ext)
            if candidate.exists() and candidate != path:
                with candidate.open("r", errors="replace") as fh:
                    nav_lines = fh.readlines()
                eph, iono, leap = _parse_nav_records(nav_lines)
                data.ephemerides.extend(eph)
                if iono and data.iono is None:
                    data.iono = iono
                if leap is not None:
                    data.leap_seconds = leap
                break
    elif file_type == "N":
        with path.open("r", errors="replace") as fh:
            nav_lines = fh.readlines()
        eph, iono, leap = _parse_nav_records(nav_lines)
        data.ephemerides = eph
        data.iono = iono
        data.leap_seconds = leap
    elif file_type == "M":
        # Merged: parse twice — once as obs, once as nav
        epochs, obs_types, approx = _parse_obs_file(path)
        data.obs_epochs = epochs
        data.obs_types = obs_types
        data.approx_position_ecef = approx
        with path.open("r", errors="replace") as fh:
            nav_lines = fh.readlines()
        eph, iono, leap = _parse_nav_records(nav_lines)
        data.ephemerides = eph
        data.iono = iono
        data.leap_seconds = leap
    else:
        raise ValueError(f"Unknown RINEX file type: {file_type!r}")

    return data

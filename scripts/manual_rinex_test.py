"""
Build a minimal RINEX 3.04-style obs file + matching nav file from synthetic
data, run the parser, and check we get sane epochs and ephemerides out.

The point is to exercise the parser end-to-end, not to produce a physically
realistic trajectory — the solver math is already validated by manual_smoke.py.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pathlib import Path
from gnss_nav.rinex import parse_rinex

tmpdir = Path("/tmp/rinex_test")
tmpdir.mkdir(exist_ok=True)

OBS = """     3.04           OBSERVATION DATA    M (MIXED)            RINEX VERSION / TYPE
test                                                        PGM / RUN BY / DATE
TEST                                                        MARKER NAME
  4641600.0000  4356600.0000  3334400.0000                  APPROX POSITION XYZ
G    2 C1C D1C                                              SYS / # / OBS TYPES
E    2 C1C D1C                                              SYS / # / OBS TYPES
  2024     6    15     0     0    0.0000000     GPS         TIME OF FIRST OBS
                                                            END OF HEADER
> 2024 06 15 00 00 00.0000000  0  3
G01  22500000.000          1234.000
G02  21000000.000           -42.000
G03  23800000.000           567.000
> 2024 06 15 00 00 01.0000000  0  3
G01  22500030.000          1234.000
G02  21000005.000           -42.000
G03  23800020.000           567.000
"""

NAV = """     3.04           N: GNSS NAV DATA    M (MIXED)           RINEX VERSION / TYPE
test                                                        PGM / RUN BY / DATE
GPSA   1.3970D-08  0.0000D+00 -5.9605D-08  5.9605D-08       IONOSPHERIC CORR
GPSB   1.1059D+05  0.0000D+00 -2.6214D+05  1.9661D+05       IONOSPHERIC CORR
                                                            END OF HEADER
G01 2024 06 15 00 00 00 0.000000000000D+00 0.000000000000D+00 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 4.272915225097D-09 0.000000000000D+00
     0.000000000000D+00 5.000000000000D-03 0.000000000000D+00 5.153720703125D+03
     5.184000000000D+05 0.000000000000D+00 0.000000000000D+00 0.000000000000D+00
     9.500000000000D-01 0.000000000000D+00 0.000000000000D+00 -7.500000000000D-09
     0.000000000000D+00 0.000000000000D+00 2.318000000000D+03 0.000000000000D+00
     2.000000000000D+00 0.000000000000D+00 -5.000000000000D-09 0.000000000000D+00
     5.184000000000D+05 4.000000000000D+00 0.000000000000D+00 0.000000000000D+00
G02 2024 06 15 00 00 00 0.000000000000D+00 0.000000000000D+00 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 4.272915225097D-09 1.570796326790D+00
     0.000000000000D+00 5.000000000000D-03 0.000000000000D+00 5.153720703125D+03
     5.184000000000D+05 0.000000000000D+00 7.853981633970D-01 0.000000000000D+00
     9.500000000000D-01 0.000000000000D+00 0.000000000000D+00 -7.500000000000D-09
     0.000000000000D+00 0.000000000000D+00 2.318000000000D+03 0.000000000000D+00
     2.000000000000D+00 0.000000000000D+00 -5.000000000000D-09 0.000000000000D+00
     5.184000000000D+05 4.000000000000D+00 0.000000000000D+00 0.000000000000D+00
G03 2024 06 15 00 00 00 0.000000000000D+00 0.000000000000D+00 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 4.272915225097D-09 3.141592653590D+00
     0.000000000000D+00 5.000000000000D-03 0.000000000000D+00 5.153720703125D+03
     5.184000000000D+05 0.000000000000D+00 1.570796326790D+00 0.000000000000D+00
     9.500000000000D-01 0.000000000000D+00 0.000000000000D+00 -7.500000000000D-09
     0.000000000000D+00 0.000000000000D+00 2.318000000000D+03 0.000000000000D+00
     2.000000000000D+00 0.000000000000D+00 -5.000000000000D-09 0.000000000000D+00
     5.184000000000D+05 4.000000000000D+00 0.000000000000D+00 0.000000000000D+00
"""

obs_path = tmpdir / "synthetic.obs"
nav_path = tmpdir / "synthetic.nav"
obs_path.write_text(OBS)
nav_path.write_text(NAV)

data = parse_rinex(obs_path)
print(f"Epochs parsed: {len(data.obs_epochs)}")
for ep in data.obs_epochs:
    print(f"  {ep.time_utc.isoformat()}: {len(ep.observations)} SVs   first={ep.observations[0].sv} PR={ep.observations[0].pseudorange:.2f} Dop={ep.observations[0].doppler}")
print(f"Ephemerides parsed: {len(data.ephemerides)}")
for e in data.ephemerides:
    print(f"  {e.sv}: toe={e.toe} sqrt_a={e.sqrt_a:.3f}")
print(f"Iono coefficients: {data.iono is not None}")
print(f"Approx position: {data.approx_position_ecef}")
print(f"Obs types: {data.obs_types}")

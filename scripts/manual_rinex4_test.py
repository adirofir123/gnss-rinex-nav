"""
Verify the parser against real RINEX 4.00 example records taken directly
from the spec (Tables A8, A12, A14). Includes:
  - v4 nav header (LEAP SECONDS mandatory, no IONOSPHERIC CORR)
  - > EPH G04 LNAV   record (Table A12)
  - > EPH E12 INAV   record (Table A14)
  - > EPH G04 CNAV   record (Table A12) — should be PARSED+SKIPPED
  - > ION G    LNAV  record (RINEX 4 way of conveying Klobuchar params)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pathlib import Path
from gnss_nav.rinex import parse_rinex

# Header + body from the spec's examples. The 60-col header labels are
# right-padded with spaces to land at column 60.
NAV_V4 = """     4.00           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
genericSW           User                20210205 000517 UTC PGM / RUN BY / DATE
    18                                                      LEAP SECONDS
                                                            END OF HEADER
> EPH G04 LNAV
G04 2019 03 14 04 00 00 1.330170780420e-04 7.275957614183e-12 0.000000000000e+00
     9.800000000000e+01-1.718750000000e+00 4.639836124941e-09 2.148941747752e+00
    -1.881271600723e-07 3.355251392350e-04 8.245930075645e-06 5.153800453186e+03
     3.600000000000e+05-1.676380634308e-08 5.171400020311e-01 1.490116119385e-08
     9.601921900531e-01 2.187187500000e+02-1.736906885738e+00-8.044977962767e-09
    -2.932264997750e-10 1.000000000000e+00 2.044000000000e+03 0.000000000000e+00
     4.000000000000e+00 6.300000000000e+01-8.847564458847e-09 8.660000000000e+02
     3.553500000000e+05 4.000000000000e+00
> EPH G04 CNAV
G04 2019 03 14 03 30 00 1.330042141490e-04 7.226219622680e-12 0.000000000000e+00
     2.001762390137e-03 6.914062500000e-01 4.625906973308e-09 1.887277537485e+00
     1.024454832077e-08 3.348654136062e-04 8.376315236092e-06 5.153800325291e+03
     2.412000000000e+05-4.656612873077e-09 5.171544951605e-01 2.328306436539e-08
     9.601927657114e-01 2.174140625000e+02-1.737767543851e+00-8.034028170143e-09
    -2.950122884460e-10-1.312310522376e-14-2.000000000000e+00 2.000000000000e+00
     0.000000000000e+00 7.000000000000e+00-8.789356797934e-09 5.000000000000e+00
    -5.820766091347e-10-6.606569513679e-09-1.178705133498e-08-1.178705133498e-08
     3.558540000000e+05 2.044000000000e+03
> EPH E12 INAV
E12 2020 09 15 00 40 00 5.605182959698e-03-1.881517164293e-11 0.000000000000e+00
     3.600000000000e+01 1.090625000000e+02 2.811188525857e-09-2.481435854929e+00
     5.209818482399e-06 1.468013506383e-04 1.532956957817e-06 5.440609727859e+03
     1.752000000000e+05-1.676380634308e-08 8.103706855689e-01 7.450580596924e-09
     9.891660140720e-01 3.219375000000e+02 5.171049929386e-01-5.815956543649e-09
     2.982267080537e-10 5.170000000000e+02 2.123000000000e+03
     3.120000000000e+00 0.000000000000e+00-1.303851604462e-08-1.280568540096e-08
     1.764340000000e+05
> ION G    LNAV
    2019 03 14 04 00 00 1.397000000000e-08 0.000000000000e+00-5.960000000000e-08
     5.960000000000e-08 1.106000000000e+05 0.000000000000e+00-2.621000000000e+05
     1.966000000000e+05 0.000000000000e+00
"""

tmpdir = Path("/tmp/rinex4_test")
tmpdir.mkdir(exist_ok=True)
nav_path = tmpdir / "test.rnx"
nav_path.write_text(NAV_V4)

data = parse_rinex(nav_path)

print(f"Ephemerides parsed: {len(data.ephemerides)}")
for e in data.ephemerides:
    print(f"  {e.sv}: toe={e.toe:.0f}  sqrt_a={e.sqrt_a:.3f}  e={e.e:.4e}  M0={e.m0:.4f}")

print(f"\nIono coefficients: {data.iono is not None}")
if data.iono:
    print(f"  alpha = {data.iono.alpha}")
    print(f"  beta  = {data.iono.beta}")

# Validate the GPS LNAV record against known spec values
gps = [e for e in data.ephemerides if e.sv == "G04"]
assert len(gps) == 1, f"Expected 1 GPS LNAV record (CNAV should be skipped), got {len(gps)}"
g = gps[0]
assert abs(g.sqrt_a - 5.153800453186e+03) < 1e-6, f"sqrt_a mismatch: {g.sqrt_a}"
assert abs(g.toe - 360000.0) < 1.0, f"toe mismatch: {g.toe}"
assert abs(g.crs - (-1.718750000000e+00)) < 1e-9, f"Crs mismatch: {g.crs}"
assert abs(g.af0 - 1.330170780420e-04) < 1e-15, f"af0 mismatch: {g.af0}"
print("\n✓ GPS LNAV ephemeris parsed correctly (CNAV correctly skipped)")

gal = [e for e in data.ephemerides if e.sv == "E12"]
assert len(gal) == 1, f"Expected 1 Galileo INAV record, got {len(gal)}"
print("✓ Galileo INAV ephemeris parsed correctly")

assert data.iono is not None, "Expected ION record to provide iono coefficients"
assert abs(data.iono.alpha[0] - 1.397e-08) < 1e-15, f"alpha0 mismatch: {data.iono.alpha[0]}"
assert abs(data.iono.beta[0] - 1.106e+05) < 1e-3, f"beta0 mismatch: {data.iono.beta[0]}"
print("✓ ION Klobuchar record parsed correctly (RINEX 4 style)")

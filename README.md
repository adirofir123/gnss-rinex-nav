# GNSS RINEX 4.0 → KML/CSV Navigation

An offline GNSS positioning algorithm that takes a **RINEX 4.0** observation/navigation file
(typically derived from Android raw GNSS measurements) and produces a **1 Hz path**
with 3D position, velocity, and UTC time as **KML + CSV** output.

The solution uses **only the RINEX file** as input — NMEA/TXT reference files are not consumed
(they exist purely for validation).

---

## Quickstart

```bash
# 1. Clone & install
git clone <your-fork-url> gnss-rinex-nav
cd gnss-rinex-nav
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# 2. Run on a RINEX 4.0 file
python -m gnss_nav.cli path/to/input.rnx --out-dir ./out

# Output:
#   ./out/trajectory.csv
#   ./out/trajectory.kml
```

Optional flags:

| Flag | Default | Meaning |
| --- | --- | --- |
| `--elevation-mask DEG` | `10` | Discard satellites below this elevation |
| `--max-iter N` | `10` | Newton iterations for the least-squares solver |
| `--systems LIST` | `G,E` | Constellations to use: G=GPS, E=Galileo, R=GLONASS, C=BeiDou |
| `--rate HZ` | `1` | Output rate (input epochs are resampled / decimated) |
| `--no-tropo` | off | Disable Saastamoinen tropospheric correction |
| `--no-iono` | off | Disable Klobuchar ionospheric correction |

---

## Algorithm overview

The pipeline mirrors the textbook GNSS single-point-positioning (SPP) flow:

```
RINEX 4.0 file
   │
   ▼
[1] Parser  ──►  observation epochs  (pseudoranges ρ, Doppler, time)
                 navigation messages (Keplerian ephemerides per SV)
                 header metadata     (iono coefficients, leap seconds)
   │
   ▼
[2] For each epoch t (1 Hz):
        ├─ Pick ephemeris closest in time per visible SV
        ├─ Compute SV position at *transmission* time  (Kepler propagation,
        │     accounting for signal travel time and Earth rotation)
        ├─ Apply SV clock correction (a0,a1,a2 + relativistic Δtᵣ)
        ├─ Apply tropo (Saastamoinen) + iono (Klobuchar) corrections
        ├─ Solve least-squares for receiver (X, Y, Z, cdt) using ≥4 SVs
        └─ Estimate velocity from Doppler measurements (preferred)
                                                          or position difference (fallback)
   │
   ▼
[3] ECEF → WGS-84 LLA, GPS time → UTC
   │
   ▼
[4] Writers: trajectory.csv  +  trajectory.kml
```

The math for step [2] follows the standard references:

- **Pseudorange equation:** ρ = ‖r_sv − r_rx‖ + c·(dt_rx − dt_sv) + I + T + ε
- **Linearization** around a guess `x₀`, then Newton-Raphson updates.
- **Reference example (4-SV closed-form):** https://mason.gmu.edu/~treid5/Math447/GPSEquations/
- **Spreadsheet walkthrough:** see the assignment link.
- **RINEX 4.0 spec:** https://files.igs.org/pub/data/format/rinex_4.00.pdf

---

## Repository layout

```
gnss-rinex-nav/
├── README.md
├── requirements.txt
├── gnss_nav/
│   ├── __init__.py
│   ├── cli.py              # Command-line entry point
│   ├── rinex.py            # RINEX 4.0 parser (obs + nav)
│   ├── ephemeris.py        # SV position/clock from broadcast ephemeris
│   ├── corrections.py      # Klobuchar iono + Saastamoinen tropo
│   ├── solver.py           # Weighted least-squares SPP + Doppler velocity
│   ├── geodesy.py          # ECEF↔LLA, time systems, constants
│   ├── trajectory.py       # Epoch loop / pipeline orchestrator
│   └── writers.py          # CSV + KML output
├── examples/
│   └── run_example.sh
├── tests/
│   └── test_geodesy.py     # Sanity tests for the deterministic bits
└── docs/
    └── ALGORITHM.md        # Detailed math derivation
```

---

## RINEX 4.0 spec conformance

The parser was verified against the worked examples in the RINEX 4.00
specification (December 2021):

- **Navigation record framing** (§5.4.1): each record starts with
  `> <REC_TYPE> <SRC> <MSG_TYPE>` where `REC_TYPE ∈ {EPH, STO, EOP, ION}`.
- **Keplerian EPH messages** propagated: GPS LNAV (Table A9), Galileo
  INAV/FNAV (Table A13), BeiDou D1/D2 (Table A21), QZSS LNAV (Table A17),
  NavIC LNAV (Table A28). Each has clock + 7 orbit lines = 8 total
  (NavIC = 7), all using `4X,4D19.12` columns.
- **Non-Keplerian EPH messages** (GPS/QZSS CNAV/CNV2, BeiDou CNV1/CNV2/CNV3,
  GLONASS FDMA) are parsed and **skipped** by the solver — they need orbital
  models (semi-major-axis-rate or numerical state-vector integration) we
  don't yet implement. The line counts are still consumed so the parser
  stays in sync across mixed-message files.
- **ION Klobuchar records** (Table A32, §5.4.11): 3 lines containing the
  8 α/β coefficients are parsed directly. The legacy v3 `IONOSPHERIC CORR`
  header line is still accepted as a fallback.
- **Mandatory `LEAP SECONDS` header** (§5.4.1, p. 8): parsed and stored.

If you give the pipeline a RINEX 4 file that includes only CNAV/CNV2
messages for some SVs, those SVs will be silently dropped from the
solution. The fix is to extend `compute_sv_state` with the CNAV
semi-major-axis-rate model (IS-GPS-200, §30.3.3.1.1). Email me if you
hit this and I'll add it.

---

## Verification

The output CSV/KML can be compared visually (Google Earth) and numerically
against the reference NMEA/TXT files distributed alongside each RINEX. Typical
SPP accuracy with broadcast ephemeris on a phone-grade chipset:

- Horizontal: **3–10 m** (open sky)
- Vertical: **5–15 m**
- Velocity (from Doppler): **0.1–0.5 m/s**

Outliers above this are usually multipath, low-elevation SVs, or sparse epochs;
tune `--elevation-mask` accordingly.

---

## Bonus: spoofing recording

If a recording from the lab's Android device contains spoofed measurements, run
the pipeline as usual — the residuals (printed when `--verbose`) will be the
diagnostic signal. Genuine SPP residuals typically sit at a few meters; spoofed
scenes often show clean, near-zero residuals despite an implausible trajectory.



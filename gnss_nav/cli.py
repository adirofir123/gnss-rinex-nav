"""
Command-line entry point.

Usage:
    python -m gnss_nav.cli INPUT.rnx --out-dir ./out
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .rinex import parse_rinex
from .trajectory import Settings, run_trajectory
from .writers import write_csv, write_kml


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="GNSS RINEX 4.0 → KML/CSV navigation")
    p.add_argument("input", type=Path, help="Path to RINEX 4.0 file (.rnx / .obs / .nav)")
    p.add_argument("--out-dir", type=Path, default=Path("./out"), help="Where to put output files")
    p.add_argument("--elevation-mask", type=float, default=10.0, help="Elevation mask in degrees (default 10)")
    p.add_argument("--max-iter", type=int, default=10, help="Newton iterations per epoch (default 10)")
    p.add_argument("--systems", type=str, default="G,E", help="Constellations: G,E,R,C (default G,E)")
    p.add_argument("--rate", type=float, default=1.0, help="Output rate in Hz (default 1)")
    p.add_argument("--no-iono", action="store_true", help="Disable Klobuchar iono correction")
    p.add_argument("--no-tropo", action="store_true", help="Disable Saastamoinen tropo correction")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv if argv is not None else sys.argv[1:])
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    log = logging.getLogger("gnss_nav")

    log.info("Parsing RINEX: %s", args.input)
    data = parse_rinex(args.input)
    log.info("  observation epochs: %d", len(data.obs_epochs))
    log.info("  ephemerides:       %d", len(data.ephemerides))
    log.info("  iono coefficients: %s", "yes" if data.iono else "no")
    if data.approx_position_ecef:
        log.info("  approx pos (ECEF): %.1f %.1f %.1f", *data.approx_position_ecef)

    if not data.obs_epochs or not data.ephemerides:
        log.error("Need both observations and ephemerides. Aborting.")
        return 2

    settings = Settings(
        elevation_mask_deg=args.elevation_mask,
        max_iter=args.max_iter,
        systems=tuple(s.strip().upper() for s in args.systems.split(",") if s.strip()),
        use_iono=not args.no_iono,
        use_tropo=not args.no_tropo,
        rate_hz=args.rate,
    )

    log.info("Running trajectory…")
    fixes = run_trajectory(data, settings)
    log.info("  produced %d fixes", len(fixes))
    if not fixes:
        log.error("No fixes produced — check elevation mask, systems, or input quality.")
        return 3

    out_dir: Path = args.out_dir
    csv_path = out_dir / "trajectory.csv"
    kml_path = out_dir / "trajectory.kml"
    write_csv(fixes, csv_path)
    write_kml(fixes, kml_path, name=args.input.stem)
    log.info("Wrote:\n  %s\n  %s", csv_path, kml_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

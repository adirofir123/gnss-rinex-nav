"""
Output writers: CSV and KML.

CSV columns:
  utc, lat_deg, lon_deg, alt_m, x_ecef, y_ecef, z_ecef,
  vx_ecef, vy_ecef, vz_ecef, speed_mps, num_sats, gdop, clock_bias_m
"""

from __future__ import annotations

import csv
from pathlib import Path
from xml.sax.saxutils import escape

from .trajectory import TrajectoryFix


def write_csv(fixes: list[TrajectoryFix], path: Path | str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "utc", "lat_deg", "lon_deg", "alt_m",
            "x_ecef", "y_ecef", "z_ecef",
            "vx_ecef", "vy_ecef", "vz_ecef", "speed_mps",
            "num_sats", "gdop", "clock_bias_m",
        ])
        for f in fixes:
            w.writerow([
                f.time_utc.isoformat().replace("+00:00", "Z"),
                f"{f.lat_deg:.9f}", f"{f.lon_deg:.9f}", f"{f.alt_m:.3f}",
                f"{f.x_ecef:.3f}", f"{f.y_ecef:.3f}", f"{f.z_ecef:.3f}",
                f"{f.vx_ecef:.4f}", f"{f.vy_ecef:.4f}", f"{f.vz_ecef:.4f}",
                f"{f.speed_mps:.4f}",
                f.num_sats, f"{f.gdop:.3f}", f"{f.clock_bias_m:.3f}",
            ])


def write_kml(fixes: list[TrajectoryFix], path: Path | str, name: str = "GNSS Trajectory") -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    coords = "\n".join(
        f"          {f.lon_deg:.9f},{f.lat_deg:.9f},{f.alt_m:.3f}"
        for f in fixes
    )

    # Per-epoch placemarks with timestamp — lets Google Earth animate the path.
    placemarks = []
    for f in fixes:
        when = f.time_utc.isoformat().replace("+00:00", "Z")
        placemarks.append(f"""    <Placemark>
      <TimeStamp><when>{when}</when></TimeStamp>
      <styleUrl>#fixStyle</styleUrl>
      <Point><coordinates>{f.lon_deg:.9f},{f.lat_deg:.9f},{f.alt_m:.3f}</coordinates></Point>
      <description><![CDATA[
        sats: {f.num_sats}, gdop: {f.gdop:.2f}, speed: {f.speed_mps:.2f} m/s
      ]]></description>
    </Placemark>""")
    placemarks_xml = "\n".join(placemarks)

    kml = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>{escape(name)}</name>
    <Style id="pathStyle">
      <LineStyle><color>ff00aaff</color><width>3</width></LineStyle>
    </Style>
    <Style id="fixStyle">
      <IconStyle>
        <scale>0.4</scale>
        <Icon><href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href></Icon>
      </IconStyle>
    </Style>
    <Placemark>
      <name>Path</name>
      <styleUrl>#pathStyle</styleUrl>
      <LineString>
        <tessellate>1</tessellate>
        <altitudeMode>absolute</altitudeMode>
        <coordinates>
{coords}
        </coordinates>
      </LineString>
    </Placemark>
{placemarks_xml}
  </Document>
</kml>
"""
    path.write_text(kml, encoding="utf-8")

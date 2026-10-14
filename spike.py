#!/usr/bin/env python3
"""Feasibility spike: which known asteroids fall inside a SPHEREx image?

Metadata only, downloads no images. Steps:
  1. IRSA SIA: list SPHEREx QR2 images covering a sky point.
  2. JPL Horizons: SPHEREx's geocentric position at the exposure time.
  3. JPL sb_ident: known asteroids predicted in that field, seen from SPHEREx.

usage: python3 spike.py [RA_DEG] [DEC_DEG] [DETECTOR 1-6]
"""
import csv, io, json, sys, urllib.parse, urllib.request

SPHEREX = "-163182"  # Horizons spacecraft ID


def get(url, params):
    return urllib.request.urlopen(url + "?" + urllib.parse.urlencode(params), timeout=120).read().decode()


def images(ra, dec):
    text = get("https://irsa.ipac.caltech.edu/SIA", {
        "COLLECTION": "spherex_qr2", "POS": f"circle {ra} {dec} 0.01", "RESPONSEFORMAT": "CSV"})
    return list(csv.DictReader(io.StringIO(text)))


def spherex_xyz(mjd):
    """Geocentric ICRF position (km) and velocity (km/s) of SPHEREx at an MJD (UTC)."""
    res = json.loads(get("https://ssd.jpl.nasa.gov/api/horizons.api", {
        "format": "json", "COMMAND": f"'{SPHEREX}'", "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@399'", "REF_PLANE": "'FRAME'", "VEC_TABLE": "'2'",
        "OUT_UNITS": "'KM-S'", "CSV_FORMAT": "'YES'", "OBJ_DATA": "'NO'",
        "TLIST": f"'{mjd + 2400000.5}'", "TLIST_TYPE": "'JD'", "TIME_TYPE": "'UT'"}))["result"]
    row = res.split("$$SOE")[1].split("$$EOE")[0].strip().split(",")
    return [float(v) for v in row[2:8]]


def sexa(deg, hours=False):
    """Degrees -> 'hh-mm-ss.ss' / 'dd-mm-ss.ss' with JPL's 'M' for minus."""
    v = deg / 15 if hours else deg
    sign = "M" if v < 0 else ""
    v = abs(v)
    d, m = int(v), int(v * 60) % 60
    s = v * 3600 - d * 3600 - m * 60
    return f"{sign}{d:02d}-{m:02d}-{s:05.2f}"


def known_asteroids(img, vmag=19):
    """Asteroids JPL predicts in the image footprint at exposure mid-time."""
    mjd = (float(img["t_min"]) + float(img["t_max"])) / 2
    xobs = ",".join(f"{v:.6f}" for v in spherex_xyz(mjd))
    jd = mjd + 2400000.5
    # ponytail: axis-aligned 3.5 deg box around the center, not the rotated
    # footprint polygon; filter by WCS pixel coords once images are downloaded.
    res = json.loads(get("https://ssd-api.jpl.nasa.gov/sb_ident.api", {
        "xobs": xobs, "obs-time": f"{jd:.6f}", "sb-kind": "a", "two-pass": "true",
        "suppress-first-pass": "true", "vmag-lim": vmag,
        "fov-ra-center": sexa(float(img["s_ra"]), hours=True), "fov-dec-center": sexa(float(img["s_dec"])),
        "fov-ra-hwidth": 1.75, "fov-dec-hwidth": 1.75}))
    return res.get("fields_second", []), res.get("data_second_pass", [])


if __name__ == "__main__":
    assert sexa(-10.5) == "M10-30-00.00" and sexa(180, hours=True) == "12-00-00.00"
    ra = float(sys.argv[1]) if len(sys.argv) > 1 else 180.0
    dec = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    det = f"SPHEREx-D{sys.argv[3] if len(sys.argv) > 3 else 1}"
    imgs = sorted((i for i in images(ra, dec) if i["energy_bandpassname"] == det), key=lambda i: float(i["t_min"]))
    print(f"{len(imgs)} {det} images cover RA {ra} Dec {dec}:")
    for i in imgs:
        print(f"  MJD {float(i['t_min']):.4f}  {i['obs_id']}  {int(i['access_estsize']) // 1024} MB")
    if not imgs:
        sys.exit()
    fields, rows = known_asteroids(imgs[0])
    print(f"\n{len(rows)} known asteroids (V<19) predicted in {imgs[0]['obs_id']}:")
    for r in sorted(rows, key=lambda r: float(r[6]))[:15]:
        print("  ", dict(zip(fields, r)))

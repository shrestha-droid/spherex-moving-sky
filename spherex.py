#!/usr/bin/env python3
"""SPHEREx data access: image search, SPHEREx position, known asteroids, remote pixels.

usage: python3 spherex.py [RA_DEG] [DEC_DEG] [DETECTOR 1-6]   # list images + known asteroids
"""
import csv, io, json, os, re, sys, time, urllib.error, urllib.parse, urllib.request

SPHEREX = "-163182"  # Horizons spacecraft ID
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def get(url, params=None, headers=None, tries=5):
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers=headers or {})
    for n in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read(), r.headers
        except (urllib.error.URLError, ConnectionError, TimeoutError):
            if n == tries - 1:
                raise
            time.sleep(2 ** n)


def cached_json(name, fn):
    path = os.path.join(CACHE, name)
    if os.path.exists(path):
        return json.load(open(path))
    os.makedirs(CACHE, exist_ok=True)
    val = fn()
    json.dump(val, open(path, "w"))
    return val


def images(ra, dec, radius=0.01):
    """SPHEREx QR2 images covering a sky point, as dicts of IRSA SIA columns."""
    body, _ = get("https://irsa.ipac.caltech.edu/SIA", {
        "COLLECTION": "spherex_qr2", "POS": f"circle {ra} {dec} {radius}", "RESPONSEFORMAT": "CSV"})
    return list(csv.DictReader(io.StringIO(body.decode())))


def spherex_xyz(mjd):
    """Geocentric ICRF position (km) and velocity (km/s) of SPHEREx at an MJD (UTC)."""
    body, _ = get("https://ssd.jpl.nasa.gov/api/horizons.api", {
        "format": "json", "COMMAND": f"'{SPHEREX}'", "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@399'", "REF_PLANE": "'FRAME'", "VEC_TABLE": "'2'",
        "OUT_UNITS": "'KM-S'", "CSV_FORMAT": "'YES'", "OBJ_DATA": "'NO'",
        "TLIST": f"'{mjd + 2400000.5}'", "TLIST_TYPE": "'JD'", "TIME_TYPE": "'UT'"})
    res = json.loads(body)["result"]
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


def parse_sexa(s, hours=False):
    """'11:59:19.30' or '+01 06\\'25.0"' -> degrees."""
    parts = [float(p) for p in re.findall(r"[\d.]+", s)]
    v = parts[0] + parts[1] / 60 + parts[2] / 3600
    v = -v if s.strip().startswith("-") else v
    return v * 15 if hours else v


def known_asteroids(ra, dec, mjd, half_width=1.75, vmag=20):
    """Asteroids JPL predicts in a box around (ra, dec) at mjd, seen from SPHEREx.

    Returns [{name, ra, dec, vmag, ra_rate, dec_rate}] with rates in arcsec/hour."""
    def fetch():
        xobs = ",".join(f"{v:.6f}" for v in spherex_xyz(mjd))
        body, _ = get("https://ssd-api.jpl.nasa.gov/sb_ident.api", {
            "xobs": xobs, "obs-time": f"{mjd + 2400000.5:.6f}", "sb-kind": "a", "two-pass": "true",
            "suppress-first-pass": "true", "vmag-lim": vmag,
            "fov-ra-center": sexa(ra, hours=True), "fov-dec-center": sexa(dec),
            "fov-ra-hwidth": half_width, "fov-dec-hwidth": half_width})
        return json.loads(body)
    res = cached_json(f"sbident_{ra:.4f}_{dec:.4f}_{mjd:.5f}_{half_width}_{vmag}.json", fetch)
    out = []
    for r in res.get("data_second_pass", []):
        out.append({"name": r[0], "ra": parse_sexa(r[1], hours=True), "dec": parse_sexa(r[2]),
                    "vmag": float(r[6]), "ra_rate": float(r[7]), "dec_rate": float(r[8])})
    return out


class RemoteFits:
    """Read pixel rectangles from a remote FITS file with HTTP range requests.

    Downloads headers plus only the requested pixels, never the whole file."""

    def __init__(self, url, until="VARIANCE"):
        self.url, self.hdus = url, {}
        name = os.path.basename(urllib.parse.urlparse(url).path)
        heads = cached_json(f"hdr_{name}.json" if until == "VARIANCE" else f"hdr_{name}_{until}.json",
                            lambda: self._read_headers(until))
        from astropy.io import fits
        for ext, (text, start) in heads.items():
            self.hdus[ext] = (fits.Header.fromstring(text), start)

    def _read_headers(self, until):
        from astropy.io import fits
        out, off = {}, 0
        while True:
            text = b""
            while True:  # headers come in 2880-byte blocks ending with an END card
                body, _ = get(self.url, headers={"Range": f"bytes={off + len(text)}-{off + len(text) + 2880 * 4 - 1}"})
                text += body
                cards = [text[i:i + 80] for i in range(0, len(text), 80)]
                end = next((i for i, c in enumerate(cards) if c.startswith(b"END ")), None)
                if end is not None:
                    break
            hlen = -(-(end + 1) * 80 // 2880) * 2880
            hdr = fits.Header.fromstring(text[:hlen].decode("ascii"))
            ext = hdr.get("EXTNAME", "PRIMARY")
            out[ext] = (text[:hlen].decode("ascii"), off + hlen)
            nbytes = abs(hdr["BITPIX"]) // 8
            for i in range(1, hdr["NAXIS"] + 1):
                nbytes *= hdr[f"NAXIS{i}"]
            if hdr["NAXIS"] == 0:
                nbytes = 0
            off += hlen + -(-nbytes // 2880) * 2880
            if ext == until:
                return out

    def wavelength(self, x, y):
        """(wavelength, bandwidth) in microns at 0-based pixel (x, y), from the WCS-WAVE lookup table.

        Needs RemoteFits(url, until="WCS-WAVE")."""
        import numpy as np
        if not hasattr(self, "_wave"):
            hdr, start = self.hdus["WCS-WAVE"]
            path = os.path.join(CACHE, "px", os.path.basename(self.url) + "_wave.npy")
            if os.path.exists(path):
                raw = np.load(path).tobytes()
            else:
                raw, _ = get(self.url, headers={"Range": f"bytes={start}-{start + hdr['NAXIS1'] - 1}"})
                os.makedirs(os.path.dirname(path), exist_ok=True)
                np.save(path, np.frombuffer(raw, np.uint8))
            n = lambda k: int(re.match(r"(\d+)", hdr[k]).group(1))  # TFORMn = '13J' etc.
            nx, ny = n("TFORM1"), n("TFORM2")
            gx = np.frombuffer(raw[:4 * nx], ">i4").astype(float) - 1  # table is 1-based
            gy = np.frombuffer(raw[4 * nx:4 * (nx + ny)], ">i4").astype(float) - 1
            vals = np.frombuffer(raw[4 * (nx + ny):4 * (nx + ny) + 8 * nx * ny], ">f4").reshape(ny, nx, 2)
            self._wave = gx, gy, vals  # FITS TDIM (2,nx,ny) -> [iy][ix][wave|bw]
        gx, gy, vals = self._wave
        fx, fy = np.interp(x, gx, np.arange(len(gx))), np.interp(y, gy, np.arange(len(gy)))
        ix, iy = min(int(fx), len(gx) - 2), min(int(fy), len(gy) - 2)
        tx, ty = fx - ix, fy - iy
        corners = [(iy, ix, (1 - tx) * (1 - ty)), (iy, ix + 1, tx * (1 - ty)), (iy + 1, ix, (1 - tx) * ty), (iy + 1, ix + 1, tx * ty)]
        good = [(vals[a, b], w) for a, b, w in corners if vals[a, b, 0] > 0]  # zeros = no filter coverage
        if not good:
            return float("nan"), float("nan")
        wsum = sum(w for _, w in good) or 1.0
        v = sum(val * w for val, w in good) / wsum if wsum > 0 else good[0][0]
        return float(v[0]), float(v[1])

    def header(self, ext="IMAGE"):
        return self.hdus[ext][0]

    def read(self, ext, x0, y0, w, h):
        """Pixels [y0:y0+h, x0:x0+w] (0-based) of an image HDU; off-image pixels are NaN."""
        import numpy as np
        path = os.path.join(CACHE, "px", f"{os.path.basename(self.url)}_{ext}_{x0}_{y0}_{w}_{h}.npy")
        if os.path.exists(path):
            return np.load(path)
        out = self._read(ext, x0, y0, w, h)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        np.save(path, out.astype(np.float32) if out.dtype.kind == "f" else out)
        return out

    def _read(self, ext, x0, y0, w, h):
        import numpy as np
        hdr, start = self.hdus[ext]
        nx, ny, bpp = hdr["NAXIS1"], hdr["NAXIS2"], abs(hdr["BITPIX"]) // 8
        dtype = {-32: ">f4", 32: ">i4", -64: ">f8", 16: ">i2"}[hdr["BITPIX"]]
        out = np.full((h, w), np.nan if hdr["BITPIX"] < 0 else 0, dtype=float if hdr["BITPIX"] < 0 else np.int64)
        xa, xb = max(x0, 0), min(x0 + w, nx)
        rows = [y for y in range(max(y0, 0), min(y0 + h, ny))]
        if xa >= xb or not rows:
            return out
        spans = [(start + (y * nx + xa) * bpp, start + (y * nx + xb) * bpp - 1) for y in rows]
        chunks = {}
        for i in range(0, len(spans), 100):  # servers cap ranges per request
            chunks.update(self._ranges(spans[i:i + 100]))
        for y, (a, _) in zip(rows, spans):
            out[y - y0, xa - x0:xb - x0] = np.frombuffer(chunks[a], dtype=dtype)
        return out

    def _ranges(self, spans):
        body, headers = get(self.url, headers={"Range": "bytes=" + ",".join(f"{a}-{b}" for a, b in spans)})
        ctype = headers.get("Content-Type", "")
        if "multipart/byteranges" not in ctype:
            return {spans[0][0]: body}
        bnd = b"--" + ctype.split("boundary=")[1].strip().encode()
        out, pos = {}, body.find(bnd)
        while True:
            hend = body.find(b"\r\n\r\n", pos)
            if hend < 0:
                return out
            m = re.search(rb"bytes (\d+)-(\d+)/", body[pos:hend], re.I)
            if not m:
                return out
            a, b = int(m.group(1)), int(m.group(2))
            out[a] = body[hend + 4:hend + 5 + b - a]
            pos = body.find(bnd, hend + 5 + b - a)


if __name__ == "__main__":
    assert sexa(-10.5) == "M10-30-00.00" and sexa(180, hours=True) == "12-00-00.00"
    assert abs(parse_sexa("11:59:19.30", hours=True) - 179.830417) < 1e-5
    assert abs(parse_sexa("-01 06'25.0\"") + 1.106944) < 1e-5
    ra = float(sys.argv[1]) if len(sys.argv) > 1 else 180.0
    dec = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    det = f"SPHEREx-D{sys.argv[3] if len(sys.argv) > 3 else 1}"
    imgs = sorted((i for i in images(ra, dec) if i["energy_bandpassname"] == det), key=lambda i: float(i["t_min"]))
    print(f"{len(imgs)} {det} images cover RA {ra} Dec {dec}:")
    for i in imgs:
        print(f"  MJD {float(i['t_min']):.4f}  {i['obs_id']}  {int(i['access_estsize']) // 1024} MB")
    if imgs:
        i = imgs[0]
        rocks = known_asteroids(float(i["s_ra"]), float(i["s_dec"]), (float(i["t_min"]) + float(i["t_max"])) / 2, vmag=19)
        print(f"\n{len(rocks)} known asteroids (V<19) predicted in {i['obs_id']}:")
        for r in sorted(rocks, key=lambda r: r["vmag"])[:15]:
            print(f"  V={r['vmag']:<5} {r['name']}")

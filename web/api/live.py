"""Live SPHEREx lookups for the web page (Vercel Python function, standard library only).

GET /api/live?op=images&ra=..&dec=..          SPHEREx exposures covering a sky position (IRSA SIA)
GET /api/live?op=cutout&url=..&ra=..&dec=..   pixels around a position, read live from the 69 MB
                                              file with HTTP range requests; returns the byte log
GET /api/live?op=known&url=..&ra=..&dec=..    known asteroids JPL predicts in that cutout, as seen
                                              from SPHEREx at the exposure time
GET /api/live?op=where&name=..&mjd=..         where JPL says a body is, seen from SPHEREx

Only NASA/IPAC IRSA and JPL hosts are ever contacted, and only SPHEREx file paths.
"""
import json, math, re, time, urllib.error, urllib.parse, urllib.request
from array import array
from http.server import BaseHTTPRequestHandler

IRSA_FILES = "https://irsa.ipac.caltech.edu/ibe/data/spherex/"
UA = {"User-Agent": "spherex-moving-sky/1.0 (NASA Space Apps 2026)"}
MAX_HALF = 31


class Log:
    def __init__(self):
        self.calls, self.t0 = [], time.time()

    def get(self, url, params=None, rng=None, timeout=25):
        if params:
            url += "?" + urllib.parse.urlencode(params)
        h = dict(UA)
        if rng:
            h["Range"] = rng
        t = time.time()
        for attempt in range(3):  # NASA servers occasionally drop a connection under load
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=timeout) as r:
                    body, status, ctype = r.read(), r.status, r.headers.get("Content-Type", "")
                break
            except (urllib.error.URLError, ConnectionError, TimeoutError):
                if attempt == 2:
                    raise
                time.sleep(0.6 * (attempt + 1))
        host = urllib.parse.urlparse(url).netloc
        self.calls.append({"host": host, "path": urllib.parse.urlparse(url).path[-60:], "range": rng,
                           "status": status, "bytes": len(body), "ms": round((time.time() - t) * 1000)})
        return body, ctype


def cards(text):
    """Minimal FITS header parser -> dict (numbers, strings, booleans)."""
    out = {}
    for i in range(0, len(text), 80):
        c = text[i:i + 80]
        if c.startswith("END "):
            break
        if c[8:10] != "= ":
            continue
        k, v = c[:8].strip(), c[10:]
        if v.lstrip().startswith("'"):
            m = re.match(r"\s*'((?:[^']|'')*)'", v)
            out[k] = m.group(1).replace("''", "'").rstrip() if m else ""
            continue
        v = v.split("/")[0].strip()
        if v in ("T", "F"):
            out[k] = v == "T"
        else:
            try:
                out[k] = int(v)
            except ValueError:
                try:
                    out[k] = float(v)
                except ValueError:
                    out[k] = v
    return out


def headers(log, url, want=("IMAGE",)):
    """Walk HDU headers with range requests until all `want` are found -> {name: (cards, data_offset)}."""
    out, off = {}, 0
    while len(out) < len(want):
        text = b""
        while b"END" + b" " * 77 not in text:
            body, _ = log.get(url, rng=f"bytes={off + len(text)}-{off + len(text) + 2880 * 4 - 1}")
            text += body
            if len(text) > 2880 * 40:
                raise ValueError("header too long")
        end = text.index(b"END" + b" " * 77)
        while end % 80:
            end = text.index(b"END" + b" " * 77, end + 1)
        hlen = -(-(end + 80) // 2880) * 2880
        h = cards(text[:hlen].decode("ascii", "replace"))
        name = h.get("EXTNAME", "PRIMARY")
        if name in want:
            out[name] = (h, off + hlen)
        n = abs(h.get("BITPIX", 8)) // 8
        for i in range(1, h.get("NAXIS", 0) + 1):
            n *= h[f"NAXIS{i}"]
        off += hlen + (-(-n // 2880) * 2880 if h.get("NAXIS", 0) else 0)
    return out


def world_to_pixel(h, ra, dec):
    """TAN-SIP sky -> 0-based pixel, using the header's inverse distortion (AP/BP)."""
    a0, d0, a, d = (math.radians(x) for x in (h["CRVAL1"], h["CRVAL2"], ra, dec))
    cosc = math.sin(d0) * math.sin(d) + math.cos(d0) * math.cos(d) * math.cos(a - a0)
    xi = math.degrees(math.cos(d) * math.sin(a - a0) / cosc)
    eta = math.degrees((math.cos(d0) * math.sin(d) - math.sin(d0) * math.cos(d) * math.cos(a - a0)) / cosc)
    c11, c12 = h.get("CDELT1", 1) * h.get("PC1_1", 1), h.get("CDELT1", 1) * h.get("PC1_2", 0)
    c21, c22 = h.get("CDELT2", 1) * h.get("PC2_1", 0), h.get("CDELT2", 1) * h.get("PC2_2", 1)
    det = c11 * c22 - c12 * c21
    U, V = (c22 * xi - c12 * eta) / det, (-c21 * xi + c11 * eta) / det
    du = sum(h.get(f"AP_{p}_{q}", 0) * U ** p * V ** q for p in range(h.get("AP_ORDER", 0) + 1) for q in range(h.get("AP_ORDER", 0) + 1 - p))
    dv = sum(h.get(f"BP_{p}_{q}", 0) * U ** p * V ** q for p in range(h.get("BP_ORDER", 0) + 1) for q in range(h.get("BP_ORDER", 0) + 1 - p))
    return h["CRPIX1"] + U + du - 1, h["CRPIX2"] + V + dv - 1


def read_rect(log, url, h, start, x0, y0, w, hh):
    """float32 pixels [y0:y0+hh, x0:x0+w] in one multi-range request."""
    nx = h["NAXIS1"]
    spans = [(start + (y * nx + x0) * 4, start + (y * nx + x0 + w) * 4 - 1) for y in range(y0, y0 + hh)]
    body, ctype = log.get(url, rng="bytes=" + ",".join(f"{a}-{b}" for a, b in spans))
    parts = {}
    if "multipart/byteranges" in ctype:
        bnd = b"--" + ctype.split("boundary=")[1].strip().encode()
        pos = body.find(bnd)
        while pos >= 0:
            hend = body.find(b"\r\n\r\n", pos)
            m = re.search(rb"bytes (\d+)-(\d+)/", body[pos:hend], re.I) if hend > 0 else None
            if not m:
                break
            a, b = int(m.group(1)), int(m.group(2))
            parts[a] = body[hend + 4:hend + 5 + b - a]
            pos = body.find(bnd, hend + 5 + b - a)
    else:
        parts[spans[0][0]] = body
    rows = []
    for a, _ in spans:
        arr = array("f", parts[a])
        arr.byteswap()  # FITS is big-endian
        rows.append(list(arr))
    return rows


def to_u8(rows):
    vals = sorted(v for r in rows for v in r if v == v)
    if not vals:
        return [0] * sum(len(r) for r in rows)
    lo, hi = vals[len(vals) // 2], vals[int(len(vals) * 0.998)]
    span = max(hi - lo, 1e-9)
    return [0 if v != v else int(255 * min(1, math.asinh(max(0, (v - lo) / span) * 10) / math.asinh(10))) for r in rows for v in r]


def sexa(deg, hours=False):
    v = deg / 15 if hours else deg
    sign, v = ("M" if v < 0 else ""), abs(v)
    d, m = int(v), int(v * 60) % 60
    return f"{sign}{d:02d}-{m:02d}-{v * 3600 - d * 3600 - m * 60:05.2f}"


def spherex_xyz(log, mjd):
    body, _ = log.get("https://ssd.jpl.nasa.gov/api/horizons.api", {
        "format": "json", "COMMAND": "'-163182'", "EPHEM_TYPE": "'VECTORS'", "CENTER": "'500@399'",
        "REF_PLANE": "'FRAME'", "VEC_TABLE": "'2'", "OUT_UNITS": "'KM-S'", "CSV_FORMAT": "'YES'",
        "OBJ_DATA": "'NO'", "TLIST": f"'{mjd + 2400000.5}'", "TLIST_TYPE": "'JD'", "TIME_TYPE": "'UT'"})
    row = json.loads(body)["result"].split("$$SOE")[1].split("$$EOE")[0].strip().split(",")
    return ",".join(f"{float(v):.6f}" for v in row[2:8])


def parse_sexa(s, hours=False):
    p = [float(x) for x in re.findall(r"[\d.]+", s)]
    v = p[0] + p[1] / 60 + p[2] / 3600
    v = -v if s.strip().startswith("-") else v
    return v * 15 if hours else v


def num(q, k, lo, hi):
    v = float(q[k][0])
    if not (lo <= v <= hi) or v != v:
        raise ValueError(f"{k} out of range")
    return v


def op_images(q, log):
    ra, dec = num(q, "ra", 0, 360), num(q, "dec", -90, 90)
    body, _ = log.get("https://irsa.ipac.caltech.edu/SIA", {
        "COLLECTION": "spherex_qr2", "POS": f"circle {ra} {dec} 0.002", "RESPONSEFORMAT": "CSV", "MAXREC": "400"})
    import csv, io
    rows = list(csv.DictReader(io.StringIO(body.decode())))
    rows.sort(key=lambda r: float(r["t_min"]))
    return {"ra": ra, "dec": dec, "count": len(rows), "images": [
        {"obs": r["obs_id"], "det": r["energy_bandpassname"][-2:], "mjd": round((float(r["t_min"]) + float(r["t_max"])) / 2, 5),
         "um": [round(float(r["em_min"]) * 1e6, 2), round(float(r["em_max"]) * 1e6, 2)], "url": r["access_url"]}
        for r in rows[:120]]}


def op_cutout(q, log):
    url = q["url"][0]
    if not url.startswith(IRSA_FILES) or not url.endswith(".fits") or ".." in url:
        raise ValueError("only SPHEREx files on IRSA are allowed")
    ra, dec = num(q, "ra", 0, 360), num(q, "dec", -90, 90)
    half = int(num(q, "half", 4, MAX_HALF)) if "half" in q else 20
    hd = headers(log, url)
    h, start = hd["IMAGE"]
    x, y = world_to_pixel(h, ra, dec)
    X, Y = int(round(x)), int(round(y))
    x0, y0 = max(0, X - half), max(0, Y - half)
    x1, y1 = min(h["NAXIS1"], X + half + 1), min(h["NAXIS2"], Y + half + 1)
    if x1 <= x0 or y1 <= y0:
        return {"inside": False, "x": x, "y": y, "log": log.calls}
    rows = read_rect(log, url, h, start, x0, y0, x1 - x0, y1 - y0)
    return {"inside": True, "x": round(x, 2), "y": round(y, 2), "x0": x0, "y0": y0, "w": x1 - x0, "h": y1 - y0,
            "stamp": to_u8(rows), "obs": h.get("OBSID"), "date": h.get("DATE-AVG"), "detector": h.get("DETECTOR"),
            "mjd": h.get("MJD-AVG"), "file_bytes": None, "log": log.calls,
            "total_ms": round((time.time() - log.t0) * 1000), "bytes": sum(c["bytes"] for c in log.calls)}


def op_known(q, log):
    url = q["url"][0]
    if not url.startswith(IRSA_FILES) or not url.endswith(".fits"):
        raise ValueError("only SPHEREx files on IRSA are allowed")
    ra, dec = num(q, "ra", 0, 360), num(q, "dec", -90, 90)
    hd = headers(log, url)
    h, _ = hd["IMAGE"]
    mjd = float(h["MJD-AVG"])
    body, _ = log.get("https://ssd-api.jpl.nasa.gov/sb_ident.api", {
        "xobs": spherex_xyz(log, mjd), "obs-time": f"{mjd + 2400000.5:.6f}", "sb-kind": "a",
        "fov-ra-center": sexa(ra, hours=True), "fov-dec-center": sexa(dec), "fov-ra-hwidth": 0.045,
        "fov-dec-hwidth": 0.045, "vmag-lim": 20})  # just the cutout (49 px x 6.15")
    res = json.loads(body)
    out = []
    for r in res.get("data_first_pass", []):
        a, d = parse_sexa(r[1], hours=True), parse_sexa(r[2])
        x, y = world_to_pixel(h, a, d)
        out.append({"name": r[0], "ra": round(a, 5), "dec": round(d, 5), "V": r[6], "x": round(x, 2), "y": round(y, 2)})
    return {"mjd": mjd, "objects": out, "log": log.calls}


def op_where(q, log):
    name = q["name"][0].strip()
    if not re.fullmatch(r"[\w .'()-]{1,40}", name):
        raise ValueError("bad name")
    mjd = num(q, "mjd", 60700, 62500)
    body, _ = log.get("https://ssd.jpl.nasa.gov/api/horizons.api", {
        "format": "json", "COMMAND": f"'{name};'", "MAKE_EPHEM": "'YES'", "EPHEM_TYPE": "'OBSERVER'",
        "CENTER": "'500@-163182'", "TLIST": f"'{mjd + 2400000.5}'", "TLIST_TYPE": "'JD'",
        "QUANTITIES": "'1,9'", "ANG_FORMAT": "'DEG'", "CSV_FORMAT": "'YES'", "OBJ_DATA": "'NO'"})
    res = json.loads(body).get("result", "")
    if "$$SOE" not in res:
        raise ValueError("JPL could not resolve that name uniquely; try its number, e.g. 872")
    p = [x.strip() for x in res.split("$$SOE")[1].split("$$EOE")[0].strip().split(",")]
    tgt = re.search(r"Target body name: (.+?)\s{2,}", res)
    return {"name": tgt.group(1) if tgt else name, "ra": float(p[3]), "dec": float(p[4]),
            "V": p[5], "mjd": mjd, "log": log.calls}


OPS = {"images": op_images, "cutout": op_cutout, "known": op_known, "where": op_where}


class handler(BaseHTTPRequestHandler):
    def do_GET(self):
        q = urllib.parse.parse_qs(urllib.parse.urlparse(self.path).query)
        log = Log()
        try:
            op = OPS[q.get("op", [""])[0]]
            code, out = 200, op(q, log)
        except (KeyError, ValueError) as e:
            code, out = 400, {"error": str(e) or "bad request", "log": log.calls}
        except Exception as e:  # upstream NASA/JPL failure
            code, out = 502, {"error": f"upstream: {e}", "log": log.calls}
        body = json.dumps(out).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Cache-Control", "public, max-age=3600" if code == 200 else "no-store")
        self.end_headers()
        self.wfile.write(body)

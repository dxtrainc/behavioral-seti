"""
RSTN 1-second solar radio flux -- acquisition.

Row 2 of the designer-priority table: radio spectral index at second cadence,
a channel with the resolving power a designer would want and no search run
against it.

THE CARRIER IS THE SPECTRAL INDEX, not the flux. RSTN reports 8 frequencies
(245, 410, 610, 1415, 2695, 4995, 8800, 15400 MHz) from each station, so the
ratio of any two is dimensionless by construction and satisfies section 3.1
without needing a unit shared with anybody. Seven are populated at Learmonth,
giving 21 pairs.

TWO STATIONS, DELIBERATELY. Learmonth and San Vito are different hardware on
opposite sides of the Earth. That is a twin gate in the sense of section 4.4: a
solar signal appears at both, while interference, a receiver calibration cycle
or a local artefact appears at one. The Oulu/Kiel gate worked for the same
reason and the GOES-16/17 gate failed because those two share a design.

Files are daily, gzipped, about 200 kB, covering the roughly ten hours the Sun
is above the horizon at each site.
"""
import os, sys, gzip, time, datetime as dt
import numpy as np, requests

BASE = ("https://www.ngdc.noaa.gov/stp/space-weather/solar-data/"
        "solar-features/solar-radio/rstn-1-second")
UA = {"User-Agent": "beacon-search/1.0 (rtg@dxtra.com; academic solar study)"}
STATIONS = {"learmonth": "apl", "san-vito": "lis", "palehua": "phf"}
MON = ["jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"]
OUT = os.path.expanduser("~/rstn")
os.makedirs(OUT, exist_ok=True)
S = requests.Session(); S.headers.update(UA)

def fetch_day(station, suf, d):
    name = "%02d%s%02d.%s.gz" % (d.day, MON[d.month-1], d.year % 100, suf)
    url = "%s/%s/%04d/%02d/%s" % (BASE, station, d.year, d.month, name)
    dst = os.path.join(OUT, "%s_%s.gz" % (station, d.strftime("%Y%m%d")))
    if os.path.exists(dst): return True
    try:
        r = S.get(url, timeout=120)
        if r.status_code != 200 or len(r.content) < 2000: return False
        open(dst, "wb").write(r.content)
        return True
    except Exception:
        return False

def main():
    y = int(sys.argv[1]) if len(sys.argv) > 1 else 2014
    m0 = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    nmon = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    got = {k: 0 for k in STATIONS}
    for station, suf in STATIONS.items():
        d = dt.date(y, m0, 1)
        end = dt.date(y + (m0+nmon-1)//13, ((m0+nmon-1) % 12) + 1, 1)
        while d < end:
            if fetch_day(station, suf, d): got[station] += 1
            d += dt.timedelta(days=1)
            time.sleep(0.15)
        print("%-12s %d days" % (station, got[station]), flush=True)
    print("total files: %d" % sum(got.values()), flush=True)

main()

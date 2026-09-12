"""
ROW 5 -- passive occulter, information in transit TIMING. SOHO/VIRGO TSI, 1 min.

§3.6 ranks this row fifth by designer preference at a level of 10 ppm dips, and
records it unsearched: "§A.2 row 27 ran daily, single-dip". This is the minute-cadence
version, and the sequence version, after Arnold (2005).

WHY THIS ARCHIVE AND NOT SPM. The SPM photometry of §4.15 carries a two-month highpass
applied upstream by the VDC. **VIRGO TSI minute does not** -- it is absolute
irradiance in W/m^2 on the WRR scale, 14,199,837 samples at 60 s from 1996 to 2023,
92.8% present. A dip search needs the low frequencies a highpass removes, so this is
the right product and SPM is the wrong one.

GATE 1, AND IT IS THE BEST CONTROL IN THIS PROGRAMME. The Sun is transited by known
planets at known times, and the depth is not fitted -- it is the ratio of two disk
areas, computable from first principles:

    Venus    (6052/696000)^2  =  75.6 ppm
    Mercury  (2440/696000)^2  =  12.3 ppm

Venus transited on 2004-06-08 and 2012-06-06; Mercury on 2003-05-07, 2006-11-08,
2016-05-09 and 2019-11-11. Both depths sit at or above the 10 ppm a designer would
set for this row. **A dip search that cannot recover Venus 2004 is not a search**, and
unlike most controls here this one has a known date, a known duration and a known
depth, none of them derived from the data.

It also calibrates the DEPTH SCALE directly: if recovered Venus comes out at 76 ppm
the photometry and the detrending are both behaving, and the limit that follows is in
real units rather than in units of the pipeline's own noise.
"""
import os, sys, json, datetime as dt
import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter

VD = os.path.expanduser("~/virgo")
TAI0 = dt.datetime(1958, 1, 1)          # TAI epoch
RSUN = 696000.0
TRANSITS = [("Venus",   "2004-06-08", (6052.0/RSUN)**2),
            ("Venus",   "2012-06-06", (6052.0/RSUN)**2),
            ("Mercury", "2003-05-07", (2440.0/RSUN)**2),
            ("Mercury", "2006-11-08", (2440.0/RSUN)**2),
            ("Mercury", "2016-05-09", (2440.0/RSUN)**2),
            ("Mercury", "2019-11-11", (2440.0/RSUN)**2)]

def load():
    d = np.asarray(fits.open(os.path.join(VD, "VIRGO_TSI_Minute.fits"))[0].data, float)
    t, v = d[:, 0], d[:, 1]              # TAI seconds, W/m^2 (WRR scale)
    return t, v

def to_dt(tai):
    return TAI0 + dt.timedelta(seconds=float(tai))

def detrend(v, win=1441):
    """relative residual against a running median ~1 day: solar variability out,
    a few-hour dip left in. Gaps are left as NaN and never interpolated."""
    ok = np.isfinite(v)
    f = np.where(ok, v, np.nan)
    med = median_filter(np.nan_to_num(f, nan=np.nanmedian(f)), size=win, mode="nearest")
    return np.where(ok, f/med - 1.0, np.nan), med

if __name__ == "__main__":
    t, v = load()
    print("VIRGO TSI minute: %d samples, %.1f%% present" % (len(v), 100*np.isfinite(v).mean()))
    print("  %s .. %s" % (to_dt(t[0]).date(), to_dt(t[-1]).date()))
    print("  median %.2f W/m^2\n" % np.nanmedian(v))

    r, med = detrend(v)
    print("  residual rms %.2f ppm per minute" % (np.nanstd(r)*1e6))

    print("\nGATE 1 -- known planetary transits against a CONTROL-DAY null")
    print("  The first version of this gate compared the deepest hour on the transit")
    print("  day against 0.4x the predicted depth. Residual rms is 62.4 ppm/min, so")
    print("  NOISE cleared that threshold and all six 'passed' -- Mercury reading")
    print("  50-70 ppm against a predicted 12.3. A control that cannot fail is not a")
    print("  control. The depth must now stand out against the same statistic measured")
    print("  on transit-free days in the same season.\n")

    def deepest_hour(day_start):
        m = (t >= day_start) & (t < day_start + 86400)
        seg = r[m]; seg = seg[np.isfinite(seg)]
        if len(seg) < 600: return np.nan
        k = 60
        return -np.convolve(seg, np.ones(k)/k, "valid").min()*1e6

    print("  %-9s %-12s %9s %9s %9s %8s %7s" % ("body","date","expected","measured","ctrl med","ctrl 99%","gate"))
    got = []
    for body, date, depth in TRANSITS:
        d0 = dt.datetime.strptime(date, "%Y-%m-%d")
        s0 = (d0 - TAI0).total_seconds()
        meas = deepest_hour(s0)
        # control days: +/-10 to +/-60 days, same season, no transit
        ctrl = []
        for off in list(range(-60, -9)) + list(range(10, 61)):
            c = deepest_hour(s0 + off*86400)
            if np.isfinite(c): ctrl.append(c)
        if not np.isfinite(meas) or len(ctrl) < 20:
            print("  %-9s %-12s %9.1f %9s %9s %8s %7s" % (body,date,depth*1e6,"no data","","","--"))
            continue
        ctrl = np.array(ctrl)
        thr = float(np.quantile(ctrl, 0.99))
        ok = meas > thr
        got.append((body, date, depth*1e6, float(meas), float(np.median(ctrl)), thr, bool(ok)))
        print("  %-9s %-12s %9.1f %9.1f %9.1f %8.1f %7s"
              % (body, date, depth*1e6, meas, np.median(ctrl), thr, "PASS" if ok else "FAIL"))

    np.save(os.path.expanduser("~/virgo_tsi_resid.npy"), r.astype(np.float32))
    json.dump([{"body": b, "date": d, "expected_ppm": e, "measured_ppm": m,
                "ctrl_median_ppm": cm, "ctrl_99_ppm": c99, "pass": o}
               for b, d, e, m, cm, c99, o in got],
              open(os.path.expanduser("~/virgo_row5_gate1.json"), "w"), indent=1)
    npass = sum(1 for g in got if g[-1])
    print("\n  %d of %d known transits recovered" % (npass, len(got)))
    print("  saved residual series for the dip search")

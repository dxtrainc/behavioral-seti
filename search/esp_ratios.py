"""
ROW 3, THE DIMENSIONLESS CARRIERS -- EUV colour ratios on SDO/EVE ESP.

search/esp_row3.py searches the four ESP light channels individually. But §3.1 is
explicit that a sender sharing no units with us can encode only in DIMENSIONLESS
quantities, so single channels are the wrong object and ratios are the right one.
This searches the six colour ratios among CH_18, CH_26, CH_30 and CH_36.

THE RATIOS COME FREE, AS THEY DID FOR VIRGO SPM. Each channel is prepared as a
relative residual r = x/trend - 1, so to first order log(x_i/x_j) = r_i - r_j and the
DIFFERENCE of two prepared channels *is* their log colour ratio. No constructed
quantity, no new normalisation.

THE THREE-CLAUSE REGULARISER APPLIES AND BITES HERE. All four bands come from ONE
instrument behind one optical path, so the ratios are not independent of each other in
the way the §3.1 rule demands, and any survivor must be read with that in mind. This is
the same trap that produced sixteen spurious survivors in the first pair sweep when
plasma beta and the Alfven Mach number were included: dimensionless is necessary and
not sufficient. It is stated here rather than discovered later.

GATE: CH_D, the dark diode. Anything appearing there is instrumental by construction.
A colour ratio cannot be formed with a dark channel that carries no signal, so the
dark channel is used as a VETO on frequency rather than as a seventh ratio.
"""
import os, sys, json, time, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import esp_row3 as E

if __name__ == "__main__":
    t0 = time.time()
    C, ok, ndays = E.build()
    print("EVE ESP colour ratios: %d days, %d samples at %.2f s, %.1f%% present"
          % (ndays, len(ok), E.DT, 100*ok.mean()), flush=True)

    R = {}
    for c in E.LIGHT + [E.DARK]:
        R[c] = E.prep(C[c], ok)
        # reported for every channel so the exclusion below is visible, not silent
        print("  %-6s residual rms %.1f ppm" % (c, np.std(R[c][ok])*1e6), flush=True)

    # dark channel spectrum, for the veto
    fd, Pd = E.spectrum(R[E.DARK])
    Rd = Pd/E.continuum(Pd)
    Nd = len(Rd); thrd = (np.median(Rd)/np.log(2.0))*np.log(Nd/0.05)
    print("\n  dark-channel veto threshold R > %.1f" % thrd, flush=True)

    print("\nBLIND NARROWBAND SEARCH -- six dimensionless colour ratios")
    allc = {}
    # CH_36 IS EXCLUDED. Its values cross zero (median 2.6e-4, minimum -5.2e-4), so the
    # relative residual x/trend - 1 diverges and its residual rms is 765,522 ppm -- 76%.
    # A ratio formed against a divergent channel is not a dimensionless carrier, it is a
    # division by something near zero, and any survivor in it would be that and nothing
    # else. Section 4.18 already excludes CH_36 from any limit; the ratios inherit the
    # exclusion. This removes three of the six pairs and leaves 18/26, 18/30 and 26/30.
    USABLE = [c for c in E.LIGHT if c != "CH_36"]
    print("\n  excluded CH_36 (residual rms 76%%, values cross zero); %d ratios from %s"
          % (len(list(itertools.combinations(USABLE, 2))), ", ".join(USABLE)), flush=True)
    for a, b in itertools.combinations(USABLE, 2):
        d = R[a] - R[b]                       # = log colour ratio, to first order
        f, P = E.spectrum(d)
        Rr = P/E.continuum(P)
        N = len(Rr); mu = np.median(Rr)/np.log(2.0); thr = mu*np.log(N/0.05)
        idx = np.where(Rr > thr)[0]
        lab = "%s/%s" % (a.replace("CH_", ""), b.replace("CH_", ""))
        print("  %-9s %d bins, threshold R > %.1f, %d above" % (lab, N, thr, len(idx)), flush=True)
        rows = []
        for i in idx[np.argsort(-Rr[idx])][:8]:
            fr = float(f[i])
            j = int(np.argmin(np.abs(fd-fr)))
            dark = float(Rd[max(0, j-2):j+3].max())
            vetoed = dark >= thrd
            rows.append((fr, float(Rr[i]), dark, vetoed))
            per = 1.0/fr
            tag = ("%.3f s" % per) if per < 3600 else (("%.2f h" % (per/3600)) if per < 86400 else ("%.2f d" % (per/86400)))
            print("      f = %.6e Hz  period %12s  R = %8.1f  dark %7.1f  %s"
                  % (fr, tag, Rr[i], dark, "VETOED (instrumental)" if vetoed else "clean"))
        if not len(idx): print("      no bins above threshold")
        allc[lab] = rows

    survivors = [(k, r) for k, rs in allc.items() for r in rs if not r[3]]
    print("\n  candidates surviving the dark-channel veto: %d" % len(survivors))
    for k, r in survivors[:12]:
        per = 1.0/r[0]
        tag = ("%.3f s" % per) if per < 3600 else (("%.2f h" % (per/3600)) if per < 86400 else ("%.2f d" % (per/86400)))
        print("      %-9s %.6e Hz  %12s  R = %.1f" % (k, r[0], tag, r[1]))
    json.dump({k: [{"f": a, "R": b, "dark": c, "vetoed": d} for a, b, c, d in v]
               for k, v in allc.items()},
              open(os.path.expanduser("~/esp_ratios.json"), "w"), indent=1)
    print("\n  %.0f s" % (time.time()-t0))

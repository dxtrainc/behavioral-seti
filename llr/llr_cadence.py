import numpy as np, glob, os
import llr_model_free as M
S = M.sessions()
dts = []
for t, v in S:
    t = np.sort(M.unwrap_day(t))
    if len(t) > 5: dts.append(np.diff(t))
d = np.concatenate(dts)
d = d[(d > 0) & (d < 600)]
print("inter-return interval: median %.3f s, mode-ish %.3f s, n=%d"
      % (np.median(d), np.median(d[d < np.percentile(d, 75)]), len(d)))
h, e = np.histogram(d, bins=np.arange(0, 60, 0.2))
k = int(np.argmax(h))
print("commonest interval: %.1f-%.1f s (%d of %d)" % (e[k], e[k+1], h[k], len(d)))
print("strongest LS peak was 18.76 s -> ratio to median interval: %.3f"
      % (18.76/np.median(d)))

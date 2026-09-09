import numpy as np, llr_model_free as M
S = M.sessions()
allt = np.concatenate([t for t, _ in S])
fs = np.linspace(1.0, 2000.0, 120000)
z = np.empty(len(fs))
for i in range(0, len(fs), 1500):
    ph = np.exp(2j*np.pi*np.outer(fs[i:i+1500], allt))
    z[i:i+1500] = np.abs(ph.mean(axis=1))
chance = 1.0/np.sqrt(len(allt)); zm = z.max()
cand = fs[z > 0.85*zm]
print("epochs %d   chance concentration %.4f   max %.4f" % (len(allt), chance, zm))
print("strongest       %9.3f Hz = %.6f s" % (fs[int(np.argmax(z))], 1.0/fs[int(np.argmax(z))]))
print("lowest in-band  %9.3f Hz = %.6f s   <- the grid the hardware fires on"
      % (cand.min(), 1.0/cand.min()))
pooled = []
for t, v in S:
    t = M.unwrap_day(t)
    if np.ptp(t) < 60: continue
    r = M.detrend(t, v); r = r[np.isfinite(r)]
    if len(r): pooled.append(r)
p = np.concatenate(pooled)
mad = 1.4826*np.median(np.abs(p-np.median(p)))
print("pooled residual  %8.1f ps = %6.1f mm one-way   n=%d" % (np.std(p)*1e12, np.std(p)*2.998e8/2*1e3, len(p)))
print("robust (MAD)     %8.1f ps = %6.1f mm" % (mad*1e12, mad*2.998e8/2*1e3))

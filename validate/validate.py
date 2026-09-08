"""Validate the GPD tail fit against direct exceedance counting.
argv: <counted.json> <tailfit.json>.  Exit 0 = PASS, 1 = FAIL."""
import json, os, sys, numpy as np
Z=json.load(open(os.path.expanduser(sys.argv[1])))
T=json.load(open(os.path.expanduser(sys.argv[2])))["results"]
zc={(r["pair"],r["form"]):r["p"] for r in Z}
tc={(r["combo"],r["form"]):r for r in T}
keys=sorted(set(zc)&set(tc))
print("matched tests: %d  (counted %d, tail-fit %d)"%(len(keys),len(zc),len(tc)))
if len(keys)<200:
    print("FAIL -- too few matched tests to validate"); sys.exit(1)
zf=1/60001.
rows=[(k,zc[k],tc[k]["p"],tc[k]["method"]) for k in keys]
def band(lo,hi,label):
    sel=[(z,t) for _,z,t,_ in rows if lo<=z<hi]
    if not sel: print("  %-30s (none)"%label); return
    lr=np.array([np.log10(t/z) for z,t in sel])
    print("  %-30s n=%4d  median %+.3f dex  90%% |dev| %.2f"%(label,len(sel),np.median(lr),np.percentile(np.abs(lr),90)))
print("\nagreement by counted-p band, log10(tailfit/counted):")
band(1e-1,1.01,"p > 0.1  (bulk)"); band(1e-2,1e-1,"0.01 - 0.1")
band(1e-3,1e-2,"0.001 - 0.01");    band(10*zf,1e-3,"10/M - 0.001")
band(0,10*zf,"p < 10/60001 (counted tail)")
meth={}
for _,_,_,m in rows: meth[m]=meth.get(m,0)+1
print("\ntail-fit method mix: %s"%meth)
res=[(z,t) for _,z,t,_ in rows if z>=10*zf]
if len(res)<100: print("FAIL -- resolved range too small"); sys.exit(1)
z_=np.array([a for a,_ in res]); t_=np.array([b for _,b in res]); lr=np.log10(t_/z_)
rk=np.corrcoef(np.argsort(np.argsort(z_)),np.argsort(np.argsort(t_)))[0,1]
print("\nRESOLVED-RANGE VERDICT (counted p >= 10/60001, n=%d):"%len(res))
print("   median bias %+.3f dex   IQR %.3f   90%% |dev| %.3f   max %.3f"
      %(np.median(lr),np.subtract(*np.percentile(lr,[75,25])),np.percentile(np.abs(lr),90),np.abs(lr).max()))
print("   rank correlation %.4f"%rk)
ok = abs(np.median(lr))<0.15 and np.percentile(np.abs(lr),90)<0.5 and rk>0.95
alpha=0.05/len(rows)
print("\nBonferroni alpha for %d tests: %.3e"%(len(rows),alpha))
for lbl,src in (("counted",[(k,z) for k,z,_,_ in rows]),("tail-fit",[(k,t) for k,_,t,_ in rows])):
    print("  %-9s survivors: %d"%(lbl,sum(1 for _,p in src if p<alpha)))
    for k,p in sorted(src,key=lambda x:x[1])[:5]:
        print("      %-44s %-14s %.3e"%(k[0][:44],k[1],p))
print("\n=== %s ==="%("PASS - tail fit validated, proceed to triples" if ok else "FAIL - do not use tail fit on triples"))
sys.exit(0 if ok else 1)

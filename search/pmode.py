"""
THE P-MODE SEARCH -- search 34.

WHY IT IS RUN ON EXIS AND NOT ON GONG OR BiSON. The canonical helioseismology
datasets are not openly downloadable right now: the BiSON portal returns 403/404 on
every data path and its only Zenodo record is a PDF; GONG sits behind an interactive
query form and serves full-disk FITS images rather than a disk-integrated series;
and PMOD's VIRGO archive is FTP-only, with ftp.pmodwrc.ch port 21 timing out from
two separate networks. LISIRD has irradiance but nothing faster than six-hourly.

What IS in hand is 1,945 days of GOES-16 EUVS at one minute and 1,121 of GOES-17,
and the earlier fast search already left a loose end: in the 4-6 minute band, Mg II
returned R = 41.6 against a threshold of 24.7, and nobody followed it up.

THE TEST THAT SETTLES IT. Solar p-modes are not one line, they are a COMB: modes of
consecutive radial order are separated by the large frequency separation, dnu = 135
uHz for the Sun, under an envelope peaking at nu_max = 3.09 mHz. So the decisive
signature is not a peak at 3.09 mHz -- lots of things could sit there -- but 135 uHz
PERIODICITY IN THE SPECTRUM ITSELF. Autocorrelate the power spectrum across the
2.5-4.0 mHz band and a real p-mode comb produces a peak at exactly 135 uHz. Nothing
instrumental has a reason to.

This matters beyond the channel. Section 16 has listed p-modes as an unsearched
public channel for four revisions, and search 33 used 1/nu_max and 1/dnu as a
candidate time base on the argument that the star's own clock needs no shared units.
If that clock is visible in EXIS at all, this is where it shows.

  gate 1  the envelope must peak near 3.09 mHz
  gate 2  the spectrum's autocorrelation must peak at 135 uHz
  gate 3  both must reproduce on GOES-17, a different spacecraft
"""
import os
import numpy as np
from scipy.ndimage import median_filter

D=os.path.expanduser("~")
CH=["irr_1175","irr_1216","irr_1335","irr_1405","MgII_EXIS"]     # EUVS-B and C only
NICE={"irr_1175":"117.5nm","irr_1216":"Ly-alpha","irr_1335":"133.5nm",
      "irr_1405":"140.5nm","MgII_EXIS":"MgII"}
NU_MAX=3.09e-3; DNU=135e-6

def load(tok):
    z=np.load(D+"/euvs1m_%s.npz"%tok); return {c:z[c] for c in CH}

def prep(x, wmin=60):
    """a 60-minute detrend, not 1440 -- p-modes live at 5 minutes and a day-long
    window leaves the whole low-frequency forest sitting under them"""
    ok=np.isfinite(x)&(x>0)
    y=np.full(len(x),np.nan); y[ok]=np.log(x[ok])
    med=np.nanmedian(y); s=1.4826*np.nanmedian(np.abs(y-med))
    y=np.clip(y,med-5*s,med+5*s)
    a=np.nan_to_num(y); m=ok.astype(np.float64)
    ca=np.concatenate([[0],np.cumsum(a)]); cm=np.concatenate([[0],np.cumsum(m)])
    h=wmin//2; i=np.arange(len(y)); lo=np.maximum(0,i-h); hi=np.minimum(len(y),i+h)
    num=ca[hi]-ca[lo]; den=cm[hi]-cm[lo]
    trend=np.where(den>5,num/np.maximum(den,1),med)
    return np.where(ok,y-trend,0.0), ok

def spec(r):
    n=len(r)//2*2; w=np.blackman(n)
    P=np.abs(np.fft.rfft(r[:n]*w))**2
    f=np.fft.rfftfreq(n,d=60.0)
    return f[1:],P[1:]

def envelope(f,P,lo=1.0e-3,hi=6.0e-3,nb=60):
    """smoothed power across the oscillation band, continuum-normalised"""
    m=(f>=lo)&(f<=hi)
    fb=f[m]; Pb=P[m]
    C=np.exp(median_filter(np.log(np.maximum(Pb,1e-300)),size=20001,mode="nearest"))
    R=Pb/C
    edges=np.linspace(lo,hi,nb+1)
    cen=0.5*(edges[1:]+edges[:-1])
    prof=np.array([R[(fb>=edges[i])&(fb<edges[i+1])].mean() for i in range(nb)])
    return fb,R,cen,prof

def comb_acf(fb,R,lo=2.5e-3,hi=4.0e-3):
    """autocorrelation of the spectrum across the mode band. A p-mode comb of
    spacing dnu puts a peak here at dnu; nothing instrumental has a reason to."""
    m=(fb>=lo)&(fb<=hi)
    x=R[m]-R[m].mean()
    df=fb[1]-fb[0]
    A=np.fft.irfft(np.abs(np.fft.rfft(x))**2, n=len(x))
    A=A[:len(A)//2]/A[0]
    lags=np.arange(len(A))*df
    return lags,A

print("P-MODE SEARCH -- GOES-R EXIS, 1-minute, EUVS-B and C\n")
out={}
for tok in ("g16","g17"):
    try: G=load(tok)
    except FileNotFoundError: continue
    print("  %s"%tok)
    for c in CH:
        r,ok=prep(G[c])
        if ok.mean()<0.3: continue
        f,P=spec(r)
        fb,R,cen,prof=envelope(f,P)
        j=np.argmax(prof)
        # gate 2: comb spacing
        lags,A=comb_acf(fb,R)
        w=(lags>80e-6)&(lags<200e-6)
        k=np.argmax(A[w]); lag_pk=lags[w][k]; a_pk=A[w][k]
        base=np.std(A[(lags>250e-6)&(lags<600e-6)])
        out[(tok,c)]=(cen[j],prof[j],lag_pk,a_pk,a_pk/base if base>0 else np.nan)
        print("     %-9s envelope peak %6.3f mHz (R=%5.2f) | ACF peak at %6.1f uHz, "
              "%.3f = %5.1f sigma"%(NICE[c],cen[j]*1e3,prof[j],lag_pk*1e6,a_pk,
                                    a_pk/base if base>0 else np.nan),flush=True)

print("\n  expected: envelope near %.2f mHz, ACF peak at %.0f uHz"%(NU_MAX*1e3,DNU*1e6))
print("\n  GATE 3 -- does any channel show BOTH gates on BOTH spacecraft?")
hit=False
for c in CH:
    a=out.get(("g16",c)); b=out.get(("g17",c))
    if not a or not b: continue
    g1=abs(a[0]-NU_MAX)<0.5e-3 and abs(b[0]-NU_MAX)<0.5e-3
    g2=abs(a[2]-DNU)<20e-6 and abs(b[2]-DNU)<20e-6 and a[4]>4 and b[4]>4
    if g1 or g2:
        hit=True
        print("     %-9s envelope %s   comb %s"%(NICE[c],"YES" if g1 else "no","YES" if g2 else "no"))
if not hit:
    print("     none -- no channel shows the p-mode envelope or the 135 uHz comb on both")

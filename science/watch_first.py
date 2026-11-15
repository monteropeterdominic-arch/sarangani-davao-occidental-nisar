import tifffile, numpy as np, zlib, math, json
from PIL import Image
import os
R=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')+'/'
def read_fp(f):
    with tifffile.TiffFile(f) as t:
        p=t.pages[0]; H,W=p.shape[:2]; S=p.samplesperpixel; tw,th=p.tilewidth,p.tilelength; fh=t.filehandle
        out=np.full((H,W,S),np.nan,np.float32); ntx=-(-W//tw)
        for i,(off,cnt) in enumerate(zip(p.dataoffsets,p.databytecounts)):
            if cnt==0: continue
            fh.seek(off); raw=np.frombuffer(zlib.decompress(fh.read(cnt)),np.uint8).reshape(th, tw*S*4)
            raw=np.cumsum(raw.reshape(th,-1,S),axis=1,dtype=np.uint64).astype(np.uint8).reshape(th,-1)  # undo byte differencing, stride = samples per pixel
            wc=tw*S; b=np.stack([raw[:,(3-k)*wc:(4-k)*wc] for k in range(4)],axis=-1)  # little-endian bytes
            v=b.copy().view(np.float32).reshape(th,tw,S)
            ty,tx=divmod(i,ntx); y0,x0=ty*th,tx*tw; hh=min(th,H-y0); ww=min(tw,W-x0)
            out[y0:y0+hh,x0:x0+ww]=v[:hh,:ww]
        return out
sens=read_fp(R+'outlook/Outlook_EVI_ElNino_sensitivity_per_degC_250m.tif')
print('sens',sens.shape, np.nanpercentile(sens,[5,50,95]))
# months Nov..Oct -> Nov-Apr = bands 0..5
s=np.nanmean(sens[:,:,0:6],axis=2)
# NISAR HV Oct aggregated to 250m grid
hv=tifffile.imread(R+'nisar/backscatter_by_month/NISAR_HV_dB_2026-10_month_30m.tif').astype('float32'); hv=np.where(hv==-32768,np.nan,hv*0.01)
H,W=hv.shape; rr,cc=np.mgrid[0:H,0:W]; x=648150+30*(cc+0.5); y=735870-30*(rr+0.5)
oc=((x-566500)//250).astype(int); orr=((729000-y)//250).astype(int)
ok=np.isfinite(hv)&(oc>=0)&(oc<s.shape[1])&(orr>=0)&(orr<s.shape[0]); idx=orr[ok]*s.shape[1]+oc[ok]
lin=10**(hv[ok]/10); sm=np.bincount(idx,lin,s.size); n=np.bincount(idx,None,s.size)
hv250=np.where(n>20,10*np.log10(np.where(n>0,sm/np.maximum(n,1),np.nan)),np.nan).reshape(s.shape)
both=np.isfinite(s)&np.isfinite(hv250)
s_q=np.nanpercentile(s[both],33.3); hv_q=np.nanpercentile(hv250[both],33.3)
model_hi=both&(s<=s_q); radar_thin=both&(hv250<=hv_q)
cls=np.zeros(s.shape,np.uint8); cls[both]=1; cls[model_hi]=2; cls[model_hi&radar_thin]=3
print('sens Nov-Apr worst third <= %.4f per degC; HV thin third <= %.2f dB'%(s_q,hv_q))
print('cells: valid',both.sum(),' El Nino hits hard',model_hi.sum(),' + radar thin plants (watch first)',(model_hi&radar_thin).sum(), ' share of land watch-first %.1f%%'%(100*(model_hi&radar_thin).sum()/both.sum()))
# chance check: expected overlap if independent = 1/9 of valid
print('expected if unrelated %.1f%%'%(100/9))
np.savez('watch.npz',cls=cls,s=s,hv250=hv250)

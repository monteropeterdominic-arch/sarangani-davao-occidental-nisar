
import json, time, os, glob, traceback, re
import numpy as np
from osgeo import gdal, ogr
from scipy import ndimage as ndi
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; ML=W+'nisar/ml30/'; O=W+'nisar/products/'; os.makedirs(O,exist_ok=True)
LOG=W+'_nisar_prod.json'; st={'step':'wait','err':None,'fin':False,'stats':{}}
def sv(): json.dump(st,open(LOG,'w'),indent=1)
sv()
try:
    while not (json.load(open(W+'_nisar_ml.json'))['fin'] and os.path.exists(W+'_h5ml.json') and json.load(open(W+'_h5ml.json')).get('fin')): time.sleep(10)
    g=json.load(open(W+'_ts.json'))['notes']['grid30']; x0,y1,Wd,Hd=g['x0'],g['y1'],g['W'],g['H']; GT=(x0,30,0,y1,0,-30)
    prj=gdal.Open(glob.glob(ML+'*.tif')[0]).GetProjection()
    M=np.load(W+'tmp/aoi_mask30.npy')
    def save(name,a,mask=True):
        a=a.astype(np.float32).copy()
        if mask: a[~M]=np.nan
        o=gdal.GetDriverByName('GTiff').Create(O+name+'.tif',Wd,Hd,1,gdal.GDT_Float32,options=['COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES'])
        o.SetGeoTransform(GT); o.SetProjection(prj); b=o.GetRasterBand(1); b.SetNoDataValue(float('nan')); b.WriteArray(a); o=None
    def lin(key):
        a=gdal.Open(ML+key+'.tif').ReadAsArray().astype(np.float32); a[~(np.isfinite(a)&(a>0))]=np.nan
        # drop zero-fill / below noise floor values
        a[10*np.log10(np.where(a>0,a,np.nan))<-35]=np.nan
        return a
    def edge_ok(a,km=3):
        v=np.isfinite(a); sm=v[::10,::10]; d=ndi.distance_transform_edt(sm)*0.3
        dd=np.kron(d,np.ones((10,10)))[:Hd,:Wd]; dd=np.pad(dd,((0,max(0,Hd-dd.shape[0])),(0,max(0,Wd-dd.shape[1]))),mode='edge')
        return v&(dd>=km)
    keys=sorted(os.path.basename(f)[:-4] for f in glob.glob(ML+'*.tif'))
    st['keys']=keys; sv()
    # ---- per-date mosaics (mean of overlapping frames, linear power) in dB
    dates=sorted(set(k[:8] for k in keys))
    for pol in ('HVHV','HHHH'):
        for d in dates:
            ks=[k for k in keys if k.startswith(d) and k.endswith(pol)]
            if not ks: continue
            st['step']=f'mosaic {d} {pol}'; sv()
            acc=np.zeros((Hd,Wd),np.float64); n=np.zeros((Hd,Wd),np.int16)
            for k in ks:
                a=lin(k); ok=edge_ok(a,1); acc+=np.where(ok,a,0); n+=ok
            m=np.where(n>0,acc/np.maximum(n,1),np.nan)
            db=10*np.log10(m); save(f'NISAR_{pol[:2]}_dB_{d}_30m',db)
            st['stats'][f'{pol[:2]}_{d}']=dict(cover_pct=round(100*float(np.isfinite(db)[M].mean()),1),median_dB=round(float(np.nanmedian(db[M])),2) if np.isfinite(db[M]).any() else None); sv()
    # ---- same-track change pairs
    pairs=[('T039A004','HVHV','20260920','20261014'),('T075D085','HVHV','20260923','20261017'),('T075D086','HHHH','20260923','20261017'),('T039A005','HHHH','20260628','20260920'),
           ('T039A004','HHHH','20260920','20261014'),('T075D085','HHHH','20260923','20261017'),('T039A005','HVHV','20260628','20261014'),('T039A005','HHHH','20260628','20261014'),('T039A005','HHHH','20260920','20261014'),('T111A005','HVHV','20260925','20261019'),('T111A005','HHHH','20260925','20261019')]
    made=[]
    for tr,pol,a_,b_ in pairs:
        ka=[k for k in keys if k.startswith(a_) and f'_{tr}_' in k and k.endswith(pol)]; kb=[k for k in keys if k.startswith(b_) and f'_{tr}_' in k and k.endswith(pol)]
        if not ka or not kb: st['stats'][f'pair {tr} {pol} {a_}-{b_}']='missing'; continue
        st['step']=f'pair {tr} {pol}'; sv()
        A=lin(ka[0]); Bb=lin(kb[0]); ok=edge_ok(A,3)&edge_ok(Bb,3)
        dB=np.where(ok,10*np.log10(Bb)-10*np.log10(A),np.nan)
        v=dB[M&np.isfinite(dB)]
        if v.size<1000: st['stats'][f'pair {tr} {pol} {a_}-{b_}']='no overlap in AOI'; continue
        med=float(np.median(v)); sig=float(1.4826*np.median(np.abs(v-med)))
        name=f'NISAR_{pol[:2]}_Change_{tr}_{a_}_to_{b_}_30m'; save(name,dB); made.append(name)
        st['stats'][name]=dict(cover_pct=round(100*float(np.isfinite(dB)[M].mean()),1),median=round(med,2),robust_sigma=round(sig,2),p1=round(float(np.percentile(v,1)),2),p99=round(float(np.percentile(v,99)),2),
            strong_drop_pct=round(100*float((v<med-3*sig).mean()),2),strong_rise_pct=round(100*float((v>med+3*sig).mean()),2)); sv()
    st['made']=made; st['step']='done'; st['fin']=True; sv()
except Exception:
    st['err']=traceback.format_exc(); sv()

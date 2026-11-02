import glob, os, re, json, traceback, numpy as np
from osgeo import gdal
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; NP=W+'nisar/products/'
ST={'fin':False}
def sv(): json.dump(ST,open(W+'_month_chg.json','w'),indent=1)
try:
    M=np.load(W+'tmp/aoi_mask30.npy')
    for pol in ('HH','HV'):
        fs=[f for f in glob.glob(NP+f'NISAR_{pol}_Change_*_2026092?_to_202610??_30m.tif')+glob.glob(NP+f'NISAR_{pol}_Change_*_2026093?_to_202610??_30m.tif')]
        acc=None
        for f in fs:
            d=gdal.Open(f); a=d.ReadAsArray().astype(np.float32)
            if acc is None: acc=np.zeros_like(a,np.float64); n=np.zeros(a.shape,np.uint8); gt=d.GetGeoTransform(); prj=d.GetProjection()
            ok=np.isfinite(a); acc[ok]+=a[ok]; n[ok]+=1
        c=np.where(n>0,acc/np.maximum(n,1),np.nan).astype(np.float32); c[~M]=np.nan
        out=NP+f'NISAR_{pol}_Change_month_2026-09_to_2026-10_30m.tif'
        o=gdal.GetDriverByName('GTiff').Create(out,c.shape[1],c.shape[0],1,gdal.GDT_Float32,['COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES'])
        o.SetGeoTransform(gt); o.SetProjection(prj); b=o.GetRasterBand(1); b.SetNoDataValue(float('nan')); b.WriteArray(c); o=None
        v=c[M&np.isfinite(c)]
        ST[pol]=dict(pairs=[os.path.basename(f) for f in fs],cover_pct=round(100*float(np.isfinite(c)[M].mean()),1),median=round(float(np.median(v)),2),pct_drop_gt2=round(100*float((v<-2).mean()),1),pct_rise_gt2=round(100*float((v>2).mean()),1)); sv()
    S2=open(W+'web_export_nisar2.py').read(); head=S2.split("\nRAMP['hv']")[0]
    H={'__name__':'h'}; exec(compile(head.replace("LOG=W+'_web_nisar2.json'","LOG=W+'_web_mchg_inner.json'"),'h','exec'),H)
    for pol in ('HH','HV'):
        H['export'](NP+f'NISAR_{pol}_Change_month_2026-09_to_2026-10_30m.tif',f'nisarmchg_{pol.lower()}','nisarmchg',f'NISAR {pol} change, September → October 2026',60,H['RAMP']['div'],-4,4,'2026-09-01','2026-10-31',unit=' dB')
    json.dump(H['manifest']['layers'],open(W+'webmap/data/mchg_layers.json','w'))
    ST['fin']=True; sv()
except Exception:
    ST['err']=traceback.format_exc()[-1500:]; sv()

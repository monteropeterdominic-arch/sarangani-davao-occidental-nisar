
import numpy as np, glob, os, json, re, traceback, time
from osgeo import gdal
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; NP=W+'nisar/products/'; ML=W+'ml/'
LOG=W+'_month.json'; ST={'step':'wait','err':None,'fin':False,'res':{}}
def sv(): json.dump(ST,open(LOG,'w'),indent=1)
sv()
try:
    while True:
        try:
            if json.load(open(W+'_web_nisar2.json'))['fin']: break
        except Exception: pass
        time.sleep(5)
    g=json.load(open(W+'_ts.json'))['notes']['grid30']; Wd,Hd=g['W'],g['H']; GT=(g['x0'],30,0,g['y1'],0,-30)
    M=np.load(W+'tmp/aoi_mask30.npy'); prj=gdal.Open(ML+'srtm_dem_30m.tif').GetProjection()
    def save(path,a):
        a=a.copy(); a[~M]=np.nan
        o=gdal.GetDriverByName('GTiff').Create(path,Wd,Hd,1,gdal.GDT_Float32,options=['COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES']); o.SetGeoTransform(GT); o.SetProjection(prj); b=o.GetRasterBand(1); b.SetNoDataValue(float('nan')); b.WriteArray(a); o=None
    MON={'2026-06':['20260628'],'2026-09':['20260920','20260923','20260925','20260930'],'2026-10':['20261014','20261017','20261019']}
    for mo,ds in MON.items():
        for pol in ('HV','HH'):
            ST['step']=f'{mo} {pol}'; sv()
            acc=np.zeros((Hd,Wd),np.float64); n=np.zeros((Hd,Wd),np.uint8)
            for d in ds:
                p=NP+f'NISAR_{pol}_dB_{d}_30m.tif'
                if not os.path.exists(p): continue
                a=gdal.Open(p).ReadAsArray().astype(np.float32); ok=np.isfinite(a); acc[ok]+=10**(a[ok]/10); n[ok]+=1; del a,ok
            c=np.where(n>0,10*np.log10(acc/np.maximum(n,1)),np.nan).astype(np.float32); save(NP+f'NISAR_{pol}_dB_{mo}_month_30m.tif',c)
            ST['res'][f'{pol} {mo}']=dict(cover_pct=round(100*float(np.isfinite(c)[M].mean()),1),passes=len(ds)); sv(); del acc,n,c
        for t in ('EVI','NDMI'):
            acc=np.zeros((Hd,Wd),np.float32); n=np.zeros((Hd,Wd),np.uint8)
            for d in ds:
                p=ML+f'pred/NISAR_estimated_{t}_{d}_30m.tif'
                if not os.path.exists(p): continue
                a=gdal.Open(p).ReadAsArray().astype(np.float32); ok=np.isfinite(a); acc[ok]+=a[ok]; n[ok]+=1; del a,ok
            c=np.where(n>0,acc/np.maximum(n,1),np.nan).astype(np.float32); save(ML+f'pred/NISAR_estimated_{t}_{mo}_month_30m.tif',c)
            ST['res'][f'ML {t} {mo}']=round(100*float(np.isfinite(c)[M].mean()),1); sv(); del acc,n,c
    # web export of monthly layers
    ST['step']='web export'; sv()
    S2=open(W+'web_export_nisar2.py').read()
    head=S2.split("\nRAMP['hv']")[0]
    H={'__name__':'h'}; exec(compile(head.replace("LOG=W+'_web_nisar2.json'","LOG=W+'_web_month.json'"),'h','exec'),H); export=H['export']; RAMP=H['RAMP']; DATA=H['DATA']; manifest=H['manifest']
    RAMP['hv']=['#33505f','#3f7f73','#73ad6c','#b8d58c','#f1f0cf']; RAMP['hh']=['#47507a','#5f7fa3','#8fb1bd','#d3dfc4','#fbf5e3']
    import calendar
    new=[]
    for mo,ds in MON.items():
        y,m=map(int,mo.split('-')); s,e=f'{mo}-01',f'{mo}-{calendar.monthrange(y,m)[1]:02d}'; nm=calendar.month_name[m]
        for pol,ramp,lo,hi in (('HV','hv',-22,-8),('HH','hh',-14,-2)):
            export(NP+f'NISAR_{pol}_dB_{mo}_month_30m.tif',f'nisarm_{pol.lower()}_{mo}','nisarm_'+pol.lower(),f'NISAR {pol} backscatter, {nm} 2026 ({len(ds)} pass{"es" if len(ds)>1 else ""})',60,RAMP[ramp],lo,hi,s,e,unit=' dB')
        for t,ramp,lo,hi,lab in (('EVI','ndvi',0,0.8,'vegetation (EVI)'),('NDMI','div',-0.2,0.6,'plant moisture (NDMI)')):
            export(ML+f'pred/NISAR_estimated_{t}_{mo}_month_30m.tif',f'mlm_{t.lower()}_{mo}','mlm_'+t.lower(),f'NISAR-estimated {lab}, {nm} 2026',60,RAMP[ramp],lo,hi,s,e)
    prev=json.load(open(DATA+'manifest.json'))
    prev['layers']=[l for l in prev['layers'] if not (l['id'].startswith('nisarm_') or l['id'].startswith('mlm_'))]+manifest['layers']
    json.dump(prev,open(DATA+'manifest.json','w'),indent=1,allow_nan=False)
    ST['step']='done'; ST['fin']=True; sv()
except Exception:
    ST['err']=traceback.format_exc(); sv()

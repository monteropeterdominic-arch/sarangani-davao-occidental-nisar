
import numpy as np, os, glob, json, re, datetime as dt, traceback, warnings, csv
from osgeo import gdal, ogr
gdal.UseExceptions(); warnings.filterwarnings('ignore')
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; D='C:/Users/PETER DOMINIC/Downloads/'; P='C:/Users/PETER DOMINIC/Pictures/'; OUT=W+'out_v2/'
LOG=W+'_ts2_v2.json'; ST={'step':'start','err':None,'fin':False}
def sv(): json.dump(ST,open(LOG,'w'))
COF=['COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES']
try:
    g=json.load(open(W+'_ts_v2.json'))['notes']['grid30']; x0=g['x0']; y1=g['y1']; x1=x0+g['W']*30; y0=y1-g['H']*30
    aoi=W+'boundary/AOI_DavaoOccidental_Sarangani_EPSG32651.gpkg'
    prj=gdal.Open(glob.glob(D+'*EVI_doy2025024*')[0]).GetProjection()
    def mask_for(gt,w,h):
        m=gdal.GetDriverByName('MEM').Create('',w,h,1,gdal.GDT_Byte); m.SetGeoTransform(gt); m.SetProjection(prj)
        v=ogr.Open(aoi); gdal.RasterizeLayer(m,[1],v.GetLayer(0),burn_values=[1],options=['ALL_TOUCHED=TRUE']); a=m.ReadAsArray()==1; v=None; return a
    def save(path,a,gt,w,h,mask):
        a=a.astype(np.float32).copy(); a[~mask]=np.nan
        o=gdal.GetDriverByName('GTiff').Create(path,w,h,1,gdal.GDT_Float32,options=COF); o.SetGeoTransform(gt); o.SetProjection(prj)
        b=o.GetRasterBand(1); b.SetNoDataValue(float('nan')); b.WriteArray(a); o.FlushCache(); o=None
    g250=(x0,250,0,y1,0,-250); W250=int(np.ceil((x1-x0)/250)); H250=int(np.ceil((y1-y0)/250)); M250=mask_for(g250,W250,H250)
    g1k=(x0,1000,0,y1,0,-1000); W1k=int(np.ceil((x1-x0)/1000)); H1k=int(np.ceil((y1-y0)/1000)); M1k=mask_for(g1k,W1k,H1k)
    def warp(f,nd):
        d=gdal.Warp('',f,format='MEM',dstSRS='EPSG:32651',outputBounds=(g1k[0],g1k[3]-H1k*1000,g1k[0]+W1k*1000,g1k[3]),width=W1k,height=H1k,resampleAlg='near',srcNodata=nd,dstNodata=nd)
        a=d.GetRasterBand(1).ReadAsArray(); d=None; return a
    lm={}
    files={}
    for f in glob.glob(D+'MOD11A2.061_LST_Day_1km_doy*.tif'):
        m=re.search(r'doy(\d{4})(\d{3})',f); files[(int(m.group(1)),int(m.group(2)))]=f
    for (y,j),f in sorted(files.items()):
        ST['step']=f'LST {y}-{j}'; sv()
        d=dt.date(y,1,1)+dt.timedelta(j-1+4); key=f'{d.year}-{d.month:02d}'
        q=warp(f.replace('_LST_Day_1km_','_QC_Day_'),None).astype(np.int32)
        v=warp(f,0).astype(np.float32)
        good=(v>0)&((q&3)<=1)&(((q>>6)&3)<=1)
        lm.setdefault(key,[]).append(np.where(good,v*0.02-273.15,np.nan))
    lst={}
    for k,L in sorted(lm.items()):
        a=np.nanmean(np.stack(L),axis=0); lst[k]=a; save(OUT+f'modis_lst/MODIS_LST_Day_{k}_1km.tif',a,g1k,W1k,H1k,M1k)
    ndvi={}
    for f in glob.glob(OUT+'modis_ndvi/MODIS_NDVI_20??-??_250m.tif'):
        k=re.search(r'NDVI_(\d{4}-\d{2})_250m',f).group(1); ndvi[k]=gdal.Open(f).ReadAsArray().astype(np.float32)
    for m in range(1,13):
        a,b=f'2026-{m:02d}',f'2025-{m:02d}'
        if a in ndvi and b in ndvi: save(OUT+f'modis_ndvi/MODIS_NDVI_Anomaly_{a}_minus_{b}_250m.tif',ndvi[a]-ndvi[b],g250,W250,H250,M250)
        if a in lst and b in lst: save(OUT+f'modis_lst/MODIS_LST_Anomaly_{a}_minus_{b}_1km.tif',lst[a]-lst[b],g1k,W1k,H1k,M1k)
    rows=[]
    for k in sorted(set(ndvi)|set(lst)):
        rows.append(dict(month=k,ndvi_mean=(round(float(np.nanmean(ndvi[k][M250])),4) if k in ndvi else None),ndvi_cover_pct=(round(100*float(np.isfinite(ndvi[k])[M250].mean()),1) if k in ndvi else None),
            lst_day_mean_C=(round(float(np.nanmean(lst[k][M1k])),2) if k in lst else None),lst_cover_pct=(round(100*float(np.isfinite(lst[k])[M1k].mean()),1) if k in lst else None)))
    with open(OUT+'AOI_monthly_means.csv','w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    ST['rows']=rows; ST['step']='done'; ST['fin']=True; sv()
except Exception:
    ST['err']=traceback.format_exc(); sv()

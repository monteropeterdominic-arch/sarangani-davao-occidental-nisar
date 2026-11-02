
import numpy as np, os, glob, json, time, traceback, re, datetime as dt, warnings
from osgeo import gdal, ogr
gdal.UseExceptions(); warnings.filterwarnings('ignore')
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; D='C:/Users/PETER DOMINIC/Downloads/'; P='C:/Users/PETER DOMINIC/Pictures/'
OUT=W+'out_v2/'; TMP=W+'tmp_v2/'
for d in (OUT,OUT+'hls',OUT+'modis_ndvi',OUT+'modis_lst',OUT+'earthquake',TMP): os.makedirs(d,exist_ok=True)
LOG=W+'_ts_v2.json'; ST={'step':'wait','err':None,'fin':False,'notes':{}}
def sv(): json.dump(ST,open(LOG,'w'),indent=1)
CO=['COMPRESS=DEFLATE','PREDICTOR=2','TILED=YES']; COF=['COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES']
def doy2date(y,j): return dt.date(y,1,1)+dt.timedelta(j-1)
def run():
  try:
    while not os.path.exists(W+'_aoi.json') or not json.load(open(W+'_aoi.json'))['fin']:
        if os.path.exists(W+'_aoi.json') and json.load(open(W+'_aoi.json'))['err']: raise RuntimeError('AOI failed')
        time.sleep(5)
    aoi=W+'boundary/AOI_DavaoOccidental_Sarangani_EPSG32651.gpkg'
    ds=ogr.Open(aoi); l=ds.GetLayer(0); xmin,xmax,ymin,ymax=l.GetExtent(); ds=None
    # ---- 30 m grid snapped to the HLS pixel grid
    ref=gdal.Open(sorted(glob.glob(D+'*EVI_doy2025024*'))[0]); rgt=ref.GetGeoTransform(); prj=ref.GetProjection()
    x0=rgt[0]+np.floor((xmin-2000-rgt[0])/30)*30; x1=rgt[0]+np.ceil((xmax+2000-rgt[0])/30)*30
    y1=rgt[3]-np.floor((rgt[3]-(ymax+2000))/30)*30; y0=rgt[3]-np.ceil((rgt[3]-(ymin-2000))/30)*30
    Wd=int(round((x1-x0)/30)); Hd=int(round((y1-y0)/30)); GT=(x0,30,0,y1,0,-30)
    ST['notes']['grid30']=dict(x0=x0,y1=y1,W=Wd,H=Hd); sv()
    def mask_for(gt,w,h):
        m=gdal.GetDriverByName('MEM').Create('',w,h,1,gdal.GDT_Byte); m.SetGeoTransform(gt); m.SetProjection(prj)
        v=ogr.Open(aoi); gdal.RasterizeLayer(m,[1],v.GetLayer(0),burn_values=[1],options=['ALL_TOUCHED=TRUE']); a=m.ReadAsArray()==1; v=None; return a
    M30=mask_for(GT,Wd,Hd); ST['notes']['aoi_px30']=int(M30.sum())
    def save(path,a,gt,w,h,dtype=gdal.GDT_Float32,nd=np.nan,mask=None):
        a=a.copy()
        if mask is not None: a[~mask]=nd
        o=gdal.GetDriverByName('GTiff').Create(path,w,h,1,dtype,options=COF if dtype==gdal.GDT_Float32 else CO)
        o.SetGeoTransform(gt); o.SetProjection(prj); b=o.GetRasterBand(1); b.SetNoDataValue(float(nd) if dtype==gdal.GDT_Float32 else nd); b.WriteArray(a); o.FlushCache(); o=None
    # ---- HLS index files
    def hls_files(idx):
        out={}
        for f in glob.glob(D+f'HLSL30_VI.020_{idx}_doy*_51N.tif'):
            m=re.search(r'doy(\d{4})(\d{3})',f); out[(int(m.group(1)),int(m.group(2)))]=f
        return dict(sorted(out.items()))
    def read_on_grid(f):
        d=gdal.Open(f); gt=d.GetGeoTransform(); b=d.GetRasterBand(1)
        a=np.full((Hd,Wd),np.nan,np.float32)
        c0=int(round((x0-gt[0])/30)); r0=int(round((gt[3]-y1)/30))
        cs,rs=max(c0,0),max(r0,0); ce=min(c0+Wd,d.RasterXSize); re_=min(r0+Hd,d.RasterYSize)
        if ce<=cs or re_<=rs: return a
        x=b.ReadAsArray(cs,rs,ce-cs,re_-rs).astype(np.float32); x[x==-19999]=np.nan; x*=0.0001
        x[(x<-1)|(x>1)]=np.nan
        a[rs-r0:re_-r0, cs-c0:ce-c0]=x; return a
    def read_rows(d,r,n):
        gt=d.GetGeoTransform(); b=d.GetRasterBand(1); a=np.full((n,Wd),np.nan,np.float32)
        c0=int(round((x0-gt[0])/30)); r0=int(round((gt[3]-y1)/30))+r
        cs,rs=max(c0,0),max(r0,0); ce=min(c0+Wd,d.RasterXSize); re_=min(r0+n,d.RasterYSize)
        if ce<=cs or re_<=rs: return a
        x=b.ReadAsArray(cs,rs,ce-cs,re_-rs).astype(np.float32); x[x==-19999]=np.nan; x*=0.0001; x[(x<-1)|(x>1)]=np.nan
        a[rs-r0:re_-r0, cs-c0:ce-c0]=x; return a
    Q={'2025Q1':((2025,1),(2025,90)),'2025Q2':((2025,91),(2025,181)),'2025Q3':((2025,182),(2025,273)),'2025Q4':((2025,274),(2025,365)),
       '2026Q1':((2026,1),(2026,90)),'2026Q2':((2026,91),(2026,181)),'2026Q3':((2026,182),(2026,273))}
    EQ={'pre2026':((2026,60),(2026,158)),'post2026':((2026,160),(2026,268)),'pre2025':((2025,60),(2025,158)),'post2025':((2025,160),(2025,268))}
    comp={}
    for idx in ('EVI','NDMI','NDWI','NBR'):
        F=hls_files(idx)
        periods=(list(Q.items()) if idx!='NBR' else [])+(list(EQ.items()) if idx in ('EVI','NBR') else [])
        for name,(a_,b_) in periods:
            ST['step']=f'HLS {idx} {name}'; sv()
            sel=[f for k,f in F.items() if a_<=k<=b_]
            if not sel: continue
            # streaming median in row strips to keep memory low
            med=np.full((Hd,Wd),np.nan,np.float32); cnt=np.zeros((Hd,Wd),np.uint8)
            dsets=[gdal.Open(f) for f in sel]
            for r in range(0,Hd,300):
                n_=min(300,Hd-r); stack=np.stack([read_rows(d_,r,n_) for d_ in dsets])
                s_=np.sort(stack,axis=0); c_=np.isfinite(s_).sum(0); lo=np.clip((c_-1)//2,0,None); hi=np.clip(c_//2,0,None)
                m_=0.5*(np.take_along_axis(s_,lo[None],0)[0]+np.take_along_axis(s_,hi[None],0)[0]); m_[c_==0]=np.nan; med[r:r+n_]=m_; cnt[r:r+n_]=c_.astype(np.uint8)
            dsets=None
            comp[(idx,name)]=med
            if name in Q:
                save(OUT+f'hls/HLS_{idx}_{name}_30m.tif',med,GT,Wd,Hd,mask=M30)
                if idx=='EVI': save(OUT+f'hls/HLS_ValidObs_{name}_30m.tif',cnt,GT,Wd,Hd,gdal.GDT_Byte,0,mask=M30)
            ST['notes'][f'{idx}_{name}']=dict(n_scenes=len(sel),cover_pct=round(100*float(np.isfinite(med)[M30].mean()),1)); sv()
        # El Nino anomaly: same quarter, 2026 minus 2025
        for q in ('Q1','Q2','Q3'):
            a=comp.get((idx,'2026'+q)); b=comp.get((idx,'2025'+q))
            if a is not None and b is not None: save(OUT+f'hls/HLS_{idx}_Anomaly_2026{q}_minus_2025{q}_30m.tif',a-b,GT,Wd,Hd,mask=M30)
        for name in Q: comp.pop((idx,name),None)
        if idx in ('NDMI','NDWI'):
            for name in EQ: comp.pop((idx,name),None)
        import gc; gc.collect()
    # ---- earthquake: 2026 change around 8 Jun minus the same-season change in 2025
    for idx in ('EVI','NBR'):
        d26=comp[(idx,'post2026')]-comp[(idx,'pre2026')]; d25=comp[(idx,'post2025')]-comp[(idx,'pre2025')]
        save(OUT+f'earthquake/HLS_{idx}_Change_2026_post_minus_pre_8Jun_30m.tif',d26,GT,Wd,Hd,mask=M30)
        save(OUT+f'earthquake/HLS_{idx}_Change_2025_sameSeason_30m.tif',d25,GT,Wd,Hd,mask=M30)
        save(OUT+f'earthquake/HLS_{idx}_Change_2026_minus_2025_seasonal_30m.tif',d26-d25,GT,Wd,Hd,mask=M30)
    for k in EQ: save(OUT+f'earthquake/HLS_EVI_{k}_30m.tif',comp[('EVI',k)],GT,Wd,Hd,mask=M30)
    np.save(TMP+'eq_evi_pre26.npy',comp[('EVI','pre2026')]); np.save(TMP+'eq_evi_post26.npy',comp[('EVI','post2026')])
    np.save(TMP+'eq_devi_excess.npy',(comp[('EVI','post2026')]-comp[('EVI','pre2026')])-(comp[('EVI','post2025')]-comp[('EVI','pre2025')]))
    json.dump(dict(GT=GT,W=Wd,H=Hd),open(TMP+'grid30.json','w'))
    comp=None
    # ---- MODIS NDVI monthly maximum-value composite (reliability 0/1), 250 m
    g250=(x0,250,0,y1,0,-250); W250=int(np.ceil((x1-x0)/250)); H250=int(np.ceil((y1-y0)/250)); M250=mask_for(g250,W250,H250)
    g1k=(x0,1000,0,y1,0,-1000); W1k=int(np.ceil((x1-x0)/1000)); H1k=int(np.ceil((y1-y0)/1000)); M1k=mask_for(g1k,W1k,H1k)
    def warp(f,gt,w,h,nd,rs='near'):
        d=gdal.Warp('',f,format='MEM',dstSRS='EPSG:32651',outputBounds=(gt[0],gt[3]-h*abs(gt[5]),gt[0]+w*gt[1],gt[3]),width=w,height=h,resampleAlg=rs,srcNodata=nd,dstNodata=nd)
        return d.GetRasterBand(1).ReadAsArray()
    nd_files={}
    for f in glob.glob(D+'MOD13Q1.061__250m_16_days_NDVI_doy*.tif'):
        m=re.search(r'doy(\d{4})(\d{3})',f); nd_files[(int(m.group(1)),int(m.group(2)))]=f
    months={}
    for (y,j),f in sorted(nd_files.items()):
        ST['step']=f'MODIS NDVI {y}-{j}'; sv()
        d=doy2date(y,j)+dt.timedelta(8); key=f'{d.year}-{d.month:02d}'   # composite mid-date
        rel=warp(f.replace('_NDVI_','_pixel_reliability_'),g250,W250,H250,-1).astype(np.int32)
        v=warp(f,g250,W250,H250,-3000).astype(np.float32); v[(v==-3000)|~np.isin(rel,[0,1])]=np.nan; v*=0.0001
        months.setdefault(key,[]).append(v)
    ndvi={}
    for k,L in sorted(months.items()):
        a=np.nanmax(np.stack(L),axis=0); ndvi[k]=a; save(OUT+f'modis_ndvi/MODIS_NDVI_{k}_250m.tif',a,g250,W250,H250,mask=M250)
    ST['step']='hls+ndvi done'; ST['fin']=True; sv(); return
    # ---- MODIS LST monthly mean of good-quality 8-day day-time LST (deg C), 1 km
    lst_files={}
    for f in glob.glob(D+'MOD11A2.061_LST_Day_1km_doy*.tif'):
        m=re.search(r'doy(\d{4})(\d{3})',f); lst_files[(int(m.group(1)),int(m.group(2)))]=f
    lm={}
    for (y,j),f in sorted(lst_files.items()):
        ST['step']=f'MODIS LST {y}-{j}'; sv()
        d=doy2date(y,j)+dt.timedelta(4); key=f'{d.year}-{d.month:02d}'
        qf=f.replace('_LST_Day_1km_','_QC_Day_')
        q=gdal.Warp('',qf,format='MEM',dstSRS='EPSG:32651',outputBounds=(g1k[0],g1k[3]-H1k*1000,g1k[0]+W1k*1000,g1k[3]),width=W1k,height=H1k,resampleAlg='near',srcNodata=None,dstNodata=255).GetRasterBand(1).ReadAsArray().astype(np.int32)
        v=warp(f,g1k,W1k,H1k,0).astype(np.float32)
        good=(v>0)&((q&3)<=1)&(((q>>6)&3)<=1)   # LST produced, LST error <= 2 K
        v=np.where(good,v*0.02-273.15,np.nan); lm.setdefault(key,[]).append(v)
    lst={}
    for k,L in sorted(lm.items()):
        a=np.nanmean(np.stack(L),axis=0); lst[k]=a; save(OUT+f'modis_lst/MODIS_LST_Day_{k}_1km.tif',a,g1k,W1k,H1k,mask=M1k)
    # ---- monthly anomalies 2026 minus 2025 (same month)
    for m in range(1,13):
        a,b=f'2026-{m:02d}',f'2025-{m:02d}'
        if a in ndvi and b in ndvi: save(OUT+f'modis_ndvi/MODIS_NDVI_Anomaly_{a}_minus_{b}_250m.tif',ndvi[a]-ndvi[b],g250,W250,H250,mask=M250)
        if a in lst and b in lst: save(OUT+f'modis_lst/MODIS_LST_Anomaly_{a}_minus_{b}_1km.tif',lst[a]-lst[b],g1k,W1k,H1k,mask=M1k)
    # area-mean table for charts
    rows=[]
    for k in sorted(set(ndvi)|set(lst)):
        rows.append(dict(month=k,ndvi_mean=(float(np.nanmean(ndvi[k][M250])) if k in ndvi else None),ndvi_cover_pct=(round(100*float(np.isfinite(ndvi[k])[M250].mean()),1) if k in ndvi else None),
                         lst_mean_C=(float(np.nanmean(lst[k][M1k])) if k in lst else None),lst_cover_pct=(round(100*float(np.isfinite(lst[k])[M1k].mean()),1) if k in lst else None)))
    json.dump(rows,open(OUT+'AOI_monthly_means.json','w'),indent=1)
    import csv
    with open(OUT+'AOI_monthly_means.csv','w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    ST['step']='done'; ST['fin']=True; sv()
  except Exception:
    ST['err']=traceback.format_exc(); sv()
run()

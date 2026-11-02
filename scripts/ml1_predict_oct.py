
import numpy as np, os, glob, re, json, time, traceback, datetime as dt, gc
from osgeo import gdal
from scipy import ndimage as ndi
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; D='C:/Users/PETER DOMINIC/Downloads/'; NP=W+'nisar/products/'; ML=W+'ml/'; os.makedirs(ML+'pred',exist_ok=True)
LOG=W+'_ml1_oct.json'; ST={'step':'start','err':None,'fin':False,'res':{}}
def sv(): json.dump(ST,open(LOG,'w'),indent=1)
try:
    g=json.load(open(W+'_ts.json'))['notes']['grid30']; x0,y1,Wd,Hd=g['x0'],g['y1'],g['W'],g['H']; GT=(x0,30,0,y1,0,-30); B=(x0,y1-Hd*30,x0+Wd*30,y1)
    prj=gdal.Open(glob.glob(NP+'*.tif')[0]).GetProjection(); M=np.load(W+'tmp/aoi_mask30.npy')
    def rd(p): return gdal.Open(p).ReadAsArray().astype(np.float32)
    ST['step']='terrain'; sv()
    if not os.path.exists(ML+'srtm_aspect_30m.tif'):
        gdal.Warp(ML+'srtm_dem_30m.tif',D+'SRTMGL1_NC.003_SRTMGL1_DEM_doy2000042000000_aid0001.tif',dstSRS='EPSG:32651',outputBounds=B,width=Wd,height=Hd,resampleAlg='bilinear',dstNodata=-32768,outputType=gdal.GDT_Float32,creationOptions=['COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES'])
        gdal.DEMProcessing(ML+'srtm_slope_30m.tif',ML+'srtm_dem_30m.tif','slope',computeEdges=True,creationOptions=['COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES'])
        gdal.DEMProcessing(ML+'srtm_aspect_30m.tif',ML+'srtm_dem_30m.tif','aspect',computeEdges=True,zeroForFlat=True,creationOptions=['COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES'])
    dem=rd(ML+'srtm_dem_30m.tif'); dem[dem<=-1000]=np.nan; slope=rd(ML+'srtm_slope_30m.tif'); slope[slope<0]=np.nan
    asp=rd(ML+'srtm_aspect_30m.tif'); asp[asp<0]=np.nan; acos=np.cos(np.radians(asp)).astype(np.float32); asin=np.sin(np.radians(asp)).astype(np.float32); del asp
    dates=sorted(re.search(r'(\d{8})',os.path.basename(f)).group(1) for f in glob.glob(NP+'NISAR_HV_dB_2026????_30m.tif'))
    dates=[d for d in dates if d!='20261019']
    ASC={'20260628':1,'20260920':1,'20260925':1,'20261014':1,'20260923':0,'20260930':0,'20261017':0,'20261019':1}
    FN=['HH_dB','HV_dB','HVminusHH_dB','RVI','HV_mean3x3','HV_std3x3','HV_mean7x7','ascending','elev','slope','aspect_cos','aspect_sin']
    NIS=list(range(8)); TER=[8,9,10,11]
    def load(d):
        hv=rd(NP+f'NISAR_HV_dB_{d}_30m.tif'); p=NP+f'NISAR_HH_dB_{d}_30m.tif'; hh=rd(p) if os.path.exists(p) else np.full_like(hv,np.nan); return hh,hv
    def feats_rows(hh,hv,d,r0,r1):
        # features for rows r0:r1 using a 3-row halo for the window statistics
        a0=max(0,r0-3); a1=min(Hd,r1+3); h=hv[a0:a1]; med=np.nanmedian(h) if np.isfinite(h).any() else -15
        hf=np.where(np.isfinite(h),h,med); m3=ndi.uniform_filter(hf,3); s3=np.sqrt(np.maximum(ndi.uniform_filter(hf*hf,3)-m3*m3,0)); m7=ndi.uniform_filter(hf,7)
        sl=slice(r0-a0,r0-a0+(r1-r0))
        HH=hh[r0:r1]; HV=hv[r0:r1]; lhh=10**(HH/10); lhv=10**(HV/10)
        F=np.stack([HH,HV,HV-HH,4*lhv/(lhh+lhv),m3[sl],s3[sl],m7[sl],np.full_like(HV,ASC.get(d,0)),dem[r0:r1],slope[r0:r1],acos[r0:r1],asin[r0:r1]],-1).astype(np.float32)
        F[~np.isfinite(HV)]=np.nan
        return F
    import joblib
    models={t:joblib.load(ML+f'model_{t}_nisar_terrain.joblib') for t in ('EVI','NDMI')}
    for d in ['20261014','20261017','20261019']:
        ST['step']=f'predict {d}'; sv()
        hh,hv=load(d)
        for tgt,m in models.items():
            P=np.full((Hd,Wd),np.nan,np.float32)
            for r0 in range(0,Hd,400):
                r1=min(Hd,r0+400); F=feats_rows(hh,hv,d,r0,r1); ok=M[r0:r1]&np.isfinite(F).all(-1)
                if ok.any(): blk=P[r0:r1]; blk[ok]=m.predict(F[ok][:,NIS+TER]); P[r0:r1]=blk
            o=gdal.GetDriverByName('GTiff').Create(ML+f'pred/NISAR_estimated_{tgt}_{d}_30m.tif',Wd,Hd,1,gdal.GDT_Float32,options=['COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES'])
            o.SetGeoTransform(GT); o.SetProjection(prj); b=o.GetRasterBand(1); b.SetNoDataValue(float('nan')); b.WriteArray(P); o=None
            ST['res'][f'{tgt}_{d}_cover_pct']=round(100*float(np.isfinite(P)[M].mean()),1); ST['res'][f'{tgt}_{d}_median']=round(float(np.nanmedian(P[M])),3) if np.isfinite(P[M]).any() else None
            del P; gc.collect()
        sv(); del hh,hv; gc.collect()
    ST['step']='done'; ST['fin']=True; sv()
except Exception:
    ST['err']=traceback.format_exc(); sv()

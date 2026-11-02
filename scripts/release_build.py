import os, glob, json, shutil, traceback, re, numpy as np
from osgeo import gdal
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; R=W+'release/'
ST={'n':0,'mb':0,'fin':False,'skip':[]}
def sv(): json.dump(ST,open(W+'_release.json','w'),indent=1)
CO=['COMPRESS=DEFLATE','ZLEVEL=9','TILED=YES','BLOCKXSIZE=512','BLOCKYSIZE=512']
def put(src,dst,scale=None,dtype='int16',desc=None):
    os.makedirs(os.path.dirname(R+dst),exist_ok=True)
    d=gdal.Open(src); nb=d.RasterCount
    if dtype=='float32':
        gdal.Translate(R+dst,d,creationOptions=CO+['PREDICTOR=3'])
    else:
        gt=gdal.GDT_Int16 if dtype=='int16' else gdal.GDT_Byte; nd=-32768 if dtype=='int16' else 255
        o=gdal.GetDriverByName('GTiff').Create(R+dst,d.RasterXSize,d.RasterYSize,nb,gt,CO+['PREDICTOR=2'])
        o.SetGeoTransform(d.GetGeoTransform()); o.SetProjection(d.GetProjection())
        for i in range(nb):
            b=d.GetRasterBand(i+1); a=b.ReadAsArray().astype(np.float64); nod=b.GetNoDataValue()
            bad=~np.isfinite(a)
            if nod is not None and np.isfinite(nod): bad|=(a==nod)
            if scale: q=np.round(a/scale)
            else: q=np.round(a)
            lo,hi=(-32767,32767) if dtype=='int16' else (0,254)
            q=np.clip(np.where(bad,0,q),lo,hi); q[bad]=nd
            ob=o.GetRasterBand(i+1); ob.WriteArray(q.astype(np.int16 if dtype=='int16' else np.uint8)); ob.SetNoDataValue(nd)
            if scale: ob.SetScale(scale); ob.SetOffset(0)
            if b.GetDescription(): ob.SetDescription(b.GetDescription())
        o=None
    ST['n']+=1; ST['mb']=round(ST['mb']+os.path.getsize(R+dst)/1e6,1); ST['cur']=dst; sv()
def warp_sin(src,dst,res,months,unit_desc):
    tmp=gdal.Warp('',src,format='MEM',dstSRS='EPSG:32651',xRes=res,yRes=res,resampleAlg='near',srcNodata=np.nan,dstNodata=np.nan,targetAlignedPixels=True)
    for i,m in enumerate(months): tmp.GetRasterBand(i+1).SetDescription(f'{m} {unit_desc}')
    os.makedirs(os.path.dirname(R+dst),exist_ok=True)
    gdal.Translate(R+dst,tmp,creationOptions=CO+['PREDICTOR=3']); ST['n']+=1; ST['mb']=round(ST['mb']+os.path.getsize(R+dst)/1e6,1); sv()
try:
    if os.path.exists(R): ST['note']='release existed'
    NP=W+'nisar/products/'; ML=W+'ml/'; O2=W+'out_v2/'; OL=W+'outlook/'
    # NISAR backscatter
    for f in sorted(glob.glob(NP+'NISAR_H?_dB_2026????_30m.tif')): put(f,'nisar/backscatter_by_date/'+os.path.basename(f),0.01)
    for f in sorted(glob.glob(NP+'NISAR_H?_dB_2026-??_month_30m.tif')):
        if '2026-06' in f: continue
        put(f,'nisar/backscatter_by_month/'+os.path.basename(f),0.01)
    for f in sorted(glob.glob(NP+'NISAR_*_Change_*_30m.tif')): put(f,'nisar/change_same_track/'+os.path.basename(f),0.01)
    # NISAR machine learning
    for f in sorted(glob.glob(ML+'pred/NISAR_estimated_*_2026-??_month_30m.tif')):
        if '2026-06' in f: continue
        put(f,'nisar_ml/vegetation_by_month/'+os.path.basename(f),0.0001)
    for f in sorted(glob.glob(ML+'pred/NISAR_estimated_*_2026????_30m.tif')): put(f,'nisar_ml/vegetation_by_date/'+os.path.basename(f),0.0001)
    put(ML+'NISAR_LandCover_2026_30m.tif','nisar_ml/landcover/NISAR_LandCover_2026_30m.tif',dtype='byte')
    put(ML+'NISAR_LandCover_confidence_pct_30m.tif','nisar_ml/landcover/NISAR_LandCover_confidence_pct_30m.tif',dtype='byte')
    put(ML+'NISAR_LandCover_Change_vs_2020_30m.tif','nisar_ml/landcover/NISAR_LandCover_Change_vs_2020_30m.tif',dtype='byte')
    put(ML+'landcover2020_8class_30m.tif','nisar_ml/landcover/LandCover2020_8class_30m.tif',dtype='byte')
    # Landsat HLS
    for f in sorted(glob.glob(O2+'hls/*.tif')):
        b=os.path.basename(f); put(f,'landsat_hls/'+b,None if 'ValidObs' in b else 0.0001,'byte' if 'ValidObs' in b else 'int16')
    # MODIS monthly
    for f in sorted(glob.glob(O2+'modis_ndvi/*.tif')): put(f,'modis/ndvi_monthly/'+os.path.basename(f),0.0001)
    for f in sorted(glob.glob(O2+'modis_lst/*.tif')): put(f,'modis/land_temperature_monthly/'+os.path.basename(f),0.01)
    # earthquake
    for f in sorted(glob.glob(O2+'earthquake/*.tif')):
        b=os.path.basename(f); put(f,'earthquake/'+b,None if 'Landslide' in b else 0.0001,'byte' if 'Landslide' in b else 'int16')
    os.makedirs(R+'earthquake',exist_ok=True)
    for f in glob.glob(O2+'earthquake/*.gpkg'): shutil.copy(f,R+'earthquake/')
    shutil.copy(W+'Possible_landslides_after_8Jun2026_ranked.csv',R+'earthquake/')
    # terrain
    put(ML+'srtm_dem_30m.tif','terrain/SRTM_elevation_m_30m.tif',0.1); put(ML+'srtm_slope_30m.tif','terrain/SRTM_slope_deg_30m.tif',0.01); put(ML+'srtm_aspect_30m.tif','terrain/SRTM_aspect_deg_30m.tif',0.01)
    # outlook (MODIS sinusoidal -> UTM 51N)
    S=json.load(open(OL+'outlook_summary.json')); months=S['evi']['months']
    for v,res in (('evi',250),('ndvi',250),('lst',1000)):
        u='degC' if v=='lst' else v.upper()
        warp_sin(OL+f'Outlook_{v}_anom_Nov2026_Oct2027.tif',f'outlook/Outlook_{v.upper()}_anomaly_Nov2026_Oct2027_{res}m.tif',res,months,f'expected {u} anomaly vs normal')
        warp_sin(OL+f'Outlook_{v}_value_Nov2026_Oct2027.tif',f'outlook/Outlook_{v.upper()}_expected_Nov2026_Oct2027_{res}m.tif',res,months,f'expected {u}')
        warp_sin(OL+f'Outlook_{v}_pbelow_Nov2026_Oct2027.tif',f'outlook/Outlook_{v.upper()}_chance_below_normal_Nov2026_Oct2027_{res}m.tif',res,months,'probability below normal')
        warp_sin(OL+f'Outlook_{v}_ONIsensitivity_perdegC.tif',f'outlook/Outlook_{v.upper()}_ElNino_sensitivity_per_degC_{res}m.tif',res,months,'anomaly per +1 C ONI')
    for f in ('Outlook_Nov2026_Oct2027_ElNino_scenario.csv','MODIS_AOI_monthly_2000_2026.csv'): shutil.copy(OL+f,R+'outlook/')
    # boundary + tables
    os.makedirs(R+'boundary',exist_ok=True)
    for f in ('AOI_DavaoOccidental_Sarangani_EPSG32651.gpkg','AOI_DavaoOccidental_Sarangani_EPSG4326.gpkg'): shutil.copy(W+'boundary/'+f,R+'boundary/')
    os.makedirs(R+'tables',exist_ok=True); shutil.copy(O2+'AOI_monthly_means.csv',R+'tables/AOI_monthly_means_2025_2026.csv')
    # scripts
    os.makedirs(R+'scripts',exist_ok=True)
    for f in ['ts_build_v2.py','ts_build2_v2.py','landslides_v2.py','nisar/ml30_build.py','nisar/products_build.py','nisar/monthly_build.py','ml/ml1_nisar_to_vegetation.py','ml/ml1_predict_oct.py','ml/ml2_nisar_landcover_change.py','outlook/modis_series.py','outlook/outlook_pixel.py','web_export.py','web_export_nisar2.py','outlook/web_export_outlook.py','release_build.py']:
        try: shutil.copy(W+f,R+'scripts/'+os.path.basename(f))
        except Exception as e: ST['skip'].append(f)
    tot=0
    for dp,dn,fn in os.walk(R):
        for f in fn: tot+=os.path.getsize(os.path.join(dp,f))
    ST['total_mb']=round(tot/1e6,1); ST['fin']=True; sv()
except Exception:
    ST['err']=traceback.format_exc()[-2000:]; sv()

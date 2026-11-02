
import numpy as np, json, os, time, traceback
from osgeo import gdal, ogr, osr
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; O2=W+'out_v2/'; ML=W+'ml/'
LOG=W+'_ls2.json'; ST={'step':'wait','err':None,'fin':False,'res':{}}
def sv(): json.dump(ST,open(LOG,'w'),indent=1)
sv()
try:
    while True:
        try:
            if json.load(open(W+'_ml1.json'))['fin']: break
        except Exception: pass
        time.sleep(10)
    ST['step']='landslides'; sv()
    g=json.load(open(W+'_ts.json'))['notes']['grid30']; GT=(g['x0'],30,0,g['y1'],0,-30); Wd,Hd=g['W'],g['H']
    M=np.load(W+'tmp/aoi_mask30.npy'); rd=lambda p: gdal.Open(p).ReadAsArray().astype(np.float32)
    pre=rd(O2+'earthquake/HLS_EVI_pre2026_30m.tif'); post=rd(O2+'earthquake/HLS_EVI_post2026_30m.tif'); ex=rd(O2+'earthquake/HLS_EVI_Change_2026_minus_2025_seasonal_30m.tif')
    slope=rd(ML+'srtm_slope_30m.tif'); dem=rd(ML+'srtm_dem_30m.tif'); prj=gdal.Open(ML+'srtm_dem_30m.tif').GetProjection()
    hasdem=(dem>-1000)&(slope>=0); obs=np.isfinite(ex)&np.isfinite(post)
    cls=np.zeros(M.shape,np.uint8); cls[M]=4; cls[M&obs&~hasdem]=5; cls[M&obs&hasdem]=1
    c2=M&obs&hasdem&(slope>=20)&(pre>=0.35)&(post<=0.20)&(ex<=-0.25); cls[c2]=2
    mem=gdal.GetDriverByName('MEM').Create('',Wd,Hd,1,gdal.GDT_Byte); mem.SetGeoTransform(GT); mem.SetProjection(prj); b=mem.GetRasterBand(1); b.WriteArray((cls==2).astype(np.uint8)); b.SetNoDataValue(0)
    gdal.SieveFilter(b,None,b,3,8); keep=b.ReadAsArray()==1; cls[(cls==2)&~keep]=1
    os.makedirs(O2+'earthquake',exist_ok=True)
    o=gdal.GetDriverByName('GTiff').Create(O2+'earthquake/Landslide_Candidates_after_8Jun2026_30m.tif',Wd,Hd,1,gdal.GDT_Byte,options=['COMPRESS=DEFLATE','TILED=YES']); o.SetGeoTransform(GT); o.SetProjection(prj); o.GetRasterBand(1).SetNoDataValue(0); o.GetRasterBand(1).WriteArray(cls); o=None
    ST['res']['classes']={str(k):int((cls==k).sum()) for k in (1,2,4,5)}; ST['res']['candidate_ha']=round(float((cls==2).sum()*0.09),1)
    base=M&obs&hasdem&(pre>=0.35); hit=cls==2
    ST['res']['rate_by_slope_per10k']={f'{a}-{b_}':round(1e4*float((hit&base&(slope>=a)&(slope<b_)).sum())/max(int((base&(slope>=a)&(slope<b_)).sum()),1),2) for a,b_ in ((0,5),(5,15),(15,25),(25,35),(35,90))}
    ST['res']['checked_pct_of_aoi']=round(100*float((M&obs&hasdem).mean()/M.mean()),1); sv()
    # polygons
    vp=O2+'earthquake/Landslide_Candidates_after_8Jun2026_polygons.gpkg'
    if os.path.exists(vp): os.remove(vp)
    m2=gdal.GetDriverByName('MEM').Create('',Wd,Hd,1,gdal.GDT_Byte); m2.SetGeoTransform(GT); m2.SetProjection(prj); m2.GetRasterBand(1).WriteArray((cls==2).astype(np.uint8)); m2.GetRasterBand(1).SetNoDataValue(0)
    vd=ogr.GetDriverByName('GPKG').CreateDataSource(vp); lyr=vd.CreateLayer('candidates',osr.SpatialReference(wkt=prj),ogr.wkbPolygon); lyr.CreateField(ogr.FieldDefn('val',ogr.OFTInteger))
    gdal.Polygonize(m2.GetRasterBand(1),m2.GetRasterBand(1),lyr,0); vd=None
    import processing
    l=processing.run('native:fieldcalculator',{'INPUT':vp,'FIELD_NAME':'area_ha','FIELD_TYPE':0,'FIELD_LENGTH':10,'FIELD_PRECISION':2,'FORMULA':'$area/10000','OUTPUT':'memory:'})['OUTPUT']
    l=processing.run('native:fieldcalculator',{'INPUT':l,'FIELD_NAME':'dist_epi_km','FIELD_TYPE':0,'FIELD_LENGTH':10,'FIELD_PRECISION':1,'FORMULA':"distance(centroid($geometry), transform(make_point(125.0561,5.5989),'EPSG:4326','EPSG:32651'))/1000",'OUTPUT':'memory:'})['OUTPUT']
    for r,pfx in ((ML+'srtm_slope_30m.tif','slope_'),(ML+'srtm_dem_30m.tif','elev_'),(O2+'earthquake/HLS_EVI_pre2026_30m.tif','evipre_'),(O2+'earthquake/HLS_EVI_post2026_30m.tif','evipost_')):
        l=processing.run('native:zonalstatisticsfb',{'INPUT':l,'INPUT_RASTER':r,'RASTER_BAND':1,'COLUMN_PREFIX':pfx,'STATISTICS':[2],'OUTPUT':'memory:'})['OUTPUT']
    l=processing.run('native:retainfields',{'INPUT':l,'FIELDS':['area_ha','dist_epi_km','slope_mean','elev_mean','evipre_mean','evipost_mean'],'OUTPUT':'memory:'})['OUTPUT']
    l=processing.run('native:joinattributesbylocation',{'INPUT':l,'JOIN':W+'boundary/AOI_DavaoOccidental_Sarangani_EPSG32651.gpkg','PREDICATE':[0],'JOIN_FIELDS':['province'],'METHOD':2,'OUTPUT':'memory:'})['OUTPUT']
    processing.run('native:savefeatures',{'INPUT':l,'OUTPUT':O2+'earthquake/Landslide_Candidates_v2.gpkg'})
    fs=list(l.getFeatures()); from collections import Counter
    ST['res']['polygons']=len(fs); ST['res']['by_province']=dict(Counter(f['province'] for f in fs)); ST['res']['poly_ha']=round(sum(f['area_ha'] for f in fs),1)
    d=np.array([f['dist_epi_km'] for f in fs]) if fs else np.array([0])
    ST['res']['dist_km_quartiles']=[round(float(x)) for x in np.percentile(d,[0,25,50,75,100])]
    gdal.VectorTranslate(W+'webmap/data/landslide_candidates.geojson',O2+'earthquake/Landslide_Candidates_v2.gpkg',format='GeoJSON',dstSRS='EPSG:4326',layerCreationOptions=['COORDINATE_PRECISION=6','RFC7946=YES'])
    ST['step']='done'; ST['fin']=True; sv()
except Exception:
    ST['err']=traceback.format_exc(); sv()

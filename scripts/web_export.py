
import numpy as np, os, glob, re, json, calendar, traceback
from osgeo import gdal, ogr, osr
from PIL import Image
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; OUT=W+'out/'; WEB=W+'webmap/'; DATA=WEB+'data/'
os.makedirs(DATA+'img',exist_ok=True)
LOG=W+'_web.json'; ST={'step':'start','n':0,'err':None,'fin':False}
def sv(): json.dump(ST,open(LOG,'w'))
RAMP={'ndvi':['#f7fcb9','#d9f0a3','#addd8e','#78c679','#31a354','#006837'],'lst':['#ffffcc','#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026'],
      'div':['#8c510a','#d8b365','#f6e8c3','#f5f5f5','#c7eae5','#5ab4ac','#01665e'],'tdiv':['#2166ac','#67a9cf','#d1e5f0','#f7f7f7','#fddbc7','#ef8a62','#b2182b'],
      'ndwi':['#f7fbff','#c6dbef','#6baed6','#2171b5','#08306b'],'elev':['#1a9850','#a6d96a','#ffffbf','#fdae61','#a0522d','#ffffff'],'slope':['#ffffcc','#fed976','#fd8d3c','#e31a1c','#800026']}
def hex2rgb(h): return [int(h[i:i+2],16) for i in (1,3,5)]
s4326=osr.SpatialReference(); s4326.ImportFromEPSG(4326); s4326.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
s3857=osr.SpatialReference(); s3857.ImportFromEPSG(3857); s3857.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
ct=osr.CoordinateTransformation(s3857,s4326)
# common 3857 extent = AOI bounds + 3 km
v=ogr.Open(W+'boundary/AOI_DavaoOccidental_Sarangani_EPSG4326.gpkg'); e=v.GetLayer(0).GetExtent(); v=None
ct2=osr.CoordinateTransformation(s4326,s3857)
ax0,ay0,_=ct2.TransformPoint(e[0],e[2]); ax1,ay1,_=ct2.TransformPoint(e[1],e[3]); ax0-=3000; ay0-=3000; ax1+=3000; ay1+=3000
def corners():
    tl=ct.TransformPoint(ax0,ay1)[:2]; tr=ct.TransformPoint(ax1,ay1)[:2]; br=ct.TransformPoint(ax1,ay0)[:2]; bl=ct.TransformPoint(ax0,ay0)[:2]
    return [[round(p[0],6),round(p[1],6)] for p in (tl,tr,br,bl)]
COR=corners()
manifest={'bounds':[[e[0],e[2]],[e[1],e[3]]],'image_corners':COR,'layers':[]}
def colorize(a,cols,lo,hi):
    t=np.clip((a-lo)/(hi-lo),0,1); k=len(cols)-1; pos=t*k; i=np.clip(np.floor(pos).astype(int),0,k-1); f=pos-i
    C=np.array([hex2rgb(c) for c in cols],np.float32)
    rgb=C[i]*(1-f[...,None])+C[i+1]*f[...,None]
    out=np.zeros(a.shape+(4,),np.uint8); ok=np.isfinite(a); out[...,:3]=rgb.astype(np.uint8); out[...,3]=np.where(ok,235,0); return out
def export(src,lid,family,label,res,cols=None,lo=None,hi=None,date0=None,date1=None,unit='',classes=None,group=None,note=None):
    w=int(round((ax1-ax0)/res)); h=int(round((ay1-ay0)/res))
    d=gdal.Warp('',src,format='MEM',dstSRS='EPSG:3857',outputBounds=(ax0,ay0,ax1,ay1),width=w,height=h,resampleAlg='near' if classes else 'average',dstNodata=np.nan if not classes else 0,outputType=gdal.GDT_Float32 if not classes else gdal.GDT_Byte)
    a=d.GetRasterBand(1).ReadAsArray(); d=None
    if classes:
        out=np.zeros(a.shape+(4,),np.uint8)
        for val,col,txt in classes:
            m=a==val; out[m,:3]=hex2rgb(col); out[m,3]=235
        leg={'type':'classes','items':[[col,txt] for val,col,txt in classes]}
    else:
        a=a.astype(np.float32); out=colorize(a,cols,lo,hi); leg={'type':'ramp','colors':cols,'min':lo,'max':hi,'unit':unit}
    fn=f'img/{lid}.png'; Image.fromarray(out,'RGBA').save(DATA+fn,optimize=True)
    item=dict(id=lid,family=family,group=group or family,label=label,file='data/'+fn,legend=leg)
    if date0: item['start']=date0; item['end']=date1
    if note: item['note']=note
    manifest['layers'].append(item); ST['n']+=1; ST['step']=lid; sv()
def mrange(k):
    y,m=map(int,k.split('-')); return f'{y}-{m:02d}-01', f'{y}-{m:02d}-{calendar.monthrange(y,m)[1]:02d}'
try:
    for f in sorted(glob.glob(OUT+'modis_ndvi/MODIS_NDVI_20??-??_250m.tif')):
        k=re.search(r'(\d{4}-\d{2})',os.path.basename(f)).group(1); export(f,'ndvi_'+k,'ndvi','NDVI '+k,250,RAMP['ndvi'],0.3,0.95,*mrange(k))
    for f in sorted(glob.glob(OUT+'modis_lst/MODIS_LST_Day_20??-??_1km.tif')):
        k=re.search(r'(\d{4}-\d{2})',os.path.basename(f)).group(1); export(f,'lst_'+k,'lst','Day land temperature '+k,1000,RAMP['lst'],20,36,*mrange(k),unit=' °C')
    for f in sorted(glob.glob(OUT+'modis_ndvi/MODIS_NDVI_Anomaly_*.tif')):
        k=re.search(r'Anomaly_(\d{4}-\d{2})',f).group(1); export(f,'ndvianom_'+k,'ndvi_anom','NDVI '+k+' vs same month 2025',250,RAMP['div'],-0.15,0.15,*mrange(k))
    for f in sorted(glob.glob(OUT+'modis_lst/MODIS_LST_Anomaly_*.tif')):
        k=re.search(r'Anomaly_(\d{4}-\d{2})',f).group(1); export(f,'lstanom_'+k,'lst_anom','Land temperature '+k+' vs same month 2025',1000,RAMP['tdiv'],-4,4,*mrange(k),unit=' °C')
    QR={'Q1':('01-01','03-31'),'Q2':('04-01','06-30'),'Q3':('07-01','09-30'),'Q4':('10-01','12-31')}
    for idx,ramp,lo,hi in (('EVI','ndvi',0,0.8),('NDMI','div',-0.2,0.6),('NDWI','ndwi',-0.8,0.3)):
        for f in sorted(glob.glob(OUT+f'hls/HLS_{idx}_20??Q?_30m.tif')):
            m=re.search(r'(\d{4})(Q\d)',f); y=m.group(1); q=m.group(2)
            export(f,f'{idx.lower()}_{y}{q}',idx.lower(),f'{idx} {y} {q}',60,RAMP[ramp],lo,hi,f'{y}-{QR[q][0]}',f'{y}-{QR[q][1]}')
    lsc=[(2,'#e31a1c','Landslide / vegetation-loss candidate'),(1,'#d9d9d9','Checked, no candidate'),(4,'#9ecae1','No usable optical pair (cloud or outside data)'),(5,'#fdd0a2','No elevation data (south of 5.9°N)')]
    export(OUT+'earthquake/Landslide_Candidates_after_8Jun2026_30m.tif','eq_classes','static','Landslide check coverage',30,classes=lsc,group='earthquake')
    export(OUT+'earthquake/HLS_EVI_Change_2026_minus_2025_seasonal_30m.tif','eq_evi_excess','static','Vegetation change after the quake, beyond normal seasonal change (EVI)',30,RAMP['div'],-0.3,0.3,group='earthquake')
    export(OUT+'earthquake/HLS_EVI_pre2026_30m.tif','eq_evi_pre','static','Vegetation before the quake (EVI, 1 Mar–7 Jun 2026)',60,RAMP['ndvi'],0,0.8,group='earthquake')
    export(OUT+'earthquake/HLS_EVI_post2026_30m.tif','eq_evi_post','static','Vegetation after the quake (EVI, 9 Jun–25 Sep 2026)',60,RAMP['ndvi'],0,0.8,group='earthquake')
    export(OUT+'terrain/Slope_AW3D30_30m.tif','slope','static','Slope (ALOS AW3D30)',60,RAMP['slope'],0,45,unit='°',group='terrain')
    # vectors -> GeoJSON
    gdal.VectorTranslate(DATA+'boundary.geojson',W+'boundary/AOI_DavaoOccidental_Sarangani_EPSG4326.gpkg',format='GeoJSON',simplifyTolerance=0.0004,layerCreationOptions=['COORDINATE_PRECISION=5','RFC7946=YES'])
    gdal.VectorTranslate(DATA+'landslide_candidates.geojson',OUT+'earthquake/Landslide_Candidates_after_8Jun2026_polygons.gpkg',format='GeoJSON',dstSRS='EPSG:4326',layerCreationOptions=['COORDINATE_PRECISION=6','RFC7946=YES'])
    import csv
    manifest['aoi_monthly']=list(csv.DictReader(open(OUT+'AOI_monthly_means.csv')))
    manifest['epicentre']={'lng':125.14699,'lat':6.02153,'label':'M7.8 epicentre, 8 Jun 2026, 7:37 am, 33 km deep'}
    json.dump(manifest,open(DATA+'manifest.json','w'),indent=1)
    tot=sum(os.path.getsize(os.path.join(DATA+'img',f)) for f in os.listdir(DATA+'img'))
    ST['step']='done'; ST['fin']=True; ST['MB']=round(tot/1e6,1); sv()
except Exception:
    ST['err']=traceback.format_exc(); sv()

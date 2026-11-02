
import numpy as np, os, glob, re, json, calendar, traceback
from osgeo import gdal, ogr, osr
from PIL import Image
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; OUT=W+'out/'; WEB=W+'webmap/'; DATA=WEB+'data/'
os.makedirs(DATA+'img',exist_ok=True)
LOG=W+'_web_nisar2.json'; ST={'step':'start','n':0,'err':None,'fin':False}
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
manifest_old=json.load(open(DATA+'manifest.json'))
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
    fn=f'img/{lid}.webp'
    m_=gdal.GetDriverByName('MEM').Create('',out.shape[1],out.shape[0],4,gdal.GDT_Byte)
    for i_ in range(4): m_.GetRasterBand(i_+1).WriteArray(out[...,i_])
    gdal.GetDriverByName('WEBP').CreateCopy(DATA+fn,m_,options=['QUALITY=88']); m_=None
    item=dict(id=lid,family=family,group=group or family,label=label,file='data/'+fn,legend=leg)
    if date0: item['start']=date0; item['end']=date1
    if note: item['note']=note
    manifest['layers'].append(item); ST['n']+=1; ST['step']=lid; sv()
def mrange(k):
    y,m=map(int,k.split('-')); return f'{y}-{m:02d}-01', f'{y}-{m:02d}-{calendar.monthrange(y,m)[1]:02d}'

RAMP['hv']=['#33505f','#3f7f73','#73ad6c','#b8d58c','#f1f0cf']; RAMP['hh']=['#47507a','#5f7fa3','#8fb1bd','#d3dfc4','#fbf5e3']
import datetime as dt
NO=W+'nisar/products/'
try:
    for pol,ramp,lo,hi in (('HV','hv',-22,-8),('HH','hh',-14,-2)):
        for f in sorted(f for f in glob.glob(NO+f'NISAR_{pol}_dB_2026????_30m.tif') if '20261019' not in f):
            d=re.search(r'(\d{8})',os.path.basename(f)).group(1); iso=f'{d[:4]}-{d[4:6]}-{d[6:]}'
            export(f,f'nisar_{pol.lower()}_{d}','nisar_'+pol.lower(),f'NISAR {pol} backscatter {iso}',60,RAMP[ramp],lo,hi,iso,iso,unit=' dB')
        export(NO+f'NISAR_{pol}_dB_SepOct2026_mean_30m.tif',f'nisar_{pol.lower()}_sepoct','nisar_'+pol.lower()+'_combined',f'NISAR {pol} backscatter, all passes 20 Sep – 17 Oct 2026 combined',60,RAMP[ramp],lo,hi,'2026-09-20','2026-10-17',unit=' dB')
    export(NO+'NISAR_passes_SepOct2026_count_30m.tif','nisar_passes','static','Number of NISAR passes, 20 Sep – 17 Oct 2026',60,classes=[(1,'#fdd0a2','1 pass'),(2,'#fdae6b','2 passes'),(3,'#f16913','3 passes'),(4,'#d94801','4 passes'),(5,'#8c2d04','5 passes')],group='nisar_cov')
    for t_,ramp_,lo_,hi_,lab_ in (('EVI','ndvi',0,0.8,'vegetation (EVI)'),('NDMI','div',-0.2,0.6,'plant moisture (NDMI)')):
        export(W+f'ml/pred/NISAR_estimated_{t_}_SepOct2026_mean_30m.tif',f'ml_{t_.lower()}_sepoct','ml_'+t_.lower()+'_combined',f'NISAR-estimated {lab_}, Sep–Oct 2026 average',60,RAMP[ramp_],lo_,hi_,'2026-09-20','2026-10-17')
    for f in sorted(glob.glob(NO+'NISAR_*_Change_*_30m.tif')):
        m=re.search(r'NISAR_(H.)_Change_T(\d{3})([AD])(\d{3})_(\d{8})_to_(\d{8})',os.path.basename(f))
        pol,tr,ad,fr,a,b=m.groups(); fa=dt.date(int(a[:4]),int(a[4:6]),int(a[6:])).strftime('%d %b'); fb=dt.date(int(b[:4]),int(b[4:6]),int(b[6:])).strftime('%d %b')
        export(f,f'nisar_chg_{pol.lower()}_{tr}{ad}{fr}_{a}',('static'),f'NISAR {pol} change, track {int(tr)} {"ascending" if ad=="A" else "descending"} frame {int(fr)}: {fa} → {fb}',60,RAMP['div'],-4,4,unit=' dB',group='nisar')
    old=[l for l in manifest_old['layers'] if not (l['id'].startswith('nisar_') or l['id'].endswith('_sepoct'))]
    for l in manifest['layers']: l['file']=l['file']
    manifest_old['layers']=old+manifest['layers']
    json.dump(manifest_old,open(DATA+'manifest.json','w'),indent=1)
    ST['step']='done'; ST['fin']=True; sv()
except Exception:
    ST['err']=traceback.format_exc(); sv()

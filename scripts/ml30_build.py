
import os, glob, json, time, re, traceback
from osgeo import gdal
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; OUTN=W+'nisar/ml30/'; os.makedirs(OUTN,exist_ok=True)
g=json.load(open(W+'_ts.json'))['notes']['grid30']; x0,y1,Wd,Hd=g['x0'],g['y1'],g['W'],g['H']; B=(x0,y1-Hd*30,x0+Wd*30,y1)
src=[]
for d in ['C:/Users/PETER DOMINIC/Pictures/NISAR_L2_GCOV_PROVISIONAL_V1_1-20261101_001247/','C:/Users/PETER DOMINIC/Pictures/','C:/Users/PETER DOMINIC/Downloads/']:
    try: fs=os.listdir(d)
    except Exception: continue
    for f in fs:
        if f.startswith('NISAR_L2_PR_GCOV') and f.endswith('_reformatted.tif'): src.append(d+f)
seen={}
for p in src:
    m=re.search(r'GCOV_(\d+)_(\d+)_([AD])_(\d+)_(\d+)_(\w{4})_\w_(\d{8})T.*?(HHHH|HVHV)_reformatted',os.path.basename(p))
    key=f'{m.group(7)}_T{m.group(2)}{m.group(3)}{m.group(4)}_{m.group(6)}_{m.group(8)}'
    if key not in seen: seen[key]=p
LOG=W+'_nisar_ml.json'; st={'todo':sorted(seen),'done':{},'skip':{},'cur':None,'fin':False}
def sv(): json.dump(st,open(LOG,'w'),indent=1)
sv()
for key in sorted(seen):
    dst=OUTN+key+'.tif'
    if os.path.exists(dst): st['done'][key]='exists'; continue
    st['cur']=key; sv(); t=time.time()
    try:
        gdal.Warp(dst+'.part',seen[key],format='GTiff',outputBounds=B,width=Wd,height=Hd,dstSRS='EPSG:32651',resampleAlg='average',srcNodata=float('nan'),dstNodata=float('nan'),
                  outputType=gdal.GDT_Float32,overviewLevel='NONE',multithread=True,warpMemoryLimit=1024,creationOptions=['COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES'])
        os.replace(dst+'.part',dst)
        d=gdal.Open(dst); b=d.GetRasterBand(1); a=b.ReadAsArray(buf_xsize=Wd//10,buf_ysize=Hd//10); import numpy as np
        v=np.isfinite(a)&(a>0); st['done'][key]=dict(sec=round(time.time()-t),valid_pct_of_grid=round(100*float(v.mean()),1),src=seen[key][:40]); d=None
        if not v.any(): os.remove(dst); st['skip'][key]='no overlap'
    except Exception as e:
        st['skip'][key]=repr(e)[:200]
    sv()
st['cur']=None; st['fin']=True; sv()

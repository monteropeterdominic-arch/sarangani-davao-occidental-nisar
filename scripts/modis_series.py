import os, re, json, time, datetime as dt, traceback
import numpy as np
from osgeo import gdal
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; D='C:/Users/PETER DOMINIC/Downloads/'; O=W+'outlook/'
LOG=W+'_outlook_modis.json'; st={'stage':'start','n':0,'fin':False}
def sv(): json.dump(st,open(LOG,'w'),indent=1)
try:
    fs=[f for f in os.listdir(D) if f.endswith('.tif') and '(1)' not in f]
    def dates(pfx):
        r={}
        for f in fs:
            if f.startswith(pfx+'_doy'):
                m=re.search(r'_doy(\d{7})',f); r[m.group(1)]=D+f
        return r
    ndvi=dates('MOD13Q1.061__250m_16_days_NDVI'); evi=dates('MOD13Q1.061__250m_16_days_EVI'); rel=dates('MOD13Q1.061__250m_16_days_pixel_reliability')
    lst=dates('MOD11A2.061_LST_Day_1km'); qc=dates('MOD11A2.061_QC_Day')
    ref=gdal.Open(ndvi[sorted(ndvi)[0]]); gt=ref.GetGeoTransform(); prj=ref.GetProjection(); X,Y=ref.RasterXSize,ref.RasterYSize
    refL=gdal.Open(lst[sorted(lst)[0]]); gtL=refL.GetGeoTransform(); XL,YL=refL.RasterXSize,refL.RasterYSize
    def mask(gt,X,Y,attr):
        m=gdal.GetDriverByName('MEM').Create('',X,Y,1,gdal.GDT_Byte); m.SetGeoTransform(gt); m.SetProjection(prj)
        for v,nm in [(1,'Davao Occidental'),(2,'Sarangani')]:
            gdal.Rasterize(m,W+'boundary/AOI_DavaoOccidental_Sarangani_EPSG4326.gpkg',burnValues=[v],where=f"province='{nm}'")
        return m.GetRasterBand(1).ReadAsArray()
    prov=mask(gt,X,Y,'fid'); provL=mask(gtL,XL,YL,'fid')
    # land cover 2020 (8 class) mode onto MODIS grid
    lcd=gdal.Warp('',W+'ml/landcover2020_8class_30m.tif',format='MEM',dstSRS=prj,outputBounds=(gt[0],gt[3]+Y*gt[5],gt[0]+X*gt[1],gt[3]),width=X,height=Y,resampleAlg='mode')
    lc=lcd.GetRasterBand(1).ReadAsArray(); lcd=None
    np.save(O+'modis250_prov.npy',prov); np.save(O+'modis1km_prov.npy',provL); np.save(O+'modis250_lc2020.npy',lc)
    json.dump(dict(gt=gt,prj=prj,X=X,Y=Y,gtL=gtL,XL=XL,YL=YL),open(O+'modis_grids.json','w'))
    def mid(doy,off): d=dt.date(int(doy[:4]),1,1)+dt.timedelta(int(doy[4:])-1+off); return d.year*100+d.month
    months=sorted({mid(k,8) for k in ndvi}); mi={m:i for i,m in enumerate(months)}; st['months']=[months[0],months[-1],len(months)]; sv()
    NV=np.lib.format.open_memmap(O+'modis_ndvi_month.npy','w+',np.int16,(len(months),Y,X)); NV[:]=-3000
    EV=np.lib.format.open_memmap(O+'modis_evi_month.npy','w+',np.int16,(len(months),Y,X)); EV[:]=-3000
    st['stage']='vi'; sv()
    for k in sorted(ndvi):
        if k not in rel or k not in evi: continue
        i=mi[mid(k,8)]
        r=gdal.Open(rel[k]).ReadAsArray(); good=(r==0)|(r==1)
        a=gdal.Open(ndvi[k]).ReadAsArray(); a=np.where(good,a,-3000).astype(np.int16)
        b=gdal.Open(evi[k]).ReadAsArray(); b=np.where(good,b,-3000).astype(np.int16)
        NV[i]=np.maximum(NV[i],a); EV[i]=np.maximum(EV[i],b); st['n']+=1
        if st['n']%25==0: st['cur']=k; sv()
    NV.flush(); EV.flush()
    lm=sorted({mid(k,4) for k in lst}); li={m:j for j,m in enumerate(lm)}
    S=np.zeros((len(lm),YL,XL),np.float32); C=np.zeros((len(lm),YL,XL),np.uint8)
    st['stage']='lst'; st['lst_months']=[lm[0],lm[-1],len(lm)]; sv()
    for k in sorted(lst):
        if k not in qc: continue
        a=gdal.Open(lst[k]).ReadAsArray().astype(np.float32); q=gdal.Open(qc[k]).ReadAsArray()
        ok=(a>0)&((q&3)<=1)
        j=li[mid(k,4)]; S[j]+=np.where(ok,a*0.02-273.15,0); C[j]+=ok
    LSTm=np.where(C>0,S/np.maximum(C,1),np.nan).astype(np.float32); np.save(O+'modis_lst_month.npy',LSTm)
    # AOI / province / class tables
    rows=[]
    for i,m in enumerate(months):
        n=NV[i].astype(np.float32); e=EV[i].astype(np.float32); vn=(n>-2000)&(prov>0)
        rec=dict(month=f'{m//100}-{m%100:02d}',ndvi=float(n[vn].mean()/1e4) if vn.any() else None,evi=float(e[vn].mean()/1e4) if vn.any() else None,valid_pct=round(100*vn.sum()/max(1,(prov>0).sum()),1))
        for p,nm in [(1,'davocc'),(2,'sarangani')]:
            v=vn&(prov==p); rec['ndvi_'+nm]=float(n[v].mean()/1e4) if v.sum()>200 else None
        for c in range(1,9):
            v=vn&(lc==c); rec[f'ndvi_lc{c}']=float(n[v].mean()/1e4) if v.sum()>100 else None
        if m in li:
            t=LSTm[li[m]]; v=np.isfinite(t)&(provL>0); rec['lst']=float(t[v].mean()) if v.sum()>50 else None; rec['lst_valid_pct']=round(100*v.sum()/max(1,(provL>0).sum()),1)
        rows.append(rec)
    json.dump(dict(months=months,lst_months=lm,rows=rows),open(O+'modis_aoi_monthly.json','w'))
    import csv
    keys=sorted({k for r in rows for k in r},key=lambda k:(k!='month',k))
    with open(O+'MODIS_AOI_monthly_2000_2026.csv','w',newline='') as f:
        w=csv.DictWriter(f,keys); w.writeheader(); [w.writerow(r) for r in rows]
    st['stage']='done'; st['fin']=True; sv()
except Exception:
    st['err']=traceback.format_exc()[-1500:]; sv()

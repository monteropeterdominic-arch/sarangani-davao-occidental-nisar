import json, time, traceback, numpy as np
from osgeo import gdal
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; O=W+'outlook/'
ST={'step':'wait','fin':False}
def sv(): json.dump(ST,open(W+'_web_outlook.json','w'),indent=1)
sv()
while True:
    try:
        a=json.load(open(W+'_outlook_px.json')); b=json.load(open(W+'_chain.json'))
        if (a.get('fin') or a.get('err')) and b.get('step')=='done': break
    except Exception: pass
    time.sleep(10)
try:
    S2=open(W+'web_export_nisar2.py').read(); head=S2.split("\nRAMP['hv']")[0]
    H={'__name__':'h'}; exec(compile(head.replace("LOG=W+'_web_nisar2.json'","LOG=W+'_web_outlook_inner.json'"),'h','exec'),H)
    export=H['export']; DATA=H['DATA']; manifest=H['manifest']; mrange=H['mrange']
    SUM=json.load(open(O+'outlook_summary.json')); months=SUM['evi']['months']
    import calendar
    def band(path,i,res):
        src=gdal.Translate('',path,format='MEM',bandList=[i+1])
        return gdal.Warp('',src,format='MEM',dstSRS='EPSG:3857',outputBounds=(H['ax0'],H['ay0'],H['ax1'],H['ay1']),xRes=res,yRes=res,resampleAlg='bilinear',srcNodata=np.nan,dstNodata=np.nan)
    DIV=['#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e']
    PB=['#01665e','#80cdc1','#f5f5f5','#dfc27d','#a6611a','#6b3a0c']
    TD=['#2166ac','#67a9cf','#d1e5f0','#f7f7f7','#fddbc7','#ef8a62','#b2182b']
    for i,mo in enumerate(months):
        s,e=mrange(mo); y,m=map(int,mo.split('-')); nm=calendar.month_name[m]
        ST['step']=mo; sv()
        export(band(O+'Outlook_evi_anom_Nov2026_Oct2027.tif',i,120),f'out_evi_{mo}','out_evi',f'Expected greenness change, {nm} {y} (EVI vs normal)',120,DIV,-0.12,0.12,s,e)
        export(band(O+'Outlook_evi_pbelow_Nov2026_Oct2027.tif',i,120),f'out_pb_{mo}','out_pb',f'Chance of below-normal greenness, {nm} {y}',120,PB,0,1,s,e)
        try: export(band(O+'Outlook_lst_anom_Nov2026_Oct2027.tif',i,250),f'out_lst_{mo}','out_lst',f'Expected daytime ground temperature vs normal, {nm} {y}',250,TD,-2,2,s,e,unit=' °C')
        except Exception as ex: ST['lst_err']=repr(ex)
    # dry-season sensitivity: mean ONI coefficient Feb-May (leads 3..6)
    d=gdal.Open(O+'Outlook_evi_ONIsensitivity_perdegC.tif'); arr=np.stack([d.GetRasterBand(k+1).ReadAsArray() for k in (3,4,5,6)])
    with np.errstate(all='ignore'): mb=np.nanmean(arr,0).astype(np.float32)
    t=gdal.GetDriverByName('MEM').Create('',d.RasterXSize,d.RasterYSize,1,gdal.GDT_Float32); t.SetGeoTransform(d.GetGeoTransform()); t.SetProjection(d.GetProjection()); t.GetRasterBand(1).WriteArray(mb); t.GetRasterBand(1).SetNoDataValue(np.nan)
    w=gdal.Warp('',t,format='MEM',dstSRS='EPSG:3857',outputBounds=(H['ax0'],H['ay0'],H['ax1'],H['ay1']),xRes=120,yRes=120,resampleAlg='bilinear',srcNodata=np.nan,dstNodata=np.nan)
    export(w,'out_sens','out_static','Greenness drop per +1 °C of El Niño, Feb–May (2000–2026)',120,DIV,-0.06,0.06)
    prev=json.load(open(DATA+'manifest.json'))
    prev['layers']=[l for l in prev['layers'] if not l['id'].startswith('out_')]+manifest['layers']
    prev['outlook']=json.loads(json.dumps(SUM).replace('NaN','null'))
    json.dump(prev,open(DATA+'manifest.json','w'),indent=1,allow_nan=False)
    ST['step']='done'; ST['fin']=True; sv()
except Exception:
    ST['err']=traceback.format_exc()[-2000:]; sv()

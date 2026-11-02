
import numpy as np, os, glob, re, json, time, traceback, datetime as dt, gc
from osgeo import gdal
from scipy import ndimage as ndi
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; D='C:/Users/PETER DOMINIC/Downloads/'; NP=W+'nisar/products/'; ML=W+'ml/'; os.makedirs(ML+'pred',exist_ok=True)
LOG=W+'_ml1.json'; ST={'step':'start','err':None,'fin':False,'res':{}}
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
    ASC={'20260628':1,'20260920':1,'20260925':1,'20261014':1,'20260923':0,'20260930':0,'20261017':0}
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
    def hls_nearest(idx,d):
        t=dt.date(int(d[:4]),int(d[4:6]),int(d[6:])); best=np.full((Hd,Wd),np.nan,np.float32); bdt=np.full((Hd,Wd),99,np.int8); n=0
        for f in glob.glob(D+f'HLSL30_VI.020_{idx}_doy2026*_51N.tif'):
            j=int(re.search(r'doy2026(\d{3})',f).group(1)); dd=abs(((dt.date(2026,1,1)+dt.timedelta(j-1))-t).days)
            if dd>8: continue
            x=gdal.Warp('',f,format='MEM',dstSRS='EPSG:32651',outputBounds=B,width=Wd,height=Hd,resampleAlg='near',srcNodata=-19999,dstNodata=-19999)
            a=x.GetRasterBand(1).ReadAsArray(); x=None
            ok=(a!=-19999)&(a>=-10000)&(a<=10000)&(dd<bdt); best[ok]=a[ok]*1e-4; bdt[ok]=dd; n+=1; del a,ok
        return (best if n else None),n
    rng=np.random.default_rng(7); BLK=5000
    rowblk=(np.arange(Hd)*30//BLK).astype(np.int32); colblk=(np.arange(Wd)*30//BLK).astype(np.int32)
    X=[];YE=[];YM=[];G=[];DDs=[]
    for d in dates:
        ST['step']=f'samples {d}'; sv()
        hh,hv=load(d); te,n1=hls_nearest('EVI',d); tm,n2=hls_nearest('NDMI',d)
        if te is None: ST['res'][f'pairs_{d}']='no HLS within 8 days'; sv(); continue
        ok=M&np.isfinite(hv)&np.isfinite(hh)&np.isfinite(te)&np.isfinite(tm)&np.isfinite(dem)
        rr,cc=np.nonzero(ok); ST['res'][f'pairs_{d}']=dict(hls_scenes=n1,matched_px=int(rr.size)); sv()
        if rr.size:
            k=rng.choice(rr.size,min(rr.size,30000),replace=False); rr=rr[k]; cc=cc[k]; o=np.argsort(rr); rr=rr[o]; cc=cc[o]
            for r0 in range(0,Hd,400):
                sel=(rr>=r0)&(rr<r0+400)
                if not sel.any(): continue
                F=feats_rows(hh,hv,d,r0,min(Hd,r0+400)); X.append(F[rr[sel]-r0,cc[sel]]); YE.append(te[rr[sel],cc[sel]]); YM.append(tm[rr[sel],cc[sel]])
                G.append(rowblk[rr[sel]]*1000+colblk[cc[sel]]); DDs.append(np.full(sel.sum(),int(d)))
        del hh,hv,te,tm,ok; gc.collect()
    X=np.concatenate(X); G=np.concatenate(G); YE=np.concatenate(YE); YM=np.concatenate(YM); DDs=np.concatenate(DDs)
    good=np.isfinite(X).all(1); X,G,YE,YM,DDs=X[good],G[good],YE[good],YM[good],DDs[good]
    np.savez_compressed(ML+'training_samples.npz',X=X,G=G,D=DDs,EVI=YE,NDMI=YM,FN=np.array(FN))
    ST['res']['n_samples']=int(X.shape[0]); ST['res']['n_blocks_5km']=int(np.unique(G).size); sv()
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
    import joblib
    def mk(): return HistGradientBoostingRegressor(max_iter=400,learning_rate=0.06,max_leaf_nodes=48,min_samples_leaf=60,l2_regularization=1.0,random_state=0)
    def cv(cols,y):
        P=np.full(y.size,np.nan)
        for tr,te_ in GroupKFold(n_splits=5).split(X,y,G):
            m=mk().fit(X[tr][:,cols],y[tr]); P[te_]=m.predict(X[te_][:,cols])
        return dict(R2=round(float(r2_score(y,P)),3),RMSE=round(float(np.sqrt(mean_squared_error(y,P))),4),MAE=round(float(mean_absolute_error(y,P)),4)),P
    models={}
    for tgt,y in (('EVI',YE),('NDMI',YM)):
        for name,cols in (('NISAR only',NIS),('Terrain only',TER),('NISAR + terrain',NIS+TER)):
            ST['step']=f'CV {tgt} {name}'; sv(); r,P=cv(cols,y); ST['res'][f'{tgt} | {name}']=r; sv()
            if name=='NISAR + terrain': np.save(ML+f'cv_pred_{tgt}.npy',P)
        models[tgt]=mk().fit(X[:,NIS+TER],y); joblib.dump(models[tgt],ML+f'model_{tgt}_nisar_terrain.joblib')
    from sklearn.inspection import permutation_importance
    tr,te_=next(GroupKFold(5).split(X,YE,G)); mm=mk().fit(X[tr][:,NIS+TER],YE[tr])
    pi=permutation_importance(mm,X[te_][:15000][:,NIS+TER],YE[te_][:15000],n_repeats=5,random_state=0)
    ST['res']['importance_EVI']=sorted([(FN[c],round(float(v),4)) for c,v in zip(NIS+TER,pi.importances_mean)],key=lambda t:-t[1]); sv()
    for d in dates:
        ST['step']=f'predict {d}'; sv()
        hh,hv=load(d)
        for tgt,m in models.items():
            P=np.full((Hd,Wd),np.nan,np.float32)
            for r0 in range(0,Hd,400):
                r1=min(Hd,r0+400); F=feats_rows(hh,hv,d,r0,r1); ok=M[r0:r1]&np.isfinite(F).all(-1)
                if ok.any(): blk=P[r0:r1]; blk[ok]=m.predict(F[ok][:,NIS+TER]); P[r0:r1]=blk
            o=gdal.GetDriverByName('GTiff').Create(ML+f'pred/NISAR_estimated_{tgt}_{d}_30m.tif',Wd,Hd,1,gdal.GDT_Float32,options=['COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES'])
            o.SetGeoTransform(GT); o.SetProjection(prj); b=o.GetRasterBand(1); b.SetNoDataValue(float('nan')); b.WriteArray(P); o=None
            if tgt=='EVI':
                te,_=hls_nearest('EVI',d); nis_ok=np.isfinite(P)
                if te is not None: ST['res'][f'gapfill_{d}']=dict(nisar_px=int(nis_ok.sum()),hls_clear_pct_of_those=round(100*float(np.isfinite(te)[nis_ok].mean()),1),cloud_gap_filled_px=int((nis_ok&~np.isfinite(te)).sum()))
                del te
            del P; gc.collect()
        sv(); del hh,hv; gc.collect()
    ST['step']='done'; ST['fin']=True; sv()
except Exception:
    ST['err']=traceback.format_exc(); sv()

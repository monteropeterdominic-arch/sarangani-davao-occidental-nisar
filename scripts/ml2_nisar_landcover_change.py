
import numpy as np, os, glob, re, json, time, traceback, gc
from osgeo import gdal, ogr
from scipy import ndimage as ndi
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; NP=W+'nisar/products/'; ML=W+'ml/'
LOG=W+'_ml2.json'; ST={'step':'start','err':None,'fin':False,'res':{}}
def sv(): json.dump(ST,open(LOG,'w'),indent=1)
try:
    g=json.load(open(W+'_ts.json'))['notes']['grid30']; x0,y1,Wd,Hd=g['x0'],g['y1'],g['W'],g['H']; GT=(x0,30,0,y1,0,-30)
    prj=gdal.Open(ML+'srtm_dem_30m.tif').GetProjection(); M=np.load(W+'tmp/aoi_mask30.npy')
    CL={1:(1,'Closed forest'),4:(2,'Open forest & mangrove'),7:(2,''),10:(3,'Brush/shrubs'),14:(4,'Grassland'),16:(5,'Annual crop'),17:(6,'Perennial crop'),20:(7,'Built-up & barren'),13:(7,''),21:(8,'Water & fishpond'),19:(8,'')}
    NAMES={v[0]:v[1] for v in CL.values() if v[1]}
    ST['step']='labels'; sv()
    mem=gdal.GetDriverByName('MEM').Create('',Wd,Hd,1,gdal.GDT_Byte); mem.SetGeoTransform(GT); mem.SetProjection(prj)
    ds=ogr.Open(W+'boundary/_lc_aoi.gpkg'); lyr=ds.GetLayer(0)
    for src,(dst,_) in CL.items():
        lyr.SetAttributeFilter(f'class_id={src}'); gdal.RasterizeLayer(mem,[1],lyr,burn_values=[dst])
    LAB=mem.ReadAsArray(); mem=None; ds=None
    o=gdal.GetDriverByName('GTiff').Create(ML+'landcover2020_8class_30m.tif',Wd,Hd,1,gdal.GDT_Byte,options=['COMPRESS=DEFLATE','TILED=YES']); o.SetGeoTransform(GT); o.SetProjection(prj); o.GetRasterBand(1).SetNoDataValue(0); o.GetRasterBand(1).WriteArray(np.where(M,LAB,0)); o=None
    # interior pixels (same label in 3x3) to avoid polygon-edge label noise
    interior=(ndi.minimum_filter(LAB,3)==ndi.maximum_filter(LAB,3))&(LAB>0)&M
    dem=gdal.Open(ML+'srtm_dem_30m.tif').ReadAsArray().astype(np.float32); dem[dem<=-1000]=np.nan
    slope=gdal.Open(ML+'srtm_slope_30m.tif').ReadAsArray().astype(np.float32); slope[slope<0]=np.nan
    dates=sorted(re.search(r'(\d{8})',os.path.basename(f)).group(1) for f in glob.glob(NP+'NISAR_HV_dB_2026????_30m.tif')); dates=[d for d in dates if d!='20261019']
    HVd=[gdal.Open(NP+f'NISAR_HV_dB_{d}_30m.tif') for d in dates]
    HHd=[gdal.Open(NP+f'NISAR_HH_dB_{d}_30m.tif') if os.path.exists(NP+f'NISAR_HH_dB_{d}_30m.tif') else None for d in dates]
    FN=['HV_mean','HV_std','HV_min','HV_max','HH_mean','HH_std','HVminusHH_mean','HVmean_std3x3','HVmean_mean7x7','n_dates','elev','slope']
    NIS=list(range(10)); ALL=list(range(12))
    def feats(r0,r1):
        a0=max(0,r0-3); a1=min(Hd,r1+3); n=a1-a0
        hv=np.stack([d.GetRasterBand(1).ReadAsArray(0,a0,Wd,n).astype(np.float32) for d in HVd])
        hh=np.stack([(d.GetRasterBand(1).ReadAsArray(0,a0,Wd,n).astype(np.float32) if d is not None else np.full((n,Wd),np.nan,np.float32)) for d in HHd])
        with np.errstate(all='ignore'):
            cnt=np.isfinite(hv).sum(0); hvm=np.nanmean(hv,0); hvs=np.nanstd(hv,0); hvmin=np.nanmin(hv,0); hvmax=np.nanmax(hv,0)
            hhm=np.nanmean(hh,0); hhs=np.nanstd(hh,0); dif=np.nanmean(hv-hh,0)
        f=np.where(np.isfinite(hvm),hvm,np.nanmedian(hvm) if np.isfinite(hvm).any() else -14)
        m3=ndi.uniform_filter(f,3); s3=np.sqrt(np.maximum(ndi.uniform_filter(f*f,3)-m3*m3,0)); m7=ndi.uniform_filter(f,7)
        sl=slice(r0-a0,r0-a0+(r1-r0))
        F=np.stack([hvm[sl],hvs[sl],hvmin[sl],hvmax[sl],hhm[sl],hhs[sl],dif[sl],s3[sl],m7[sl],cnt[sl].astype(np.float32),dem[r0:r1],slope[r0:r1]],-1)
        bad=(cnt[sl]<2)|~np.isfinite(hhm[sl])
        F[bad]=np.nan
        return F
    # ---- sample
    rng=np.random.default_rng(11); per={c:[] for c in NAMES}
    for c in NAMES:
        rr,cc=np.nonzero(interior&(LAB==c)); k=rng.choice(rr.size,min(rr.size,60000),replace=False); per[c]=(rr[k],cc[k])
    RR=np.concatenate([per[c][0] for c in NAMES]); CC=np.concatenate([per[c][1] for c in NAMES]); o=np.argsort(RR); RR,CC=RR[o],CC[o]
    X=[];Y=[];G=[]
    for r0 in range(0,Hd,400):
        ST['step']=f'features rows {r0}'; sv()
        s=(RR>=r0)&(RR<r0+400)
        if not s.any(): continue
        F=feats(r0,min(Hd,r0+400)); X.append(F[RR[s]-r0,CC[s]]); Y.append(LAB[RR[s],CC[s]]); G.append((RR[s]*30//5000)*1000+(CC[s]*30//5000))
        del F
    X=np.concatenate(X); Y=np.concatenate(Y); G=np.concatenate(G); ok=np.isfinite(X[:,NIS]).all(1)&np.isfinite(X[:,10:]).all(1); X,Y,G=X[ok],Y[ok],G[ok]
    # cap per class to 15000 after filtering
    keep=np.concatenate([rng.permutation(np.flatnonzero(Y==c))[:15000] for c in NAMES]); X,Y,G=X[keep],Y[keep],G[keep]
    ST['res']['samples_per_class']={NAMES[c]:int((Y==c).sum()) for c in NAMES}; sv()
    np.savez_compressed(ML+'lc_training_samples.npz',X=X,Y=Y,G=G,FN=np.array(FN))
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, cohen_kappa_score
    import joblib
    def mk(): return HistGradientBoostingClassifier(max_iter=300,learning_rate=0.08,max_leaf_nodes=48,min_samples_leaf=40,l2_regularization=1.0,class_weight='balanced',random_state=0)
    for name,cols in (('NISAR only',NIS),('NISAR + terrain',ALL)):
        ST['step']=f'CV {name}'; sv(); P=np.zeros_like(Y)
        for tr,te in GroupKFold(5).split(X,Y,G):
            P[te]=mk().fit(X[tr][:,cols],Y[tr]).predict(X[te][:,cols])
        cl=sorted(NAMES)
        ST['res'][name]=dict(overall_accuracy=round(float(accuracy_score(Y,P)),3),kappa=round(float(cohen_kappa_score(Y,P)),3),macro_F1=round(float(f1_score(Y,P,average='macro')),3),
            F1_per_class={NAMES[c]:round(float(v),2) for c,v in zip(cl,f1_score(Y,P,average=None,labels=cl))},confusion=confusion_matrix(Y,P,labels=cl).tolist(),labels=[NAMES[c] for c in cl]); sv()
    model=mk().fit(X[:,ALL],Y); joblib.dump(model,ML+'model_landcover_nisar_terrain.joblib')
    # ---- predict map + change vs 2020
    ST['step']='predict'; sv()
    PRED=np.zeros((Hd,Wd),np.uint8); PMAX=np.zeros((Hd,Wd),np.uint8); CHG=np.zeros((Hd,Wd),np.uint8)
    classes=model.classes_
    for r0 in range(0,Hd,300):
        r1=min(Hd,r0+300); F=feats(r0,r1); ok=M[r0:r1]&np.isfinite(F).all(-1)
        if ok.any():
            pr=model.predict_proba(F[ok]); c=classes[pr.argmax(1)]; p=pr.max(1)
            a=PRED[r0:r1]; a[ok]=c; PRED[r0:r1]=a; b=PMAX[r0:r1]; b[ok]=np.round(p*100).astype(np.uint8); PMAX[r0:r1]=b
        del F
    FOREST={1,2}; VEG={1,2,3,4,5,6}
    lab=np.where(M,LAB,0); conf=PMAX>=60; diff=(PRED>0)&(lab>0)&(PRED!=lab)&conf
    CHG[(PRED>0)&(lab>0)]=1                                        # 1 agrees with 2020 (or low confidence)
    CHG[diff]=5                                                    # 5 other class change
    CHG[diff&np.isin(lab,list(FOREST))&~np.isin(PRED,list(FOREST))]=2  # 2 forest -> non-forest
    CHG[diff&np.isin(lab,list(VEG))&(PRED==7)]=3                   # 3 vegetation -> built-up/bare
    CHG[diff&(PRED==8)&(lab!=8)]=4                                 # 4 land -> water
    for nm,a in (('NISAR_LandCover_2026_30m.tif',PRED),('NISAR_LandCover_confidence_pct_30m.tif',PMAX),('NISAR_LandCover_Change_vs_2020_30m.tif',CHG)):
        o=gdal.GetDriverByName('GTiff').Create(ML+nm,Wd,Hd,1,gdal.GDT_Byte,options=['COMPRESS=DEFLATE','TILED=YES']); o.SetGeoTransform(GT); o.SetProjection(prj); o.GetRasterBand(1).SetNoDataValue(0)
        o.GetRasterBand(1).WriteArray(np.where(M,a,0)); o=None
    tot=int(((PRED>0)&M).sum())
    ST['res']['mapped_km2']=round(tot*0.0009,1); ST['res']['change_km2']={k:round(float((CHG==v).sum())*0.0009,1) for k,v in (('agrees',1),('forest_to_nonforest',2),('veg_to_builtup_bare',3),('land_to_water',4),('other',5))}
    ST['res']['landslide_overlap']=None
    try:
        ls=gdal.Open(W+'out_v2/earthquake/Landslide_Candidates_after_8Jun2026_30m.tif').ReadAsArray()==2
        ST['res']['landslide_overlap']=dict(candidate_px=int(ls.sum()),with_nisar=int((ls&(PRED>0)).sum()),flagged_as_change=int((ls&(CHG>=2)).sum()),pct=round(100*float((ls&(CHG>=2)).sum())/max(int((ls&(PRED>0)).sum()),1),1),
            background_change_pct=round(100*float(((CHG>=2)&M).sum())/max(tot,1),1))
    except Exception as e: ST['res']['landslide_overlap']=repr(e)
    ST['step']='done'; ST['fin']=True; sv()
except Exception:
    ST['err']=traceback.format_exc(); sv()

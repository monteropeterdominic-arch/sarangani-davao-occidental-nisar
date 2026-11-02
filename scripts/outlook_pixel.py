# 12-month El Nino-conditioned outlook per pixel (MODIS EVI/NDVI 250 m, LST 1 km)
# anomaly(Y,lead k) = a_k + b_k*ONI_NDJ(Y) + c_k*init(Y) ; init = Sep-Oct anomaly of year Y
import json, traceback, math, numpy as np
from osgeo import gdal, osr
gdal.UseExceptions()
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; O=W+'outlook/'; LOG=W+'_outlook_px.json'; st={'stage':'start'}
def sv(): json.dump(st,open(LOG,'w'),indent=1)
ONI={2000:-0.78,2001:-0.32,2002:1.08,2003:0.29,2004:0.71,2005:-0.82,2006:0.88,2007:-1.72,2008:-0.79,2009:1.5,2010:-1.48,2011:-0.92,2012:-0.03,2013:-0.11,2014:0.8,2015:2.59,2016:-0.37,2017:-0.76,2018:1.05,2019:0.75,2020:-1.06,2021:-0.83,2022:-0.7,2023:1.99,2024:-0.43,2025:-0.6}
SCEN={'central':2.5,'low':2.0,'high':3.0}; YRS=list(range(2000,2026)); HOLD=[2002,2009,2015,2023]
G=json.load(open(O+'modis_grids.json')); MJ=json.load(open(O+'modis_aoi_monthly.json')); months=MJ['months']
def idx(y,m): v=y*100+m; return months.index(v) if v in months else None
def run(name,stack,scale,ndv,gt,X,Y,prov,strip):
    T=len(months); mon=np.array([m%100 for m in months]); tyr=np.array([(m//100)+(m%100-0.5)/12 for m in months])
    cl_sel=np.array([2001<=m//100<=2025 for m in months])
    leads=[(Y0,m) for Y0,m in [(0,11),(0,12)]+[(1,mm) for mm in range(1,11)]]
    fut=[(2026+dy,m) for dy,m in leads]
    out={k:np.full((12,Y,X),np.nan,np.float32) for k in ['anom','value','pbelow','sd']}
    beta=np.full((12,Y,X),np.nan,np.float32); r2=np.full((12,Y,X),np.nan,np.float32)
    hc={h:{'pred':np.zeros(12),'obs':np.zeros(12),'n':np.zeros(12),'sxy':np.zeros(12),'sxx':np.zeros(12),'syy':np.zeros(12),'sx':np.zeros(12),'sy':np.zeros(12)} for h in HOLD}
    for r0 in range(0,Y,strip):
        r1=min(Y,r0+strip); a=np.array(stack[:,r0:r1,:],dtype=np.float32)
        a=np.where(a==ndv,np.nan,a*scale) if ndv is not None else a
        a[:, ~(prov[r0:r1]>0)]=np.nan
        clim=np.full((13,)+a.shape[1:],np.nan,np.float32)
        for m in range(1,13):
            s=(mon==m)&cl_sel
            with np.errstate(all='ignore'): clim[m]=np.nanmedian(a[s],axis=0)
        an=a-clim[mon]
        # per-pixel linear trend of anomalies
        tt=(tyr-2013)[:,None,None]; ok=np.isfinite(an); n=ok.sum(0)
        with np.errstate(all='ignore'):
            tm=np.where(ok,tt,0).sum(0)/n; am=np.nansum(an,0)/n
            sl=np.nansum(np.where(ok,(tt-tm)*(an-am),0),0)/np.nansum(np.where(ok,(tt-tm)**2,0),0)
        sl=np.nan_to_num(sl); ic=np.nan_to_num(am-sl*tm)
        an=an-(ic+sl*tt)
        def init(Y0):
            i9,i10=idx(Y0,9),idx(Y0,10); v=[an[i] for i in (i9,i10) if i is not None]
            with np.errstate(all='ignore'): x=np.nanmean(np.stack(v),0)
            return np.nan_to_num(x)
        INIT={Y0:init(Y0) for Y0 in YRS+[2026]}
        sdi=np.nanstd(np.stack([INIT[y] for y in YRS]),0)+1e-6
        for k,(dy,m) in enumerate(leads):
            ys=[];xs=[]
            for Y0 in YRS:
                i=idx(Y0+dy,m); ys.append(an[i] if i is not None else np.full(an.shape[1:],np.nan)); xs.append(Y0)
            Yk=np.stack(ys); wk=np.isfinite(Yk); Yk=np.nan_to_num(Yk)
            on=np.array([ONI[y] for y in xs])[:,None,None]*np.ones_like(Yk); iv=np.stack([INIT[y] for y in xs])/sdi
            def fit(wmask):
                Xs=[np.ones_like(Yk),on,iv]; A=np.zeros(Yk.shape[1:]+(3,3)); b=np.zeros(Yk.shape[1:]+(3,))
                for p in range(3):
                    b[...,p]=(wmask*Xs[p]*Yk).sum(0)
                    for q in range(3): A[...,p,q]=(wmask*Xs[p]*Xs[q]).sum(0)
                A[...,1,1]+=1.0; A[...,2,2]+=1.0   # ridge
                A[...,0,0]+=1e-6
                return np.linalg.solve(A,b[...,None])[...,0]
            cf=fit(wk.astype(np.float32)); nn=wk.sum(0)
            pr=cf[...,0][None]+cf[...,1][None]*on+cf[...,2][None]*iv
            res=np.where(wk,Yk-pr,0); sd=np.sqrt((res**2).sum(0)/np.maximum(nn-3,1))
            vy=np.where(wk,(Yk-(Yk*wk).sum(0)/np.maximum(nn,1))**2,0).sum(0)
            with np.errstate(all='ignore'): R2=1-(res**2).sum(0)/vy
            good=nn>=12
            fi=INIT[2026]/sdi
            p=cf[...,0]+cf[...,1]*SCEN['central']+cf[...,2]*fi
            tfut=(2026+dy)+(m-0.5)/12-2013
            val=clim[m]+ic+sl*tfut+p
            out['anom'][k,r0:r1]=np.where(good,p,np.nan); out['value'][k,r0:r1]=np.where(good,val,np.nan); out['sd'][k,r0:r1]=np.where(good,sd,np.nan)
            z=-p/np.maximum(sd,1e-6); out['pbelow'][k,r0:r1]=np.where(good,0.5*(1+_erf(z/math.sqrt(2))),np.nan)
            beta[k,r0:r1]=np.where(good,cf[...,1],np.nan); r2[k,r0:r1]=np.where(good,R2,np.nan)
            for h in HOLD:
                j=YRS.index(h); w2=wk.astype(np.float32).copy(); w2[j]=0
                c2=fit(w2); ph=c2[...,0]+c2[...,1]*ONI[h]+c2[...,2]*iv[j]; ob=Yk[j]; v=wk[j]&good
                H=hc[h]; H['pred'][k]+=ph[v].sum(); H['obs'][k]+=ob[v].sum(); H['n'][k]+=v.sum()
                H['sx'][k]+=ph[v].sum(); H['sy'][k]+=ob[v].sum(); H['sxy'][k]+=(ph[v]*ob[v]).sum(); H['sxx'][k]+=(ph[v]**2).sum(); H['syy'][k]+=(ob[v]**2).sum()
        st[name]=r1; sv()
    return out,beta,r2,hc,fut
def _erf(x):
    # vectorised erf (Abramowitz-Stegun 7.1.26)
    s=np.sign(x); x=np.abs(x); t=1/(1+0.3275911*x)
    y=1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-0.284496736)*t+0.254829592)*t*np.exp(-x*x)
    return s*y
def save(arr,gt,X,Y,path):
    d=gdal.GetDriverByName('GTiff').Create(path,X,Y,arr.shape[0],gdal.GDT_Float32,['COMPRESS=DEFLATE','PREDICTOR=3'])
    d.SetGeoTransform(gt); d.SetProjection(G['prj'])
    for i in range(arr.shape[0]): b=d.GetRasterBand(i+1); b.WriteArray(arr[i]); b.SetNoDataValue(float('nan'))
    d=None
try:
    summ={}
    prov=np.load(O+'modis250_prov.npy'); provL=np.load(O+'modis1km_prov.npy')
    jobs=[('evi',np.load(O+'modis_evi_month.npy',mmap_mode='r'),1e-4,-3000,G['gt'],G['X'],G['Y'],prov,40),
          ('ndvi',np.load(O+'modis_ndvi_month.npy',mmap_mode='r'),1e-4,-3000,G['gt'],G['X'],G['Y'],prov,40),
          ('lst',np.load(O+'modis_lst_month.npy',mmap_mode='r'),1.0,None,G['gtL'],G['XL'],G['YL'],provL,155)]
    for name,stk,sc,nd,gt,X,Y,pv,strip in jobs:
        st['stage']=name; sv()
        out,beta,r2,hc,fut=run(name,stk,sc,nd,gt,X,Y,pv,strip)
        for k in out: save(out[k],gt,X,Y,O+f'Outlook_{name}_{k}_Nov2026_Oct2027.tif')
        save(beta,gt,X,Y,O+f'Outlook_{name}_ONIsensitivity_perdegC.tif'); save(r2,gt,X,Y,O+f'Outlook_{name}_R2.tif')
        m=pv>0; S={'months':[f'{y}-{mm:02d}' for y,mm in fut]}
        for k in ['anom','value','pbelow','sd']: S[k]=[float(np.nanmean(out[k][i][m])) for i in range(12)]
        S['beta']=[float(np.nanmean(beta[i][m])) for i in range(12)]; S['r2']=[float(np.nanmedian(r2[i][m])) for i in range(12)]
        S['frac_pbelow_gt_0.66']=[float(np.nanmean((out['pbelow'][i][m]>0.66))) for i in range(12)]
        for pnum,pn in [(1,'davocc'),(2,'sarangani')]:
            mm=pv==pnum; S['anom_'+pn]=[float(np.nanmean(out['anom'][i][mm])) for i in range(12)]
        hh={}
        for h,H in hc.items():
            n=np.maximum(H['n'],1); cov=H['sxy']/n-(H['sx']/n)*(H['sy']/n); vx=H['sxx']/n-(H['sx']/n)**2; vy=H['syy']/n-(H['sy']/n)**2
            hh[h]={'pred_mean':(H['pred']/n).tolist(),'obs_mean':(H['obs']/n).tolist(),'pixel_r':(cov/np.sqrt(np.maximum(vx*vy,1e-12))).tolist(),'n':H['n'].tolist()}
        S['hindcast']=hh; summ[name]=S
        json.dump(summ,open(O+'outlook_summary.json','w'))
    st['stage']='done'; st['fin']=True; sv()
except Exception:
    st['err']=traceback.format_exc()[-2000:]; sv()

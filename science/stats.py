import tifffile, numpy as np, json
import os
R=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')+'/'
rd=lambda f,s=0.01: (lambda a: np.where(a==-32768,np.nan,a*s))(tifffile.imread(R+f).astype('float32'))
hv=rd('nisar/backscatter_by_month/NISAR_HV_dB_2026-10_month_30m.tif')
slope=rd('terrain/SRTM_slope_deg_30m.tif'); asp=rd('terrain/SRTM_aspect_deg_30m.tif')
cls=tifffile.imread(R+'earthquake/Landslide_Candidates_after_8Jun2026_30m.tif')
lc=tifffile.imread(R+'nisar_ml/landcover/LandCover2020_8class_30m.tif')
ch=rd('nisar/change_by_month/NISAR_HH_Change_month_2026-09_to_2026-10_30m.tif')
cand=(cls==2); steep=(cls==1)&(slope>=20); land=np.isin(cls,[1,2,4,5])
names={1:'Closed forest',2:'Open forest & mangrove',3:'Brush/shrubs',4:'Grassland',5:'Annual crop',6:'Perennial crop',7:'Built-up & barren',8:'Water & fishpond'}
out={}
a=hv[cand]; b=hv[steep]; a=a[np.isfinite(a)]; b=b[np.isfinite(b)]
out['hv_cand']=np.round(a,2).tolist(); rng=np.random.default_rng(0); out['hv_steep_sample']=np.round(rng.choice(b,60000,replace=False),2).tolist()
out['hv_median']={'cand':float(np.median(a)),'steep':float(np.median(b)),'p25_steep':float(np.percentile(b,25)),'n_cand_px':int(a.size),'n_steep_px':int(b.size)}
# bootstrap CI of median difference
d=[np.median(rng.choice(a,a.size))-np.median(rng.choice(b,20000)) for _ in range(400)]
out['hv_median_diff_ci95']=[float(np.percentile(d,2.5)),float(np.percentile(d,97.5))]
strat=[]
for k,n in names.items():
    x=hv[cand&(lc==k)]; y=hv[steep&(lc==k)]; x=x[np.isfinite(x)]; y=y[np.isfinite(y)]
    if x.size>50: strat.append({'class':n,'n_cand_px':int(x.size),'cand':round(float(np.median(x)),2),'steep':round(float(np.median(y)),2),'diff':round(float(np.median(x)-np.median(y)),2)})
out['by_landcover']=strat
sb=[]
for lo,hi in [(20,28),(28,35),(35,45)]:
    x=hv[cand&(slope>=lo)&(slope<hi)]; y=hv[(cls==1)&(slope>=lo)&(slope<hi)]; x=x[np.isfinite(x)]; y=y[np.isfinite(y)]
    sb.append({'band':f'{lo}–{hi}°','n_cand_px':int(x.size),'cand':round(float(np.median(x)),2),'steep':round(float(np.median(y)),2),'diff':round(float(np.median(x)-np.median(y)),2)})
out['by_slope']=sb
f=lambda m: float(np.nanmean((asp[m]>=0)&(asp[m]<180)))
out['east_facing']={'cand':round(f(cand),3),'steep':round(f(steep),3)}
out['slope_median']={'cand':float(np.nanmedian(slope[cand])),'steep':float(np.nanmedian(slope[steep]))}
# El Nino: HH change by landcover
v=ch[land]; v=v[np.isfinite(v)]
out['hh_change']={'median':float(np.median(v)),'wetter_gt1':float((v>1).mean()),'drier_lt_m1':float((v<-1).mean()),'n':int(v.size)}
out['hh_change_hist']=np.histogram(np.clip(v,-6,6),bins=48,range=(-6,6))[0].tolist()
bylc=[]
for k,n in names.items():
    x=ch[land&(lc==k)]; x=x[np.isfinite(x)]
    if x.size>1000: bylc.append({'class':n,'median':round(float(np.median(x)),2),'wetter':round(float((x>1).mean()),3),'drier':round(float((x<-1).mean()),3),'area_km2':round(x.size*0.0009,0)})
out['hh_by_landcover']=bylc
json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'stats.json'),'w'))
print(json.dumps({k:v for k,v in out.items() if not k.startswith('hv_cand') and k not in('hv_steep_sample','hh_change_hist')},indent=1))

import json, os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
S=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'stats.json'))); M=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'inputs_summary.json')))
TEAL='#01665e'; TEAL2='#5ab4ac'; BROWN='#8c510a'; SAND='#dfc27d'; RED='#a50f15'; GOLD='#b8915a'; INK='#2a2620'; GREY='#9a9488'; BG='#fbf8f1'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.edgecolor':'#cfc6b4','axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':BG,'axes.facecolor':BG,'savefig.facecolor':BG,'axes.titleweight':'bold','axes.titlesize':12.5})
import textwrap
def finish(fig,title,sub,src,name):
    fig.suptitle(title,x=0.012,ha='left',fontsize=15,fontweight='bold',color=INK,y=0.985)
    w=int(fig.get_size_inches()[0]*11.8)
    fig.text(0.012,0.93,textwrap.fill(sub,w),ha='left',va='top',fontsize=10.5,color='#5b5446',linespacing=1.3)
    fig.text(0.012,0.012,'Source: '+src,ha='left',fontsize=8.5,color='#7a7263')
    fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),'figures',f'fig_{name}.png'),dpi=170); plt.close(fig)

# F1 clouds vs NISAR
ml=M['ml']; dates=[('28 Jun','20260628'),('20 Sep','20260920'),('23 Sep','20260923'),('25 Sep','20260925'),('30 Sep','20260930'),('17 Oct','20261017')]
clear=[ml['gapfill_'+d]['hls_clear_pct_of_those'] for _,d in dates]
fig,ax=plt.subplots(figsize=(10,4.6)); fig.subplots_adjust(top=0.74,bottom=0.12,left=0.08,right=0.98)
x=np.arange(len(dates)); ax.bar(x-0.2,[100]*len(x),0.38,color=TEAL,label='NISAR radar (sees through cloud)')
ax.bar(x+0.2,clear,0.38,color=SAND,label='Landsat HLS, clear view of the same ground (±8 days)')
for i,v in enumerate(clear): ax.text(i+0.2,v+2,f'{v:.1f}%',ha='center',fontsize=10,color=BROWN,fontweight='bold')
ax.set_xticks(x,[d for d,_ in dates]); ax.set_ylabel('% of NISAR footprint observed'); ax.set_ylim(0,112); ax.legend(frameon=False,loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=9.5)
finish(fig,'1 · Why NISAR: clouds hid the ground from optical satellites','On NISAR pass dates in 2026, Landsat could see only 0.5–62% of the same ground. NISAR saw all of it.','NISAR L2 GCOV footprints; NASA HLS L30 v2.0 Fmask, matched within 8 days (team data release 1505bae)','1_clouds_vs_nisar')

# F2 NISAR coverage timeline
cs=M['nisar_change_stats']; dd=['20260628','20260920','20260923','20260925','20260930','20261014','20261017','20261019']
lab=['28 Jun','20 Sep','23 Sep','25 Sep','30 Sep','14 Oct','17 Oct','19 Oct']; hh=[cs['HH_'+d]['cover_pct'] for d in dd]
fig,ax=plt.subplots(figsize=(10,4.2)); fig.subplots_adjust(top=0.78,bottom=0.16,left=0.08,right=0.98)
ax.bar(range(8),hh,color=TEAL2); ax.bar([9,10],[100,100],color=TEAL)
for i,v in enumerate(hh): ax.text(i,v+2,f'{v:.0f}%',ha='center',fontsize=9.5)
for i,v in zip([9,10],[100,100]): ax.text(i,v+2,'100%',ha='center',fontsize=9.5,fontweight='bold')
ax.set_xticks(list(range(8))+[9,10],lab+['Sep\nmonthly','Oct\nmonthly'],fontsize=9.5); ax.set_ylim(0,115); ax.set_ylabel('% of both provinces covered')
ax.axvline(8.0,color='#cfc6b4',lw=1,ls='--')
finish(fig,'2 · NISAR data used: 8 pass dates → 2 complete monthly maps','Single passes cover 19–89% of the area; averaging every pass in a month (in power) gives 100% coverage for September and October 2026.','NISAR L2 GCOV HH, provisional (ASF DAAC); team data release 1505bae','2_nisar_coverage')

# F3 landslide check
a=np.array(S['hv_cand']); b=np.array(S['hv_steep_sample']); md=S['hv_median']; ci=S['hv_median_diff_ci95']
fig,(ax,ax2)=plt.subplots(1,2,figsize=(12,5.1),gridspec_kw={'width_ratios':[1.25,1]}); fig.subplots_adjust(top=0.78,bottom=0.15,left=0.07,right=0.98,wspace=0.32)
bins=np.linspace(-26,-2,49); ax.hist(b,bins,density=True,color=GREY,alpha=.55,label='Other steep land (≥20°), no slide')
ax.hist(a,bins,density=True,color=RED,alpha=.7,label='Possible landslides (Landsat, 165 ha)')
ax.axvline(md['steep'],color='#55504a',lw=1.5,ls='--'); ax.axvline(md['cand'],color=RED,lw=1.5,ls='--')
ax.text(md['cand']-0.3,ax.get_ylim()[1]*0.55,f"median {md['cand']:.1f} dB",ha='right',color=RED,fontsize=9.5)
ax.text(md['steep']+0.3,ax.get_ylim()[1]*0.55,f"median {md['steep']:.1f} dB",ha='left',color='#55504a',fontsize=9.5)
ax.set_xlabel('NISAR HV backscatter, Oct 2026 (dB) — lower = less leaf/branch volume'); ax.set_ylabel('density'); ax.legend(frameon=False,fontsize=9,loc='upper left')
ax.set_title(f'Difference of medians {md["cand"]-md["steep"]:.2f} dB (95% CI {ci[0]:.2f} to {ci[1]:.2f})',fontsize=10.5,loc='left')
rows=[(x['class'],x['diff']) for x in S['by_landcover']]+[('Slope '+x['band'],x['diff']) for x in S['by_slope']]
y=np.arange(len(rows))[::-1]; ax2.barh(y,[r[1] for r in rows],color=[TEAL if i<len(S['by_landcover']) else GOLD for i in range(len(rows))])
for yi,(n,v) in zip(y,rows): ax2.text(-0.08,yi,f'{v:.1f} dB',va='center',ha='right',fontsize=9.5,color='white',fontweight='bold')
ax2.set_yticks(y,[r[0] for r in rows],fontsize=9.5); ax2.axvline(0,color=INK,lw=.8); ax2.set_xlabel('landslide minus same-class steep land (dB)')
ax2.set_title('Control: same land type / same steepness',fontsize=10.5,loc='left')
finish(fig,'3 · NISAR independently confirms the earthquake landslides','Spots Landsat flagged after the 8 June M7.8 quake are ~2.9 dB darker in NISAR HV — in every land-cover class and slope band. 194 of 305 patches (64%) fall in the darkest quarter of steep land (chance: 25%).','NISAR L2 GCOV HV Oct 2026; HLS L30 landslide candidates; SRTM slope; NAMRIA 2020 land cover. Bootstrap CI, 400 resamples.','3_landslide_check')

# F4 sensor comparison
fig,(ax,ax2)=plt.subplots(1,2,figsize=(12,4.6)); fig.subplots_adjust(top=0.78,bottom=0.14,left=0.07,right=0.98,wspace=0.3)
m=['Looks bare in NISAR HV\n(darkest quarter)','NISAR land-cover\nchange flagged']; s1=[64,27.9]; s2=[25,11.9]; x=np.arange(2)
ax.bar(x-0.2,s1,0.38,color=RED,label='Possible landslide spots'); ax.bar(x+0.2,s2,0.38,color=GREY,label='Other steep land')
for i,(p,q) in enumerate(zip(s1,s2)): ax.text(i-0.2,p+1.5,f'{p:g}%',ha='center',fontweight='bold',color=RED); ax.text(i+0.2,q+1.5,f'{q:g}%',ha='center',color='#55504a')
ax.set_xticks(x,m); ax.set_ylim(0,78); ax.set_ylabel('% of area / patches'); ax.legend(frameon=False,fontsize=9); ax.set_title('NISAR L-band (24 cm): separates them',loc='left',fontsize=10.5)
t=['Track 163\n(3→10 Jun)','Track 69\n(4→16 Jun)']; a1=[0.024,0.000]; a2=[0.025,-0.003]
ax2.bar(x-0.2,a1,0.38,color=RED,label='Possible landslide spots'); ax2.bar(x+0.2,a2,0.38,color=GREY,label='Other steep land')
for i,(p,q) in enumerate(zip(a1,a2)): ax2.text(i-0.2,p+0.002,f'{p:.3f}',ha='center',fontsize=9.5); ax2.text(i+0.2,max(q,0)+0.002,f'{q:.3f}',ha='center',fontsize=9.5)
ax2.set_xticks(x,t); ax2.set_ylim(-0.01,0.05); ax2.axhline(0,color=INK,lw=.8); ax2.set_ylabel('median coherence drop across the quake')
ax2.set_title('Sentinel-1 C-band (5.5 cm): no difference',loc='left',fontsize=10.5)
finish(fig,'4 · Longer radar waves see what shorter ones miss','Over dense tropical vegetation, Sentinel-1 coherence was already at noise level before the quake (>0.5 at only 9% and 2% of spots), so it could not tell landslides apart. NISAR could.','NISAR L2 GCOV; Copernicus Sentinel-1 interferograms via ASF HyP3 (A. Pelicano), with pre-quake control pairs','4_sensor_comparison')

# F5 El Nino rain + HH change
r=M['rain']; obs=r['observed_2025_2026']; ks=list(obs.keys()); clim=r['clim_by_month']
cl=[clim[int(k[5:7])-1] for k in ks]; pct=[obs[k]/c*100 for k,c in zip(ks,cl)]
fig,(ax,ax2)=plt.subplots(1,2,figsize=(12.5,5.1),gridspec_kw={'width_ratios':[1.5,1]}); fig.subplots_adjust(top=0.76,bottom=0.17,left=0.06,right=0.98,wspace=0.32)
cols=[BROWN if p<60 else (SAND if p<90 else TEAL2) for p in pct]
ax.bar(range(len(ks)),[obs[k] for k in ks],color=cols); ax.plot(range(len(ks)),cl,color=INK,lw=1.5,marker='o',ms=3,label='1998–2025 normal')
ax.set_xticks(range(len(ks)),[k[2:4]+'-'+k[5:] for k in ks],rotation=90,fontsize=8.5); ax.set_ylabel('rainfall (mm/month)')
for i,k in enumerate(ks):
    if k in('2026-08','2026-09','2026-10'): ax.text(i,obs[k]+8,f'{pct[i]:.0f}%',ha='center',fontsize=9,fontweight='bold',color=BROWN if pct[i]<60 else INK)
ax.legend(handles=[Patch(color=BROWN,label='<60% of normal'),Patch(color=SAND,label='60–90%'),Patch(color=TEAL2,label='≥90%'),plt.Line2D([],[],color=INK,marker='o',ms=3,label='1998–2025 normal')],frameon=False,fontsize=8.5,ncol=2,loc='upper right')
ax.set_title('GPM IMERG rainfall over both provinces',loc='left',fontsize=10.5)
hb=sorted(S['hh_by_landcover'],key=lambda x:x['wetter']); y=np.arange(len(hb))
ax2.barh(y,[x['wetter']*100 for x in hb],color=TEAL,label='wetter (> +1 dB)'); ax2.barh(y,[-x['drier']*100 for x in hb],color=BROWN,label='drier (< −1 dB)')
for yi,x_ in zip(y,hb): ax2.text(x_['wetter']*100+1,yi,f"{x_['wetter']*100:.0f}%",va='center',fontsize=9)
ax2.set_yticks(y,[x['class'] for x in hb],fontsize=9.5); ax2.axvline(0,color=INK,lw=.8); ax2.set_xlabel('% of class area, Sep → Oct 2026'); ax2.set_title('NISAR HH change by land type\nteal = wetter (> +1 dB) · brown = drier (< −1 dB)',loc='left',fontsize=10.5)
finish(fig,'5 · NISAR saw the soil wet up after the El Niño dry spell','Aug–Sep 2026 rain was 17–44% of normal; October recovered to 79%. NISAR HH rose by >1 dB on 32% of the land (10% fell) — most on exposed soils: annual crops 48%, grassland 40%, closed forest 28%.','NASA GPM IMERG V07 (GIOVANNI); NISAR L2 GCOV HH monthly composites; NAMRIA 2020 land cover','5_elnino_nisar')

# F6 outlook
o=M['outlook']['evi']; mo=r['months']; lab=[m_[2:4]+'-'+m_[5:] for m_ in mo]
fig,(ax,ax2)=plt.subplots(1,2,figsize=(12.5,4.7)); fig.subplots_adjust(top=0.78,bottom=0.14,left=0.06,right=0.97,wspace=0.3)
ax.fill_between(range(12),r['p10'],r['p90'],color=SAND,alpha=.5,label='10–90% range (ONI +2.0 to +3.0)'); ax.plot(range(12),r['central'],color=BROWN,lw=2,marker='o',ms=4,label='expected (ONI +2.5)')
ax.plot(range(12),r['clim'],color=INK,lw=1.3,ls='--',label='normal'); ax.set_xticks(range(12),lab,rotation=90,fontsize=8.5); ax.set_ylabel('rainfall (mm/month)'); ax.legend(frameon=False,fontsize=8.5)
ax.set_title('Rain: Nov–Apr expected at 41% of normal (341 vs 825 mm)',loc='left',fontsize=10.5)
ax2.bar(range(12),np.array(o['anom'])*100,color=[BROWN if v<-0.03 else SAND for v in o['anom']]); ax2.set_ylabel('area-mean EVI anomaly (×100)')
ax3=ax2.twinx(); ax3.plot(range(12),np.array(o['pbelow'])*100,color=RED,marker='o',ms=4,lw=1.8); ax3.set_ylabel('chance below normal (%)',color=RED); ax3.set_ylim(0,100); ax3.spines['top'].set_visible(False)
ax2.set_xticks(range(12),lab,rotation=90,fontsize=8.5); ax2.axhline(0,color=INK,lw=.8)
ax2.set_title('EVI anomaly (bars) and chance below normal (line):\nbrowning peaks Apr 2027, 78% chance',loc='left',fontsize=10.5)
finish(fig,'6 · The next 12 months under a strong El Niño (scenario, not a forecast)','Per-pixel model on 26 years of MODIS: anomaly = a + b·ONI(NDJ) + c·(Sep–Oct anomaly). Rainfall from a 1998–2025 GPM regression on ONI.','NASA MODIS MOD13Q1/MOD11A2 v6.1 (AppEEARS), GPM IMERG V07, NOAA CPC ONI (8 Oct 2026 outlook)','6_outlook')

# F7 hindcast + AI model
hc=r['hindcast']; yrs=list(hc.keys())
fig,(ax,ax2)=plt.subplots(1,2,figsize=(12.5,4.9)); fig.subplots_adjust(top=0.72,bottom=0.14,left=0.07,right=0.98,wspace=0.3)
x=np.arange(len(yrs)); ax.bar(x-0.2,[hc[y]['pred_NovApr_mm'] for y in yrs],0.38,color=GOLD,label='predicted (left out of training)'); ax.bar(x+0.2,[hc[y]['obs_NovApr_mm'] for y in yrs],0.38,color=TEAL,label='observed')
ax.axhline(825,color=INK,ls='--',lw=1.2); ax.text(4.45,840,'normal 825 mm',ha='right',fontsize=9)
ax.set_xticks(x,[f'{y}–{int(y[2:])+1:02d}\n(El Niño)' if y!='1998' else '1998–99\n(La Niña)' for y in yrs],fontsize=9); ax.set_ylabel('Nov–Apr rainfall (mm)'); ax.legend(frameon=False,fontsize=8.5,loc='upper right')
cv=r['cv_NovApr']; ax.set_title(f"Leave-one-year-out: r = {cv['r']}, MAE {cv['mae_mm']} mm (n = {cv['n_years']} yrs)",loc='left',fontsize=10.5)
k=['NISAR only','Terrain only','NISAR + terrain']; e=[M['ml']['EVI | '+i]['R2'] for i in k]; n=[M['ml']['NDMI | '+i]['R2'] for i in k]
xx=np.arange(3); ax2.bar(xx-0.2,e,0.38,color=TEAL,label='EVI (greenness)'); ax2.bar(xx+0.2,n,0.38,color=TEAL2,label='NDMI (leaf water)')
for i,(p,q) in enumerate(zip(e,n)): ax2.text(i-0.2,p+0.01,f'{p:.2f}',ha='center',fontsize=9.5); ax2.text(i+0.2,q+0.01,f'{q:.2f}',ha='center',fontsize=9.5)
lcm=M['ml_lc']; ax2.set_xticks(xx,k); ax2.set_ylabel('R² (spatial-block cross-validation)'); ax2.set_ylim(0,0.5); ax2.legend(frameon=False,fontsize=8.5,loc='upper left')
ax2.set_title(f"NISAR → vegetation AI, 167,858 px\nscored on 210 held-out 5 km blocks",loc='left',fontsize=10.5)
finish(fig,'7 · Validation: tested on years and places the models never saw','Left: El Niño rainfall model hindcasts. Right: NISAR → vegetation model, scored on held-out 5 km blocks (no spatial leakage). NISAR adds most of the skill; terrain alone is weak. NISAR land cover vs the 2020 map: accuracy 63%, κ 0.57.','GPM IMERG V07; NOAA ONI; NISAR L2 GCOV; HLS L30; SRTM; NAMRIA 2020 land cover (reference for κ)','7_validation')
print('ok')

import os, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, textwrap
BG='#fbf8f1'; INK='#2a2620'; HEAD='#3f5a4a'
def table(name,title,cols,rows,widths,src):
    wrap=[[textwrap.fill(str(c),w) for c,w in zip(r,widths)] for r in rows]
    lines=[max(c.count('\n')+1 for c in r) for r in wrap]
    H=1.1+0.27*sum(lines)+0.35*len(rows); W=sum(widths)*0.085+0.4
    fig=plt.figure(figsize=(W,H)); fig.patch.set_facecolor(BG); ax=fig.add_axes([0.012,0.35/H,0.976,1-1.0/H]); ax.axis('off')
    fig.text(0.012,1-0.18/H,title,fontsize=14,fontweight='bold',color=INK,va='top')
    tb=ax.table(cellText=wrap,colLabels=cols,colWidths=[w/sum(widths) for w in widths],bbox=[0,0,1,1],cellLoc='left')
    tb.auto_set_font_size(False); tb.set_fontsize(9.6)
    for (i,j),c in tb.get_celld().items():
        c.set_edgecolor('#d9cfbb'); c.set_linewidth(0.6)
        if i==0: c.set_facecolor(HEAD); c.get_text().set_color('white'); c.get_text().set_fontweight('bold'); pass
        else: c.set_facecolor('#fffdf8' if i%2 else '#f4eee1'); pass
        c.PAD=0.03
    fig.text(0.012,0.012,'Source: '+src,fontsize=8,color='#7a7263')
    fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),'figures',f'table_{name}.png'),dpi=170); plt.close(fig)

table('1_data','Table 1 · Data used (all free and open)',
 ['Dataset','Sensor / agency','Resolution','Dates used','What we used it for'],
 [['NISAR L2 GCOV (HH, HV), provisional','NISAR L-band SAR, NASA–ISRO (ASF DAAC)','30 m','8 dates: 28 Jun – 19 Oct 2026; Sep & Oct monthly composites','Core input: vegetation volume (HV), soil wetness (HH), Sep→Oct change, AI vegetation & land cover, landslide check, watch list'],
  ['HLS L30 v2.0 (EVI, NBR, NDMI, NDWI)','Landsat 8/9 OLI, NASA/USGS (LP DAAC)','30 m','Jan 2025 – Oct 2026','Landslide detection (pre/post 8 Jun vs 2025), AI training labels, quarterly El Niño series'],
  ['MODIS MOD13Q1 / MOD11A2 v6.1','Terra MODIS, NASA (AppEEARS)','250 m / 1 km','2000 – Oct 2026','Vegetation & land-temperature anomalies; 12-month El Niño outlook model'],
  ['GPM IMERG V07','GPM constellation, NASA (GIOVANNI)','0.1°','1998 – Oct 2026','Rainfall anomalies; ONI → Nov–Apr rain regression and hindcasts'],
  ['SRTM GL1','Space Shuttle radar, NASA','30 m','2000','Slope, aspect, elevation: landslide filter, AI features, stratified controls'],
  ['Sentinel-1 interferograms (7 pairs, 2 tracks)','Copernicus C-band SAR (ESA), ASF HyP3','HyP3 InSAR','pairs spanning 8 Jun 2026, plus pre-quake control pairs','Independent radar cross-check of the quake (by A. Pelicano)'],
  ['Oceanic Niño Index','NOAA CPC','monthly','1998 – Oct 2026 (+ 8 Oct 2026 outlook)','El Niño strength driving the outlook (ONI +2.5 °C scenario)'],
  ['2020 regional land cover (8 classes)','NAMRIA','30 m','2020','AI training reference, land-cover strata for controls, province boundary']],
 [26,28,10,24,40],'team data release github.com/monteropeterdominic-arch/sarangani-davao-occidental-nisar (commit 1505bae)')

table('2_hypotheses','Table 2 · Questions, tests and results (scientific method)',
 ['Hypothesis','Test (with control)','Result','Verdict'],
 [['H1. NISAR L-band can see vegetation loss at the quake landslides through clouds','Compare NISAR HV (Oct 2026) at 305 Landsat-flagged patches vs other steep land; repeat within each land-cover class and slope band; bootstrap CI','−2.86 dB (95% CI −3.02 to −2.66); −2.2 to −3.3 dB in every class/band; 194/305 patches (64%) in darkest quarter vs 25% expected','Supported'],
  ['H2. Shorter C-band radar (Sentinel-1) can do the same','Coherence drop across the quake at spots vs steep land, two tracks, with pre-quake control pairs','0.024 vs 0.025 and 0.000 vs −0.003; pre-quake coherence >0.5 at only 9% / 2% of spots','Not supported (signal at noise level)'],
  ['H3. NISAR HH tracks soil re-wetting after the El Niño dry spell','Sep→Oct HH change vs GPM rain (Sep 44% → Oct 79% of normal), by land-cover class','Median +0.46 dB; 32% of land wetter vs 10% drier; strongest on annual crops (48%) and grassland (40%), weakest under forest','Supported (consistent direction and physics)'],
  ['H4. El Niño strength predicts dry-season rain','GPM Nov–Apr rain vs ONI(NDJ), 1998–2025; leave-one-year-out','r = −0.81; LOYO r = 0.83, MAE 138 mm (n = 28); 2002, 2015 hindcasts within 10%','Supported; 2009 and 2023 over-predicted the drought'],
  ['H5. NISAR + terrain can estimate greenness under clouds','Gradient-boosted trees, 167,858 matched pixels, 210 held-out 5 km blocks','EVI R² 0.37 (NISAR only 0.24, terrain only 0.12); land cover κ 0.57','Partly supported: good for patterns, not single fields']],
 [30,40,46,22],'team analysis of NISAR, HLS, GPM, MODIS, SRTM and Sentinel-1 (data release 1505bae)')

table('3_limits','Table 3 · Uncertainty and limitations (what we do NOT claim)',
 ['Limitation','Why it matters','How we handled it'],
 [['NISAR data is provisional and starts 19 Jun 2026','No NISAR image before the 8 Jun quake, so no NISAR before/after','Used NISAR to check optical results after the quake; stated in app and docs'],
  ['Single passes cover 19–89% of the area','Per-date maps have gaps','Monthly power-mean composites (100% cover); gaps left blank, never filled with guesses'],
  ['Landslide patches are not field-checked','Rain, clearing or harvest can mimic slides','Called "possible landslides"; matched PhilSA estimate (165 vs ~138 ha); controls by class and slope'],
  ['AI vegetation map R² 0.37','Field-level values are uncertain','Shown as patterns; spatial-block CV to avoid leakage'],
  ['Outlook assumes ONI +2.5 °C (range 2.0–3.0)','Future El Niño may differ','Labelled a scenario with a 10–90% range and hindcasts, not a forecast'],
  ['Watch-first map combines model + radar','Priority guide, not a damage map','Updated monthly as new NISAR passes arrive (12-day repeat)']],
 [34,36,48],'team analysis; see docs and data release README')
print('ok')

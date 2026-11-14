# Cross-check: do the optical landslide spots lose Sentinel-1 coherence after the quake more than other steep land?
# And how often does NISAR flag the same places as changed land? Run in the QGIS Python console (GDAL + numpy).
import numpy as np
from osgeo import gdal
W='C:/Users/PETER DOMINIC/Documents/DavOcc_Sarangani/'; G=W+'s1_quake/geotiff/'
M=np.load(W+'tmp/aoi_mask30.npy'); rd=lambda p: gdal.Open(p).ReadAsArray()
ls=rd(W+'out_v2/earthquake/Landslide_Candidates_after_8Jun2026_30m.tif')   # 2 = possible landslide, 1 = checked, no slide
slope=rd(W+'ml/srtm_slope_30m.tif')
def coh(name):
    a=rd(G+name+'_coherence.tif').astype(np.float32); a[a==-32768]=np.nan; return a*1e-4
spots=M&(ls==2); steep=M&(ls==1)&(slope>=20)
for ctl,qk in [('t163_20260528_20260603_control','t163_20260603_20260610_quake'),('t069_20260522_20260603_control','t069_20260604_20260616_quake')]:
    c0,c1=coh(ctl),coh(qk); drop=c0-c1
    for nm,m in [('spots',spots),('steep no-slide',steep)]:
        ok=m&np.isfinite(drop)
        print(ctl[:4],nm,'median drop %.3f'%np.median(drop[ok]),'| coherence>0.5 before quake %.0f%%'%(100*(c0[ok]>0.5).mean()))
ch=rd(W+'ml/NISAR_LandCover_Change_vs_2020_30m.tif')   # 1 agree, 2-5 possible change
for nm,m in [('spots',spots),('steep no-slide',steep)]:
    v=ch[m&(ch>=1)&(ch<=5)]; print('NISAR',nm,'flagged change %.1f%%'%(100*(v>=2).mean()))

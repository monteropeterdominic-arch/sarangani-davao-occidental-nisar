import sqlite3, struct, numpy as np, tifffile, json, math
from PIL import Image, ImageDraw
import os
R=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')+'/'
rd=lambda f,s=0.01: (lambda a: np.where(a==-32768,np.nan,a*s))(tifffile.imread(R+f).astype('float32'))
hv=rd('nisar/backscatter_by_month/NISAR_HV_dB_2026-10_month_30m.tif'); slope=rd('terrain/SRTM_slope_deg_30m.tif')
cls=tifffile.imread(R+'earthquake/Landslide_Candidates_after_8Jun2026_30m.tif'); lcc=tifffile.imread(R+'nisar_ml/landcover/NISAR_LandCover_Change_vs_2020_30m.tif')
ref=(cls==1)&(slope>=20); p25=np.nanpercentile(hv[ref],25); med=np.nanmedian(hv[ref])
X0,Y0,P=648150.0,735870.0,30.0
def wkb_poly(b):
    # GPKG header
    flags=b[3]; env=(flags>>1)&7; hl=8+[0,32,48,48,64][env]; w=b[hl:]
    bo='<' if w[0]==1 else '>'; typ=struct.unpack(bo+'I',w[1:5])[0]; nr=struct.unpack(bo+'I',w[5:9])[0]; o=9; rings=[]
    for _ in range(nr):
        n=struct.unpack(bo+'I',w[o:o+4])[0]; o+=4; pts=struct.unpack(bo+'%dd'%(2*n),w[o:o+16*n]); o+=16*n; rings.append(list(zip(pts[0::2],pts[1::2])))
    return rings
def utm2ll(x,y,zone=51):
    a=6378137.0; f=1/298.257223563; e2=f*(2-f); k0=0.9996; x-=500000.0
    ep2=e2/(1-e2); M=y/k0; mu=M/(a*(1-e2/4-3*e2**2/64-5*e2**3/256)); e1=(1-math.sqrt(1-e2))/(1+math.sqrt(1-e2))
    p1=mu+(3*e1/2-27*e1**3/32)*math.sin(2*mu)+(21*e1**2/16-55*e1**4/32)*math.sin(4*mu)+(151*e1**3/96)*math.sin(6*mu)+(1097*e1**4/512)*math.sin(8*mu)
    C1=ep2*math.cos(p1)**2; T1=math.tan(p1)**2; N1=a/math.sqrt(1-e2*math.sin(p1)**2); R1=a*(1-e2)/(1-e2*math.sin(p1)**2)**1.5; D=x/(N1*k0)
    lat=p1-(N1*math.tan(p1)/R1)*(D**2/2-(5+3*T1+10*C1-4*C1**2-9*ep2)*D**4/24+(61+90*T1+298*C1+45*T1**2-252*ep2-3*C1**2)*D**6/720)
    lon=(D-(1+2*T1+C1)*D**3/6+(5-2*C1+28*T1-3*C1**2+8*ep2+24*T1**2)*D**5/120)/math.cos(p1)
    return round(math.degrees(lon)+(zone-1)*6-180+3,6), round(math.degrees(lat),6)
db=sqlite3.connect(R+'earthquake/Landslide_Candidates_after_8Jun2026_polygons.gpkg')
feats=[]; stats={'bare':0,'veg':0,'nodata':0}
for fid,g,val in db.execute('select fid,geom,val from candidates'):
    rings=wkb_poly(g); ext=rings[0]
    cols=[(x-X0)/P for x,y in ext]; rows=[(Y0-y)/P for x,y in ext]
    c0,c1=int(min(cols))-1,int(max(cols))+2; r0,r1=int(min(rows))-1,int(max(rows))+2
    im=Image.new('L',(c1-c0,r1-r0),0); ImageDraw.Draw(im).polygon([(c-c0,r-r0) for c,r in zip(cols,rows)],fill=1)
    mk=np.array(im,bool); sub=hv[r0:r1,c0:c1][mk]; sub=sub[np.isfinite(sub)]
    if sub.size==0:  # tiny polygon: centroid pixel 3x3
        rc,cc=int(np.mean(rows)),int(np.mean(cols)); sub=hv[rc-1:rc+2,cc-1:cc+2].ravel(); sub=sub[np.isfinite(sub)]
    m=float(np.median(sub)) if sub.size else None
    lcm=lcc[r0:r1,c0:c1][mk] if mk.sum() else np.array([])
    tag='nodata' if m is None else ('bare' if m<p25 else 'veg'); stats[tag]+=1
    cx=np.mean([x for x,y in ext]); cy=np.mean([y for x,y in ext]); lon,lat=utm2ll(cx,cy)
    feats.append({'type':'Feature','properties':{'id':fid,'nisar_hv_db':None if m is None else round(m,2),'nisar':tag,'nisar_lc_changed':bool(((lcm>1)&(lcm!=255)).any()) if lcm.size else False},'geometry':{'type':'Point','coordinates':[lon,lat]}})
json.dump({'type':'FeatureCollection','properties':{'steep_ref_hv_p25_db':round(float(p25),2),'steep_ref_hv_median_db':round(float(med),2)},'features':feats},open('landslides_nisar_check.geojson','w'))
print(stats, 'p25',round(float(p25),2), 'lc changed', sum(f['properties']['nisar_lc_changed'] for f in feats))
print(feats[0])

import numpy as np
import re
import Ngl
from sys import argv

wktype = "eps"
schemes = ["Ah21","Ah42","Ah63"]
resolutions = ["N256_dt075"]

path = argv[1]
te = []
en = []
for scheme in schemes:
    for resolution in resolutions:
        TE = [] 
        EN = []
        fname = path+"/swm_"+resolution+"_"+scheme+".out"
        fd = open(fname,"rb")
        for line in fd.readlines():
            a = re.findall("TE =\s*(\d+\.\d+)",line)
            if a: TE.append(float(a[0]))
            a = re.findall("Enstrophy =\s*(\d+\.\d+)",line)
            if a: EN.append(float(a[0]))
        TE = np.array(TE)
        TE = TE / TE[0]-1.0
        EN = np.array(EN)
        EN = EN / EN[0]-1.0
        te.append(TE)
        en.append(EN)
te = np.array(te)
en = np.array(en)
Nt_te = te.shape[1]
Nt_en = en.shape[1]
wk_res = Ngl.Resources()
wk_res.wkOrientation="Landscape"
wks = Ngl.open_wks(wktype, "barotropic_invariants",wk_res)
res = Ngl.Resources()
res.xyLineColors = ["red","green","blue","orange"]
res.xyLineThicknessF = 2
res.tmXBMode = "Explicit"
res.tmXBLabels = [0,10,20,30]
res.tmXBValues=[0,240,480,719]
res.tiXAxisString = "days"
res.nglFrame = False
res.nglDraw = False

print(te[:,-1])
print(en[:,-1])

lg_res = Ngl.Resources()
lg_res.vpWidthF = 0.12
lg_res.vpHeightF = 0.11
lg_res.lgLineThicknessF  = 2.0
lg_res.lgLineColors = ["red","green","blue"]
lg_res.lgDashIndexes = [0,0,0,0]
lg_res.lgLabelFontHeightF = 0.013
lg_res.lgPerimOn = True
lg1 = Ngl.legend_ndc(wks,len(schemes),schemes,0.1,0.45,lg_res)
#lg2 = Ngl.legend_ndc(wks,len(schemes),["Ah"+scheme for scheme in schemes],0.6,0.45,lg_res)

res.tiYAxisString = "relative error 10~S~-6"
res.tiMainString = "Total energy conservation"
plot_te = Ngl.xy(wks,np.arange(0,Nt_te),te*1e6,res)
res.tiYAxisString = "relative error 10~S~-2"
res.tiMainString = "Potential enstrophy conservation"
plot_en = Ngl.xy(wks,np.arange(0,Nt_en),en*1e2,res)
pres = Ngl.Resources()
#pres.nglPanelXWhiteSpacePercent=5.
#pres.nglMaximize=True
#pres.nglPaperOrientation="Landscape"
Ngl.panel(wks,(plot_te,plot_en),(1,2),pres)
Ngl.delete_wks(wks)
Ngl.end()

import numpy as np
import re
import Ngl

def plot_te_en(label,path,schemes,resolutions,wktype):

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
    wks = Ngl.open_wks(wktype, label+"_te_en",wk_res)
    res = Ngl.Resources()
    res.xyLineColors = ["red","green","blue","orange"]
    res.xyLineThicknessF = 2
    res.tmXBMode = "Explicit"
    if(Nt_te <= 360):
        res.tmXBLabels = [0,5,10,15]
        res.tmXBValues=[0,120,240,359]
    else:
        res.tmXBLabels = [0,10,20,30]
        res.tmXBValues=[0,240,480,719]
    res.tiXAxisString = "days"
    res.nglFrame = False
    res.nglDraw = False

    print(te[:,-1])
    print(en[:,-1])

    te_max = np.max(np.abs(te))
    te_pow = int(np.floor(np.log(te_max)/np.log(10)))
    en_max = np.max(np.abs(en))
    en_pow = int(np.floor(np.log(en_max)/np.log(10)))
    print "pows", [te_pow,en_pow]

    lg_res = Ngl.Resources()
    lg_res.vpWidthF = 0.12
    lg_res.vpHeightF = 0.11
    lg_res.lgLineThicknessF  = 2.0
    lg_res.lgLineColors = ["red","green","blue","orange"]
    lg_res.lgDashIndexes = [0,0,0,0]
    lg_res.lgLabelFontHeightF = 0.013
    lg_res.lgPerimOn = True
 
    scheme_labels = [scheme+"N~B~c~N~="+str(int(resolution[1:4])) for scheme in schemes for resolution in resolutions]
    lg1 = Ngl.legend_ndc(wks,len(schemes)*len(resolutions),scheme_labels,0.1,0.45,lg_res)
    #lg2 = Ngl.legend_ndc(wks,len(schemes),["Ah"+scheme for scheme in schemes],0.6,0.45,lg_res)

    res.tiYAxisString = "relative error 10~S~"+str(te_pow)
    res.tiMainString = "Total energy conservation"
    plot_te = Ngl.xy(wks,np.arange(0,Nt_te),te*10**(-te_pow),res)
    res.tiYAxisString = "relative error 10~S~"+str(en_pow)
    res.tiMainString = "Potential enstrophy conservation"
    plot_en = Ngl.xy(wks,np.arange(0,Nt_en),en*10**(-en_pow),res)
    print np.min(en*10**(-en_pow))
    pres = Ngl.Resources()
    #pres.nglPanelXWhiteSpacePercent=5.
    #pres.nglMaximize=True
    #pres.nglPaperOrientation="Landscape"
    Ngl.panel(wks,(plot_te,plot_en),(1,2),pres)

    Ngl.delete_wks(wks)

import numpy as np
import re
import Ngl
from sys import argv

path = argv[1]

def read_err(fname):
    fd = open(fname,"rb")
    l2_h = []
    linf_h = []
    l2_u = []
    linf_u = []
    for line in fd.readlines():
        errs = re.findall("\d+\.\d+E.\d+",line)
        l2_h.append(float(errs[0]))
        linf_h.append(float(errs[1]))
        l2_u.append(float(errs[2]))
        linf_u.append(float(errs[3]))
    return np.array(l2_h), np.array(linf_h),  np.array(l2_u),  np.array(linf_u)

schemes = ["21","42","43","63"]
path = path+"/errors_"
resolutions=["N032_dt600","N064_dt300","N128_dt150","N256_dt075"]
resources = Ngl.Resources()
resources.trYLog = True
resources.trYMinF = 1e-5
resources.trYMaxF = 10.0
resources.tmYMajorGrid = True
resources.tmXMajorGrid = True
for scheme in schemes:
    l2_errs = []
    for res in resolutions:
        fname = path+res+"_Ah"+scheme+".txt"
        l2_h, linf_h, l2, linf = read_err(fname)
        l2_errs.append(l2)
    wks = Ngl.open_wks("png", "l2_t_ah"+scheme)
    plot = Ngl.xy(wks,range(1,241),np.array(l2_errs),resources)
    Ngl.delete_wks(wks)

#l2_conv = []
#linf_conv = []
#schemes = ["21","42","43","63"]
#path = path+"/errors_"
#resolutions=["N032_dt600","N064_dt300","N128_dt150","N256_dt075"]
#
#for scheme in schemes:
#    for res in resolutions:
#        fname = path+res+"_Ah"+scheme+".txt"
#        l2, linf = read_err(fname)
#        l2_conv.append(np.max(l2))
#        linf_conv.append(np.max(linf))
#l2_order = []
#linf_order=[]
#x = np.array([32,64,128,256])
#l2_0 = {2: 1e-3, 3: 2e-4, 4: 6e-5}
#linf_0 = {2: 3e-3, 3: 5e-4, 4: 2e-4}
#for order in [2,3,4]:
#    for i in range(len(x)):
#        l2_order.append(l2_0[order]*(1.0*x[0]/x[i])**order)
#        linf_order.append(linf_0[order]*(1.0*x[0]/x[i])**order)
#
#l2 = np.array(l2_conv).reshape((len(schemes),len(resolutions)))
#linf = np.array(linf_conv).reshape((len(schemes),len(resolutions)))
#l2_ord = np.array(l2_order).reshape((3,len(resolutions)))
#linf_ord = np.array(linf_order).reshape((3,len(resolutions)))
#
#wks = Ngl.open_wks("png", "ts2_conv")
#res = Ngl.Resources()
#res.trXLog = True
#res.trYLog = True
#res.tiXAxisString="N~B~c"
#res.tmXBMode = "Explicit"
#res.tmXBValues = x
#res.tmXBLabels = x
#res.trXMinF = 19
#res.trXMaxF = 170
#res.nglDraw = False
#res.nglFrame = False
#res.xyLineThicknessF  = 2
#res.xyLineColors=["red","green","blue","orange"]
#res.xyMarkLineMode="MarkLines"
#res.xyMarker = 16
##res.xyLabelMode = "Custom"
##res.xyExplicitLabels = ['Ah21', 'Ah42', 'Ah63']
##res.xyLineLabelFontHeightF = 0.015
#
#res.tiYAxisString="~F10~l~B~2~N~ ~F21~error"
#plot1 =  Ngl.xy(wks,x,l2,res)
#res.tiYAxisString="~F10~l~B~~F34~%~N~~F21~ error"
#plot2 =  Ngl.xy(wks,x,linf,res)
#
#lres = Ngl.Resources()
#text_label_orders=['2nd order', '3rd order', '4th order']
#for i in range(3):
#    lres.gsLineLabelString = text_label_orders[i]
#    Ngl.add_polyline(wks, plot1, x, l2_ord[i,:],lres)
#    Ngl.add_polyline(wks, plot2, x, linf_ord[i,:],lres)
#
#lg_res = Ngl.Resources()
#lg_res.vpWidthF = 0.12
#lg_res.vpHeightF = 0.11
#lg_res.lgLineThicknessF  = 2.0
#lg_res.lgLineColors = ["red","green","blue","yellow"]
#lg_res.lgDashIndexes = [0,0,0,0]
#lg_res.lgLabelFontHeightF = 0.013
#lg_res.lgPerimOn = True
#lg1 = Ngl.legend_ndc(wks,len(schemes),["Ah"+scheme for scheme in schemes],0.1,0.45,lg_res)
#lg2 = Ngl.legend_ndc(wks,len(schemes),["Ah"+scheme for scheme in schemes],0.6,0.45,lg_res)
#
#Ngl.panel(wks,(plot1,plot2),(1,2))

import numpy as np
import re
import Ngl
from sys import argv
from add_cubed_sphere import add_cubed_sphere

path = argv[1]

schemes = ["Ah21","Ah42","Ah63"]
cn_res = Ngl.Resources()
cn_res.cnFillOn = True
cn_res.cnLinesOn = False
cn_res.cnLineLabelsOn = False
cn_res.mpCenterLonF = 180.0
cn_res.mpGridAndLimbOn = False
cn_res.lbOrientation="Horizontal"
cn_res.lbLabelFontHeightF = 0.015
cn_res.tiMainFontHeightF = 0.02
cn_res.mpGeophysicalLineColor = "Transparent"
cn_res.mpGreatCircleLinesOn = True

cn_res.nglDraw = False
cn_res.nglFrame = False

overlay_res = Ngl.Resources()
overlay_res.nglDraw = False
overlay_res.nglFrame = False
overlay_res.cnLineThicknessF = 2.0
overlay_res.cnInfoLabelOn = False
# overlay_res.cnGridA

for N in [20,40]:#[20,40,80,160]:
    Nlon = 4*N
    Nlat = 2*N+1
    cn_res.sfXArray = np.linspace(0.0,360.0,Nlon,endpoint=False)
    cn_res.sfYArray = np.linspace(-90.0,90.0,Nlat)
    overlay_res.sfXArray = np.linspace(0.0,360.0,Nlon,endpoint=False)
    overlay_res.sfYArray = np.linspace(-90.0,90.0,Nlat)
    wksres = Ngl.Resources()
    wksres.wkOrientation="portrait"
    wks = Ngl.open_wks("eps", "ts2_h_error_N{:03d}".format(N),wksres)
    Ngl.define_colormap(wks,"GMT_polar")
    plots = []
    for scheme in schemes:
        fd = open("h_N"+"{:03d}_".format(N)+scheme+".dat","rb")
        h0 = np.fromfile(fd,count=Nlon*Nlat,dtype=np.float32).reshape(Nlat,Nlon)
        fd.seek(4*Nlon*Nlat*10,0)
        h1 = np.fromfile(fd,count=Nlon*Nlat,dtype=np.float32).reshape(Nlat,Nlon)
        fd.close()
        cn_res.tiMainString = scheme
        plots.append(Ngl.contour_map(wks,h1-h0,cn_res))
        overlay_plot = Ngl.contour(wks, h0, overlay_res)
        Ngl.overlay(plots[-1], overlay_plot)
        add_cubed_sphere(wks,plots[-1])

    textres = Ngl.Resources()
    textres.txFontHeightF = 0.020
    #Ngl.text_ndc(wks,"Height error at day 10 [m]",0.5,.97,textres)
    pres = Ngl.Resources()
    Ngl.panel(wks,plots,(3,1),pres)
    Ngl.delete_wks(wks)

def read_err(fname):
    fd = open(fname,"rb")
    l2_h = []
    linf_h = []
    for line in fd.readlines():
        errs = re.findall("\d+\.\d+E.\d+",line)
        l2_h.append(float(errs[0]))
        linf_h.append(float(errs[1]))
    return np.array(l2_h), np.array(linf_h)

l221, linf = read_err(path+"/errors_N160_dt100_Ah21.txt")
l242, linf = read_err(path+"/errors_N160_dt100_Ah42.txt")
l263, linf = read_err(path+"/errors_N160_dt100_Ah63.txt")


wks = Ngl.open_wks("eps", "l2_h_t")

res = Ngl.Resources()
res.trYLog = True

plot = Ngl.xy(wks,range(1,241),np.array([l221,l242,l263]),res)

Ngl.delete_wks(wks)

l2_conv = []
linf_conv = []
schemes = ["21","42","63"]
path = path+"/errors_"
resolutions=["N020_dt800","N040_dt400","N080_dt200","N160_dt100"]

for scheme in schemes:
    for res in resolutions:
        fname = path+res+"_Ah"+scheme+".txt"
        l2, linf = read_err(fname)
        l2_conv.append(np.max(l2))
        linf_conv.append(np.max(linf))
l2_order = []
linf_order=[]
x = np.array([20,40,80,160])
l2_0 = {2: 1e-3, 3: 2e-4, 4: 2e-5}
linf_0 = {2: 3e-3, 3: 5e-4, 4: 0.7e-4}
for order in [2,3,4]:
    for i in range(len(x)):
        l2_order.append(l2_0[order]*(1.0*x[0]/x[i])**order)
        linf_order.append(linf_0[order]*(1.0*x[0]/x[i])**order)

l2 = np.array(l2_conv).reshape((len(schemes),len(resolutions)))
linf = np.array(linf_conv).reshape((len(schemes),len(resolutions)))
l2_ord = np.array(l2_order).reshape((3,len(resolutions)))
linf_ord = np.array(linf_order).reshape((3,len(resolutions)))

wks = Ngl.open_wks("eps", "ts2_conv")
res = Ngl.Resources()
res.trXLog = True
res.trYLog = True
res.tiXAxisString="N~B~c"
res.tmXBMode = "Explicit"
res.tmXBValues = x
res.tmXBLabels = x
res.trXMinF = 19
res.trXMaxF = 170
res.nglDraw = False
res.nglFrame = False
res.xyLineThicknessF  = 2
res.xyLineColors=["red","green","blue","orange"]
res.xyMarkLineMode="MarkLines"
res.xyMarker = 16
#res.xyLabelMode = "Custom"
#res.xyExplicitLabels = ['Ah21', 'Ah42', 'Ah63']
#res.xyLineLabelFontHeightF = 0.015

res.tiYAxisString="~F10~l~B~2~N~ ~F21~error"
plot1 =  Ngl.xy(wks,x,l2,res)
res.tiYAxisString="~F10~l~B~~F34~%~N~~F21~ error"
plot2 =  Ngl.xy(wks,x,linf,res)

lres = Ngl.Resources()
text_label_orders=['2nd order', '3rd order', '4th order']
for i in range(3):
    lres.gsLineLabelString = text_label_orders[i]
    Ngl.add_polyline(wks, plot1, x, l2_ord[i,:],lres)
    Ngl.add_polyline(wks, plot2, x, linf_ord[i,:],lres)

lg_res = Ngl.Resources()
lg_res.vpWidthF = 0.12
lg_res.vpHeightF = 0.11
lg_res.lgLineThicknessF  = 2.0
lg_res.lgLineColors = ["red","green","blue"]
lg_res.lgDashIndexes = [0,0,0,0]
lg_res.lgLabelFontHeightF = 0.013
lg_res.lgPerimOn = True
lg1 = Ngl.legend_ndc(wks,len(schemes),["Ah"+scheme for scheme in schemes],0.1,0.45,lg_res)
lg2 = Ngl.legend_ndc(wks,len(schemes),["Ah"+scheme for scheme in schemes],0.6,0.45,lg_res)

Ngl.panel(wks,(plot1,plot2),(1,2))

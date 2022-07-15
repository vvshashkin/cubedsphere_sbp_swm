import numpy as np
import re
import Ngl
from sys import argv
from add_cubed_sphere import add_cubed_sphere

from plot_te_en_cons import plot_te_en

wktype = "eps"
path = argv[1]
schemes = ["Ah21","Ah42","Ah43","Ah63"]

plot_te_en("ts5",path,schemes,["N040_dt400"],wktype)

cn_res = Ngl.Resources()
cn_res.cnFillOn = True
#cn_res.cnLinesOn = False
cn_res.cnLineLabelsOn = False
cn_res.mpCenterLonF = 180.0
cn_res.mpGridAndLimbOn = False
cn_res.lbOrientation="Horizontal"
cn_res.lbLabelFontHeightF = 0.015
cn_res.tiMainFontHeightF = 0.02
cn_res.mpGeophysicalLineColor="Transparent"
cn_res.mpGreatCircleLinesOn = True

cn_res.nglDraw = False
cn_res.nglFrame = False

orog_nlat = 181
orog_nlon = 360
orog_lon = np.linspace(0.0,2.0*np.pi,orog_nlon,endpoint=False)
orog_lat = np.linspace(-0.5*np.pi,0.5*np.pi,orog_nlat)
orog = np.empty((orog_nlat,orog_nlon),dtype = np.float32)
orog_phi0 = np.pi / 6.0
orog_lam0 = 0.5*np.pi
orog_r = np.pi / 9.0
orog_h = 2000.0
for j in range(orog_nlat):
    r = np.arccos(np.sin(orog_phi0)*np.sin(orog_lat[j])+\
                  np.cos(orog_phi0)*np.cos(orog_lat[j])*np.cos(orog_lam0-orog_lon))
    orog[j,:] = orog_h*(1-np.minimum(orog_r,r)/orog_r)


overlay_res = Ngl.Resources()
overlay_res.nglDraw = False
overlay_res.nglFrame = False
overlay_res.cnLineThicknessF = 2.0
overlay_res.sfXArray = np.linspace(0.0,360.0,orog_nlon,endpoint=False)
overlay_res.sfYArray = np.linspace(-90.0,90.0,orog_nlat)
overlay_res.cnLineLabelsOn = False
overlay_res.cnInfoLabelOn = False
overlay_res.cnLevelSelectionMode = "ExplicitLevels"
overlay_res.cnLevels = [0.5]

for N in [20,40,80,160]:#[32,48,64,96]:
    Nlon = 4*N
    Nlat = 2*N+1
    cn_res.sfXArray = np.linspace(0.0,360.0,Nlon,endpoint=False)
    cn_res.sfYArray = np.linspace(-90.0,90.0,Nlat)
    for scheme in schemes:
        wks = Ngl.open_wks(wktype, "ts5_h_"+scheme+"_N{:03d}".format(N))
        Ngl.define_colormap(wks,"rainbow+white")
        fd = open(path+"/h_N"+"{:03d}_".format(N)+scheme+".dat","rb")
        fd.seek(4*Nlon*Nlat*15,0)
        h15 = np.fromfile(fd,count=Nlon*Nlat,dtype=np.float32).reshape(Nlat,Nlon)
        fd.close()
        fd = open(path+"/curl_N"+"{:03d}_".format(N)+scheme+".dat","rb")
        fd.seek(4*Nlon*Nlat*15,0)
        z15 = np.fromfile(fd,count=Nlon*Nlat,dtype=np.float32).reshape(Nlat,Nlon)
        fd.close()
        cn_res.tiMainString = "total height field [m] "+scheme + " N~B~c~N~="+str(N)
        cn_res.cnLevelSelectionMode = "ExplicitLevels"
        #cn_res.cnLevels = [5100,5200,5300,5400,5500,5600,5700,5800,5900]
        cn_res.cnLevels = [5050,5100,5150,5200,5250,5300,5350,5400,5450,5500,5550,5600,5650,5700,5750,5800,5850,5900,5950]
        cn_res.cnFillColors = list(np.linspace(33,237,21,dtype=int))#len(cn_res.cnLevels)+1))
        plot1 = Ngl.contour_map(wks,h15,cn_res)
        print(np.max(orog),np.min(orog))
        plot2 = Ngl.contour(wks,orog,overlay_res)
        add_cubed_sphere(wks,plot1)

        cn_res.tiMainString = "relative vorticity field 10~S~-5~N~[s] "+scheme + " N~B~c~N~="+str(N)
        cn_res.cnLevels = [-3,-2.5,-2,-1.5,-1,-0.5,0.5,1.,1.5,2,2.5,3,3.5,4,4.5,5]
        cn_res.cnFillColors = [17,33,49,65,81,97,238,177,183,189,195,201,207,213,219,225,231]
        plot3 = Ngl.contour_map(wks,z15*1e5,cn_res)
        plot4 = Ngl.contour(wks,orog,overlay_res)
        add_cubed_sphere(wks,plot3)
        Ngl.overlay(plot1,plot2)
        Ngl.overlay(plot3,plot4)

        Ngl.panel(wks,(plot1,plot3),(1,2))
        #Ngl.draw(plot1)
        #Ngl.frame(wks)
        Ngl.delete_wks(wks)


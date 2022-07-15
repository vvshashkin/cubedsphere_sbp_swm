import numpy as np
import re
import Ngl
from sys import argv

path = argv[1]
if len(argv)>2 :
    stswm_sol_path = argv[2]
else:
    stswm_sol_path = None
schemes = ["Ah21","Ah42","Ah63"]
#schemes = ["Ah42"]
t1 = 401
t2 = 2400
Nc = 96
Nlon_stswm = 388
Nlat_stswm = 194
wktype = "eps"


def add_cubed_sphere(wks,plot):
    clat = 180.0*np.arcsin(1.0/np.sqrt(3))/np.pi
    line_res = Ngl.Resources()
    line_res.gsLineThicknessF = 1.0
    line_res.gsLineDashPattern = 2
    l1 = Ngl.add_polyline(wks,plot,(45,45,315,315,45),(clat,-clat,-clat,clat,clat),line_res)
    l2 = Ngl.add_polyline(wks,plot,(135,135,45,45,135),(clat,-clat,-clat,clat,clat),line_res)
    l3 = Ngl.add_polyline(wks,plot,(225,225,135,135,225),(clat,-clat,-clat,clat,clat),line_res)
    l4 = Ngl.add_polyline(wks,plot,(315,315,225,225,315),(clat,-clat,-clat,clat,clat),line_res)
    return (l1,l2,l3,l4)

def get_mean_and_std(fname,t1,t2,nx,ny):
    Nt = t2-t1+1.0
    f_mean  = np.zeros((ny,nx),dtype=np.float64)
    f_mean2 = np.zeros((ny,nx),dtype=np.float64)
    fd = open(fname,"rb")
    fd.seek(4*t1*nx*ny,0)
    for i in range(t2-t1+1):
        f = np.fromfile(fd,count=nx*ny,dtype=np.float32).reshape((ny,nx))
        f_mean  = f_mean+f/Nt
        f_mean2 = f_mean2+f**2/Nt

    fd.close()
    return f_mean, np.sqrt(f_mean2-f_mean**2)

def plot_Eldred(scheme, N, path,t1,t2):

    Nlon = 4*N
    Nlat = 2*N+1
    suffix = "_N{:03d}_".format(N)+scheme
    
    cn_res = Ngl.Resources()
    cn_res.cnFillOn = True
    cn_res.cnLinesOn = False
    #cn_res.cnLineThicknessF = 0.5
    cn_res.cnLineLabelsOn = False
    cn_res.mpCenterLonF = 180.0
    cn_res.mpGridAndLimbOn = False
    cn_res.mpGeophysicalLineColor = "Transparent"
    cn_res.mpGreatCircleLinesOn = True
    cn_res.lbLabelBarOn = True
    cn_res.lbBoxMinorExtentF = 0.12
    cn_res.lbOrientation="Horizontal"
    cn_res.lbLabelFontHeightF = 0.015
    cn_res.tiMainFontHeightF = 0.01
    cn_res.sfXArray = np.linspace(0.0,360.0,Nlon,endpoint=False)
    cn_res.sfYArray = np.linspace(-90.0,90.0,Nlat)
    cn_res.tiMainFontHeightF = 0.02
    cn_res.nglDraw = False
    cn_res.nglFrame = False
   
   
    wkres = Ngl.Resources()
    wkres.nglMaximize = True
    wkres.wkColorMap = "cmp_b2r"
    wks = Ngl.open_wks(wktype, "Eldred_curl_N{:03d}_".format(N)+scheme,wkres)

    cn_res.cnLevelSelectionMode = "ExplicitLevels"
    cn_res.cnLevels = [-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3,3.5,4]
    cn_res.cnFillColors = [2,6,10,14,18,22,26,30,0,37,41,45,49,53,57,61,65]

    fd = open(path+"/curl"+suffix+".dat","rb")
    fd.seek(4*t2*Nlon*Nlat,0)
    z = np.fromfile(fd,count=Nlon*Nlat,dtype=np.float32).reshape((Nlat,Nlon))
    fd.close()
    plot = Ngl.contour_map(wks,z*1e5,cn_res)
    cs = add_cubed_sphere(wks,plot)
    Ngl.draw(plot)
    Ngl.frame(wks)
    Ngl.delete_wks(wks)


    z_mean, z_std     = get_mean_and_std(path+"/curl"+suffix+".dat", t1, t2, Nlon, Nlat)
    div_mean, div_std = get_mean_and_std(path+"/div"+suffix+".dat", t1, t2, Nlon, Nlat)
    h_mean, h_std     = get_mean_and_std(path+"/h"+suffix+".dat", t1, t2, Nlon, Nlat)

    wks = Ngl.open_wks(wktype, "Eldred_N{:03d}_".format(N)+scheme,wkres)
   
    cn_res.tiMainString="mean relative vorticity 10~S~-5~N~ [1/s]" 
    cn_res.cnLevels = [-2,-1.75,-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25,1.5,1.75,2]
    cn_res.cnFillColors = [2,6,10,14,18,22,26,30,0,37,41,45,49,53,57,61,65]
    plot1 = Ngl.contour_map(wks,z_mean*1e5,cn_res)
    cs1 = add_cubed_sphere(wks,plot1)

    cn_res.tiMainString="relative vorticity std 10~S~-5~N~ [1/s]" 
    cn_res.cnLevels = [0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
    cn_res.cnFillColors = [2,8,14,20,26,32,38,45,52,59,65]
    plot2 = Ngl.contour_map(wks,z_std*1e5,cn_res)
    cs2 = add_cubed_sphere(wks,plot2)

    cn_res.tiMainString="mean divergence 10~S~-8~N~ [1/s]" 
    cn_res.cnLevels = [-2,-1.75,-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25,1.5,1.75,2]
    cn_res.cnFillColors = [2,6,10,14,18,22,26,30,0,37,41,45,49,53,57,61,65]
    plot3 = Ngl.contour_map(wks,div_mean*1e8,cn_res)
    #cs3 = add_cubed_sphere(wks,plot3)

    cn_res.tiMainString="divergence std 10~S~-8~N~ [1/s]" 
    cn_res.cnLevels = [4,8,12,16,20,24,28,32,36,40,44]
    cn_res.cnFillColors = [2,8,14,20,26,32,38,44,50,56,62,65]
    plot4 = Ngl.contour_map(wks,div_std*1e8,cn_res)
    cs4 = add_cubed_sphere(wks,plot4)

    cn_res.tiMainString="mean height [m]" 
    cn_res.cnLevels = [550,600,650,700,750,800,850,900,950,1000,1050,1100,1150]
    cn_res.cnFillColors = [2,7,12,17,22,27,32,37,42,47,52,56,61,65]
    plot5 = Ngl.contour_map(wks,h_mean,cn_res)
    cs5 = add_cubed_sphere(wks,plot5)

    cn_res.tiMainString="height std [m]" 
    cn_res.cnLevels = [4,8,12,16,20,24,28,32,36,40,44]
    cn_res.cnFillColors = [2,8,14,20,26,32,38,44,50,56,62,65]
    plot6 = Ngl.contour_map(wks,h_std,cn_res)
    cs6 = add_cubed_sphere(wks,plot6)

    pres = Ngl.Resources()
    Ngl.panel(wks,(plot5,plot6,plot1,plot2,plot3,plot4),(3,2),pres)
    Ngl.delete_wks(wks)

    #wks = Ngl.open_wks(wktype, "xy_N{:03d}_".format(N)+scheme,wkres)
    #xy_res = Ngl.Resources()
    #xy_res.nglDraw = False
    #xy_res.nglFrame = False
    #plot1 = Ngl.xy(wks,cn_res.sfYArray,np.mean(h_std,axis=1),xy_res)
    #plot2 = Ngl.xy(wks,cn_res.sfYArray,np.mean(z_std*1e5,axis=1),xy_res)
    #Ngl.panel(wks,(plot1,plot2),(1,2),pres)
    #Ngl.delete_wks(wks)
    return np.mean(z_std,axis=1), np.mean(div_std,axis=1)

z_stds = np.empty((0,2*Nc+1))
div_stds = np.empty((0,2*Nc+1))
for scheme in schemes:
    z_std, div_std = plot_Eldred(scheme,Nc,path,t1,t2)
    z_stds = np.vstack((z_stds,z_std))
    div_stds = np.vstack((div_stds,div_std))
    print div_std.shape

if(stswm_sol_path):
    stswm_lat = np.fromfile(stswm_sol_path+"/lats_stswm_128.txt",sep=",")
    z_mean, z_std     = get_mean_and_std(stswm_sol_path+"/Eldred_z.dat", t1, t2, Nlon_stswm, Nlat_stswm)
    d_mean, d_std     = get_mean_and_std(stswm_sol_path+"/Eldred_d.dat", t1, t2, Nlon_stswm, Nlat_stswm)
    z_std_stswm = np.mean(z_std,axis=1)
    d_std_stswm = np.mean(d_std,axis=1)

wk_res = Ngl.Resources()
wk_res.wkOrientation="Landscape"
wks = Ngl.open_wks(wktype, "z_div_std".format(Nc),wk_res)
xy_res = Ngl.Resources()
xy_res.nglDraw = False
xy_res.nglFrame = False
xy_res.xyLineColors=["red","green","blue"]
xy_res.xyLineThicknessF=2.0

x = np.linspace(-90.0,90.0,2*Nc+1)

xy_res.tiMainString="divergence std 10~S~-8 ~N~ [1/s]"
plot1 = Ngl.xy(wks,x,div_stds*1e8,xy_res)
if(stswm_sol_path):
    line_res = Ngl.Resources()
    line_res.gsLineColor="orange"
    line_res.gsLineThicknessF = 2.0
    line1 = Ngl.add_polyline(wks,plot1,stswm_lat,d_std_stswm*1e8,line_res)
xy_res.tiMainString="relative vorticity std 10~S~-5 ~N~ [1/s]"
plot2 = Ngl.xy(wks,x,z_stds*1e5,xy_res)
if(stswm_sol_path): line2 = Ngl.add_polyline(wks,plot2,stswm_lat,z_std_stswm*1e5,line_res)

lg_res = Ngl.Resources()
lg_res.vpWidthF = 0.12
lg_res.vpHeightF = 0.13
lg_res.lgLineThicknessF  = 2.0
lg_res.lgLineColors = ["red","green","blue"]
if(stswm_sol_path): lg_res.lgLineColors = ["red","green","blue","orange"]
lg_res.lgDashIndexes = [0,0,0,0,0]
lg_res.lgLabelFontHeightF = 0.013
lg_res.lgPerimOn = True
if(stswm_sol_path):
    lg1 = Ngl.legend_ndc(wks,4,["Ah21","Ah42","Ah63","T128"],0.2,0.67,lg_res)
else :
    lg1 = Ngl.legend_ndc(wks,3,["Ah21","Ah42","Ah63"],0.2,0.67,lg_res)
#lg2 = Ngl.legend_ndc(wks,4,["Ah21","Ah42","Ah43","Ah63"],0.72,0.65,lg_res)

pres = Ngl.Resources()
pres.nglPanelXWhiteSpacePercent = 5.0
Ngl.panel(wks,(plot1,plot2),(1,2),pres)
Ngl.delete_wks(wks)

import Ngl
import numpy as np

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

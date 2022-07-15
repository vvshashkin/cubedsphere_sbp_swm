module sph_coords_mod

use const_mod, only : pi

implicit none

contains

pure subroutine sph2cart(lam, phi, x, y, z)

    real(kind=8), intent(in)  :: lam, phi
    real(kind=8), intent(out) :: x, y, z

    x = cos(phi)*cos(lam)
    y = cos(phi)*sin(lam)
    z = sin(phi)

end subroutine sph2cart

pure subroutine cart2sph(x, y, z, lam, phi)

    real(kind=8), intent(in)  :: x, y, z
    real(kind=8), intent(out) :: lam, phi

    lam = atan2(y,x)
    phi = atan2(z,sqrt(x**2 + y**2))

    lam = lam+pi*(1.0_8-sign(1.0_8, lam))

end subroutine cart2sph

subroutine sph2cart_vec(lam, phi, V_lam, V_phi, Vx, Vy, Vz)

    real(kind=8), intent(in)  :: lam, phi, V_lam, V_phi
    real(kind=8), intent(out) :: Vx, Vy, Vz

    Vx = -sin(lam)*V_lam - sin(phi)*cos(lam)*V_phi
    Vy =  cos(lam)*V_lam - sin(phi)*sin(lam)*V_phi
    Vz =                            cos(phi)*V_phi

end subroutine sph2cart_vec

subroutine cart2sph_vec(lam, phi, Vx, Vy, Vz, V_lam, V_phi)

    real(kind=8), intent(in)  :: lam, phi, Vx, Vy, Vz
    real(kind=8), intent(out) :: V_lam, V_phi

    V_lam =          -sin(lam)*Vx +          cos(lam)*Vy
    V_phi = -cos(lam)*sin(phi)*Vx - sin(lam)*sin(phi)*Vy + cos(phi)*Vz

end subroutine cart2sph_vec

subroutine rotate_3D_y(x, y, z, theta, xn, yn, zn)

    real(kind=8), intent(in)  :: x, y, z, theta
    real(kind=8), intent(out) :: xn, yn, zn

    xn =  cos(theta)*x - sin(theta)*z
    zn =  sin(theta)*x + cos(theta)*z
    yn = y

end subroutine rotate_3D_y

subroutine sph_rotate_phi(lam, phi, alpha, lam_n, phi_n)

    real(kind=8), intent(in)  :: lam,   phi, alpha
    real(kind=8), intent(out) :: lam_n, phi_n

    real(kind=8) :: x, y, z, xn, yn, zn

    call sph2cart(lam, phi, x, y, z)
    call rotate_3D_y(x, y, z, -alpha, xn, yn, zn)
    call cart2sph(xn, yn, zn, lam_n, phi_n)

end subroutine sph_rotate_phi

end module sph_coords_mod

module barotropic_instability_u_mod

implicit none

contains

function barotropic_instability_u(u0, phi) result(u)
    use const_mod, only : pi

    real(kind=8), intent(in) :: u0, phi
    real(kind=8)             :: u

    real(kind=8), parameter :: phi0 = pi/7d0, phi1 = .5d0*pi-phi0

    u = u0/exp(-4._8/(phi1-phi0)**2) * exp(1._8/min((phi-phi0)*(phi-phi1),-1e-2))

end function barotropic_instability_u

end module barotropic_instability_u_mod

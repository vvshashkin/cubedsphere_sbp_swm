module test_vertical_profiles_mod

use abstract_vertical_profile_mod, only : vertical_profile_t
use const_N_profile_mod,           only : const_N_profile_t
use isothermal_profile_mod,        only : isothermal_profile_t

implicit none

real(kind=8), parameter :: tolerance1 = 1e-14
real(kind=8), parameter :: tolerance2 = 1e-8

contains

subroutine test_vertical_profiles()

    call test_N_const()
    call test_isothermal()

end subroutine

subroutine test_N_const
    real(kind=8), parameter :: N=0.01, p0 = 101325._8, t0=297._8, grav = 10.0_8 !test if works with non-standard values
    real(kind=8), parameter :: h(4) = [0.0_8, 1e3_8, 10e3_8, 30e3_8]
    real(kind=8) :: theta(size(h)), d_theta(size(h))
    real(kind=8) err

    class(vertical_profile_t), allocatable :: profile


    profile = const_N_profile_t(N=N,grav=grav)

    err = check_surface(profile)
    if(err > tolerance1) then
        print *, "const_N profile surface checks failed", err
        stop
    end if
    theta   = profile%calc_theta(h,t0,p0)
    d_theta = profile%calc_dtheta_dz(h,t0,p0)
    err = maxval(abs(grav*d_theta/theta-N**2))
    if(err > tolerance1) then
        print *, "const_N profile N2 check failed", err
        stop
    end if
    err = check_hydristatic_equilibrium(profile,grav,h)
    if(err > tolerance1) then
        print *, "const_N profile hydrostatic equilibrium check failed", err
        stop
    end if
    print *, "const_N_profile test passed"
end subroutine test_N_const

subroutine test_isothermal
    real(kind=8), parameter :: p0 = 101325._8, t0=297._8, grav = 10.0_8 !test if works with non-standard values
    real(kind=8), parameter :: h(4) = [0.0_8, 1e3_8, 10e3_8, 30e3_8]
    real(kind=8) :: theta(size(h)), d_theta(size(h))
    real(kind=8) err

    class(vertical_profile_t), allocatable :: profile

    profile = isothermal_profile_t(grav=grav,p0=101325.0_8)

    err = check_surface(profile)
    if(err > tolerance1) then
        print *, "isothermal profile surface checks failed", err
        stop
    end if

    err = check_hydristatic_equilibrium(profile,grav,h)
    if(err > tolerance1) then
        print *, "isothermal profile hydrostatic equilibrium check failed", err
        stop
    end if
    print *, "isothermal_profile test passed"
end subroutine test_isothermal

function check_surface(profile) result(err)
    class(vertical_profile_t), intent(in) :: profile
    real(kind=8) :: err

    real(kind=8), parameter :: p0 = 101325._8, t0=297._8

    err = abs(t0 - profile%calc_temp(0.0_8,t0,p0))     + &
          abs(t0 - profile%calc_theta(0.0_8,t0,p0))    + &
          abs(p0 - profile%calc_pressure(0.0_8,t0,p0)) + &
          abs(1.0_8 - profile%calc_ExnerP(0.0_8,t0,p0))
end function check_surface

function check_hydristatic_equilibrium(profile, grav, h) result(err)
    use const_mod, only : Cp, rgaz

    class(vertical_profile_t), intent(in) :: profile
    real(kind=8), intent(in) :: grav, h(:)
    real(kind=8) :: err

    real(kind=8), parameter :: p0 = 101325._8, t0=297._8
    real(kind=8), parameter :: dh = 0.5_8

    real(kind=8) :: theta(size(h)), d_P(size(h)), pres(size(h)), temp(size(h))

    theta = profile%calc_theta(h, t0, p0)
    d_P = profile%calc_dExnerP_dz(h, t0, p0)

    err = maxval(abs(Cp*theta*d_P+grav))

end function check_hydristatic_equilibrium

end module test_vertical_profiles_mod

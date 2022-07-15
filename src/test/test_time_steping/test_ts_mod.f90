module test_rk4_mod

implicit none

contains

subroutine test_rk4()
    use stvec_mod,               only: stvec_t
    use stvec_iomega_mod,        only: stvec_iomega_t, init_stvec_iomega
    use domain_mod,              only: domain_t
    use operator_iomega_mod,     only: operator_iomega_t, init_iomega_operator
    use timescheme_mod,          only: timescheme_t
    use timescheme_factory_mod,  only: create_timescheme
    !use explicit_Eul1_mod,       only: init_explicit_Eul1_ts
!    use rk4_mod,                  only: rk4_t, init_rk4
!    use exp_taylor_mod,           only: exp_taylor_t, init_exp_taylor
!    use exp_krylov_mod,           only: exp_krylov_t, init_exp_krylov
    use const_mod,                only: pi, Day24h_sec
    type(stvec_iomega_t) v1
    class(stvec_t), allocatable :: v2, v3, v4
    type(domain_t)  :: domain
    type(operator_iomega_t), allocatable :: oper
    class(timescheme_t), allocatable :: explicit_Eul1, rk4
!    type(parameters_iomega_t), allocatable :: model_params
!    type(explicit_Eul1_t) ts_exEul
!    type(rk4_t) ts_rk4
!    type(exp_taylor_t) ts_exp_taylor
!    class(timescheme_abstract_t), allocatable :: ts_exp_krylov
    integer, parameter :: N = 10
    real(kind=8),    parameter :: dt = 6._8*3600._8
    complex(kind=8)  omega(N), ftrue1(N), ftrue_rk4(N), ftrue_exp(N)
    real(kind=8),    parameter :: tolerance = 1e-12_8
    real(kind=8),    parameter :: tolerance_exp = 1e-10_8

    integer i

    do i = 1,N
        omega(i) = cmplx(1.0_8 / (Day24h_sec*i), 2._8*pi*i / Day24h_sec)
        !Analytical responces of Eul1 & RK4 scheme
        ftrue1(i) = 1._8+dt*omega(i)
        ftrue_rk4(i) = 1._8+dt*omega(i)+0.5_8*(dt*omega(i))**2      + &
                                              (dt*omega(i))**3/6._8 + &
                                              (dt*omega(i))**4/24._8
        ftrue_exp(i) = exp(dt*omega(i))
    end do

    call init_stvec_iomega(v1,N)
    call init_iomega_operator(oper,omega)

!    v1%f(1:N) = 1._8
!    call v1%copy_to(v2)
!    call v1%copy_to(v3)
!    call v3%create_similar(v2)
!
!    call oper%apply_to(v2,v3,domain)
!    call v3%update(-dt,v2,domain)
!    call oper%solve(v2,v3,dt,domain)
!
!    select type(v2)
!    class is (stvec_iomega_t)
!        !print *, v2%f(:)-v1%f(:)
!        !print *, v1%f(:)
!    end select

    v1%f(1:N) = 1._8
    call create_timescheme(explicit_Eul1,v1,"explicit_Eul1")
    call explicit_Eul1%step(v1, oper, domain, dt)
    call ifpassed(v1%f, ftrue1,    tolerance, "Explicit Eulerian")

    v1%f(1:N) = 1._8
    call create_timescheme(rk4,v1,"rk4")
    call rk4%step(v1, oper, domain, dt)
    call ifpassed(v1%f, ftrue_rk4,    tolerance, "RK4")

!
!    ts_rk4 = init_rk4(oper, v2)
!    call ts_rk4%step(v2, model_params, dt)
!
!    ts_exp_taylor = init_exp_taylor(oper, v3, tolerance)
!    call ts_exp_taylor%step(v3, model_params, dt)
!
!    call init_exp_krylov(ts_exp_krylov, oper, v3, 2*N)
!    call ts_exp_krylov%step(v4, model_params, dt)
!
!    call ifpassed(v2%f, ftrue_rk4, tolerance, "RK4")
!    call ifpassed(v3%f, ftrue_exp, tolerance_exp, "EXP_Taylor")
!    call ifpassed(v4%f, ftrue_exp, tolerance, "EXP_Krylov")
!
end subroutine test_rk4

subroutine ifpassed(f, ftrue, tolerance, scheme_name)
    complex(kind=8), intent(in) :: f(:), ftrue(:)
    real(kind=8),    intent(in) :: tolerance
    character(*),    intent(in) :: scheme_name

    real(kind=8) err
    integer(kind=4) N

    N = size(f)

    err = sum( abs( f(1:N)-ftrue(1:N)))
    if(err < tolerance) then
        print *, scheme_name //" scheme test passed"
    else
        print *, scheme_name //" scheme test failed, error = ", err
    end if
end subroutine ifpassed

end module test_rk4_mod

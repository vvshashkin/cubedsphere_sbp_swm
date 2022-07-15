module test_metric_class_mod

implicit none

private
public :: test_metric_class

real(kind=8), parameter :: test_tolerace = 8e-15_8

contains

subroutine test_metric_class(topology_type,metric_type)
    use metric_mod,           only : metric_t
    use metric_factory_mod,   only : create_metric
    use config_metric_mod,    only : config_metric_t
    use topology_mod,         only : topology_t
    use topology_factory_mod, only : create_topology

    character(len=*), intent(in) :: topology_type, metric_type

    class(topology_t), allocatable :: topology
    class(metric_t),   allocatable :: metric
    real(kind=8) a1(4), a2(4), a3(4)
    real(kind=8) b1(3), b2(3), b3(4)
    real(kind=8) r(3)
    real(kind=8) Q(6), QI(6), Jac
    real(kind=8), parameter :: h_surf = 1000.0_8, h_top = 30000.0_8
    real(kind=8), parameter :: dhdalpha = 0.5_8, dhdbeta = -0.2_8
    real(kind=8), parameter :: stepx = 0.1
    real(kind=8), parameter :: stepy = 0.11
    real(kind=8), parameter :: stepz = 0.2
    integer(kind=4), parameter :: Nx = 1/stepx
    integer(kind=4), parameter :: Ny = 1/stepy
    integer(kind=4), parameter :: Nz = 1/stepz
    integer(kind=4) npanels, panel_ind, i, j, k
    real(kind=8) alpha, beta, eta, alpha0, beta0, da, db, detQ, h
    logical :: is_correct
    type(config_metric_t) :: config_metric

    topology = create_topology(topology_type)

    call config_metric%set_defaults()
    config_metric%vertical_scale = h_top
    call create_metric(metric, topology, metric_type,config_metric)

    npanels = topology%npanels
    alpha0 = metric%alpha0
    beta0 = metric%beta0
    da = metric%alpha1-metric%alpha0
    db = metric%beta1-metric%beta0

    is_correct = .true.

    do panel_ind = 1, npanels
        do k=0,Nz
            eta = k*stepz
            do j=0,Nx
                do i = 0,Ny
                    alpha = alpha0+i*stepx*da
                    beta = beta0+j*stepy*db
                    a1 = metric%calculate_a1(panel_ind, alpha, beta, eta, h_surf, dhdalpha, h_top)
                    a2 = metric%calculate_a2(panel_ind, alpha, beta, eta, h_surf, dhdbeta, h_top)
                    a3 = metric%calculate_a3(panel_ind, alpha, beta, eta, h_surf, h_top)
                    b1 = metric%calculate_b1(panel_ind, alpha, beta, eta, h_surf, dhdalpha, dhdbeta, h_top)
                    b2 = metric%calculate_b2(panel_ind, alpha, beta, eta, h_surf, dhdalpha, dhdbeta, h_top)
                    b3 = metric%calculate_b3(panel_ind, alpha, beta, eta, h_surf, dhdalpha, dhdbeta, h_top)
                    Q  = metric%calculate_Q(panel_ind,alpha,beta, eta, h_surf, dhdalpha, dhdbeta, h_top)
                    QI = metric%calculate_QI(panel_ind,alpha,beta, eta, h_surf, dhdalpha, dhdbeta, h_top)
                    Jac  = metric%calculate_J(panel_ind,alpha,beta, eta, h_surf, h_top)
                    call check("a1 is not normal to b2",sum(a1(1:3)*b2(1:3)), 0.0_8, test_tolerace, is_correct)
                    call check("a2 is not normal to b1",sum(a2(1:3)*b1(1:3)), 0.0_8, test_tolerace, is_correct)
                    call check("Q(1) /= a1*a1",sum(a1(1:4)*a1(1:4)), Q(1), test_tolerace, is_correct)
                    call check("Q(2) /= a1*a2",sum(a1(1:4)*a2(1:4)), Q(2), test_tolerace, is_correct)
                    call check("Q(3) /= a2*a2",sum(a2(1:4)*a2(1:4)), Q(3), test_tolerace, is_correct)
                    call check("Q(4) /= a1*a3",sum(a1(1:4)*a3(1:4)), Q(4), test_tolerace, is_correct)
                    call check("Q(5) /= a2*a3",sum(a2(1:4)*a3(1:4)), Q(5), test_tolerace, is_correct)
                    call check("Q(6) /= a3*a3",sum(a3(1:4)*a3(1:4)), Q(6), test_tolerace, is_correct)

                    call check("QI(1) /= b1*b1",sum(b1(1:3)*b1(1:3)), QI(1), test_tolerace, is_correct)
                    call check("QI(2) /= b1*b2",sum(b1(1:3)*b2(1:3)), QI(2), test_tolerace, is_correct)
                    call check("QI(3) /= b2*b2",sum(b2(1:3)*b2(1:3)), QI(3), test_tolerace, is_correct)
                    call check("QI(4) /= b1*b3",sum(b1(1:3)*b3(1:3)), QI(4), test_tolerace, is_correct)
                    call check("QI(5) /= b2*b3",sum(b2(1:3)*b3(1:3)), QI(5), test_tolerace, is_correct)
                    call check("QI(6) /= b3*b3",sum(b3(1:4)*b3(1:4)), QI(6), test_tolerace, is_correct)

                    detQ = Q(1)*Q(3)*Q(6)+2.0_8*Q(2)*Q(5)*Q(4)-Q(4)**2*Q(3)-Q(1)*Q(5)**2-Q(2)**2*Q(6)
                    call check("J**2 /= det(Q)",Jac**2,detQ,test_tolerace, is_correct)
                    detQ = Qi(1)*Qi(3)*Qi(6)+2.0_8*Qi(2)*Qi(5)*Qi(4)-Qi(4)**2*Qi(3)-Qi(1)*Qi(5)**2-Qi(2)**2*Qi(6)
                    call check("J**(-2) /= det(Qi)",1.0_8/Jac**2,detQ,test_tolerace, is_correct)


                    if(metric_type == "ecs") then
                        r = metric%calculate_r(panel_ind,alpha,beta,eta, h_surf, h_top)
                        call check("a1 is not normal to r",sum(a1(1:3)*r(1:3)), 0.0_8, test_tolerace, is_correct)
                        call check("a2 is not normal to r",sum(a2(1:3)*r(1:3)), 0.0_8, test_tolerace, is_correct)
                        call check("b1 is not normal to r",sum(b1(1:3)*r(1:3)), 0.0_8, test_tolerace, is_correct)
                        call check("b2 is not normal to r",sum(b2(1:3)*r(1:3)), 0.0_8, test_tolerace, is_correct)
                    end if
                    if(metric_type == "shallow_atmosphere_metric" .and. k == 0) then
                        h = metric%calculate_h(panel_ind,alpha,beta,0.0_8,h_surf,h_top)
                        call check("h(eta=0) /= h_surf",h, h_surf, test_tolerace, is_correct)
                        h = metric%calculate_h(panel_ind,alpha,beta,1.0_8,h_surf,h_top)
                        call check("h(eta=1) /= h_top",h, h_top, test_tolerace, is_correct)
                    end if

                    if(.not. is_correct) then
                        print *, "metric_class test (" // metric_type//") failed"
                        return
                    end if
                end do
            end do
        end do
    end do
    print *, "metric_class test (" // metric_type//") passed"
end subroutine test_metric_class

subroutine check(message,a,b,tol, is_correct)
    !check if abs(a-b) < tol and prints message otherwise
    character(len=*), intent(in)    :: message
    real(kind=8),     intent(in)    :: a, b
    real(kind=8),     intent(in)    :: tol
    logical,          intent(inout) :: is_correct

    logical is_correct_local

    is_correct_local = abs(a-b) < tol
    !print *, message, is_correct, is_correct_local
    is_correct = is_correct .and. is_correct_local

    if(.not. is_correct_local) then
        print *, message,", error: ", a-b
    end if

end

end module test_metric_class_mod

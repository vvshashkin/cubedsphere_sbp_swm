module test_quadrature_mod

use parcomm_mod,             only : parcomm_global
use abstract_quadrature_mod, only : quadrature_t
use quadrature_factory_mod,  only : create_quadrature
use domain_mod,              only : domain_t
use domain_factory_mod,      only : create_domain
use grid_field_mod,          only : grid_field_t
use grid_field_factory_mod,  only : create_grid_field
use mesh_mod,                only : mesh_t

use const_mod, only : pi

implicit none

contains

subroutine test_quadrature(quadrature_name, max_order_exact, staggering, points_type)
    character(len=*), intent(in) :: quadrature_name
    integer(kind=4),  intent(in) :: max_order_exact
    character(len=*), intent(in) :: staggering, points_type

    class(quadrature_t), allocatable :: quadrature
    type(domain_t), target           :: domain
    type(grid_field_t)               :: test_fun
    type(mesh_t), pointer            :: mesh

    real(kind=8) :: exact_quadrature, num_quadrature
    real(kind=8), parameter :: tolerance = 1e-12
    integer(kind=4)    :: order, t, i, j, k, p
    integer, parameter :: N=20, nz = 1
    real(kind=8) :: a, b


    call create_domain(domain, "cube", staggering, N, nz)

    select case(points_type)
    case('o')
        mesh => domain%mesh_o
    case('x')
        mesh => domain%mesh_x
    case('y')
        mesh => domain%mesh_y
    case('xy')
        mesh => domain%mesh_xy
    case default
        call parcomm_global%abort("unknown points type:" // points_type)
    end select

    call create_quadrature(quadrature, quadrature_name, mesh)
    call create_grid_field(test_fun, 0, 0, mesh)

    do order = 0, max_order_exact
        exact_quadrature = 0.0_8
        do t=mesh%ts, mesh%te
            do k = mesh%tile(t)%ks, mesh%tile(t)%ke
                do j = mesh%tile(t)%js, mesh%tile(t)%je
                    b = mesh%tile(t)%get_beta(j)
                    do i = mesh%tile(t)%is, mesh%tile(t)%ie
                        a = mesh%tile(t)%get_alpha(i)
                        test_fun%tile(t)%p(i,j,k) = 1.0_8/mesh%tile(t)%J(i,j,k)*a**order*b**order
                    end do
                end do
            end do
        end do
        p = (order+1)
        exact_quadrature = 6._8*((0.25_8*pi)**p - (-0.25_8*pi)**p)**2 / p**2
        num_quadrature = quadrature%mass(test_fun,mesh,domain%parcomm)

        if(abs(exact_quadrature - num_quadrature) > tolerance) then
            print *, quadrature_name //" quadrature test failed at order", &
                     order, "exact integral=", exact_quadrature, &
                            "calculated=integral=", num_quadrature, &
                            exact_quadrature-num_quadrature
            return
        end if
    end do

    call parcomm_global%print(quadrature_name //" quadrature test passed")

end subroutine test_quadrature

end module test_quadrature_mod

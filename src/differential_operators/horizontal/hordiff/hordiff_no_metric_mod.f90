module hordiff_no_metric_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t
use mesh_mod,       only : mesh_t
use halo_mod,       only : halo_t

use abstract_hordiff_mod,   only : hordiff_operator_t
use abstract_laplace_mod,   only : laplace_operator_t

implicit none

type, public, extends(hordiff_operator_t) :: hordiff_no_metric_t
    integer(kind=4) :: diff_order
    real(kind=8)    :: diff_coeff

    class(laplace_operator_t), allocatable :: laplace_op
    type(grid_field_t) :: f_tend_inter

contains
    procedure, public :: calc_diff
end type hordiff_no_metric_t

contains

subroutine calc_diff(this, f_tend, f, mesh, domain)

    use vec_math_mod, only : divide_by_J_self

    class(hordiff_no_metric_t), intent(inout) :: this
    type(grid_field_t),         intent(inout) :: f_tend
    type(grid_field_t),         intent(inout) :: f
    type(mesh_t),               intent(in)    :: mesh
    type(domain_t),             intent(in)    :: domain

    integer(kind=4) :: p

    real(kind=8) :: coeff

    call this%f_tend_inter%assign(1.0_8,f,mesh)
    call this%laplace_op%calc_laplace(this%f_tend_inter,this%f_tend_inter,domain)

    do p= 2, this%diff_order
        call this%laplace_op%calc_laplace(this%f_tend_inter,this%f_tend_inter,domain)
    end do
    coeff = (-1.0_8)**(this%diff_order+1)*this%diff_coeff**(2*this%diff_order)
    call f_tend%assign(coeff, this%f_tend_inter, mesh)

    ! call divide_by_J_self(f_tend, mesh)

end subroutine calc_diff

end module hordiff_no_metric_mod

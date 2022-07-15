module hordiff_scalar_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t
use mesh_mod,       only : mesh_t
use halo_mod,       only : halo_t

use abstract_hordiff_mod,   only : hordiff_operator_t
use abstract_div_mod,       only : div_operator_t
use abstract_grad_mod,      only : grad_operator_t
use abstract_co2contra_mod, only : co2contra_operator_t
use abstract_laplace_mod,   only : laplace_operator_t

implicit none

type, public, extends(hordiff_operator_t) :: hordiff_scalar_t
    integer(kind=4) :: diff_order
    real(kind=8)    :: diff_coeff

    class(laplace_operator_t), allocatable :: laplace_op
    type(grid_field_t) :: f_tend_inter
    ! class(div_operator_t),       allocatable :: div_op
    ! class(grad_operator_t),      allocatable :: grad_op
    ! class(co2contra_operator_t), allocatable :: co2contra_op
    !
    ! class(halo_t), allocatable :: edge_sync
    !
    ! type(grid_field_t) :: div, gxt, gyt, gx, gy
contains
    procedure, public :: calc_diff
end type hordiff_scalar_t

contains

subroutine calc_diff(this, f_tend, f, mesh, domain)

    class(hordiff_scalar_t), intent(inout) :: this
    type(grid_field_t),      intent(inout) :: f_tend
    type(grid_field_t),      intent(inout) :: f
    type(mesh_t),            intent(in)    :: mesh
    type(domain_t),          intent(in)    :: domain

    integer(kind=4) :: p

    real(kind=8) :: coeff

    call this%f_tend_inter%assign(1.0_8,f,mesh)
    call this%laplace_op%calc_laplace(this%f_tend_inter,this%f_tend_inter,domain)

    do p= 2, this%diff_order
        call this%laplace_op%calc_laplace(this%f_tend_inter,this%f_tend_inter,domain)
    end do
    coeff = (-1.0_8)**(this%diff_order+1)*this%diff_coeff**(2*this%diff_order)
    call f_tend%assign(coeff, this%f_tend_inter, mesh)

    ! call this%grad_op%calc_grad(this%gx, this%gy, f, domain)
    !
    ! call this%co2contra_op%transform(this%gxt, this%gyt, this%gx, this%gy, domain)
    !
    ! call this%div_op%calc_div(this%div, this%gxt, this%gyt, domain)
    !
    ! do p = 2, this%diff_order
    !     call this%grad_op%calc_grad(this%gx, this%gy, this%div, domain)
    !
    !     call this%co2contra_op%transform(this%gxt, this%gyt, this%gx, this%gy, domain)
    !
    !     call this%div_op%calc_div(this%div, this%gxt, this%gyt, domain)
    ! end do
    !
    ! coeff = (-1.0_8)**(this%diff_order+1)*this%diff_coeff**(2*this%diff_order)
    !
    ! if(allocated(this%edge_sync)) then
    !     call this%edge_sync%get_halo_scalar(this%div, domain, 1)
    ! end if
    !
    ! call f_tend%assign(coeff, this%div, mesh)

end subroutine calc_diff

end module hordiff_scalar_mod

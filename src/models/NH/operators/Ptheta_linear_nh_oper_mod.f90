module Ptheta_linear_nh_oper_mod

use operator_mod,                   only : operator_t
use grid_field_mod,                 only : grid_field_t
use abstract_grad_3d_mod,           only : grad_3d_operator_t
use abstract_div_3d_mod,            only : div_3d_operator_t
use abstract_co2contra_mod,         only : co2contra_operator_t
use abstract_interpolators3d_mod,   only : interpolator_w2uv_t
use abstract_vertical_operator_mod, only : vertical_operator_t
use stvec_mod,                      only : stvec_t
use stvec_nh_mod,                   only : stvec_nh_t
use parcomm_mod,                    only : parcomm_global
use domain_mod,                     only : domain_t

implicit none

type, extends(operator_t) :: Ptheta_linear_nh_operator_t

    class(grad_3d_operator_t),   allocatable :: grad_op
    class(div_3d_operator_t),    allocatable :: div_op
    class(co2contra_operator_t), allocatable :: co2contra_op
    class(interpolator_w2uv_t),  allocatable :: w2uv_op
    class(vertical_operator_t),  allocatable :: w2p_op

    type(grid_field_t) :: P0, dP0_dz, dP0_dz_w
    type(grid_field_t) :: theta0, theta0_u, theta0_v, dtheta0_dz
    type(grid_field_t) :: w, wp
    type(grid_field_t) :: grad_x, grad_y, grad_z
    type(grid_field_t) :: div3
    !type(grid_field_t) :: grad_x_contra, grad_y_contra

    contains
    procedure, public :: apply

end type Ptheta_linear_nh_operator_t

contains

subroutine apply(this, vout, vin, domain)

    use const_mod, only : rgaz, Cv, Cp

    class(Ptheta_linear_nh_operator_t), intent(inout) :: this
    class(stvec_t),                     intent(inout) :: vout !inout to enable preallocated vectors
    class(stvec_t),                     intent(inout) :: vin
    type(domain_t),                     intent(in)    :: domain

    select type(vout)
    class is (stvec_nh_t)
    select type(vin)
    class is (stvec_nh_t)
        ! call vout%u%assign(0.0_8, domain%mesh_u)
        ! call vout%v%assign(0.0_8, domain%mesh_v)
        ! call vout%eta_dot%assign(0.0_8, domain%mesh_w)
        ! call vout%theta%assign(0.0_8, domain%mesh_w)
        ! call vout%P%assign(0.0_8,domain%mesh_p)

        call this%grad_op%calc_grad(this%grad_x, this%grad_y, this%grad_z,vin%P,domain)
        call this%grad_x%assign_prod(-Cp,this%grad_x,this%theta0_u,domain%mesh_u)
        call this%grad_y%assign_prod(-Cp,this%grad_y,this%theta0_v,domain%mesh_v)
        call this%grad_z%assign_prod(-Cp,this%grad_z,this%theta0,domain%mesh_w)
        call this%co2contra_op%transform(vout%u,vout%v,this%grad_x, this%grad_y, domain)

        call vout%eta_dot%assign_prod(-Cp,vin%theta,this%dP0_dz_w,domain%mesh_w)
        call vout%eta_dot%update(1.0_8, this%grad_z, domain%mesh_w)

        call vout%theta%assign_prod(-1.0_8,this%dtheta0_dz,vin%eta_dot,domain%mesh_w)

        call this%w2p_op%apply(this%wp,vin%eta_dot,domain)
        call vout%P%assign_prod(-1.0_8, this%wp, this%dP0_dz,domain%mesh_p)
        call this%div_op%calc_div(this%div3,vin%u,vin%v,vin%eta_dot,domain)
        call this%div3%assign_prod(-rgaz/Cv,this%P0,this%div3,domain%mesh_p)
        call vout%P%update(1.0_8,this%div3,domain%mesh_p)

        vout%model_time = 1.0_8 != dt / dt

        call apply_0boundary_conds(vout%eta_dot,domain%mesh_w)
    class default
        call parcomm_global%abort("Ptheta_linear_nh_oper_t, vin type error")
    end select
    class default
        call parcomm_global%abort("Ptheta_linear_nh_oper_t, vout type error")
    end select


end subroutine apply

subroutine apply_0boundary_conds(eta_dot, mesh)

    use mesh_mod, only : mesh_t

    type(grid_field_t), intent(inout) :: eta_dot
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t, i, j, ks, ke, is,ie, js, je

    do t = mesh%ts, mesh%te
        is = mesh%tile(t)%is; ie = mesh%tile(t)%ie
        js = mesh%tile(t)%js; je = mesh%tile(t)%je
        ks = mesh%tile(t)%ks; ke = mesh%tile(t)%ke

        do j = js, je
            do i=is,ie
                eta_dot%tile(t)%p(i,j,ks) = 0.0_8
                eta_dot%tile(t)%p(i,j,ke) = 0.0_8
            end do
        end do
    end do
end subroutine apply_0boundary_conds

end module Ptheta_linear_nh_oper_mod

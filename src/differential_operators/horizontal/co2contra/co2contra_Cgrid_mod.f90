module co2contra_Cgrid_mod

use abstract_co2contra_mod,          only : co2contra_operator_t
use grid_field_mod,                  only : grid_field_t, tile_field_t
use mesh_mod,                        only : tile_mesh_t
use domain_mod,                      only : domain_t
use exchange_abstract_mod,           only : exchange_t
use parcomm_mod,                     only : parcomm_global
use halo_mod,                        only : halo_vec_t
use sbp_operator_mod,                only : sbp_operator_t
use abstract_interpolators2d_mod,    only : interpolator2d_vec2vec_t

implicit none

type, extends(co2contra_operator_t), public :: co2contra_c_sbp_t
    character(len=:),      allocatable :: operator_name
    class(exchange_t),     allocatable :: exchange_inner
    class(sbp_operator_t), allocatable :: sbp_interp_v2h, sbp_interp_h2v
contains
    procedure :: transform => transform_co2contra_c_sbp
end type co2contra_c_sbp_t

type, extends(co2contra_operator_t), public :: co2contra_c_sbp_new_t
    character(len=:), allocatable :: operator_name
    class(interpolator2d_vec2vec_t), allocatable  :: interp_h2v_op, interp_v2h_op
    type(grid_field_t)                            :: uh, vh
    contains
        procedure :: transform    => transform_co2contra_c_sbp_new
        procedure :: transform2co => transform_contra2co_c_sbp_new
end type co2contra_c_sbp_new_t

contains

subroutine transform_co2contra_c_sbp(this, u_contra, v_contra, u_cov, v_cov, domain)
    class(co2contra_c_sbp_t), intent(inout) :: this
    type(domain_t),           intent(in)    :: domain
    type(grid_field_t),       intent(inout) :: u_cov, v_cov
    !output:
    type(grid_field_t),       intent(inout) :: u_contra, v_contra

    integer(kind=4) :: t

    call this%exchange_inner%do_vec(u_cov, v_cov, domain%parcomm)

    select case(this%operator_name)
    case ("co2contra_c_sbp21")
        do t=domain%mesh_o%ts, domain%mesh_o%te
            call transform_co2contra_c_sbp21_tile(u_contra%tile(t), v_contra%tile(t),&
                                                  u_cov%tile(t), v_cov%tile(t),      &
                                                  domain%mesh_u%tile(t), domain%mesh_v%tile(t), &
                                                  domain%mesh_o%tile(t))
        end do
    case("co2contra_c_sbp42")
        do t=domain%mesh_o%ts, domain%mesh_o%te
            call transform_co2contra_c_sbp42_tile(u_contra%tile(t), v_contra%tile(t),           &
                                                  u_cov%tile(t), v_cov%tile(t),                 &
                                                  this%sbp_interp_h2v, this%sbp_interp_v2h,     &
                                                  domain%mesh_u%tile(t), domain%mesh_v%tile(t), &
                                                  domain%mesh_o%tile(t))
        end do
    case default
        call parcomm_global%abort("unknown co2contra_c_sbp operator "// this%operator_name)
    end select
end subroutine transform_co2contra_c_sbp

subroutine transform_co2contra_c_sbp_new(this, u_contra, v_contra, u_cov, v_cov, domain)
    class(co2contra_c_sbp_new_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u_cov, v_cov
    !output:
    type(grid_field_t),           intent(inout) :: u_contra, v_contra

    integer(kind=4) :: t

    call this%interp_v2h_op%interp2d_vec2vec(this%uh, this%vh, u_cov, v_cov, domain)

    do t = domain%partition%ts, domain%partition%te
        call start_co2contra_transform_at_h_tile(this%uh%tile(t), this%vh%tile(t), &
                                       domain%mesh_p%tile(t))
    end do

    !u_contra contains v part at u, v_contra contains u part at v
    call this%interp_h2v_op%interp2d_vec2vec(u_contra, v_contra, this%uh, this%vh, domain)

    do t = domain%partition%ts, domain%partition%te
        call finalize_co2contra_transform_at_uv_tile(u_contra%tile(t), v_contra%tile(t), &
                   u_cov%tile(t), v_cov%tile(t), &
                domain%mesh_u%tile(t), domain%mesh_v%tile(t))
    end do

end subroutine transform_co2contra_c_sbp_new
subroutine start_co2contra_transform_at_h_tile(uh, vh, mesh_o)

    type(tile_field_t),     intent(inout) :: uh, vh
    type(tile_mesh_t),      intent(in)    :: mesh_o

    real(kind=8)    :: u, v
    integer(kind=4) :: i, j, k

    do k = mesh_o%ks,mesh_o%ke
        do j = mesh_o%js, mesh_o%je
            do i = mesh_o%is, mesh_o%ie
                u = uh%p(i,j,k)*mesh_o%J(i,j,k)*mesh_o%Qi(2,i,j,k)
                v = vh%p(i,j,k)*mesh_o%J(i,j,k)*mesh_o%Qi(2,i,j,k)
                !Change components for (u,v) to (v,u) interpolation
                uh%p(i,j,k) = v
                vh%p(i,j,k) = u
            end do
        end do
    end do
end subroutine start_co2contra_transform_at_h_tile
subroutine finalize_co2contra_transform_at_uv_tile(u_contra, v_contra, u_cov, v_cov, mesh_u, mesh_v)

    type(tile_field_t),     intent(inout) :: u_contra, v_contra, u_cov, v_cov
    type(tile_mesh_t),      intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k

    do k = mesh_u%ks,mesh_u%ke
        do j = mesh_u%js, mesh_u%je
            do i = mesh_u%is, mesh_u%ie
                u_contra%p(i,j,k) = mesh_u%Qi(1,i,j,k)*u_cov%p(i,j,k)+u_contra%p(i,j,k)/mesh_u%J(i,j,k)
            end do
        end do
    end do
    do k = mesh_v%ks,mesh_v%ke
        do j = mesh_v%js, mesh_v%je
            do i = mesh_v%is, mesh_v%ie
                v_contra%p(i,j,k) = mesh_v%Qi(3,i,j,k)*v_cov%p(i,j,k)+v_contra%p(i,j,k)/mesh_v%J(i,j,k)
            end do
        end do
    end do

end subroutine finalize_co2contra_transform_at_uv_tile
subroutine transform_contra2co_c_sbp_new(this, u_cov, v_cov, u_contra, v_contra, domain)
    class(co2contra_c_sbp_new_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u_contra, v_contra
    !output:
    type(grid_field_t),           intent(inout) :: u_cov, v_cov

    integer(kind=4) :: t

    call this%interp_v2h_op%interp2d_vec2vec(this%uh, this%vh, u_contra, v_contra, domain)

    do t = domain%partition%ts, domain%partition%te
        call start_contra2co_transform_at_h_tile(this%uh%tile(t), this%vh%tile(t),&
                                                 domain%mesh_p%tile(t))
    end do

    !u_cov contains v part at u, v_cov contains u part at v
    call this%interp_h2v_op%interp2d_vec2vec(u_cov, v_cov, this%uh, this%vh, domain)

    do t = domain%partition%ts, domain%partition%te
        call finalize_contra2co_transform_at_uv_tile(u_cov%tile(t), v_cov%tile(t), &
                   u_contra%tile(t), v_contra%tile(t), &
                   domain%mesh_u%tile(t), domain%mesh_v%tile(t))
    end do

end subroutine transform_contra2co_c_sbp_new
subroutine start_contra2co_transform_at_h_tile(uh, vh, mesh_o)

    type(tile_field_t),     intent(inout) :: uh, vh
    type(tile_mesh_t),      intent(in)    :: mesh_o

    real(kind=8)    :: u, v
    integer(kind=4) :: i, j, k

    do k = mesh_o%ks,mesh_o%ke
        do j = mesh_o%js, mesh_o%je
            do i = mesh_o%is, mesh_o%ie
                u = uh%p(i,j,k)*mesh_o%J(i,j,k)*mesh_o%Q(2,i,j,k)
                v = vh%p(i,j,k)*mesh_o%J(i,j,k)*mesh_o%Q(2,i,j,k)
                !Change components for (u,v) to (v,u) interpolation
                uh%p(i,j,k) = v
                vh%p(i,j,k) = u
            end do
        end do
    end do
end subroutine start_contra2co_transform_at_h_tile
subroutine finalize_contra2co_transform_at_uv_tile(u_cov, v_cov, u_contra, v_contra, mesh_u, mesh_v)

    type(tile_field_t),     intent(inout) :: u_contra, v_contra, u_cov, v_cov
    type(tile_mesh_t),      intent(in)    :: mesh_u, mesh_v

    integer(kind=4) :: i, j, k

    do k = mesh_u%ks,mesh_u%ke
        do j = mesh_u%js, mesh_u%je
            do i = mesh_u%is, mesh_u%ie
                u_cov%p(i,j,k) = mesh_u%Q(1,i,j,k)*u_contra%p(i,j,k)+u_cov%p(i,j,k)/mesh_u%J(i,j,k)
            end do
        end do
    end do
    do k = mesh_v%ks,mesh_v%ke
        do j = mesh_v%js, mesh_v%je
            do i = mesh_v%is, mesh_v%ie
                v_cov%p(i,j,k) = mesh_v%Q(3,i,j,k)*v_contra%p(i,j,k)+v_cov%p(i,j,k)/mesh_v%J(i,j,k)
            end do
        end do
    end do

end subroutine finalize_contra2co_transform_at_uv_tile
! subroutine transform_co2contra_c_halo(this, u_contra, v_contra, u_cov, v_cov, domain)
!     class(co2contra_c_halo_t), intent(inout) :: this
!     type(domain_t),            intent(in)    :: domain
!     type(grid_field_t),        intent(inout) :: u_cov, v_cov
!     !output:
!     type(grid_field_t),        intent(inout) :: u_contra, v_contra
!
!     integer(kind=4) :: t
!
!   call parcomm_global%abort("this subroutine is not working, no C-halo procedure for covariant vectors")
!
!     call this%halo%get_halo_vector(u_cov, v_cov, domain, 1)
!
!     do t=domain%mesh_o%ts, domain%mesh_o%te
!         call transform_co2contra_c_halo_tile(u_contra%tile(t), v_contra%tile(t),&
!                                               u_cov%tile(t), v_cov%tile(t),      &
!                                               domain%mesh_u%tile(t), domain%mesh_v%tile(t), &
!                                               domain%mesh_o%tile(t))
!     end do
!
! end subroutine transform_co2contra_c_halo

subroutine transform_co2contra_c_sbp21_tile(u_contra, v_contra, u_cov, v_cov, mesh_u, mesh_v, mesh_o)

    type(tile_field_t), intent(in)    :: u_cov, v_cov
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v, mesh_o
    !output
    type(tile_field_t), intent(inout) :: u_contra, v_contra

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke
    real(kind=8)    :: u_at_o(mesh_u%is:mesh_u%ie,mesh_u%js-1:mesh_u%je)
    real(kind=8)    :: v_at_o(mesh_u%is-1:mesh_u%ie,mesh_u%js:mesh_u%je)
    real(kind=8)    :: v_at_u, u_at_v

    ks = mesh_o%ks; ke = mesh_o%ke

    do k=ks, ke

        is = mesh_u%is; ie = mesh_u%ie
        js = mesh_u%js; je = mesh_u%je

        do j=js,je
            do i=max(is-1,1),min(ie,mesh_v%nx)
                v_at_o(i,j) = 0.5_8*(v_cov%p(i,j,k)+v_cov%p(i,j+1,k))*mesh_o%J(i,j,k)*mesh_o%Qi(2,i,j,k)
            end do

            if(is == 1) then
                u_contra%p(1,j,k) = mesh_u%Qi(1,1,j,k)*u_cov%p(1,j,k)+v_at_o(1,j) / mesh_u%J(1,j,k)
            end if
            do i = max(is,2), min(ie,mesh_v%nx)
                v_at_u = 0.5_8*(v_at_o(i,j)+v_at_o(i-1,j))
                u_contra%p(i,j,k) = mesh_u%Qi(1,i,j,k)*u_cov%p(i,j,k)+v_at_u / mesh_u%J(i,j,k)
            end do
            if(ie == mesh_u%nx) then
                u_contra%p(ie,j,k) = mesh_u%Qi(1,ie,j,k)*u_cov%p(ie,j,k)+v_at_o(ie-1,j) / mesh_u%J(ie,j,k)
            end if
        end do

        is = mesh_v%is; ie = mesh_v%ie
        js = mesh_v%js; je = mesh_v%je

        do j=max(js-1,1),min(je,mesh_u%ny)
            do i = is,ie
                u_at_o(i,j) = 0.5_8*(u_cov%p(i+1,j,k)+u_cov%p(i,j,k))*mesh_o%J(i,j,k)*mesh_o%Qi(2,i,j,k)
            end do
        end do

        if(js == 1) then
            do i=is, ie
                v_contra%p(i,1,k) = mesh_v%Qi(3,i,1,k)*v_cov%p(i,1,k)+u_at_o(i,1) / mesh_v%J(i,1,k)
            end do
        end if
        do j = max(js,2), min(je,mesh_u%ny)
            do i = is, ie
                u_at_v = 0.5_8*(u_at_o(i,j)+u_at_o(i,j-1))
                v_contra%p(i,j,k) = mesh_v%Qi(3,i,j,k)*v_cov%p(i,j,k)+u_at_v / mesh_v%J(i,j,k)
            end do
        end do
        if(je == mesh_v%ny) then
            do i = is, ie
                v_contra%p(i,je,k) = mesh_v%Qi(3,i,je,k)*v_cov%p(i,je,k)+u_at_o(i,je-1) / mesh_v%J(i,je,k)
            end do
        end if

    end do

end subroutine transform_co2contra_c_sbp21_tile

subroutine transform_co2contra_c_sbp42_tile(gx, gy, dx, dy,                 &
                                            sbp_interp_h2v, sbp_interp_v2h, &
                                            mesh_x, mesh_y, mesh_o)

    use tile_mod, only : tile_t

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(inout) :: dx, dy
    class(sbp_operator_t),  intent(in)    :: sbp_interp_v2h, sbp_interp_h2v
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_o

    integer(kind=4) :: i, j, k, l, nsrc, ntarg
    integer(kind=4) :: is, ie, js, je, ks, ke
    !indices for grad components:
    integer(kind=4) :: isx, iex, jsx, jex
    integer(kind=4) :: isy, iey, jsy, jey
    !indices for grad components interpolated to p-points:
    integer(kind=4) :: ispx, iepx, jspx, jepx
    integer(kind=4) :: ispy, iepy, jspy, jepy

    integer(kind=4), parameter :: hw = 3 !maximum of interpolation stencil displacement
    real(kind=8)    :: dx_at_p(mesh_y%is:mesh_y%ie,mesh_y%js-hw:mesh_y%je+hw,1)
    real(kind=8)    :: dy_at_p(mesh_x%is-hw:mesh_x%ie+hw,mesh_x%js:mesh_x%je,1)

    type(tile_t) :: dy_at_p_bounds, dy_at_p_work, gx_work, gx_bounds
    type(tile_t) :: dx_at_p_bounds, dx_at_p_work, gy_work, gy_bounds

    ks = mesh_o%ks; ke = mesh_o%ke

    !bounding box for required values of gx
    isx = mesh_x%is; iex = mesh_x%ie
    jsx = mesh_x%js; jex = mesh_x%je
    !bounding box for needed values of dy_at_p:
    ispy = interp_p_stencil_start(isx,mesh_o%nx)
    iepy = interp_p_stencil_end  (iex,mesh_o%nx)
    jspy = jsx
    jepy = jex

    !bounding box for required values of gy
    isy = mesh_y%is; iey = mesh_y%ie
    jsy = mesh_y%js; jey = mesh_y%je
    !bounding box for needed values of dx_at_p:
    ispx = isy
    iepx = iey
    jspx = interp_p_stencil_start(jsx,mesh_o%ny)
    jepx = interp_p_stencil_end  (jex,mesh_o%ny)

    dy_at_p_bounds = tile_t(is = isx-hw,ie = iex+hw, js = jsx,  je = jex,  ks=1, ke=1)
    dy_at_p_work   = tile_t(is = ispy,  ie = iepy,   js = jspy, je = jepy, ks=1, ke=1)

    dx_at_p_bounds = tile_t(is = isy, ie = iey, js = jsx-hw, je = jex+hw,  ks=1, ke=1)
    dx_at_p_work   = tile_t(is = ispx,ie = iepx,js = jspx,   je = jepx,    ks=1, ke=1)

    gx_bounds = tile_t(is = gx%is, ie = gx%ie, js = gx%js, je = gx%je, ks = gx%ks, ke = gx%ke)
    gx_work   = tile_t(is = isx,   ie = iex,   js = jsx,   je = jex,   ks = 1,ke = 1)

    gy_bounds = tile_t(is = gy%is, ie = gy%ie, js = gy%js, je = gy%je, ks = gy%ks, ke = gy%ke)
    gy_work   = tile_t(is = isy,   ie = iey,   js = jsy,   je = jey,   ks = 1,ke = 1)

    do k = ks,ke

        !interpolate dy to p-points
        ntarg = mesh_o%ny
        nsrc  = mesh_o%ny+1

        dy_at_p_bounds%ks = k; dy_at_p_bounds%ke = k
        dy_at_p_work%ks   = k; dy_at_p_work%ke   = k
        call sbp_interp_v2h%apply(dy_at_p, dy_at_p_work, dy_at_p_bounds, mesh_o%ny, 'y', dy)

        !mutiplicate by metric terms in p-points
        do j = jspy,jepy
            do i = ispy, iepy
                dy_at_p(i,j,1) = mesh_o%J(i,j,k)*mesh_o%Qi(2,i,j,k)*dy_at_p(i,j,1)
            end do
        end do
        !interpolate to u-points
        gx_work%ks = k; gx_work%ke = k
        call sbp_interp_h2v%apply(gx%p, gx_work, gx_bounds, mesh_x%nx, 'x', dy_at_p, dy_at_p_bounds)

        dx_at_p_bounds%ks = k; dx_at_p_bounds%ke = k
        dx_at_p_work%ks   = k; dx_at_p_work%ke   = k
        call sbp_interp_v2h%apply(dx_at_p, dx_at_p_work, dx_at_p_bounds, mesh_o%nx, 'x', dx)

        !multiplicate by metric terms in p-points
        do j = jspx,jepx
            do i = ispx, iepx
                dx_at_p(i,j,1) = mesh_o%J(i,j,k)*mesh_o%Qi(2,i,j,k)*dx_at_p(i,j,1)
            end do
        end do

        !interpolate dx to v-points
        gy_work%ks = k; gy_work%ke = k
        call sbp_interp_h2v%apply(gy%p, gy_work, gy_bounds, mesh_y%ny, 'y', dx_at_p, dx_at_p_bounds)

        do j=jsx,jex
            do i=isx,iex
                gx%p(i,j,k) = mesh_x%Qi(1,i,j,k)*dx%p(i,j,k)+gx%p(i,j,k)/mesh_x%J(i,j,k)
            end do
        end do

        do j=jsy,jey
            do i=isy,iey
                gy%p(i,j,k) = mesh_y%Qi(3,i,j,k)*dy%p(i,j,k)+gy%p(i,j,k)/mesh_y%J(i,j,k)
            end do
        end do
    end do

    contains
    !from u/v points to p points
    integer(kind=4) function interp_v_stencil_start(i, n) result(is)
        integer(kind=4) :: i, n
        if(i <= 4) then
            is = 1
        else
            is = min(i-1,n-3)
        end if
    end

    integer(kind=4) function interp_v_stencil_end(i, n) result(ie)
        integer(kind=4) :: i, n
        if(i <= n-5) then
            ie = max(4,i+2)
        else
            ie = n
        end if
    end

    !from p-points to u/v points
    integer(kind=4) function interp_p_stencil_start(i, n) result(is)
        integer(kind=4) :: i, n
        if(i <= 4) then
            is = 1
        else
            is = min(i-2,n-3)
        end if
    end

    integer(kind=4) function interp_p_stencil_end(i, n) result(ie)
        integer(kind=4) :: i, n
        if(i <= n-3) then
            ie = max(4,i+1)
        else
            ie = n
        end if
    end
end subroutine transform_co2contra_c_sbp42_tile

! subroutine transform_co2contra_c_halo_tile(u_contra, v_contra, u_cov, v_cov, mesh_u, mesh_v, mesh_o)
!
!     type(tile_field_t), intent(in)    :: u_cov, v_cov
!     type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v, mesh_o
!     !output
!     type(tile_field_t), intent(inout) :: u_contra, v_contra
!
!     integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke
!     real(kind=8)    :: u_at_v, v_at_u
!
!     u_contra%p = 0.0
!     v_contra%p = 0.0
!
!     ks = mesh_o%ks; ke = mesh_o%ke
!
!     do k = ks, ke
!         is = mesh_u%is; ie = mesh_u%ie
!         js = mesh_u%js; je = mesh_u%je
!         do j = js, je
!             do i = is, ie
!                 v_at_u = 0.25_8*(v_cov%p(i,j,k)+v_cov%p(i-1,j,k)+v_cov%p(i,j+1,k)+v_cov%p(i-1,j+1,k))
!                 u_contra%p(i,j,k) = u_cov%p(i,j,k)!mesh_u%Qi(1,i,j)*u_cov%p(i,j,k)+mesh_u%Qi(2,i,j)*v_at_u
!             end do
!         end do
!
!         is = mesh_v%is; ie = mesh_v%ie
!         js = mesh_v%js; je = mesh_v%je
!         do j = js, je
!             do i = is, ie
!                 u_at_v = 0.25_8*(u_cov%p(i,j,k)+u_cov%p(i,j-1,k)+u_cov%p(i+1,j,k)+u_cov%p(i+1,j-1,k))
!                 v_contra%p(i,j,k) = v_cov%p(i,j,k)!mesh_v%Qi(3,i,j)*v_cov%p(i,j,k)+mesh_v%Qi(2,i,j)*u_at_v
!             end do
!         end do
!     end do
!
! end subroutine transform_co2contra_c_halo_tile

end module co2contra_Cgrid_mod

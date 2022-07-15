module massflux_Cgrid_mod

use abstract_massflux_mod, only : massflux_operator_t
use grid_field_mod,        only : grid_field_t
use domain_mod,            only : domain_t
use halo_mod,              only : halo_t, halo_vec_t
use parcomm_mod,           only : parcomm_global
use sbp_operator_mod,      only : sbp_operator_t

implicit none

type, extends(massflux_operator_t), public :: massflux_chalo_t

    integer(kind=4)                :: order
    class(halo_t),     allocatable :: halo
    class(halo_vec_t), allocatable :: halo_flux

    contains

    procedure :: calc_massflux => calc_c2_massflux

end type massflux_chalo_t

type, extends(massflux_operator_t), public :: massflux_c_up4_t

    class(halo_t),     allocatable :: halo
    class(halo_vec_t), allocatable :: halo_flux

    contains

    procedure :: calc_massflux => calc_up4_massflux

end type massflux_c_up4_t

type, extends(massflux_operator_t), public :: massflux_c_sbp21_t

    contains

    procedure :: calc_massflux => calc_c_sbp21_massflux

end type massflux_c_sbp21_t

type, extends(massflux_operator_t), public :: massflux_c_sbp42_t

    class(sbp_operator_t), allocatable :: sbp_interp_h2v

    contains

    procedure :: calc_massflux => calc_c_sbp42_massflux

end type massflux_c_sbp42_t

contains

subroutine calc_c2_massflux(this, fx, fy, f, u, v, domain)
    class(massflux_chalo_t), intent(inout) :: this
    type(domain_t),          intent(in)    :: domain
    type(grid_field_t),      intent(inout) :: f, u, v
    !output:
    type(grid_field_t),      intent(inout) :: fx, fy

    integer(kind=4) :: t


    select case(this%order)
    case(2)
        call this%halo%get_halo_scalar(f,domain,1)
        do t = domain%mesh_p%ts, domain%mesh_p%te
            call calc_c2_massflux_tile(fx%tile(t), fy%tile(t), &
                                       f%tile(t), u%tile(t), v%tile(t), &
                                       domain%mesh_x%tile(t), domain%mesh_y%tile(t))
        end do
    case(4)
        call this%halo%get_halo_scalar(f,domain,2)
        do t = domain%mesh_p%ts, domain%mesh_p%te
            call calc_c4_massflux_tile(fx%tile(t), fy%tile(t), &
                                       f%tile(t), u%tile(t), v%tile(t), &
                                       domain%mesh_x%tile(t), domain%mesh_y%tile(t))
        end do
    case default
        call parcomm_global%abort("massflux_chalo_t is currently implemented only for orders=2,4")
    end select
    call this%halo_flux%get_halo_vector(fx,fy,domain,1)

end subroutine calc_c2_massflux

subroutine calc_up4_massflux(this, fx, fy, f, u, v, domain)
    class(massflux_c_up4_t), intent(inout) :: this
    type(domain_t),          intent(in)    :: domain
    type(grid_field_t),      intent(inout) :: f, u, v
    !output:
    type(grid_field_t),      intent(inout) :: fx, fy

    integer(kind=4) :: t

    call this%halo%get_halo_scalar(f,domain,3)
    do t = domain%mesh_p%ts, domain%mesh_p%te
        call calc_up4_massflux_tile(fx%tile(t), fy%tile(t), &
                                    f%tile(t), u%tile(t), v%tile(t), &
                                    domain%mesh_x%tile(t), domain%mesh_y%tile(t))
    end do
    !call this%halo_flux%get_halo_vector(fx,fy,domain,1)

end subroutine calc_up4_massflux

subroutine calc_c2_massflux_tile(fx,fy,f,u,v,mesh_x,mesh_y)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_x, mesh_y

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke

    ks = mesh_x%ks; ke = mesh_x%ke

    do k=ks,ke
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je

        do j=js,je; do i=is,ie
            fx%p(i,j,k) = u%p(i,j,k)*0.5_8*(f%p(i,j,k)+f%p(i-1,j,k))
        end do; end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        do j=js,je; do i=is,ie
            fy%p(i,j,k) = v%p(i,j,k)*0.5_8*(f%p(i,j,k)+f%p(i,j-1,k))
        end do; end do

    end do
end subroutine calc_c2_massflux_tile

subroutine calc_c4_massflux_tile(fx,fy,f,u,v,mesh_x,mesh_y)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_x, mesh_y

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke

    ks = mesh_x%ks; ke = mesh_x%ke

    do k=ks,ke
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je

        do j=js,je; do i=is,ie
            fx%p(i,j,k) = u%p(i,j,k)*(-f%p(i-2,j,k)+7._8*f%p(i-1,j,k)+7._8*f%p(i,j,k)-f%p(i+1,j,k))/12._8
        end do; end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        do j=js,je; do i=is,ie
            fy%p(i,j,k) = v%p(i,j,k)*(-f%p(i,j-2,k)+7._8*f%p(i,j-1,k)+7._8*f%p(i,j,k)-f%p(i,j+1,k))/12._8
        end do; end do

    end do
end subroutine calc_c4_massflux_tile

subroutine calc_up4_massflux_tile(fx,fy,f,u,v,mesh_x,mesh_y)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_x, mesh_y

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke
    real(kind=8)    :: zl, zr

    ks = mesh_x%ks; ke = mesh_x%ke

    do k=ks,ke
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je

        do j=js,je; do i=is,ie
            zl = .5_8+sign(.5_8,u%p(i,j,k))
            zr = 1._8-zl
            fx%p(i,j,k) = u%p(i,j,k)*( &
                          zl*(5._8*f%p(i,j,k)+15._8*f%p(i-1,j,k)-&
                                              5._8*f%p(i-2,j,k)+f%p(i-3,j,k))+&
                          zr*(5._8*f%p(i-1,j,k)+15._8*f%p(i,j,k)-&
                                              5._8*f%p(i+1,j,k)+f%p(i+2,j,k)))/16._8
        end do; end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        do j=js,je; do i=is,ie
            zl = .5_8+sign(.5_8,v%p(i,j,k))
            zr = 1._8-zl
            fy%p(i,j,k) = v%p(i,j,k)*( &
                            zl*(5._8*f%p(i,j,k)+15._8*f%p(i,j-1,k)-&
                                              5._8*f%p(i,j-2,k)+f%p(i,j-3,k))+&
                            zr*(5._8*f%p(i,j-1,k)+15._8*f%p(i,j,k)-&
                                             5._8*f%p(i,j+1,k)+f%p(i,j+2,k)))/16._8
        end do; end do

    end do
end subroutine calc_up4_massflux_tile

subroutine calc_c_sbp21_massflux(this, fx, fy, f, u, v, domain)

    class(massflux_c_sbp21_t), intent(inout) :: this
    type(domain_t),            intent(in)    :: domain
    type(grid_field_t),        intent(inout) :: f, u, v
    !output:
    type(grid_field_t),        intent(inout) :: fx, fy

    integer(kind=4) :: t

    do t = domain%mesh_p%ts, domain%mesh_p%te

        call calc_c_sbp21_massflux_tile(fx%tile(t), fy%tile(t), &
                                        f%tile(t), u%tile(t), v%tile(t), &
                                        domain%mesh_x%tile(t), domain%mesh_y%tile(t),&
                                        domain%mesh_o%tile(t))

    end do

end subroutine calc_c_sbp21_massflux

subroutine calc_c_sbp21_massflux_tile(fx,fy,f,u,v,mesh_x,mesh_y, mesh_o)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_x, mesh_y, mesh_o

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke

    ks = mesh_x%ks; ke = mesh_x%ke

    do k=ks,ke
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je

        do j=js,je
            if(is == 1) fx%p(1,j,k) = u%p(1,j,k)*f%p(1,j,k)
            do i=max(is,2),min(ie,mesh_o%nx)
                fx%p(i,j,k) = u%p(i,j,k)*0.5_8*(f%p(i,j,k)+f%p(i-1,j,k))
            end do
            if(ie == mesh_o%nx+1) fx%p(mesh_o%nx+1,j,k) = u%p(mesh_o%nx+1,j,k)*f%p(mesh_o%nx,j,k)
        end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        if(js == 1) then
            do i=is,ie
                fy%p(i,1,k) = v%p(i,1,k)*f%p(i,1,k)
            end do
        end if

        do j=max(js,2),min(je,mesh_o%ny); do i=is,ie
            fy%p(i,j,k) = v%p(i,j,k)*0.5_8*(f%p(i,j,k)+f%p(i,j-1,k))
        end do; end do

        if(je == mesh_o%ny+1) then
            do i=is,ie
                fy%p(i,mesh_o%ny+1,k) = v%p(i,mesh_o%ny+1,k)*f%p(i,mesh_o%ny,k)
            end do
        end if

    end do
end subroutine calc_c_sbp21_massflux_tile

subroutine calc_c_sbp42_massflux(this, fx, fy, f, u, v, domain)

    class(massflux_c_sbp42_t), intent(inout) :: this
    type(domain_t),            intent(in)    :: domain
    type(grid_field_t),        intent(inout) :: f, u, v
    !output:
    type(grid_field_t),        intent(inout) :: fx, fy

    integer(kind=4) :: t

    do t = domain%mesh_p%ts, domain%mesh_p%te

        call calc_c_sbp42_massflux_tile(fx%tile(t), fy%tile(t), &
                                        f%tile(t), u%tile(t), v%tile(t), &
                                        this%sbp_interp_h2v, &
                                        domain%mesh_x%tile(t), domain%mesh_y%tile(t),&
                                        domain%mesh_o%tile(t))

    end do

end subroutine calc_c_sbp42_massflux

subroutine calc_c_sbp42_massflux_tile(fx,fy,f,u,v,sbp_interp_h2v,mesh_x,mesh_y, mesh_o)

    use tile_mod,   only : tile_t

    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t),    intent(inout) :: fx, fy
    type(tile_field_t),    intent(in)    :: f, u, v
    class(sbp_operator_t), intent(in)    :: sbp_interp_h2v
    type(tile_mesh_t),     intent(in)    :: mesh_x, mesh_y, mesh_o

    integer(kind=4) :: k, ks, ke
    integer(kind=4) :: i, j
    integer(kind=4) :: is, ie, js, je

    real(kind=8)    :: fp(f%is:f%ie,f%js:f%je,1)
    type(tile_t)    :: fp_tile, fx_tile, fx_work_tile,  fy_tile, fy_work_tile

    fp_tile = tile_t(is = f%is, ie = f%ie, js = f%js, je = f%je, ks = 1, ke = 1)
    fx_tile = tile_t(is = fx%is, ie = fx%ie, js = fx%js, je = fx%je, ks = fx%ks, ke = fx%ke)
    fx_work_tile = tile_t(is = mesh_x%is, ie = mesh_x%ie, js = mesh_x%js, je = mesh_x%je, &
                          ks = 1, ke = 1)
    fy_tile = tile_t(is = fy%is, ie = fy%ie, js = fy%js, je = fy%je, ks = fy%ks, ke = fy%ke)
    fy_work_tile = tile_t(is = mesh_y%is, ie = mesh_y%ie, js = mesh_y%js, je = mesh_y%je, &
                          ks = 1, ke = 1)

    do k = mesh_o%ks, mesh_o%ke

        is = mesh_o%is; ie = mesh_o%ie
        js = mesh_o%js; je = mesh_o%je

        do j=js,je
            do i=is,ie
                fp(i,j,1) = f%p(i,j,k)!mesh_o%G(i,j)*f%p(i,j,k)
            end do
        end do

        fp_tile%ks = k; fp_tile%ke = k
        call sbp_interp_h2v%apply(fx%p,fx_work_tile,fx_tile, mesh_o%nx+1, 'x', fp, fp_tile)

        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je
        do j=js, je
            do i=is, ie
                fx%p(i,j,k) = u%p(i,j,k)*fx%p(i,j,k)! / mesh_x%G(i,j)
            end do
        end do

        call sbp_interp_h2v%apply(fy%p,fy_work_tile,fy_tile, mesh_o%ny+1, 'y', fp, fp_tile)

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je
        do j=js, je
            do i=is, ie
                fy%p(i,j,k) = v%p(i,j,k)*fy%p(i,j,k)! / mesh_y%G(i,j)
            end do
        end do
    end do

end subroutine calc_c_sbp42_massflux_tile

end module massflux_Cgrid_mod

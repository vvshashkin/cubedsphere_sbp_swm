module halo_C_default_mod

use halo_mod,          only : halo_t, halo_vec_t
use exchange_halo_mod, only : exchange_t

implicit none

type, extends(halo_vec_t) :: halo_C_vec_default_t

    class(exchange_t), allocatable  :: exch_halo
    integer(kind=4)                 :: ts_u, te_u, ts_v, te_v
    logical, allocatable            :: is_right_edge(:), is_left_edge(:)
    logical, allocatable            :: is_top_edge(:), is_bottom_edge(:)

    contains

    procedure :: get_halo_vector => get_C_default_vector_halo

end type

contains

subroutine get_C_default_vector_halo(this,u,v,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(halo_C_vec_default_t),  intent(inout) :: this
    class(grid_field_t),          intent(inout) :: u, v
    type(domain_t),               intent(in)    :: domain
    integer(kind=4),              intent(in)    :: halo_width

    integer(kind=4) :: t,is,ie,js,je,ks,ke,i,j,nx,ny

    call this%exch_halo%do_vec(u, v, domain%parcomm)

    nx = domain%partition%nh
    ny = domain%partition%nh

    do t = this%ts_u, this%te_u
        js = domain%mesh_u%tile(t)%js
        je = domain%mesh_u%tile(t)%je
        ks = domain%mesh_u%tile(t)%ks
        ke = domain%mesh_u%tile(t)%ke

        if(this%is_left_edge(t)) then
            u%tile(t)%p(1,js:je,ks:ke) = 0.5_8*(u%tile(t)%p(1,js:je,ks:ke)+u%tile(t)%p(0,js:je,ks:ke))
            do i=1,halo_width
                u%tile(t)%p(1-i,js:je,ks:ke) = u%tile(t)%p(-i,js:je,ks:ke)
            end do
        end if
        if(this%is_right_edge(t)) then
            u%tile(t)%p(nx+1,js:je,ks:ke) = 0.5_8*(u%tile(t)%p(nx+1,js:je,ks:ke)+u%tile(t)%p(nx+2,js:je,ks:ke))
            do i=1,halo_width
                u%tile(t)%p(nx+1+i,js:je,ks:ke) = u%tile(t)%p(nx+i+2,js:je,ks:ke)
            end do
        end if
    end do

    do t = this%ts_v, this%te_v
        is = domain%mesh_v%tile(t)%is
        ie = domain%mesh_v%tile(t)%ie
        ks = domain%mesh_v%tile(t)%ks
        ke = domain%mesh_v%tile(t)%ke

        if(this%is_bottom_edge(t)) then
            v%tile(t)%p(is:ie,1,ks:ke) = 0.5_8*(v%tile(t)%p(is:ie,1,ks:ke)+v%tile(t)%p(is:ie,0,ks:ke))
            do i=1,halo_width
                v%tile(t)%p(is:ie,1-i,ks:ke) = v%tile(t)%p(is:ie,-i,ks:ke)
            end do
        end if
        if(this%is_top_edge(t)) then
            v%tile(t)%p(is:ie,ny+1,ks:ke) = 0.5_8*(v%tile(t)%p(is:ie,ny+1,ks:ke)+v%tile(t)%p(is:ie,ny+2,ks:ke))
            do i=1,halo_width
                v%tile(t)%p(is:ie,ny+1+i,ks:ke) = v%tile(t)%p(is:ie,ny+2+i,ks:ke)
            end do
        end if
    end do

end subroutine get_C_default_vector_halo

end module halo_C_default_mod

!Initialize weights to syncronize Ah-grid vector components across edges of cubed sphere
!Sync process is following:
!1) transform components from neoghbouring face
! to the coordinate system of current face
!2) make half sums with components from current face at the same point
!The precomputed weights storage structure and subroutine performing this algorithm
!is in ecs_halo_Ah_vec_sync_mod
module ecs_halo_Ah_vec_sync_factory_mod

use halo_mod,                 only : halo_vec_t
use ecs_halo_Ah_vec_sync_mod, only : ecs_halo_Ah_vec_sync_t
use domain_mod,               only : domain_t
use exchange_factory_mod,     only : create_xy_points_halo_exchange
use ecs_metric_mod,           only : ecs_b1_proto, ecs_b2_proto, ecs_a1_proto, ecs_a2_proto
use parcomm_mod,              only : parcomm_global

implicit none

contains

subroutine create_ecs_Ah_vec_sync(halo_out,domain,halo_width,components_type)

    use const_mod,                only : pi

    class(halo_vec_t), allocatable, intent(out) :: halo_out
    class(domain_t),                intent(in)  :: domain
    integer(kind=4),                intent(in)  :: halo_width
    character(len=*),               intent(in)  :: components_type

    type(ecs_halo_Ah_vec_sync_t), allocatable :: halo
    integer(kind=4), parameter    :: halo_width_edges=1
    integer(kind=4) :: t, is, ie, js, je, i, j
    real(kind=8)    :: alpha, beta, a(3), b(3)

    allocate(halo)
    halo%exch_edges =  &
              create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                           domain%topology,  halo_width_edges, 'cross')
    halo%ts = domain%mesh_xy%ts
    halo%te = domain%mesh_xy%te

    allocate(halo%tile(halo%ts:halo%te))

    halo%components_type = components_type

    select case(components_type)
    case("contravariant")

    do t=halo%ts,halo%te
        is = domain%mesh_xy%tile(t)%is
        ie = domain%mesh_xy%tile(t)%ie
        js = domain%mesh_xy%tile(t)%js
        je = domain%mesh_xy%tile(t)%je

        if(js == 1) then
            allocate(halo%tile(t)%qb(is:ie))
            do i = is,ie
                alpha = -0.25_8*pi+(i-1)*domain%mesh_xy%tile(t)%hx
                b(1:3) = ecs_b1_proto(alpha,-0.25_8*pi)
                a(1:3) = ecs_a2_proto(alpha, 0.25_8*pi)
                a(1:3) = [a(1),-a(3),a(2)]
                halo%tile(t)%qb(i) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(je == domain%mesh_xy%tile(t)%ny) then
            allocate(halo%tile(t)%qt(is:ie))
            do i = is,ie
                alpha = -0.25_8*pi+(i-1)*domain%mesh_xy%tile(t)%hx
                b(1:3) = ecs_b1_proto(alpha, 0.25_8*pi)
                a(1:3) = ecs_a2_proto(alpha,-0.25_8*pi)
                a(1:3) = [a(1),a(3),-a(2)]
                halo%tile(t)%qt(i) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(is == 1) then
            allocate(halo%tile(t)%ql(js:je))
            do j = js,je
                beta = -0.25_8*pi+(j-1)*domain%mesh_xy%tile(t)%hx
                b(1:3) = ecs_b2_proto(-0.25_8*pi, beta)
                a(1:3) = ecs_a1_proto( 0.25_8*pi, beta)
                a(1:3) = [-a(3),a(2),a(1)]
                halo%tile(t)%ql(j) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(ie == domain%mesh_xy%tile(t)%nx) then
            allocate(halo%tile(t)%qr(js:je))
            do j = js,je
                beta = -0.25_8*pi+(j-1)*domain%mesh_xy%tile(t)%hx
                b(1:3) = ecs_b2_proto( 0.25_8*pi, beta)
                a(1:3) = ecs_a1_proto(-0.25_8*pi, beta)
                a(1:3) = [a(3),a(2),-a(1)]
                halo%tile(t)%qr(j) = sum(a(1:3)*b(1:3))
            end do
        end if
    end do

    case("covariant")
    do t=domain%mesh_xy%ts, domain%mesh_xy%te
        is = domain%mesh_xy%tile(t)%is
        ie = domain%mesh_xy%tile(t)%ie
        js = domain%mesh_xy%tile(t)%js
        je = domain%mesh_xy%tile(t)%je

        if(js == 1) then
            allocate(halo%tile(t)%qb(is:ie))
            do i = is,ie
                alpha = -0.25_8*pi+(i-1)*domain%mesh_xy%tile(t)%hx
                a(1:3) = ecs_a2_proto(alpha,-0.25_8*pi)
                b(1:3) = ecs_b1_proto(alpha, 0.25_8*pi)
                b(1:3) = [b(1),-b(3),b(2)]
                halo%tile(t)%qb(i) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(je == domain%mesh_xy%tile(t)%ny) then
            allocate(halo%tile(t)%qt(is:ie))
            do i = is,ie
                alpha = -0.25_8*pi+(i-1)*domain%mesh_xy%tile(t)%hx
                a(1:3) = ecs_a2_proto(alpha, 0.25_8*pi)
                b(1:3) = ecs_b1_proto(alpha,-0.25_8*pi)
                b(1:3) = [b(1),b(3),-b(2)]
                halo%tile(t)%qt(i) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(is == 1) then
            allocate(halo%tile(t)%ql(js:je))
            do j = js,je
                beta = -0.25_8*pi+(j-1)*domain%mesh_xy%tile(t)%hx
                a(1:3) = ecs_a1_proto(-0.25_8*pi, beta)
                b(1:3) = ecs_b2_proto( 0.25_8*pi, beta)
                b(1:3) = [-b(3),b(2),b(1)]
                halo%tile(t)%ql(j) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(ie == domain%mesh_xy%tile(t)%nx) then
            allocate(halo%tile(t)%qr(js:je))
            do j = js,je
                beta = -0.25_8*pi+(j-1)*domain%mesh_xy%tile(t)%hx
                a(1:3) = ecs_a1_proto( 0.25_8*pi, beta)
                b(1:3) = ecs_b2_proto(-0.25_8*pi, beta)
                b(1:3) = [b(3),b(2),-b(1)]
                halo%tile(t)%qr(j) = sum(a(1:3)*b(1:3))
            end do
        end if
    end do

    case default
        call parcomm_global%abort("create_ecs_Ah_vec_sync: unknown components type: "//&
                                  components_type)
    end select

    call move_alloc(halo,halo_out)
end subroutine create_ecs_Ah_vec_sync

end module ecs_halo_Ah_vec_sync_factory_mod

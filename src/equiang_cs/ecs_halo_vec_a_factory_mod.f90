!Initialization routine for halo_vec object on eq.cubsph grid
!no staggering:
!1)save information about scalar halo-procedure (interpolation of
!  individual components)
!2)precompute vector-component transformation matrices from source panel
!  to target panel
module ecs_halo_vec_a_factory_mod
use ecs_halo_mod,       only : ecs_halo_t
use ecs_halo_vec_a_mod, only : ecs_halo_vec_t, ecs_tile_halo_vec_t
use parcomm_mod,        only : parcomm_global

implicit none

integer, parameter :: corner_halo_width = 5!minimum halo-width to compute 2x2 corner-halo-areas

private
public   :: create_ecs_A_vec_halo_procedure

contains

subroutine create_ecs_A_vec_halo_procedure(halo_out,domain,halo_width)
    use halo_mod,               only : halo_vec_t
    use domain_mod,             only : domain_t
    use exchange_factory_mod,   only : create_o_points_halo_exchange

    class(halo_vec_t), allocatable, intent(out) :: halo_out
    class(domain_t),                intent(in)  :: domain
    integer(kind=4),                intent(in)  :: halo_width

    !locals
    type(ecs_halo_vec_t), allocatable :: halo
    integer(kind=4)      :: ex_halo_width = 8
    integer(kind=4)      :: ts, te, is,ie, js, je, nh, t
    real(kind=8)         :: hx

    allocate(halo)
    ts = domain%partition%ts
    te = domain%partition%te
    nh = domain%partition%nh
    halo%ts = ts
    halo%te = te
    allocate(halo%tile(ts:te))

    halo%exch_halo = create_o_points_halo_exchange(domain%partition, domain%parcomm, domain%topology, &
                                                 ex_halo_width, 'full')
    do t=ts,te
        hx = domain%mesh_p%tile(t)%hx
        call domain%partition%tiles_o%tile(t)%getind(is,ie,js,je)
        call init_ecs_tile_halo_vect(halo%tile(t),domain%partition%panel_map(t), &
                                     is,ie,js,je,nh,halo_width,hx)
    end do

    call move_alloc(halo, halo_out)

end

subroutine init_ecs_tile_halo_vect(tile_halo, panel_ind,is,ie,js,je, &
                              nx,halo_width,hx)

    use ecs_halo_factory_mod, only : init_ecs_tile_halo

    integer(kind=4),   intent(in) :: panel_ind,is,ie,js,je,nx
    integer(kind=4),   intent(in) :: halo_width
    real(kind=8),      intent(in) :: hx

    type(ecs_tile_halo_vec_t), intent(inout) :: tile_halo

    !locals
    integer(kind=4) ish, ieh, jsh, jeh
    integer(kind=4) cor_hw

    allocate(tile_halo%scalar_halo)
    call init_ecs_tile_halo(tile_halo%scalar_halo,is,ie,js,je,nx,halo_width,hx)

    tile_halo%n = nx
    tile_halo%lhalo = tile_halo%scalar_halo%lhalo
    tile_halo%panel_ind = panel_ind
    tile_halo%halo_width = max(halo_width,corner_halo_width)
    tile_halo%corner_halo_width = corner_halo_width
    cor_hw = corner_halo_width

    if(tile_halo%lhalo(1) .or. tile_halo%lhalo(2)) then
        ish = max(1,is-halo_width)
        ish = max(1,minval(tile_halo%scalar_halo%indx(ish,1:halo_width)-1))
        ieh = min(nx,ie+halo_width)
        ieh = min(nx,maxval(tile_halo%scalar_halo%indx(ieh,1:halo_width))+2)
        tile_halo%ish = ish; tile_halo%ieh = ieh
    end if
    if(tile_halo%lhalo(3) .or. tile_halo%lhalo(4)) then
        jsh = max(1,js-halo_width)
        jsh = max(1,minval(tile_halo%scalar_halo%indy(jsh,1:halo_width)-1))
        jeh = min(nx,je+halo_width)
        jeh = min(nx,maxval(tile_halo%scalar_halo%indy(jeh,1:halo_width))+2)
        tile_halo%jsh = jsh; tile_halo%jeh = jeh
    end if

    if(tile_halo%lhalo(1)) then
        allocate(tile_halo%TM1(4,ish:ieh,1:max(halo_width,cor_hw)))
        call init_transform_matrix(tile_halo%TM1,ish,ieh,max(halo_width,cor_hw), hx, 1)
    end if
    if(tile_halo%lhalo(2)) then
        allocate(tile_halo%TM2(4,ish:ieh,1:max(halo_width,cor_hw)))
        call init_transform_matrix(tile_halo%TM2,ish,ieh,max(halo_width,cor_hw), hx, 2)
    end if
    if(tile_halo%lhalo(3)) then
        allocate(tile_halo%TM3(4,jsh:jeh,1:max(halo_width,cor_hw)))
        call init_transform_matrix(tile_halo%TM3,jsh,jeh,max(halo_width,cor_hw), hx, 3)
    end if
    if(tile_halo%lhalo(4)) then
        allocate(tile_halo%TM4(4,jsh:jeh,1:max(halo_width,cor_hw)))
        call init_transform_matrix(tile_halo%TM4,jsh,jeh,max(halo_width,cor_hw), hx, 4)
    end if
end subroutine init_ecs_tile_halo_vect

subroutine init_transform_matrix(TM, i1, i2, hw, hx, edge_num)
    use const_mod,        only : pi
    use ecs_metric_mod,   only : ecs_a1_proto, ecs_a2_proto, &
                                 ecs_b1_proto, ecs_b2_proto, &
                                 ecs_point_r_proto

    integer(kind=4), intent(in)  :: i1, i2, hw
    real(kind=8),    intent(in)  :: hx
    integer(kind=4), intent(in)  :: edge_num
    real(kind=8),    intent(out) :: TM(4,i1:i2,1:hw)

    !locals
    integer(kind=4) i, j
    real(kind=8) a1s(3), a2s(3), b1t(3), b2t(3), v(3), r(3)
    real(kind=8) alpha, beta, alpha_s, beta_s, alpha_t, beta_t
    integer(kind=4), parameter :: rot_matrix(3,3,4) = reshape(&
                                        [[[ 1, 0, 0],[ 0, 0,-1],[ 0, 1, 0]],&
                                         [[ 1, 0, 0],[ 0, 0, 1],[ 0,-1, 0]],&
                                         [[ 0, 0,-1],[ 0, 1, 0],[ 1, 0, 0]],&
                                         [[ 0, 0, 1],[ 0, 1, 0],[-1, 0, 0]]],&
                                         [3,3,4])
    integer(kind=4), parameter :: ab_shift(4) = [0,0,1,1]
    integer(kind=4), parameter :: b_sign(4)   = [1,-1,1,-1]

    do j=1, hw
        do i=i1,i2
            alpha = -0.25_8*pi+(i-0.5_8)*hx
            beta  =  b_sign(edge_num)*(0.25_8*pi-(j-0.5_8)*hx)
            alpha_s = (1-ab_shift(edge_num))*alpha+ab_shift(edge_num)*beta
            beta_s  = ab_shift(edge_num)*alpha+(1-ab_shift(edge_num))*beta
            !covariant basis vectors at source panel
            v   = ecs_a1_proto(alpha_s, beta_s)
            a1s(1) = rot_matrix(1,1,edge_num)*v(1)+rot_matrix(2,1,edge_num)*v(2)+rot_matrix(3,1,edge_num)*v(3)
            a1s(2) = rot_matrix(1,2,edge_num)*v(1)+rot_matrix(2,2,edge_num)*v(2)+rot_matrix(3,2,edge_num)*v(3)
            a1s(3) = rot_matrix(1,3,edge_num)*v(1)+rot_matrix(2,3,edge_num)*v(2)+rot_matrix(3,3,edge_num)*v(3)
            v   = ecs_a2_proto(alpha_s, beta_s)
            a2s(1) = rot_matrix(1,1,edge_num)*v(1)+rot_matrix(2,1,edge_num)*v(2)+rot_matrix(3,1,edge_num)*v(3)
            a2s(2) = rot_matrix(1,2,edge_num)*v(1)+rot_matrix(2,2,edge_num)*v(2)+rot_matrix(3,2,edge_num)*v(3)
            a2s(3) = rot_matrix(1,3,edge_num)*v(1)+rot_matrix(2,3,edge_num)*v(2)+rot_matrix(3,3,edge_num)*v(3)
            v   = ecs_point_r_proto(alpha_s, beta_s)
            r(1) = rot_matrix(1,1,edge_num)*v(1)+rot_matrix(2,1,edge_num)*v(2)+rot_matrix(3,1,edge_num)*v(3)
            r(2) = rot_matrix(1,2,edge_num)*v(1)+rot_matrix(2,2,edge_num)*v(2)+rot_matrix(3,2,edge_num)*v(3)
            r(3) = rot_matrix(1,3,edge_num)*v(1)+rot_matrix(2,3,edge_num)*v(2)+rot_matrix(3,3,edge_num)*v(3)

            alpha_t = atan(r(1)/r(3))
            beta_t  = atan(r(2)/r(3))
            b1t = ecs_b1_proto(alpha_t, beta_t)
            b2t = ecs_b2_proto(alpha_t, beta_t)
            !transform matrix <- dot-products of target_contra_vec*source_cov_vec
            !TM = |1  2|
            !     |3  4|
            TM(1,i,j) = sum(a1s(1:3)*b1t(1:3))
            TM(2,i,j) = sum(a2s(1:3)*b1t(1:3))
            TM(3,i,j) = sum(a1s(1:3)*b2t(1:3))
            TM(4,i,j) = sum(a2s(1:3)*b2t(1:3))
        end do
    end do

end subroutine init_transform_matrix

end module ecs_halo_vec_a_factory_mod

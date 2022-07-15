!Object to interpolate vector contravariant component values to virtual points beyond face edge
!C-staggered grid.

module ecs_halo_vec_c_mod

use halo_mod,          only : halo_vec_t
use exchange_halo_mod, only : exchange_t

implicit none

type, extends(halo_vec_t) :: ecs_halo_vec_c_t

    integer(kind=4)                         :: ts, te
    integer(kind=4)                         :: exchange_width
    class(exchange_t),          allocatable :: exch_halo
    type(ecs_tile_halo_cvec_t), allocatable :: tile(:)

    contains

    procedure :: get_halo_vector => get_ecs_c_vector_halo

end type

type, extends(halo_vec_t) :: ecs_halo_vec_c_cov_t

    integer(kind=4)                         :: ts, te
    integer(kind=4)                         :: exchange_width
    class(exchange_t),          allocatable :: exch_halo
    type(ecs_tile_halo_cvec_cov_t), allocatable :: tile(:)

    contains

    procedure :: get_halo_vector => get_ecs_c_vector_cov_halo

end type

type ecs_tile_halo_cvec_t
    logical is_left_edge, is_right_edge, is_bottom_edge, is_top_edge
    !interpolation weights and stencils for the vector component normal to edge
    real(kind=8),    allocatable :: wn_left(:,:,:), wn_right(:,:,:), wn_bottom(:,:,:), wn_top(:,:,:)
    integer(kind=4), allocatable :: in_left(:,:), in_right(:,:), in_bottom(:,:), in_top(:,:)
    !interpolation weights and stencils for the vector component tangential to edge
    real(kind=8),    allocatable :: wt_left(:,:,:), wt_right(:,:,:), wt_bottom(:,:,:), wt_top(:,:,:)
    integer(kind=4), allocatable :: it_left(:,:), it_right(:,:), it_bottom(:,:), it_top(:,:)
    !transformation from source panel contravariant to target panel contravariant
    real(kind=8),    allocatable :: qtt_left(:,:), qtt_right(:,:), qtt_bottom(:,:), qtt_top(:,:)
    real(kind=8),    allocatable :: qtn_left(:,:), qtn_right(:,:), qtn_bottom(:,:), qtn_top(:,:)
    real(kind=8),    allocatable :: wtn_left(:,:,:), wtn_right(:,:,:), wtn_bottom(:,:,:), wtn_top(:,:,:)
    integer(kind=4), allocatable :: itn_left(:,:), itn_right(:,:), itn_bottom(:,:), itn_top(:,:)
    contains

    procedure :: interpv => interp_ecs_tile_halo_cvec
end type ecs_tile_halo_cvec_t

type ecs_tile_halo_cvec_cov_t
    logical is_left_edge, is_right_edge, is_bottom_edge, is_top_edge
    !interpolation weights and stencils for the vector component normal to edge
    real(kind=8),    allocatable :: wn_left(:,:,:), wn_right(:,:,:), wn_bottom(:,:,:), wn_top(:,:,:)
    integer(kind=4), allocatable :: in_left(:,:), in_right(:,:), in_bottom(:,:), in_top(:,:)
    !interpolation weights and stencils for the vector component tangential to edge
    real(kind=8),    allocatable :: wt_left(:,:,:), wt_right(:,:,:), wt_bottom(:,:,:), wt_top(:,:,:)
    integer(kind=4), allocatable :: it_left(:,:), it_right(:,:), it_bottom(:,:), it_top(:,:)
    real(kind=8),    allocatable :: ct_left(:,:), ct_right(:,:), ct_bottom(:,:), ct_top(:,:)
    !transformation from source panel contravariant to target panel contravariant
    real(kind=8),    allocatable :: qnt_left(:,:), qnt_right(:,:), qnt_bottom(:,:), qnt_top(:,:)
    real(kind=8),    allocatable :: qnn_left(:,:), qnn_right(:,:), qnn_bottom(:,:), qnn_top(:,:)
    real(kind=8),    allocatable :: wnt_left(:,:,:), wnt_right(:,:,:), wnt_bottom(:,:,:), wnt_top(:,:,:)
    integer(kind=4), allocatable :: int_left(:,:), int_right(:,:), int_bottom(:,:), int_top(:,:)
    contains

    procedure :: interpv => interp_ecs_tile_halo_cvec_cov
end type ecs_tile_halo_cvec_cov_t

contains

subroutine get_ecs_c_vector_halo(this,u,v,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(ecs_halo_vec_c_t),  intent(inout) :: this
    class(grid_field_t),      intent(inout) :: u,v
    type(domain_t),           intent(in)    :: domain
    integer(kind=4),          intent(in)    :: halo_width

    integer(kind=4) t

    call this%exch_halo%do_vec(u,v, domain%parcomm)

    do t=this%ts,this%te
        call this%tile(t)%interpv(u%tile(t),v%tile(t),domain%partition%tiles_x%tile(t), &
                                  domain%partition%tiles_y%tile(t),halo_width,this%exchange_width,t)
    end do
end subroutine get_ecs_c_vector_halo

subroutine get_ecs_c_vector_cov_halo(this,u,v,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(ecs_halo_vec_c_cov_t),  intent(inout) :: this
    class(grid_field_t),          intent(inout) :: u,v
    type(domain_t),               intent(in)    :: domain
    integer(kind=4),              intent(in)    :: halo_width

    integer(kind=4) t

    call this%exch_halo%do_vec(u,v, domain%parcomm)

    do t=this%ts,this%te
        call this%tile(t)%interpv(u%tile(t),v%tile(t),domain%partition%tiles_x%tile(t), &
                                  domain%partition%tiles_y%tile(t),halo_width,this%exchange_width,t)
    end do
end subroutine get_ecs_c_vector_cov_halo

subroutine interp_ecs_tile_halo_cvec(this,u,v,tile_u,tile_v,halo_width,exchange_width,t)
    use grid_field_mod, only: tile_field_t
    use tile_mod,       only: tile_t

    class(ecs_tile_halo_cvec_t), intent(in)    :: this
    type(tile_field_t),          intent(inout) :: u, v
    type(tile_t),                intent(in)    :: tile_u, tile_v
    integer(kind=4),             intent(in)    :: halo_width, exchange_width,t

    integer(kind=4) :: i,is,ie,js,je,ks,ke,wst,wend
    integer(kind=4) :: is1,ie1,js1,je1
    integer(kind=4) :: jsv,jev,isv,iev

    is = tile_u%is
    ie = tile_u%ie
    ie1 = tile_v%ie
    js = tile_u%js
    je = tile_u%je
    js1 = tile_v%js
    je1 = tile_v%je
    ks = tile_u%ks
    ke = tile_u%ke
    jsv = js-exchange_width
    jev = je+exchange_width

    if(this%is_left_edge) then
        u%p(1,js:je,ks:ke) = 0.5_8*(u%p(1,js:je,ks:ke)+u%p(0,js:je,ks:ke))
        do i=1,halo_width+1
            u%p(1-i,jsv:jev,ks:ke) = u%p(-i,jsv:jev,ks:ke)
        end do
        wst  = lbound(this%wt_left,1)
        wend = ubound(this%wt_left,1)
        call interp_edge(v,this%wt_left,this%it_left,wst,wend,js1,je1,ks,ke,halo_width,ie1,'left')
        wst  = lbound(this%wtn_left,1)
        wend = ubound(this%wtn_left,1)
        call tangent_comp_transform(v,u,this%wtn_left,this%itn_left,this%qtt_left, this%qtn_left, &
                                    wst,wend,js1,je1,ks,ke,halo_width,ie1,'left')
        wst  = lbound(this%wn_left,1)
        wend = ubound(this%wn_left,1)
        call interp_edge(u,this%wn_left,this%in_left,wst,wend,js,je,ks,ke,halo_width,ie,'left')
    end if
    if(this%is_right_edge) then
        u%p(ie,js:je,ks:ke) = 0.5_8*(u%p(ie,js:je,ks:ke)+u%p(ie+1,js:je,ks:ke))
        do i=1,halo_width+1
            u%p(ie+i,jsv:jev,ks:ke) = u%p(ie+i+1,jsv:jev,ks:ke)
        end do
        wst  = lbound(this%wt_right,1)
        wend = ubound(this%wt_right,1)
        call interp_edge(v,this%wt_right,this%it_right,wst,wend,js1,je1,ks,ke,halo_width,ie1,'right')
        wst  = lbound(this%wtn_right,1)
        wend = ubound(this%wtn_right,1)
        call tangent_comp_transform(v,u,this%wtn_right,this%itn_right,this%qtt_right, this%qtn_right, &
                                    wst,wend,js1,je1,ks,ke,halo_width,ie1,'right')
        wst  = lbound(this%wn_right,1)
        wend = ubound(this%wn_right,1)
        call interp_edge(u,this%wn_right,this%in_right,wst,wend,js,je,ks,ke,halo_width,ie,'right')
    end if

    is = tile_v%is
    ie = tile_v%ie
    is1 = tile_u%is
    ie1 = tile_u%ie
    js = tile_v%js
    je = tile_v%je
    je1 = tile_u%je
    ks = tile_v%ks
    ke = tile_v%ke
    isv = is-exchange_width
    iev = ie+exchange_width

    if(this%is_bottom_edge) then
        v%p(is:ie,1,ks:ke) = 0.5_8*(v%p(is:ie,1,ks:ke)+v%p(is:ie,0,ks:ke))
        do i=1,halo_width+1
            v%p(isv:iev,1-i,ks:ke) = v%p(isv:iev,-i,ks:ke)
        end do
        wst  = lbound(this%wt_bottom,1)
        wend = ubound(this%wt_bottom,1)
        call interp_edge(u,this%wt_bottom,this%it_bottom,wst,wend,is1,ie1,ks,ke,halo_width,je1,'bottom')
        wst  = lbound(this%wtn_bottom,1)
        wend = ubound(this%wtn_bottom,1)
        call tangent_comp_transform(u,v,this%wtn_bottom,this%itn_bottom,this%qtt_bottom, this%qtn_bottom, &
                                    wst,wend,is1,ie1,ks,ke,halo_width,je1,'bottom')
        wst  = lbound(this%wn_bottom,1)
        wend = ubound(this%wn_bottom,1)
        call interp_edge(v,this%wn_bottom,this%in_bottom,wst,wend,is,ie,ks,ke,halo_width,je,'bottom')
    end if
    if(this%is_top_edge) then
        v%p(is:ie,je,ks:ke) = 0.5_8*(v%p(is:ie,je,ks:ke)+v%p(is:ie,je+1,ks:ke))
        do i=1,halo_width+1
            v%p(isv:iev,je+i,ks:ke) = v%p(isv:iev,je+i+1,ks:ke)
        end do
        wst  = lbound(this%wt_top,1)
        wend = ubound(this%wt_top,1)
        call interp_edge(u,this%wt_top,this%it_top,wst,wend,is1,ie1,ks,ke,halo_width,je1,'top')
        wst  = lbound(this%wtn_top,1)
        wend = ubound(this%wtn_top,1)
        call tangent_comp_transform(u,v,this%wtn_top,this%itn_top,this%qtt_top,this%qtn_top, &
                                    wst,wend,is1,ie1,ks,ke,halo_width,je1,'top')
        wst  = lbound(this%wn_top,1)
        wend = ubound(this%wn_top,1)
        call interp_edge(v,this%wn_top,this%in_top,wst,wend,is,ie,ks,ke,halo_width,je,'top')
    end if

end subroutine interp_ecs_tile_halo_cvec

subroutine interp_ecs_tile_halo_cvec_cov(this,u,v,tile_u,tile_v,halo_width,exchange_width,t)
    use grid_field_mod, only: tile_field_t
    use tile_mod,       only: tile_t

    class(ecs_tile_halo_cvec_cov_t), intent(in)    :: this
    type(tile_field_t),          intent(inout) :: u, v
    type(tile_t),                intent(in)    :: tile_u, tile_v
    integer(kind=4),             intent(in)    :: halo_width, exchange_width,t

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke,wst,wend
    integer(kind=4) :: is1,ie1,js1,je1
    integer(kind=4) :: jsv,jev,isv,iev

    is = tile_u%is
    ie = tile_u%ie
    ie1 = tile_v%ie
    js = tile_u%js
    je = tile_u%je
    js1 = tile_v%js
    je1 = tile_v%je
    ks = tile_u%ks
    ke = tile_u%ke
    jsv = js-exchange_width
    jev = je+exchange_width

    if(this%is_left_edge) then
        do i=1,halo_width+1
            u%p(1-i,jsv:jev,ks:ke) = u%p(-i,jsv:jev,ks:ke)
        end do
        wst  = lbound(this%wn_left,1)
        wend = ubound(this%wn_left,1)
        call interp_edge(u,this%wn_left,this%in_left,wst,wend,js,je,ks,ke,halo_width,ie,'left')
        wst  = lbound(this%wnt_left,1)
        wend = ubound(this%wnt_left,1)
        call normal_comp_transform(u,v,this%wnt_left,this%int_left,this%qnt_left, this%qnn_left, &
                                    wst,wend,js1,je1,ks,ke,halo_width,ie1,'left')

        wst  = lbound(this%wt_left,1)
        wend = ubound(this%wt_left,1)
        call interp_edge(v,this%wt_left,this%it_left,wst,wend,js1,je1,ks,ke,halo_width,ie1,'left')
        do k=ks,ke
            do i=1, halo_width
                do j=js1,je1
                    v%p(1-i,j,k) = v%p(1-i,j,k)*this%ct_left(j,i)
                end do
            end do
        end do
    end if
    if(this%is_right_edge) then
        do i=1,halo_width+1
            u%p(ie+i,jsv:jev,ks:ke) = u%p(ie+i+1,jsv:jev,ks:ke)
        end do
        wst  = lbound(this%wn_right,1)
        wend = ubound(this%wn_right,1)
        call interp_edge(u,this%wn_right,this%in_right,wst,wend,js,je,ks,ke,halo_width,ie,'right')
        wst  = lbound(this%wnt_right,1)
        wend = ubound(this%wnt_right,1)
        call normal_comp_transform(u,v,this%wnt_right,this%int_right,this%qnt_right, this%qnn_right, &
                                   wst,wend,js1,je1,ks,ke,halo_width,ie1,'right')

        wst  = lbound(this%wt_right,1)
        wend = ubound(this%wt_right,1)
        call interp_edge(v,this%wt_right,this%it_right,wst,wend,js1,je1,ks,ke,halo_width,ie1,'right')
        do k=ks,ke
            do i=1, halo_width
                do j=js1,je1
                    v%p(ie1+i,j,k) = v%p(ie1+i,j,k)*this%ct_right(j,i)
                end do
            end do
        end do
    end if

    is = tile_v%is
    ie = tile_v%ie
    is1 = tile_u%is
    ie1 = tile_u%ie
    js = tile_v%js
    je = tile_v%je
    je1 = tile_u%je
    ks = tile_v%ks
    ke = tile_v%ke
    isv = is-exchange_width
    iev = ie+exchange_width

    if(this%is_bottom_edge) then
        do i=1,halo_width+1
            v%p(isv:iev,1-i,ks:ke) = v%p(isv:iev,-i,ks:ke)
        end do
        wst  = lbound(this%wn_bottom,1)
        wend = ubound(this%wn_bottom,1)
        call interp_edge(v,this%wn_bottom,this%in_bottom,wst,wend,is,ie,ks,ke,halo_width,je,'bottom')
        wst  = lbound(this%wnt_bottom,1)
        wend = ubound(this%wnt_bottom,1)
        call normal_comp_transform(v,u,this%wnt_bottom,this%int_bottom,this%qnt_bottom, this%qnn_bottom, &
                                   wst,wend,is1,ie1,ks,ke,halo_width,je1,'bottom')

        wst  = lbound(this%wt_bottom,1)
        wend = ubound(this%wt_bottom,1)
        call interp_edge(u,this%wt_bottom,this%it_bottom,wst,wend,is1,ie1,ks,ke,halo_width,je1,'bottom')
        do k=ks,ke
            do j=1, halo_width
                do i=is1,ie1
                    u%p(i,1-j,k) = u%p(i,1-j,k)*this%ct_bottom(i,j)
                end do
            end do
        end do
    end if
    if(this%is_top_edge) then
        do i=1,halo_width+1
            v%p(isv:iev,je+i,ks:ke) = v%p(isv:iev,je+i+1,ks:ke)
        end do
        wst  = lbound(this%wn_top,1)
        wend = ubound(this%wn_top,1)
        call interp_edge(v,this%wn_top,this%in_top,wst,wend,is,ie,ks,ke,halo_width,je,'top')
        wst  = lbound(this%wnt_top,1)
        wend = ubound(this%wnt_top,1)
        call normal_comp_transform(v,u,this%wnt_top,this%int_top,this%qnt_top,this%qnn_top, &
                                   wst,wend,is1,ie1,ks,ke,halo_width,je1,'top')
        wst  = lbound(this%wt_top,1)
        wend = ubound(this%wt_top,1)
        call interp_edge(u,this%wt_top,this%it_top,wst,wend,is1,ie1,ks,ke,halo_width,je1,'top')
        do k=ks,ke
            do j=1, halo_width
                do i=is1,ie1
                    u%p(i,je1+j,k) = u%p(i,je1+j,k)*this%ct_top(i,j)
                end do
            end do
        end do
    end if

end subroutine interp_ecs_tile_halo_cvec_cov

subroutine interp_edge(f, w, ind, wst, wend, si, ei, ks, ke, halo_width, nh, edge)
    use grid_field_mod, only : tile_field_t
    use parcomm_mod,    only : parcomm_global

    type(tile_field_t), intent(inout) :: f
    integer(kind=4),    intent(in)    :: wst, wend
    integer(kind=4),    intent(in)    :: si, ei, ks, ke, halo_width, nh
    real(kind=8),       intent(in)    :: w(wst:wend,si:ei,1:halo_width)
    integer(kind=4),    intent(in)    :: ind(si:ei,1:halo_width)
    character(len=*),   intent(in)    :: edge

    integer(kind=4) :: i, j, k
    real(kind=8)    :: f1(si:ei)

    do k=ks,ke
        do j=1,halo_width
            if(edge=="left") then
                do i=si,ei
                    f1(i) = sum(w(wst:wend,i,j)*f%p(1-j,ind(i,j)+wst:ind(i,j)+wend,k))
                end do
                f%p(1-j,si:ei,k) = f1(si:ei)
            elseif(edge=="right") then
                do i=si,ei
                    f1(i) = sum(w(wst:wend,i,j)*f%p(nh+j,ind(i,j)+wst:ind(i,j)+wend,k))
                end do
                f%p(nh+j,si:ei,k) = f1(si:ei)
            elseif(edge=="bottom") then
                do i=si,ei
                    f1(i) = sum(w(wst:wend,i,j)*f%p(ind(i,j)+wst:ind(i,j)+wend,1-j,k))
                end do
                f%p(si:ei,1-j,k) = f1(si:ei)
            elseif(edge=="top") then
                do i=si,ei
                    f1(i) = sum(w(wst:wend,i,j)*f%p(ind(i,j)+wst:ind(i,j)+wend,nh+j,k))
                end do
                f%p(si:ei,nh+j,k) = f1(si:ei)
            else
                call parcomm_global%abort("ecs_halo_vec_c_mod: interp_edge, unknown edge: "//edge)
            end if
        end do
    end do

end subroutine interp_edge

subroutine tangent_comp_transform(ut, un, w, ind, qtt, qtn, wst, wend, si, ei, ks, ke, halo_width, nh, edge)
    use grid_field_mod, only : tile_field_t
    use parcomm_mod,    only : parcomm_global

    type(tile_field_t), intent(inout) :: ut !contravariant velocity components: tangent& normal to edge
    type(tile_field_t), intent(in)    :: un !contravariant velocity components: tangent& normal to edge
    integer(kind=4),    intent(in)    :: wst, wend
    integer(kind=4),    intent(in)    :: si, ei, ks, ke, halo_width, nh
    real(kind=8),       intent(in)    :: w(wst:wend,si:ei,1:halo_width)
    integer(kind=4),    intent(in)    :: ind(si:ei,1:halo_width)
    real(kind=8),       intent(in)    :: qtt(si:ei,1:halo_width), qtn(si:ei,1:halo_width)
    character(len=*),   intent(in)    :: edge

    integer(kind=4) :: i, j, k
    real(kind=8)    :: un1, un2, un3, un4, un0

    do k=ks,ke
        select case(edge)
        case("left")
            do i=si,ei
                un1 = sum(w(wst:wend,i,1)*un%p(1,ind(i,1)+wst:ind(i,1)+wend,k))
                un2 = sum(w(wst:wend,i,1)*un%p(0,ind(i,1)+wst:ind(i,1)+wend,k))
                un3 = sum(w(wst:wend,i,1)*un%p(-1,ind(i,1)+wst:ind(i,1)+wend,k))
                un4 = sum(w(wst:wend,i,1)*un%p(-2,ind(i,1)+wst:ind(i,1)+wend,k))
                un0 = (un4-5.0_8*un3+15.0_8*un2+5.0_8*un1)/16.0_8
                ut%p(0,i,k) = qtt(i,1)*ut%p(0,i,k)-qtn(i,1)*un0
            end do
            do j=2,halo_width
                do i=si,ei
                    un1 = sum(w(wst:wend,i,j)*un%p(3-j,ind(i,j)+wst:ind(i,j)+wend,k))
                    un2 = sum(w(wst:wend,i,j)*un%p(2-j,ind(i,j)+wst:ind(i,j)+wend,k))
                    un3 = sum(w(wst:wend,i,j)*un%p(1-j,ind(i,j)+wst:ind(i,j)+wend,k))
                    un4 = sum(w(wst:wend,i,j)*un%p(-j,ind(i,j)+wst:ind(i,j)+wend,k))
                    un0 = (-un1+9.0_8*un2+9.0_8*un3-un4)/16.0_8
                    ut%p(1-j,i,k) = qtt(i,j)*ut%p(1-j,i,k)-qtn(i,j)*un0
                end do
            end do
        case("right")
            do i=si,ei
                un1 = sum(w(wst:wend,i,1)*un%p(nh+1,ind(i,1)+wst:ind(i,1)+wend,k))
                un2 = sum(w(wst:wend,i,1)*un%p(nh+2,ind(i,1)+wst:ind(i,1)+wend,k))
                un3 = sum(w(wst:wend,i,1)*un%p(nh+3,ind(i,1)+wst:ind(i,1)+wend,k))
                un4 = sum(w(wst:wend,i,1)*un%p(nh+4,ind(i,1)+wst:ind(i,1)+wend,k))
                un0 = (un4-5.0_8*un3+15.0_8*un2+5.0_8*un1)/16.0_8
                ut%p(nh+1,i,k) = qtt(i,1)*ut%p(nh+1,i,k)+qtn(i,1)*un0
            end do
            do j=2,halo_width
                do i=si,ei
                    un1 = sum(w(wst:wend,i,j)*un%p(nh+j-1,ind(i,j)+wst:ind(i,j)+wend,k))
                    un2 = sum(w(wst:wend,i,j)*un%p(nh+j  ,ind(i,j)+wst:ind(i,j)+wend,k))
                    un3 = sum(w(wst:wend,i,j)*un%p(nh+j+1,ind(i,j)+wst:ind(i,j)+wend,k))
                    un4 = sum(w(wst:wend,i,j)*un%p(nh+j+2,ind(i,j)+wst:ind(i,j)+wend,k))
                    un0 = (-un1+9.0_8*un2+9.0_8*un3-un4)/16.0_8
                    ut%p(nh+j,i,k) = qtt(i,j)*ut%p(nh+j,i,k)+qtn(i,j)*un0
                end do
            end do
        case("bottom")
            do i=si,ei
                un1 = sum(w(wst:wend,i,1)*un%p(ind(i,1)+wst:ind(i,1)+wend,1,k))
                un2 = sum(w(wst:wend,i,1)*un%p(ind(i,1)+wst:ind(i,1)+wend,0,k))
                un3 = sum(w(wst:wend,i,1)*un%p(ind(i,1)+wst:ind(i,1)+wend,-1,k))
                un4 = sum(w(wst:wend,i,1)*un%p(ind(i,1)+wst:ind(i,1)+wend,-2,k))
                un0 = (un4-5.0_8*un3+15.0_8*un2+5.0_8*un1)/16.0_8
                ut%p(i,0,k) = qtt(i,1)*ut%p(i,0,k)-qtn(i,1)*un0
            end do
            do j=2,halo_width
                do i=si,ei
                    un1 = sum(w(wst:wend,i,j)*un%p(ind(i,j)+wst:ind(i,j)+wend,3-j,k))
                    un2 = sum(w(wst:wend,i,j)*un%p(ind(i,j)+wst:ind(i,j)+wend,2-j,k))
                    un3 = sum(w(wst:wend,i,j)*un%p(ind(i,j)+wst:ind(i,j)+wend,1-j,k))
                    un4 = sum(w(wst:wend,i,j)*un%p(ind(i,j)+wst:ind(i,j)+wend,-j,k))
                    un0 = (-un1+9.0_8*un2+9.0_8*un3-un4)/16.0_8
                    ut%p(i,1-j,k) = qtt(i,j)*ut%p(i,1-j,k)-qtn(i,j)*un0
                end do
            end do
        case("top")
            do i=si,ei
                un1 = sum(w(wst:wend,i,1)*un%p(ind(i,1)+wst:ind(i,1)+wend,nh+1,k))
                un2 = sum(w(wst:wend,i,1)*un%p(ind(i,1)+wst:ind(i,1)+wend,nh+2,k))
                un3 = sum(w(wst:wend,i,1)*un%p(ind(i,1)+wst:ind(i,1)+wend,nh+3,k))
                un4 = sum(w(wst:wend,i,1)*un%p(ind(i,1)+wst:ind(i,1)+wend,nh+4,k))
                un0 = (un4-5.0_8*un3+15.0_8*un2+5.0_8*un1)/16.0_8
                ut%p(i,nh+1,k) = qtt(i,1)*ut%p(i,nh+1,k)+qtn(i,1)*un0
            end do
            do j=2,halo_width
                do i=si,ei
                    un1 = sum(w(wst:wend,i,j)*un%p(ind(i,j)+wst:ind(i,j)+wend,nh+j-1,k))
                    un2 = sum(w(wst:wend,i,j)*un%p(ind(i,j)+wst:ind(i,j)+wend,nh+j  ,k))
                    un3 = sum(w(wst:wend,i,j)*un%p(ind(i,j)+wst:ind(i,j)+wend,nh+j+1,k))
                    un4 = sum(w(wst:wend,i,j)*un%p(ind(i,j)+wst:ind(i,j)+wend,nh+j+2,k))
                    un0 = (-un1+9.0_8*un2+9.0_8*un3-un4)/16.0_8
                    ut%p(i,nh+j,k) = qtt(i,j)*ut%p(i,nh+j,k)+qtn(i,j)*un0
                end do
            end do
        case default
        call parcomm_global%abort("ecs_halo_vec_c_mod, tangent_comp_transform - "//&
                                  "unkniown edge type: "//edge)
        end select
    end do
end subroutine tangent_comp_transform

subroutine normal_comp_transform(un, ut, w, ind, qnt, qnn, wst, wend, si, ei, ks, ke, halo_width, nh, edge)
    use grid_field_mod, only : tile_field_t
    use parcomm_mod,    only : parcomm_global

    type(tile_field_t), intent(inout) :: un !contravariant velocity components: tangent& normal to edge
    type(tile_field_t), intent(in)    :: ut !contravariant velocity components: tangent& normal to edge
    integer(kind=4),    intent(in)    :: wst, wend
    integer(kind=4),    intent(in)    :: si, ei, ks, ke, halo_width, nh
    real(kind=8),       intent(in)    :: w(wst:wend,si:ei,1:halo_width)
    integer(kind=4),    intent(in)    :: ind(si:ei,1:halo_width)
    real(kind=8),       intent(in)    :: qnt(si:ei,1:halo_width), qnn(si:ei,1:halo_width)
    character(len=*),   intent(in)    :: edge

    integer(kind=4) :: i, j, k
    real(kind=8)    :: ut1, ut2, ut3, ut4, ut0

    do k=ks,ke
        select case(edge)
        case("left")
            do i=si,ei
                ut1 = sum(w(wst:wend,i,1)*ut%p(0,ind(i,1)+wst:ind(i,1)+wend,k))
                ut2 = sum(w(wst:wend,i,1)*ut%p(-1,ind(i,1)+wst:ind(i,1)+wend,k))
                ut3 = sum(w(wst:wend,i,1)*ut%p(-2,ind(i,1)+wst:ind(i,1)+wend,k))
                ut4 = sum(w(wst:wend,i,1)*ut%p(-3,ind(i,1)+wst:ind(i,1)+wend,k))
                ut0 = (ut4-5.0_8*ut3+15.0_8*ut2+5.0_8*ut1)/16.0_8
                un%p(0,i,k) = qnn(i,1)*un%p(0,i,k)-qnt(i,1)*ut0
            end do
            do j=2,halo_width
                do i=si,ei
                    ut1 = sum(w(wst:wend,i,j)*ut%p(2-j,ind(i,j)+wst:ind(i,j)+wend,k))
                    ut2 = sum(w(wst:wend,i,j)*ut%p(1-j,ind(i,j)+wst:ind(i,j)+wend,k))
                    ut3 = sum(w(wst:wend,i,j)*ut%p(-j,ind(i,j)+wst:ind(i,j)+wend,k))
                    ut4 = sum(w(wst:wend,i,j)*ut%p(-j-1,ind(i,j)+wst:ind(i,j)+wend,k))
                    ut0 = (-ut1+9.0_8*ut2+9.0_8*ut3-ut4)/16.0_8
                    un%p(1-j,i,k) = qnn(i,j)*un%p(1-j,i,k)-qnt(i,j)*ut0
                end do
            end do
        case("right")
            do i=si,ei
                ut1 = sum(w(wst:wend,i,1)*ut%p(nh+1,ind(i,1)+wst:ind(i,1)+wend,k))
                ut2 = sum(w(wst:wend,i,1)*ut%p(nh+2,ind(i,1)+wst:ind(i,1)+wend,k))
                ut3 = sum(w(wst:wend,i,1)*ut%p(nh+3,ind(i,1)+wst:ind(i,1)+wend,k))
                ut4 = sum(w(wst:wend,i,1)*ut%p(nh+4,ind(i,1)+wst:ind(i,1)+wend,k))
                ut0 = (ut4-5.0_8*ut3+15.0_8*ut2+5.0_8*ut1)/16.0_8
                un%p(nh+2,i,k) = qnn(i,1)*un%p(nh+2,i,k)+qnt(i,1)*ut0
            end do
            do j=2,halo_width
                do i=si,ei
                    ut1 = sum(w(wst:wend,i,j)*ut%p(nh+j-1  ,ind(i,j)+wst:ind(i,j)+wend,k))
                    ut2 = sum(w(wst:wend,i,j)*ut%p(nh+j  ,ind(i,j)+wst:ind(i,j)+wend,k))
                    ut3 = sum(w(wst:wend,i,j)*ut%p(nh+j+1,ind(i,j)+wst:ind(i,j)+wend,k))
                    ut4 = sum(w(wst:wend,i,j)*ut%p(nh+j+2,ind(i,j)+wst:ind(i,j)+wend,k))
                    ut0 = (-ut1+9.0_8*ut2+9.0_8*ut3-ut4)/16.0_8
                    un%p(nh+j+1,i,k) = qnn(i,j)*un%p(nh+j+1,i,k)+qnt(i,j)*ut0
                end do
            end do
        case("bottom")
            do i=si,ei
                ut1 = sum(w(wst:wend,i,1)*ut%p(ind(i,1)+wst:ind(i,1)+wend,0,k))
                ut2 = sum(w(wst:wend,i,1)*ut%p(ind(i,1)+wst:ind(i,1)+wend,-1,k))
                ut3 = sum(w(wst:wend,i,1)*ut%p(ind(i,1)+wst:ind(i,1)+wend,-2,k))
                ut4 = sum(w(wst:wend,i,1)*ut%p(ind(i,1)+wst:ind(i,1)+wend,-3,k))
                ut0 = (ut4-5.0_8*ut3+15.0_8*ut2+5.0_8*ut1)/16.0_8
                un%p(i,0,k) = qnn(i,1)*un%p(i,0,k)-qnt(i,1)*ut0
            end do
            do j=2,halo_width
                do i=si,ei
                    ut1 = sum(w(wst:wend,i,j)*ut%p(ind(i,j)+wst:ind(i,j)+wend,2-j,k))
                    ut2 = sum(w(wst:wend,i,j)*ut%p(ind(i,j)+wst:ind(i,j)+wend,1-j,k))
                    ut3 = sum(w(wst:wend,i,j)*ut%p(ind(i,j)+wst:ind(i,j)+wend,-j,k))
                    ut4 = sum(w(wst:wend,i,j)*ut%p(ind(i,j)+wst:ind(i,j)+wend,-j-1,k))
                    ut0 = (-ut1+9.0_8*ut2+9.0_8*ut3-ut4)/16.0_8
                    un%p(i,-j+1,k) = qnn(i,j)*un%p(i,-j+1,k)-qnt(i,j)*ut0
                end do
            end do
        case("top")
            do i=si,ei
                ut1 = sum(w(wst:wend,i,1)*ut%p(ind(i,1)+wst:ind(i,1)+wend,nh+1,k))
                ut2 = sum(w(wst:wend,i,1)*ut%p(ind(i,1)+wst:ind(i,1)+wend,nh+2,k))
                ut3 = sum(w(wst:wend,i,1)*ut%p(ind(i,1)+wst:ind(i,1)+wend,nh+3,k))
                ut4 = sum(w(wst:wend,i,1)*ut%p(ind(i,1)+wst:ind(i,1)+wend,nh+4,k))
                ut0 =(ut4-5.0_8*ut3+15.0_8*ut2+5.0_8*ut1)/16.0_8
                un%p(i,nh+2,k) = qnn(i,1)*un%p(i,nh+2,k)+qnt(i,1)*ut0
            end do
            do j=2,halo_width
                do i=si,ei
                    ut1 = sum(w(wst:wend,i,j)*ut%p(ind(i,j)+wst:ind(i,j)+wend,nh+j-1,k))
                    ut2 = sum(w(wst:wend,i,j)*ut%p(ind(i,j)+wst:ind(i,j)+wend,nh+j  ,k))
                    ut3 = sum(w(wst:wend,i,j)*ut%p(ind(i,j)+wst:ind(i,j)+wend,nh+j+1,k))
                    ut4 = sum(w(wst:wend,i,j)*ut%p(ind(i,j)+wst:ind(i,j)+wend,nh+j+2,k))
                    ut0 = (-ut1+9.0_8*ut2+9.0_8*ut3-ut4)/16.0_8
                    un%p(i,nh+j+1,k) = qnn(i,j)*un%p(i,nh+j+1,k)+qnt(i,j)*ut0
                end do
            end do
        case default
        call parcomm_global%abort("ecs_halo_vec_c_mod, tangent_comp_transform - "//&
                                  "unkniown edge type: "//edge)
        end select
    end do
end subroutine normal_comp_transform

end module ecs_halo_vec_c_mod

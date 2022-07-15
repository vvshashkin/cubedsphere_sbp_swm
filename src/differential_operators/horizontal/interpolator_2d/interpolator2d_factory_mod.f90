module interpolator2d_factory_mod

use abstract_interpolators2d_mod,    only : interpolator2d_scalar2vec_t, &
                                            interpolator2d_vec2vec_t
use interpolator_p2uv_sbp_C_mod,     only : interpolator_p2uv_sbp_C_t,   &
                                            interpolator_pvec2uv_sbp_C_t,&
                                            interpolator_pvec2uv_sbp_Ch_t
use interpolator_uv2p_mod,           only : interpolator2d_uv2p_sbp_C_t
use interpolator_q2uv_sbp_mod,       only : interpolator_q2uv_sbp_Ch_t
use interpolator_uv2q_sbp_mod,       only : interpolator_uv2q_sbp_Ch_t
use domain_mod,                      only : domain_t
use sbp_factory_mod,                 only : create_sbp_operator
use exchange_factory_mod,            only : create_o_points_halo_exchange,  &
                                            create_xy_points_halo_exchange, &
                                            create_symmetric_halo_vec_exchange_C
use parcomm_mod,                     only : parcomm_global

implicit none

contains

subroutine create_scalar2vec_interpolator2d(interpolator2d, interp2d_name, domain)
    class(interpolator2d_scalar2vec_t), allocatable, intent(out) :: interpolator2d
    character(len=*),                                intent(in)  :: interp2d_name
    type(domain_t),                                  intent(in)  :: domain

    select case(interp2d_name)
    case("interp2d_p2uv_C_sbp21")
        call create_p2uv_sbp_interpolator(interpolator2d,"W21_stagered_interp_c2i",&
                                         domain,halo_width=1)
    case("interp2d_p2uv_C_sbp42")
        call create_p2uv_sbp_interpolator(interpolator2d,"W42_stagered_interp_c2i",&
                                         domain,halo_width=3)
    case default
        call parcomm_global%abort("create_scalar2vec_interpolator2d, unknown interpolator name: "//&
                                  interp2d_name)
    end select

end subroutine create_scalar2vec_interpolator2d

subroutine create_vec2vec_interpolator2d(interpolator2d, interp2d_name, domain)
    class(interpolator2d_vec2vec_t), allocatable,    intent(out) :: interpolator2d
    character(len=*),                                intent(in)  :: interp2d_name
    type(domain_t),                                  intent(in)  :: domain

    select case(interp2d_name)
    case("interp2d_pvec2uv_C_sbp21")
        call create_pvec2uv_sbp_interpolator(interpolator2d,"W21_stagered_interp_c2i",&
                                         domain,halo_width=1)
    case("interp2d_pvec2uv_C_sbp42")
        call create_pvec2uv_sbp_interpolator(interpolator2d,"W42_stagered_interp_c2i",&
                                         domain,halo_width=3)
    case("interp2d_uv2pvec_C_sbp21")
        call create_uv2pvec_sbp_interpolator(interpolator2d,"W21_stagered_interp_i2c",&
                                             domain,halo_width=1)
    case("interp2d_uv2pvec_C_sbp42")
        call create_uv2pvec_sbp_interpolator(interpolator2d,"W42_stagered_interp_i2c",&
                                             domain,halo_width=3)
    case("interp2d_pvec2uv_Ch_sbp21")
        call create_pvec2uv_Ch_sbp_interpolator(interpolator2d,"W21_stagered_interp_i2c",&
                                                domain,halo_width=1)
    case("interp2d_pvec2uv_Ch_sbp42")
        call create_pvec2uv_Ch_sbp_interpolator(interpolator2d,"W42_stagered_interp_i2c",&
                                                domain,halo_width=3)
    case("interp2d_qvec2uv_Ch_sbp21")
        call create_qvec2uv_Ch_sbp_interpolator(interpolator2d,"W21_stagered_interp_c2i",&
                                                domain,halo_width=1)
    case("interp2d_qvec2uv_Ch_sbp42")
        call create_qvec2uv_Ch_sbp_interpolator(interpolator2d,"W42_stagered_interp_c2i",&
                                                domain,halo_width=3)
    case("interp2d_uv2qvec_Ch_sbp21")
        call create_uv2qvec_Ch_sbp_interpolator(interpolator2d,"W21_stagered_interp_i2c",&
                                                domain,halo_width=1)
    case("interp2d_uv2qvec_Ch_sbp42")
        call create_uv2qvec_Ch_sbp_interpolator(interpolator2d,"W42_stagered_interp_i2c",&
                                                domain,halo_width=3)
    case default
        call parcomm_global%abort("create_vec2vec_interpolator2d, unknown interpolator name: "//&
                                  interp2d_name)
    end select

end subroutine create_vec2vec_interpolator2d

subroutine create_p2uv_sbp_interpolator(interpolator_p2uv, sbp_p2uv_interp_name, &
                                                             domain, halo_width)

    class(interpolator2d_scalar2vec_t), allocatable, intent(out) :: interpolator_p2uv
    character(len=*),                                intent(in)  :: sbp_p2uv_interp_name
    type(domain_t),                                  intent(in)  :: domain
    integer(kind=4),                                 intent(in)  :: halo_width

    type(interpolator_p2uv_sbp_C_t), allocatable :: interpolator_p2uv_sbp

    allocate(interpolator_p2uv_sbp)

    interpolator_p2uv_sbp%sbp_interp_p2uv = create_sbp_operator(sbp_p2uv_interp_name)

    interpolator_p2uv_sbp%exchange = &
          create_o_points_halo_exchange(domain%partition, domain%parcomm, &
                                      domain%topology, halo_width, 'full')

    call move_alloc(interpolator_p2uv_sbp, interpolator_p2uv)
end subroutine create_p2uv_sbp_interpolator

subroutine create_pvec2uv_sbp_interpolator(interpolator_p2uv, sbp_p2uv_interp_name, &
                                                             domain, halo_width)

    class(interpolator2d_vec2vec_t), allocatable, intent(out) :: interpolator_p2uv
    character(len=*),                             intent(in)  :: sbp_p2uv_interp_name
    type(domain_t),                               intent(in)  :: domain
    integer(kind=4),                              intent(in)  :: halo_width

    type(interpolator_pvec2uv_sbp_C_t), allocatable :: interpolator_p2uv_sbp

    allocate(interpolator_p2uv_sbp)

    interpolator_p2uv_sbp%sbp_interp_p2uv = create_sbp_operator(sbp_p2uv_interp_name)

    interpolator_p2uv_sbp%exchange = &
          create_o_points_halo_exchange(domain%partition, domain%parcomm, &
                                      domain%topology, halo_width, 'full')

    call move_alloc(interpolator_p2uv_sbp, interpolator_p2uv)
end subroutine create_pvec2uv_sbp_interpolator

subroutine create_uv2pvec_sbp_interpolator(interpolator_uv2p, sbp_uv2p_interp_name, &
                                                             domain, halo_width)

    class(interpolator2d_vec2vec_t), allocatable, intent(out) :: interpolator_uv2p
    character(len=*),                             intent(in)  :: sbp_uv2p_interp_name
    type(domain_t),                               intent(in)  :: domain
    integer(kind=4),                              intent(in)  :: halo_width

    type(interpolator2d_uv2p_sbp_C_t), allocatable :: interpolator_uv2p_sbp

    allocate(interpolator_uv2p_sbp)

    interpolator_uv2p_sbp%sbp_interp_uv2p = create_sbp_operator(sbp_uv2p_interp_name)

    interpolator_uv2p_sbp%exchange = &
         create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                              domain%topology, halo_width, 'full')

    call move_alloc(interpolator_uv2p_sbp, interpolator_uv2p)
end subroutine create_uv2pvec_sbp_interpolator

subroutine create_pvec2uv_Ch_sbp_interpolator(interpolator_p2uv, sbp_p2uv_interp_name, &
                                                             domain, halo_width)

    class(interpolator2d_vec2vec_t), allocatable, intent(out) :: interpolator_p2uv
    character(len=*),                             intent(in)  :: sbp_p2uv_interp_name
    type(domain_t),                               intent(in)  :: domain
    integer(kind=4),                              intent(in)  :: halo_width

    type(interpolator_pvec2uv_sbp_Ch_t), allocatable :: interpolator_p2uv_sbp

    allocate(interpolator_p2uv_sbp)

    interpolator_p2uv_sbp%sbp_interp_p2uv = create_sbp_operator(sbp_p2uv_interp_name)

    interpolator_p2uv_sbp%exchange = &
          create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                       domain%topology, halo_width, 'full')

    call move_alloc(interpolator_p2uv_sbp, interpolator_p2uv)
end subroutine create_pvec2uv_Ch_sbp_interpolator

subroutine create_qvec2uv_Ch_sbp_interpolator(interpolator_q2uv, sbp_q2uv_interp_name, &
                                                             domain, halo_width)

    class(interpolator2d_vec2vec_t), allocatable, intent(out) :: interpolator_q2uv
    character(len=*),                             intent(in)  :: sbp_q2uv_interp_name
    type(domain_t),                               intent(in)  :: domain
    integer(kind=4),                              intent(in)  :: halo_width

    type(interpolator_q2uv_sbp_Ch_t), allocatable :: interpolator_q2uv_sbp

    allocate(interpolator_q2uv_sbp)

    interpolator_q2uv_sbp%sbp_interp_q2uv = create_sbp_operator(sbp_q2uv_interp_name)

    interpolator_q2uv_sbp%exchange = &
          create_o_points_halo_exchange(domain%partition, domain%parcomm, &
                                       domain%topology, halo_width, 'full')

    call move_alloc(interpolator_q2uv_sbp, interpolator_q2uv)
end subroutine create_qvec2uv_Ch_sbp_interpolator

subroutine create_uv2qvec_Ch_sbp_interpolator(interpolator_uv2q, sbp_uv2q_interp_name, &
                                                             domain, halo_width)

    class(interpolator2d_vec2vec_t), allocatable, intent(out) :: interpolator_uv2q
    character(len=*),                             intent(in)  :: sbp_uv2q_interp_name
    type(domain_t),                               intent(in)  :: domain
    integer(kind=4),                              intent(in)  :: halo_width

    type(interpolator_uv2q_sbp_Ch_t), allocatable :: interpolator_uv2q_sbp

    allocate(interpolator_uv2q_sbp)

    interpolator_uv2q_sbp%sbp_interp_uv2q = create_sbp_operator(sbp_uv2q_interp_name)

    interpolator_uv2q_sbp%exchange = &
          create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                              domain%topology, halo_width, 'full')

    call move_alloc(interpolator_uv2q_sbp, interpolator_uv2q)
end subroutine create_uv2qvec_Ch_sbp_interpolator

end module interpolator2d_factory_mod

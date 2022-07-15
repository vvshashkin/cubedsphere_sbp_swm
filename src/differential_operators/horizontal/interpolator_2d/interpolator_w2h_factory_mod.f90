module interpolator_w2h_factory_mod

use interpolator_w2h_mod,         only : interpolator_w2h_t
use domain_mod,                   only : domain_t
use sbp_factory_mod,              only : create_sbp_operator
use exchange_factory_mod,         only : create_xy_points_halo_exchange
use grid_field_factory_mod,       only : create_grid_field
use interpolator2d_factory_mod,   only : create_vec2vec_interpolator2d
use parcomm_mod,                  only : parcomm_global

implicit none

contains

subroutine create_w2h_interpolator(interpolator_w2h, sbp_i2c_interp_name, domain)

    type(interpolator_w2h_t), intent(out) :: interpolator_w2h
    character(len=*),         intent(in)  :: sbp_i2c_interp_name
    type(domain_t),           intent(in)  :: domain

    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 3

    interpolator_w2h%sbp_interp_w2v = create_sbp_operator(sbp_i2c_interp_name)

    interpolator_w2h%exchange = &
          create_xy_points_halo_exchange(domain%partition, domain%parcomm, &
                                       domain%topology, halo_width, 'full')

    call create_grid_field(interpolator_w2h%fu, halo_width+1, 0, domain%mesh_u)
    call create_grid_field(interpolator_w2h%fv, halo_width+1, 0, domain%mesh_v)
    call create_grid_field(interpolator_w2h%fuh,           0, 0, domain%mesh_p)
    call create_grid_field(interpolator_w2h%fvh,           0, 0, domain%mesh_p)

    select case(sbp_i2c_interp_name)
    case ("W21_stagered_interp_i2c")
        call create_vec2vec_interpolator2d(interpolator_w2h%interp_v2h_op, "interp2d_uv2pvec_C_sbp21", domain)
    case ("W42_stagered_interp_i2c")
        call create_vec2vec_interpolator2d(interpolator_w2h%interp_v2h_op, "interp2d_uv2pvec_C_sbp42", domain)
    case default
        call parcomm_global%abort("create_w2h_interpolator, unknown sbp_i2c_interp_name:"//&
                                 sbp_i2c_interp_name)
    end select

end subroutine create_w2h_interpolator

end module interpolator_w2h_factory_mod

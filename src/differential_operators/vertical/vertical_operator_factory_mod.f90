module vertical_operator_factory_mod

use abstract_vertical_operator_mod, only : vertical_operator_t
use parcomm_mod,                    only : parcomm_global

implicit none

contains

subroutine create_vertical_operator(vertical_op, vertical_op_name)

    use identity_vertical_operator_mod, only : identity_vertical_operator_t

    class(vertical_operator_t), allocatable, intent(out) :: vertical_op
    character(len=*), intent(in) :: vertical_op_name

    select case(vertical_op_name)
    case("eta_diff_p2w_sbp21")
        call create_sbp_vertical_op(vertical_op, "D21_staggered_c2i","w",is_diff=.true.)
    case("eta_diff_p2w_sbp42")
        call create_sbp_vertical_op(vertical_op, "D42_staggered_c2i","w",is_diff=.true.)
    case("eta_diff_sbp21")
        call create_sbp_vertical_op(vertical_op, "d21","w",is_diff=.true.)
    case("eta_diff_sbp42")
        call create_sbp_vertical_op(vertical_op, "d42","w",is_diff=.true.)
    case("eta_diff_sbp63")
        call create_sbp_vertical_op(vertical_op, "d63","w",is_diff=.true.)
    case("eta_diff_w2p_sbp21")
        call create_sbp_vertical_op(vertical_op, "D21_staggered_i2c","p",is_diff=.true.)
    case("eta_diff_w2p_sbp42")
        call create_sbp_vertical_op(vertical_op, "D42_staggered_i2c","p",is_diff=.true.)
    case("vertical_interp_p2w_sbp21")
        call create_sbp_vertical_op(vertical_op, "W21_stagered_interp_c2i","w",is_diff=.false.)
    case("vertical_interp_p2w_sbp42")
        call create_sbp_vertical_op(vertical_op, "W42_stagered_interp_c2i","w",is_diff=.false.)
    case("vertical_interp_w2p_sbp21")
        call create_sbp_vertical_op(vertical_op, "W21_stagered_interp_i2c","p",is_diff=.false.)
    case("vertical_interp_w2p_sbp42")
        call create_sbp_vertical_op(vertical_op, "W42_stagered_interp_i2c","p",is_diff=.false.)
    case("identity_p")
        vertical_op = identity_vertical_operator_t(points_type = "p")
    case default
        call parcomm_global%abort("create_vertical_operator error, unknown operator "//&
                                  vertical_op_name)
    end select
end subroutine create_vertical_operator

subroutine create_sbp_vertical_op(vertical_op, sbp_op_name, target_grid_name, is_diff)

    use sbp_vertical_operator_mod, only : sbp_vertical_op_t
    use sbp_factory_mod,           only : create_sbp_operator

    class(vertical_operator_t), allocatable, intent(out) :: vertical_op
    character(len=*), intent(in) :: sbp_op_name, target_grid_name
    logical,          intent(in) :: is_diff

    type(sbp_vertical_op_t), allocatable :: sbp_vertical_op

    allocate(sbp_vertical_op)
    sbp_vertical_op%sbp_op = create_sbp_operator(sbp_op_name)
    sbp_vertical_op%tagret_grid = target_grid_name
    sbp_vertical_op%is_diff = is_diff

    call move_alloc(sbp_vertical_op, vertical_op)
end subroutine create_sbp_vertical_op

end module vertical_operator_factory_mod

module sbp_factory_mod

use sbp_operator_mod,             only : sbp_operator_t
use parcomm_mod,                  only : parcomm_global

use sbp_operators_collection_mod, only : Q21, lastnonzeroQ21, Da2_in, Da2_inshift, &
                                         Q42, lastnonzeroQ42, Da4_in, Da4_inshift, &
                                         Q43, lastnonzeroQ43,                      &
                                         Q63, lastnonzeroQ63, Da6_in, Da6_inshift, &
                                         W21_staggered_i2c, W21_staggered_i2c_last_nonzero, &
                                         W21_staggered_i2c_in_shift,                        &
                                         W21_staggered_c2i, W21_staggered_c2i_last_nonzero, &
                                         W21_staggered_c2i_in_shift, W21_staggered_in,      &
                                         W42_staggered_i2c_noopt, W42_staggered_i2c_opt1, &
                                         W42_staggered_i2c_opt2, W42_staggered_i2c_last_nonzero, &
                                         W42_staggered_c2i_noopt, W42_staggered_c2i_opt1, &
                                         W42_staggered_c2i_opt2, W42_staggered_c2i_last_nonzero, &
                                         W42_staggered_in, W42_staggered_c2i_in_shift, &
                                         W42_staggered_c2i_in_shift,  W42_staggered_i2c_in_shift, &
                                         D42_staggered_c2i, D42_staggered_c2i_last_nonzero, &
                                         D42_staggered_in, D42_staggered_c2i_in_shift, &
                                         D42_boundary_proj, D42_A_centers, D42_A_interfaces, &
                                         D42_staggered_i2c, D42_staggered_i2c_last_nonzero, &
                                         D42_staggered_i2c_in_shift, &
                                         D42_boundary_proj, D42_A_centers, D42_A_interfaces, &
                                         D21_staggered_c2i, D21_staggered_c2i_last_nonzero, &
                                         D21_staggered_in, D21_staggered_c2i_in_shift, &
                                         D21_boundary_proj, D21_A_centers, D21_A_interfaces, &
                                         D21_staggered_i2c, D21_staggered_i2c_last_nonzero, &
                                         D21_staggered_i2c_in_shift, &
                                         D2_21_in, D2_21_shift, D2_21_edge, lastnonzeroD2_21_edge, &
                                         D2_21_boundary_proj, Q21_A

implicit none

contains

function create_sbp_operator(sbp_operator_name) result(sbp_op)
    character(len=*), intent(in) :: sbp_operator_name
    type(sbp_operator_t)         :: sbp_op

    !the sign of right corner part of matrix
    !should be -1 for derivatives, 1 for interpolations
    real(kind=8)                 :: right_side_sign
    !size of edge matrix block
    integer(kind=4)              :: n_edge, nw_edge, ne
    !

    !Initialization of Left-side and inner matrix coefficients
    !Right side matrix block is initialized below for all options
    !using symmetry
    if(sbp_operator_name == "d21") then
        sbp_op%W_edge_l    = Q21
        sbp_op%edge_last_l = lastnonzeroQ21
        sbp_op%W_in        = Da2_in
        sbp_op%in_shift    = Da2_inshift
        sbp_op%dnx = 0
        right_side_sign = -1.0_8
    else if(sbp_operator_name == "d2_21") then
            sbp_op%W_edge_l    = D2_21_edge
            sbp_op%edge_last_l = lastnonzeroD2_21_edge
            sbp_op%W_in        = D2_21_in
            sbp_op%in_shift    = D2_21_shift
            sbp_op%dnx = 0
            right_side_sign = 1
            sbp_op%proj_operator_l = D2_21_boundary_proj
            sbp_op%Al_in = Q21_A
            sbp_op%Al_out = Q21_A
    else if(sbp_operator_name == "d42") then
        sbp_op%W_edge_l    = Q42
        sbp_op%edge_last_l = lastnonzeroQ42
        sbp_op%W_in        = Da4_in
        sbp_op%in_shift    = Da4_inshift
        sbp_op%dnx = 0
        right_side_sign = -1.0_8
    else if(sbp_operator_name == "d43") then
        sbp_op%W_edge_l    = Q43
        sbp_op%edge_last_l = lastnonzeroQ43
        sbp_op%W_in        = Da4_in
        sbp_op%in_shift    = Da4_inshift
        sbp_op%dnx = 0
        right_side_sign = -1.0_8
    else if(sbp_operator_name == "d63") then
        sbp_op%W_edge_l    = Q63
        sbp_op%edge_last_l = lastnonzeroQ63
        sbp_op%W_in        = Da6_in
        sbp_op%in_shift    = Da6_inshift
        sbp_op%dnx = 0
        right_side_sign = -1.0_8
    else if(sbp_operator_name == "W21_stagered_interp_c2i") then
        sbp_op%W_edge_l    = W21_staggered_c2i
        sbp_op%edge_last_l = W21_staggered_c2i_last_nonzero
        sbp_op%W_in        = W21_staggered_in
        sbp_op%in_shift    = W21_staggered_c2i_in_shift
        sbp_op%dnx = 1
        right_side_sign = 1.0_8
    else if(sbp_operator_name == "W21_stagered_interp_i2c") then
        sbp_op%W_edge_l    = W21_staggered_i2c
        sbp_op%edge_last_l = W21_staggered_i2c_last_nonzero
        sbp_op%W_in        = W21_staggered_in
        sbp_op%in_shift    = W21_staggered_i2c_in_shift
        sbp_op%dnx = -1
        right_side_sign = 1.0_8
    else if(sbp_operator_name == "W42_stagered_interp_c2i") then
        sbp_op%W_edge_l    = W42_staggered_c2i_opt1
        sbp_op%edge_last_l = W42_staggered_c2i_last_nonzero
        sbp_op%W_in        = W42_staggered_in
        sbp_op%in_shift    = W42_staggered_c2i_in_shift
        sbp_op%dnx = 1
        right_side_sign = 1.0_8
    else if(sbp_operator_name == "W42_stagered_interp_i2c") then
        sbp_op%W_edge_l    = W42_staggered_i2c_opt1
        sbp_op%edge_last_l = W42_staggered_i2c_last_nonzero
        sbp_op%W_in        = W42_staggered_in
        sbp_op%in_shift    = W42_staggered_i2c_in_shift
        sbp_op%dnx = -1
        right_side_sign = 1.0_8
    else if(sbp_operator_name == "D21_staggered_c2i") then
        sbp_op%W_edge_l    = D21_staggered_c2i
        sbp_op%edge_last_l = D21_staggered_c2i_last_nonzero
        sbp_op%W_in        = D21_staggered_in
        sbp_op%in_shift    = D21_staggered_c2i_in_shift
        sbp_op%dnx = 1
        right_side_sign =-1.0_8
        sbp_op%proj_operator_l = D21_boundary_proj
        sbp_op%Al_in = D21_A_centers
        sbp_op%Al_out = D21_A_interfaces
    else if(sbp_operator_name == "D21_staggered_i2c") then
        sbp_op%W_edge_l    = D21_staggered_i2c
        sbp_op%edge_last_l = D21_staggered_i2c_last_nonzero
        sbp_op%W_in        = D21_staggered_in
        sbp_op%in_shift    = D21_staggered_i2c_in_shift
        sbp_op%dnx =-1
        right_side_sign =-1.0_8
        sbp_op%proj_operator_l = D21_boundary_proj
        sbp_op%Al_in = D21_A_interfaces
        sbp_op%Al_out = D21_A_centers
    else if(sbp_operator_name == "D42_staggered_c2i") then
        sbp_op%W_edge_l    = D42_staggered_c2i
        sbp_op%edge_last_l = D42_staggered_c2i_last_nonzero
        sbp_op%W_in        = D42_staggered_in
        sbp_op%in_shift    = D42_staggered_c2i_in_shift
        sbp_op%dnx = 1
        right_side_sign =-1.0_8
        sbp_op%proj_operator_l = D42_boundary_proj
        sbp_op%Al_in = D42_A_centers
        sbp_op%Al_out = D42_A_interfaces
    else if(sbp_operator_name == "D42_staggered_i2c") then
        sbp_op%W_edge_l    = D42_staggered_i2c
        sbp_op%edge_last_l = D42_staggered_i2c_last_nonzero
        sbp_op%W_in        = D42_staggered_in
        sbp_op%in_shift    = D42_staggered_i2c_in_shift
        sbp_op%dnx =-1
        right_side_sign =-1.0_8
        sbp_op%proj_operator_l = D42_boundary_proj
        sbp_op%Al_in = D42_A_interfaces
        sbp_op%Al_out = D42_A_centers
    else
        call parcomm_global%abort("sbp_factory_mod, unknown sbp operator name: "// sbp_operator_name)
    end if

    nw_edge = size(sbp_op%W_edge_l,1)
    n_edge  = size(sbp_op%W_edge_l,2)
    sbp_op%W_edge_r = right_side_sign * sbp_op%W_edge_l(nw_edge:1:-1,n_edge:1:-1)
    sbp_op%edge_first_shift_r = -sbp_op%edge_last_l(n_edge:1:-1)+1
    if(allocated(sbp_op%Al_in)) then
        ne = size(sbp_op%Al_in)
        sbp_op%Ar_in = sbp_op%Al_in(ne:1:-1)
    end if
    if(allocated(sbp_op%Al_out)) then
        ne = size(sbp_op%Al_out)
        sbp_op%Ar_out = sbp_op%Al_out(ne:1:-1)
    end if
    if(allocated(sbp_op%proj_operator_l)) then
        ne = size(sbp_op%proj_operator_l)
        sbp_op%proj_operator_r = sbp_op%proj_operator_l(ne:1:-1)
    end if

end function create_sbp_operator

end module sbp_factory_mod

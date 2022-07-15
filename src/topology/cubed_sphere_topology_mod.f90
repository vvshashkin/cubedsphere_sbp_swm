module cubed_sphere_topology_mod

use topology_mod, only: topology_t

type, extends(topology_t) :: cubed_sphere_topology_t

contains
    procedure :: init => init_cs_topology
    procedure :: transform_tile_coords
    procedure :: transform_index
    procedure :: find_basis_orientation
    procedure :: get_true_panel_coords
end type cubed_sphere_topology_t

!!             ______
!!            |      :
!!            |  6   :
!!            |......:
!!   ......    ......   ......    ......
!!  :      |  :      | :      |  :      |
!!  :   4  |  :   1  | :   2  |  :   3  |
!!  :______|  :______| :______|  :______|
!!             ______
!!            |      :
!!            |  5   :
!!            |...,,,:
!!

integer(kind=4), parameter :: ex(3,6) = reshape( (/  (/-1, 0,  0/), &
                                                     (/ 0, 0, -1/), &
                                                     (/ 1, 0,  0/), &
                                                     (/ 0, 0,  1/), &
                                                     (/ 1, 0,  0/), &
                                                     (/ 1, 0,  0/)  /),  (/3 ,6/) )

integer(kind=4), parameter :: ey(3,6) = reshape( (/  (/0, 1, 0/), &
                                                     (/0, 1, 0/), &
                                                     (/0, 1, 0/), &
                                                     (/0, 1, 0/), &
                                                     (/0, 0, 1/), &
                                                     (/0, 0,-1/)    /),  (/3 ,6/) )

integer(kind=4), parameter :: n(3,6) = reshape( (/   (/ 0, 0,  1/), &
                                                     (/-1, 0,  0/), &
                                                     (/ 0, 0, -1/), &
                                                     (/ 1, 0,  0/), &
                                                     (/ 0, 1,  0/), &
                                                     (/ 0,-1,  0/)    /),  (/3 ,6/) )

integer(kind=4), parameter ::  r(3,6) = reshape( (/  (/1, 0, 0/), &
                                                     (/1, 0, 1/), &
                                                     (/0, 0, 1/), &
                                                     (/0, 0, 0/), &
                                                     (/0, 0, 0/), &
                                                     (/0, 1, 1/)    /),  (/3 ,6/) )
contains

subroutine init_cs_topology(this)
    class(cubed_sphere_topology_t), intent(inout) :: this

    this%npanels = 6
    this%ex = ex
    this%ey = ey
    this%n  = n
    this%r  = r

end

subroutine calc_xyz_cords(panel_ind, i, j, npoints, ix, iy, iz)

    integer(kind=4), intent(in)  :: panel_ind, i, j, npoints
    integer(kind=4), intent(out) :: ix, iy, iz

    ix = (i-1)*ex(1, panel_ind) + (j-1)*ey(1, panel_ind) + (npoints-1)*r(1,panel_ind) + 1
    iy = (i-1)*ex(2, panel_ind) + (j-1)*ey(2, panel_ind) + (npoints-1)*r(2,panel_ind) + 1
    iz = (i-1)*ex(3, panel_ind) + (j-1)*ey(3, panel_ind) + (npoints-1)*r(3,panel_ind) + 1

end subroutine calc_xyz_cords

subroutine transform_index(this, pn_source, i_source, j_source, &
                                 pn_target, i_target, j_target, &
                                 Nx, Ny)

    !transform indices at source panel to indices at target panel

    class(cubed_sphere_topology_t), intent(in) :: this

    integer(kind=4),  intent(in)  :: pn_source, pn_target ! index of source and target panels
    integer(kind=4),  intent(in)  :: i_source , j_source
    integer(kind=4),  intent(out) :: i_target , j_target
    integer(kind=4),  intent(in)  :: nx, ny !number of points along dimensions

    integer(kind=4) :: ex_local(1:3), ey_local(1:3), r_local(1:3), shift(2)
    integer(kind=4) :: ex_ex, ex_ey, ex_r, ey_ex, ey_ey, ey_r
    logical :: is_opposite, is_right, is_left, is_bottom, is_top


    if (pn_source == pn_target) then
        i_target = i_source
        j_target = j_source
        return
    end if

    is_right    = (dot_product(this%ex(1:3,pn_target), this%n(1:3,pn_source)) == -1)
    is_left     = (dot_product(this%ex(1:3,pn_target), this%n(1:3,pn_source)) ==  1)
    is_bottom   = (dot_product(this%ey(1:3,pn_target), this%n(1:3,pn_source)) ==  1)
    is_top      = (dot_product(this%ey(1:3,pn_target), this%n(1:3,pn_source)) == -1)
    is_opposite = (dot_product( this%n(1:3,pn_target), this%n(1:3,pn_source)) == -1)

    if (is_right) then
        ex_local = this%n (1:3, pn_target)
        ey_local = this%ey(1:3, pn_target)
        r_local  = this%r (1:3, pn_target) + this%ex(1:3, pn_target)
        shift = [nx,0]
    else if (is_left) then
        ex_local = -this%n (1:3, pn_target)
        ey_local =  this%ey(1:3, pn_target)
        r_local  =  this%r (1:3, pn_target) - ex_local(1:3)
        shift = [-nx,0]
    else if (is_bottom) then
        ex_local = this%ex(1:3, pn_target)
        ey_local = -this%n(1:3, pn_target)
        r_local  =  this%r(1:3, pn_target) + this%n(1:3, pn_target)
        shift = [0,-ny]
    else if (is_top) then
        ex_local =  this%ex(1:3, pn_target)
        ey_local =   this%n(1:3, pn_target)
        r_local  =   this%r(1:3, pn_target) +this%ey(1:3,pn_target)
        shift = [0, ny]
    else if (is_opposite) then
        ex_local = -this%ex(1:3, pn_target)
        ey_local =  this%ey(1:3, pn_target)
        r_local  =  this%r (1:3, pn_target) + this%ex(1:3, pn_target) + this%n(1:3, pn_target)
        shift = [2*nx,0]
    end if

    ex_ex = dot_product(this%ex(1:3, pn_source), ex_local(1:3))
    ex_ey = dot_product(this%ey(1:3, pn_source), ex_local(1:3))
    ex_r  = dot_product( this%r(1:3, pn_source) - r_local(1:3), ex_local(1:3))

    ey_ex = dot_product(this%ex(1:3, pn_source), ey_local(1:3))
    ey_ey = dot_product(this%ey(1:3, pn_source), ey_local(1:3))
    ey_r  = dot_product( this%r(1:3, pn_source) - r_local(1:3), ey_local(1:3))

    i_target = (i_source-1)*ex_ex + (j_source-1)*ex_ey + (nx-1)*ex_r + 1 + shift(1)
    j_target = (i_source-1)*ey_ex + (j_source-1)*ey_ey + (ny-1)*ey_r + 1 + shift(2)

end subroutine transform_index
subroutine transform_tile_coords(this, pn_source, source_tile, &
                                       pn_target, target_tile, &
                                       nx, ny)

    use tile_mod, only : tile_t

    !transform tile coords at source panel to coords at target panel
    class(cubed_sphere_topology_t), intent(in) :: this

    integer(kind=4), intent(in)  :: pn_source, pn_target ! index of source and target panels
    type(tile_t),    intent(in)  :: source_tile ! tile coords at source panel
    type(tile_t),    intent(out) :: target_tile ! tile coords at target panel
    integer(kind=4), intent(in)  :: nx, ny !number of points along dimensions

    integer(kind=4) is, ie, js, je

    call this%transform_index(pn_source, source_tile%is, source_tile%js, &
                              pn_target, is, js, Nx, Ny)

    call this%transform_index(pn_source, source_tile%ie, source_tile%je, &
                              pn_target, ie, je, Nx, Ny)

    target_tile%is = min(is,ie)
    target_tile%ie = max(is,ie)
    target_tile%js = min(js,je)
    target_tile%je = max(js,je)

    target_tile%ks = source_tile%ks
    target_tile%ke = source_tile%ke

end subroutine transform_tile_coords
subroutine find_basis_orientation(this, pn_source, pn_target, i_step, j_step, first_dim_index )

    !transform tile coords at source panel to coords at target panel
    class(cubed_sphere_topology_t), intent(in) :: this
    integer(kind=4), intent(in)  :: pn_source, pn_target ! index of source and target panels

    integer(kind=4),  intent(out), optional :: i_step, j_step  !Determines descending or ascending order of loop
    character(len=1), intent(out), optional :: first_dim_index !Determines first loop index in horizontal direction

    integer(kind=4) :: ex_local(1:3), ey_local(1:3)
    integer(kind=4) :: ex_ex, ex_ey, ey_ex, ey_ey
    logical :: is_opposite, is_right, is_left, is_bottom, is_top

    if (pn_source == pn_target) then
        first_dim_index = 'i'
        i_step = 1
        j_step = 1
        return
    end if

    is_right    = (dot_product(this%ex(1:3,pn_target), this%n(1:3,pn_source)) == -1)
    is_left     = (dot_product(this%ex(1:3,pn_target), this%n(1:3,pn_source)) ==  1)
    is_bottom   = (dot_product(this%ey(1:3,pn_target), this%n(1:3,pn_source)) ==  1)
    is_top      = (dot_product(this%ey(1:3,pn_target), this%n(1:3,pn_source)) == -1)
    is_opposite = (dot_product( this%n(1:3,pn_target), this%n(1:3,pn_source)) == -1)

    if (is_right) then
        ex_local = this%n (1:3, pn_target)
        ey_local = this%ey(1:3, pn_target)
    else if (is_left) then
        ex_local = -this%n (1:3, pn_target)
        ey_local =  this%ey(1:3, pn_target)
    else if (is_bottom) then
        ex_local =   this%ex(1:3, pn_target)
        ey_local =  -this%n (1:3, pn_target)
    else if (is_top) then
        ex_local =  this%ex(1:3, pn_target)
        ey_local =   this%n(1:3, pn_target)
    else if (is_opposite) then
        ex_local = -this%ex(1:3, pn_target)
        ey_local =  this%ey(1:3, pn_target)
    end if

    ex_ex = dot_product(this%ex(1:3, pn_source), ex_local(1:3))
    ex_ey = dot_product(this%ey(1:3, pn_source), ex_local(1:3))

    ey_ex = dot_product(this%ex(1:3, pn_source), ey_local(1:3))
    ey_ey = dot_product(this%ey(1:3, pn_source), ey_local(1:3))

    if ( ex_ex /= 0 ) then
        first_dim_index = 'i'
        i_step = ex_ex
        j_step = ey_ey
    else
        first_dim_index = 'j'
        i_step = ex_ey
        j_step = ey_ex
    end if

end subroutine find_basis_orientation
subroutine get_true_panel_coords(this, i_in, j_in, pn_in, nx, ny, i_out, j_out, pn_out)

    class(cubed_sphere_topology_t), intent(in) :: this

    integer(kind=4), intent(in)  :: i_in,  j_in,  pn_in, nx, ny
    integer(kind=4), intent(out) :: i_out, j_out, pn_out

    logical :: is_x_left, is_x_right, is_x_center, &
               is_y_top, is_y_bottom, is_y_center, &
               is_right, is_left, is_top, is_bottom, is_center

   integer(kind=4)  :: i_step, j_step
   character(len=1) :: first_dim_index

    is_x_left   = (i_in >= 1-nx .and. i_in <= 0)
    is_x_center = (i_in >= 1    .and. i_in <= nx)
    is_x_right  = (i_in >= 1+nx .and. i_in <= 2*nx)

    is_y_bottom = (j_in >= 1-ny .and. j_in <= 0)
    is_y_center = (j_in >= 1    .and. j_in <= ny)
    is_y_top    = (j_in >= 1+ny .and. j_in <= 2*ny)

    is_right  = (is_x_right  .and. is_y_center)
    is_left   = (is_x_left   .and. is_y_center)
    is_top    = (is_x_center .and. is_y_top   )
    is_bottom = (is_x_center .and. is_y_bottom)
    is_center = (is_x_center .and. is_y_center)

    do pn_out = 1, 6
        if (is_right  .and. dot_product(this%ex(1:3,pn_in), this%n(1:3,pn_out)) == -1) exit
        if (is_left   .and. dot_product(this%ex(1:3,pn_in), this%n(1:3,pn_out)) ==  1) exit
        if (is_bottom .and. dot_product(this%ey(1:3,pn_in), this%n(1:3,pn_out)) ==  1) exit
        if (is_top    .and. dot_product(this%ey(1:3,pn_in), this%n(1:3,pn_out)) == -1) exit
        if (is_center .and. pn_in == pn_out) exit
    end do

    if (pn_out == 7) then
        print*, 'Something wrong in get_true_panel_coords! Stop!'
        stop
    end if

    call this%find_basis_orientation(pn_in, pn_out, i_step, j_step, first_dim_index)

    if (first_dim_index=='i') then
        call this%transform_index(pn_in, i_in, j_in, pn_out, i_out, j_out, nx, ny)
    else
        call this%transform_index(pn_in, i_in, j_in, pn_out, i_out, j_out, ny, nx)
    end if
end subroutine get_true_panel_coords
end module cubed_sphere_topology_mod

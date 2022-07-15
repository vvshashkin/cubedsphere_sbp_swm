module topology_mod

use tile_mod, only: tile_t

implicit none

type, abstract :: topology_t
    integer(kind=4) npanels
    integer(kind=4), allocatable :: ex(:,:)
    integer(kind=4), allocatable :: ey(:,:)
    integer(kind=4), allocatable :: n(:,:)
    integer(kind=4), allocatable :: r(:,:)

contains
    procedure(init_topology_i),          deferred :: init
    procedure(transform_tile_coords_i),  deferred :: transform_tile_coords
    procedure(transform_index_i),        deferred :: transform_index
    procedure(find_basis_orientation_i), deferred :: find_basis_orientation
    procedure(get_true_panel_coords_i),  deferred :: get_true_panel_coords
end type topology_t

abstract interface
    subroutine init_topology_i(this)
        import topology_t
        class(topology_t), intent(inout) :: this
    end subroutine init_topology_i
    subroutine transform_tile_coords_i(this, pn_source, source_tile, &
                                             pn_target, target_tile, &
                                             nx, ny)
        import topology_t, tile_t
        !transform tile coords at source panel to coords at target panel
        class(topology_t), intent(in)  :: this
        integer(kind=4),   intent(in)  :: pn_source, pn_target ! index of source and target panels
        type(tile_t),      intent(in)  :: source_tile ! tile coords at source panel
        type(tile_t),      intent(out) :: target_tile ! tile coords at target panel

        integer(kind=4),  intent(in)            :: nx, ny !number of points along dimensions
    end subroutine transform_tile_coords_i

    subroutine transform_index_i(this, pn_source, i_source, j_source, &
                                       pn_target, i_target, j_target, &
                                       Nx, Ny)
        import topology_t
        class(topology_t), intent(in)            :: this
        integer(kind=4),  intent(in)  :: pn_source, pn_target ! index of source and target panels
        integer(kind=4),  intent(in)  :: i_source , j_source
        integer(kind=4),  intent(out) :: i_target , j_target
        integer(kind=4),  intent(in)  :: nx, ny !number of points along dimensions

    end subroutine transform_index_i

    subroutine find_basis_orientation_i(this, pn_source, pn_target, i_step, j_step, first_dim_index )
        import topology_t
        import tile_t
        class(topology_t), intent(in)            :: this
        integer(kind=4),   intent(in)            :: pn_source, pn_target ! index of source and target panels
        integer(kind=4),   intent(out), optional :: i_step, j_step  !Determines descending or ascending order of loop
        character(len=1),  intent(out), optional :: first_dim_index !Determines first loop index in horizontal direction
    end subroutine find_basis_orientation_i
    subroutine get_true_panel_coords_i(this, i_in, j_in, pn_in, nx, ny, i_out, j_out, pn_out)

        import topology_t
        class(topology_t), intent(in) :: this
        integer(kind=4), intent(in)   :: i_in,  j_in,  pn_in, nx, ny
        integer(kind=4), intent(out)  :: i_out, j_out, pn_out

    end subroutine get_true_panel_coords_i

end interface

!!             ......
!!            |      :
!!            |  6   :
!!            |______:
!!   ......    ......   ......    ......
!!  |      :  |      : |      :  |      :
!!  |   4  :  |   1  : |   2  :  |   3  :
!!  |______:  |______: |______:  |______:
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
!INNER normal to cube face
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

subroutine calc_xyz_cords(panel_ind, i, j, npoints, ix, iy, iz)

    integer(kind=4), intent(in)  :: panel_ind, i, j, npoints
    integer(kind=4), intent(out) :: ix, iy, iz

    ix = (i-1)*ex(1, panel_ind) + (j-1)*ey(1, panel_ind) + (npoints-1)*r(1,panel_ind) + 1
    iy = (i-1)*ex(2, panel_ind) + (j-1)*ey(2, panel_ind) + (npoints-1)*r(2,panel_ind) + 1
    iz = (i-1)*ex(3, panel_ind) + (j-1)*ey(3, panel_ind) + (npoints-1)*r(3,panel_ind) + 1

end subroutine calc_xyz_cords

end module topology_mod

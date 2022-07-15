module tiles_mod

use parcomm_mod, only : parcomm_global
use tile_mod,    only : tile_t

implicit none

type, public :: tiles_t
    integer(kind=4) :: Nx, Ny, Nz, Nt
    type(tile_t), allocatable :: tile(:)
contains
    procedure, public :: init
end type tiles_t

contains

subroutine init(this, Nt, Nz, Ny, Nx, tile)
    class(tiles_t),  intent(inout) :: this
    integer(kind=4), intent(in)    :: Nt, Nx, Ny, Nz
    type(tile_t), optional, intent(in) :: tile(1:Nt)

    this%Nt = Nt
    this%Nz = Nz
    this%Ny = Ny
    this%Nx = Nx

    allocate(this%tile(1:Nt))

    if (present(tile)) then
        this%tile(1:Nt) = tile(1:Nt)
    end if

end subroutine init

end module tiles_mod

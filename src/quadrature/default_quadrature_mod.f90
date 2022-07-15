module default_quadrature_mod

use parcomm_mod,    only : parcomm_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t
use abstract_quadrature_mod, only : tile_quadrature_t

use mpi

implicit none

type, extends(tile_quadrature_t) :: default_tile_quadrature_t
    contains
        procedure, public :: mass => calc_tile_mass_default
        procedure, public :: dot => calc_tile_dot_default
end type default_tile_quadrature_t

contains

function calc_tile_mass_default(this, f, mesh) result(out)

    class(default_tile_quadrature_t), intent(in) :: this
    type(tile_field_t), intent(in) :: f
    type(tile_mesh_t),  intent(in) :: mesh
    real(kind=8)                   :: out

    integer(kind=4) :: k, j, i

    out = 0.0_8

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                out = out + f%p(i,j,k)*mesh%J(i,j,k)*mesh%hx*mesh%hy
            end do
        end do
    end do

end function

function calc_tile_dot_default(this, f1, f2, mesh) result(out)

    class(default_tile_quadrature_t), intent(in) :: this
    type(tile_field_t), intent(in) :: f1, f2
    type(tile_mesh_t),  intent(in) :: mesh
    real(kind=8)                   :: out

    integer(kind=4) :: k, j, i

    out = 0.0_8

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                out = out + f1%p(i,j,k)*f2%p(i,j,k)*mesh%J(i,j,k)*mesh%hx*mesh%hy
            end do
        end do
    end do

end function

end module default_quadrature_mod

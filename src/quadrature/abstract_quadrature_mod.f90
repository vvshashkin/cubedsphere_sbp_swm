module abstract_quadrature_mod

use parcomm_mod,    only : parcomm_global, parcomm_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t
use parcomm_mod,    only : parcomm_t

use mpi

implicit none

type :: quadrature_t
    class(tile_quadrature_t), allocatable :: tile(:)
    contains
        procedure, public :: mass => calc_mass
        procedure, public :: dot => calc_dot
end type quadrature_t

type, abstract :: tile_quadrature_t
    contains
        procedure, public :: mass => calc_tile_mass_not_implemented
        procedure, public :: dot => calc_tile_dot_not_implemented
end type tile_quadrature_t

contains

function calc_mass(this, f, mesh, parcomm) result(mass)
    class(quadrature_t),         intent(in) :: this
    type(grid_field_t),          intent(in) :: f
    type(mesh_t),                intent(in) :: mesh
    type(parcomm_t),             intent(in) :: parcomm

    real(kind=8) :: mass

    real(kind=8) :: mass_loc
    integer(kind=4) :: t, err

    mass_loc = 0.0_8

    do t = mesh%ts, mesh%te
        mass_loc = mass_loc + this%tile(t)%mass(f%tile(t), mesh%tile(t))
    end do
    call mpi_allreduce(mass_loc, mass, 1, mpi_double, mpi_sum, parcomm%comm_w, err)
end

function calc_dot(this, f1, f2, mesh, parcomm) result(dot_prod)
    class(quadrature_t),         intent(in) :: this
    type(grid_field_t),          intent(in) :: f1, f2
    type(mesh_t),                intent(in) :: mesh
    type(parcomm_t),             intent(in) :: parcomm

    real(kind=8) :: dot_prod

    real(kind=8) :: dot_prod_loc
    integer(kind=4) :: t, err

    dot_prod_loc = 0.0_8

    do t = mesh%ts, mesh%te
        dot_prod_loc = dot_prod_loc + this%tile(t)%dot(f1%tile(t), f2%tile(t), mesh%tile(t))
    end do
    call mpi_allreduce(dot_prod_loc, dot_prod, 1, mpi_double, mpi_sum, parcomm%comm_w, err)
end

function calc_tile_mass_not_implemented(this, f, mesh) result(mass)
    class(tile_quadrature_t), intent(in) :: this
    type(tile_field_t),       intent(in) :: f
    type(tile_mesh_t),        intent(in) :: mesh

    real(kind=8) :: mass

    mass = 0.0_8

    call parcomm_global%abort("tried to use not implemented function mass of quadrature, subtype of mesh_quadrature_t class")

end

function calc_tile_dot_not_implemented(this, f1, f2, mesh) result(mass)
    class(tile_quadrature_t), intent(in) :: this
    type(tile_field_t),       intent(in) :: f1, f2
    type(tile_mesh_t),        intent(in) :: mesh

    real(kind=8) :: mass

    mass = 0.0_8

    call parcomm_global%abort("tried to use not implemented function dot of quadrature, subtype of mesh_quadrature_t class")

end

end module abstract_quadrature_mod

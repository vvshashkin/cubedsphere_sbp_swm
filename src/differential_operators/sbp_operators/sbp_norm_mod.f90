module sbp_norm_mod

use grid_field_mod,         only : grid_field_t, tile_field_t
use grid_field_factory_mod, only : create_grid_field
use parcomm_mod,            only : parcomm_t
use mesh_mod,               only : mesh_t, tile_mesh_t
use mpi

implicit none

contains

function calc_mass(f, Ax, Ay, mesh, parcomm) result(out)

    type(grid_field_t),           intent(in)    :: f
    real(kind=8),                 intent(in)    :: Ax(:), Ay(:)
    type(mesh_t),                 intent(in)    :: mesh
    type(parcomm_t),              intent(in)    :: parcomm

    real(kind=8) :: out
    real(kind=8) :: out_loc

    type(grid_field_t) :: mass_mtrx

    integer(kind=4) :: t, mesh_nx, mesh_ny, err

    out = 0.0_8
    out_loc = 0.0_8
    call create_grid_field(mass_mtrx, 0, 0, mesh)

    mesh_nx = -huge(1)
    mesh_ny = -huge(1)

    do t = mesh%ts, mesh%te
        mesh_nx = max(mesh_nx, mesh%tile(t)%ie)
        mesh_ny = max(mesh_ny, mesh%tile(t)%je)
    end do

    do t = mesh%ts, mesh%te
        call calc_mass_matrix_tile(mass_mtrx%tile(t), Ax, Ay, mesh%tile(t), mesh_nx, mesh_ny)
    end do

    do t = mesh%ts, mesh%te
        out_loc = out_loc + calc_mass_tile(f%tile(t), mass_mtrx%tile(t), mesh%tile(t))
    end do

    call mpi_allreduce(out_loc, out, 1, mpi_double, mpi_sum, parcomm%comm_w, err)

end function calc_mass

subroutine calc_mass_matrix_tile(mass_mtrx, Ax, Ay, mesh, mesh_nx, mesh_ny)

    type(tile_field_t), intent(inout) :: mass_mtrx
    type(tile_mesh_t),  intent(in)    :: mesh
    integer(kind=4),    intent(in)    :: mesh_nx, mesh_ny
    real(kind=8),       intent(in)    :: Ax(:), Ay(:)


    integer(kind=4) :: nx, ny, i, j
    real(kind=8) :: wi, wj


    nx = size(Ax, dim=1)
    ny = size(Ay, dim=1)

    do j = mesh%js, mesh%je
        wj = 1.0_8
        if (j<=ny .and. j>=1) wj = Ay(j)
        if (j<=mesh_ny .and. j>=mesh_ny-ny+1) wj = Ay(mesh_ny-j+1)
        do i = mesh%is, mesh%ie
            wi = 1.0_8
            if (i<=nx .and. i>=1) wi = Ax(i)
            if (i<=mesh_nx .and. i>=mesh_nx-nx+1) wi = Ax(mesh_nx-i+1)

            mass_mtrx%p(i,j,1) = wi*wj

        end do
    end do


end subroutine calc_mass_matrix_tile
function calc_mass_tile(f, mass_mtrx, mesh) result(out)

    type(tile_field_t), intent(in) :: f, mass_mtrx
    type(tile_mesh_t),  intent(in) :: mesh
    real(kind=8)                   :: out

    integer(kind=4) :: k, j, i

    out = 0.0_8

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                out = out + mass_mtrx%p(i,j,1)*f%p(i,j,k)*mesh%G(i,j)*mesh%hx*mesh%hy
            end do
        end do
    end do

end function calc_mass_tile
end module sbp_norm_mod

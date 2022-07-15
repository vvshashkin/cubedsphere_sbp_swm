module grid_field_factory_mod

use grid_field_mod, only : grid_field_t

implicit none

contains

subroutine create_grid_field(grid_field, halo_width_xy, halo_width_z, mesh)

    use mesh_mod,       only : mesh_t

    type(grid_field_t), intent(out) :: grid_field
    integer(kind=4),    intent(in)  :: halo_width_xy, halo_width_z

    type(mesh_t),  intent(in)  :: mesh

    integer(kind=4) :: t

    allocate(grid_field%tile(mesh%ts:mesh%te))

    grid_field%ts = mesh%ts
    grid_field%te = mesh%te

    do t = mesh%ts, mesh%te
        call grid_field%tile(t)%init(mesh%tile(t)%is-halo_width_xy, mesh%tile(t)%ie+halo_width_xy, &
                                     mesh%tile(t)%js-halo_width_xy, mesh%tile(t)%je+halo_width_xy, &
                                     mesh%tile(t)%ks-halo_width_z,  mesh%tile(t)%ke+halo_width_z)
    end do

end subroutine create_grid_field
subroutine create_grid_field_global(grid_field, halo_width_xy, halo_width_z, tiles)

    use tiles_mod, only : tiles_t

    type(grid_field_t), intent(out) :: grid_field
    integer(kind=4),    intent(in)  :: halo_width_xy, halo_width_z
    type(tiles_t),      intent(in)  :: tiles

    integer(kind=4) :: t, ts, te
    integer(kind=4) :: is, ie, js, je, ks, ke

    allocate(grid_field%tile(1:tiles%Nt))

    grid_field%ts = 1
    grid_field%te = tiles%Nt

    do t = 1, tiles%Nt
        call tiles%tile(t)%getind(is, ie, js, je, ks, ke)
        call grid_field%tile(t)%init(is-halo_width_xy, ie+halo_width_xy, &
                                     js-halo_width_xy, je+halo_width_xy, &
                                     ks-halo_width_z,  ke+halo_width_z)
    end do

end subroutine create_grid_field_global
end module grid_field_factory_mod

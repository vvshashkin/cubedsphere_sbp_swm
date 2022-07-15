module stvec_advection_factory_mod

implicit none

contains

subroutine create_stvec_advection(stvec, domain, halo_width_xy, halo_width_z)

    use domain_mod,             only : domain_t
    use stvec_mod,              only : stvec_t
    use stvec_advection_mod,    only : stvec_advection_t
    use grid_field_factory_mod, only : create_grid_field

    class(stvec_t), allocatable, intent(out) :: stvec
    type(domain_t),              intent(in)  :: domain
    integer(kind=4),             intent(in)  :: halo_width_xy, halo_width_z

    type(stvec_advection_t), allocatable :: stvec_advection

    allocate(stvec_advection)

    call create_grid_field(stvec_advection%h, halo_width_xy, halo_width_z, domain%mesh_p)
    call create_grid_field(stvec_advection%u, halo_width_xy, halo_width_z, domain%mesh_u)
    call create_grid_field(stvec_advection%v, halo_width_xy, halo_width_z, domain%mesh_v)

    call move_alloc(stvec_advection, stvec)

end subroutine create_stvec_advection

end module stvec_advection_factory_mod

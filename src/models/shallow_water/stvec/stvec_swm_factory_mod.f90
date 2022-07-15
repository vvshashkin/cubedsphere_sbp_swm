module stvec_swm_factory_mod

implicit none

contains

subroutine create_stvec_swm(stvec, domain, halo_width_xy, halo_width_z)

    use domain_mod,             only : domain_t
    use stvec_mod,              only : stvec_t
    use stvec_swm_mod,          only : stvec_swm_t
    use grid_field_factory_mod, only : create_grid_field

    class(stvec_t), allocatable, intent(out) :: stvec
    type(domain_t),              intent(in)  :: domain
    integer(kind=4),             intent(in)  :: halo_width_xy, halo_width_z

    type(stvec_swm_t), allocatable :: stvec_swm

    allocate(stvec_swm)

    call create_grid_field(stvec_swm%h, halo_width_xy, halo_width_z, domain%mesh_p)
    call create_grid_field(stvec_swm%u, halo_width_xy, halo_width_z, domain%mesh_u)
    call create_grid_field(stvec_swm%v, halo_width_xy, halo_width_z, domain%mesh_v)

    call move_alloc(stvec_swm, stvec)

end subroutine create_stvec_swm

end module stvec_swm_factory_mod

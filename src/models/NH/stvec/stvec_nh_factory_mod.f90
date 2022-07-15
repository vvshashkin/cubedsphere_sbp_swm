module stvec_nh_factory_mod

implicit none

contains

subroutine create_stvec_nh(stvec, domain, halo_width_xy, halo_width_z)

    use domain_mod,             only : domain_t
    use stvec_mod,              only : stvec_t
    use stvec_nh_mod,           only : stvec_nh_t
    use grid_field_factory_mod, only : create_grid_field

    class(stvec_t), allocatable, intent(out) :: stvec
    type(domain_t),              intent(in)  :: domain
    integer(kind=4),             intent(in)  :: halo_width_xy, halo_width_z

    type(stvec_nh_t), allocatable :: stvec_nh

    allocate(stvec_nh)

    call create_grid_field(stvec_nh%u, halo_width_xy, halo_width_z, domain%mesh_u)
    call create_grid_field(stvec_nh%v, halo_width_xy, halo_width_z, domain%mesh_v)
    call create_grid_field(stvec_nh%eta_dot, halo_width_xy, halo_width_z, domain%mesh_w)
    call create_grid_field(stvec_nh%theta, halo_width_xy, halo_width_z, domain%mesh_w)
    call create_grid_field(stvec_nh%P, halo_width_xy, halo_width_z, domain%mesh_P)

    stvec_nh%model_time = 0.0_8

    call move_alloc(stvec_nh, stvec)

end subroutine create_stvec_nh

end module stvec_nh_factory_mod

module domain_factory_mod

use domain_mod,                only : domain_t
use topology_factory_mod,      only : create_topology
use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
use metric_mod,                only : metric_t
use metric_factory_mod,        only : create_metric, create_metric_by_config
use config_domain_mod,         only : config_domain_t
use mpi

implicit none

generic :: create_domain => create_domain_by_config, create_domain_by_arguments

contains

subroutine create_domain_by_arguments(domain, topology_type, staggering_type, nh, nz, &
                                      parcomm)

    use parcomm_mod,       only: parcomm_global, parcomm_t
    use config_domain_mod, only: config_domain_t

    type(domain_t),   intent(out) :: domain
    character(len=*), intent(in)  :: topology_type
    character(len=*), intent(in)  :: staggering_type
    integer(kind=4),  intent(in)  :: nh, nz
    type(parcomm_t),  optional, intent(in) :: parcomm

    type(config_domain_t) :: config_domain

    config_domain%N  = nh
    config_domain%Nz = nz
    config_domain%staggering_type     = staggering_type
    config_domain%vertical_staggering = "None"
    config_domain%metric_type         = "ecs"
    config_domain%topology_type       = topology_type
    config_domain%h_top               = 1.0_8
    call config_domain%config_metric%set_defaults()

    call create_domain_by_config(domain,config_domain,parcomm)

end subroutine create_domain_by_arguments

subroutine create_domain_by_config(domain, config, parcomm)

    use mesh_factory_mod,    only : create_mesh
    use parcomm_factory_mod, only : create_parcomm
    use parcomm_mod,         only : parcomm_global, parcomm_t

    type(domain_t),             intent(out) :: domain
    type(config_domain_t),      intent(in)  :: config
    type(parcomm_t),  optional, intent(in)  :: parcomm

    !class(metric_t), allocatable  :: metric
    integer(kind=4) :: halo_width
    integer(kind=4) :: mpi_comm_local

    !have to be passed as an argument in future
    halo_width = 8

    domain%topology = create_topology(config%topology_type)

    call create_metric_by_config(domain%metric, domain%topology, &
                                 config%metric_type, config%config_metric)

    domain%horizontal_staggering = config%staggering_type

    mpi_comm_local = parcomm_global%comm_w
    if(present(parcomm)) then
        domain%parcomm = parcomm
    else
        !call create_parcomm(parcomm_glocal%comm_w, domain%parcomm)
        domain%parcomm = parcomm_global
    end if

    call domain%partition%init(config%N, config%Nz, max(1,domain%parcomm%np/6), &
                               domain%parcomm%myid, domain%parcomm%Np,          &
                               config%staggering_type, strategy = 'default')

    !WORKAROUND vertical staggering
    call create_mesh(domain%mesh_o,  domain%partition, domain%metric, halo_width, &
                                     config%h_top, 'o', config%vertical_staggering)
    call create_mesh(domain%mesh_x,  domain%partition, domain%metric, halo_width, &
                                     config%h_top, 'x', config%vertical_staggering)
    call create_mesh(domain%mesh_y,  domain%partition, domain%metric, halo_width, &
                                     config%h_top, 'y', config%vertical_staggering)
    call create_mesh(domain%mesh_xy, domain%partition, domain%metric, halo_width, &
                                     config%h_top, 'xy',config%vertical_staggering)

    if (config%vertical_staggering ==  "CharneyPhilips") then
        call create_mesh(domain%mesh_z,  domain%partition, domain%metric, halo_width, &
                                         config%h_top, 'z', config%vertical_staggering)
        call create_mesh(domain%mesh_xyz,domain%partition, domain%metric, halo_width, &
                                         config%h_top, 'xyz',config%vertical_staggering)
    end if

    select case(config%staggering_type)
    case ('A')
        domain%mesh_p = domain%mesh_o
        domain%mesh_u = domain%mesh_o
        domain%mesh_v = domain%mesh_o
        domain%mesh_q = domain%mesh_o
    case ('Ah') !all degrees of freedom at corner points
        domain%mesh_p = domain%mesh_xy
        domain%mesh_u = domain%mesh_xy
        domain%mesh_v = domain%mesh_xy
        domain%mesh_q = domain%mesh_xy
    case ('C')
        domain%mesh_p = domain%mesh_o
        domain%mesh_u = domain%mesh_x
        domain%mesh_v = domain%mesh_y
        domain%mesh_q = domain%mesh_xy
    case ('Ch')
        domain%mesh_p = domain%mesh_xy
        domain%mesh_u = domain%mesh_y
        domain%mesh_v = domain%mesh_x
        domain%mesh_q = domain%mesh_o
    case default
        call parcomm_global%abort("domain_factory_mod, unknown staggering type: "//config%staggering_type)
    end select
    select case(config%vertical_staggering)
    case("None")
        domain%mesh_w = domain%mesh_p
    case("CharneyPhilips")
            select case(config%staggering_type)
            case ('A','C')
                domain%mesh_w = domain%mesh_z
            case ('Ah','Ch')
                domain%mesh_w = domain%mesh_xyz
            case default
                call parcomm_global%abort("domain_factory_mod, unknown staggering type: "//config%staggering_type)
            end select
    case default
        call parcomm_global%abort("domain_factory_mod, unknown vertical staggering type: "//config%vertical_staggering)
    end select
end subroutine create_domain_by_config

end module domain_factory_mod

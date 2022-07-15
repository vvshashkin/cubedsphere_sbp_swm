module metric_factory_mod

use metric_mod,        only : metric_t
use topology_mod,      only : topology_t
use parcomm_mod,       only : parcomm_global
use config_metric_mod, only : config_metric_t

implicit none

generic :: create_metric => create_metric_by_param, create_metric_by_config

contains

subroutine create_metric_by_config(metric, topology, metric_type, config)
    use ecs_metric_mod,            only : ecs_metric_t
    use ecs_metric_factory_mod,    only : create_ecs_metric

    class(metric_t), allocatable, intent(out) :: metric
    class(topology_t),            intent(in)  :: topology
    character(len=*),             intent(in)  :: metric_type
    type(config_metric_t),        intent(in)  :: config

    select case(metric_type)
    case("ecs")
        call create_ecs_metric(metric, topology, &
         config%scale, config%omega, config%rotation_matrix, config%rotation_axis)
    case("shallow_atmosphere_metric")
        call create_shallow_atmosphere_metric(metric,topology,metric_type,config)
    case default
        call parcomm_global%abort("Unknown metric_type var " // metric_type // &
                                                    " in metric_factory_mod")
    end select
end subroutine create_metric_by_config

subroutine create_metric_by_param(metric, topology, metric_type)
    use ecs_metric_mod,            only : ecs_metric_t
    use ecs_metric_factory_mod,    only : create_ecs_metric

    class(topology_t),            intent(in)  :: topology
    character(len=*),             intent(in)  :: metric_type
    class(metric_t), allocatable, intent(out) :: metric

    type(config_metric_t) :: config

    call config%set_defaults()
    call create_metric_by_config(metric,topology,metric_type,config)

end subroutine create_metric_by_param

subroutine create_shallow_atmosphere_metric(metric, topology, metric_type, config)

    use shallow_atm_metric_mod, only : shallow_atm_metric_t
    use vertical_transform_factory_mod, only : create_vertical_transform

    class(metric_t), allocatable, intent(out) :: metric
    class(topology_t),            intent(in)  :: topology
    character(len=*),             intent(in)  :: metric_type
    type(config_metric_t),        intent(in)  :: config

    type(shallow_atm_metric_t), allocatable :: metric_shallow_atm

    allocate(metric_shallow_atm)

    call create_metric_by_config(metric_shallow_atm%metric_2d, topology, config%metric_2d_type, config)

    metric_shallow_atm%scale = config%scale
    metric_shallow_atm%vertical_scale = config%vertical_scale
    metric_shallow_atm%omega = config%omega
    metric_shallow_atm%rotation_axis = config%rotation_axis
    metric_shallow_atm%rotation_matrix = config%rotation_matrix
    metric_shallow_atm%alpha0 = metric_shallow_atm%metric_2d%alpha0
    metric_shallow_atm%alpha1 = metric_shallow_atm%metric_2d%alpha1
    metric_shallow_atm%beta0 = metric_shallow_atm%metric_2d%beta0
    metric_shallow_atm%beta1 = metric_shallow_atm%metric_2d%beta1

    metric_shallow_atm%vertical_transform = &
                    create_vertical_transform(config%vertical_transform_name)

    call move_alloc(metric_shallow_atm, metric)

end subroutine create_shallow_atmosphere_metric

end module metric_factory_mod

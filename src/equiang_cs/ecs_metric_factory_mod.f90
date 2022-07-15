module ecs_metric_factory_mod

use parcomm_mod,  only : parcomm_global

implicit none

contains

subroutine create_ecs_metric(metric, topology, sphere_r, sphere_omega, rotation_matrix, rotation_axis)

    use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
    use topology_mod,              only : topology_t
    use metric_mod,                only : metric_t
    use ecs_metric_mod,            only : ecs_metric_t
    use const_mod,                 only : pi

    class(topology_t),             intent(in)  :: topology
    class(metric_t), allocatable,  intent(out) :: metric
    real(kind=8), optional,        intent(in)  :: sphere_r, sphere_omega
    real(kind=8), optional,        intent(in)  :: rotation_matrix(3,3), rotation_axis(3)

    type(ecs_metric_t), allocatable :: ecs_metric
    real(kind=8) :: local_rotation_matrix(3, 3), local_rotation_axis(3)
    real(kind=8) :: local_sphere_r, local_sphere_omega

    local_rotation_matrix = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    if (present(rotation_matrix)) local_rotation_matrix = rotation_matrix

    local_rotation_axis = [0, 0, 1]
    if (present(rotation_axis)) local_rotation_axis = rotation_axis

    local_sphere_omega = 1.0_8
    if (present(sphere_omega)) local_sphere_omega = sphere_omega

    local_sphere_r = 1.0_8
    if (present(sphere_r)) local_sphere_r = sphere_r

    allocate(ecs_metric)

    ecs_metric%alpha0 =-0.25*pi
    ecs_metric%beta0  =-0.25*pi
    ecs_metric%alpha1 = 0.25*pi
    ecs_metric%beta1  = 0.25*pi

    ecs_metric%scale = local_sphere_r
    ecs_metric%vertical_scale = 1.0_8
    ecs_metric%omega = local_sphere_omega

    ecs_metric%rotation_matrix = local_rotation_matrix
    ecs_metric%rotation_axis   = local_rotation_axis

    select type(topology)
    class is (cubed_sphere_topology_t)
        ecs_metric%topology = topology
    class default
        call parcomm_global%abort("Wrong topology class in create_ecs_metric!")
    end select

    call move_alloc(ecs_metric, metric)

end

end module ecs_metric_factory_mod

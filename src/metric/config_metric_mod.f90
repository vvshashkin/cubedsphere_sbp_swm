module config_metric_mod

use config_mod, only : config_t

implicit none

type, public, extends(config_t) :: config_metric_t

    integer(kind=4) :: N, Nz

    real(kind=8) :: scale = 1.0_8
    real(kind=8) :: vertical_scale = 1.0_8
    real(kind=8) :: omega = 1.0_8
    real(kind=8) :: rotation_matrix(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    real(kind=8) :: rotation_axis(3)     = [0, 0, 1]
    character(:), allocatable :: metric_2d_type
    character(:), allocatable :: vertical_transform_name
contains
    procedure, public :: parse
    procedure, public :: set_defaults
end type config_metric_t

contains

subroutine set_defaults(this)
    class(config_metric_t), intent(inout) :: this

    this%scale = 1.0_8
    this%vertical_scale = 1.0_8
    this%omega = 1.0_8
    this%rotation_matrix(1:3,1:3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    this%rotation_axis(1:3)     = [0, 0, 1]
    this%metric_2d_type = "ecs"
    this%vertical_transform_name = "vertical_transform_default"
end subroutine set_defaults

subroutine parse(this, config_string)

    use parcomm_mod, only : parcomm_global

    class(config_metric_t), intent(inout) :: this
    character(len=*),       intent(in)    :: config_string

    real(kind=8) :: scale = 1.0_8
    real(kind=8) :: omega = 1.0_8
    real(kind=8) :: rotation_matrix(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    real(kind=8) :: rotation_axis(3)     = [0, 0, 1]
    character(len=255) :: metric_2d_type = "ecs"
    character(len=255) :: vertical_transform_name = "vertical_transform_default"

    namelist /metric/ scale, omega, rotation_matrix, rotation_axis, &
                      metric_2d_type, vertical_transform_name

    read(config_string, metric)

    this%scale = scale
    this%omega = omega
    this%rotation_matrix = rotation_matrix
    this%rotation_axis(1:3) = rotation_axis(1:3) / sqrt(sum(rotation_axis(1:3)**2))
    this%metric_2d_type = trim(metric_2d_type)
    this%vertical_transform_name = trim(vertical_transform_name)

end subroutine parse

end module config_metric_mod

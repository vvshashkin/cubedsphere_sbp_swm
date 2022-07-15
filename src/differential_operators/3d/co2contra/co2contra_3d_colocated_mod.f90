module co2contra_3d_colocated_mod

use abstract_co2contra_3d_mod, only : co2contra_3d_operator_t
use grid_field_mod,            only : grid_field_t
use domain_mod,                only : domain_t

implicit none

type, extends(co2contra_3d_operator_t), public :: co2contra_3d_colocated_t
contains
    procedure :: transform    => transform_co2contra_3d_colocated
    procedure :: transform2co => transform_contra2co_3d_colocated
end type co2contra_3d_colocated_t

contains

subroutine transform_co2contra_3d_colocated(this, u_contra, v_contra, w_contra, u_cov, v_cov, w_cov, domain)
    class(co2contra_3d_colocated_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u_cov, v_cov, w_cov
    !output:
    type(grid_field_t),           intent(inout) :: u_contra, v_contra, w_contra

    integer(kind=4) :: t

    do t=domain%mesh_p%ts, domain%mesh_p%te
        call transform_co2contra_3d_colocated_tile(                   &
             u_contra%tile(t), v_contra%tile(t), w_contra%tile(t), &
             u_cov%tile(t),    v_cov%tile(t),    w_cov%tile(t),    &
             domain%mesh_p%tile(t))
    end do

end subroutine transform_co2contra_3d_colocated

subroutine transform_co2contra_3d_colocated_tile(u_contra, v_contra, w_contra, u_cov, v_cov, w_cov, mesh)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(in)    :: u_cov, v_cov, w_cov
    type(tile_mesh_t),  intent(in)    :: mesh
    !output
    type(tile_field_t), intent(inout) :: u_contra, v_contra, w_contra

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                u_contra%p(i,j,k) = mesh%Qi(1,i,j,k)*u_cov%p(i,j,k) + &
                                    mesh%Qi(2,i,j,k)*v_cov%p(i,j,k) + &
                                    mesh%Qi(4,i,j,k)*w_cov%p(i,j,k)

                v_contra%p(i,j,k) = mesh%Qi(2,i,j,k)*u_cov%p(i,j,k) + &
                                    mesh%Qi(3,i,j,k)*v_cov%p(i,j,k) + &
                                    mesh%Qi(5,i,j,k)*w_cov%p(i,j,k)

                w_contra%p(i,j,k) = mesh%Qi(4,i,j,k)*u_cov%p(i,j,k) + &
                                    mesh%Qi(5,i,j,k)*v_cov%p(i,j,k) + &
                                    mesh%Qi(6,i,j,k)*w_cov%p(i,j,k)
            end do
        end do
    end do

end subroutine transform_co2contra_3d_colocated_tile

subroutine transform_contra2co_3d_colocated(this, u_cov, v_cov, w_cov, u_contra, v_contra, w_contra, domain)
    class(co2contra_3d_colocated_t), intent(inout) :: this
    type(domain_t),               intent(in)    :: domain
    type(grid_field_t),           intent(inout) :: u_contra, v_contra, w_contra
    !output:
    type(grid_field_t),           intent(inout) :: u_cov, v_cov, w_cov

    integer(kind=4) :: t

    do t=domain%mesh_p%ts, domain%mesh_p%te
        call transform_contra2co_3d_colocated_tile(                 &
              u_cov%tile(t),    v_cov%tile(t),    w_cov%tile(t),    &
              u_contra%tile(t), v_contra%tile(t), w_contra%tile(t), &
              domain%mesh_u%tile(t) )
    end do

end subroutine transform_contra2co_3d_colocated

subroutine transform_contra2co_3d_colocated_tile(u_cov, v_cov, w_cov, u_contra, v_contra, w_contra, mesh)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(in)    :: u_contra, v_contra, w_contra
    type(tile_mesh_t),  intent(in)    :: mesh
    !output
    type(tile_field_t), intent(inout) :: u_cov, v_cov, w_cov

    integer(kind=4) :: i, j, k

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                u_cov%p(i,j,k) = mesh%Q(1,i,j,k)*u_contra%p(i,j,k) + &
                                 mesh%Q(2,i,j,k)*v_contra%p(i,j,k) + &
                                 mesh%Q(4,i,j,k)*w_contra%p(i,j,k)

                v_cov%p(i,j,k) = mesh%Q(2,i,j,k)*u_contra%p(i,j,k) + &
                                 mesh%Q(3,i,j,k)*v_contra%p(i,j,k) + &
                                 mesh%Q(5,i,j,k)*w_contra%p(i,j,k)

                w_cov%p(i,j,k) = mesh%Q(4,i,j,k)*u_contra%p(i,j,k) + &
                                 mesh%Q(5,i,j,k)*v_contra%p(i,j,k) + &
                                 mesh%Q(6,i,j,k)*w_contra%p(i,j,k)
            end do
        end do
    end do

end subroutine transform_contra2co_3d_colocated_tile

end module co2contra_3d_colocated_mod

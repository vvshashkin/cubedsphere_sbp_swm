module test_ecs_halo_c_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

implicit none

private
public   :: test_ecs_cvec_halo

contains

subroutine test_ecs_cvec_halo(halo_procedure, components_type)

    use domain_factory_mod,         only : create_domain
    use grid_field_factory_mod,     only : create_grid_field
    use halo_mod,                   only : halo_t, halo_vec_t
    use halo_factory_mod,           only : create_halo_procedure, create_vector_halo_procedure
    use test_fields_mod,            only : set_vector_test_field, solid_rot => solid_rotation_field_generator, &
                                           cross_polar => cross_polar_flow_generator

    character(len=*), intent(in) :: halo_procedure, components_type

    integer(kind=4), parameter         :: nh=64, nz=3, halo_width=3, ex_halo_width=8

    type(domain_t)             :: domain
    type(grid_field_t)         :: u_test, v_test
    type(grid_field_t)         :: u_true, v_true
    class(halo_vec_t), allocatable :: halo_vec

    real(kind=8) inface_err, cross_edge_err
    real(kind=8) inface_err_max, cross_edge_err_max
    real(kind=8) inface_corner_err, inface_corner_err_max
    real(kind=8) inedge_corner_err, inedge_corner_err_max
    real(kind=8) halo_corner_err, halo_corner_err_max

    logical is_test_passed
    real(kind=8), parameter :: tolerance = 3.0e-7_8
    integer(kind=4) :: t, j, i

    call create_domain(domain, "cube", 'C', nh, nz)

    call create_vector_halo_procedure(halo_vec,domain,halo_width,halo_procedure)
    !call create_vector_halo_procedure(halo_vec,domain,halo_width,"C_vec_default")

    call create_grid_field(u_test, ex_halo_width, 0, domain%mesh_u)
    call create_grid_field(u_true, ex_halo_width, 0, domain%mesh_u)
    call create_grid_field(v_test, ex_halo_width, 0, domain%mesh_v)
    call create_grid_field(v_true, ex_halo_width, 0, domain%mesh_v)

    call set_vector_test_field(u_test,v_test,solid_rot, domain%mesh_u, domain%mesh_v, &
                               0, components_type, 0.0_8)
    call set_vector_test_field(u_true,v_true,solid_rot, domain%mesh_u, domain%mesh_v, &
                               halo_width, components_type, 0.0_8)

    call domain%parcomm%print('equiangular cubed-sphere C-grid halo-zone interpolation test')

    call halo_vec%get_halo_vector(u_test,v_test,domain,halo_width)

    call halo_err(inface_err, inface_err_max, cross_edge_err, cross_edge_err_max, &
                  inface_corner_err, inface_corner_err_max, inedge_corner_err,    &
                  inedge_corner_err_max, halo_corner_err, halo_corner_err_max,    &
                  domain%partition%nh+1, domain%partition%nh, halo_width,         &
                  domain%mesh_u, u_test, u_true)

    if(domain%parcomm%myid == 0) then
        print *, "halo errors: U"
 !       print *, "inface mean abs, inface max, cross edge mean abs, cross edge max"
 !       print '(4e15.7)', inface_err, inface_err_max, cross_edge_err, cross_edge_err_max
 !       print *, "inface-corner mean abs, inface-corner max, inedge-corner mean abs, inedge-corner max, halo-corner mean abs, halo-corner max"
 !       print '(6e15.7)', inface_corner_err, inface_corner_err_max, inedge_corner_err, inedge_corner_err_max,&
 !                       halo_corner_err, halo_corner_err_max

        print *, cross_edge_err, cross_edge_err_max
        is_test_passed = max(inedge_corner_err_max, halo_corner_err_max, &
                             cross_edge_err_max) < tolerance
    end if

    call halo_err(inface_err, inface_err_max, cross_edge_err, cross_edge_err_max, &
                  inface_corner_err, inface_corner_err_max, inedge_corner_err,    &
                  inedge_corner_err_max, halo_corner_err, halo_corner_err_max,    &
                  domain%partition%nh, domain%partition%nh+1, halo_width,         &
                  domain%mesh_v, v_test, v_true)

    if(domain%parcomm%myid == 0) then
        print *, "halo errors: V"
!        print *, "inface mean abs, inface max, cross edge mean abs, cross edge max"
!        print '(4e15.7)', inface_err, inface_err_max, cross_edge_err, cross_edge_err_max
!        print *, "inface-corner mean abs, inface-corner max, inedge-corner mean abs, inedge-corner max, halo-corner mean abs, halo-corner max"
!        print '(6e15.7)', inface_corner_err, inface_corner_err_max, inedge_corner_err, inedge_corner_err_max,&
!                        halo_corner_err, halo_corner_err_max

        print *, cross_edge_err, cross_edge_err_max
        is_test_passed = max(inedge_corner_err_max, halo_corner_err_max, &
                             cross_edge_err_max) < tolerance
    end if

    if(domain%parcomm%myid == 0) then
        if (is_test_passed) then
            print *, "halo test passed"
        else
            print *, "halo test failed"
        end if
    end if

    !do t=1,6
    !    print *, "face ", t
    !    do j=1, 4
    !        print '(i2,7f15.7)', j, u_test%tile(t)%p(0:6,j,2)
    !        print '(i2,7f15.7)', j, u_true%tile(t)%p(0:6,j,2)
    !    end do
    !end do
    ! do i=1, nh
    !     print '(i3,3F15.7)', i, u_test%tile(2)%p(nh+2,i,1)-u_true%tile(2)%p(nh+2,i,1),&
    !                             u_test%tile(2)%p(nh+3,i,1)-u_true%tile(2)%p(nh+3,i,1)
    ! end do
    !print *, "Err1", maxval(abs(v_test%tile(1)%p(1:nh,-2,1)-v_true%tile(1)%p(1:nh,-2,1)))
    !print *, "Err1", maxval(abs(v_test%tile(1)%p(1:nh,-1,1)-v_true%tile(1)%p(1:nh,-1,1)))
    !print *, "Err1", maxval(abs(v_test%tile(1)%p(1:nh,0,1)-v_true%tile(1)%p(1:nh,0,1)))
    !print *, "Err2", maxval(abs(v_test%tile(1)%p(1:nh,nh+2,1)-v_true%tile(1)%p(1:nh,nh+2,1)))
    !print *, "Err3", maxval(abs(v_test%tile(1)%p(1:nh,nh+3,1)-v_true%tile(1)%p(1:nh,nh+3,1)))
    !print *, "Err4", maxval(abs(v_test%tile(1)%p(1:nh,nh+4,1)-v_true%tile(1)%p(1:nh,nh+4,1)))
    ! print *, "V1", v_test%tile(1)%p(nh+2,17,1), v_true%tile(1)%p(nh+2,17,1)
    ! print *, "Err1", maxval(abs(u_test%tile(2)%p(1:nh+1,-2,1)-u_true%tile(2)%p(1:nh+1,-2,1)))
    ! print *, "Err1", maxval(abs(u_test%tile(2)%p(1:nh+1,-1,1)-u_true%tile(2)%p(1:nh+1,-1,1)))
    ! print *, "Err1", maxval(abs(u_test%tile(2)%p(1:nh+1,0,1)-u_true%tile(2)%p(1:nh+1,0,1)))
    ! print *, "Err2", maxval(abs(u_test%tile(2)%p(1:nh+1,nh+1,1)-u_true%tile(2)%p(1:nh+1,nh+1,1)))
    ! print *, "Err3", maxval(abs(u_test%tile(2)%p(1:nh+1,nh+2,1)-u_true%tile(2)%p(1:nh+1,nh+2,1)))
    ! print *, "Err4", maxval(abs(u_test%tile(2)%p(1:nh+1,nh+3,1)-u_true%tile(2)%p(1:nh+1,nh+3,1)))
    !do j=domain%mesh_u%tile(3)%is,domain%mesh_u%tile(3)%ie
    !    !print '(i3,3E15.7)', j, v_test%tile(1)%p(nh+2,j,1),v_true%tile(1)%p(nh+2,j,1),&
    !    !                        v_test%tile(1)%p(nh+2,j,1)-v_true%tile(1)%p(nh+2,j,1)
    !    print '(i3,3E15.7)', j, u_test%tile(3)%p(j,nh+1,1),u_true%tile(3)%p(j,nh+1,1),&
    !                            u_test%tile(3)%p(j,nh+1,1)-u_true%tile(3)%p(j,nh+1,1)
    !end do
end subroutine test_ecs_cvec_halo

subroutine halo_err(gl_inface_err, gl_inface_err_max, gl_cross_edge_err, gl_cross_edge_err_max, &
                    gl_inface_corner_err, gl_inface_corner_err_max, gl_inedge_corner_err,       &
                    gl_inedge_corner_err_max, gl_halo_corner_err, gl_halo_corner_err_max,       &
                    nx, ny, halo_width, mesh, f1, f2)

use mpi
use mesh_mod, only : mesh_t

real(kind=8),          intent(out) :: gl_inface_err, gl_cross_edge_err
real(kind=8),          intent(out) :: gl_inface_err_max, gl_cross_edge_err_max
real(kind=8),          intent(out) :: gl_inface_corner_err, gl_inface_corner_err_max
real(kind=8),          intent(out) :: gl_inedge_corner_err, gl_inedge_corner_err_max
real(kind=8),          intent(out) :: gl_halo_corner_err, gl_halo_corner_err_max
integer(kind=4),       intent(in)  :: nx, ny, halo_width
type(mesh_t),          intent(in)  :: mesh
type(grid_field_t),    intent(in)  :: f1, f2
!locals
real(kind=8)      err, inface_err, cross_edge_err
real(kind=8)      err_max, inface_err_max, cross_edge_err_max
real(kind=8)      inface_corner_err, inface_corner_err_max
real(kind=8)      inedge_corner_err, inedge_corner_err_max
real(kind=8)      halo_corner_err, halo_corner_err_max
integer(kind=4)   num_inedge_corners, gl_num_inedge_corners
integer(kind=4)   is, ie, js, je, ks, ke, klev, ind, ierr

inface_err = 0._8;          inface_err_max = 0._8
cross_edge_err = 0._8;      cross_edge_err_max = 0._8
inface_corner_err = 0._8;   inface_corner_err_max = 0._8
inedge_corner_err = 0._8;   inedge_corner_err_max = 0._8
halo_corner_err = 0._8;     halo_corner_err_max = 0._8
num_inedge_corners = 0

do ind = mesh%ts, mesh%te

    is = mesh%tile(ind)%is; ie = mesh%tile(ind)%ie
    js = mesh%tile(ind)%js; je = mesh%tile(ind)%je
    ks = mesh%tile(ind)%ks; ke = mesh%tile(ind)%ke
    klev = ke-ks+1
    !tile edge errors
    err     = sum(abs(f1%tile(ind)%p(is-halo_width:is-1,js:je,ks:ke)-f2%tile(ind)%p(is-halo_width:is-1,js:je,ks:ke)))/ny
    err_max = maxval(abs(f1%tile(ind)%p(is-halo_width:is-1,js:je,ks:ke)-f2%tile(ind)%p(is-halo_width:is-1,js:je,ks:ke)))
    if(is == 1) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     =    sum(abs(f1%tile(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)-f2%tile(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)))/ny
    err_max = maxval(abs(f1%tile(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)-f2%tile(ind)%p(ie+1:ie+halo_width,js:je,ks:ke)))
    if(ie == nx) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     =    sum(abs(f1%tile(ind)%p(is:ie,js-halo_width:js-1,ks:ke)-f2%tile(ind)%p(is:ie,js-halo_width:js-1,ks:ke)))/nx
    err_max = maxval(abs(f1%tile(ind)%p(is:ie,js-halo_width:js-1,ks:ke)-f2%tile(ind)%p(is:ie,js-halo_width:js-1,ks:ke)))
    if(js == 1) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err
        inface_err_max = max(inface_err_max,err_max)
    end if

    err     =    sum(abs(f1%tile(ind)%p(is:ie,je+1:je+halo_width,ks:ke)-f2%tile(ind)%p(is:ie,je+1:je+halo_width,ks:ke)))/nx
    err_max = maxval(abs(f1%tile(ind)%p(is:ie,je+1:je+halo_width,ks:ke)-f2%tile(ind)%p(is:ie,je+1:je+halo_width,ks:ke)))
    if(je == ny) then
        cross_edge_err     = cross_edge_err+err
        cross_edge_err_max = max(cross_edge_err_max,err_max)
    else
        inface_err     = inface_err+err
        inface_err_max = max(inface_err_max,err_max)
    end if

 !   err     = sum(abs(f1%tile(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke)-f2%tile(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke))) / halo_width**2
 !   err_max = maxval(abs(f1%tile(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke)-f2%tile(ind)%p(is-halo_width:is-1,js-halo_width:js-1,ks:ke)))
 !   if(is == 1 .and. js == 1) then
 !       halo_corner_err = halo_corner_err+err
 !       halo_corner_err_max = max(halo_corner_err_max,err_max)
 !   else if(is == 1 .or. js == 1) then
 !       inedge_corner_err = inedge_corner_err+err
 !       inedge_corner_err_max = max(inedge_corner_err_max,err_max)
 !       num_inedge_corners = num_inedge_corners+1
 !   else
 !       inface_corner_err = inface_corner_err+err
 !       inface_corner_err_max = max(inface_corner_err_max,err_max)
 !   end if

 !   err     = sum(abs(f1%tile(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke)-f2%tile(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke))) / halo_width**2
 !   err_max = maxval(abs(f1%tile(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke)-f2%tile(ind)%p(ie+1:ie+halo_width,js-halo_width:js-1,ks:ke)))
 !   if(ie == nx .and. js == 1) then
 !       halo_corner_err = halo_corner_err+err
 !       halo_corner_err_max = max(halo_corner_err_max,err_max)
 !   else if(ie == nx .or. js == 1) then
 !       inedge_corner_err = inedge_corner_err+err
 !       inedge_corner_err_max = max(inedge_corner_err_max,err_max)
 !       num_inedge_corners = num_inedge_corners+1
 !   else
 !       inface_corner_err = inface_corner_err+err
 !       inface_corner_err_max = max(inface_corner_err_max,err_max)
 !   end if

!    err     = sum(abs(f1%tile(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke)-f2%tile(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke))) / halo_width**2
!    err_max = maxval(abs(f1%tile(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke)-f2%tile(ind)%p(ie+1:ie+halo_width,je+1:je+halo_width,ks:ke)))
!    if(ie == nx .and. je == ny) then
!        halo_corner_err = halo_corner_err+err
!        halo_corner_err_max = max(halo_corner_err_max,err_max)
!    else if(ie == nx .or. je == ny) then
!        inedge_corner_err = inedge_corner_err+err
!        inedge_corner_err_max = max(inedge_corner_err_max,err_max)
!        num_inedge_corners = num_inedge_corners+1
!    else
!        inface_corner_err = inface_corner_err+err
!    end if
!        inface_corner_err_max = max(inface_corner_err_max,err_max)!
!
!    err     = sum(abs(f1%tile(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke)-f2%tile(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke))) / halo_width**2
!    err_max = maxval(abs(f1%tile(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke)-f2%tile(ind)%p(is-halo_width:is-1,je+1:je+halo_width,ks:ke)))
!    if(is == 1 .and. je == ny) then
!        halo_corner_err = halo_corner_err+err
!        halo_corner_err_max = max(halo_corner_err_max,err_max)
!    else if(is == 1 .or. je == ny) then
!        inedge_corner_err = inedge_corner_err+err
!        inedge_corner_err_max = max(inedge_corner_err_max,err_max)
!        num_inedge_corners = num_inedge_corners+1
!    else
!        inface_corner_err = inface_corner_err+err
!        inface_corner_err_max = max(inface_corner_err_max,err_max)
!    end if
    !print *, "error", ind, cross_edge_err_max
end do

call mpi_allreduce(cross_edge_err, gl_cross_edge_err, 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(inface_err,     gl_inface_err    , 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(cross_edge_err_max, gl_cross_edge_err_max, 1, mpi_double, mpi_max, mpi_comm_world, ierr)
call mpi_allreduce(inface_err_max,     gl_inface_err_max    , 1, mpi_double, mpi_max, mpi_comm_world, ierr)

gl_cross_edge_err = gl_cross_edge_err/(6*4*klev*halo_width) !scale to 1 element error

call mpi_allreduce(inface_corner_err,     gl_inface_corner_err    , 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(inface_corner_err_max,     gl_inface_corner_err_max    , 1, mpi_double, mpi_max, mpi_comm_world, ierr)
call mpi_allreduce(inedge_corner_err,     gl_inedge_corner_err    , 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(inedge_corner_err_max,     gl_inedge_corner_err_max    , 1, mpi_double, mpi_max, mpi_comm_world, ierr)
call mpi_allreduce(num_inedge_corners,     gl_num_inedge_corners    , 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(halo_corner_err,     gl_halo_corner_err    , 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(halo_corner_err_max,     gl_halo_corner_err_max    , 1, mpi_double, mpi_max, mpi_comm_world, ierr)

gl_inedge_corner_err = gl_inedge_corner_err / max(gl_num_inedge_corners,1) / klev
gl_halo_corner_err = gl_halo_corner_err / 48 / klev

end subroutine halo_err

end module test_ecs_halo_c_mod

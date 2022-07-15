module test_latlon_output_mod

use domain_mod,             only : domain_t
use domain_factory_mod,     only : create_domain
use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field
use outputer_abstract_mod,  only : outputer_t, outputer_vector_t
use outputer_factory_mod,   only : create_latlon_outputer,&
                                   create_latlon_vec_outputer
use parcomm_mod,            only : parcomm_global

implicit none

contains

subroutine test_latlon_output(staggering, scalar_grid,irec)

    use test_fields_mod,        only : set_scalar_test_field, &
                                       xyz_f=>xyz_scalar_field_generator

    character(len=*), intent(in) :: staggering, scalar_grid
    integer(kind=4),  intent(in) :: irec

    class(outputer_t), allocatable  :: outputer
    type(domain_t)                  :: domain
    type(grid_field_t)              :: f

    character(*), parameter    :: file_name = "h.dat"
    integer(kind=4), parameter :: nh = 32, nz = 3
    integer(kind=4), parameter :: nlon=120, nlat=61

    call create_domain(domain, "cube", staggering, nh, nz)

    if(scalar_grid == "A") then
        call create_grid_field(f, 0, 0, domain%mesh_o)
        call set_scalar_test_field(f,xyz_f,domain%mesh_o,0)
    else if(scalar_grid == "Ah") then
        call create_grid_field(f, 0, 0, domain%mesh_xy)
        call set_scalar_test_field(f,xyz_f,domain%mesh_xy,0)
    else
        call parcomm_global%abort("test_latlon_outputter_mod, "//&
                                  "unsupported type of scalar point:"//&
                                  scalar_grid)
    end if

    call create_latlon_outputer(outputer, nlat, nlon, scalar_grid, domain)
    call outputer%write(f,domain,"h.dat",irec)

    if(domain%parcomm%myid == 0) then
        print *, "Test latlon outputer, staggering: "//staggering //&
                " scalar  grid: "// scalar_grid//" passed"
    end if

end subroutine test_latlon_output

subroutine test_latlon_vec_output(staggering,components_type,irec)

    use test_fields_mod,      only : set_vector_test_field, &
                                     vec_field_gen=>solid_rotation_field_generator

    character(len=*), intent(in) :: staggering, components_type
    integer(kind=4),  intent(in) :: irec

    class(outputer_vector_t), allocatable  :: outputer
    type(domain_t)                         :: domain
    type(grid_field_t)                     :: u, v

    character(*), parameter    :: file_name = "uv.dat"
    integer(kind=4), parameter :: nh = 32, nz = 3
    integer(kind=4), parameter :: nlon=120, nlat=61

    call create_domain(domain, "cube", staggering, nh, nz)

    call create_grid_field(u, 0, 0, domain%mesh_u)
    call create_grid_field(v, 0, 0, domain%mesh_v)
    call set_vector_test_field(u,v,vec_field_gen,domain%mesh_u,domain%mesh_v,0,components_type)

    call create_latlon_vec_outputer(outputer, nlat, nlon, staggering, &
                                   components_type, domain)
    call outputer%write(u,v,domain,"u.dat","v.dat",irec)

    if(domain%parcomm%myid == 0) then
        print *, "Test vector latlon outputer, staggering: "//staggering //&
                ", components type: "// components_type//" passed"
    end if

end subroutine test_latlon_vec_output

end module test_latlon_output_mod

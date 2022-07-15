module test_regrid_mod

use abstract_regrid_mod,    only : regrid_t, regrid_vec_t
use regrid_factory_mod,     only : create_latlon_regrid, create_latlon_vector_regrid
use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field
use domain_mod,             only : domain_t
use domain_factory_mod,     only : create_domain

use const_mod,              only : pi

implicit none

contains

subroutine test_regrid(staggering)
    use test_fields_mod,        only : set_scalar_test_field, &
                                       xyz_f=>xyz_scalar_field_generator

    character(len=*), intent(in)    :: staggering

    class(regrid_t),    allocatable :: regrid
    type(grid_field_t)              :: f
    type(domain_t)                  :: domain
    integer(kind=4),    parameter   :: N=32, nz=3
    integer(kind=4),    parameter   :: Nlon = 4*N, Nlat=2*N+1

    real(kind=8),       parameter   :: hlon = 2._8*pi / real(Nlon,8)
    real(kind=8),       parameter   :: hlat = pi / (Nlat-1)

    real(kind=8), dimension(Nlon,Nlat)    :: x_latlon, y_latlon, z_latlon
    real(kind=8), dimension(Nlon,Nlat,nz) :: ptest, ptest_cub, ptrue
    real(kind=8) :: lon, lat
    integer(kind=4) :: i,j

    print *, "scalar latlon regrid test, staggering: "//staggering

    call create_domain(domain, "cube", staggering, N, nz)
    call create_grid_field(f,0,0,domain%mesh_p)
    call set_scalar_test_field(f,xyz_f,domain%mesh_p,0)

    do j=1, Nlat
        lat = -0.5_8*pi+(j-1)*hlat
        do i=1, Nlon
            lon = (i-1)*hlon
            x_latlon(i,j) = cos(lat)*cos(lon)
            y_latlon(i,j) = cos(lat)*sin(lon)
            z_latlon(i,j) = sin(lat)
        end do
    end do
    call xyz_f%get_scalar_field(ptrue,Nlat*Nlon,nz,x_latlon,y_latlon,z_latlon)

    call create_latlon_regrid(regrid,domain,Nlat=Nlat,Nlon=Nlon,interp_type="linear",&
                              scalar_grid_type=staggering)
    call regrid%do_regrid(ptest,f,domain)

    call create_latlon_regrid(regrid,domain,Nlat=Nlat,Nlon=Nlon,interp_type="cubic",&
                              scalar_grid_type=staggering)
    call regrid%do_regrid(ptest_cub,f,domain)

    print *, "Interpolation error (linear)", maxval(abs(ptest-ptrue))
    print *, "Interpolation error (cubic) ", maxval(abs(ptest_cub-ptrue))

end subroutine test_regrid

subroutine test_regrid_vec(staggering,components_type)
    use test_fields_mod,      only : set_vector_test_field, &
                                     vec_field_gen=>solid_rotation_field_generator
    use latlon_functions_mod, only : get_latlon_uv=>transform_cartesian_hor_vector_to_latlon

    character(len=*), intent(in)    :: staggering, components_type

    class(regrid_vec_t), allocatable :: regrid
    type(grid_field_t)               :: u, v
    type(domain_t)                   :: domain
    integer(kind=4),     parameter   :: N=32, nz=3
    integer(kind=4),     parameter   :: Nlon = 4*N, Nlat=2*N+1

    real(kind=8),        parameter   :: hlon = 2._8*pi / real(Nlon,8)
    real(kind=8),        parameter   :: hlat = pi / (Nlat-1)

    real(kind=8), dimension(Nlon,Nlat)    :: lat, lon, x_latlon, y_latlon, z_latlon
    real(kind=8), dimension(Nlon,Nlat,nz) :: utest, utest_cub, utrue
    real(kind=8), dimension(Nlon,Nlat,nz) :: vtest, vtest_cub, vtrue
    real(kind=8), dimension(Nlon,Nlat,nz) :: vx, vy, vz
    integer(kind=4) :: i,j

    print *, "vector latlon regrid test, staggering: "//staggering//&
             ", components type: "//components_type

    call create_domain(domain, "cube", staggering, N, nz)
    call create_grid_field(u,0,0,domain%mesh_u)
    call create_grid_field(v,0,0,domain%mesh_v)
    call set_vector_test_field(u,v,vec_field_gen,domain%mesh_u,domain%mesh_v,0,components_type)

    do j=1, Nlat
        do i=1, Nlon
            lat(i,j) = -0.5_8*pi+(j-1)*hlat
            lon(i,j) = (i-1)*hlon
            x_latlon(i,j) = cos(lat(i,j))*cos(lon(i,j))
            y_latlon(i,j) = cos(lat(i,j))*sin(lon(i,j))
            z_latlon(i,j) = sin(lat(i,j))
        end do
    end do
    call vec_field_gen%get_vector_field(vx,vy,vz,Nlat*Nlon,nz,x_latlon,y_latlon,z_latlon)
    call get_latlon_uv(utrue,vtrue,Nlat*Nlon,nz,vx,vy,vz,lat,lon)

    call create_latlon_vector_regrid(regrid,domain,Nlat=Nlat,Nlon=Nlon,interp_type="linear",&
                                     vector_grid_type=staggering,components_type=components_type)
    call regrid%do_regrid(utest,vtest,u,v,domain)

    call create_latlon_vector_regrid(regrid,domain,Nlat=Nlat,Nlon=Nlon,interp_type="cubic",&
                                     vector_grid_type=staggering,components_type=components_type)
    call regrid%do_regrid(utest_cub,vtest_cub,u,v,domain)

    print *, "Interpolation error (linear)", maxval(abs(utest-utrue)),maxval(abs(vtest-vtrue))
    print *, "Interpolation error (cubic) ", maxval(abs(utest_cub-utrue)),maxval(abs(vtest_cub-vtrue))

end subroutine test_regrid_vec



end module test_regrid_mod

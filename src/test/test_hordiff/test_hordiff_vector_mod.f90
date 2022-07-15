module test_hordiff_vector_mod

use domain_mod,                     only : domain_t
use domain_factory_mod,             only : create_domain

use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field

use abstract_hordiff_mod, only : hordiff_operator_t
use hordiff_factory_mod,  only : create_hordiff_operator

use test_fields_mod, only : random_vector_field_generator_t, set_vector_test_field, solid_rotation_t

use outputer_abstract_mod, only : outputer_t, outputer_vector_t
use outputer_factory_mod,  only : create_master_paneled_outputer,&
                                           create_latlon_outputer, create_latlon_vec_outputer
use halo_mod,               only : halo_vec_t
use halo_factory_mod,       only : create_vector_halo_procedure
implicit none

contains

subroutine hordiff_vector_test()

    class(hordiff_operator_t), allocatable :: diff_op
    class(outputer_vector_t),  allocatable :: outputer_vec
    class(outputer_t),         allocatable :: outputer, outputer_u, outputer_v
    class(halo_vec_t),         allocatable :: Ah_sync
    type(domain_t) :: domain
    type(grid_field_t) :: u, v, u_tend, v_tend

    type(random_vector_field_generator_t) :: field

    integer(kind=4), parameter :: N = 32, nz = 1, halo_width = 10
    real(kind=8),    parameter :: diff_coeff = 0.3_8

    character(len=2), parameter :: staggering = "Ah"

    integer(kind=4) :: it


    call create_domain(domain, "cube", trim(staggering), N, nz)

    call create_grid_field(u,      5, 0, domain%mesh_u)
    call create_grid_field(v,      5, 0, domain%mesh_v)
    call create_grid_field(u_tend, 5, 0, domain%mesh_u)
    call create_grid_field(v_tend, 5, 0, domain%mesh_v)

    call create_hordiff_operator(diff_op, "hordiff_vec_xyz_Ah_sbp_42_narrow", diff_coeff, domain)

    call set_vector_test_field(u, v, field, domain%mesh_u, domain%mesh_v, 0, "covariant")

    if(trim(staggering)=="Ah") then
        call create_vector_halo_procedure(Ah_sync, domain, 1, "ecs_Ah_vec_sync_covariant")
        call Ah_sync%get_halo_vector(u, v, domain, 1)
    end if

    if(trim(staggering) == "Ah") then
        call create_latlon_outputer(outputer,          2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "Ah", &
                                   "covariant", domain)
    else if(trim(staggering) == "C") then
        call create_latlon_outputer(outputer,          2*domain%partition%Nh+1, 4*domain%partition%Nh, "A", domain)
        call create_latlon_vec_outputer(outputer_vec,  2*domain%partition%Nh+1, 4*domain%partition%Nh, "C", &
                                       "covariant", domain)
    end if

    call create_master_paneled_outputer(outputer_u, "u", domain)
    call create_master_paneled_outputer(outputer_v, "v", domain)

    do it = 1,400

        if (domain%parcomm%myid==0) print*, "it = ", it

        ! call outputer%write(h, domain, 'h_diff.dat', it)
        call outputer_vec%write(u, v, domain, 'u_diff.dat', 'v_diff.dat', it)
        call outputer_u%write(u, domain, 'u_diff_mp.dat', it)
        call outputer_u%write(v, domain, 'v_diff_mp.dat', it)

        call diff_op%calc_diff_vec(u_tend, v_tend, u, v, domain)
        call u%update(diff_coeff, u_tend, domain%mesh_u)
        call v%update(diff_coeff, v_tend, domain%mesh_v)
    end do

end subroutine hordiff_vector_test

end module test_hordiff_vector_mod

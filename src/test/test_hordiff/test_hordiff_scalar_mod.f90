module test_hordiff_scalar_mod

use domain_mod,                     only : domain_t
use domain_factory_mod,             only : create_domain

use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field

use abstract_hordiff_mod, only : hordiff_operator_t
use hordiff_factory_mod,  only : create_hordiff_operator

use test_fields_mod, only : random_scalar_field_generator_t, set_scalar_test_field

use outputer_abstract_mod, only : outputer_t
use outputer_factory_mod,  only : create_master_paneled_outputer,&
                                           create_latlon_outputer
use halo_mod,               only : halo_t
use halo_factory_mod,       only : create_halo_procedure

use abstract_quadrature_mod, only : quadrature_t
use quadrature_factory_mod,  only : create_quadrature

implicit none

contains

subroutine hordiff_scalar_test()

    class(hordiff_operator_t), allocatable :: diff_op
    class(outputer_t),   allocatable :: outputer, outputer_mp
    class(halo_t),       allocatable :: Ah_sync
    class(quadrature_t), allocatable :: quadrature
    type(domain_t) :: domain
    type(grid_field_t) :: h, h_tend

    type(random_scalar_field_generator_t) :: field

    integer(kind=4), parameter :: N = 32, nz = 1, halo_width = 10
    real(kind=8),    parameter :: diff_coeff = 0.25_8

    character(len=2), parameter :: staggering = "Ah"

    integer(kind=4) :: it
    real(kind=8) :: mass, l2_norm


    call create_domain(domain, "cube", trim(staggering), N, nz)

    call create_grid_field(h,      halo_width, 0, domain%mesh_p)
    call create_grid_field(h_tend, halo_width, 0, domain%mesh_p)

    call create_hordiff_operator(diff_op, "hordiff_scalar_Ah_sbp_63_narrow", diff_coeff, domain)

    call set_scalar_test_field(h, field, domain%mesh_p, 0)

    if(trim(staggering)=="Ah") then
        call create_halo_procedure(Ah_sync, domain, 1, "Ah_scalar_sync")
        call Ah_sync%get_halo_scalar(h, domain, 1)
    end if

    call create_latlon_outputer(outputer, 2*N+1, 4*N, "A", domain)

    call create_master_paneled_outputer(outputer_mp, "p", domain)

    call create_quadrature(quadrature, "SBP_Ah63_quadrature", domain%mesh_p)

    do it = 1, 600

        if (domain%parcomm%myid==0) print*, it

        call outputer%write(h, domain, 'h_diff.dat', it)

        call outputer_mp%write(h, domain, 'h_diff_mp.dat', it)

        call diff_op%calc_diff(h_tend, h, domain%mesh_p, domain)
        call h%update(diff_coeff, h_tend, domain%mesh_p)

        mass = quadrature%mass(h, domain%mesh_p, domain%parcomm)
        l2_norm = sqrt(quadrature%dot(h, h, domain%mesh_p, domain%parcomm))

        if (domain%parcomm%myid==0) print*, "mass", mass
        if (domain%parcomm%myid==0) print*, "l2 norm", l2_norm

    end do

end subroutine hordiff_scalar_test

end module test_hordiff_scalar_mod

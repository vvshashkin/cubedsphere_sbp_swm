module nh_model_mod

use domain_mod,     only : domain_t
use timescheme_mod, only : timescheme_t
use stvec_mod,      only : stvec_t
use operator_mod,   only : operator_t

use stvec_nh_mod, only : stvec_nh_t

use abstract_postprocessing_mod, only : postprocessing_t
use mpi

implicit none

type nh_model_t

    type(domain_t) :: domain

    class(timescheme_t),     allocatable :: timescheme
    class(stvec_t),          allocatable :: stvec
    class(operator_t),       allocatable :: operator
    class(postprocessing_t), allocatable :: postproc

    real(kind=8) :: dt
    real(kind=8) :: tau_write
    real(kind=8) :: tau_diagnostics
    real(kind=8) :: simulation_time

    contains
    procedure, public :: run

end type nh_model_t

contains

subroutine run(this)
    class(nh_model_t), intent(inout) :: this

    integer(kind=4) :: nstep, nzap
    integer(kind=4) :: istep
    real(kind=8)    :: tmax, pmax, wmax
    real(kind=8)    :: time

    nstep = int(this%simulation_time / this%dt)
    nzap  = int(this%tau_write / this%dt)

    call this%postproc%write_fields(1,this%stvec,this%domain)

    do istep = 1, nstep
        time = mpi_wtime()
        call this%timescheme%step(this%stvec, this%operator, this%domain, this%dt)
        if(mod(istep,nzap) == 0) then
            call this%postproc%write_fields(istep/nzap+1,this%stvec,this%domain)
            if(this%domain%parcomm%myid == 0) &
            print *, "wr", istep/nzap+1
        end if
        select type(stvec=>this%stvec)
        type is (stvec_nh_t)
            tmax = stvec%theta%maxabs(this%domain%mesh_w,this%domain%parcomm)
            wmax = stvec%eta_dot%maxabs(this%domain%mesh_w,this%domain%parcomm)
            pmax = stvec%P%maxabs(this%domain%mesh_p,this%domain%parcomm)
            if(this%domain%parcomm%myid == 0) print *, istep, "w_maxabs", wmax, "t_maxabs", tmax, "pmaxabs", pmax, "time=", mpi_wtime()-time
            !print *, istep, "w maxabs", stvec%eta_dot%maxabs(this%domain%mesh_w,this%domain%parcomm) &
            !              , "T maxabs", stvec%theta%maxabs(this%domain%mesh_w,this%domain%parcomm)   &
            !              , "P maxabs", stvec%P%maxabs(this%domain%mesh_p,this%domain%parcomm)
        class default
            print *, istep, "***"
        end select
    end do

end subroutine run

end module nh_model_mod

module timescheme_factory_mod

use stvec_mod,         only : stvec_t
use timescheme_mod,    only : timescheme_t
use parcomm_mod,       only : parcomm_global
use explicit_Eul1_mod, only : explicit_Eul1_t
use rk4_mod,           only : rk4_t

implicit none

contains

subroutine create_timescheme(timescheme, v, timescheme_name)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),   intent(in) :: v !example of model state-vector
    character(len=*), intent(in) :: timescheme_name

    select case(timescheme_name)
    case("explicit_Eul1")
        call create_explicit_Eul1_timescheme(timescheme, v)
    case("rk4")
        call create_rk4_timescheme(timescheme, v)
    case default
        call parcomm_global%abort("unknown timescheme_name in create_timescheme: "// &
                                  timescheme_name)
    end select

end subroutine create_timescheme

subroutine create_explicit_Eul1_timescheme(timescheme, v)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in)  :: v! model state-vector example

    type(explicit_Eul1_t), allocatable :: Eul1_timescheme

    allocate(Eul1_timescheme)
    call v%create_similar(Eul1_timescheme%tendency)
    call move_alloc(Eul1_timescheme, timescheme)

end subroutine create_explicit_Eul1_timescheme

subroutine create_rk4_timescheme(timescheme, v)

    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in)  :: v !example of model state vector

    type(rk4_t), allocatable :: rk4

    allocate(rk4)

    !preallocate additional state vectors
    call v%create_similar(rk4%k1)
    call v%create_similar(rk4%k2)
    call v%create_similar(rk4%k3)
    call v%create_similar(rk4%k4)
    call v%create_similar(rk4%y)

    call move_alloc(rk4, timescheme)
end subroutine create_rk4_timescheme

end module timescheme_factory_mod

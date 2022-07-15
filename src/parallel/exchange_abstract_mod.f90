module exchange_abstract_mod

use grid_field_mod, only : grid_field_t
use parcomm_mod,    only : parcomm_t

implicit none

type, abstract, public :: exchange_t
contains
    procedure(exchange_procedure),     deferred :: do
    procedure(exchange_vec_procedure), deferred :: do_vec
end type exchange_t

abstract interface
    subroutine exchange_procedure(this, f, parcomm)
        import exchange_t, grid_field_t, parcomm_t
        class(exchange_t),  intent(inout) :: this
        type(parcomm_t),    intent(in)    :: parcomm
        type(grid_field_t), intent(inout) :: f
    end subroutine
    subroutine exchange_vec_procedure(this, u, v, parcomm)
        import exchange_t, grid_field_t, parcomm_t
        class(exchange_t),  intent(inout) :: this
        type(parcomm_t),    intent(in)    :: parcomm
        type(grid_field_t), intent(inout) :: u, v
    end subroutine
end interface

contains

end module exchange_abstract_mod

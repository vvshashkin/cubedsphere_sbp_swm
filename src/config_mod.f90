module config_mod

implicit none

type, public :: config_t
contains
    procedure, public :: parse
end type config_t

contains

subroutine parse(this, config_string)

    class(config_t),  intent(inout) :: this
    character(len=*), intent(in)    :: config_string

end subroutine parse

end module config_mod

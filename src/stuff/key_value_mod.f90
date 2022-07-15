module key_value_mod

implicit none

type string_t
    character(:), allocatable :: str
end type string_t

type key_value_r8_t
    type(string_t), allocatable :: keys(:)
    real(kind=8),   allocatable :: values(:)

    contains
    procedure, public :: print => print_key_value_r8
end type key_value_r8_t

contains

subroutine print_key_value_r8(this, key_value_delim, line_separator)
    class(key_value_r8_t), intent(in) :: this
    character(len=*), optional, intent(in) :: line_separator
    character(len=*), optional, intent(in) :: key_value_delim

    character(len=:), allocatable :: line_sep_loc, key_value_delim_loc
    character(len=:), allocatable :: result_line
    character(len=256) :: buff
    integer(kind=4) :: i

    line_sep_loc = " "
    key_value_delim_loc = " = "

    result_line = ""

    do i=1,size(this%keys)
        write(buff,*), this%keys(i)%str, key_value_delim_loc, this%values(i), line_sep_loc
        result_line = result_line // trim(buff)
    end do

    i = len(result_line)
    write(buff,"(A,I8.8,A)") "(A",i,")"
    print trim(buff), result_line

end subroutine print_key_value_r8

end module key_value_mod

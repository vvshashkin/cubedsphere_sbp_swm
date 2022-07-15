module cmd_args_mod

    type cmd_arg_t
        character(:), allocatable :: str
    end type cmd_arg_t

contains
    subroutine get_cmd_args(arg_list,arg_list_len)
        type(cmd_arg_t), allocatable, intent(out) :: arg_list(:)
        integer(kind=4), optional,    intent(out) :: arg_list_len

        integer(kind=4) i, nargs, arg_len

        nargs = command_argument_count()+1 !c style indexing -> fortran indexing
                                           !(0-th arg is program name)
        if(present(arg_list_len)) arg_list_len = nargs

        allocate(arg_list(nargs))

        do i=1, nargs
            call get_command_argument(i-1, length=arg_len)
            allocate(character(arg_len) :: arg_list(i)%str)
            call get_command_argument(i-1, arg_list(i)%str)
        end do

    end subroutine get_cmd_args

end module cmd_args_mod

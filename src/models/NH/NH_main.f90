program NH_main

    use parcomm_mod,  only : init_global_parallel_enviroment,   &
                             deinit_global_parallel_enviroment, &
                             parcomm_global
    use cmd_args_mod, only : cmd_arg_t, get_cmd_args
    use namelist_read_mod, only : read_namelist_as_str

    use nh_model_mod,         only : nh_model_t
    use config_nh_model_mod,  only : config_nh_model_t
    use nh_model_factory_mod, only : create_nh_model

    implicit none

    type(cmd_arg_t),  allocatable :: cmd_args(:)
    integer(kind=4) nargs
    character(len=:), allocatable :: namelist_file
    character(:), allocatable :: namelist_string
    logical :: read_namelist

    type(config_nh_model_t) :: config
    type(nh_model_t) :: nh_model

    call init_global_parallel_enviroment()

    call get_cmd_args(cmd_args, nargs)

    read_namelist = .false.
    if(nargs >= 2) then
        namelist_file = cmd_args(2)%str
        call read_namelist_as_str(namelist_string, namelist_file, parcomm_global%myid)
        read_namelist = .true.
    end if

    if(.not. read_namelist) then
        call parcomm_global%abort("NH_main usage: NH_main namelist=namelist_file_name")
    end if

    call config%parse(namelist_string)
    call create_nh_model(nh_model,config)

    call nh_model%run()

    call deinit_global_parallel_enviroment()

end program NH_main

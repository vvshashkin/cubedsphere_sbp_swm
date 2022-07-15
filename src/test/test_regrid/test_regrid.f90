program main

    use test_regrid_mod, only : test_regrid, test_regrid_vec
    use parcomm_mod,     only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_regrid('A')
    call test_regrid('Ah')

    call test_regrid_vec('Ah','covariant')
    call test_regrid_vec('Ah','contravariant')
    call test_regrid_vec('A','covariant')
    call test_regrid_vec('A','contravariant')
    call test_regrid_vec('C','covariant')
    call test_regrid_vec('C','contravariant')

    call deinit_global_parallel_enviroment()

end

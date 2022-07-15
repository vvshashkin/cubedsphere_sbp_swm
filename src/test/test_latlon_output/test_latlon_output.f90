program latlon_output_test

    use test_latlon_output_mod,  only : test_latlon_output, test_latlon_vec_output
    use parcomm_mod,             only : init_global_parallel_enviroment, &
                                        deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_latlon_output(staggering="A", scalar_grid="A",irec=1)
    call test_latlon_output(staggering="Ah",scalar_grid="Ah",irec=2)
    call test_latlon_vec_output(staggering="Ah", components_type="contravariant",irec=1)
    call test_latlon_vec_output(staggering="A", components_type="contravariant",irec=2)
    call test_latlon_vec_output(staggering="C", components_type="covariant",irec=3)

    call deinit_global_parallel_enviroment()

end program latlon_output_test

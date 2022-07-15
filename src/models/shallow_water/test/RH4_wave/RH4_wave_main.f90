program RH4_wave_test

use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment, parcomm_global

use RH4_wave_mod, only : run_RH4_wave

call init_global_parallel_enviroment()


call run_RH4_wave()


call deinit_global_parallel_enviroment()

end program RH4_wave_test

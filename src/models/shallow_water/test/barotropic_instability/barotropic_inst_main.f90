program barotropic_inst_test

use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment, parcomm_global

use barotropic_inst_mod, only : run_barotropic_inst

call init_global_parallel_enviroment()


call run_barotropic_inst()


call deinit_global_parallel_enviroment()

end program barotropic_inst_test

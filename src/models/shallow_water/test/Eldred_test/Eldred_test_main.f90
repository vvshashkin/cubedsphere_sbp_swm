program Eldred_test

use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment, parcomm_global

use Eldred_test_mod, only : run_Eldred_test

call init_global_parallel_enviroment()


call run_Eldred_test()


call deinit_global_parallel_enviroment()

end program Eldred_test

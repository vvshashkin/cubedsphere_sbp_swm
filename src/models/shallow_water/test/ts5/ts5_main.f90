program ts5_test

use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment, parcomm_global

use ts5_mod, only : run_ts5

call init_global_parallel_enviroment()

call run_ts5()

call deinit_global_parallel_enviroment()

end program ts5_test

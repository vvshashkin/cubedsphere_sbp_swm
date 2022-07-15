program ts2_test

use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment, parcomm_global

use ts2_mod, only : run_ts2

call init_global_parallel_enviroment()

call run_ts2()

call deinit_global_parallel_enviroment()

end program ts2_test

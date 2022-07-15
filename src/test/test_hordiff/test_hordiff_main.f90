program horidff_test

use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment, parcomm_global

use test_hordiff_scalar_mod, only : hordiff_scalar_test
use test_hordiff_vector_mod, only : hordiff_vector_test

call init_global_parallel_enviroment()


call hordiff_scalar_test()
! call hordiff_vector_test()


call deinit_global_parallel_enviroment()

end program horidff_test

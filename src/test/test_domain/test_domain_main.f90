program domain_test

use test_domain_mod, only : test_domain
use parcomm_mod,     only : init_global_parallel_enviroment, &
                            deinit_global_parallel_enviroment

implicit none

integer(kind=4) :: ierr

    call init_global_parallel_enviroment()

    call test_domain()

    call deinit_global_parallel_enviroment()

end program domain_test

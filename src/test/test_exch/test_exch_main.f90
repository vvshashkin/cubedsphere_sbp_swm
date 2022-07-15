program main

use test_mod,    only : test_A_halo_exchange, test_Ah_halo_exchange, &
                        test_halo_vec_C_exchange, test_gather_exchange, &
                        test_halo_vec_Ch_exchange
use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_A_halo_exchange()
    ! call test_Ah_halo_exchange()
     call test_halo_vec_Ch_exchange()
    call test_halo_vec_C_exchange()
    call test_gather_exchange()

    call deinit_global_parallel_enviroment()

end

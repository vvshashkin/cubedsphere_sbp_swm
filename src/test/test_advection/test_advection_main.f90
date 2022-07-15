program advection_test

use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment, parcomm_global

use test_solid_rotation_mod, only : test_solid_rotation

call init_global_parallel_enviroment()


call test_solid_rotation()


call deinit_global_parallel_enviroment()

end program advection_test

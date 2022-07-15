program test_lap
!prints maximum/minimum real and imag parts of div*grad operator eigen values

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment
use test_diffops_mod, only: test_laplace_spectre

use cmd_args_mod,     only: cmd_arg_t, get_cmd_args

implicit none

type(cmd_arg_t), allocatable :: cmd_args(:)
integer(kind=4)              :: nargs, iarg
character(:), allocatable    :: div_op_name
character(:), allocatable    :: grad_op_name
character(:), allocatable    :: co2contra_name
character(:), allocatable    :: staggering

div_op_name    = "divergence_ah_c_sbp21"
grad_op_name   = "gradient_ah_c21_sbp_ecs"
co2contra_name = "co2contra_ah_c_sbp21"
staggering     = "Ah_C"

call get_cmd_args(cmd_args,nargs)

if(nargs >=2) div_op_name = cmd_args(2)%str
if(nargs >=3) grad_op_name = cmd_args(3)%str
if(nargs >=4) co2contra_name = cmd_args(4)%str
if(nargs >=5) staggering = cmd_args(5)%str

call init_global_parallel_enviroment()

call test_laplace_spectre(div_op_name, grad_op_name, co2contra_name, staggering)

call deinit_global_parallel_enviroment()

end program test_lap

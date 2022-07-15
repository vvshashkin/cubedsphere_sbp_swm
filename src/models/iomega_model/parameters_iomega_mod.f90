module parameters_iomega_mod

use container_abstract_mod, only: model_parameters_abstract_t

implicit none

type, extends(model_parameters_abstract_t) :: parameters_iomega_t
    integer(kind=4)              :: N     !vector dimension
    complex(kind=8), allocatable :: omega(:) !eigen values, imag == oscillation frequency,
                                             !              real == amplification/decay
end type parameters_iomega_t

contains

function init_iomega_params(omega) result(params)
    type(parameters_iomega_t)       :: params
    complex(kind=8),    intent(in)  :: omega(:)

    params%N = size(omega)
    params%omega = omega

end function init_iomega_params

end module parameters_iomega_mod

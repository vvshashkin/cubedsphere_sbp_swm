module container_abstract_mod

implicit none

!Empty abstract container-types for model state & model parameters
    type, abstract :: state_abstract_t

    end type state_abstract_t

    type, abstract :: model_parameters_abstract_t

    end type model_parameters_abstract_t

end module container_abstract_mod

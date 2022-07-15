module halo_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

implicit none

type, abstract :: halo_t
    contains
    procedure(get_halo_scalar),  deferred, public :: get_halo_scalar  !scalar halo procedure
end type halo_t

type, abstract :: halo_vec_t
    contains
    procedure(get_halo_vector),  deferred, public :: get_halo_vector  !vector halo procedure
end type halo_vec_t

interface
    subroutine get_halo_scalar(this,f,domain,halo_width)
        import halo_t, grid_field_t, domain_t
        class(halo_t),       intent(inout) :: this
        class(grid_field_t), intent(inout) :: f
        type(domain_t),      intent(in)    :: domain
        integer(kind=4),     intent(in)    :: halo_width
    end subroutine get_halo_scalar

    subroutine get_halo_vector(this,u,v,domain,halo_width)
        import halo_vec_t, grid_field_t, domain_t
        class(halo_vec_t),   intent(inout) :: this
        class(grid_field_t), intent(inout) :: u, v
        type(domain_t),      intent(in)    :: domain
        integer(kind=4),     intent(in)    :: halo_width
    end subroutine get_halo_vector
end interface

end module halo_mod

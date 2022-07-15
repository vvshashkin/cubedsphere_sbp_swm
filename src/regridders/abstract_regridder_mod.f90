module abstract_regrid_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t


implicit none



type, abstract :: regrid_t
    contains
        procedure(do_regrid_scalar_interface), deferred :: do_regrid
end type regrid_t

type, abstract :: regrid_vec_t
    contains
        procedure(do_regrid_vector_interface), deferred :: do_regrid
end type regrid_vec_t

interface
    subroutine do_regrid_scalar_interface(this,fout,f,domain)
        import grid_field_t, regrid_t, domain_t
        class(regrid_t),     intent(inout) :: this
        real(kind=8),        intent(inout) :: fout(:,:,:)
        type(grid_field_t),  intent(in)    :: f
        type(domain_t),      intent(in)    :: domain
    end subroutine do_regrid_scalar_interface

    subroutine do_regrid_vector_interface(this,uout,vout,u,v,domain)
        import grid_field_t, regrid_vec_t, domain_t
        class(regrid_vec_t), intent(inout) :: this
        real(kind=8),        intent(inout) :: uout(:,:,:), vout(:,:,:)
        type(grid_field_t),  intent(in)    :: u, v
        type(domain_t),      intent(in)    :: domain
    end subroutine do_regrid_vector_interface
end interface

end module abstract_regrid_mod

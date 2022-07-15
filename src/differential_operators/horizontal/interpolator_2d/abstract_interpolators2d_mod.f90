module abstract_interpolators2d_mod

use domain_mod,     only : domain_t
use grid_field_mod, only : grid_field_t

implicit none

type, abstract :: interpolator2d_scalar2vec_t
contains
    procedure(interp_p2uv), deferred :: interp2d_scalar2vec
end type interpolator2d_scalar2vec_t

type, abstract :: interpolator2d_vec2vec_t
contains
    procedure(interp_uv2uv), deferred :: interp2d_vec2vec
end type interpolator2d_vec2vec_t

type, abstract :: interpolator2d_scalar2scalar_t
contains
    procedure(interp_p2p), deferred :: interp2d_p2p
end type interpolator2d_scalar2scalar_t

abstract interface
    subroutine interp_p2uv(this, u, v, p, domain)
        import domain_t, grid_field_t, interpolator2d_scalar2vec_t
        class(interpolator2d_scalar2vec_t), intent(inout) :: this
        type(grid_field_t),                 intent(inout) :: p, u, v
        type(domain_t),                     intent(in)    :: domain
    end subroutine interp_p2uv
    subroutine interp_uv2uv(this, u, v, u_source, v_source, domain)
        import domain_t, grid_field_t, interpolator2d_vec2vec_t
        class(interpolator2d_vec2vec_t), intent(inout) :: this
        type(grid_field_t),              intent(inout) :: u_source, v_source, u, v
        type(domain_t),                  intent(in)    :: domain
    end subroutine interp_uv2uv
    subroutine interp_p2p(this, f_out, f_in, domain)
        import domain_t, grid_field_t, interpolator2d_scalar2scalar_t
        class(interpolator2d_scalar2scalar_t), intent(inout) :: this
        type(grid_field_t),              intent(inout) :: f_out, f_in
        type(domain_t),                  intent(in)    :: domain
    end subroutine interp_p2p
end interface

end module abstract_interpolators2d_mod

module test_fields_mod
use grid_field_mod, only : grid_field_t
use mesh_mod,       only : mesh_t

implicit none

private
public :: set_scalar_test_field, set_vector_test_field, set_perp_vector_test_field
public :: xyz_scalar_field_generator_t, xyz_scalar_field_generator
public :: solid_rotation_field_generator_t, solid_rotation_field_generator
public :: solid_rotation_vecadv_tend_t, solid_rotation_vecadv_tend
public :: xyz_grad_generator_t, xyz_grad_generator
public :: cross_polar_flow_generator_t, cross_polar_flow_generator
public :: cross_polar_flow_div_generator_t, cross_polar_flow_div_generator
public :: cross_polar_flow_zeta_generator_t, cross_polar_flow_zeta_generator
public :: random_vector_field_generator_t, random_vector_field_generator
public :: VSH_curl_free_10_generator_t, VSH_curl_free_10_generator
public :: zero_scalar_field_generator_t, zero_scalar_field_generator
public :: random_scalar_field_generator_t, random_scalar_field_generator

public :: gaussian_hill_scalar_field_generator_t, gaussian_hill_scalar_field_generator

public :: coriolis_force_field_generator_t!, coriolis_force_field_generator
public :: solid_rotation_t

public :: solid_rotated_scalar_field_t

public :: ts2_height_generator_t
public :: rh4_wave_height_generator_t, rh4_wave_wind_generator_t
public :: barotropic_instability_wind_generator_t
public :: barotropic_instability_height_generator_t
public :: Eldred_test_height_generator_t, Eldred_test_height_generator
public :: Eldred_test_wind_generator_t, Eldred_test_wind_generator
public :: ts5_orography_generator_t

public :: KE_scalar_field_t, Ylm2_field_generator_t

!!!!!!!!!!!!!Abstract scalar and vector fields generators
type, public, abstract :: scalar_field_generator_t
    real(kind=8) :: time = 0.0_8
contains
    procedure(get_scalar_field), deferred :: get_scalar_field
end type scalar_field_generator_t

type, public, abstract :: vector_field_generator_t
    real(kind=8) :: time = 0.0_8
contains
    procedure(get_vector_field), deferred :: get_vector_field
end type vector_field_generator_t

!!!!!!!!!!!!Specific fields
type, extends(scalar_field_generator_t) :: xyz_scalar_field_generator_t
contains
    procedure :: get_scalar_field => generate_xyz_scalar_field
end type xyz_scalar_field_generator_t

type, extends(scalar_field_generator_t) :: random_scalar_field_generator_t
contains
    procedure :: get_scalar_field => generate_random_scalar_field
end type random_scalar_field_generator_t

type, extends(scalar_field_generator_t) :: gaussian_hill_scalar_field_generator_t
contains
    procedure :: get_scalar_field => generate_gaussian_hill_scalar_field
end type gaussian_hill_scalar_field_generator_t

type, extends(vector_field_generator_t) :: solid_rotation_field_generator_t
contains
    procedure :: get_vector_field => generate_solid_rotation_vector_field
end type solid_rotation_field_generator_t

type, extends(vector_field_generator_t) :: solid_rotation_vecadv_tend_t
contains
    procedure :: get_vector_field => generate_solid_rotation_vecadv_tend
end type solid_rotation_vecadv_tend_t

type, extends(vector_field_generator_t) :: coriolis_force_field_generator_t
    class(vector_field_generator_t), allocatable :: input_field
contains
    procedure :: get_vector_field => generate_coriolis_force_field
end type coriolis_force_field_generator_t

type, extends(vector_field_generator_t) :: xyz_grad_generator_t
contains
    procedure :: get_vector_field => generate_xyz_grad_field
end type xyz_grad_generator_t

type, extends(vector_field_generator_t) :: cross_polar_flow_generator_t
contains
    procedure :: get_vector_field => generate_cross_polar_flow
end type cross_polar_flow_generator_t

type, extends(scalar_field_generator_t) :: cross_polar_flow_div_generator_t
contains
    procedure :: get_scalar_field => generate_cross_polar_flow_div
end type cross_polar_flow_div_generator_t

type, extends(scalar_field_generator_t) :: cross_polar_flow_zeta_generator_t
contains
    procedure :: get_scalar_field => generate_cross_polar_flow_zeta
end type cross_polar_flow_zeta_generator_t

type, extends(vector_field_generator_t) :: random_vector_field_generator_t
contains
    procedure :: get_vector_field => generate_random_vector_field
end type random_vector_field_generator_t

type, extends(vector_field_generator_t) :: VSH_curl_free_10_generator_t
contains
    procedure :: get_vector_field => generate_VSH_curl_free_10
end type VSH_curl_free_10_generator_t

type, extends(scalar_field_generator_t) :: zero_scalar_field_generator_t
contains
    procedure :: get_scalar_field => generate_zero_scalar_field
end type zero_scalar_field_generator_t

type, extends(scalar_field_generator_t) :: Ylm2_field_generator_t
contains
    procedure :: get_scalar_field => generate_Ylm2_scalar_field
end type Ylm2_field_generator_t

type, extends(vector_field_generator_t) :: solid_rotation_t
    real(kind=8) :: alpha = 0.0_8 !rotation axis angle
    real(kind=8) :: axis(3) = [0.0_8,0.0_8,1.0_8]
    real(kind=8) :: u0 = 1.0_8    ! equator speed
contains
    procedure :: get_vector_field => gen_solid_rotation_vec_field
end type solid_rotation_t

type, extends(scalar_field_generator_t) :: solid_rotated_scalar_field_t
    type(solid_rotation_t),          allocatable :: solid_rotation_vector_field
    class(scalar_field_generator_t), allocatable :: scalar_field
contains
    procedure :: get_scalar_field => generate_solid_rotated_scalar_field
end type solid_rotated_scalar_field_t

type, extends(scalar_field_generator_t) :: KE_scalar_field_t
    class(vector_field_generator_t), allocatable :: vector_field
contains
    procedure :: get_scalar_field => generate_KE_scalar_field
end type KE_scalar_field_t

type, extends(scalar_field_generator_t) :: ts2_height_generator_t
    real(kind=8) :: h_mean  = 1.0_8
    real(kind=8) :: alpha   = 0.0_8 !rotation axis angle
    real(kind=8) :: axis(3) = [0.0_8, 0.0_8, 1.0_8]
    real(kind=8) :: u0      = 1.0_8 !equator speed
    real(kind=8) :: omega   = 1.0_8 !angular velocity of the sphere
    real(kind=8) :: a       = 1.0_8 !radii of the sphere
    real(kind=8) :: grav    = 1.0_8 !gravity acceleration
contains
    procedure :: get_scalar_field => generate_ts2_height_field
end type ts2_height_generator_t

type, extends(scalar_field_generator_t) :: rh4_wave_height_generator_t
    real(kind=8) :: h_mean = 1.0_8
    real(kind=8) :: omega  = 1.0_8 !angular velocity of the sphere
    real(kind=8) :: a      = 1.0_8 !radii of the sphere
    real(kind=8) :: grav   = 1.0_8 !gravity acceleration
contains
    procedure :: get_scalar_field => generate_rh4_wave_height_field
end type rh4_wave_height_generator_t
type, extends(vector_field_generator_t) :: rh4_wave_wind_generator_t
    real(kind=8) :: omega  = 1.0_8 !angular velocity of the sphere
    real(kind=8) :: a      = 1.0_8 !radii of the sphere
contains
    procedure :: get_vector_field => generate_rh4_wave_wind_field
end type rh4_wave_wind_generator_t

type, extends(scalar_field_generator_t) :: barotropic_instability_height_generator_t
    real(kind=8)    :: H0
    real(kind=8)    :: H_north  !height at North pole
    integer(kind=4) :: Nq  !number of pre-computed latitudes of zonal symmetric h
    real(kind=8), allocatable :: H_zonal(:)
    real(kind=8)    :: dphi !step between precomputed latitudes
    real(kind=8)    :: h_pert !height of perturbation
contains
    procedure :: get_scalar_field => generate_barotropic_instability_height
end type barotropic_instability_height_generator_t

type, extends(vector_field_generator_t) :: barotropic_instability_wind_generator_t
    real(kind=8) :: u0 = 1.0_8
contains
    procedure :: get_vector_field => generate_barotropic_instability_wind
end type barotropic_instability_wind_generator_t

type, extends(scalar_field_generator_t) :: Eldred_test_height_generator_t
    real(kind=8) :: h_mean  = 1e3_8
    real(kind=8) :: delta_h = 8e3_8 / 9.80616_8
contains
    procedure :: get_scalar_field => generate_Eldred_test_height
end type Eldred_test_height_generator_t

type, extends(vector_field_generator_t) :: Eldred_test_wind_generator_t
    real(kind=8)    :: u0 = 120.0_8
    integer(kind=4) :: m  = 12 !number of periods pole to pole
contains
    procedure :: get_vector_field => generate_Eldred_test_wind
end type Eldred_test_wind_generator_t

type, extends(scalar_field_generator_t) :: ts5_orography_generator_t
    real(kind=8) :: r_mount = 1.0
    real(kind=8) :: h_mount = 0.0_8
    real(kind=8) :: x_mount = 0.0_8
    real(kind=8) :: y_mount = 0.5_8*sqrt(3.0_8)
    real(kind=8) :: z_mount = 0.5_8
contains
    procedure :: get_scalar_field => generate_ts5_orography_field
end type ts5_orography_generator_t

!!!field generator instances
type(xyz_scalar_field_generator_t)     :: xyz_scalar_field_generator
type(solid_rotation_field_generator_t) :: solid_rotation_field_generator
type(solid_rotation_vecadv_tend_t)     :: solid_rotation_vecadv_tend
type(xyz_grad_generator_t)             :: xyz_grad_generator
type(cross_polar_flow_generator_t)     :: cross_polar_flow_generator
type(cross_polar_flow_div_generator_t) :: cross_polar_flow_div_generator
type(cross_polar_flow_zeta_generator_t) :: cross_polar_flow_zeta_generator
type(random_vector_field_generator_t)  :: random_vector_field_generator
type(random_scalar_field_generator_t)  :: random_scalar_field_generator
type(VSH_curl_free_10_generator_t)     :: VSH_curl_free_10_generator
type(zero_scalar_field_generator_t)    :: zero_scalar_field_generator
type(gaussian_hill_scalar_field_generator_t) :: gaussian_hill_scalar_field_generator
type(Eldred_test_height_generator_t)   :: Eldred_test_height_generator
type(Eldred_test_wind_generator_t)     :: Eldred_test_wind_generator
! type(coriolis_force_field_generator_t) :: coriolis_force_field_generator = &


abstract interface
    subroutine get_scalar_field(this,f,npts,nlev,x,y,z)
        import scalar_field_generator_t
        class(scalar_field_generator_t), intent(in)  :: this
        integer(kind=4),                 intent(in)  :: npts, nlev
        real(kind=8),                    intent(in)  :: x(npts), y(npts), z(npts)
        real(kind=8),                    intent(out) :: f(npts,nlev)
    end subroutine get_scalar_field
    subroutine get_vector_field(this,vx,vy,vz,npts,nlev,x,y,z)
        import vector_field_generator_t
        class(vector_field_generator_t),    intent(in)  :: this
        integer(kind=4),                    intent(in)  :: npts, nlev
        real(kind=8),                       intent(in)  :: x(npts), y(npts), z(npts)
        real(kind=8), dimension(npts,nlev), intent(out) :: vx, vy, vz
    end subroutine get_vector_field
end interface

contains

subroutine set_scalar_test_field(f, generator, mesh, halo_width, fill_value)
    type(grid_field_t),              intent(inout) :: f
    class(scalar_field_generator_t), intent(in)    :: generator
    type(mesh_t),                    intent(in)    :: mesh
    integer(kind=4),                 intent(in)    :: halo_width
    real(kind=8),          optional, intent(in)    :: fill_value

    !locals
    integer(kind=4) t, isv, iev, jsv, jev, npoints, klev, ks
    real(kind=8), allocatable :: p(:,:,:)
    real(kind=8), dimension(:,:), allocatable :: px, py, pz

    do t = mesh%ts, mesh%te
        isv = mesh%tile(t)%is-halo_width
        iev = mesh%tile(t)%ie+halo_width
        jsv = mesh%tile(t)%js-halo_width
        jev = mesh%tile(t)%je+halo_width
        ks  = mesh%tile(t)%ks
        npoints = (iev-isv+1)*(jev-jsv+1)
        klev = mesh%tile(t)%ke-mesh%tile(t)%ks+1
        if(present(fill_value)) f%tile(t)%p(:,:,:) = -fill_value
        !suppress temporary array warning (made tmp arrays manually)
        allocate(p(isv:iev,jsv:jev,1:klev))
        allocate(px(isv:iev,jsv:jev),py(isv:iev,jsv:jev),pz(isv:iev,jsv:jev))
        !WORKAROUND: ks
        px(isv:iev,jsv:jev) = mesh%tile(t)%rx(isv:iev,jsv:jev,ks)
        py(isv:iev,jsv:jev) = mesh%tile(t)%ry(isv:iev,jsv:jev,ks)
        pz(isv:iev,jsv:jev) = mesh%tile(t)%rz(isv:iev,jsv:jev,ks)
        call generator%get_scalar_field(p, npoints,klev,px,py,pz)
        f%tile(t)%p(isv:iev,jsv:jev,1:klev) = p(isv:iev,jsv:jev,1:klev)
        deallocate(p,px,py,pz)
    end do
end subroutine set_scalar_test_field

subroutine set_vector_test_field(u, v, generator, mesh_u, mesh_v, halo_width, components_type, fill_value)

    use parcomm_mod, only : parcomm_global

    type(grid_field_t),              intent(inout) :: u, v
    class(vector_field_generator_t), intent(in)    :: generator
    type(mesh_t),                    intent(in)    :: mesh_u, mesh_v
    integer(kind=4),                 intent(in)    :: halo_width
    character(len=*),                intent(in)    :: components_type
    real(kind=8),          optional, intent(in)    :: fill_value

    !locals
    integer(kind=4) :: t, ncomponents

     do t = mesh_u%ts, mesh_u%te
        if(components_type=="contravariant") then
            ncomponents = size(mesh_u%tile(t)%b1,1)
            call set_vector_test_field_1tile_1comp(u%tile(t),generator,mesh_u%tile(t),mesh_u%tile(t)%b1, &
                                                   ncomponents, halo_width,fill_value)
            ncomponents = size(mesh_u%tile(t)%b2,1)
            call set_vector_test_field_1tile_1comp(v%tile(t),generator,mesh_v%tile(t),mesh_v%tile(t)%b2, &
                                                   ncomponents, halo_width,fill_value)
        else if (components_type=="covariant") then
            ncomponents = size(mesh_u%tile(t)%a1,1)
            call set_vector_test_field_1tile_1comp(u%tile(t),generator,mesh_u%tile(t),mesh_u%tile(t)%a1, &
                                                   ncomponents, halo_width,fill_value)
            ncomponents = size(mesh_u%tile(t)%a2,1)
            call set_vector_test_field_1tile_1comp(v%tile(t),generator,mesh_v%tile(t),mesh_v%tile(t)%a2, &
                                                   ncomponents, halo_width,fill_value)
        else
            call parcomm_global%abort("test_functions_mod, set_vector_test_field: unknown components type "// &
                                      components_type//" use covariant or contravariant")
        end if

     end do

end subroutine set_vector_test_field

subroutine set_vector_test_field_1tile_1comp(u,generator,mesh,bvec,ncomp,halo_width,fill_value)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t),              intent(inout) :: u
    class(vector_field_generator_t), intent(in)    :: generator
    type(tile_mesh_t),               intent(in)    :: mesh
    real(kind=8),                    intent(in)    :: bvec(ncomp,mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width, &
                                                                 mesh%js-mesh%halo_width:mesh%je+mesh%halo_width, &
                                                                 mesh%ks:mesh%ke)
    integer(kind=4),                 intent(in)    :: ncomp
    integer(kind=4),                 intent(in)    :: halo_width
    real(kind=8),          optional, intent(in)    :: fill_value

    !locals
    integer(kind=4) :: isv, iev, jsv, jev, i, j, k, k1, klev, npoints, ks
    real(kind=8), dimension(:,:,:), allocatable ::  vx, vy, vz
    real(kind=8), dimension(:,:), allocatable   :: px, py, pz

    if(present(fill_value)) u%p(:,:,:) = fill_value

    isv = mesh%is-halo_width
    iev = mesh%ie+halo_width
    jsv = mesh%js-halo_width
    jev = mesh%je+halo_width
    klev = mesh%ke-mesh%ks+1
    ks = mesh%ks

    allocate(vx(isv:iev,jsv:jev,1:klev),vy(isv:iev,jsv:jev,1:klev),vz(isv:iev,jsv:jev,1:klev))
    allocate(px(isv:iev,jsv:jev),py(isv:iev,jsv:jev),pz(isv:iev,jsv:jev))
    px(isv:iev,jsv:jev) = mesh%rx(isv:iev,jsv:jev,ks)
    py(isv:iev,jsv:jev) = mesh%ry(isv:iev,jsv:jev,ks)
    pz(isv:iev,jsv:jev) = mesh%rz(isv:iev,jsv:jev,ks)
    npoints = (iev-isv+1)*(jev-jsv+1)
    call generator%get_vector_field(vx,vy,vz,npoints,klev,px,py,pz)

    do k=mesh%ks, mesh%ke
        k1 = k-mesh%ks+1
        do j = jsv, jev
            do i = isv, iev
            u%p(i,j,k) = sum([vx(i,j,k1),vy(i,j,k1),vz(i,j,k1)]*bvec(1:3,i,j,ks))
            end do
        end do
    end do
    deallocate(vx,vy,vz)
    deallocate(px,py,pz)
end subroutine set_vector_test_field_1tile_1comp

subroutine set_perp_vector_test_field(u, v, generator, mesh_u, mesh_v, halo_width, components_type, fill_value)

    use parcomm_mod, only : parcomm_global

    type(grid_field_t),              intent(inout) :: u, v
    class(vector_field_generator_t), intent(in)    :: generator
    type(mesh_t),                    intent(in)    :: mesh_u, mesh_v
    integer(kind=4),                 intent(in)    :: halo_width
    character(len=*),                intent(in)    :: components_type
    real(kind=8),          optional, intent(in)    :: fill_value

    !locals
    integer(kind=4) :: t
    integer(kind=4) :: ncomponents

     do t = mesh_u%ts, mesh_u%te
        if(components_type=="contravariant") then
            ncomponents = size(mesh_u%tile(t)%b1,1)
            call set_perp_vector_test_field_1tile_1comp(u%tile(t),generator,mesh_u%tile(t),mesh_u%tile(t)%b1, &
                                                   ncomponents, halo_width,fill_value)
            ncomponents = size(mesh_u%tile(t)%b2,1)
            call set_perp_vector_test_field_1tile_1comp(v%tile(t),generator,mesh_v%tile(t),mesh_v%tile(t)%b2, &
                                                   ncomponents, halo_width,fill_value)
        else if (components_type=="covariant") then
            ncomponents = size(mesh_u%tile(t)%a1,1)
            call set_perp_vector_test_field_1tile_1comp(u%tile(t),generator,mesh_u%tile(t),mesh_u%tile(t)%a1, &
                                                   ncomponents, halo_width,fill_value)
            ncomponents = size(mesh_u%tile(t)%a2,1)
            call set_perp_vector_test_field_1tile_1comp(v%tile(t),generator,mesh_v%tile(t),mesh_v%tile(t)%a2, &
                                                   ncomponents, halo_width,fill_value)
        else
            call parcomm_global%abort("test_functions_mod, set_vector_test_field: unknown components type "// &
                                      components_type//" use covariant or contravariant")
        end if

     end do

end subroutine set_perp_vector_test_field

subroutine set_perp_vector_test_field_1tile_1comp(u,generator,mesh,bvec,ncomp,halo_width,fill_value)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t),              intent(inout) :: u
    class(vector_field_generator_t), intent(in)    :: generator
    type(tile_mesh_t),               intent(in)    :: mesh
    real(kind=8),                    intent(in)    :: bvec(ncomp,mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width, &
                                                                 mesh%js-mesh%halo_width:mesh%je+mesh%halo_width, &
                                                                 mesh%ks:mesh%ke)
    integer(kind=4),                 intent(in)    :: ncomp
    integer(kind=4),                 intent(in)    :: halo_width
    real(kind=8),          optional, intent(in)    :: fill_value

    !locals
    integer(kind=4) :: isv, iev, jsv, jev, i, j, k, k1, klev, npoints, ks
    real(kind=8), dimension(:,:,:), allocatable ::  vx, vy, vz
    real(kind=8), dimension(:,:), allocatable   :: px, py, pz
    real(kind=8) :: nx, ny, nz

    if(present(fill_value)) u%p(:,:,:) = fill_value

    isv = mesh%is-halo_width
    iev = mesh%ie+halo_width
    jsv = mesh%js-halo_width
    jev = mesh%je+halo_width
    ks  = mesh%ks
    klev = mesh%ke-mesh%ks+1

    allocate(vx(isv:iev,jsv:jev,1:klev),vy(isv:iev,jsv:jev,1:klev),vz(isv:iev,jsv:jev,1:klev))
    allocate(px(isv:iev,jsv:jev),py(isv:iev,jsv:jev),pz(isv:iev,jsv:jev))
    px(isv:iev,jsv:jev) = mesh%rx(isv:iev,jsv:jev,ks)
    py(isv:iev,jsv:jev) = mesh%ry(isv:iev,jsv:jev,ks)
    pz(isv:iev,jsv:jev) = mesh%rz(isv:iev,jsv:jev,ks)
    npoints = (iev-isv+1)*(jev-jsv+1)
    call generator%get_vector_field(vx,vy,vz,npoints,klev,px,py,pz)

    do k=mesh%ks, mesh%ke
        k1 = k-mesh%ks+1
        do j = jsv, jev
            do i = isv, iev
            nx = px(i, j) / sqrt(px(i,j)**2+py(i,j)**2+pz(i,j)**2)
            ny = py(i, j) / sqrt(px(i,j)**2+py(i,j)**2+pz(i,j)**2)
            nz = pz(i, j) / sqrt(px(i,j)**2+py(i,j)**2+pz(i,j)**2)
            !u%p(i,j,k) = sum([vx(i,j,k1),vy(i,j,k1),vz(i,j,k1)]*bvec(1:3,i,j))
            u%p(i,j,k) = sum([ny*vz(i,j,k)-nz*vy(i,j,k), &
                             -nx*vz(i,j,k)+nz*vx(i,j,k), &
                              nx*vy(i,j,k)-ny*vx(i,j,k)]*bvec(1:3,i,j,ks))
            end do
        end do
    end do
    deallocate(vx,vy,vz)
    deallocate(px,py,pz)
end subroutine set_perp_vector_test_field_1tile_1comp

subroutine generate_xyz_scalar_field(this,f,npts,nlev,x,y,z)
    class(xyz_scalar_field_generator_t),  intent(in) :: this
    integer(kind=4), intent(in)                  :: npts, nlev
    real(kind=8),    intent(in)                  :: x(npts), y(npts), z(npts)
    real(kind=8),    intent(out)                 :: f(npts,nlev)

    integer(kind=4) :: i, k, k1
    real(kind=8)    :: M(3,3) = reshape([1.0_8, 0.0_8, 0.0_8, &
                                        0.0_8, 1.0_8, 0.0_8,  &
                                        0.0_8, 0.0_8, 1.0_8],[3,3])

    do k = 1, nlev
        k1 = mod(k-1,3)+1 !1,2,3,1,2,3,1,2 etc
        do i=1, npts
            f(i,k) = M(1,k1)*x(i)+M(2,k1)*y(i)+M(3,k1)*z(i) !k1=1:x, k1=2:y, k1=3:z
        end do
    end do
end subroutine generate_xyz_scalar_field

subroutine generate_gaussian_hill_scalar_field(this, f, npts, nlev, x, y, z)

    use const_mod,      only : pi
    use sph_coords_mod, only : cart2sph

    class(gaussian_hill_scalar_field_generator_t),  intent(in)  :: this
    integer(kind=4),                                intent(in)  :: npts, nlev
    real(kind=8),                                   intent(in)  :: x(npts), y(npts), z(npts)
    real(kind=8),                                   intent(out) :: f(npts,nlev)

    integer(kind=4) :: i, k
    real(kind=8)    :: lam, phi, r, gamma

    gamma = 0.0_8

    do k = 1, nlev
        do i=1, npts
            call cart2sph(x(i), y(i), z(i), lam, phi)
            r = acos(sin(phi)*sin(gamma)+cos(phi)*cos(gamma)*(cos(lam)))
            f(i,k) = exp(-(4.0_8*r)**2)
        end do
    end do
end subroutine generate_gaussian_hill_scalar_field

subroutine generate_solid_rotation_vector_field(this,vx,vy,vz,npts,nlev,x,y,z)
    class(solid_rotation_field_generator_t),  intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k, k1
    real(kind=8)    :: axis(3,3) = reshape([1.0_8, 0.0_8, 0.0_8,  & !rotation axes for different levs
                                            0.0_8, 1.0_8, 0.0_8,  &
                                            0.0_8, 0.0_8, 1.0_8],[3,3])

    do k = 1, nlev
        k1 = mod(k-1,3)+1  !1,2,3,1,2,3,1,2 etc
        do i=1, npts
            vx(i,k) = axis(2,k1)*z(i)-axis(3,k1)*y(i)
            vy(i,k) =-axis(1,k1)*z(i)+axis(3,k1)*x(i)
            vz(i,k) = axis(1,k1)*y(i)-axis(2,k1)*x(i)
        end do
    end do
end subroutine generate_solid_rotation_vector_field

subroutine generate_solid_rotation_vecadv_tend(this,vx,vy,vz,npts,nlev,x,y,z)
    class(solid_rotation_vecadv_tend_t),      intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k, k1
    real(kind=8)    :: ux, uy, uz, ra
    real(kind=8)    :: axis(3,3) = reshape([1.0_8, 0.0_8, 0.0_8,  & !rotation axes for different levs
                                            0.0_8, 1.0_8, 0.0_8,  &
                                            0.0_8, 0.0_8, 1.0_8],[3,3])

    do k = 1, nlev
        k1 = mod(k-1,3)+1  !1,2,3,1,2,3,1,2 etc
        do i=1, npts
            ux = axis(2,k1)*z(i)-axis(3,k1)*y(i)
            uy =-axis(1,k1)*z(i)+axis(3,k1)*x(i)
            uz = axis(1,k1)*y(i)-axis(2,k1)*x(i)
            ra = axis(1,k1)*x(i)+axis(2,k1)*y(i)+axis(3,k1)*z(i)
            vx(i,k) =-(y(i)*uz-z(i)*uy)*ra
            vy(i,k) = (x(i)*uz-z(i)*ux)*ra
            vz(i,k) =-(x(i)*uy-y(i)*ux)*ra
        end do
    end do
end subroutine generate_solid_rotation_vecadv_tend

subroutine gen_solid_rotation_vec_field(this, vx, vy, vz, npts, nlev, x, y, z)

    use sph_coords_mod, only : rotate_3D_y, cart2sph, sph2cart_vec

    class(solid_rotation_t),                  intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k
    real(kind=8) :: lam, phi, v_lam, v_phi, xn, yn, zn, vxn, vyn, vzn
    real(kind=8) :: ax, ay, az, amod

    amod = sqrt(sum(this%axis(1:3)**2))
    ax = this%axis(1) / amod
    ay = this%axis(2) / amod
    az = this%axis(3) / amod

    do k = 1, nlev
        do i=1, npts
            xn = x(i) / sqrt(x(i)**2+y(i)**2+z(i)**2)
            yn = y(i) / sqrt(x(i)**2+y(i)**2+z(i)**2)
            zn = z(i) / sqrt(x(i)**2+y(i)**2+z(i)**2)
            vx(i,k) = this%u0*(ay*zn-az*yn)
            vy(i,k) =-this%u0*(ax*zn-az*xn)
            vz(i,k) = this%u0*(ax*yn-ay*xn)
            ! !translate north pole from [0, pi/2] to [0, pi/2-alpha]
            ! call rotate_3D_y(x(i), y(i), z(i), -this%alpha, xn, yn, zn)
            ! call cart2sph(xn, yn, zn, lam, phi)
            !
            ! v_lam = this%u0*cos(phi)
            ! v_phi = 0.0_8
            !
            ! call sph2cart_vec(lam, phi, v_lam, v_phi, vxn, vyn, vzn)
            ! call rotate_3D_y(vxn, vyn, vzn, this%alpha, vx(i,k), vy(i,k), vz(i,k))
        end do
    end do
end subroutine gen_solid_rotation_vec_field

subroutine generate_coriolis_force_field(this, vx, vy, vz, npts, nlev, x, y, z)

    use sph_coords_mod, only : cart2sph, sph2cart_vec, cart2sph_vec

    class(coriolis_force_field_generator_t),  intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    real(kind=8)    :: v_lam, v_phi, lam, phi, pcori, f_lam, f_phi
    integer(kind=4) :: i, k

    call this%input_field%get_vector_field(vx, vy, vz, npts, nlev, x, y, z)

    do k = 1, nlev
        do i=1, npts
            call cart2sph(x(i), y(i), z(i), lam, phi)
            call cart2sph_vec(lam, phi, vx(i,k), vy(i,k), vz(i,k), v_lam, v_phi)
            pcori = 2*sin(phi)
            f_lam =  pcori*v_phi
            f_phi = -pcori*v_lam
            call sph2cart_vec(lam, phi, f_lam, f_phi, vx(i,k), vy(i,k), vz(i,k))
        end do
    end do
end subroutine generate_coriolis_force_field

subroutine generate_xyz_grad_field(this,vx,vy,vz,npts,nlev,x,y,z)
    class(xyz_grad_generator_t),              intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k, k1
    real(kind=8)    :: e(3,3) = reshape([1.0_8, 0.0_8, 0.0_8,  & ! level1: f=x,
                                         0.0_8, 1.0_8, 0.0_8,  & ! level2: f=y
                                         0.0_8, 0.0_8, 1.0_8],[3,3]) !level 3: f=z
    !(nabla)_h f = (nabla)_3d f- n *(n, nabla_3d f)
    !(nabla)_3d f = e(k)
    real(kind=8)    :: n(3) !normal to the surface
    real(kind=8)    :: nne(3) !n *(n, nabla_3d f)
    real(kind=8)    :: nabla_h(3)

    do k = 1, nlev
        k1 = mod(k-1,3)+1  !1,2,3,1,2,3,1,2 etc
        do i=1, npts
            n(1:3) = [x(i),y(i),z(i)] / sqrt(x(i)**2+y(i)**2+z(i)**2)
            nne(1:3) = n(1:3)*sum(n(1:3)*e(1:3,k))
            nabla_h(1:3) = e(1:3,k)-nne(1:3)
            vx(i,k) = nabla_h(1)
            vy(i,k) = nabla_h(2)
            vz(i,k) = nabla_h(3)
        end do
    end do
end subroutine generate_xyz_grad_field

subroutine generate_cross_polar_flow(this,vx,vy,vz,npts,nlev,x,y,z)
    use latlon_functions_mod, only : sin_phi, cos_phi, sin_lam, cos_lam

    class(cross_polar_flow_generator_t),      intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k, k1
    integer(kind=4) :: k3
    real(kind=8)    :: x1, y1, z1
    real(kind=8)    :: vx1, vy1, vz1
    real(kind=8)    :: cphi, sphi, clam, slam
    real(kind=8)    :: u, v
    real(kind=8)    :: M(3,3) = reshape([1.0_8, 0.0_8, 0.0_8,  & !rotation matrix
                                         0.0_8, 1.0_8, 0.0_8,  &
                                         0.0_8, 0.0_8, 1.0_8],[3,3])

    k3(k) = mod(k-1,3)+1

    do k = 1, nlev
        do i=1, npts
            !permutate coordinates k=1: (x,y,z), k=2: (y,z,x), k=3: (z,x,y)
            !x1 = sum(M(1:3,k3(k))  *[x(i),y(i),z(i)])
            !y1 = sum(M(1:3,k3(k+1))*[x(i),y(i),z(i)])
            !z1 = sum(M(1:3,k3(k+2))*[x(i),y(i),z(i)])
            x1 = x(i); y1 = y(i); z1 = z(i)
            cphi = cos_phi(x1,y1,z1)
            sphi = sin_phi(x1,y1,z1)
            clam = cos_lam(x1,y1,z1)
            slam = sin_lam(x1,y1,z1)
            u = sphi*(sphi**2-3.0_8*cphi**2)*slam-0.5_8*cphi
            v = sphi**2*clam
            vx1 =-slam*u-sphi*clam*v
            vy1 = clam*u-sphi*slam*v
            vz1 = cphi*v
            vx(i,k) = vx1; vy(i,k) = vy1; vz(i,k) = vz1
            !permute vector components back
            !vx(i,k) = sum(M(k3(k),  1:3)  *[vx1,vy1,vz1])
            !vy(i,k) = sum(M(k3(k+1),1:3)  *[vx1,vy1,vz1])
            !vz(i,k) = sum(M(k3(k+2),1:3)  *[vx1,vy1,vz1])
        end do
    end do
end subroutine generate_cross_polar_flow

subroutine generate_cross_polar_flow_div(this,f,npts,nlev,x,y,z)
    import scalar_field_generator_t
    class(cross_polar_flow_div_generator_t),  intent(in) :: this
    integer(kind=4), intent(in)                  :: npts, nlev
    real(kind=8),    intent(in)                  :: x(npts), y(npts), z(npts)
    real(kind=8),    intent(out)                 :: f(npts,nlev)

    integer(kind=4) :: i, k

    do k = 1, nlev
        do i=1, npts
            f(i,k) = -x(i)*z(i) / (x(i)**2+y(i)**2+z(i)**2)
        end do
    end do
end subroutine generate_cross_polar_flow_div

subroutine generate_cross_polar_flow_zeta(this,f,npts,nlev,x,y,z)
    import scalar_field_generator_t
    class(cross_polar_flow_zeta_generator_t),  intent(in) :: this
    integer(kind=4), intent(in)                  :: npts, nlev
    real(kind=8),    intent(in)                  :: x(npts), y(npts), z(npts)
    real(kind=8),    intent(out)                 :: f(npts,nlev)

    integer(kind=4) :: i, k

    do k = 1, nlev
        do i=1, npts
            f(i,k) = (16._8*(1-z(i)**2)*y(i)-13._8*y(i)-z(i)) / (x(i)**2+y(i)**2+z(i)**2)
        end do
    end do
end subroutine generate_cross_polar_flow_zeta

subroutine generate_random_vector_field(this,vx,vy,vz,npts,nlev,x,y,z)
    use latlon_functions_mod, only : sin_phi, cos_phi, sin_lam, cos_lam

    class(random_vector_field_generator_t),   intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k
    real(kind=8)    :: proj

    call random_number(vx)
    call random_number(vy)
    call random_number(vz)

    vx = vx-0.5_8
    vy = vy-0.5_8
    vz = vz-0.5_8

    !ensure the vectors are tangent to the sphere
!    do k = 1, nlev
!        do i=1, npts
!            proj = (x(i)*vx(i)+y(i)*vy(i)+z(i)*vz(i)) / (x(i)**2+y(i)**2+z(i)**2)
!            vx(i) = vx(i)-proj*x(i)
!            vy(i) = vy(i)-proj*y(i)
!            vz(i) = vz(i)-proj*z(i)
!        end do
!    end do
end subroutine generate_random_vector_field
subroutine generate_VSH_curl_free_10(this, vx, vy, vz, npts, nlev, x, y, z)

    use sph_coords_mod, only : cart2sph, sph2cart_vec

    class(VSH_curl_free_10_generator_t),  intent(in)  :: this
    integer(kind=4),                      intent(in)  :: npts, nlev
    real(kind=8), dimension(npts),        intent(in)  :: x, y, z
    real(kind=8), dimension(npts, nlev),  intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k
    real(kind=8)    :: lam, phi, v_lam, v_phi

    do k = 1, nlev
        do i = 1, npts
            call cart2sph(x(i), y(i), z(i), lam, phi)
            v_lam = 0.0_8
            v_phi = -cos(phi)
            call sph2cart_vec(lam, phi, v_lam, v_phi, vx(i,k), vy(i,k), vz(i,k))
        end do
    end do

end subroutine generate_VSH_curl_free_10
subroutine generate_solid_rotated_scalar_field(this, f, npts, nlev, x, y, z)

    use sph_coords_mod, only : rotate_3D_y, cart2sph, sph2cart_vec, sph2cart

    class(solid_rotated_scalar_field_t), intent(in)  :: this
    integer(kind=4),                     intent(in)  :: npts, nlev
    real(kind=8),                        intent(in)  :: x(npts), y(npts), z(npts)
    real(kind=8),                        intent(out) :: f(npts,nlev)

    integer(kind=4) :: i, k
    real(kind=8) :: lam, phi, v_lam, v_phi, xn(npts), yn(npts), zn(npts)
    real(kind=8) :: alpha, u0

    alpha = this%solid_rotation_vector_field%alpha
    u0    = this%solid_rotation_vector_field%u0

    do i=1, npts
        !translate north pole from [0, pi/2] to [0, pi/2-alpha]
        call rotate_3D_y(x(i), y(i), z(i), -alpha, xn(i), yn(i), zn(i))
        call cart2sph(xn(i), yn(i), zn(i), lam, phi)

        lam =  lam - u0*this%time

        call sph2cart(lam, phi, xn(i), yn(i), zn(i))
    end do

        call this%scalar_field%get_scalar_field(f, npts, nlev, xn, yn, zn)

end subroutine generate_solid_rotated_scalar_field
subroutine generate_KE_scalar_field(this, f, npts, nlev, x, y, z)

    use sph_coords_mod, only : rotate_3D_y, cart2sph, sph2cart_vec, sph2cart

    class(KE_scalar_field_t), intent(in)  :: this
    integer(kind=4),          intent(in)  :: npts, nlev
    real(kind=8),             intent(in)  :: x(npts), y(npts), z(npts)
    real(kind=8),             intent(out) :: f(npts,nlev)

    real(kind=8):: vx(npts,nlev), vy(npts,nlev), vz(npts,nlev)

    integer(kind=4) :: i, k

    call this%vector_field%get_vector_field(vx, vy,vz, npts, nlev, x, y, z)

    do k = 1, nlev
        do i = 1, npts
            f(i,k) = 0.5_8*(vx(i,k)**2+vy(i,k)**2+vz(i,k)**2)
        end do
    end do

end subroutine generate_KE_scalar_field
subroutine generate_ts2_height_field(this, f, npts, nlev, x, y, z)

    use sph_coords_mod, only : rotate_3D_y, cart2sph, sph2cart_vec, sph2cart

    class(ts2_height_generator_t), intent(in)  :: this
    integer(kind=4),               intent(in)  :: npts, nlev
    real(kind=8),                  intent(in)  :: x(npts), y(npts), z(npts)
    real(kind=8),                  intent(out) :: f(npts,nlev)

    integer(kind=4) :: i, k
    real(kind=8) :: lam, phi, v_lam, v_phi, xn(npts), yn(npts), zn(npts)
    real(kind=8) :: alpha, u0, h_mean, a, omega, grav, axis(3)

    alpha     = this%alpha
    axis(1:3) = this%axis(1:3) / sqrt(sum(this%axis(1:3)**2))
    u0        = this%u0
    h_mean    = this%h_mean
    omega     = this%omega
    a         = this%a
    grav      = this%grav
    do k = 1, nlev
        do i=1, npts
            ! !translate north pole from [0, pi/2] to [0, pi/2-alpha]
            ! call rotate_3D_y(x(i), y(i), z(i), -alpha, xn(i), yn(i), zn(i))
            ! call cart2sph(xn(i), yn(i), zn(i), lam, phi)
            xn(1) = x(i) / sqrt(x(i)**2+y(i)**2+z(i)**2)
            yn(1) = y(i) / sqrt(x(i)**2+y(i)**2+z(i)**2)
            zn(1) = z(i) / sqrt(x(i)**2+y(i)**2+z(i)**2)
            phi = asin(xn(1)*axis(1)+yn(1)*axis(2)+zn(1)*axis(3))
            f(i,k) = h_mean - (a*omega*u0+0.5_8*u0**2)/grav*sin(phi)**2
        end do
    end do

end subroutine generate_ts2_height_field
subroutine generate_rh4_wave_wind_field(this, vx, vy, vz, npts, nlev, x, y, z)

    use sph_coords_mod, only : rotate_3D_y, cart2sph, sph2cart_vec

    class(rh4_wave_wind_generator_t),         intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k
    real(kind=8) :: lam, phi, v_lam, v_phi, xn, yn, zn, vxn, vyn, vzn
    real(kind=8) :: a, omega
    real(kind=8) :: rW, rK, R

    omega  = this%omega
    a      = this%a

    rW = 7.848*10.0_8**(-6)
    rK = 7.848*10.0_8**(-6)
    R  = 4

    do k = 1, nlev
        do i=1, npts
            !translate north pole from [0, pi/2] to [0, pi/2-alpha]
            ! call rotate_3D_y(x(i), y(i), z(i), -this%alpha, xn, yn, zn)
            ! call cart2sph(xn, yn, zn, lam, phi)

            call cart2sph(x(i), y(i), z(i), lam, phi)

            v_lam = a*(rW*cos(phi)+ &
                      +rK*(cos(phi)**(R-1)) &
                      *(R*sin(phi)**2-cos(phi)**2)*cos(R*lam))
            v_phi = -a*rK*R*(cos(phi)**(R-1))*sin(phi)*sin(R*lam)

            call sph2cart_vec(lam, phi, v_lam, v_phi, vx(i,k), vy(i,k), vz(i,k))
        end do
    end do
end subroutine generate_rh4_wave_wind_field
subroutine generate_rh4_wave_height_field(this, f, npts, nlev, x, y, z)

    use sph_coords_mod, only : rotate_3D_y, cart2sph, sph2cart_vec, sph2cart

    class(rh4_wave_height_generator_t), intent(in)  :: this
    integer(kind=4),                    intent(in)  :: npts, nlev
    real(kind=8),                       intent(in)  :: x(npts), y(npts), z(npts)
    real(kind=8),                       intent(out) :: f(npts,nlev)

    integer(kind=4) :: i, k
    real(kind=8) :: lam, phi, v_lam, v_phi, xn(npts), yn(npts), zn(npts)
    real(kind=8) :: h_mean, a, omega, grav
    real(kind=8) :: rA, rB, rC, rW, rK, R

    h_mean = this%h_mean
    omega  = this%omega
    a      = this%a
    grav   = this%grav

    rW = 7.848*10.0_8**(-6)
    rK = 7.848*10.0_8**(-6)
    R  = 4

    do k = 1, nlev
        do i=1, npts
            !translate north pole from [0, pi/2] to [0, pi/2-alpha]
            ! call rotate_3D_y(x(i), y(i), z(i), -alpha, xn(i), yn(i), zn(i))
            ! call cart2sph(xn(i), yn(i), zn(i), lam, phi)

            call cart2sph(x(i), y(i), z(i), lam, phi)

            rA = rW/2*(2*omega+rW)*cos(phi)**2 + &
               + 0.25*rK**2*(cos(phi)**(2*R))*( (R+1)*cos(phi)**2 + &
                 (2*R**2-R-2) -2*R**2*cos(phi)**(-2) )
            rB = 2*((omega+rW)*rK/(R+1)/(R+2))*(cos(phi)**R)*( (R**2+2*R+2) &
               - (R+1)**2*cos(phi)**2 )
            rC = 0.25*rK**2*(cos(phi)**(2*R))*( (R+1)*cos(phi)**2-(R+2) )

            f(i,k) = h_mean + a**2/grav*(rA+rB*cos(R*lam)+rC*cos(2*R*lam))
        end do
    end do

end subroutine generate_rh4_wave_height_field

subroutine generate_barotropic_instability_wind(this, vx, vy, vz, npts, nlev, x, y, z)

    use const_mod, only : pi
    use barotropic_instability_u_mod, only : u_fun=>barotropic_instability_u

    class(barotropic_instability_wind_generator_t),         intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k
    real(kind=8)    :: phi, u, ex, ey
    real(kind=8), parameter :: phi0 = pi/7d0, phi1 = .5d0*pi-phi0

    do k = 1, nlev
        do i=1, npts
            phi = asin(z(i) / sqrt(x(i)**2+y(i)**2+z(i)**2))
            u = u_fun(this%u0, phi)
            ex = -y(i) / max(sqrt(x(i)**2+y(i)**2),1e-14)
            ey =  x(i) / max(sqrt(x(i)**2+y(i)**2),1e-14)

            vx(i,k) = ex*u
            vy(i,k) = ey*u
            vz(i,k) = 0.0_8
        end do
    end do
end subroutine generate_barotropic_instability_wind

subroutine generate_barotropic_instability_height(this, f, npts, nlev, x, y, z)

    use const_mod, only : pi
    use sph_coords_mod, only : cart2sph

    class(barotropic_instability_height_generator_t), intent(in)  :: this
    integer(kind=4),                    intent(in)  :: npts, nlev
    real(kind=8),                       intent(in)  :: x(npts), y(npts), z(npts)
    real(kind=8),                       intent(out) :: f(npts,nlev)

    integer(kind=4) :: i, k, indy
    real(kind=8), parameter :: phi0 = pi/7d0, phi1 = .5d0*pi-phi0
    real(kind=8),  parameter :: lat_diam = 1.0_8/15.0_8, lon_diam = 1.0_8/3.0_8
    real(kind=8) :: phi, lam, zdy, w(-1:2)


    do k = 1, nlev
        do i=1, npts
            call cart2sph(x(i), y(i), z(i), lam, phi)
            if(phi<=phi0) then
                f(i,k) = this%H0
            else if(phi >= phi1) then
                f(i,k) = this%H_north
            else
                zdy = (phi-phi0) / this%dphi
                indy = floor(zdy)
                zdy = zdy - indy
                w(-1) =-zdy*(zdy-1._8)*(zdy-2._8) / 6._8
                w(0) = (zdy+1._8)*(zdy-1._8)*(zdy-2._8) / 2._8
                w(1) =-(zdy+1._8)*zdy*(zdy-2._8) / 2._8
                w(2) = (zdy+1._8)*zdy*(zdy-1._8) / 6._8
                ! f(i,k) = this%H_zonal(indy) + (this%H_zonal(indy+1)-this%H_zonal(indy))*zdy
                f(i,k) = sum(this%H_zonal(indy-1:indy+2)*w(-1:2))
            end if
            if(lam>pi) lam = lam -2.0_8*pi
            f(i,k) = f(i,k) + this%h_pert*cos(phi)* &
                                 exp(-(lam/lon_diam)**2)*exp(-((phi-.25_8*pi)/lat_diam)**2)
        end do
    end do

end subroutine generate_barotropic_instability_height
subroutine generate_zero_scalar_field(this, f, npts, nlev, x, y, z)

    class(zero_scalar_field_generator_t),  intent(in) :: this
    integer(kind=4), intent(in)                       :: npts, nlev
    real(kind=8),    intent(in)                       :: x(npts), y(npts), z(npts)
    real(kind=8),    intent(out)                      :: f(npts,nlev)

    integer(kind=4) :: i, k

    do k = 1, nlev
        do i=1, npts
            f(i,k) = 0.0_8
        end do
    end do

end subroutine generate_zero_scalar_field

subroutine generate_random_scalar_field(this, f, npts, nlev, x, y, z)

    class(random_scalar_field_generator_t),  intent(in) :: this
    integer(kind=4), intent(in)                         :: npts, nlev
    real(kind=8),    intent(in)                         :: x(npts), y(npts), z(npts)
    real(kind=8),    intent(out)                        :: f(npts,nlev)

    integer(kind=4) :: i, k

    call random_number(f)

end subroutine generate_random_scalar_field

subroutine generate_Eldred_test_height(this, f, npts, nlev, x, y, z)

    class(Eldred_test_height_generator_t), intent(in)   :: this
    integer(kind=4),                       intent(in)   :: npts, nlev
    real(kind=8),                          intent(in)   :: x(npts), y(npts), z(npts)
    real(kind=8),                          intent(out)  :: f(npts,nlev)

    integer(kind=4) :: i, k

    do k = 1, nlev
        do i=1, npts
            f(i,k) = this%h_mean+this%delta_h*(1.0_8/3.0_8 - z(i)**2)
        end do
    end do

end subroutine generate_Eldred_test_height

subroutine generate_Eldred_test_wind(this, vx, vy, vz, npts, nlev, x, y, z)

    use const_mod, only : pi

    class(Eldred_test_wind_generator_t),      intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k
    real(kind=8)    :: phi, u_div_cos
    ! real(kind=8), parameter :: phi0 = pi/7d0, phi1 = .5d0*pi-phi0
    !
    do k = 1, nlev
        do i=1, npts
            phi = asin(z(i) / sqrt(x(i)**2+y(i)**2+z(i)**2))
            !u / cos(phi):
            u_div_cos = 4.0_8*this%u0*sin(phi)**2*cos(phi)*(sin(this%m*phi)**2-0.5_8)
            vx(i,k) =-u_div_cos*y(i)
            vy(i,k) = u_div_cos*x(i)
            vz(i,k) = 0.0_8
        end do
    end do
end subroutine generate_Eldred_test_wind

subroutine generate_ts5_orography_field(this, f, npts, nlev, x, y, z)


    class(ts5_orography_generator_t), intent(in)  :: this
    integer(kind=4),                  intent(in)  :: npts, nlev
    real(kind=8),                     intent(in)  :: x(npts), y(npts), z(npts)
    real(kind=8),                     intent(out) :: f(npts,nlev)

    integer(kind=4) :: i, k
    real(kind=8)    :: dist

    do k = 1, nlev
        do i=1, npts
            dist = acos(x(i)*this%x_mount+y(i)*this%y_mount+z(i)*this%z_mount)
            f(i,k) = this%h_mount*max(0.0_8,1.0-dist/this%r_mount)
        end do
    end do

end subroutine generate_ts5_orography_field

subroutine generate_Ylm2_scalar_field(this, f, npts, nlev, x, y, z)

    use sph_coords_mod, only : cart2sph

    class(Ylm2_field_generator_t), intent(in)  :: this
    integer(kind=4),               intent(in)  :: npts, nlev
    real(kind=8),                  intent(in)  :: x(npts), y(npts), z(npts)
    real(kind=8),                  intent(out) :: f(npts,nlev)

    integer(kind=4) :: i, k
    real(kind=8)    :: lam, phi

    do k = 1, nlev
        do i=1, npts
            call cart2sph(x(i), y(i), z(i), lam, phi)
            if(mod(k,3) == 1) then
                f(i,k) = cos(2.0_8*lam)*cos(phi)**2
            else if(mod(k,3) == 2) then
                f(i,k) = cos(lam)*cos(phi)*sin(phi)
            else
                f(i,k) = 3.0_8*sin(phi)**2 -1.0_8
            end if
            ! f(i,k) = cos(phi)*cos(lam)
        end do
    end do

end subroutine generate_Ylm2_scalar_field

end module test_fields_mod

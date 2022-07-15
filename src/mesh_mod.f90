module mesh_mod
implicit none

type, public :: mesh_t
    type(tile_mesh_t), allocatable :: tile(:)
    real(kind=8) :: scale !projection scale factor
    real(kind=8) :: vertical_scale !vertical projection scale factor
    real(kind=8) :: omega !angular speed
    real(kind=8), allocatable :: rotation_axis(:) !angular velocity direction
    integer(kind=4) :: ts, te
contains
end type mesh_t

type, public :: tile_mesh_t

    integer(kind=4) :: is, ie, js, je, ks, ke
    integer(kind=4) :: nx, ny, nz                !global grid dimension
!    integer(kind=4) :: panel_ind
    integer(kind=4) :: halo_width

    real(kind=8), allocatable    :: rx(:,:,:), ry(:,:,:), rz(:,:,:) !cartesian coordinates of mesh points
    real(kind=8), allocatable    :: h(:,:,:)               !mean sea level height
    real(kind=8), allocatable    :: a1(:,:,:,:), a2(:,:,:,:)   !cartesian coordinates of covariant vecs at mesh points
    real(kind=8), allocatable    :: a3(:,:,:,:)
    real(kind=8), allocatable    :: b1(:,:,:,:), b2(:,:,:,:)   !cartesian coordinates of contravariant vecs at mesh points
    real(kind=8), allocatable    :: b3(:,:,:,:)
    real(kind=8), allocatable    :: Q(:,:,:,:)             !metric tensor at mesh-points
    real(kind=8), allocatable    :: QI(:,:,:,:)            !inverse metric tensor at mesh-points
    real(kind=8), allocatable    :: J(:,:,:)               !sqrt of metric tensor det at mesh-points
    real(kind=8), allocatable    :: G(:,:,:,:,:,:)         !Christofel symbols (i j,k,ix,iy,iz)
    real(kind=8)                 :: hx, hy                 !horizontal grid step
    real(kind=8)                 :: hz                     !vertical grid step
    real(kind=8)                 :: shift_i, shift_j       !determines shift of the first grid point from the boundary
    real(kind=8)                 :: shift_k                !determines shift of the first grid point from the boundary
    real(kind=8)                 :: alpha_0, beta_0        !determines coord start

    character(:), allocatable    :: points_type
contains

    procedure, public :: init => init_tile_mesh
    procedure, public :: get_alpha, get_beta, get_eta
end type tile_mesh_t

contains

subroutine init_tile_mesh(this, is, ie, js, je, ks, ke, halo_width)

    class(tile_mesh_t), intent(out) :: this
    integer(kind=4),    intent(in)  :: is, ie, js, je, ks, ke
    integer(kind=4),    intent(in)  :: halo_width

    integer(kind=4) :: isv, iev, jsv, jev

    isv = is - halo_width; iev = ie + halo_width
    jsv = js - halo_width; jev = je + halo_width

    allocate(this%rx(isv:iev,jsv:jev,ks:ke))
    allocate(this%ry(isv:iev,jsv:jev,ks:ke))
    allocate(this%rz(isv:iev,jsv:jev,ks:ke))
    allocate(this%h(isv:iev,jsv:jev,ks:ke))
    allocate(this%a1(4,isv:iev,jsv:jev,ks:ke)) !4-th component for shallow atmosphere approximation
    allocate(this%a2(4,isv:iev,jsv:jev,ks:ke))
    allocate(this%a3(4,isv:iev,jsv:jev,ks:ke))
    allocate(this%b1(3,isv:iev,jsv:jev,ks:ke))
    allocate(this%b2(3,isv:iev,jsv:jev,ks:ke))
    allocate(this%b3(4,isv:iev,jsv:jev,ks:ke))
    allocate(this%Q (6,isv:iev,jsv:jev,ks:ke)) !6 elements of 3x3 matrix are stored due to symmetricity
    allocate(this%QI(6,isv:iev,jsv:jev,ks:ke)) ! -'-'-
    allocate(this%J(isv:iev,jsv:jev,ks:ke))
    allocate(this%G(3,3,3,isv:iev,jsv:jev,ks:ke))

    this%is = is
    this%ie = ie
    this%js = js
    this%je = je
    this%ks = ks
    this%ke = ke

    this%halo_width = halo_width

end subroutine init_tile_mesh

pure function get_alpha(this, i) result(alpha)

    class(tile_mesh_t), intent(in) :: this
    integer(kind=4),    intent(in) :: i

    real(kind=8) :: alpha

    alpha = this%alpha_0 + (i-1+this%shift_i)*this%hx

end function get_alpha

pure function get_beta(this, j) result(beta)

    class(tile_mesh_t), intent(in) :: this
    integer(kind=4),    intent(in) :: j

    real(kind=8) :: beta

    beta = this%beta_0 + (j-1+this%shift_j)*this%hy

end function get_beta

pure function get_eta(this, k) result(eta)

    class(tile_mesh_t), intent(in) :: this
    integer(kind=4),    intent(in) :: k

    real(kind=8) :: eta

    eta = (k-1+this%shift_k)*this%hz

end function get_eta

end module mesh_mod

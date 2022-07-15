module latlon_regrid_mod

use abstract_regrid_mod,    only : regrid_t, regrid_vec_t
use domain_mod,             only : domain_t
use grid_field_mod,         only : grid_field_t
use halo_mod,               only : halo_t, halo_vec_t
use parcomm_mod,            only : parcomm_global

implicit none

type, public :: latlon_interp_t
    integer(kind=4)               :: Nlat, Nlon
    integer(kind=4)               :: ks, ke
    character(len=:), allocatable :: interp_type
    integer(kind=4)               :: stencil_width
    real(kind=8),     allocatable :: wx(:,:,:), wy(:,:,:)
    integer(kind=4),  allocatable :: tile_ind(:,:)
    integer(kind=4),  allocatable :: x_ind(:,:)
    integer(kind=4),  allocatable :: y_ind(:,:)

    contains
    procedure, public :: interpolate => latlon_interpolate
end type latlon_interp_t

type, public, extends(regrid_t) :: latlon_regrid_t
    integer(kind=4) :: Nlon, Nlat
    character(len=:), allocatable :: scalar_grid_type
    type(latlon_interp_t)         :: interp_scalar
    class(halo_t), allocatable    :: scalar_halo
    type(grid_field_t)            :: work_field

    contains
    procedure :: do_regrid     => do_latlon_scalar_regrid
end type latlon_regrid_t

type, public, extends(regrid_vec_t) :: latlon_regrid_vec_t
    integer(kind=4) :: Nlon, Nlat
    character(len=:),  allocatable :: vector_grid_type
    character(len=:),  allocatable :: components_type
    type(latlon_interp_t)         :: interp_components
    !vector transformation from curvilinear grid uv to latlon grid uv:
    real(kind=8),      allocatable :: u2u(:,:), u2v(:,:), v2u(:,:), v2v(:,:)
    class(halo_vec_t), allocatable :: halo
    type(grid_field_t)             :: work_u, work_v
    integer(kind=4)                :: ks, ke

    contains
    procedure          :: do_regrid => do_latlon_vector_regrid
    procedure, private :: transform_components => transform_components_to_latlon_uv
end type latlon_regrid_vec_t

contains

subroutine do_latlon_scalar_regrid(this,fout,f,domain)
    class(latlon_regrid_t),    intent(inout) :: this
    real(kind=8),              intent(inout) :: fout(:,:,:)
    type(grid_field_t),        intent(in)    :: f
    type(domain_t),            intent(in)    :: domain

    if(this%scalar_grid_type == "Ah") then
        call this%interp_scalar%interpolate(fout,f)
    else if(this%scalar_grid_type == "A") then
        call this%work_field%assign(f,domain%mesh_o)
        call this%scalar_halo%get_halo_scalar(this%work_field,&
                                     domain,this%interp_scalar%stencil_width/2)
        call this%interp_scalar%interpolate(fout,this%work_field)
    else
        call parcomm_global%abort(this%scalar_grid_type//&
                                  " scalar grid type is not supported in latlon regrid,"//&
                                  " use A or Ah")
    end if

end subroutine do_latlon_scalar_regrid

subroutine do_latlon_vector_regrid(this,uout,vout,u,v,domain)
    class(latlon_regrid_vec_t), intent(inout) :: this
    real(kind=8),               intent(inout) :: uout(:,:,:), vout(:,:,:)
    type(grid_field_t),         intent(in)    :: u, v
    type(domain_t),             intent(in)    :: domain

    integer(kind=4) :: i, j, k
    real(kind=8)    :: zu, zv

    if(this%vector_grid_type == "Ah") then
        call this%interp_components%interpolate(uout,u)
        call this%interp_components%interpolate(vout,v)
    else if(this%vector_grid_type == "A" .or. &
            this%vector_grid_type == "C") then

        if(this%vector_grid_type == "A") then
            call this%work_u%assign(u,domain%mesh_o)
            call this%work_v%assign(v,domain%mesh_o)
        else if(this%vector_grid_type == "C") then
            !WORKAROUND, interpolation of C-vector components to A-points
            call interpolate_C_vec_to_A_serial4(this%work_u, this%work_v, u, v, &
                                                                   domain%mesh_o)
        end if
        !WORKAROUND, need standard subroutine here
        if(this%components_type == "covariant") &
            call transform_cov_to_contra_non_staggered(this%work_u,this%work_v,domain%mesh_o)

        call this%halo%get_halo_vector(this%work_u, this%work_v, &
                                       domain,this%interp_components%stencil_width/2)
        call this%interp_components%interpolate(uout,this%work_u)
        call this%interp_components%interpolate(vout,this%work_v)
    else
        call parcomm_global%abort(this%vector_grid_type//&
                                  " scalar grid type is not supported in latlon regrid,"//&
                                  " use A or Ah")
    end if

    call this%transform_components(uout,vout)

end subroutine do_latlon_vector_regrid

subroutine transform_components_to_latlon_uv(this, uout,vout)

    class(latlon_regrid_vec_t), intent(in) :: this

    real(kind=8), intent(inout) :: uout(1:this%Nlon,1:this%Nlat,this%ks:this%ke),&
                                    vout(1:this%Nlon,1:this%Nlat,this%ks:this%ke)

    real(kind=8)    :: zu, zv
    integer(kind=4) :: i, j, k

    do k=this%ks,this%ke
        do j=1,this%Nlat
            do i=1,this%Nlon
                zu = uout(i,j,k)
                zv = vout(i,j,k)
                uout(i,j,k) = this%u2u(i,j)*zu+this%v2u(i,j)*zv
                vout(i,j,k) = this%u2v(i,j)*zu+this%v2v(i,j)*zv
            end do
        end do
    end do

end subroutine transform_components_to_latlon_uv

subroutine latlon_interpolate(this,fout,f)
    class(latlon_interp_t), intent(in)    :: this
    real(kind=8),           intent(inout) :: fout(this%Nlon,this%Nlat,this%ks:this%ke)
    type(grid_field_t),     intent(in)    :: f

    integer(kind=4) :: i, j, k, i1, j1
    real(kind=8)    :: p

    do k=this%ks, this%ke
        do j = 1, this%Nlat
            do i = 1, this%Nlon

                fout(i,j,k) = 0.0_8
                do j1=1, this%stencil_width
                    do i1 = 1, this%stencil_width
                        p = f%tile(this%tile_ind(i,j))%p(this%x_ind(i,j)+i1-1,this%y_ind(i,j)+j1-1,k)
                        fout(i,j,k) = fout(i,j,k)+this%wx(i1,i,j)*this%wy(j1,i,j)*p
                    end do
                end do

            end do
        end do
    end do
end subroutine latlon_interpolate

subroutine transform_cov_to_contra_non_staggered(u,v,mesh)
    use mesh_mod, only : mesh_t

    type(grid_field_t), intent(inout) :: u, v
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: i, j, k, t
    real(kind=8)    :: zu, zv

    do t=mesh%ts, mesh%te
        do k=mesh%tile(t)%ks,mesh%tile(t)%ke
            do j=mesh%tile(t)%js,mesh%tile(t)%je
                do i=mesh%tile(t)%is,mesh%tile(t)%ie
                    zu = u%tile(t)%p(i,j,k)
                    zv = v%tile(t)%p(i,j,k)
                    u%tile(t)%p(i,j,k) = mesh%tile(t)%Qi(1,i,j,k)*zu+ &
                                         mesh%tile(t)%Qi(2,i,j,k)*zv
                    v%tile(t)%p(i,j,k) = mesh%tile(t)%Qi(2,i,j,k)*zu+ &
                                         mesh%tile(t)%Qi(3,i,j,k)*zv
                end do
            end do
        end do
    end do
end subroutine transform_cov_to_contra_non_staggered

subroutine interpolate_C_vec_to_A_serial4(uout, vout, u, v, mesh_o)
    use mesh_mod, only : mesh_t

    type(grid_field_t), intent(inout) :: uout, vout
    type(grid_field_t), intent(in)    :: u, v
    type(mesh_t),       intent(in)    :: mesh_o

    integer(kind=8) i, j, k, t, n
    !Interpolates components of vector at C-staggered grid
    !to cell centers (A-grid) using 4th order interpolation
    !Serial case is assumed: we have all needed values (indices from 1 to n) in hand

    do t=mesh_o%ts,mesh_o%te
    do k=mesh_o%tile(t)%ks, mesh_o%tile(t)%ke
        do j=mesh_o%tile(t)%js, mesh_o%tile(t)%je
            do i = 1,mesh_o%tile(t)%nx
                uout%tile(t)%p(i,j,k) = 0.5_8*(u%tile(t)%p(i,j,k)+u%tile(t)%p(i+1,j,k))
            end do
        end do

        do j=1, mesh_o%tile(t)%ny
            do i = mesh_o%tile(t)%is,mesh_o%tile(t)%ie
                vout%tile(t)%p(i,j,k) = 0.5_8*(v%tile(t)%p(i,j,k)+v%tile(t)%p(i,j+1,k))
            end do
        end do

    end do
    end do

    do t=mesh_o%ts, mesh_o%te

    do k=mesh_o%tile(t)%ks, mesh_o%tile(t)%ke
        do j=mesh_o%tile(t)%js, mesh_o%tile(t)%je
            uout%tile(t)%p(1,j,k) = (5._8*u%tile(t)%p(1,j,k)+15._8*u%tile(t)%p(2,j,k)-&
                                     5._8*u%tile(t)%p(3,j,k)+u%tile(t)%p(4,j,k))/16._8
            do i = 2,mesh_o%tile(t)%nx-1
                uout%tile(t)%p(i,j,k) = (-u%tile(t)%p(i-1,j,k)+9._8*u%tile(t)%p(i,j,k)+&
                                      9._8*u%tile(t)%p(i+1,j,k)-u%tile(t)%p(i+2,j,k))/16._8
            end do
            n = mesh_o%tile(t)%nx
            uout%tile(t)%p(n,j,k) = (5._8*u%tile(t)%p(n+1,j,k)+15._8*u%tile(t)%p(n,j,k)-&
                                     5._8*u%tile(t)%p(n-1,j,k)+u%tile(t)%p(n-2,j,k))/16._8
        end do

        do i= mesh_o%tile(t)%is, mesh_o%tile(t)%ie
            vout%tile(t)%p(i,1,k) = (5._8*v%tile(t)%p(i,1,k)+15._8*v%tile(t)%p(i,2,k)-&
                                     5._8*v%tile(t)%p(i,3,k)+v%tile(t)%p(i,4,k))/16._8
        end do
        do j=2, mesh_o%tile(t)%ny-1
            do i = mesh_o%tile(t)%is,mesh_o%tile(t)%ie
                vout%tile(t)%p(i,j,k) = (-v%tile(t)%p(i,j-1,k)+9._8*v%tile(t)%p(i,j,k)+&
                                         9._8*v%tile(t)%p(i,j+1,k)-v%tile(t)%p(i,j+2,k))/16._8
            end do
        end do
        n = mesh_o%tile(t)%ny
        do i= mesh_o%tile(t)%is, mesh_o%tile(t)%ie
            vout%tile(t)%p(i,n,k) = (5._8*v%tile(t)%p(i,n+1,k)+15._8*v%tile(t)%p(i,n,k)-&
                                     5._8*v%tile(t)%p(i,n-1,k)+v%tile(t)%p(i,n-2,k))/16._8
        end do
    end do
    end do

end subroutine interpolate_C_vec_to_A_serial4

end module latlon_regrid_mod

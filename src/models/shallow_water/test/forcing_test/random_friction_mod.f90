module random_friction_mod

use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field
use mesh_mod,               only : mesh_t
use domain_mod,             only : domain_t

implicit none

type random_friction_t
    integer(kind=4)           :: l !spherical harmonic number
    real(kind=8)              :: tau !decorrelation time
    real(kind=8)              :: sigma_amp ! amplitude
    real(kind=8), allocatable :: sigma_lm(:) !friction coefficient expansion in spherical harmonics
    type(grid_field_t), allocatable :: ylm_u(:), ylm_v(:)
    contains
        procedure :: apply_update_forcing => random_friction_apply_update
end type

type random_scalar_t
    integer(kind=4)           :: l !spherical harmonic number
    real(kind=8)              :: tau !decorrelation time
    real(kind=8)              :: amp ! amplitude
    real(kind=8), allocatable :: plm(:) !field expansion in spherical harmonics
    type(grid_field_t), allocatable :: ylm(:)
    contains
        procedure :: apply_update_forcing => random_scalar_apply_update
end type

contains

subroutine initialize_random_friction(random_friction, l, tau, sigma_amp, domain)
    type(random_friction_t), intent(out) :: random_friction
    integer(kind=4),         intent(in)  :: l
    real(kind=8),            intent(in)  :: tau, sigma_amp
    type(domain_t),          intent(in)  :: domain

    integer(kind=4) :: m

    random_friction%l         = l
    random_friction%tau       = tau
    random_friction%sigma_amp = sigma_amp

    allocate(random_friction%sigma_lm(-l:l))
    random_friction%sigma_lm(-l:l) = 0.0_8
    allocate(random_friction%ylm_u(-l:l), random_friction%ylm_v(-l:l))
    do m=-l,l
        call create_grid_field(random_friction%ylm_u(m),0,0,domain%mesh_u)
        call create_grid_field(random_friction%ylm_v(m),0,0,domain%mesh_v)
    end do

    call generate_spherical_functions_l(random_friction%ylm_u,domain%mesh_u,l)
    call generate_spherical_functions_l(random_friction%ylm_v,domain%mesh_v,l)

end subroutine initialize_random_friction

subroutine initialize_random_scalar(random_scalar, l, tau, amp, domain)
    type(random_scalar_t),   intent(out) :: random_scalar
    integer(kind=4),         intent(in)  :: l
    real(kind=8),            intent(in)  :: tau, amp
    type(domain_t),          intent(in)  :: domain

    integer(kind=4) :: m

    random_scalar%l   = l
    random_scalar%tau = tau
    random_scalar%amp = amp

    allocate(random_scalar%plm(-l:l))
    call random_number(random_scalar%plm(-l:l))
    random_scalar%plm(-l:l) = random_scalar%plm(-l:l)-0.5_8
    allocate(random_scalar%ylm(-l:l))
    do m=-l,l
        call create_grid_field(random_scalar%ylm(m),0,0,domain%mesh_p)
    end do

    call generate_spherical_functions_l(random_scalar%ylm,domain%mesh_p,l)

end subroutine initialize_random_scalar

subroutine generate_spherical_functions_l(f,mesh,l)

    use const_mod, only : pi

    type(grid_field_t), intent(inout) :: f(-l:l)
    type(mesh_t),       intent(in)    :: mesh
    integer(kind=4),    intent(in)    :: l

    integer(kind=4) :: t, is, ie, js, je,ks,ke, nlat
    integer(kind=8) :: i,j,k,m,ip
    real(kind=8), allocatable :: aplg(:,:)
    real(kind=8) :: lam

    do t = mesh%ts, mesh%te
        is = mesh%tile(t)%is; ie = mesh%tile(t)%ie
        js = mesh%tile(t)%js; je = mesh%tile(t)%je
        ks = mesh%tile(t)%ks; ke = mesh%tile(t)%ke
        nlat = (ie-is+1)*(je-js+1)
        allocate(aplg(nlat,0:l))
        call generate_associated_legendre_l(aplg,l,nlat,mesh%tile(t)%rz(is:ie,js:je,1))

        m = 0
        do k=ks,ke
            do j=js, je
                do i=is, ie
                    ip = (j-js)*(ie-is+1)+i-is+1
                    f(m)%tile(t)%p(i,j,k) = aplg(ip,m)/sqrt(2.0_8*pi)
                end do
            end do
        end do
        do m=1,l
            do k=ks,ke
                do j=js, je
                    do i=is, ie
                        lam = atan2(mesh%tile(t)%ry(i,j,1),mesh%tile(t)%rx(i,j,1))
                        ip = (j-js)*(ie-is+1)+i-is+1
                        f( m)%tile(t)%p(i,j,k) = aplg(ip,m)*sin(m*lam) / sqrt(2.0_8*pi)
                        f(-m)%tile(t)%p(i,j,k) = aplg(ip,m)*cos(m*lam) / sqrt(2.0_8*pi)
                    end do
                end do
            end do
        end do
        deallocate(aplg)
    end do
end subroutine generate_spherical_functions_l

subroutine generate_associated_legendre_l(p,l,nlat,sinlat)
    integer(kind=4), intent(in)  :: l, nlat
    real(kind=8),    intent(in)  :: sinlat(nlat)
    real(kind=8),    intent(out) :: p(nlat,0:l)

    real(kind=8) :: aplg((l+1)*(l+2)/2,nlat), c1((l+1)*(l+2)/2), c2((l+1)*(l+2)/2)
    real(kind=8) :: y
    integer(kind=4) :: ip, m, k, j, ind

    ind(l,m) = l*(l+1)/2+m+1
    do m=0,l
        do k=m,l
            ip = ind(k,m)
            if(k==m) then !P_mm except P_00
                c1(ip) =-sqrt((2._8*m+1._8)/(2._8*max(m,1)))
            else if(k==m+1) then !P_(m+1)m
                c1(ip) = sqrt(2._8*m+3._8)
            else
                !normalization mutiplicators:
                c1(ip) = sqrt((2._8*k+1._8)*(2._8*k-1._8)/(1._8*(k+m)*(k-m)))
                c2(ip) =-sqrt((2._8*k+1._8)*(k-m-1._8)*(k+m-1._8)/((2._8*k-3._8)*(k+m)*(k-m)))
            end if
        end do
    end do
    do j=1,nlat
        y = sinlat(j)
        do m=0,l
            if(m==0) then ! P_m=0_m=0 also
                ip = ind(0,0)
                aplg(ip,j) = sqrt(0.5_8)
            else !Pmm, m>0
                ip = ind(m,m)
                aplg(ip,j) =c1(ip)*sqrt(1.0_8-y**2)*aplg(ind(m-1,m-1),j)
            end if
            if(m<l) then !Pm+1,m
                ip = ind(m+1,m)
                aplg(ip,j) = c1(ip)*y*aplg(ind(m,m),j)
            end if
            do k=m+2,l ! regular case Plm via Pl-1,m and Pl-2,m
                ip = ind(k,m)
                aplg(ip,j) = c1(ip)*y*aplg(ind(k-1,m),j)+c2(ip)*aplg(ind(k-2,m),j)
            end do
        end do
    end do

    do j=1,nlat
        do m=0,l
            p(j,m) = aplg(ind(l,m),j)
        end do
    end do

end subroutine generate_associated_legendre_l

subroutine random_friction_apply_update(this, u_tend, v_tend, u, v, domain, dt)
    class(random_friction_t), intent(inout) :: this
    type(grid_field_t),       intent(inout) :: u_tend, v_tend
    type(grid_field_t),       intent(inout) :: u, v
    type(domain_t),           intent(in)    :: domain
    real(kind=8),             intent(in)    :: dt

    real(kind=8)    :: sigma_lm(-this%l:this%l), sigma_min
    integer(kind=4) :: l, m

    l = this%l
    call random_number(sigma_lm(-l:l))
    sigma_lm(-l:l) = sigma_lm(-l:l)-0.5_8

    this%sigma_lm(-l:l) = (1.0_8-dt/this%tau)*this%sigma_lm(-l:l)+&
                                    dt/this%tau*sigma_lm(-l:l)
    this%sigma_lm(0) = sqrt(0.5_8)*this%sigma_lm(0)
    this%sigma_lm(-l:l) = this%sigma_lm(-l:l) / sqrt(sum(this%sigma_lm(-l:l)**2))

    call u_tend%assign(0.0_8, domain%mesh_u)
    call v_tend%assign(0.0_8, domain%mesh_v)

    do m=-l,l
        call u_tend%update(this%sigma_lm(m),this%ylm_u(m), domain%mesh_u)
        call v_tend%update(this%sigma_lm(m),this%ylm_v(m), domain%mesh_v)
    end do

    !non-negativity of friction coefficients:
    do m=-l,l
        sigma_min = min(u_tend%minimum(domain%mesh_u,domain%parcomm),&
                        v_tend%minimum(domain%mesh_v,domain%parcomm), 0.0_8)
        call u_tend%update(-sigma_min,domain%mesh_u)
        call v_tend%update(-sigma_min,domain%mesh_v)
    end do

    call u_tend%assign_prod(-0.5_8*this%sigma_amp,u_tend,u,domain%mesh_u)
    call v_tend%assign_prod(-0.5_8*this%sigma_amp,v_tend,v,domain%mesh_v)
end subroutine

subroutine random_scalar_apply_update(this, h, domain, dt)

    use mpi

    class(random_scalar_t),   intent(inout) :: this
    type(grid_field_t),       intent(inout) :: h
    type(domain_t),           intent(in)    :: domain
    real(kind=8),             intent(in)    :: dt

    real(kind=8)    :: plm(-this%l:this%l), a
    integer(kind=4) :: l, m, ierr

    l = this%l
    if(domain%parcomm%myid == 0) then
        call random_number(plm(-l:l))
        plm(-l:l) = plm(-l:l)-0.5_8

        !this%plm(-l:l) = (1.0_8-dt/this%tau)*this%plm(-l:l)+&
        !                                             dt/this%tau*plm(-l:l)
        !this%plm(0) = sqrt(0.5_8)*this%plm(0)
        !this%plm(-l:l) = this%plm(-l:l) / sqrt(sum(this%plm(-l:l)**2))
        a = 1.0_8-dt/this%tau
        this%plm(-l:l) = a*this%plm(-l:l)+sqrt(1.0_8-a**2)*plm(-l:l)
    end if
    call mpi_bcast(this%plm,2*l+1,MPI_DOUBLE,0,domain%parcomm%comm_w,ierr)

    call h%assign(sqrt(0.5_8)*this%plm(0), this%ylm(0), domain%mesh_p)
    do m=1,l
        call h%update(this%plm( m),this%ylm( m), domain%mesh_p)
        call h%update(this%plm(-m),this%ylm(-m), domain%mesh_p)
    end do
    call h%assign(this%amp/sqrt(2*l+1.0_8),h,domain%mesh_p)

end subroutine

end module random_friction_mod

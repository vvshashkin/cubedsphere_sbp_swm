module solid_rotation3d_therm_mod

use test_fields_3d_mod, only : scalar_field3d_t
use mesh_mod,           only : tile_mesh_t
use grid_field_mod,     only : tile_field_t
use const_mod,          only : Cp, rgaz, kappa

implicit none

private
public :: solid_rotation_theta_Nb_t, solid_rotation_PExner_Nb_t, &
          solid_rotation_theta_isoT_t, solid_rotation_PExner_isoT_t

type, extends(scalar_field3d_t) :: solid_rotation_theta_Nb_t
    real(kind=8) :: U0, omega, sphere_rad, Nb, grav
    real(kind=8) :: alpha = 0.0_8
    real(kind=8) :: Teq = 300._8
    real(kind=8) :: peq = 1e5_8
    real(kind=8) :: pref = 1e5_8
    contains
    procedure :: get_scalar_field_tile => get_theta_Nb_tile
end type solid_rotation_theta_Nb_t

type, extends(scalar_field3d_t) :: solid_rotation_PExner_Nb_t
    real(kind=8) :: U0, omega, sphere_rad, Nb, grav
    real(kind=8) :: alpha = 0.0_8
    real(kind=8) :: Teq = 300._8
    real(kind=8) :: peq = 1e5_8
    real(kind=8) :: pref = 1e5_8
    contains
    procedure :: get_scalar_field_tile => get_PExner_Nb_tile
end type solid_rotation_PExner_Nb_t

type, extends(scalar_field3d_t) :: solid_rotation_theta_isoT_t
    real(kind=8) :: U0, omega, sphere_rad, grav
    real(kind=8) :: alpha = 0.0_8
    real(kind=8) :: T0 = 300._8
    real(kind=8) :: p0 = 93000._8
    real(kind=8) :: pref = 1e5_8
    contains
    procedure :: get_scalar_field_tile => get_theta_isoT_tile
end type solid_rotation_theta_isoT_t

type, extends(scalar_field3d_t) :: solid_rotation_PExner_isoT_t
    real(kind=8) :: U0, omega, sphere_rad, grav
    real(kind=8) :: alpha = 0.0_8
    real(kind=8) :: T0 = 300._8
    real(kind=8) :: p0 = 93000._8
    real(kind=8) :: pref = 1e5_8
    contains
    procedure :: get_scalar_field_tile => get_PExner_isoT_tile
end type solid_rotation_PExner_isoT_t

contains

subroutine get_theta_Nb_tile(this,f,mesh,halo_width)
    class(solid_rotation_theta_Nb_t),    intent(in)    :: this
    type(tile_field_t),                  intent(inout) :: f
    type(tile_mesh_t),                   intent(in)    :: mesh
    integer(kind=4),                     intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: phi, T, P

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                phi = asin(mesh%rx(i,j,k)*sin(this%alpha)+mesh%rz(i,j,k)*cos(this%alpha))
                T = T_Nb(this%u0,this%omega,this%sphere_rad,this%Nb,&
                         this%grav,this%Teq,phi,mesh%h(i,j,k))
                P = PExner_Nb(this%u0,this%omega,this%sphere_rad,this%Nb,&
                              this%grav,this%Teq,this%peq,this%pref,phi,mesh%h(i,j,k))
                f%p(i,j,k) = T / P
            end do
        end do
    end do

end subroutine get_theta_Nb_tile

subroutine get_PExner_Nb_tile(this,f,mesh,halo_width)
    class(solid_rotation_PExner_Nb_t),   intent(in)    :: this
    type(tile_field_t),                  intent(inout) :: f
    type(tile_mesh_t),                   intent(in)    :: mesh
    integer(kind=4),                     intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: phi

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                phi = asin(mesh%rx(i,j,k)*sin(this%alpha)+mesh%rz(i,j,k)*cos(this%alpha))
                f%p(i,j,k) = PExner_Nb(this%u0,this%omega,this%sphere_rad,this%Nb,&
                                       this%grav,this%Teq,this%peq,this%pref,phi,mesh%h(i,j,k))
            end do
        end do
    end do

end subroutine get_PExner_Nb_tile

subroutine get_theta_isoT_tile(this,f,mesh,halo_width)
    class(solid_rotation_theta_isoT_t),    intent(in)    :: this
    type(tile_field_t),                    intent(inout) :: f
    type(tile_mesh_t),                     intent(in)    :: mesh
    integer(kind=4),                       intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: phi, ps

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                phi = asin(mesh%rx(i,j,k)*sin(this%alpha)+mesh%rz(i,j,k)*cos(this%alpha))
                ps = p_surf_isoT(this%u0,this%omega,this%sphere_rad,this%T0,&
                                 this%p0,this%grav,phi)
                f%p(i,j,k) =  this%T0*(this%pref/ps)**kappa*exp(this%grav*mesh%h(i,j,k) / (Cp*this%T0))
            end do
        end do
    end do

end subroutine get_theta_isoT_tile

subroutine get_PExner_isoT_tile(this,f,mesh,halo_width)
    class(solid_rotation_PExner_isoT_t),   intent(in)    :: this
    type(tile_field_t),                    intent(inout) :: f
    type(tile_mesh_t),                     intent(in)    :: mesh
    integer(kind=4),                       intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: phi, ps

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                phi = asin(mesh%rx(i,j,k)*sin(this%alpha)+mesh%rz(i,j,k)*cos(this%alpha))
                ps = p_surf_isoT(this%u0,this%omega,this%sphere_rad,this%T0,&
                                 this%p0,this%grav,phi)
                f%p(i,j,k) =  (ps/this%pref)**kappa*exp(-this%grav*mesh%h(i,j,k) / (Cp*this%T0))
            end do
        end do
    end do

end subroutine get_PExner_isoT_tile

real(kind=8) pure function T_surf_Nb(u0,omega,a,Nb,grav,Teq,phi) result(ts)
    real(kind=8), intent(in) ::  u0,omega,a,Nb,grav,phi,Teq

    real(kind=8) :: G

    G = grav**2 / (Nb**2*Cp)
    ts = G+(Teq-G)*exp(-0.25_8*u0*Nb**2/grav**2 * (u0+2.0_8*omega*a)*(cos(2.0_8*phi)-1.0_8))
end function

real(kind=8) pure function T_Nb(u0,omega,a,Nb,grav,Teq,phi,h) result(T)
    real(kind=8), intent(in) ::  u0,omega,a,Nb,grav,phi,Teq,h

    real(kind=8) :: G, exponent

    exponent = exp(Nb**2*h/grav)
    G = grav**2 / (Nb**2*Cp)
    T = G*(1._8-exponent)+T_surf_Nb(u0,omega,a,Nb,grav,Teq,phi)*exponent
end function

real(kind=8) pure function PExner_Nb(u0,omega,a,Nb,grav,Teq,peq,pref,phi,h) result(P)
    real(kind=8), intent(in) ::  u0,omega,a,Nb,grav,phi,Teq,peq,pref,h

    real(kind=8) :: G

    G = grav**2 / (Nb**2*Cp)
    P = (peq/pref)**kappa*(T_nb(u0,omega,a,Nb,grav,Teq,phi,h)/Teq)* &
        exp(0.25_8*u0/(G*Cp)*(u0+2.0_8*omega*a)*(cos(2.0_8*phi)-1.0_8) - grav*h/(G*Cp))
end function

real(kind=8) pure function p_surf_isoT(u0,omega,a,T0,p0,grav,phi) result(ps)
    real(kind=8), intent(in) ::  u0,omega,a,T0,p0,grav,phi

    real(kind=8) :: Nb2

    Nb2 = grav**2 / (Cp*T0)
    ps = p0*exp(-0.5_8*Nb2*u0/(grav**2*kappa)*(u0+2.0_8*omega*a)*(sin(phi)**2-1.0_8))
end function

end module solid_rotation3d_therm_mod

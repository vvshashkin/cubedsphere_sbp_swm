module test_grid_field_mod

implicit none

private
public :: test_grid_field

!real(kind=8), parameter :: a=0.1_8,b=1000.0_8,c=123.456_8

abstract interface
    real(kind=8) function fxyzh(x,y,z,h) result(f)
        real(kind=8), intent(in) ::  x,y,z,h
    end function fxyzh
end interface

contains

subroutine test_grid_field()
    use grid_field_mod,         only : grid_field_t
    use domain_mod,             only : domain_t
    use domain_factory_mod,     only : create_domain
    use grid_field_factory_mod, only : create_grid_field

    type(domain_t)     :: domain
    type(grid_field_t) :: f1, f2, f3, f4

    integer(kind=4)  :: nh=100, nz=10, halo_width=10
    character(len=1) :: hor_grid_type = 'C'

    integer(kind=4)  :: t, i, j, k
    real(kind=8)     :: err

    logical :: is_passed

    call create_domain(domain, "cube", hor_grid_type, nh, nz)

    call create_grid_field(f1, halo_width, 0, domain%mesh_p)
    f2 = f1%create_similar()
    f3 = f1%create_similar()

    call init_grid_field_fxyz(f1,domain%mesh_p,fx)
    call init_grid_field_fxyz(f2,domain%mesh_p,fy)
    call init_grid_field_fxyz(f3,domain%mesh_p,fz)

    f4 = f1%copy()
    call f4%update(-1.0_8, f1, domain%mesh_p)

    is_passed = (f4%algebraic_norm2(domain%mesh_p,domain%parcomm)==0.0_8)

    call f4%assign(0.0_8,domain%mesh_p)
    call f4%update(10.0_8,f2,domain%mesh_p)
    call f4%update(0.1_8,f3,domain%mesh_p)
    call f1%assign(-10.0_8,f2,-0.1_8,f3,domain%mesh_p)
    call f4%assign(1.0_8,f4,1.0_8,f1,domain%mesh_p)

    is_passed = is_passed .and. (f4%algebraic_norm2(domain%mesh_p,domain%parcomm)==0.0_8)

    call f4%assign(20.0_8,f2,0.1_8,f3,domain%mesh_p)
    call f4%update(2.0_8,f1,0.1_8,f3,domain%mesh_p)

    is_passed = is_passed .and. (f4%algebraic_norm2(domain%mesh_p,domain%parcomm)<1e-16_8*6*nh**2)

    is_passed = is_passed .and. &
               abs(1.0_8-f2%algebraic_norm2(domain%mesh_p,domain%parcomm)**2/f2%algebraic_dot(f2,domain%mesh_p,domain%parcomm))<5e-16_8

    call f4%assign(2.0_8,f2,1.5_8,f3,domain%mesh_p)
    call f1%assign(2.0_8,f2,0.5_8,f3,domain%mesh_p)
    !f4 = f1+f3
    !||f4||^2 = ||f1||^2+||f3||^2+(f1,f3):
    err = 1.0_8 - &
          (f1%algebraic_norm2(domain%mesh_p,domain%parcomm)**2+ &
           f3%algebraic_norm2(domain%mesh_p,domain%parcomm)**2+ &
           f3%algebraic_dot(f1,domain%mesh_p,domain%parcomm)+   &
           f1%algebraic_dot(f3,domain%mesh_p,domain%parcomm))/ &
           f4%algebraic_norm2(domain%mesh_p,domain%parcomm)**2

    is_passed = is_passed .and. abs(err) < 1e-14


    if(is_passed) then
        call domain%parcomm%print("grid_field test passed")
    else
        call domain%parcomm%print("grid_field test failed")
    end if

end subroutine test_grid_field

real(kind=8) function fx(x,y,z,h) result(f)
    real(kind=8), intent(in) ::  x,y,z,h
    f = (h+1.0_8)*x
end
real(kind=8) function fy(x,y,z,h) result(f)
    real(kind=8), intent(in) ::  x,y,z,h
    f = (h+1.0_8)*y
end
real(kind=8) function fz(x,y,z,h) result(f)
    real(kind=8), intent(in) ::  x,y,z,h
    f = (h+1.0_8)*z
end

subroutine init_grid_field_fxyz(gf, mesh, f)
    use grid_field_mod,         only : grid_field_t
    use mesh_mod,               only : mesh_t

    type(grid_field_t), intent(inout) :: gf
    type(mesh_t),       intent(in)    :: mesh
    procedure(fxyzh)                  :: f

    integer(kind=4) :: t,i,j,k

    do t=mesh%ts,mesh%te
        do k = mesh%tile(t)%ks, mesh%tile(t)%ke
            do j = mesh%tile(t)%js, mesh%tile(t)%je
                do i = mesh%tile(t)%is, mesh%tile(t)%ie
                    gf%tile(t)%p(i,j,k) = f(mesh%tile(t)%rx(i,j,k), mesh%tile(t)%ry(i,j,k), mesh%tile(t)%rz(i,j,k),1.0_8*k)
                end do
            end do
        end do
    end do
end subroutine init_grid_field_fxyz

end module test_grid_field_mod

module test_halo_mod

implicit none

private
public   :: test_halo

contains

subroutine test_halo()

    use mpi
    use domain_mod,             only : domain_t
    use domain_factory_mod,     only : create_domain
    use halo_mod,               only : halo_t, halo_vec_t
    use halo_factory_mod,       only : create_halo_procedure, create_vector_halo_procedure
    use grid_field_mod,         only : grid_field_t
    use grid_field_factory_mod, only : create_grid_field

    type(domain_t) domain
    class(halo_t),     allocatable :: halo
    class(halo_vec_t), allocatable :: halo_vec

    integer(kind=4), parameter :: nh = 4,nz=11,halo_width=1

    type(grid_field_t) :: f1, f2, u1, v1, u2, v2
    integer(kind=4) :: t, i, j, k
    integer(kind=8) :: ierr

    real(kind=8) err, gl_err

    call create_domain(domain, "cube", 'A', nh, nz)
    call domain%parcomm%print("cubed sphere default A grid scalar halo test")

    call create_halo_procedure(halo,domain,halo_width,"A_default")

    call create_grid_field(f1, halo_width, 0, domain%mesh_p)
    call create_grid_field(f2, halo_width, 0, domain%mesh_p)

    call init_field_xz(f1,domain,0,halo_width)
    call init_field_xz(f2,domain,halo_width,halo_width)

    call halo%get_halo_scalar(f1,domain,halo_width)

    err = compare_grid_functions(domain,halo_width,f1,f2)
    call mpi_allreduce(err, gl_err, 1, mpi_double, mpi_max, domain%parcomm%comm_w, ierr)

    if(gl_err<1e-16) then
        call domain%parcomm%print( "A-default halo test passed")
    else
        call domain%parcomm%print( "A-default halo test failed, error = ")
    end if

    call create_vector_halo_procedure(halo_vec,domain,halo_width,"A_vec_default")

    call create_grid_field(u1, halo_width, 0, domain%mesh_u)
    call create_grid_field(v1, halo_width, 0, domain%mesh_v)
    call create_grid_field(u2, halo_width, 0, domain%mesh_u)
    call create_grid_field(v2, halo_width, 0, domain%mesh_v)

    call init_vec_field(u1,v1,domain,0,halo_width)
    call init_vec_field(u2,v2,domain,halo_width,halo_width)

    call halo_vec%get_halo_vector(u1,v1,domain,halo_width)

    err = compare_grid_functions(domain,halo_width,u1,u2)+compare_grid_functions(domain,halo_width,v1,v2)
    call mpi_allreduce(err, gl_err, 1, mpi_double, mpi_max, domain%parcomm%comm_w, ierr)

    if(gl_err<1e-16) then
        call domain%parcomm%print( "A-default vector halo test passed")
    else
        call domain%parcomm%print( "A-default vector halo test failed, error = ")
    end if

end subroutine test_halo

subroutine init_field_xz(f,domain,halo_width,halo_widthf)
    use domain_mod,         only : domain_t
    use grid_field_mod,     only : grid_field_t

    type(grid_field_t), intent(inout) :: f
    type(domain_t),     intent(in)    :: domain
    integer(kind=4),    intent(in)    :: halo_width,halo_widthf


    integer(kind=4) t, i, j, k
    integer(kind=4) p
    integer(kind=4) n(3), ex(3), ey(3)
    integer(kind=4) nh,nz

    real(kind=8) xyz(3), nh2
    real(kind=8) q(6)

    nh=domain%partition%nh
    nh2 = 0.5_8*nh
    nz=domain%partition%nz

    do t = domain%mesh_p%ts, domain%mesh_p%te
        p = domain%partition%panel_map(t)
        n(1:3)  = domain%topology%n(1:3,p)
        ex(1:3) = domain%topology%ex(1:3,p)
        ey(1:3) = domain%topology%ey(1:3,p)
        do j = domain%mesh_p%tile(t)%js-halo_width, domain%mesh_p%tile(t)%je+halo_width
            do i = domain%mesh_p%tile(t)%is-halo_width, domain%mesh_p%tile(t)%ie+halo_width
                xyz(1:3) = nh2*n(1:3)+(i-nh2-0.5)*ex(1:3)+(j-nh2-0.5)*ey(1:3)
                q(1) = testfun(nh2,halo_widthf, xyz(1),xyz(2),xyz(3))
                q(2) = testfun(nh2,halo_widthf, xyz(1),xyz(3),xyz(2))
                q(3) = testfun(nh2,halo_widthf, xyz(2),xyz(1),xyz(3))
                q(4) = testfun(nh2,halo_widthf, xyz(2),xyz(3),xyz(1))
                q(5) = testfun(nh2,halo_widthf, xyz(3),xyz(1),xyz(2))
                q(6) = testfun(nh2,halo_widthf, xyz(3),xyz(2),xyz(1))
                do k=domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
                    f%tile(t)%p(i,j,k) = q(mod(k,6)+1)
                end do
            end do
        end do
    end do
contains
    real(kind=8) function testfun(n,h,x,y,z) result(f)
        real(kind=8),    intent(in) :: n
        integer(kind=4), intent(in) :: h
        real(kind=8),    intent(in) :: x,y,z

        if(y >-n+h .and. y <= n-h) then
            f = 10*(n+0.5+(abs(x)-abs(z)))+(y+n-0.5)
        else
            f = 0.0_8
        end if
    end
end subroutine init_field_xz

real(kind=8) function compare_grid_functions(domain,halo_width,f1,f2) result(err)
    use domain_mod,         only : domain_t
    use grid_field_mod,     only : grid_field_t

    type(domain_t),     intent(in)    :: domain
    integer(kind=4),    intent(in)    :: halo_width
    type(grid_field_t), intent(inout) :: f1,f2

    integer(kind=4) t,i,j,k,is,ie

    err = 0.0
    do t = domain%mesh_p%ts, domain%mesh_p%te
        do k = domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
            do j = domain%mesh_p%tile(t)%js-halo_width, domain%mesh_p%tile(t)%je+halo_width
                is = domain%mesh_p%tile(t)%is-halo_width
                ie = domain%mesh_p%tile(t)%ie+halo_width
                !exclude halo corners at panel corners
                if(j < 1 .or. j > domain%partition%nh) then
                    is = max(is,1)
                    ie = min(ie,domain%partition%nh)
                end if

                do i = is,ie
                    err = max(err,abs(f1%tile(t)%p(i,j,k)-f2%tile(t)%p(i,j,k)))
                end do
            end do
        end do
    end do

end function compare_grid_functions

subroutine init_vec_field(u,v,domain,halo_width,halo_widthf)
    use domain_mod,         only : domain_t
    use grid_field_mod,     only : grid_field_t

    type(grid_field_t), intent(inout) :: u,v
    type(domain_t),     intent(in)    :: domain
    integer(kind=4),    intent(in)    :: halo_width,halo_widthf


    integer(kind=4) t, i, j, k
    integer(kind=4) p
    integer(kind=4) n(3), ex(3), ey(3)
    integer(kind=4) nh,nz

    real(kind=8) xyz(3), nh2
    real(kind=8) q(3,6)

    nh=domain%partition%nh
    nh2 = 0.5_8*nh
    nz=domain%partition%nz

    do t = domain%mesh_p%ts, domain%mesh_p%te
        p = domain%partition%panel_map(t)
        n(1:3)  = domain%topology%n(1:3,p)
        ex(1:3) = domain%topology%ex(1:3,p)
        ey(1:3) = domain%topology%ey(1:3,p)
        do j = domain%mesh_p%tile(t)%js-halo_width, domain%mesh_p%tile(t)%je+halo_width
            do i = domain%mesh_p%tile(t)%is-halo_width, domain%mesh_p%tile(t)%ie+halo_width
                xyz(1:3) = nh2*n(1:3)+(i-nh2-0.5)*ex(1:3)+(j-nh2-0.5)*ey(1:3)
                q(1:3,1) = testfun(nh2,halo_widthf, xyz,'x')
                q(1:3,2) = testfun(nh2,halo_widthf, xyz,'y')
                q(1:3,3) = testfun(nh2,halo_widthf, xyz,'z')
                q(1:3,4) = testfun(nh2,halo_widthf, xyz,'xy')
                q(1:3,5) = testfun(nh2,halo_widthf, xyz,'xz')
                q(1:3,6) = testfun(nh2,halo_widthf, xyz,'yz')
                do k=domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
                    u%tile(t)%p(i,j,k) = sum(q(1:3,mod(k,6)+1)*ex(1:3))
                    v%tile(t)%p(i,j,k) = sum(q(1:3,mod(k,6)+1)*ey(1:3))
                end do
            end do
        end do
    end do
contains
    function testfun(n,h,xyz,axis) result(v)
        real(kind=8),    intent(in) :: n
        integer(kind=4), intent(in) :: h
        real(kind=8),    intent(in) :: xyz(3)
        character(len=*),intent(in) :: axis
        !output:
        real(kind=8) , target :: v(3)
        !local
        integer(kind=4) i
        real(kind=8), parameter :: margin = 1e-16
        real(kind=8) x,y,z
        real(kind=8) vx,vy,vz
        real(kind=8) magnitude

        x = xyz(1); y=xyz(2); z=xyz(3)
        v(1:3) = 0

        do i=1,len(axis)
            select case(axis(i:i))
            case('x')
                x = xyz(2)
                y = xyz(3)
                z = xyz(1)
            case('y')
                x =-xyz(1)
                y = xyz(3)
                z = xyz(2)
            case('z')
                x = xyz(1)
                y = xyz(2)
                z = xyz(3)
            end select

            vx = 0.0_8; vy = 0.0_8
            if(z >-n+h .and. z <= n-h) then
                magnitude = z+2._8*n
                if(abs(x-n) < margin .or. abs(x+n) < margin) then
                    vy = sign(magnitude,x)
                else if(abs(y-n) < margin .or. abs(y+n) < margin) then
                    vx = sign(magnitude,-y)
                end if
            end if

            select case(axis(i:i))
            case('x')
                v(2) = v(2)+vx
                v(3) = v(3)+vy
            case('y')
                v(1) = v(1)-vx
                v(3) = v(3)+vy
            case('z')
                v(1) = v(1)+vx
                v(2) = v(2)+vy
            end select
        end do
    end
end subroutine init_vec_field

end module test_halo_mod

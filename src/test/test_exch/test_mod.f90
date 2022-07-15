module test_mod

use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field, create_grid_field_global

implicit none

contains

subroutine test_A_halo_exchange()

    use mpi

    use domain_mod,             only : domain_t
    use domain_factory_mod,     only : create_domain
    use exchange_abstract_mod,  only : exchange_t
    use exchange_halo_mod,      only : exchange_2D_halo_t
    use exchange_factory_mod,   only : create_o_points_halo_exchange

    type(domain_t) :: domain

    class(exchange_t), allocatable     :: exch_halo
    type(grid_field_t)                 :: f, u, v
    integer(kind=4)                    :: nh=50, nz=5, halo_width=25
    integer(kind=4)                    :: ierr, code

    integer(kind=4) :: ts, te
    integer(kind=4) :: ind, t, i, j, k, err_sum, gl_err_sum

    integer(kind=4)  :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number
    integer(kind=4)  :: i_out, j_out, i_step, j_step, i1, j1, p1, pn, is, ie
    character(len=1) :: first_dim_index

    real(kind=8) :: u_remote, v_remote, f_remote

    logical :: is_succeed

#define __fun1(i,j,k,pn) ((pn-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k)
#define __fun2(i,j,k,pn) ((pn-1)*nh*nh*nz + nh*nh*(k-1) + nh*(j-1) + i)

    call create_domain(domain, "cube", 'A', nh, nz)

    call domain%parcomm%print('Running A-grid halo exchange test!')

    call create_grid_field(f, halo_width, 0, domain%mesh_p)


    ts = domain%partition%ts
    te = domain%partition%te

    do t = ts, te
        f%tile(t)%p = huge(1.0_8)
        pn = domain%partition%panel_map(t)
        do k = domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
            do j = domain%mesh_p%tile(t)%js, domain%mesh_p%tile(t)%je
                do i = domain%mesh_p%tile(t)%is, domain%mesh_p%tile(t)%ie
                    f%tile(t)%p(i,j,k) =  __fun1(i,j,k,pn)
                end do
            end do
        end do
    end do

    !Init exchange
    exch_halo = create_o_points_halo_exchange( &
                domain%partition, domain%parcomm, domain%topology,  halo_width, 'full')

    !Perform exchange
    call exch_halo%do(f, domain%parcomm)

    err_sum = 0

    do t = ts, te
        pn = domain%partition%panel_map(t)
        do k = domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
            do j = domain%mesh_p%tile(t)%js-halo_width, domain%mesh_p%tile(t)%je+halo_width
                is = domain%mesh_p%tile(t)%is-halo_width
                ie = domain%mesh_p%tile(t)%ie+halo_width
                if (j<1 .or. j>nh) then
                    is = max(is,1)
                    ie = min(ie,domain%partition%nh)
                end if
                do i = is, ie
                    call domain%topology%get_true_panel_coords(i, j,pn, nh, nh, i1, j1, p1)
                    f_remote = __fun1(i1,j1,k,p1)
                    err_sum = err_sum + abs(int(f%tile(t)%p(i,j,k) - f_remote))
                end do
            end do
        end do
    end do

    call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)

    if (gl_err_sum==0) then
        call domain%parcomm%print('Scalar exchange test passed!')
    else
        call domain%parcomm%abort('Scalar exchange test failed! Error! Abort!')
    end if

    call create_grid_field(u, halo_width, 0, domain%mesh_u)
    call create_grid_field(v, halo_width, 0, domain%mesh_v)


    do t = ts, te
        pn = domain%partition%panel_map(t)
        u%tile(t)%p = huge(1.0_8)
        do k = domain%mesh_u%tile(t)%ks, domain%mesh_u%tile(t)%ke
            do j = domain%mesh_u%tile(t)%js, domain%mesh_u%tile(t)%je
                do i = domain%mesh_u%tile(t)%is, domain%mesh_u%tile(t)%ie
                    u%tile(t)%p(i,j,k) =  __fun1(i,j,k,pn)
                end do
            end do
        end do
    end do
    do t = ts, te
        pn = domain%partition%panel_map(t)
        v%tile(t)%p = huge(1.0_8)
        do k = domain%mesh_v%tile(t)%ks, domain%mesh_v%tile(t)%ke
            do j = domain%mesh_v%tile(t)%js, domain%mesh_v%tile(t)%je
                do i = domain%mesh_v%tile(t)%is, domain%mesh_v%tile(t)%ie
                    v%tile(t)%p(i,j,k) =  __fun2(i,j,k,pn)
                end do
            end do
        end do
    end do

    call exch_halo%do_vec(u, v, domain%parcomm)

    err_sum = 0

!In case of A grid mesh_p = mesh_u = mesh_v
    do t = ts, te
        pn = domain%partition%panel_map(t)
        do k = domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
            do j = domain%mesh_p%tile(t)%js-halo_width, domain%mesh_p%tile(t)%je+halo_width
                is = domain%mesh_p%tile(t)%is-halo_width
                ie = domain%mesh_p%tile(t)%ie+halo_width
                if (j<1 .or. j>nh) then
                    is = max(is,1)
                    ie = min(ie,domain%partition%nh)
                end if
                do i = is, ie
                    call domain%topology%get_true_panel_coords(i, j, pn, nh, nh, i1, j1, p1)
                    call domain%topology%find_basis_orientation(pn, p1, i_step, j_step, first_dim_index)
                    if (first_dim_index == "i") then
                        u_remote = i_step*__fun1(i1,j1,k,p1)
                        v_remote = j_step*__fun2(i1,j1,k,p1)
                    else
                        u_remote = j_step*__fun2(i1,j1,k,p1)
                        v_remote = i_step*__fun1(i1,j1,k,p1)
                    end if
                    err_sum = err_sum + abs(int(u%tile(t)%p(i,j,k) - u_remote)) + abs(int(v%tile(t)%p(i,j,k) - v_remote))
                end do
            end do
        end do
    end do

    call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, domain%parcomm%comm_w, ierr)

    if (gl_err_sum==0) then
        call domain%parcomm%print('Vector exchange test passed!')
    else
        call domain%parcomm%abort('Vector exchange test not passed! Error! Abort!')
    end if

#undef __fun1
#undef __fun2

end subroutine test_A_halo_exchange

subroutine test_Ah_halo_exchange()

    use mpi

    use domain_mod,             only : domain_t
    use domain_factory_mod,     only : create_domain
    use exchange_abstract_mod,  only : exchange_t
    use exchange_halo_mod,      only : exchange_2D_halo_t
    use exchange_factory_mod,   only : create_xy_points_halo_exchange
    use test_fields_mod,        only : set_vector_test_field, set_scalar_test_field, &
                                       xyz_f => xyz_scalar_field_generator

    type(domain_t) :: domain

    class(exchange_t), allocatable     :: exch_halo
    type(grid_field_t)                 :: f, u, v
    integer(kind=4)                    :: nh=50, nz=5, halo_width=25
    integer(kind=4)                    :: ierr, code

    call create_domain(domain, "cube", 'Ah', nh, nz)
    call domain%parcomm%print('Running Ah-grid halo exchange test!')
    exch_halo = create_xy_points_halo_exchange( &
                    domain%partition, domain%parcomm, domain%topology,  halo_width, 'full')

    call create_grid_field(f, halo_width, 0, domain%mesh_xy)
    call set_scalar_test_field(f,xyz_f,domain%mesh_xy,0)

    call exch_halo%do(f, domain%parcomm)

    print *, f%tile(1)%p(101,1,1)
    print *, f%tile(1)%p(101,0,1)
    print *, f%tile(1)%p(102,1,1)

end subroutine test_Ah_halo_exchange

subroutine test_halo_vec_C_exchange()

    use mpi

    use domain_mod,             only : domain_t
    use domain_factory_mod,     only : create_domain
    use exchange_abstract_mod,  only : exchange_t
    use exchange_halo_C_mod,    only : exchange_2D_halo_C_t
    use exchange_factory_mod,   only : create_symmetric_halo_vec_exchange_C

    type(domain_t) :: domain

    class(exchange_t), allocatable :: exch_halo
    type(grid_field_t)             :: f, u, v
    integer(kind=4)                :: nh=50, nz=5, halo_width=20
    integer(kind=4)                :: ierr

    integer(kind=4) :: ts, te
    integer(kind=4) :: ind, t, i, j, k, err_sum, gl_err_sum

    integer(kind=4)  :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number
    integer(kind=4)  :: i_out, j_out, i_step, j_step, pn, p1, i1, j1, l_halo, r_halo, is, ie
    character(len=1) :: first_dim_index

    real(kind=8) :: u_remote, v_remote

#define __fun1(i,j,k,pn) ((pn-1)*nh*(nh+1)*nz + nz*(nh+1)*(j-1) + nz*(i-1) + k)
#define __fun2(i,j,k,pn) ((pn-1)*nh*(nh+1)*nz + nh*(nh+1)*(k-1) + (nh+1)*(j-1) + i)

    call create_domain(domain, "cube", 'C', nh, nz)

    call domain%parcomm%print('Running halo_vec_C_exchange test!')

    call create_grid_field(u, halo_width+1, 0, domain%mesh_u)
    call create_grid_field(v, halo_width+1, 0, domain%mesh_v)

    ts = domain%partition%ts
    te = domain%partition%te

    do t = ts, te
        pn = domain%partition%panel_map(t)
        u%tile(t)%p = huge(1.0_8)
        do k = domain%mesh_u%tile(t)%ks, domain%mesh_u%tile(t)%ke
            do j = domain%mesh_u%tile(t)%js, domain%mesh_u%tile(t)%je
                do i = domain%mesh_u%tile(t)%is, domain%mesh_u%tile(t)%ie
                    u%tile(t)%p(i,j,k) =  __fun1(i,j,k,pn)
                end do
            end do
        end do
        v%tile(t)%p = huge(1.0_8)
        do k = domain%mesh_v%tile(t)%ks, domain%mesh_v%tile(t)%ke
            do j = domain%mesh_v%tile(t)%js, domain%mesh_v%tile(t)%je
                do i = domain%mesh_v%tile(t)%is, domain%mesh_v%tile(t)%ie
                    v%tile(t)%p(i,j,k) =  __fun2(i,j,k,pn)
                end do
            end do
        end do
    end do

    !Init exchange

    exch_halo = create_symmetric_halo_vec_exchange_C( &
                domain%partition, domain%parcomm, domain%topology, halo_width, 'full')

    !Perform exchange
    call exch_halo%do_vec(u, v, domain%parcomm)

    err_sum = 0

    do t = ts, te
        pn = domain%partition%panel_map(t)
        l_halo = halo_width
        r_halo = halo_width
        if (domain%mesh_u%tile(t)%is ==    1) l_halo = halo_width+1
        if (domain%mesh_u%tile(t)%ie == nh+1) r_halo = halo_width+1
        do k = domain%mesh_u%tile(t)%ks, domain%mesh_u%tile(t)%ke
            do j = domain%mesh_u%tile(t)%js-halo_width, domain%mesh_u%tile(t)%je+halo_width
                is = domain%mesh_p%tile(t)%is-l_halo
                ie = domain%mesh_p%tile(t)%ie+r_halo
                if (j<1 .or. j>nh) then
                    is = max(is,1)
                    ie = min(ie,nh+1)
                end if
                do i = is, ie
                    call domain%topology%get_true_panel_coords(i, j, pn, nh+1, nh, i1, j1, p1)
                    call domain%topology%find_basis_orientation(pn, p1, i_step, j_step, first_dim_index)
                    if (first_dim_index == "i") then
                        u_remote = i_step*__fun1(i1,j1,k,p1)
                    else
                        u_remote = j_step*__fun2(i1,j1,k,p1)
                    end if
                    err_sum = err_sum + abs(int(u%tile(t)%p(i,j,k) - u_remote))
                end do
            end do
        end do
    end do

    do t = ts, te
        pn = domain%partition%panel_map(t)
        l_halo = halo_width
        r_halo = halo_width
        if (domain%mesh_v%tile(t)%js ==    1) l_halo = halo_width+1
        if (domain%mesh_v%tile(t)%je == nh+1) r_halo = halo_width+1
        do k = domain%mesh_v%tile(t)%ks, domain%mesh_v%tile(t)%ke
            do j = domain%mesh_v%tile(t)%js-l_halo, domain%mesh_v%tile(t)%je+r_halo
                is = domain%mesh_p%tile(t)%is-halo_width
                ie = domain%mesh_p%tile(t)%ie+halo_width
                if (j<1 .or. j>nh+1) then
                    is = max(is,1)
                    ie = min(ie,nh)
                end if
                do i = is, ie
                    call domain%topology%get_true_panel_coords(i, j, pn, nh, nh+1, i1, j1, p1)
                    call domain%topology%find_basis_orientation(pn, p1, i_step, j_step, first_dim_index)
                    if (first_dim_index == "i") then
                        v_remote = j_step*__fun2(i1,j1,k,p1)
                    else
                        v_remote = i_step*__fun1(i1,j1,k,p1)
                    err_sum = err_sum + abs(int(v%tile(t)%p(i,j,k) - v_remote))
                    end if
                end do
            end do
        end do
    end do

    call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, domain%parcomm%comm_w, ierr)

    if (gl_err_sum==0) then
        call domain%parcomm%print('Test passed!')
    else
        call domain%parcomm%abort('Test not passed! Error! Abort!')
    end if

#undef __fun1
#undef __fun2

end subroutine test_halo_vec_C_exchange
subroutine test_halo_vec_Ch_exchange()

    use mpi

    use domain_mod,            only : domain_t
    use domain_factory_mod,    only : create_domain
    use exchange_abstract_mod, only : exchange_t
    use exchange_halo_Ch_mod,  only : exchange_2D_halo_Ch_t
    use exchange_factory_mod,  only : create_symmetric_halo_vec_exchange_Ch
    use mesh_mod,              only : mesh_t

    type(domain_t), target :: domain

    class(exchange_t), allocatable :: exch_halo
    type(grid_field_t)             :: f, u, v
    integer(kind=4)                :: nh=50, nz=5, halo_width=20
    integer(kind=4)                :: ierr

    integer(kind=4) :: ts, te
    integer(kind=4) :: ind, t, i, j, k, err_sum, gl_err_sum

    integer(kind=4)  :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number
    integer(kind=4)  :: i_out, j_out, i_step, j_step, pn, p1, i1, j1, l_halo, r_halo, is, ie
    character(len=1) :: first_dim_index

    real(kind=8) :: u_remote, v_remote

    type(mesh_t), pointer :: mesh_u, mesh_v

#define __fun1(i,j,k,pn) ((pn-1)*nh*(nh+1)*nz + nz*(nh+1)*(j-1) + nz*(i-1) + k)
#define __fun2(i,j,k,pn) ((pn-1)*nh*(nh+1)*nz + nh*(nh+1)*(k-1) + (nh+1)*(j-1) + i)

    call create_domain(domain, "cube", 'C', nh, nz)

    mesh_u => domain%mesh_y
    mesh_v => domain%mesh_x

    call domain%parcomm%print('Running halo_vec_Ch_exchange test!')

    call create_grid_field(u, halo_width+1, 0, mesh_u)
    call create_grid_field(v, halo_width+1, 0, mesh_v)

    ts = domain%partition%ts
    te = domain%partition%te

    do t = ts, te
        pn = domain%partition%panel_map(t)
        u%tile(t)%p = huge(1.0_8)
        do k = mesh_u%tile(t)%ks, mesh_u%tile(t)%ke
            do j = mesh_u%tile(t)%js, mesh_u%tile(t)%je
                do i = mesh_u%tile(t)%is, mesh_u%tile(t)%ie
                    u%tile(t)%p(i,j,k) =  __fun1(i,j,k,pn)
                end do
            end do
        end do
        v%tile(t)%p = huge(1.0_8)
        do k = mesh_v%tile(t)%ks, mesh_v%tile(t)%ke
            do j = mesh_v%tile(t)%js, mesh_v%tile(t)%je
                do i = mesh_v%tile(t)%is, mesh_v%tile(t)%ie
                    v%tile(t)%p(i,j,k) =  __fun2(i,j,k,pn)
                end do
            end do
        end do
    end do

    !Init exchange

    exch_halo = create_symmetric_halo_vec_exchange_Ch( &
                domain%partition, domain%parcomm, domain%topology, halo_width, 'full')

    !Perform exchange
    call exch_halo%do_vec(u, v, domain%parcomm)

    err_sum = 0

    do t = ts, te
        pn = domain%partition%panel_map(t)
        l_halo = halo_width
        r_halo = halo_width
        if (mesh_v%tile(t)%is ==    1) l_halo = halo_width+1
        if (mesh_v%tile(t)%ie == nh+1) r_halo = halo_width+1
        do k = mesh_v%tile(t)%ks, mesh_v%tile(t)%ke
            do j = mesh_v%tile(t)%js-halo_width, mesh_v%tile(t)%je+halo_width
                is = domain%mesh_p%tile(t)%is-l_halo
                ie = domain%mesh_p%tile(t)%ie+r_halo
                if (j<1 .or. j>nh) then
                    is = max(is,1)
                    ie = min(ie,nh+1)
                end if
                do i = is, ie
                    call domain%topology%get_true_panel_coords(i, j, pn, nh+1, nh, i1, j1, p1)
                    call domain%topology%find_basis_orientation(pn, p1, i_step, j_step, first_dim_index)
                    if (first_dim_index == "i") then
                        v_remote = j_step*__fun2(i1,j1,k,p1)
                    else
                        v_remote = i_step*__fun1(i1,j1,k,p1)
                    end if
                    err_sum = err_sum + abs(int(v%tile(t)%p(i,j,k) - v_remote))
                end do
            end do
        end do
    end do

    do t = ts, te
        pn = domain%partition%panel_map(t)
        l_halo = halo_width
        r_halo = halo_width
        if (mesh_u%tile(t)%js ==    1) l_halo = halo_width+1
        if (mesh_u%tile(t)%je == nh+1) r_halo = halo_width+1
        do k = mesh_u%tile(t)%ks, mesh_u%tile(t)%ke
            do j = mesh_u%tile(t)%js-l_halo, mesh_u%tile(t)%je+r_halo
                is = domain%mesh_p%tile(t)%is-halo_width
                ie = domain%mesh_p%tile(t)%ie+halo_width
                if (j<1 .or. j>nh+1) then
                    is = max(is,1)
                    ie = min(ie,nh)
                end if
                do i = is, ie
                    call domain%topology%get_true_panel_coords(i, j, pn, nh, nh+1, i1, j1, p1)
                    call domain%topology%find_basis_orientation(pn, p1, i_step, j_step, first_dim_index)
                    if (first_dim_index == "i") then
                        u_remote = i_step*__fun1(i1,j1,k,p1)
                    else
                        u_remote = j_step*__fun2(i1,j1,k,p1)
                    err_sum = err_sum + abs(int(u%tile(t)%p(i,j,k) - u_remote))
                    end if
                end do
            end do
        end do
    end do

    call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, domain%parcomm%comm_w, ierr)

    if (gl_err_sum==0) then
        call domain%parcomm%print('Test passed!')
    else
        call domain%parcomm%abort('Test not passed! Error! Abort!')
    end if

#undef __fun1
#undef __fun2

end subroutine test_halo_vec_Ch_exchange
subroutine test_halo_u_exchange()

    use mpi

    use domain_mod,             only : domain_t
    use domain_factory_mod,     only : create_domain
    use exchange_abstract_mod,  only : exchange_t
    use exchange_halo_C_mod,    only : exchange_2D_halo_C_t
    use exchange_factory_mod,   only : create_symm_halo_vec_exchange_U_points

    type(domain_t) :: domain

    class(exchange_t), allocatable :: exch_halo
    type(grid_field_t)             :: f, u, v
    integer(kind=4)                :: nh=50, nz=5, halo_width=25
    integer(kind=4)                :: ierr

    integer(kind=4) :: ts, te
    integer(kind=4) :: ind, t, i, j, k, err_sum, gl_err_sum

    integer(kind=4)  :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number
    integer(kind=4)  :: i_out, j_out, i_step, j_step, pn, p1, i1, j1, l_halo, r_halo, is, ie
    character(len=1) :: first_dim_index

    real(kind=8) :: u_remote, v_remote

#define __fun1(i,j,k,pn) ((pn-1)*nh*(nh+1)*nz + nz*(nh+1)*(j-1) + nz*(i-1) + k)
#define __fun2(i,j,k,pn) ((pn-1)*nh*(nh+1)*nz + nh*(nh+1)*(k-1) + (nh+1)*(j-1) + i)

    call create_domain(domain, "cube", 'C', nh, nz)

    call domain%parcomm%print('Running u_halo_C_exchange test!')

    call create_grid_field(u, halo_width+1, 0, domain%mesh_u)
    call create_grid_field(v, halo_width+1, 0, domain%mesh_v)

    ts = domain%partition%ts
    te = domain%partition%te

    do t = ts, te
        pn = domain%partition%panel_map(t)
        u%tile(t)%p = huge(1.0_8)
        do k = domain%mesh_u%tile(t)%ks, domain%mesh_u%tile(t)%ke
            do j = domain%mesh_u%tile(t)%js, domain%mesh_u%tile(t)%je
                do i = domain%mesh_u%tile(t)%is, domain%mesh_u%tile(t)%ie
                    u%tile(t)%p(i,j,k) =  __fun1(i,j,k,pn)
                end do
            end do
        end do
        v%tile(t)%p = huge(1.0_8)
        do k = domain%mesh_v%tile(t)%ks, domain%mesh_v%tile(t)%ke
            do j = domain%mesh_v%tile(t)%js, domain%mesh_v%tile(t)%je
                do i = domain%mesh_v%tile(t)%is, domain%mesh_v%tile(t)%ie
                    v%tile(t)%p(i,j,k) =  __fun2(i,j,k,pn)
                end do
            end do
        end do
    end do

    !Init exchange

    exch_halo = create_symm_halo_vec_exchange_U_points( &
                domain%partition, domain%parcomm, domain%topology, halo_width, 'full')

    !Perform exchange
    call exch_halo%do_vec(u, v, domain%parcomm)

    err_sum = 0

    do t = ts, te
        pn = domain%partition%panel_map(t)
        l_halo = halo_width
        r_halo = halo_width
        if (domain%mesh_u%tile(t)%is ==    1) l_halo = halo_width+1
        if (domain%mesh_u%tile(t)%ie == nh+1) r_halo = halo_width+1
        do k = domain%mesh_u%tile(t)%ks, domain%mesh_u%tile(t)%ke
            do j = domain%mesh_u%tile(t)%js-halo_width, domain%mesh_u%tile(t)%je+halo_width
                is = domain%mesh_p%tile(t)%is-l_halo
                ie = domain%mesh_p%tile(t)%ie+r_halo
                if (j<1 .or. j>nh) then
                    is = max(is,1)
                    ie = min(ie,nh+1)
                end if
                do i = is, ie
                    call domain%topology%get_true_panel_coords(i, j, pn, nh+1, nh, i1, j1, p1)
                    call domain%topology%find_basis_orientation(pn, p1, i_step, j_step, first_dim_index)
                    if (first_dim_index == "i") then
                        u_remote = i_step*__fun1(i1,j1,k,p1)
                    else
                        u_remote = j_step*__fun2(i1,j1,k,p1)
                    end if
                    err_sum = err_sum + abs(int(u%tile(t)%p(i,j,k) - u_remote))
                end do
            end do
        end do
    end do

    call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, domain%parcomm%comm_w, ierr)

    if (gl_err_sum==0) then
        call domain%parcomm%print('Test passed!')
    else
        call domain%parcomm%abort('Test not passed! Error! Abort!')
    end if

#undef __fun1
#undef __fun2

end subroutine test_halo_u_exchange
subroutine test_halo_v_exchange()

    use mpi

    use domain_mod,             only : domain_t
    use domain_factory_mod,     only : create_domain
    use exchange_abstract_mod,  only : exchange_t
    use exchange_halo_C_mod,    only : exchange_2D_halo_C_t
    use exchange_factory_mod,   only : create_symm_halo_vec_exchange_V_points

    type(domain_t) :: domain

    class(exchange_t), allocatable :: exch_halo
    type(grid_field_t)             :: f, u, v
    integer(kind=4)                :: nh=50, nz=5, halo_width=25
    integer(kind=4)                :: ierr

    integer(kind=4) :: ts, te
    integer(kind=4) :: ind, t, i, j, k, err_sum, gl_err_sum

    integer(kind=4)  :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number
    integer(kind=4)  :: i_out, j_out, i_step, j_step, pn, p1, i1, j1, l_halo, r_halo, is, ie
    character(len=1) :: first_dim_index

    real(kind=8) :: u_remote, v_remote

#define __fun1(i,j,k,pn) ((pn-1)*nh*(nh+1)*nz + nz*(nh+1)*(j-1) + nz*(i-1) + k)
#define __fun2(i,j,k,pn) ((pn-1)*nh*(nh+1)*nz + nh*(nh+1)*(k-1) + (nh+1)*(j-1) + i)

    call create_domain(domain, "cube", 'C', nh, nz)

    call domain%parcomm%print('Running v_halo_C_exchange test!')

    call create_grid_field(u, halo_width+1, 0, domain%mesh_u)
    call create_grid_field(v, halo_width+1, 0, domain%mesh_v)

    ts = domain%partition%ts
    te = domain%partition%te

    do t = ts, te
        pn = domain%partition%panel_map(t)
        u%tile(t)%p = huge(1.0_8)
        do k = domain%mesh_u%tile(t)%ks, domain%mesh_u%tile(t)%ke
            do j = domain%mesh_u%tile(t)%js, domain%mesh_u%tile(t)%je
                do i = domain%mesh_u%tile(t)%is, domain%mesh_u%tile(t)%ie
                    u%tile(t)%p(i,j,k) =  __fun1(i,j,k,pn)
                end do
            end do
        end do
        v%tile(t)%p = huge(1.0_8)
        do k = domain%mesh_v%tile(t)%ks, domain%mesh_v%tile(t)%ke
            do j = domain%mesh_v%tile(t)%js, domain%mesh_v%tile(t)%je
                do i = domain%mesh_v%tile(t)%is, domain%mesh_v%tile(t)%ie
                    v%tile(t)%p(i,j,k) =  __fun2(i,j,k,pn)
                end do
            end do
        end do
    end do

    !Init exchange

    exch_halo = create_symm_halo_vec_exchange_V_points( &
                domain%partition, domain%parcomm, domain%topology, halo_width, 'full')

    !Perform exchange
    call exch_halo%do_vec(u, v, domain%parcomm)

    err_sum = 0

    do t = ts, te
        pn = domain%partition%panel_map(t)
        l_halo = halo_width
        r_halo = halo_width
        if (domain%mesh_v%tile(t)%js ==    1) l_halo = halo_width+1
        if (domain%mesh_v%tile(t)%je == nh+1) r_halo = halo_width+1
        do k = domain%mesh_v%tile(t)%ks, domain%mesh_v%tile(t)%ke
            do j = domain%mesh_v%tile(t)%js-l_halo, domain%mesh_v%tile(t)%je+r_halo
                is = domain%mesh_p%tile(t)%is-halo_width
                ie = domain%mesh_p%tile(t)%ie+halo_width
                if (j<1 .or. j>nh+1) then
                    is = max(is,1)
                    ie = min(ie,nh)
                end if
                do i = is, ie
                    call domain%topology%get_true_panel_coords(i, j, pn, nh, nh+1, i1, j1, p1)
                    call domain%topology%find_basis_orientation(pn, p1, i_step, j_step, first_dim_index)
                    if (first_dim_index == "i") then
                        v_remote = j_step*__fun2(i1,j1,k,p1)
                    else
                        v_remote = i_step*__fun1(i1,j1,k,p1)
                    err_sum = err_sum + abs(int(v%tile(t)%p(i,j,k) - v_remote))
                    end if
                end do
            end do
        end do
    end do

    call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, domain%parcomm%comm_w, ierr)

    if (gl_err_sum==0) then
        call domain%parcomm%print('Test passed!')
    else
        call domain%parcomm%abort('Test not passed! Error! Abort!')
    end if

end subroutine test_halo_v_exchange

subroutine test_gather_exchange()

    use domain_mod,            only : domain_t
    use domain_factory_mod,    only : create_domain
    use exchange_abstract_mod, only : exchange_t
    use exchange_factory_mod,  only : create_gather_exchange

    type(domain_t) :: domain

    class(exchange_t), allocatable :: exch_gather
    type(grid_field_t)             :: f

    integer(kind=4) :: nh=50, nz=5, halo_width=10
    integer(kind=4) :: ierr, code
    integer(kind=4) :: master_id = 0

    integer(kind=4) :: ts, te
    integer(kind=4) :: i, j, k, t, N_tiles, pn
    integer(kind=4) :: is, ie, js, je, ks, ke

    real(kind=8) :: err_sum

    integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number

#define __fun(i,j,k,pn) ((pn-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k)

    call create_domain(domain, "cube", 'A', nh, nz)

    call domain%parcomm%print('Running gather exchange test!')

    !find start and end index of tiles belonging to the current proccesor
    ts = domain%partition%ts
    te = domain%partition%te

    !Init arrays

    N_tiles = domain%partition%num_tiles*domain%partition%num_panels

    if (domain%parcomm%myid == master_id) then
        call create_grid_field_global(f, 0, 0, domain%partition%tiles_p)
        do t = 1, domain%partition%Nt
            f%tile(t)%p = huge(1.0_8)
        end do
    else
        call create_grid_field(f, halo_width, 0, domain%mesh_p)
    end if

    do t = ts, te
        f%tile(t)%p = huge(1.0_8)
        pn = domain%partition%panel_map(t)
        do k = domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
            do j = domain%mesh_p%tile(t)%js, domain%mesh_p%tile(t)%je
                do i = domain%mesh_p%tile(t)%is, domain%mesh_p%tile(t)%ie
                    f%tile(t)%p(i,j,k) =  __fun(i,j,k,pn)
                end do
            end do
        end do
    end do

    !Init exchange

    call create_gather_exchange(exch_gather, 'p', domain%parcomm, domain%partition, master_id)

    !Perform exchange
    call exch_gather%do(f, domain%parcomm)

    if (domain%parcomm%myid == master_id) then

        err_sum = 0

        do t = 1, domain%partition%Nt
            pn = domain%partition%panel_map(t)
            call domain%partition%tiles_p%tile(t)%getind(is, ie, js, je, ks, ke)
            do k = ks, ke
                do j = js, je
                    do i = is, ie
                        err_sum = err_sum + abs(f%tile(t)%p(i,j,k) - __fun(i,j,k,pn))
                    end do
                end do
            end do
        end do

        if (int(err_sum)==0) then
            print*, 'Test passed!'
        else
            call domain%parcomm%abort('Test not passed! Error! Abort!')
            stop
        end if

    end if

#undef __fun

end subroutine test_gather_exchange


end module test_mod

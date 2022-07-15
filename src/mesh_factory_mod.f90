module mesh_factory_mod

use mesh_mod,    only : mesh_t
use tiles_mod,   only : tiles_t
use parcomm_mod, only : parcomm_global

implicit none

contains

subroutine create_mesh(mesh, partition, metric, halo_width, h_top, points_type, vertical_staggering)

    use partition_mod, only : partition_t
    use metric_mod,    only : metric_t

    type(partition_t), intent(in)  :: partition
    class(metric_t),   intent(in)  :: metric
    type(mesh_t),      intent(out) :: mesh

    integer(kind=4),  intent(in)   :: halo_width
    real(kind=8),     intent(in)   :: h_top
    character(len=*), intent(in)   :: points_type, vertical_staggering

    integer(kind=4) :: t, pind, i, j, k, ts, te, is, ie, js, je, ks, ke, nh, nx, ny, nz
    real(kind=8) :: alpha, beta, eta, hx, hy, hz, vec(3)

    type(tiles_t) :: tiles

    real(kind=8) :: shift_i, shift_j, shift_k

    call partition%get_tiles(points_type, tiles)

    ! define horizontal grid parameters
    shift_i = 0.5_8
    shift_j = 0.5_8
    select case(points_type)
    case("o", "z")
    case("x");         shift_i = 0.0_8
    case("y");         shift_j = 0.0_8
    case("xy", "xyz"); shift_i = 0.0_8
                       shift_j = 0.0_8
    case default
        call parcomm_global%abort("mesh factory error! Wrong points_type: "// points_type)
    end select

    nh = partition%nh
    hx = (metric%alpha1 - metric%alpha0)/real(nh,8)
    hy = (metric%beta1  - metric%beta0 )/real(nh,8)

    ! define vertical grid parameters
    shift_k = 0.0_8
    select case(vertical_staggering)
    case("None")
        if (partition%Nz == 1) then
            !call parcomm_global%print("Single layer model. Setting hz = 1")
            hz = 1.0_8
        else
            hz = 1.0_8/(partition%Nz-1)
        end if
        select case(points_type)
        case("z, xyz")
            call parcomm_global%abort( &
            "Wrong combination of points_type and vertical staggering in mesh factory "//&
                                      points_type//" "//vertical_staggering)
        end select
    case("CharneyPhilips")
        hz = 1.0_8/partition%Nz
        select case(points_type)
        case("o", "x", "y", "xy"); shift_k = 0.5_8
        end select
    end select

    ts = partition%ts
    te = partition%te

    allocate(mesh%tile(ts:te))
    mesh%ts = ts
    mesh%te = te

    mesh%scale          = metric%scale
    mesh%vertical_scale = metric%vertical_scale
    mesh%omega          = metric%omega
    mesh%rotation_axis  = metric%rotation_axis

    do t = ts, te

        call tiles%tile(t)%getind(is, ie, js, je, ks, ke)
        pind = partition%panel_map(t)

        call mesh%tile(t)%init(is, ie, js, je, ks, ke, halo_width)
        mesh%tile(t)%points_type = points_type

        mesh%tile(t)%nx = tiles%Nx
        mesh%tile(t)%ny = tiles%Ny
        mesh%tile(t)%nz = tiles%Nz

        mesh%tile(t)%hx = hx
        mesh%tile(t)%hy = hy
        mesh%tile(t)%hz = hz

        mesh%tile(t)%shift_i = shift_i
        mesh%tile(t)%shift_j = shift_j
        mesh%tile(t)%shift_k = shift_k

        mesh%tile(t)%alpha_0 = metric%alpha0
        mesh%tile(t)%beta_0  = metric%beta0

        do k = ks, ke
            eta = mesh%tile(t)%get_eta(k)
            do j = js - halo_width, je + halo_width
                beta = mesh%tile(t)%get_beta(j)
                do i = is - halo_width, ie + halo_width
                    alpha = mesh%tile(t)%get_alpha(i)
                    vec(1:3) = metric%calculate_r(pind, alpha, beta)

                    mesh%tile(t)%rx(i,j,k) = vec(1)
                    mesh%tile(t)%ry(i,j,k) = vec(2)
                    mesh%tile(t)%rz(i,j,k) = vec(3)
                    mesh%tile(t)%h(i,j,k)  = metric%calculate_h(pind,alpha,beta,eta,0.0_8,h_top)

                    mesh%tile(t)%Q (1:6,i,j,k) = metric%calculate_Q(pind, alpha, beta)
                    mesh%tile(t)%Qi(1:6,i,j,k) = metric%calculate_Qi(pind, alpha, beta)

                    mesh%tile(t)%J(i,j,k)= metric%calculate_J(pind, alpha, beta)

                    mesh%tile(t)%G(1:3,1:3,1:3,i,j,k) = metric%calculate_G(pind,alpha,beta)

                    mesh%tile(t)%a1(1:4,i,j,k) = metric%calculate_a1(pind,alpha,beta)
                    mesh%tile(t)%a2(1:4,i,j,k) = metric%calculate_a2(pind,alpha,beta)
                    mesh%tile(t)%a3(1:4,i,j,k) = metric%calculate_a3(pind,alpha,beta)
                    mesh%tile(t)%b1(1:3,i,j,k) = metric%calculate_b1(pind,alpha,beta)
                    mesh%tile(t)%b2(1:3,i,j,k) = metric%calculate_b2(pind,alpha,beta)
                    mesh%tile(t)%b3(1:4,i,j,k) = metric%calculate_b3(pind,alpha,beta)
                end do
            end do
        end do

    end do


end subroutine create_mesh


end module mesh_factory_mod

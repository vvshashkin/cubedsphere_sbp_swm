module partition_mod
use tiles_mod,   only : tiles_t
use parcomm_mod, only : parcomm_global
implicit none

type, public :: partition_t
    type(tiles_t) :: tiles_o, tiles_x,  tiles_y,  tiles_xy, tiles_z, tiles_xyz, &
                     tiles_u, tiles_v, tiles_p
    integer(kind=4), allocatable :: proc_map(:) !determine belonging of the tile to the specific processor
    integer(kind=4), allocatable :: panel_map(:)!determine belonging of the tile to the specific panel
    integer(kind=4)              :: Nh, Nz      !number of grid points in x/y, z direction for the one panel
    integer(kind=4)              :: num_tiles   !number of tiles at each panel in the partition
    integer(kind=4)              :: num_panels  !
    integer(kind=4)              :: Nt          ! total number of tiles
    integer(kind=4)              :: ts, te      ! start and end index of the tiles belonging to the specific processor
contains
    procedure, public :: init
    procedure, public :: get_tiles
    ! procedure, public :: write_to_txt
    ! procedure, public :: write_to_txt_3d
end type partition_t

contains

subroutine init(this, Nh, Nz, num_tiles, myid, Np, staggering_type, strategy)
    class(partition_t), intent(inout) :: this
    integer(kind=4),    intent(in)    :: Nh        ! num of points in x and y direction at each panel
    integer(kind=4),    intent(in)    :: Nz        ! num of points in z direction at each panel
    integer(kind=4),    intent(in)    :: num_tiles ! num of tiles at each panel
    integer(kind=4),    intent(in)    :: myid, Np  ! myid and num of processors
    character(*),       intent(in)    :: staggering_type, strategy

    integer(kind=4) :: num_panels ! this must bs argument of the routine

    integer(kind=4) :: Nt, t
    integer(kind=4), allocatable :: is(:), ie(:), js(:), je(:), ks(:), ke(:)

!We assume that number of tiles at each panel are the same
    num_panels = 6
    Nt = num_panels*num_tiles

    this%num_panels = num_panels ! need to modify

    ! allocate(this%tile(Nt))
    allocate(this%proc_map(Nt))
    allocate(this%panel_map(Nt))
    this%num_tiles = num_tiles
    this%Nt = Nt
    this%Nh = Nh
    this%Nz = Nz

    call this%tiles_o%init(Nt, Nz, Nh, Nh)

    if (Np>Nt) then
        call parcomm_global%abort('Error!!! Number of tiles is less than number of processors!!!')
    end if
    select case (strategy)
    case ('default')
        call default_strategy(this, Nh, Nz, Np)
    case default
        call parcomm_global%abort('Wrong strategy in partition_mod.f90: '// strategy)
    end select

    this%ts = findloc(this%proc_map, myid, dim=1)
    this%te = findloc(this%proc_map, myid, back = .true., dim=1)

    call this%tiles_x%init (Nt, Nz, Nh,   Nh+1, this%tiles_o%tile)
    call this%tiles_y%init (Nt, Nz, Nh+1, Nh,   this%tiles_o%tile)
    call this%tiles_xy%init(Nt, Nz, Nh+1, Nh+1, this%tiles_o%tile)

    call this%tiles_z%init  (Nt, Nz+1, Nh,   Nh,   this%tiles_o%tile)
    call this%tiles_xyz%init(Nt, Nz+1, Nh+1, Nh+1, this%tiles_o%tile)

    do t = 1, Nt
        if (this%tiles_x%tile(t)%ie == nh) this%tiles_x%tile(t)%ie = nh+1
        if (this%tiles_y%tile(t)%je == nh) this%tiles_y%tile(t)%je = nh+1

        if (this%tiles_xy%tile(t)%ie == nh) this%tiles_xy%tile(t)%ie = nh+1
        if (this%tiles_xy%tile(t)%je == nh) this%tiles_xy%tile(t)%je = nh+1

        if (this%tiles_xyz%tile(t)%ie == nh) this%tiles_xyz%tile(t)%ie = nh+1
        if (this%tiles_xyz%tile(t)%je == nh) this%tiles_xyz%tile(t)%je = nh+1

        if(this%tiles_z%tile(t)%ke  == Nz)  this%tiles_z%tile(t)%ke  = Nz+1
        if(this%tiles_xyz%tile(t)%ke == Nz) this%tiles_xyz%tile(t)%ke = Nz+1
    end do

    if (staggering_type == 'A') then
        this%tiles_u = this%tiles_o
        this%tiles_v = this%tiles_o
        this%tiles_p = this%tiles_o
    else if (staggering_type == 'Ah') then
        this%tiles_u = this%tiles_xy
        this%tiles_v = this%tiles_xy
        this%tiles_p = this%tiles_xy
    else if (staggering_type == 'C') then
        this%tiles_u = this%tiles_x
        this%tiles_v = this%tiles_y
        this%tiles_p = this%tiles_o
    else if (staggering_type == 'Ch') then
        this%tiles_u = this%tiles_y
        this%tiles_v = this%tiles_x
        this%tiles_p = this%tiles_xy
    else
        call parcomm_global%abort('Unknown staggering_type in partition initialization!')
    end if

    ! if(myid == 0) call this%write_to_txt('partition.txt')

end subroutine init

subroutine default_strategy(partition, Nh, Nz, Np)
    type(partition_t), intent(inout) :: partition
    integer(kind=4),    intent(in)   :: Nh, Nz, Np

    integer(kind=4) :: npx_npy(2), wt(np), s(0:np)
    integer(kind=4), allocatable :: wx(:), wy(:)
    integer(kind=4) :: ind, i, j, panel_ind
    integer(kind=4) :: is, ie, js, je, ks, ke

    npx_npy = get_factors(partition%num_tiles)

    allocate(wx(npx_npy(1)), wy(npx_npy(2)))

    call partition_1d(npx_npy(1), Nh, wx)
    call partition_1d(npx_npy(2), Nh, wy)

    ind = 0
    do panel_ind = 1, partition%num_panels
        do j = 1, npx_npy(2)
            do i = 1, npx_npy(1)

                ind = ind+1

                if (j==1) then
                    js = 1
                    je = wy(j)
                else
                    js = partition%tiles_o%tile(ind-npx_npy(1))%js + wy(j-1)
                    je = partition%tiles_o%tile(ind-npx_npy(1))%je + wy(j  )
                end if

                if (i==1) then
                    is = 1
                    ie = wx(i)
                else
                    is = partition%tiles_o%tile(ind-1)%is + wx(i-1)
                    ie = partition%tiles_o%tile(ind-1)%ie + wx(i  )
                endif

                ks = 1
                ke = Nz

                partition%panel_map(ind) = panel_ind
                call partition%tiles_o%tile(ind)%init(is, ie, js, je, ks, ke)

            end do
        end do
    end do

    call partition_1d(Np, partition%Nt, wt)

    s(0) = 0
    do ind = 1, Np
        s(ind) = s(ind-1) + wt(ind)
        partition%proc_map(s(ind-1)+1:s(ind)) = ind-1
    end do

end subroutine default_strategy
! subroutine write_to_txt(this, filename)
!     class(partition_t), intent(in) :: this
!     character(*),       intent(in) :: filename
!
!     integer(kind=4) :: ind
!
!     open(120, file=filename)
!
!     do ind = 1, size(this%tile,1)
!         write(120,*) this%tile(ind)%panel_number, this%tile(ind)%is, this%tile(ind)%ie &
!                                            , this%tile(ind)%js, this%tile(ind)%je, this%proc_map(ind)
!     end do
!
!     close(120)
!
! end subroutine
! subroutine write_to_txt_3d(this, filename)
!
!     use topology_mod, only : ex, ey, r
!
!     class(partition_t), intent(in) :: this
!     character(*),       intent(in) :: filename
!
!     integer(kind=4) :: ind, xyz_s(1:3), xyz_e(1:3)
!
!     open(120, file=filename)
!
!     do ind = 1, size(this%tile,1)
!
!         xyz_s(1:3) = r(1:3,this%tile(ind)%panel_number)*this%Nh + (this%tile(ind)%is-1)*ex(1:3, this%tile(ind)%panel_number) + &
!                                                              (this%tile(ind)%js-1)*ey(1:3, this%tile(ind)%panel_number)
!         xyz_e(1:3) = r(1:3,this%tile(ind)%panel_number)*this%Nh + (this%tile(ind)%ie-1)*ex(1:3, this%tile(ind)%panel_number) + &
!                                                              (this%tile(ind)%je-1)*ey(1:3, this%tile(ind)%panel_number)
!
!         write(120,'(8I4)') this%tile(ind)%panel_number, xyz_s(1), xyz_e(1), xyz_s(2), xyz_e(2), xyz_s(3), xyz_e(3), this%proc_map(ind)
!     end do
!
!     close(120)
!
! end subroutine write_to_txt_3d

subroutine partition_1d(Np, N, work)
    integer(kind=4), intent(in)  :: Np, N
    integer(kind=4), intent(out) :: work(Np)

    real(kind=8) :: mean_work, deficite

    mean_work = N/Np
    work(1:Np)=mean_work
    deficite = mod(N,Np)
    work(1:deficite) = work(1:deficite)+1

end subroutine partition_1d

function get_factors(n) result(factors)
    integer(kind=4), intent(in)  :: n
    integer(kind=4) :: factors(2)

    integer(kind=4) :: f1

    f1 = int(sqrt(real(n,8)))

    do
        if (mod(n,f1)==0) exit
        f1 = f1-1
    end do

    factors = (/n/f1, f1/)

end function get_factors
subroutine get_tiles(this, points_type, tiles)

    use parcomm_mod, only : parcomm_global

    class(partition_t), intent(in)  :: this
    character(len=*),   intent(in)  :: points_type
    type(tiles_t),      intent(out) :: tiles

    select case(points_type)
    case('o')
        tiles = this%tiles_o
    case('x')
        tiles = this%tiles_x
    case('y')
        tiles = this%tiles_y
    case('xy')
        tiles = this%tiles_xy
    case('z')
        tiles = this%tiles_z
    case('xyz')
        tiles = this%tiles_xyz
    case('p')
        tiles = this%tiles_p
    case('u')
        tiles = this%tiles_u
    case('v')
        tiles = this%tiles_v
    case default
        call parcomm_global%abort("Wrong points_type in get_tiles: "//points_type)
    end select

end subroutine get_tiles
end module partition_mod

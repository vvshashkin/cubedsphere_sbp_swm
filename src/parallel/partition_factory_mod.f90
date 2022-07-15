module partition_factory_mod

use partition_mod, only : partition_t

implicit none

contains

subroutine create_partition(partition, Nh, Nz, num_panels, myid, Np, staggering_type, strategy)

    type(partition_t),  intent(out) :: partition
    integer(kind=4),    intent(in)  :: Nh         ! num of points in hor direction
    integer(kind=4),    intent(in)  :: Nz         ! num of points in z direction
    integer(kind=4),    intent(in)  :: num_panels !
    integer(kind=4),    intent(in)  :: myid, Np   ! myid and num of processors
    character(*),       intent(in)  :: staggering_type, strategy

    integer(kind=4) :: t, Ntiles

    partition%num_panels = num_panels

    allocate(this%tile_o(num_panels*num_tiles))
    allocate(this%tile_x(num_panels*num_tiles))
    allocate(this%tile_y(num_panels*num_tiles))
    allocate(this%tile_xy(num_panels*num_tiles))
    allocate(this%tile(num_panels*num_tiles))
    allocate(this%proc_map(num_panels*num_tiles))
    allocate(this%panel_map(num_panels*num_tiles))
    this%num_tiles = num_tiles
    this%Nh = Nh
    this%Nz = Nz

end subroutine create_partition

end module partition_factory_mod

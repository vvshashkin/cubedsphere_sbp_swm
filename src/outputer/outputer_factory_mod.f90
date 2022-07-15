module outputer_factory_mod

use outputer_abstract_mod, only : outputer_t, outputer_vector_t
use domain_mod,            only : domain_t
use parcomm_mod,           only : parcomm_global

implicit none

contains

subroutine create_master_paneled_outputer(outputer, points_type, domain, master_id)

    use master_paneled_outputer_mod, only : master_paneled_outputer_t
    use grid_field_factory_mod,      only : create_grid_field_global
    use exchange_factory_mod,        only : create_gather_exchange

    class(outputer_t), allocatable, intent(out) :: outputer
    character(len=*),               intent(in)  :: points_type
    type(domain_t),                 intent(in)  :: domain
    integer(kind=4),  optional,     intent(in)  :: master_id

    type(master_paneled_outputer_t), allocatable :: master_outputer
    integer(kind=4) :: i, master_id_loc

    master_id_loc = 0
    if (present(master_id)) master_id_loc = master_id

    allocate(master_outputer)

    if (domain%parcomm%myid == master_id_loc) then
        call domain%partition%get_tiles(points_type, master_outputer%tiles)
        call create_grid_field_global(master_outputer%buf, 0, 0, master_outputer%tiles)

    end if

    master_outputer%master_id = master_id_loc

    call create_gather_exchange(master_outputer%gather_exch, points_type, &
                                domain%parcomm, domain%partition, master_id_loc)

    call move_alloc(master_outputer, outputer)

end subroutine create_master_paneled_outputer

subroutine create_latlon_outputer(outputer, Nlat, Nlon, scalar_grid_type, domain, &
                                  master_id, is_z_interfaces)

    use latlon_outputer_mod,         only : latlon_outputer_t
    use grid_field_factory_mod,      only : create_grid_field_global, &
                                            create_grid_field
    use exchange_factory_mod,        only : create_gather_exchange
    use domain_factory_mod,          only : create_domain
    use parcomm_mod,                 only : parcomm_t
    use parcomm_factory_mod,         only : create_group_parcomm
    use regrid_factory_mod,          only : create_latlon_regrid

    class(outputer_t), allocatable, intent(out) :: outputer
    integer(kind=4),                intent(in)  :: Nlon, Nlat
    character(len=*),               intent(in)  :: scalar_grid_type
    type(domain_t),                 intent(in)  :: domain
    integer(kind=4),  optional,     intent(in)  :: master_id
    logical,          optional,     intent(in)  :: is_z_interfaces

    type(latlon_outputer_t), allocatable :: latlon_outputer
    integer(kind=4) :: i, master_id_loc
    character(len=:), allocatable :: points_type
    type(parcomm_t) :: master_parcomm
    logical :: is_z_interfaces_loc

    master_id_loc = 0
    if (present(master_id)) master_id_loc = master_id

    is_z_interfaces_loc = .false.
    if(present(is_z_interfaces)) is_z_interfaces_loc = is_z_interfaces

    allocate(latlon_outputer)

    call create_group_parcomm(master_parcomm, domain%parcomm, [master_id_loc])

    if(domain%parcomm%myid == master_id_loc) then
        !A - staggering is WORKAROUND
        call create_domain(latlon_outputer%regrid_domain, "cube", 'A', &
                           domain%partition%nh, 1, parcomm=master_parcomm)

        call create_latlon_regrid(latlon_outputer%regrid, &
                                  latlon_outputer%regrid_domain,&
                                  Nlat=Nlat,Nlon=Nlon,&
                                  interp_type="cubic",&
                                  scalar_grid_type=scalar_grid_type)

        if(scalar_grid_type == "A") then
            call create_grid_field(latlon_outputer%regrid_work, 0, 0, &
                                   latlon_outputer%regrid_domain%mesh_o)
        else if(scalar_grid_type == "Ah") then
            call create_grid_field(latlon_outputer%regrid_work, 0, 0, &
                                   latlon_outputer%regrid_domain%mesh_xy)
        end if
        latlon_outputer%Nlat = Nlat
        latlon_outputer%Nlon = Nlon
        allocate(latlon_outputer%buffer(Nlon,Nlat,1))
    end if

    if(is_z_interfaces_loc) then
        if(scalar_grid_type == 'A') then
            latlon_outputer%tiles = domain%partition%tiles_z
            points_type = "z"
        else if(scalar_grid_type == 'Ah') then
            latlon_outputer%tiles = domain%partition%tiles_xyz
            points_type = "xyz"
        else
            call parcomm_global%abort("latlon outputer supports only A and Ah scalar grids")
        end if
    else
        if(scalar_grid_type == 'A') then
            latlon_outputer%tiles = domain%partition%tiles_o
            points_type = "o"
        else if(scalar_grid_type == 'Ah') then
            latlon_outputer%tiles = domain%partition%tiles_xy
            points_type = "xy"
        else
            call parcomm_global%abort("latlon outputer supports only A and Ah scalar grids")
        end if
    end if

    if (domain%parcomm%myid == master_id_loc) then
        call create_grid_field_global(latlon_outputer%exchange_buf, 0, 0,  latlon_outputer%tiles)
    end if

    latlon_outputer%master_id = master_id_loc

    call create_gather_exchange(latlon_outputer%gather_exch, points_type, &
                                domain%parcomm, domain%partition, master_id_loc)

    call move_alloc(latlon_outputer, outputer)

end subroutine create_latlon_outputer

subroutine create_latlon_vec_outputer(outputer, Nlat, Nlon, staggering, &
                                      components_type, domain, master_id)

    use latlon_outputer_mod,         only : latlon_vec_outputer_t
    use grid_field_factory_mod,      only : create_grid_field_global, &
                                            create_grid_field
    use exchange_factory_mod,        only : create_gather_exchange
    use domain_factory_mod,          only : create_domain
    use parcomm_mod,                 only : parcomm_t
    use parcomm_factory_mod,         only : create_group_parcomm
    use regrid_factory_mod,          only : create_latlon_vector_regrid

    class(outputer_vector_t), allocatable, intent(out) :: outputer
    integer(kind=4),                       intent(in)  :: Nlon, Nlat
    character(len=*),                      intent(in)  :: staggering
    character(len=*),                      intent(in)  :: components_type
    type(domain_t),                        intent(in)  :: domain
    integer(kind=4),  optional,            intent(in)  :: master_id

    type(latlon_vec_outputer_t), allocatable :: latlon_vec_outputer
    integer(kind=4) :: i, master_id_loc
    type(parcomm_t) :: master_parcomm

    master_id_loc = 0
    if (present(master_id)) master_id_loc = master_id

    allocate(latlon_vec_outputer)

    call create_group_parcomm(master_parcomm, domain%parcomm, [master_id_loc])

    if(domain%parcomm%myid == master_id_loc) then
        call create_domain(latlon_vec_outputer%regrid_domain, "cube", staggering, &
                           domain%partition%nh, nz=1, parcomm=master_parcomm)

        call create_latlon_vector_regrid(latlon_vec_outputer%regrid, &
                                         latlon_vec_outputer%regrid_domain,&
                                         Nlat=Nlat,Nlon=Nlon,&
                                         interp_type="cubic",&
                                         vector_grid_type=staggering, &
                                         components_type=components_type)

        call create_grid_field(latlon_vec_outputer%regrid_work_u, 0, 0, &
                               latlon_vec_outputer%regrid_domain%mesh_u)
        call create_grid_field(latlon_vec_outputer%regrid_work_v, 0, 0, &
                               latlon_vec_outputer%regrid_domain%mesh_v)

        latlon_vec_outputer%Nlat = Nlat
        latlon_vec_outputer%Nlon = Nlon
        allocate(latlon_vec_outputer%ulatlon(Nlon,Nlat,1))
        allocate(latlon_vec_outputer%vlatlon(Nlon,Nlat,1))
    end if

    latlon_vec_outputer%tiles_u = domain%partition%tiles_u
    latlon_vec_outputer%tiles_v = domain%partition%tiles_v

    if (domain%parcomm%myid == master_id_loc) then
        call create_grid_field_global(latlon_vec_outputer%exchange_buf_u, 0, 0, &
                                      latlon_vec_outputer%tiles_u)
        if(staggering == "A" .or. &
           staggering == "Ah") then
            latlon_vec_outputer%exchange_buf_v = latlon_vec_outputer%exchange_buf_u
        else
            call create_grid_field_global(latlon_vec_outputer%exchange_buf_v, 0, 0,&
                                          latlon_vec_outputer%tiles_v)
        end if
    end if

    latlon_vec_outputer%master_id = master_id_loc

    call create_gather_exchange(latlon_vec_outputer%gather_exch_u, "u", &
                                domain%parcomm, domain%partition, master_id_loc)
    if(staggering == "A" .or. staggering == "Ah") then
        latlon_vec_outputer%gather_exch_v = latlon_vec_outputer%gather_exch_u
    else
        call create_gather_exchange(latlon_vec_outputer%gather_exch_v, "v", &
                                    domain%parcomm, domain%partition, master_id_loc)
    end if

    call move_alloc(latlon_vec_outputer, outputer)

end subroutine create_latlon_vec_outputer


! function create_mpi_paneled_outputer(partition) result(outputer)
!
!     use mpi_paneled_outputer_mod, only : mpi_paneled_outputer_t
!     use partition_mod,            only : partition_t
!     use mpi
!
!     type(partition_t), intent(in)  :: partition
!     type(mpi_paneled_outputer_t)   :: outputer
!
!     integer(kind=mpi_offset_kind) :: disp
!     integer(kind=4) :: type_local(partition%ts:partition%te)
!     integer(kind=4) :: type_global(partition%ts:partition%te)
!     integer(kind=4) :: type_domain
!     integer(kind=4) :: dim, gsize(4), lsize(4), start(4), size
!     integer(kind=4) :: myid, ierr
!     integer(kind=4) :: t, is, ie, js, je, ks, ke
!
!     call mpi_comm_rank(mpi_comm_world, myid, ierr)
!
!     dim = 4
!
!     gsize = [partition%Nh, partition%Nh, partition%num_panels, partition%Nz]
!     lsize = [partition%Nh, partition%Nh, partition%num_panels, partition%Nz]
!     start = [0,0,0,0]
!
!     call MPI_Type_create_subarray(dim, gsize, lsize, start, MPI_ORDER_FORTRAN, mpi_real4, type_domain, ierr)
!     call MPI_Type_commit(type_domain, ierr)
!     call mpi_type_size(type_domain, size, ierr)
!     call mpi_type_free(type_domain, ierr)
!
!     disp = int(size,kind=mpi_offset_kind)
!
!     do t = partition%ts, partition%te
!
!         is = partition%tile(t)%is; ie = partition%tile(t)%ie
!         js = partition%tile(t)%js; je = partition%tile(t)%je
!         ks = partition%tile(t)%ks; ke = partition%tile(t)%ke
!
!         gsize = [partition%Nh, partition%Nh, partition%num_panels, partition%Nz]
!         lsize = [ie-is+1, je-js+1, 1, ke-ks+1]
!         start = [is-1, js-1,partition%tile(t)%panel_number-1, ks-1]
!
!         call MPI_Type_create_subarray(dim, gsize, lsize, start, MPI_ORDER_FORTRAN, mpi_real4, type_global(t), ierr)
!         call MPI_Type_commit(type_global(t), ierr)
!
!         gsize = [ie-is+1, je-js+1, 1, ke-ks+1]
!         lsize = [ie-is+1, je-js+1, 1, ke-ks+1]
!         start = [0, 0, 0, 0]
!
!         call MPI_Type_create_subarray(dim, gsize, lsize, start, MPI_ORDER_FORTRAN, mpi_real4, type_local(t), ierr)
!         call MPI_Type_commit(type_local(t), ierr)
!
!     end do
!
!     outputer = mpi_paneled_outputer_t(record_disp = disp,      &
!                                       type_local = type_local, &
!                                       type_global = type_global )
!
! end function create_mpi_paneled_outputer

end module outputer_factory_mod

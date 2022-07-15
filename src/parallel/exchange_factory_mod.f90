module exchange_factory_mod

use parcomm_mod,   only : parcomm_t
use partition_mod, only : partition_t
use topology_mod,  only : topology_t
use tile_mod,      only : tile_t
use tiles_mod,     only : tiles_t

implicit none

contains

function create_o_points_halo_exchange(partition, parcomm, topology, halo_width, halo_type) result(exchange)

    use exchange_halo_mod, only : exchange_2D_halo_t

    type(partition_t), intent(in)    :: partition
    class(topology_t), intent(in) :: topology
    integer(kind=4),   intent(in) :: halo_width
    character(*),      intent(in) :: halo_type
    type(parcomm_t),   intent(in) :: parcomm

    type(exchange_2D_halo_t) :: exchange

    exchange = create_o_z_points_halo_exchange(partition%tiles_o, partition, parcomm, topology, &
                                             halo_width, halo_type)
end function create_o_points_halo_exchange

function create_z_points_halo_exchange(partition, parcomm, topology, halo_width, halo_type) result(exchange)

    use exchange_halo_mod, only : exchange_2D_halo_t

    type(partition_t), intent(in)    :: partition
    class(topology_t), intent(in) :: topology
    integer(kind=4),   intent(in) :: halo_width
    character(*),      intent(in) :: halo_type
    type(parcomm_t),   intent(in) :: parcomm

    type(exchange_2D_halo_t) :: exchange

    exchange = create_o_z_points_halo_exchange(partition%tiles_z, partition, parcomm, topology, &
                                             halo_width, halo_type)
end function create_z_points_halo_exchange

function create_o_z_points_halo_exchange(tiles, partition, parcomm, topology, halo_width, &
                                       halo_type) result(exchange)

    use exchange_halo_mod, only : exchange_2D_halo_t

    type(partition_t), intent(in) :: partition
    type(tiles_t), target, intent(in) :: tiles
    class(topology_t), intent(in) :: topology
    integer(kind=4),   intent(in) :: halo_width
    character(*),      intent(in) :: halo_type
    type(parcomm_t),   intent(in) :: parcomm

    type(exchange_2D_halo_t) :: exchange

    type(tile_t), pointer :: local_tile, remote_tile

    type(tile_t),     dimension(partition%Nt**2) :: send_tile, recv_tile
    integer(kind=4),  dimension(partition%Nt**2) :: send_to_proc_id, recv_from_proc_id,   &
                                                    recv_to_tile_ind, send_from_tile_ind, &
                                                    send_tag, recv_tag
    integer(kind=4),  dimension(partition%Nt**2) :: send_pts_num, recv_pts_num
    integer(kind=4),  dimension(partition%Nt**2) :: send_i_step, send_j_step
    character(len=1), dimension(partition%Nt**2) :: first_dim_index

    type(tile_t) :: temp_tile, halo_tile

    integer(kind=4) :: local_ind, remote_ind, exch_num, ind,ts ,te, rnum, snum
    integer(kind=4) :: local_panel, remote_panel
    logical :: is_intersection

    if (halo_width > partition%Nh) then
        write(*,*) 'Error! Halo is too wide! Abort!'
        stop
    end if

    rnum = 0
    snum = 0

    do local_ind = partition%ts, partition%te

        local_tile => tiles%tile(local_ind)
        local_panel = partition%panel_map(local_ind)

        do remote_ind = 1, partition%Nt

            if (remote_ind == local_ind) cycle

            remote_tile => tiles%tile(remote_ind)
            remote_panel = partition%panel_map(remote_ind)

            call topology%transform_tile_coords(remote_panel, remote_tile, &
                                                local_panel,  temp_tile,   &
                                                tiles%Nx, tiles%Ny)

            !recv exchange
            call find_tiles_halo_intersection(local_tile, halo_width, halo_width, &
                          temp_tile, halo_type, recv_tile(rnum+1), is_intersection)

            if (is_intersection) then
                rnum = rnum + 1

                recv_to_tile_ind(rnum) = local_ind
                recv_pts_num(rnum) = recv_tile(rnum)%get_points_num()
                recv_from_proc_id(rnum) = partition%proc_map(remote_ind)
                recv_tag(rnum) = partition%Nt*(remote_ind-1) + local_ind
            end if

            !send exchange
            call topology%find_basis_orientation(remote_panel, local_panel, &
                                                  send_i_step(snum+1), send_j_step(snum+1), &
                                                  first_dim_index(snum+1))
            call find_tiles_halo_intersection(temp_tile, halo_width, halo_width, &
                        local_tile, halo_type, send_tile(snum+1), is_intersection)

            if (is_intersection) then
                snum = snum + 1
                send_from_tile_ind(snum) = local_ind
                send_pts_num(snum) = send_tile(snum)%get_points_num()
                send_to_proc_id(snum) = partition%proc_map(remote_ind)
                send_tag(snum) = partition%Nt*(local_ind-1) + remote_ind
            end if
        end do
    end do

    exchange%send_tile          = send_tile(1:snum)
    exchange%send_points_num    = send_pts_num(1:snum)
    exchange%send_from_tile_ind = send_from_tile_ind(1:snum)
    exchange%send_to_proc_id    = send_to_proc_id(1:snum)
    exchange%send_tag           = send_tag(1:snum)
    exchange%send_i_step        = send_i_step(1:snum)
    exchange%send_j_step        = send_j_step(1:snum)
    exchange%first_dim_index    = first_dim_index(1:snum)
    exchange%send_number        = snum
    exchange%recv_tile          = recv_tile(1:rnum)
    exchange%recv_to_tile_ind   = recv_to_tile_ind(1:rnum)
    exchange%recv_from_proc_id  = recv_from_proc_id(1:rnum)
    exchange%recv_points_num    = recv_pts_num(1:rnum)
    exchange%recv_tag           = recv_tag(1:rnum)
    exchange%recv_number        = rnum

    allocate(exchange%send_buff(snum)   , exchange%recv_buff(rnum)   )
    allocate(exchange%mpi_send_req(snum), exchange%mpi_recv_req(rnum))

    ! points num multiplied by 2 for vec excahnges
    do ind = 1, rnum
        call exchange%recv_buff(ind)%init(2*recv_pts_num(ind))
    end do
    do ind = 1, snum
        call exchange%send_buff(ind)%init(2*send_pts_num(ind))
    end do

end function create_o_z_points_halo_exchange
function create_xy_points_halo_exchange(partition, parcomm, topology, halo_width, halo_type) result(exchange)

    use exchange_halo_mod, only : exchange_2D_halo_t

    type(partition_t), intent(in)    :: partition
    class(topology_t), intent(in) :: topology
    integer(kind=4),   intent(in) :: halo_width
    character(*),      intent(in) :: halo_type
    type(parcomm_t),   intent(in) :: parcomm

    type(exchange_2D_halo_t) :: exchange

    exchange = create_xy_xyz_points_halo_exchange(partition%tiles_xy, partition, parcomm, topology, &
                                             halo_width, halo_type)
end function create_xy_points_halo_exchange

function create_xyz_points_halo_exchange(partition, parcomm, topology, halo_width, halo_type) result(exchange)

    use exchange_halo_mod, only : exchange_2D_halo_t

    type(partition_t), intent(in)    :: partition
    class(topology_t), intent(in) :: topology
    integer(kind=4),   intent(in) :: halo_width
    character(*),      intent(in) :: halo_type
    type(parcomm_t),   intent(in) :: parcomm

    type(exchange_2D_halo_t) :: exchange

    exchange = create_xy_xyz_points_halo_exchange(partition%tiles_xyz, partition, parcomm, topology, &
                                             halo_width, halo_type)
end function create_xyz_points_halo_exchange
function create_xy_xyz_points_halo_exchange(tiles, partition, parcomm, topology, &
                                            halo_width, halo_type) result(exchange)

    use exchange_halo_mod, only : exchange_2D_halo_t

    type(tiles_t), target, intent(in) :: tiles
    type(partition_t),     intent(in) :: partition
    class(topology_t),     intent(in) :: topology
    integer(kind=4),       intent(in) :: halo_width
    character(*),          intent(in) :: halo_type
    type(parcomm_t),       intent(in) :: parcomm

    type(exchange_2D_halo_t) :: exchange

    type(tile_t), pointer :: local_tile, remote_tile

    type(tile_t),     dimension(partition%Nt**2) :: send_tile, recv_tile
    integer(kind=4),  dimension(partition%Nt**2) :: send_to_proc_id, recv_from_proc_id,   &
                                                    recv_to_tile_ind, send_from_tile_ind, &
                                                    send_tag, recv_tag
    integer(kind=4),  dimension(partition%Nt**2) :: send_pts_num, recv_pts_num
    integer(kind=4),  dimension(partition%Nt**2) :: send_i_step, send_j_step
    character(len=1), dimension(partition%Nt**2) :: first_dim_index

    type(tile_t) :: temp_tile, halo_tile

    integer(kind=4) :: local_ind, remote_ind, exch_num, ind,ts ,te, rnum, snum
    integer(kind=4) :: local_panel, remote_panel
    logical :: is_intersection

    if (halo_width > partition%Nh) then
        write(*,*) 'Error! Halo is too wide! Abort!'
        stop
    end if

    rnum = 0
    snum = 0

    do local_ind = partition%ts, partition%te

        local_tile => tiles%tile(local_ind)
        local_panel = partition%panel_map(local_ind)

        do remote_ind = 1, partition%Nt

            if (remote_ind == local_ind) cycle

            remote_tile => tiles%tile(remote_ind)
            remote_panel = partition%panel_map(remote_ind)

            call topology%transform_tile_coords(remote_panel, remote_tile, &
                                                local_panel,  temp_tile,   &
                                                tiles%Nx, tiles%Ny)

            !recv exchange
            call find_tiles_halo_intersection(local_tile, halo_width, halo_width, &
                          temp_tile, halo_type, recv_tile(rnum+1), is_intersection)

            if (is_intersection) then
                rnum = rnum + 1

                recv_to_tile_ind(rnum) = local_ind
                recv_pts_num(rnum) = recv_tile(rnum)%get_points_num()
                recv_from_proc_id(rnum) = partition%proc_map(remote_ind)
                recv_tag(rnum) = partition%Nt*(remote_ind-1) + local_ind
            end if

            !send exchange
            call topology%find_basis_orientation(remote_panel, local_panel, &
                                                  send_i_step(snum+1), send_j_step(snum+1), &
                                                  first_dim_index(snum+1))
            call find_tiles_halo_intersection(temp_tile, halo_width, halo_width, &
                        local_tile, halo_type, send_tile(snum+1), is_intersection)

            if (is_intersection) then
                snum = snum + 1
                send_from_tile_ind(snum) = local_ind
                send_pts_num(snum) = send_tile(snum)%get_points_num()
                send_to_proc_id(snum) = partition%proc_map(remote_ind)
                send_tag(snum) = partition%Nt*(local_ind-1) + remote_ind
            end if
        end do
    end do

    exchange%send_tile          = send_tile(1:snum)
    exchange%send_points_num    = send_pts_num(1:snum)
    exchange%send_from_tile_ind = send_from_tile_ind(1:snum)
    exchange%send_to_proc_id    = send_to_proc_id(1:snum)
    exchange%send_tag           = send_tag(1:snum)
    exchange%send_i_step        = send_i_step(1:snum)
    exchange%send_j_step        = send_j_step(1:snum)
    exchange%first_dim_index    = first_dim_index(1:snum)
    exchange%send_number        = snum
    exchange%recv_tile          = recv_tile(1:rnum)
    exchange%recv_to_tile_ind   = recv_to_tile_ind(1:rnum)
    exchange%recv_from_proc_id  = recv_from_proc_id(1:rnum)
    exchange%recv_points_num    = recv_pts_num(1:rnum)
    exchange%recv_tag           = recv_tag(1:rnum)
    exchange%recv_number        = rnum

    allocate(exchange%send_buff(snum)   , exchange%recv_buff(rnum)   )
    allocate(exchange%mpi_send_req(snum), exchange%mpi_recv_req(rnum))

    ! points num multiplied by 2 for vec excahnges
    do ind = 1, rnum
        call exchange%recv_buff(ind)%init(2*recv_pts_num(ind))
    end do
    do ind = 1, snum
        call exchange%send_buff(ind)%init(2*send_pts_num(ind))
    end do

end function create_xy_xyz_points_halo_exchange
function create_symmetric_halo_vec_exchange_Ch(partition, parcomm, topology, halo_width, halo_type) result(exchange)

    use exchange_halo_Ch_mod, only : exchange_2D_halo_Ch_t

    type(partition_t), target, intent(in) :: partition
    class(topology_t),         intent(in) :: topology
    integer(kind=4),           intent(in) :: halo_width
    character(*),              intent(in) :: halo_type
    type(parcomm_t),           intent(in) :: parcomm

    type(exchange_2D_halo_Ch_t) :: exchange

    exchange%exch_u = create_symm_halo_vec_exchange_U_points(partition, parcomm, topology, halo_width, halo_type)
    exchange%exch_v = create_symm_halo_vec_exchange_V_points(partition, parcomm, topology, halo_width, halo_type)

end function create_symmetric_halo_vec_exchange_Ch

function create_symmetric_halo_vec_exchange_C(partition, parcomm, topology, halo_width, halo_type) result(exchange)

    use exchange_halo_C_mod, only : exchange_2D_halo_C_t

    type(partition_t), target, intent(in)    :: partition
    class(topology_t), intent(in) :: topology
    integer(kind=4), intent(in) :: halo_width
    character(*),    intent(in) :: halo_type
    type(parcomm_t), intent(in) :: parcomm

    type(exchange_2D_halo_C_t) :: exchange

    exchange%exch_u = create_symm_halo_vec_exchange_U_points(partition, parcomm, topology, halo_width, halo_type)
    exchange%exch_v = create_symm_halo_vec_exchange_V_points(partition, parcomm, topology, halo_width, halo_type)

end function create_symmetric_halo_vec_exchange_C

function create_symm_halo_vec_exchange_U_points(partition, parcomm, topology, halo_width, halo_type) result(exchange)

    use exchange_halo_mod, only : exchange_2d_halo_t

    type(partition_t), target, intent(in)    :: partition
    class(topology_t), intent(in) :: topology
    integer(kind=4),  intent(in) :: halo_width
    character(*),     intent(in) :: halo_type
    type(parcomm_t),  intent(in) :: parcomm

    type(exchange_2d_halo_t) :: exchange

    type(tile_t),     dimension(partition%Nt**2) :: send_tile, recv_tile
    integer(kind=4),  dimension(partition%Nt**2) :: send_to_proc_id,  recv_from_proc_id,  &
                                                    recv_to_tile_ind, send_from_tile_ind, &
                                                    send_tag, recv_tag
    integer(kind=4),  dimension(partition%Nt**2) :: send_pts_num, recv_pts_num
    integer(kind=4),  dimension(partition%Nt**2) :: send_i_step, send_j_step
    character(len=1), dimension(partition%Nt**2) :: first_dim_index

    type(tile_t), pointer :: local_tile, remote_tile

    type(tile_t) :: temp_tile, halo_tile

    integer(kind=4) :: local_ind, remote_ind, exch_num, ind,ts ,te, rnum, snum
    integer(kind=4) :: local_panel, remote_panel, halo_width_x, halo_width_y
    logical :: is_intersection

    if (halo_width > partition%Nh) then
        write(*,*) 'Error! Halo is too wide! Abort!'
        stop
    end if

    rnum = 0
    snum = 0

    do local_ind = partition%ts, partition%te

        local_panel = partition%panel_map(local_ind)

        do remote_ind = 1, partition%Nt

            if (remote_ind == local_ind) cycle

            halo_width_x = halo_width
            halo_width_y = halo_width

            !recv part
            local_tile => partition%tiles_x%tile(local_ind)

            remote_panel = partition%panel_map(remote_ind)

            call topology%find_basis_orientation(remote_panel, local_panel,       &
                                        send_i_step(snum+1), send_j_step(snum+1), &
                                        first_dim_index(snum+1))

            if (first_dim_index(snum+1) == 'i') then
            !recv from u points to u points
                remote_tile => partition%tiles_x%tile(remote_ind)
            else if (first_dim_index(snum+1) == 'j') then
            !recv from v points to u points
                remote_tile => partition%tiles_y%tile(remote_ind)
            end if

            call topology%transform_tile_coords(remote_panel, remote_tile,     &
                                                local_panel,  temp_tile,       &
                                        partition%tiles_x%Nx, partition%tiles_x%Ny)

            !recv exchange. +1 to x halo_width to account for edge points
            if (remote_panel /= local_panel) halo_width_x = halo_width+1

            call find_tiles_halo_intersection(local_tile, halo_width_x, halo_width_y, &
                            temp_tile, halo_type, recv_tile(rnum+1), is_intersection)

            if (is_intersection) then
                rnum = rnum + 1
                recv_to_tile_ind(rnum) = local_ind
                recv_pts_num(rnum) = recv_tile(rnum)%get_points_num()
                recv_from_proc_id(rnum) = partition%proc_map(remote_ind)
                recv_tag(rnum) = partition%Nt*(remote_ind-1) + local_ind
            end if
            !end of recv part

            halo_width_x = halo_width
            halo_width_y = halo_width

            if (first_dim_index(snum+1) == 'i') then
                !send from u points to u points!local/remote tiles already assigned and transformed
                !send exchange. +1 to x halo_width to account for edge points
                if (remote_panel /= local_panel) halo_width_x = halo_width+1
                call find_tiles_halo_intersection(temp_tile, halo_width_x, halo_width_y, &
                              local_tile, halo_type, send_tile(snum+1), is_intersection)
            else if (first_dim_index(snum+1) == 'j') then
                !send from v points to u points
                local_tile  => partition%tiles_y%tile(local_ind)
                remote_tile => partition%tiles_x%tile(remote_ind)
                call topology%transform_tile_coords(remote_panel, remote_tile, &
                                                    local_panel, temp_tile,     &
                                                    partition%tiles_y%Nx, partition%tiles_y%Ny)
                !send exchange. +1 to y halo_width to account for edge points
                if (remote_panel /= local_panel) halo_width_y = halo_width+1
                call find_tiles_halo_intersection(temp_tile, halo_width_x, halo_width_y, &
                              local_tile, halo_type, send_tile(snum+1), is_intersection)
            end if

            if (is_intersection) then
                snum = snum + 1
                send_from_tile_ind(snum) = local_ind
                send_pts_num(snum) = send_tile(snum)%get_points_num()
                send_to_proc_id(snum) = partition%proc_map(remote_ind)
                send_tag(snum) = partition%Nt*(local_ind-1) + remote_ind
            end if
        end do
    end do

    exchange%send_tile          = send_tile(1:snum)
    exchange%send_points_num    = send_pts_num(1:snum)
    exchange%send_from_tile_ind = send_from_tile_ind(1:snum)
    exchange%send_to_proc_id    = send_to_proc_id(1:snum)
    exchange%send_tag           = send_tag(1:snum)
    exchange%send_i_step        = send_i_step(1:snum)
    exchange%send_j_step        = send_j_step(1:snum)
    exchange%first_dim_index    = first_dim_index(1:snum)
    exchange%send_number        = snum
    exchange%recv_tile          = recv_tile(1:rnum)
    exchange%recv_to_tile_ind   = recv_to_tile_ind(1:rnum)
    exchange%recv_from_proc_id  = recv_from_proc_id(1:rnum)
    exchange%recv_points_num    = recv_pts_num(1:rnum)
    exchange%recv_tag           = recv_tag(1:rnum)
    exchange%recv_number        = rnum

    allocate(exchange%send_buff(snum)   , exchange%recv_buff(rnum)   )
    allocate(exchange%mpi_send_req(snum), exchange%mpi_recv_req(rnum))

    do ind = 1, rnum
        call exchange%recv_buff(ind)%init(recv_pts_num(ind))
    end do
    do ind = 1, snum
        call exchange%send_buff(ind)%init(send_pts_num(ind))
    end do

end function create_symm_halo_vec_exchange_U_points
function create_symm_halo_vec_exchange_V_points(partition, parcomm, topology, halo_width, halo_type) result(exchange)

    use exchange_halo_mod, only : exchange_2d_halo_t

    type(partition_t), target, intent(in)    :: partition
    class(topology_t), intent(in) :: topology
    integer(kind=4), intent(in)   :: halo_width
    character(*),    intent(in)   :: halo_type
    type(parcomm_t), intent(in)   :: parcomm

    type(exchange_2d_halo_t) :: exchange

    type(tile_t),     dimension(partition%Nt**2) :: send_tile, recv_tile
    integer(kind=4),  dimension(partition%Nt**2) :: send_to_proc_id, recv_from_proc_id,   &
                                                    recv_to_tile_ind, send_from_tile_ind, &
                                                    send_tag, recv_tag
    integer(kind=4),  dimension(partition%Nt**2) :: send_pts_num, recv_pts_num
    integer(kind=4),  dimension(partition%Nt**2) :: send_i_step, send_j_step
    character(len=1), dimension(partition%Nt**2) :: first_dim_index

    type(tile_t), pointer :: local_tile, remote_tile

    type(tile_t) :: temp_tile, halo_tile

    integer(kind=4) :: local_ind, remote_ind, exch_num, ind,ts ,te, rnum, snum
    integer(kind=4) :: local_panel, remote_panel, halo_width_x, halo_width_y
    logical :: is_intersection

    if (halo_width > partition%Nh) then
        write(*,*) 'Error! Halo is too wide! Abort!'
        stop
    end if

    rnum = 0
    snum = 0

    do local_ind = partition%ts, partition%te

        local_panel = partition%panel_map(local_ind)

        do remote_ind = 1, partition%Nt

            if (remote_ind == local_ind) cycle

            halo_width_x = halo_width
            halo_width_y = halo_width

            !recv part
            local_tile => partition%tiles_y%tile(local_ind)

            remote_panel = partition%panel_map(remote_ind)

            call topology%find_basis_orientation(remote_panel, local_panel, &
                                        send_i_step(snum+1), send_j_step(snum+1),&
                                        first_dim_index(snum+1))

            if (first_dim_index(snum+1) == 'i') then
            !recv from v points to v points
                remote_tile => partition%tiles_y%tile(remote_ind)
            else if (first_dim_index(snum+1) == 'j') then
            !recv from u points to v points
                remote_tile => partition%tiles_x%tile(remote_ind)
            end if

            call topology%transform_tile_coords(remote_panel, remote_tile,     &
                                                local_panel, temp_tile,        &
                                                partition%tiles_y%Nx, partition%tiles_y%Ny)

            !recv exchange. +1 to y halo_width to account for edge points
            if (remote_panel /= local_panel) halo_width_y = halo_width+1

            call find_tiles_halo_intersection(local_tile, halo_width_x, halo_width_y, &
                            temp_tile, halo_type, recv_tile(rnum+1), is_intersection)

            if (is_intersection) then
                rnum = rnum + 1
                recv_to_tile_ind(rnum) = local_ind
                recv_pts_num(rnum) = recv_tile(rnum)%get_points_num()
                recv_from_proc_id(rnum) = partition%proc_map(remote_ind)
                recv_tag(rnum) = partition%Nt*(remote_ind-1) + local_ind
            end if
            !end of recv part

            halo_width_x = halo_width
            halo_width_y = halo_width

            if (first_dim_index(snum+1) == 'i') then
                !send from v points to v points!local/remote tiles already assigned and transformed
                !send exchange. +1 to y halo_width to account for edge points
                if (remote_panel /= local_panel) halo_width_y = halo_width+1
                call find_tiles_halo_intersection(temp_tile, halo_width_x, halo_width_y, &
                              local_tile, halo_type, send_tile(snum+1), is_intersection)
            else if (first_dim_index(snum+1) == 'j') then
                !send from u points to v points
                local_tile  => partition%tiles_x%tile(local_ind)
                remote_tile => partition%tiles_y%tile(remote_ind)
                call topology%transform_tile_coords(remote_panel, remote_tile,   &
                                                    local_panel, temp_tile,     &
                                                    partition%tiles_x%Nx, partition%tiles_x%Ny)
                !send exchange. +1 to x halo_width to account for edge points
                if (remote_panel /= local_panel) halo_width_x = halo_width+1
                call find_tiles_halo_intersection(temp_tile, halo_width_x, halo_width_y, &
                              local_tile, halo_type, send_tile(snum+1), is_intersection)
            end if

            if (is_intersection) then
                snum = snum + 1
                send_from_tile_ind(snum) = local_ind
                send_pts_num(snum) = send_tile(snum)%get_points_num()
                send_to_proc_id(snum) = partition%proc_map(remote_ind)
                send_tag(snum) = partition%Nt*(local_ind-1) + remote_ind
            end if
        end do
    end do

    exchange%send_tile          = send_tile(1:snum)
    exchange%send_points_num    = send_pts_num(1:snum)
    exchange%send_from_tile_ind = send_from_tile_ind(1:snum)
    exchange%send_to_proc_id    = send_to_proc_id(1:snum)
    exchange%send_tag           = send_tag(1:snum)
    exchange%send_i_step        = send_i_step(1:snum)
    exchange%send_j_step        = send_j_step(1:snum)
    exchange%first_dim_index    = first_dim_index(1:snum)
    exchange%send_number        = snum
    exchange%recv_tile          = recv_tile(1:rnum)
    exchange%recv_to_tile_ind   = recv_to_tile_ind(1:rnum)
    exchange%recv_from_proc_id  = recv_from_proc_id(1:rnum)
    exchange%recv_points_num    = recv_pts_num(1:rnum)
    exchange%recv_tag           = recv_tag(1:rnum)
    exchange%recv_number        = rnum

    allocate(exchange%send_buff(snum)   , exchange%recv_buff(rnum)   )
    allocate(exchange%mpi_send_req(snum), exchange%mpi_recv_req(rnum))

    do ind = 1, rnum
        call exchange%recv_buff(ind)%init(recv_pts_num(ind))
    end do
    do ind = 1, snum
        call exchange%send_buff(ind)%init(send_pts_num(ind))
    end do

end function create_symm_halo_vec_exchange_V_points
subroutine find_tiles_halo_intersection(tile1, halo_width1_x, halo_width1_y, tile2, halo_type, out_tile, is_intersect)

    type(tile_t),     intent(in) :: tile1, tile2
    integer(kind=4),  intent(in) :: halo_width1_x, halo_width1_y
    character(len=*), intent(in) :: halo_type

    type(tile_t),    intent(out) :: out_tile
    logical,         intent(out) :: is_intersect

    is_intersect = .true.

    if (halo_type=='full') then
        out_tile%ks = tile1%ks
        out_tile%ke = tile1%ke

        out_tile%is = max(tile1%is-halo_width1_x, tile2%is)
        out_tile%ie = min(tile1%ie+halo_width1_x, tile2%ie)

        out_tile%js = max(tile1%js-halo_width1_y, tile2%js)
        out_tile%je = min(tile1%je+halo_width1_y, tile2%je)

        if ((out_tile%is<=out_tile%ie) .and. (out_tile%js<=out_tile%je)) return

    else if (halo_type == 'cross') then

        out_tile%ks = tile1%ks
        out_tile%ke = tile1%ke

        out_tile%is = max(tile1%is-halo_width1_x, tile2%is)
        out_tile%ie = min(tile1%ie+halo_width1_x, tile2%ie)

        out_tile%js = max(tile1%js, tile2%js)
        out_tile%je = min(tile1%je, tile2%je)

        if ((out_tile%is<=out_tile%ie) .and. (out_tile%js<=out_tile%je)) return

        out_tile%is = max(tile1%is, tile2%is)
        out_tile%ie = min(tile1%ie, tile2%ie)

        out_tile%js = max(tile1%js-halo_width1_y, tile2%js)
        out_tile%je = min(tile1%je+halo_width1_y, tile2%je)

        if ((out_tile%is<=out_tile%ie) .and. (out_tile%js<=out_tile%je)) return
    else
        print*, 'Error! Wrong halo_type!'
        stop
    end if

    is_intersect = .false.

end subroutine find_tiles_halo_intersection

subroutine create_gather_exchange(exchange, points_type, parcomm, partition, master_id)

    use exchange_abstract_mod, only : exchange_t
    use exchange_gather_mod,   only : exchange_gather_t
    use parcomm_mod,           only : parcomm_t

    class(exchange_t), allocatable, intent(out) :: exchange
    character(len=*),               intent(in)  :: points_type
    type(parcomm_t),                intent(in)  :: parcomm
    type(partition_t),              intent(in)  :: partition
    integer(kind=4),   optional,    intent(in)  :: master_id

    type(exchange_gather_t), allocatable :: gather_exchange

    type(tiles_t) :: tiles

    integer(kind=4),  dimension(partition%Nt) :: &
                      recv_is, recv_ie, recv_js, recv_je, recv_ks, recv_ke, &
                      send_is, send_ie, send_js, send_je, send_ks, send_ke

    integer(kind=4),  dimension(partition%Nt) :: &
                      recv_from_proc_id , recv_to_tile_ind , &
                      send_from_tile_ind,  send_tag, recv_tag

    integer(kind=4),  dimension(partition%Nt) :: &
                      send_pts_num, recv_pts_num

    integer(kind=4) :: t, send_exch_num, recv_exch_num, ts, te
    integer(kind=4) :: master_id_loc

    master_id_loc = 0
    if (present(master_id)) master_id_loc = master_id

    recv_exch_num = 0
    send_exch_num = 0

    call partition%get_tiles(points_type, tiles)

    if (parcomm%myid == master_id_loc) then

        do t = 1, partition%Nt

            if (partition%proc_map(t) == parcomm%myid) cycle

            recv_exch_num = recv_exch_num + 1

            recv_is(recv_exch_num) = tiles%tile(t)%is
            recv_ie(recv_exch_num) = tiles%tile(t)%ie

            recv_js(recv_exch_num) = tiles%tile(t)%js
            recv_je(recv_exch_num) = tiles%tile(t)%je

            recv_ks(recv_exch_num) = tiles%tile(t)%ks
            recv_ke(recv_exch_num) = tiles%tile(t)%ke

            recv_pts_num(recv_exch_num) = (recv_ie(recv_exch_num) - recv_is(recv_exch_num) + 1)* &
                                          (recv_je(recv_exch_num) - recv_js(recv_exch_num) + 1)* &
                                          (recv_ke(recv_exch_num) - recv_ks(recv_exch_num) + 1)

            recv_to_tile_ind  (recv_exch_num) = t
            recv_from_proc_id(recv_exch_num)  = partition%proc_map(t)
            recv_tag(recv_exch_num) = t

        end do

    else

        do t = partition%ts, partition%te

            send_exch_num = send_exch_num + 1

            send_is(send_exch_num) = tiles%tile(t)%is
            send_ie(send_exch_num) = tiles%tile(t)%ie

            send_js(send_exch_num) = tiles%tile(t)%js
            send_je(send_exch_num) = tiles%tile(t)%je

            send_ks(send_exch_num) = tiles%tile(t)%ks
            send_ke(send_exch_num) = tiles%tile(t)%ke

            send_pts_num(send_exch_num) = (send_ie(send_exch_num) - send_is(send_exch_num) + 1)* &
                                          (send_je(send_exch_num) - send_js(send_exch_num) + 1)* &
                                          (send_ke(send_exch_num) - send_ks(send_exch_num) + 1)

            send_from_tile_ind(send_exch_num) = t
            send_tag(send_exch_num) = t

        end do
    end if

    allocate(gather_exchange,  source = exchange_gather_t(            &
            send_is            = send_is(1:send_exch_num),            &
            send_ie            = send_ie(1:send_exch_num),            &
            send_js            = send_js(1:send_exch_num),            &
            send_je            = send_je(1:send_exch_num),            &
            send_ks            = send_ks(1:send_exch_num),            &
            send_ke            = send_ke(1:send_exch_num),            &
            recv_is            = recv_is(1:recv_exch_num),            &
            recv_ie            = recv_ie(1:recv_exch_num),            &
            recv_js            = recv_js(1:recv_exch_num),            &
            recv_je            = recv_je(1:recv_exch_num),            &
            recv_ks            = recv_ks(1:recv_exch_num),            &
            recv_ke            = recv_ke(1:recv_exch_num),            &
            recv_to_tile_ind   = recv_to_tile_ind(1:recv_exch_num),   &
            send_from_tile_ind = send_from_tile_ind(1:send_exch_num), &
            recv_from_proc_id  = recv_from_proc_id(1:recv_exch_num),  &
            send_points_num    = send_pts_num(1:send_exch_num),       &
            recv_points_num    = recv_pts_num(1:recv_exch_num),       &
            send_tag           = send_tag(1:send_exch_num),           &
            recv_tag           = recv_tag(1:recv_exch_num),           &
            send_number        = send_exch_num,                       &
            recv_number        = recv_exch_num,                       &
            master_id          = master_id))

    allocate(gather_exchange%send_buff(1:send_exch_num)   , gather_exchange%recv_buff(1:recv_exch_num)   )
    allocate(gather_exchange%mpi_send_req(1:send_exch_num), gather_exchange%mpi_recv_req(1:recv_exch_num))
    do t = 1, recv_exch_num
        call gather_exchange%recv_buff(t)%init(recv_pts_num(t))
    end do
    do t = 1, send_exch_num
        call gather_exchange%send_buff(t)%init(send_pts_num(t))
    end do

    call move_alloc(gather_exchange, exchange)

end subroutine create_gather_exchange

end module exchange_factory_mod

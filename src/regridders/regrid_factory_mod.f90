module regrid_factory_mod

use parcomm_mod,            only : parcomm_global

implicit none

contains

subroutine create_latlon_regrid(regrid_out, domain, Nlon, Nlat,interp_type, &
                                scalar_grid_type)

    use abstract_regrid_mod,    only : regrid_t
    use latlon_regrid_mod,      only : latlon_regrid_t, latlon_interp_t
    use domain_mod,             only : domain_t
    use halo_factory_mod,       only : create_halo_procedure
    use grid_field_factory_mod, only : create_grid_field
    use tile_mod,               only : tile_t

    class(regrid_t), allocatable, intent(out) :: regrid_out
    type(domain_t),               intent(in)  :: domain
    integer(kind=4),              intent(in)  :: Nlon, Nlat
    character(len=*),             intent(in)  :: interp_type
    character(len=*),             intent(in)  :: scalar_grid_type

    class(latlon_regrid_t), allocatable :: regrid
    type(tile_t) :: stencil_bounds
    integer(kind=4), parameter :: A_halo_width = 2, A_ex_halo_width = 8

    allocate(regrid)

    regrid%Nlon = Nlon
    regrid%Nlat = Nlat
    regrid%scalar_grid_type = scalar_grid_type

    !Check:
    if(domain%partition%ts /= 1 .or. &
       domain%partition%te /= domain%partition%num_tiles*domain%partition%num_panels) then
       call parcomm_global%abort("cannot create latlon interpolator for distributed domains")
    end if

    if(regrid%scalar_grid_type == "Ah") then
        call stencil_bounds%init(is=1,ie=domain%mesh_o%tile(1)%nx+1,&
                                 js=1,je=domain%mesh_o%tile(1)%ny+1,&
                                 ks=domain%mesh_o%tile(1)%ks,ke=domain%mesh_o%tile(1)%ke)
        call create_latlon_interp(regrid%interp_scalar,domain,domain%mesh_xy, &
                                  stencil_bounds, Nlat, Nlon, interp_type)
    else if(regrid%scalar_grid_type == "A") then
        call create_halo_procedure(regrid%scalar_halo,domain,A_halo_width,"ECS_O")
        call create_grid_field(regrid%work_field,A_ex_halo_width,0,domain%mesh_o)
        call stencil_bounds%init(is=1-A_halo_width,ie=domain%mesh_o%tile(1)%nx+A_halo_width,&
                                 js=1-A_halo_width,je=domain%mesh_o%tile(1)%ny+A_halo_width,&
                                 ks=domain%mesh_o%tile(1)%ks,ke=domain%mesh_o%tile(1)%ke)
        call create_latlon_interp(regrid%interp_scalar,domain,domain%mesh_o, &
                                  stencil_bounds, Nlat, Nlon, interp_type)
    else
       call parcomm_global%abort(scalar_grid_type//&
                                 " scalar grid type is not supported in latlon regrid,"//&
                                 " use A or Ah")
    end if

    call move_alloc(regrid,regrid_out)
end subroutine create_latlon_regrid

subroutine create_latlon_vector_regrid(regrid_out, domain, Nlon, Nlat, interp_type, &
                                       vector_grid_type,components_type)

    use abstract_regrid_mod,    only : regrid_vec_t
    use latlon_regrid_mod,      only : latlon_regrid_vec_t, latlon_interp_t
    use domain_mod,             only : domain_t
    use halo_factory_mod,       only : create_vector_halo_procedure
    use grid_field_factory_mod, only : create_grid_field
    use tile_mod,               only : tile_t

    use const_mod,              only : pi

    class(regrid_vec_t), allocatable, intent(out) :: regrid_out
    type(domain_t),                   intent(in)  :: domain
    integer(kind=4),                  intent(in)  :: Nlon, Nlat
    character(len=*),                 intent(in)  :: interp_type
    character(len=*),                 intent(in)  :: vector_grid_type
    character(len=*),                 intent(in)  :: components_type

    class(latlon_regrid_vec_t), allocatable :: regrid
    type(tile_t) :: stencil_bounds
    integer(kind=4), parameter :: A_halo_width = 2, A_ex_halo_width = 8
    real(kind=8) :: hlon, hlat, lon, lat
    real(kind=8) :: r(3), iv(3), jv(3), a1(4), a2(4), b1(3), b2(3)
    real(kind=8) :: alpha, beta
    integer(kind=4) :: i, j, panel_ind

    !Check:
    if(domain%partition%ts /= 1 .or. &
       domain%partition%te /= domain%partition%num_tiles*domain%partition%num_panels) then
       call parcomm_global%abort("cannot create latlon interpolator for distributed domains")
    end if

    allocate(regrid)

    regrid%Nlon = Nlon
    regrid%Nlat = Nlat
    regrid%vector_grid_type = vector_grid_type
    regrid%components_type = components_type
    regrid%ks = domain%mesh_o%tile(domain%mesh_o%ts)%ks
    regrid%ke = domain%mesh_o%tile(domain%mesh_o%ts)%ke

    allocate(regrid%u2u(Nlon,Nlat),regrid%u2v(Nlon,Nlat),&
             regrid%v2u(Nlon,Nlat),regrid%v2v(Nlon,Nlat))

    if(regrid%vector_grid_type == "Ah") then
        call stencil_bounds%init(is=1,ie=domain%mesh_o%tile(1)%nx+1,&
                                 js=1,je=domain%mesh_o%tile(1)%ny+1,&
                                 ks=domain%mesh_o%tile(1)%ks,ke=domain%mesh_o%tile(1)%ke)
        call create_latlon_interp(regrid%interp_components,domain,domain%mesh_xy, &
                                  stencil_bounds, Nlat, Nlon, interp_type)
    else if(regrid%vector_grid_type == "A" .or. &
            regrid%vector_grid_type == "C") then
        call create_vector_halo_procedure(regrid%halo,domain,A_halo_width,"ecs_A_vec")
        call create_grid_field(regrid%work_u,A_ex_halo_width,0,domain%mesh_o)
        call create_grid_field(regrid%work_v,A_ex_halo_width,0,domain%mesh_o)
        call stencil_bounds%init(is=1-A_halo_width,ie=domain%mesh_o%tile(1)%nx+A_halo_width,&
                                 js=1-A_halo_width,je=domain%mesh_o%tile(1)%ny+A_halo_width,&
                                 ks=domain%mesh_o%tile(1)%ks,ke=domain%mesh_o%tile(1)%ke)
        call create_latlon_interp(regrid%interp_components,domain,domain%mesh_o, &
                                  stencil_bounds, Nlat, Nlon, interp_type)
    else
       call parcomm_global%abort(vector_grid_type//&
                                 " scalar grid type is not supported in latlon vector regrid,"//&
                                 " use A, Ah or C")
    end if

    hlon = 2._8*pi / real(Nlon,8)
    hlat = pi / (Nlat-1.0_8)

    !Init transform from curvilinear grid components to
    if(components_type == "contravariant" .or. &
       vector_grid_type == "A" .or. vector_grid_type == "C") then
       !Only contravariant halo-procedures are currently implemented for A and C grids
       !Thus vector components will be (anyway) transformed to contravariant
        do j = 1,Nlat
            lat = -0.5_8*pi+hlat*(j-1)
            do i = 1,Nlon
                lon = hlon*(i-1)

                r(1:3)  = [cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)]
                iv(1:3) = [-sin(lon),cos(lon),0.0_8]
                jv(1:3) = [-sin(lat)*cos(lon),-sin(lat)*sin(lon),cos(lat)]
                call domain%metric%transform_cartesian_to_native(panel_ind, alpha, beta, r)
                a1(1:4) = domain%metric%calculate_a1(panel_ind,alpha,beta)
                a2(1:4) = domain%metric%calculate_a2(panel_ind,alpha,beta)
                regrid%u2u(i,j) = sum(iv(1:3)*a1(1:3))
                regrid%u2v(i,j) = sum(jv(1:3)*a1(1:3))
                regrid%v2u(i,j) = sum(iv(1:3)*a2(1:3))
                regrid%v2v(i,j) = sum(jv(1:3)*a2(1:3))
            end do
        end do
    else if(components_type == "covariant") then
        do j = 1,Nlat
            lat = -0.5_8*pi+hlat*(j-1)
            do i = 1,Nlon
                lon = hlon*(i-1)

                r(1:3)  = [cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)]
                iv(1:3) = [-sin(lon),cos(lon),0.0_8]
                jv(1:3) = [-sin(lat)*cos(lon),-sin(lat)*sin(lon),cos(lat)]
                call domain%metric%transform_cartesian_to_native(panel_ind, alpha, beta, r)
                b1(1:3) = domain%metric%calculate_b1(panel_ind,alpha,beta)
                b2(1:3) = domain%metric%calculate_b2(panel_ind,alpha,beta)
                regrid%u2u(i,j) = sum(iv(1:3)*b1(1:3))
                regrid%u2v(i,j) = sum(jv(1:3)*b1(1:3))
                regrid%v2u(i,j) = sum(iv(1:3)*b2(1:3))
                regrid%v2v(i,j) = sum(jv(1:3)*b2(1:3))
            end do
        end do
    else
        call parcomm_global%abort("regrid_factory_mod, create_latlon_vector_regrid" // &
                                  " - unknown vector components type:" // components_type)
    end if
    call move_alloc(regrid,regrid_out)
end subroutine create_latlon_vector_regrid

subroutine create_latlon_interp(interp,domain,mesh,stencil_bounds,Nlat,Nlon,interp_type)
    use latlon_regrid_mod,      only : latlon_interp_t
    use domain_mod,             only : domain_t
    use tile_mod,               only : tile_t
    use mesh_mod,               only : mesh_t

    use const_mod,              only : pi

    type(latlon_interp_t),  intent(out) :: interp
    type(domain_t),         intent(in)  :: domain
    type(mesh_t),           intent(in)  :: mesh
    type(tile_t),           intent(in)  :: stencil_bounds
    integer(kind=4),        intent(in)  :: Nlat, Nlon
    character(len=*),       intent(in)  :: interp_type

    real(kind=8)    :: lon, lat, hlon, hlat, r(3)
    real(kind=8)    :: alpha, beta
    integer(kind=4) :: panel_ind, i, j, ipanel, t, nx, ny
    integer(kind=4) :: iright
    real(kind=8)    :: hx, hy, alpha0, beta0, shift_alpha, shift_beta
    real(kind=8)    :: zdx, zdy

    interp%Nlon = Nlon
    interp%Nlat = Nlat
    interp%interp_type = interp_type
    interp%ks = mesh%tile(mesh%ts)%ks
    interp%ke = mesh%tile(mesh%ts)%ke

    select case(interp_type)
    case('linear')
        interp%stencil_width = 2
    case('cubic')
        interp%stencil_width = 4
    case default
        call parcomm_global%abort("regrid_factory_mod, create_latlon_interp"//&
                                  " - unknown interpolation type: "//interp_type //&
                                  ", use linear or cubic")
    end select

    allocate(interp%wx(interp%stencil_width,Nlon,Nlat))
    allocate(interp%wy(interp%stencil_width,Nlon,Nlat))
    allocate(interp%tile_ind(Nlon,Nlat))
    allocate(interp%x_ind(Nlon,Nlat),interp%y_ind(Nlon,Nlat))

    hx = mesh%tile(mesh%ts)%hx
    hy = mesh%tile(mesh%ts)%hy
    alpha0 = mesh%tile(mesh%ts)%alpha_0
    beta0  = mesh%tile(mesh%ts)%beta_0
    shift_alpha = mesh%tile(mesh%ts)%shift_i
    shift_beta  = mesh%tile(mesh%ts)%shift_j
    nx = mesh%tile(mesh%ts)%nx
    ny = mesh%tile(mesh%ts)%ny

    hlon = 2._8*pi / real(Nlon,8)
    hlat = pi / (Nlat-1.0_8)

    do j = 1,Nlat
        lat = -0.5_8*pi+hlat*(j-1)
        do i = 1,Nlon
            lon = hlon*(i-1)
            r(1) = cos(lat)*cos(lon)
            r(2) = cos(lat)*sin(lon)
            r(3) = sin(lat)

            call domain%metric%transform_cartesian_to_native(panel_ind, alpha, beta, r)

            zdx = (alpha-alpha0)/hx+1-shift_alpha
            zdy = (beta -beta0) /hx+1-shift_beta
            !find leftmost point of stencil that lies in stencil bounds
            !we assume that (stencil_bounds%ie-stencil_bounds%is+1) > stencil_width
            !i.e. domain is wide enough
            iright = min(floor(zdx)+interp%stencil_width/2,stencil_bounds%ie)
            interp%x_ind(i,j) = max(iright-interp%stencil_width+1,stencil_bounds%is)
            iright = min(floor(zdy)+interp%stencil_width/2,stencil_bounds%je)
            interp%y_ind(i,j) = max(iright-interp%stencil_width+1,stencil_bounds%js)
            !save normalized displacement of point for future computation of interpolation weights
            interp%wx(1,i,j)  = zdx - interp%x_ind(i,j)
            interp%wy(1,i,j)  = zdy - interp%y_ind(i,j)
            !find tile
            do t = domain%partition%ts, domain%partition%te
                if(domain%partition%panel_map(t) == panel_ind   .and.&
                   mesh%tile(t)%is <= max(interp%x_ind(i,j),1)  .and.&
                   mesh%tile(t)%ie >= min(interp%x_ind(i,j),nx) .and.&
                   mesh%tile(t)%js <= max(interp%y_ind(i,j),1)  .and.&
                   mesh%tile(t)%je >= min(interp%y_ind(i,j),ny))  then
                   interp%tile_ind(i,j) = t
               end if
            end do
        end do
    end do

    if(interp_type == 'linear') then
        do j=1, Nlat
            do i=1, Nlon
                zdx = interp%wx(1,i,j)
                interp%wx(1,i,j) = 1._8-zdx
                interp%wx(2,i,j) = zdx
                zdy = interp%wy(1,i,j)
                interp%wy(1,i,j) = 1._8-zdy
                interp%wy(2,i,j) = zdy
            end do
        end do
    else if(interp_type == 'cubic') then
        do j=1, Nlat
            do i=1, Nlon
                zdx = interp%wx(1,i,j)
                zdy = interp%wy(1,i,j)
                interp%wx(1,i,j) =-(zdx-1._8)*(zdx-2._8)*(zdx-3._8) / 6._8
                interp%wx(2,i,j) = zdx*(zdx-2._8)*(zdx-3._8) / 2._8
                interp%wx(3,i,j) =-zdx*(zdx-1._8)*(zdx-3._8) / 2._8
                interp%wx(4,i,j) = zdx*(zdx-1._8)*(zdx-2._8) / 6._8
                interp%wy(1,i,j) =-(zdy-1._8)*(zdy-2._8)*(zdy-3._8) / 6._8
                interp%wy(2,i,j) = zdy*(zdy-2._8)*(zdy-3._8) / 2._8
                interp%wy(3,i,j) =-zdy*(zdy-1._8)*(zdy-3._8) / 2._8
                interp%wy(4,i,j) = zdy*(zdy-1._8)*(zdy-2._8) / 6._8
            end do
        end do
    end if
end

end module regrid_factory_mod

module sbp_operator_mod

use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global

implicit none

!Implements application of 1D SBP operator (difference or interpolation)

type sbp_operator_t
    real(kind=8), allocatable    :: W_edge_l(:,:) !operator weights at left edge
    real(kind=8), allocatable    :: W_edge_r(:,:) !operator weights at right edge
    integer(kind=4), allocatable :: edge_last_l(:)  !last nonzero weight in W_edge_l
    integer(kind=4), allocatable :: edge_first_shift_r(:) !first nonzero weight in W_edge_r (minus n)
    real(kind=8), allocatable    :: W_in(:)         !inner stencil weights
    integer(kind=4)              :: in_shift        !shift of inner stencil first element with respect to
                                                    !the index of calculated element
    real(kind=8), allocatable    :: Al_out(:), Al_in(:) !Mass matrices at output and input grids left side
    real(kind=8), allocatable    :: Ar_out(:), Ar_in(:) !Mass matrices at output and input grids right side
    integer(kind=4)              :: dnx  ! difference of the dimensions of output and input grids
                                         ! in the direction of operator application
                                         ! nx_output - nx_input
    real(kind=8), allocatable    :: proj_operator_l(:), proj_operator_r(:)

    contains
    !Grid field in -> array out
    procedure, public :: apply_gf_to_array      => apply_sbp_gf_to_array
    !grid field in (+tile) -> grid field out
    procedure, public :: apply_gf_to_gf         => apply_sbp_gf_to_gf
    !array in -> array out
    procedure, public :: apply_array_to_array   => apply_sbp_array_to_array
    generic :: apply => apply_gf_to_array, apply_array_to_array, apply_gf_to_gf

    procedure, public:: apply_z_grid_field => apply_sbp_z_gridfield
    procedure, public:: apply_z_tile_field => apply_sbp_z_tilefield
    generic :: apply_z => apply_z_grid_field, apply_z_tile_field

    procedure, public :: add_penalty_gf
    procedure, public :: add_penalty_array
    generic :: add_penalty => add_penalty_gf, add_penalty_array
end type sbp_operator_t

contains

subroutine apply_sbp_gf_to_gf(this, f_out, work_tile, &
                                                nx_out_grid, direction, f_in)
    class(sbp_operator_t), intent(in) :: this
    !work_tile = indices where operator action will be calculated
    type(tile_t),          intent(in)  :: work_tile
    integer(kind=4),       intent(in)  :: nx_out_grid
    character(len=*),      intent(in)  :: direction
    type(tile_field_t),    intent(in)  :: f_in
    !output:
    type(tile_field_t),    intent(inout) :: f_out


    integer(kind=4) :: nx_in_grid, k, f_is, f_ie, f_js, f_je

    nx_in_grid = nx_out_grid-this%dnx

    do k=work_tile%ks, work_tile%ke
        f_is = f_in%is; f_ie = f_in%ie
        f_js = f_in%js; f_je = f_in%je

        call sbp_apply(f_out%p(f_out%is:f_out%ie,f_out%js:f_out%je,k),&
                       f_out%is,  f_out%ie,  f_out%js,  f_out%je,     &
                       work_tile%is, work_tile%ie, work_tile%js, work_tile%je,    &
                       f_in%p(f_is:f_ie,f_js:f_je,k), f_is, f_ie, f_js, f_je,     &
                       nx_in_grid, nx_out_grid,                                   &
                       this%W_edge_l, this%edge_last_l,                           &
                       this%W_edge_r, this%edge_first_shift_r, &
                       this%W_in, this%in_shift, direction)
    end do

end subroutine apply_sbp_gf_to_gf

subroutine apply_sbp_gf_to_array(this, a_out, work_tile, out_tile, &
                                                nx_out_grid, direction, f_in)
    class(sbp_operator_t), intent(in) :: this
    !work_tile = indices where operator action will be calculated
    !bounding_tile = output arrat bounds
    type(tile_t),          intent(in) :: work_tile, out_tile
    integer(kind=4),       intent(in) :: nx_out_grid
    character(len=*),      intent(in) :: direction
    type(tile_field_t),    intent(in) :: f_in
    !output:
    real(kind=8), intent(out) :: a_out(out_tile%is:out_tile%ie, &
                                       out_tile%js:out_tile%je, &
                                       out_tile%ks:out_tile%ke)


    integer(kind=4) :: nx_in_grid, k, f_is, f_ie, f_js, f_je

    nx_in_grid = nx_out_grid-this%dnx

    do k=work_tile%ks, work_tile%ke
        f_is = f_in%is; f_ie = f_in%ie
        f_js = f_in%js; f_je = f_in%je

        call sbp_apply(a_out(out_tile%is:out_tile%ie, out_tile%js:out_tile%je, k),&
                       out_tile%is,  out_tile%ie,  out_tile%js,  out_tile%je,     &
                       work_tile%is, work_tile%ie, work_tile%js, work_tile%je,    &
                       f_in%p(f_is:f_ie,f_js:f_je,k), f_is, f_ie, f_js, f_je,     &
                       nx_in_grid, nx_out_grid,                                   &
                       this%W_edge_l, this%edge_last_l,                           &
                       this%W_edge_r, this%edge_first_shift_r, &
                       this%W_in, this%in_shift, direction)
    end do

end subroutine apply_sbp_gf_to_array

subroutine apply_sbp_array_to_array(this, a_out, work_tile, out_tile, &
                                      nx_out_grid, direction, a_in, in_tile)
    class(sbp_operator_t), intent(in) :: this
    !work_tile = indices where operator action will be calculated
    !bounding_tile = output arrat bounds
    type(tile_t),          intent(in) :: work_tile, out_tile
    integer(kind=4),       intent(in) :: nx_out_grid
    character(len=*),      intent(in) :: direction
    type(tile_t),          intent(in) :: in_tile
    real(kind=8),          intent(in) :: a_in(in_tile%is:in_tile%ie, &
                                              in_tile%js:in_tile%je, &
                                              in_tile%ks:in_tile%ke)
    !output:
    real(kind=8), intent(out) :: a_out(out_tile%is:out_tile%ie, &
                                       out_tile%js:out_tile%je, &
                                       out_tile%ks:out_tile%ke)


    integer(kind=4) :: nx_in_grid, k, f_is, f_ie, f_js, f_je

    nx_in_grid = nx_out_grid-this%dnx

    do k=work_tile%ks, work_tile%ke
        f_is = in_tile%is; f_ie = in_tile%ie
        f_js = in_tile%js; f_je = in_tile%je

        call sbp_apply(a_out(out_tile%is:out_tile%ie, out_tile%js:out_tile%je, k),&
                       out_tile%is,  out_tile%ie,  out_tile%js,  out_tile%je,     &
                       work_tile%is, work_tile%ie, work_tile%js, work_tile%je,    &
                       a_in(f_is:f_ie,f_js:f_je,k), f_is, f_ie, f_js, f_je,     &
                       nx_in_grid, nx_out_grid,                                   &
                       this%W_edge_l, this%edge_last_l,                           &
                       this%W_edge_r, this%edge_first_shift_r, &
                       this%W_in, this%in_shift, direction)
    end do

end subroutine apply_sbp_array_to_array

subroutine sbp_apply(df,is1,ie1,js1,je1,is,ie,js,je,                          &
                     f,is0,ie0,js0,je0,nsrc,ntarg,                            &
                     Ql,l_last_nonzero,Qr,r_first_nonzero_shift,Qin,inshift,direction)
    !dimensions of input function
    real(kind=8),     intent(in)    :: Ql(:,:)
    integer(kind=4),  intent(in)    :: l_last_nonzero(:)
    real(kind=8),     intent(in)    :: Qr(:,:)
    integer(kind=4),  intent(in)    :: r_first_nonzero_shift(:)
    real(kind=8),     intent(in)    :: Qin(:)
    integer(kind=4),  intent(in)    :: inshift
    integer(kind=4),  intent(in)    :: is0, ie0, js0, je0
    !global mesh dimension in x:
    integer(kind=4),  intent(in)    :: nsrc, ntarg
    !output function dimensions:
    integer(kind=4),  intent(in)    :: is1,ie1,js1,je1
    integer(kind=4),  intent(in)    :: is, ie, js, je

    real(kind=8),     intent(in)    :: f(is0:ie0,js0:je0)
    character(len=*), intent(in)    :: direction

    real(kind=8),     intent(inout) :: df(is1:ie1,js1:je1)

    integer(kind=4) :: i, ir, j, jr
    integer(kind=4) :: mx, my, nwst, nwr
    integer(kind=4) :: l_edge_size, r_edge_size, inwidth, inend

    l_edge_size = size(Ql,2)
    r_edge_size = size(Qr,2)
    nwr = size(Qr,1)
    inwidth = size(Qin)
    inend = inwidth+inshift-1

    select case(direction)
    case('x')

    do j=js,je
        do i = max(is,1),min(l_edge_size,ie)
            mx = l_last_nonzero(i)
            df(i,j) = sum(Ql(1:mx,i)*f(1:mx,j))
        end do
        do i = max(l_edge_size+1,is), min(ie,ntarg-r_edge_size)
            df(i,j) = sum(Qin(1:inwidth)*f(i+inshift:i+inend,j))
        end do
        do i=max(is,ntarg-r_edge_size+1),min(ie,ntarg)
            ir = i - (ntarg-r_edge_size)
            mx = nsrc+r_first_nonzero_shift(ir)
            nwst = nwr+r_first_nonzero_shift(ir)
            df(i,j) = sum(Qr(nwst:nwr,ir)*f(mx:nsrc,j))
        end do
    end do

    case('y')

    do j = max(js,1),min(l_edge_size,je)
        do i=is,ie
            my = l_last_nonzero(j)
            df(i,j) = sum(Ql(1:my,j)*f(i,1:my))
        end do
    end do
    do j = max(l_edge_size+1,js), min(je,ntarg-r_edge_size)
        do i=is,ie
            df(i,j) = sum(Qin(1:inwidth)*f(i,j+inshift:j+inend))
        end do
    end do
    do j=max(js,ntarg-r_edge_size+1),min(je,ntarg)
        do i=is,ie
            jr = j - (ntarg-r_edge_size)
            my = nsrc+r_first_nonzero_shift(jr)
            nwst = nwr+r_first_nonzero_shift(jr)
            df(i,j) = sum(Qr(nwst:nwr,jr)*f(i,my:nsrc))
        end do
    end do

    case default
    call parcomm_global%abort("unknown direction in sbp_apply: "//direction)
    end select
end subroutine sbp_apply

subroutine apply_sbp_z_gridfield(this, f_out, f_in, mesh_out)
    class(sbp_operator_t), intent(in)  :: this
    !work_tile = indices where operator action will be calculated
    type(grid_field_t),    intent(in)  :: f_in
    type(mesh_t),          intent(in)  :: mesh_out
    !output:
    type(grid_field_t),    intent(inout) :: f_out

    integer(kind=4) :: t

    do t = mesh_out%ts, mesh_out%te
        call apply_sbp_z_tilefield(this,f_out%tile(t), f_in%tile(t), mesh_out%tile(t))
    end do
end subroutine apply_sbp_z_gridfield

subroutine apply_sbp_z_tilefield(this, f_out, f_in, mesh_out)
    class(sbp_operator_t), intent(in) :: this
    !work_tile = indices where operator action will be calculated
    type(tile_field_t),    intent(in)  :: f_in
    type(tile_mesh_t),     intent(in)  :: mesh_out
    !output:
    type(tile_field_t),    intent(inout) :: f_out

    integer(kind=4) i, j, k, is, ie, js, je, ks, ke, m, kr
    integer(kind=4) :: l_edge_size, r_edge_size, nwr, nwst, inwidth, inend
    integer(kind=4) :: nz_out, nz_in

    l_edge_size = size(this%W_edge_l,2)
    r_edge_size = size(this%W_edge_r,2)
    nwr = size(this%W_edge_r,1)
    inwidth = size(this%W_in)
    inend = inwidth+this%in_shift-1
    nz_out = mesh_out%nz
    nz_in = nz_out - this%dnx

    is = mesh_out%is; ie = mesh_out%ie
    js = mesh_out%js; je = mesh_out%je
    ks = mesh_out%ks; ke = mesh_out%ke

    do k = max(ks,1),min(l_edge_size,ke)
        do j = js, je;  do i = is, ie
            m = this%edge_last_l(k)
            f_out%p(i,j,k) = sum(this%W_edge_l(1:m,k)*f_in%p(i,j,1:m))
        end do; end do
    end do
    do k = max(l_edge_size+1,ks), min(ke,nz_out-r_edge_size)
        do j = js, je;  do i = is, ie
            f_out%p(i,j,k) = sum(this%W_in(1:inwidth)*f_in%p(i,j,k+this%in_shift:k+inend))
        end do; end do
    end do
    do k = max(ks,nz_out-r_edge_size+1),min(ke,nz_out)
        do j = js, je;  do i = is, ie
            kr = k - (nz_out-r_edge_size)
            m  = nz_in+this%edge_first_shift_r(kr)
            nwst = nwr+this%edge_first_shift_r(kr)
            f_out%p(i,j,k) = sum(this%W_edge_r(nwst:nwr,kr)*f_in%p(i,j,m:nz_in))
        end do; end do
    end do

end subroutine apply_sbp_z_tilefield

subroutine add_penalty_gf(this, f_out, work_tile, nx_out_grid, direction, &
                                                             penalty_type,f_in)
    class(sbp_operator_t), intent(in) :: this
    !work_tile = indices where operator action will be calculated
    type(tile_t),          intent(in)  :: work_tile
    integer(kind=4),       intent(in)  :: nx_out_grid
    character(len=*),      intent(in)  :: direction, penalty_type
    type(tile_field_t),    intent(in)  :: f_in
    !output:
    type(tile_field_t),    intent(inout) :: f_out

    type(tile_t) :: out_tile

    out_tile = tile_t(is=f_out%is,ie=f_out%ie,js=f_out%js,je=f_out%je,&
                      ks=f_out%ks,ke=f_out%ke)

    if(penalty_type == "at_center") then
        call add_penalty_center(this, f_out%p, work_tile, out_tile, nx_out_grid, direction, f_in)
    else if(penalty_type == "at_interface") then
        call add_penalty_interface(this, f_out%p, work_tile, out_tile, nx_out_grid, direction, f_in)
    else if(penalty_type == "at_interface_new_temp") then
        call add_penalty_interface_new_temp(this, f_out%p, work_tile, out_tile, nx_out_grid, direction, f_in)
    else
        call parcomm_global%abort("sbp operator mod, add penalty, unknown penalty type:"//&
                                  penalty_type)
    end if

end subroutine add_penalty_gf

subroutine add_penalty_array(this, f_out, work_tile, out_tile, nx_out_grid, direction, &
                                                             penalty_type,f_in)
    class(sbp_operator_t), intent(in) :: this
    !work_tile = indices where operator action will be calculated
    type(tile_t),          intent(in)  :: work_tile, out_tile
    integer(kind=4),       intent(in)  :: nx_out_grid
    character(len=*),      intent(in)  :: direction, penalty_type
    type(tile_field_t),    intent(in)  :: f_in
    !output:
    real(kind=8),          intent(inout) :: f_out(out_tile%is:out_tile%ie,&
                                                  out_tile%js:out_tile%je,&
                                                  out_tile%ks:out_tile%ke)

    if(penalty_type == "at_center") then
        call add_penalty_center(this, f_out, work_tile, out_tile, nx_out_grid, direction, f_in)
    else if(penalty_type == "at_interface") then
        call add_penalty_interface(this, f_out, work_tile, out_tile, nx_out_grid, direction, f_in)
    else
        call parcomm_global%abort("sbp operator mod, add penalty, unknown penalty type:"//&
                                  penalty_type)
    end if

end subroutine add_penalty_array

subroutine add_penalty_center(this, f_out, work_tile, out_tile, nx_out_grid, direction,f_in)
    class(sbp_operator_t), intent(in) :: this
    !work_tile = indices where operator action will be calculated
    type(tile_t),          intent(in)  :: out_tile, work_tile
    integer(kind=4),       intent(in)  :: nx_out_grid
    character(len=*),      intent(in)  :: direction
    type(tile_field_t),    intent(in)  :: f_in
    !output:
    !type(tile_field_t),    intent(inout) :: f_out
    real(kind=8),          intent(inout) :: f_out(out_tile%is:out_tile%ie,&
                                                  out_tile%js:out_tile%je,&
                                                  out_tile%ks:out_tile%ke)

    integer(kind=4) :: nx_in_grid, is, ie, js, je, ks, ke
    integer(kind=4) :: i, j, k, ne_left, ne_right, n_Ar
    real(kind=8)    :: penalty_l(size(this%proj_operator_l))
    real(kind=8)    :: penalty_r(size(this%proj_operator_r))
    real(kind=8)    :: df

    nx_in_grid = nx_out_grid-this%dnx
    ne_left  = size(this%proj_operator_l)
    ne_right = size(this%proj_operator_r)

    is = work_tile%is; ie = work_tile%ie
    js = work_tile%js; je = work_tile%je
    ks = work_tile%ks; ke = work_tile%ke

    penalty_l(1:ne_left) = this%proj_operator_l(1:ne_left) / this%Al_out(1:ne_left)
    n_Ar = size(this%Ar_out)
    penalty_r(1:ne_right) = this%proj_operator_r(1:ne_right) / this%Ar_out(n_Ar-ne_right+1:n_Ar)

    if(direction == "x") then
        do k=ks, ke
            do j = js, je
                if(is <= ne_left) df = f_in%p(1,j,k)-f_in%p(0,j,k)
                do i=max(1,is),min(ie,ne_left)
                    !f_out%p(i,j,k) = f_out%p(i,j,k)+0.5_8*penalty_l(i)*df
                    f_out(i,j,k) = f_out(i,j,k)+0.5_8*penalty_l(i)*df
                end do
                if(ie >= nx_out_grid-ne_right+1) &
                    df = f_in%p(nx_in_grid+1,j,k)-f_in%p(nx_in_grid,j,k)
                do i=max(nx_out_grid-ne_right+1,is),min(ie,nx_out_grid)
                    !f_out%p(i,j,k) = f_out%p(i,j,k)+0.5_8*penalty_r(i-nx_out_grid+ne_right)*df
                    f_out(i,j,k) = f_out(i,j,k)+0.5_8*penalty_r(i-nx_out_grid+ne_right)*df
                end do
            end do
        end do
    else if(direction == "y") then
        do k=ks, ke
            do j=max(1,js), min(ie, ne_left)
                do i=is, ie
                    df = f_in%p(i,1,k)-f_in%p(i,0,k)
                    ! f_out%p(i,j,k) = f_out%p(i,j,k)+0.5_8*penalty_l(j)*df
                    f_out(i,j,k) = f_out(i,j,k)+0.5_8*penalty_l(j)*df
                end do
            end do
            do j=max(nx_out_grid-ne_right+1,js), min(je, nx_out_grid)
                do i=is, ie
                    df = f_in%p(i,nx_in_grid+1,k)-f_in%p(i,nx_in_grid,k)
                    ! f_out%p(i,j,k) = f_out%p(i,j,k)+0.5_8*penalty_r(j-nx_out_grid+ne_right)*df
                    f_out(i,j,k) = f_out(i,j,k)+0.5_8*penalty_r(j-nx_out_grid+ne_right)*df
                end do
            end do
        end do
    else
        call parcomm_global%abort("sbp operator mod, add penalty_centers, unknown direction:"//&
                                  direction)
    end if

end subroutine add_penalty_center

subroutine add_penalty_interface(this, f_out, work_tile, out_tile, nx_out_grid, direction,f_in)
    class(sbp_operator_t), intent(in) :: this
    !work_tile = indices where operator action will be calculated
    type(tile_t),          intent(in)  :: work_tile, out_tile
    integer(kind=4),       intent(in)  :: nx_out_grid
    character(len=*),      intent(in)  :: direction
    type(tile_field_t),    intent(in)  :: f_in
    !output:
    real(kind=8),          intent(inout) :: f_out(out_tile%is:out_tile%ie,&
                                                  out_tile%js:out_tile%je,&
                                                  out_tile%ks:out_tile%ke)

    integer(kind=4) :: nx_in_grid, is, ie, js, je, ks, ke
    integer(kind=4) :: i, j, k, ne_left, ne_right, n_Ar
    real(kind=8)    :: penalty_l(size(this%proj_operator_l))
    real(kind=8)    :: penalty_r(size(this%proj_operator_r))
    real(kind=8)    :: df

    nx_in_grid = nx_out_grid-this%dnx
    ne_left  = size(this%proj_operator_l)
    ne_right = size(this%proj_operator_r)

    is = work_tile%is; ie = work_tile%ie
    js = work_tile%js; je = work_tile%je
    ks = work_tile%ks; ke = work_tile%ke

    penalty_l(1:ne_left) = this%proj_operator_l(1:ne_left) / this%Al_out(1)
    n_Ar = size(this%Ar_out)
    penalty_r(1:ne_right) = this%proj_operator_r(1:ne_right) / this%Ar_out(n_Ar)

    if(direction == "x") then
        do k=ks, ke
            if(is == 1) then; do j = js, je
                df = 0.0_8
                do i=1, ne_left
                    df = df+penalty_l(i)*(f_in%p(i,j,k)-f_in%p(1-i,j,k))
                end do
                f_out(1,j,k) = f_out(1,j,k)+0.5_8*df
            end do; end if
            if(ie == nx_out_grid) then; do j = js, je
                df = 0.0_8
                do i=1, ne_right
                    df = df+penalty_r(ne_right-i+1)*(f_in%p(nx_in_grid+i,j,k)- &
                                                     f_in%p(nx_in_grid-i+1,j,k))
                end do
                f_out(ie,j,k) = f_out(ie,j,k)+0.5_8*df
            end do; end if
        end do
    else if(direction == "y") then
        do k=ks, ke
            if(js == 1) then
                do i=is, ie
                    df = 0.0_8
                    do j=1, ne_left
                        df = df+penalty_l(j)*(f_in%p(i,j,k)-f_in%p(i,1-j,k))
                    end do
                    f_out(i,1,k) = f_out(i,1,k)+0.5_8*df
                end do
            end if
            if(je == nx_out_grid) then
                do i=is, ie
                    df = 0.0_8
                    do j=1, ne_right
                        df = df+penalty_r(ne_right-j+1)*(f_in%p(i,nx_in_grid+j,k)- &
                                                         f_in%p(i,nx_in_grid-j+1,k))
                    end do
                    f_out(i,je,k) = f_out(i,je,k)+0.5_8*df
                end do
            end if
        end do
    else
        call parcomm_global%abort("sbp operator mod, add penalty_centers, unknown direction:"//&
                                  direction)
    end if

end subroutine add_penalty_interface
subroutine add_penalty_interface_new_temp(this, f_out, work_tile, out_tile, nx_out_grid, direction,f_in)
    class(sbp_operator_t), intent(in) :: this
    !work_tile = indices where operator action will be calculated
    type(tile_t),          intent(in)  :: work_tile, out_tile
    integer(kind=4),       intent(in)  :: nx_out_grid
    character(len=*),      intent(in)  :: direction
    type(tile_field_t),    intent(in)  :: f_in
    !output:
    real(kind=8),          intent(inout) :: f_out(out_tile%is:out_tile%ie,&
                                                  out_tile%js:out_tile%je,&
                                                  out_tile%ks:out_tile%ke)

    integer(kind=4) :: nx_in_grid, is, ie, js, je, ks, ke
    integer(kind=4) :: i, j, k, ne_left, ne_right, n_Ar
    real(kind=8)    :: penalty_l(size(this%proj_operator_l))
    real(kind=8)    :: penalty_r(size(this%proj_operator_r))
    real(kind=8)    :: df

    nx_in_grid = nx_out_grid-this%dnx
    ne_left  = size(this%proj_operator_l)
    ne_right = size(this%proj_operator_r)

    is = work_tile%is; ie = work_tile%ie
    js = work_tile%js; je = work_tile%je
    ks = work_tile%ks; ke = work_tile%ke

    penalty_l(1:ne_left) = this%proj_operator_l(1:ne_left) / this%Al_out(1)
    n_Ar = size(this%Ar_out)
    penalty_r(1:ne_right) = this%proj_operator_r(1:ne_right) / this%Ar_out(n_Ar)

    if(direction == "x") then
        do k=ks, ke
            if(is == 1) then; do j = js, je
                df = 0.0_8
                do i=1, ne_left
                    df = df-penalty_l(i)*(f_in%p(i,j,k)+f_in%p(1-i,j,k))
                end do
                f_out(1,j,k) = f_out(1,j,k)+0.5_8*df
            end do; end if
            if(ie == nx_out_grid) then; do j = js, je
                df = 0.0_8
                do i=1, ne_right
                    df = df-penalty_r(ne_right-i+1)*(f_in%p(nx_in_grid+i,j,k)+ &
                                                     f_in%p(nx_in_grid-i+1,j,k))
                end do
                f_out(ie,j,k) = f_out(ie,j,k)+0.5_8*df
            end do; end if
        end do
    else if(direction == "y") then
        do k=ks, ke
            if(js == 1) then
                do i=is, ie
                    df = 0.0_8
                    do j=1, ne_left
                        df = df-penalty_l(j)*(f_in%p(i,j,k)+f_in%p(i,1-j,k))
                    end do
                    f_out(i,1,k) = f_out(i,1,k)+0.5_8*df
                end do
            end if
            if(je == nx_out_grid) then
                do i=is, ie
                    df = 0.0_8
                    do j=1, ne_right
                        df = df-penalty_r(ne_right-j+1)*(f_in%p(i,nx_in_grid+j,k)+ &
                                                         f_in%p(i,nx_in_grid-j+1,k))
                    end do
                    f_out(i,je,k) = f_out(i,je,k)+0.5_8*df
                end do
            end if
        end do
    else
        call parcomm_global%abort("sbp operator mod, add penalty_centers, unknown direction:"//&
                                  direction)
    end if

end subroutine add_penalty_interface_new_temp
end module sbp_operator_mod

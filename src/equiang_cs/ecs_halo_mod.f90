!Module to interpolate values to virtual points beyond face edge
!Edge interpolations are treated as one following case:
!    __________
!   / source  /
!  / /^b     /
! / +-->a   /
!/_________/     ^z
!|         |     |  /^y
!|  ^b     |     | /
!|  | targ |     |/
!|  +-->a  |     +------->x
!|_________|

!Edge and corners numeration:
!4__________3
!| edge2    |
!|e        e|
!|d        d|
!|g        g|
!|3        4|
!1----------2
!   edge1

!corner procedures are treated as corner N1

module ecs_halo_mod

use halo_mod,          only : halo_t
use exchange_halo_mod, only : exchange_t
use parcomm_mod,       only : parcomm_global

implicit none

!type-container for precomputed interpolation weights to have
!multiple halo-interpolators for multiple grids (for geometric MG-solver)

type, extends(halo_t) :: ecs_halo_t

    integer(kind=4)                    :: ts, te
    class(exchange_t),     allocatable :: exch_halo
    type(ecs_tile_halo_t), allocatable :: tile(:)
    logical :: is_z_interfaces

    contains

    procedure :: get_halo_scalar => get_ecs_halo

end type

type ecs_tile_halo_t
  integer(kind=4) n          !number of real grid points along cubed-sphere face edge
  integer(kind=4) halo_width !number of rows in halo-zone
  integer(kind=4) wsx, wex     !starting and ending index of interpolation weights x-edges
  integer(kind=4) wsy, wey     !starting and ending index of interpolation weights y-edges
  logical lhalo(4), lcorn(4) !flags to perform haloprocedures at specific edge and corner
  real(kind=8)   , allocatable :: wx(:,:,:) !interpolation weights for edges along x-axis
  integer(kind=4), allocatable :: indx(:,:) !interpolation stencil base point index for edges along x-axis
                                           !   s-----s-t----s------s
                                           !   base--^
  real(kind=8)   , allocatable :: wy(:,:,:) !interpolation weights, edges along y-axis
  integer(kind=4), allocatable :: indy(:,:) !interpolation stencil base point index, edgws along x-axis

  contains
  procedure, public :: interp  => ext_ecs_tile_halo

end type ecs_tile_halo_t

contains


subroutine get_ecs_halo(this,f,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(ecs_halo_t),        intent(inout) :: this
    class(grid_field_t),      intent(inout) :: f
    type(domain_t),           intent(in)    :: domain
    integer(kind=4),          intent(in)    :: halo_width

    integer(kind=4) t

    call this%exch_halo%do(f, domain%parcomm)

    do t=this%ts,this%te
        if(this%is_z_interfaces) then
            call this%tile(t)%interp(f%tile(t),domain%partition%tiles_z%tile(t),halo_width)
        else
            call this%tile(t)%interp(f%tile(t),domain%partition%tiles_o%tile(t),halo_width)
        end if
    end do
end subroutine get_ecs_halo

subroutine ext_ecs_tile_halo(this, f, tile, halo_width)
!interpolate source face values at halo zones to target face virtual points
use grid_field_mod, only : tile_field_t
use tile_mod,       only : tile_t

class(ecs_tile_halo_t), intent(in)    :: this
type(tile_field_t),     intent(inout) :: f
type(tile_t),           intent(in)    :: tile
integer(kind=4),        intent(in)    :: halo_width

!local
integer(kind=4) n, thw, ihw
integer(kind=4) is, ie, js, je
integer(kind=4) nvi, nvj, nvk
integer(kind=4) isv, iev, jsv,jev, ksv, kev
integer(kind=4) klev
integer(kind=4) i,j,k
logical lhalo(4) !halo-procedure at edge
logical lcorn(4) !corner halo-procedure  for numeration of edges and corners see below
logical lfail_hw, lfail_corn, lfail_halo_long
real(kind=8) zf_csp(6,f%ks:f%ke,4) !values at corner-points and near-corner points@first halo-row
!local real(kind=8) zbufc(1-f%nvi:f%nvi,1-f%nvj:f%nvj,f%ks-f%nvk:f%ke+f%nvk) !store values for corner procedure

!short names for needed params
n = this%n
thw = this%halo_width !dimensions of initialized weights, theoretically max hw
ihw = halo_width      !width of user-requested halo-interpolations

if(thw < ihw) call parcomm_global%abort("requested halo width is greater than " // &
                                          "initialized maximum halo_width")

is = tile%is;      ie = tile%ie
js = tile%js;      je = tile%je
nvi = tile%is-f%is;    nvj = tile%js-f%js; nvk = tile%ks-f%ks

isv = is-nvi;   iev = ie+nvi
jsv = js-nvj;   jev = je+nvj
ksv = tile%ks-nvk; kev = tile%ke+nvk

klev = kev-ksv+1

lhalo = this%lhalo
lcorn = this%lcorn

if(lcorn(1).or.lcorn(2).or.lcorn(3).or.lcorn(4)) then
    !interpolate f%p to special points in corners, store obtained values in zf_csp
    !zf_csp -> f%p after edge halo procedure
    call ecs_halo_corners(zf_csp,                                       & !special points interpolated values
                          f%p,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                          this%wx,this%wsx,this%wex,thw,                & !weight array and its bounds
                          lcorn,n,ihw)                                  !operating parameters
end if

if(lhalo(1).or.lhalo(2)) then
    call ecs_halo_edges_x(f%p,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                          this%wx,this%wsx,this%wex,thw,this%indx,      & !weight array and its bounds
                          lhalo,n,ihw)                                  !operating parameters
end if
if(lhalo(3).or.lhalo(4)) then
    call ecs_halo_edges_y(f%p,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                          this%wy,this%wsy,this%wey,thw,this%indy,      & !weight array and its bounds
                          lhalo,n,ihw)                                  !operating parameters
end if

call ecs_halo_corners_put(f%p,isv,iev,jsv,jev,klev,zf_csp,n,lcorn,ihw)

end subroutine ext_ecs_tile_halo

!corner procedures for halo-zones
!currently interpolates only left/right values in first halo row
!planned: halo-corner interpolation
subroutine ecs_halo_corners(pfcsp,                                       &
                            pf,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                            w,ws,we,whw,                              & !weight array and its bounds
                            lcorn,n,hw)                                    !global operating parameters

!output
real(kind=8),   intent(out) :: pfcsp(6,klev,4) !values at corner-points and near-corner points in first halo-row
!pfcsp placement for corner#1:
!  ha-|real points
!  lo |
!  __1|____
!   43|2 halo
!   65|
!  halo
!  corner
!input
!function values:
integer(kind=4), intent(in) :: is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj
real(kind=8),    intent(in) :: pf(isv:iev,jsv:jev,klev)
!interpolation weights
integer(kind=4), intent(in) :: ws,we,whw
real(kind=8),    intent(in) :: w(-1:2,ws:we,whw) !interpolation weights
!operting parameters
logical,         intent(in) :: lcorn(4)
integer(kind=4), intent(in) :: n, hw
!local
real(kind=8)    zbufc(1-nvi:nvi,1-nvj:nvj,klev) !store values for corner procedure
integer(kind=4) k, j, i, icor
logical lfail_corn
integer(kind=4) ist(4),jst(4),jdir(4),idir(4)
real(kind=8) zpw(-1:2,-1:2,whw)

ist = [0,n+1,n+1,0  ]
idir= [1, -1, -1, 1]
jst = [0, 0, n+1,n+1]
jdir= [1, 1,  -1, -1]

!check if we have enough data for corner-halo-procedure
lfail_corn = (lcorn(1) .and. (isv> -2 .or. iev<  3 .or. jsv> -2 .or. jev<  3)) .or.& !corner 1
             (lcorn(2) .and. (isv>n-2 .or. iev<n+3 .or. jsv> -2 .or. jev<  3)) .or.& !corner 2
             (lcorn(3) .and. (isv>n-2 .or. iev<n+3 .or. jsv>n-2 .or. jev<n+3)) .or.& !corner 3
             (lcorn(4) .and. (isv> -2 .or. iev<  3 .or. jsv>n-2 .or. jev<n+3))       !corner 4
if(lfail_corn) call parcomm_global%abort("Must be at least 3 halo points for interpolation at corners")

do icor = 1,4
    if(lcorn(icor)) then
        do k=1, klev
            do j=1-nvj,nvj
                do i=1-nvi,nvi
                    zbufc(i,j,k) = pf(ist(icor)+idir(icor)*i,jst(icor)+jdir(icor)*j,k)
                end do
            end do
        end do
        !adjust weights to the case of corner 1
        if(idir(icor) == 1) then
            zpw(-1:2,-1:2,1:whw) = w(-1:2,-1:2,1:whw)
        else
            zpw(-1:2,-1:2,1:whw) = w(2:-1:-1,n+2:n-1:-1,1:whw)
        end if
        call ecs_halo_1corner(pfcsp(1:6,1:klev,icor),zbufc,nvi,nvj,klev,zpw,-1,2,whw,hw)
    end if
end do

end subroutine ecs_halo_corners

subroutine ecs_halo_1corner(pfcsp,fc,nvi,nvj,klev,pw,ws,we,whw,hw)
!output
real(kind=8),   intent(out) :: pfcsp(6,klev)
!input:
real(kind=8),    intent(in) :: fc(1-nvi:nvi,1-nvj:nvj,klev)
integer(kind=4), intent(in) :: nvi, nvj, klev
integer(kind=4), intent(in) :: ws, we, whw, hw
real(kind=8),    intent(in) :: pw(-1:2,ws:we,whw)
!locals:
integer(kind=4) k
real(kind=8) zpw(-1:2)
real(kind=8) zx, zy, zz
real(kind=8) zw, zw2, zw3
real(kind=8) zfy1, zfy2, zff1, zff2, zff11, zff22

!pfcsp = 0._8

zpw= pw(:,1,1)
zw = pw(-1,1,1);   zw2 = zw*zw;   zw3 = zw2*zw

do k=1, klev
    zx = zpw(0)*fc(0,1,k)+zpw(1)*fc(0, 2,k)+zpw(2)*fc(0, 3,k)
    zy = zpw(0)*fc(1,0,k)+zpw(1)*fc(1,-1,k)+zpw(2)*fc(1,-2,k)
    zz = zpw(0)*fc(1,1,k)+zpw(1)*fc(2, 1,k)+zpw(2)*fc(3, 1,k)
    pfcsp(1,k) = (zx+zw*zy+zw2*zz)/(1-zw3)
    zfy1       = (zy+zw*zz+zw2*zx)/(1-zw3)

    zx = zpw(0)*fc(1,0,k)+zpw(1)*fc( 2,0,k)+zpw(2)*fc( 3,0,k)
    zy = zpw(0)*fc(0,1,k)+zpw(1)*fc(-1,1,k)+zpw(2)*fc(-2,1,k)
    zz = zpw(0)*fc(1,1,k)+zpw(1)*fc( 1,2,k)+zpw(2)*fc( 1,3,k)
    pfcsp(2,k) = (zx+zw*zy+zw2*zz)/(1-zw3)
    zfy2       = (zy+zw*zz+zw2*zx)/(1-zw3)

    zff1 = pw(-1,1,2)*fc(2,0,k)+pw(0,1,2)*fc(2,-1,k)+pw(1,1,2)*fc(2,-2,k)+pw(2,1,2)*fc(2,-3,k)
    zff2 = pw(-1,1,2)*fc(0,2,k)+pw(0,1,2)*fc(-1,2,k)+pw(1,1,2)*fc(-2,2,k)+pw(2,1,2)*fc(-3,2,k)
    pfcsp(3,k) = 0.5_8*(pw(-1,0,1)*zff1+pw(0,0,1)*zfy1+pw(1,0,1)*fc(0,1,k)+pw(2,0,1)*fc(0,2,k)+&
                        pw(-1,0,1)*zff2+pw(0,0,1)*zfy2+pw(1,0,1)*fc(1,0,k)+pw(2,0,1)*fc(2,0,k))
    if(hw>1) then
        zff1  = pw(-1,2,1)*fc(1,0,k)+pw(0,2,1)*fc(1,-1,k)+pw(1,2,1)*fc(1,-2,k)+pw(2,2,1)*fc(1,-3,k)
        zff11 = pw(-1,2,2)*fc(2,-1,k)+pw(0,2,2)*fc(2,-2,k)+pw(1,2,2)*fc(2,-3,k)+pw(2,2,2)*fc(2,-4,k)
        zff2  = pw(-1,2,1)*fc(0,1,k)+pw(0,2,1)*fc(-1,1,k)+pw(1,2,1)*fc(-2,1,k)+pw(2,2,1)*fc(-3,1,k)
        zff22 = pw(-1,2,2)*fc(-1,2,k)+pw(0,2,2)*fc(-2,2,k)+pw(1,2,2)*fc(-3,2,k)+pw(2,2,2)*fc(-4,2,k)

        pfcsp(4,k) = pw(-1,0,2)*zff1+pw(0,0,2)*fc(-1,1,k)+pw(1,0,2)*fc(-1,2,k)+pw(2,0,2)*fc(-1,3,k)
        pfcsp(5,k) = pw(-1,0,2)*zff2+pw(0,0,2)*fc(1,-1,k)+pw(1,0,2)*fc(2,-1,k)+pw(2,0,2)*fc(3,-1,k)

        pfcsp(6,k) = 0.5_8*(pw(-1,-1,2)*zff11+pw(0,-1,2)*zff1+pw(1,-1,2)*fc(-1,1,k)+pw(2,-1,2)*fc(-1,2,k)+&
                            pw(-1,-1,2)*zff22+pw(0,-1,2)*zff2+pw(1,-1,2)*fc(1,-1,k)+pw(2,-1,2)*fc(2,-1,k))
    end if
end do

end subroutine ecs_halo_1corner

subroutine ecs_halo_corners_put(pf,isv,iev,jsv,jev,klev,pfcsp,n,lcorn,hw)
!output
real(kind=8), intent(inout) :: pf(isv:iev,jsv:jev,klev) !output field
!input
integer(kind=4), intent(in) :: isv,iev,jsv,jev,klev   !dimensions of pf
real(kind=8),    intent(in) :: pfcsp(6,klev,4)          !previously interpolated values at 'special' corner points
integer(kind=4), intent(in) :: n, hw
logical,         intent(in) :: lcorn(4)                 !if pf corner is cubedsphere corner

if(lcorn(1)) then
    pf( 0, 1,1:klev) = pfcsp(1,1:klev,1)
    pf( 1, 0,1:klev) = pfcsp(2,1:klev,1)
    pf( 0, 0,1:klev) = pfcsp(3,1:klev,1)
    if(hw>1) then
        pf(-1, 0,1:klev) = pfcsp(4,1:klev,1)
        pf( 0,-1,1:klev) = pfcsp(5,1:klev,1)
        pf(-1,-1,1:klev) = pfcsp(6,1:klev,1)
    end if
end if
if(lcorn(2)) then
    pf(n+1, 1,1:klev) = pfcsp(1,1:klev,2)
    pf(n  , 0,1:klev) = pfcsp(2,1:klev,2)
    pf(n+1, 0,1:klev) = pfcsp(3,1:klev,2)
    if(hw>1) then
        pf(n+2, 0,1:klev) = pfcsp(4,1:klev,2)
        pf(n+1,-1,1:klev) = pfcsp(5,1:klev,2)
        pf(n+2,-1,1:klev) = pfcsp(6,1:klev,2)
    end if
end if
if(lcorn(3)) then
    pf(n+1,n  ,1:klev) = pfcsp(1,1:klev,3)
    pf(n  ,n+1,1:klev) = pfcsp(2,1:klev,3)
    pf(n+1,n+1,1:klev) = pfcsp(3,1:klev,3)
    if(hw>1) then
        pf(n+2,n+1,1:klev) = pfcsp(4,1:klev,3)
        pf(n+1,n+2,1:klev) = pfcsp(5,1:klev,3)
        pf(n+2,n+2,1:klev) = pfcsp(6,1:klev,3)
    end if
end if
if(lcorn(4)) then
    pf( 0,n  ,1:klev) = pfcsp(1,1:klev,4)
    pf( 1,n+1,1:klev) = pfcsp(2,1:klev,4)
    pf( 0,n+1,1:klev) = pfcsp(3,1:klev,4)
    if(hw>1) then
        pf(-1,n+1,1:klev) = pfcsp(4,1:klev,4)
        pf( 0,n+2,1:klev) = pfcsp(5,1:klev,4)
        pf(-1,n+2,1:klev) = pfcsp(6,1:klev,4)
    end if
end if
end subroutine ecs_halo_corners_put

subroutine ecs_halo_edges_x(pf,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                            wx,wsx,wex,whw,indx,                         & !weights and indices  arrays and their bounds
                            lhalo,n,hw)                                    !global operating parameters

!in-output
real(kind=8), intent(inout) :: pf(isv:iev,jsv:jev,klev)
!input
integer(kind=4), intent(in) :: is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj
!interpolation weights
integer(kind=4), intent(in) :: wsx,wex,whw
real(kind=8),    intent(in) :: wx(-1:2,wsx:wex,whw) !interpolation weights
integer(kind=4), intent(in) :: indx(wsx:wex,whw)
!operting parameters
logical,         intent(in) :: lhalo(4)
integer(kind=4), intent(in) :: n, hw
!local
integer(kind=4) k, j, i
integer(kind=4) ish, ieh, jsh, jeh
logical lfail_hw
logical lfail_halo_long
real(kind=8) zbufx(isv:iev, hw, klev)

!Internal checks to avoid (at least some part of) out-of-array errors
!check if we have enough data for edge halo-procedure
lfail_hw = nvj< hw
if(lfail_hw) call parcomm_global%abort("cubed sphere halo_width > grid_function_t halo width, can't continue")
!check if we have enough data along edges to perform interpolations
ish = max(1,is-hw); ieh = min(n,ie+hw)
lfail_halo_long = .false.
lfail_halo_long = (minval(indx(ish,1:hw))-1<isv .or. maxval(indx(ieh,1:hw))+2>iev)
if(lfail_halo_long) call parcomm_global%abort("not enough points along cub.sph edge to perform halo-interpolations (nvi=5 needed)")

!calculate values at corner special points
!store them in sepparate arrays to not to spoil original f%p

!Halo-procedures along edges
if(lhalo(1)) then
    zbufx(isv:iev,1:hw,1:klev) = pf(isv:iev, 0:1-hw:-1,1:klev)
    call ecs_ext_halo_1e(zbufx,isv,iev,is,ie,hw,klev,wx,wsx,wex,whw,indx,n)
    pf(isv:iev,1-hw:0,1:klev) = zbufx(isv:iev,hw:1:-1,1:klev)
end if
if(lhalo(2)) then
    zbufx(isv:iev,1:hw,1:klev) = pf(isv:iev, n+1:n+hw, 1:klev)
    call ecs_ext_halo_1e(zbufx,isv,iev,is,ie,hw,klev,wx,wsx,wex,whw,indx,n)
    pf(isv:iev, n+1:n+hw, 1:klev) = zbufx(isv:iev,1:hw,1:klev)
end if

end subroutine ecs_halo_edges_x

subroutine ecs_halo_edges_y(pf,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                            wy,wsy,wey,whw,indy,                         & !weights and indices  arrays and their bounds
                            lhalo,n,hw)                                    !global operating parameters

!in-output
real(kind=8), intent(inout) :: pf(isv:iev,jsv:jev,klev)
!input
integer(kind=4), intent(in) :: is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj
!interpolation weights
integer(kind=4), intent(in) :: wsy,wey,whw
real(kind=8),    intent(in) :: wy(-1:2,wsy:wey,whw) !interpolation weights
integer(kind=4), intent(in) :: indy(wsy:wey,whw)
!operting parameters
logical,         intent(in) :: lhalo(4)
integer(kind=4), intent(in) :: n, hw
!local
integer(kind=4) k, j, i
integer(kind=4) jsh, jeh
logical lfail_hw
logical lfail_halo_long
real(kind=8) zbufy(jsv:jev, hw, klev)

!Internal checks to avoid (at least some part of) out-of-array errors
!check if we have enough data for edge halo-procedure
lfail_hw = nvj< hw
if(lfail_hw) call parcomm_global%abort("cubed sphere halo_width > grid_function_t halo width, can't continue")
!check if we have enough data along edges to perform interpolations
jsh = max(1,js-hw); jeh = min(n,je+hw)
lfail_halo_long = (minval(indy(jsh,1:hw))-1<jsv .or. maxval(indy(jeh,1:hw))+2>jev)
if(lfail_halo_long) call parcomm_global%abort("not enough points along cub.sph edge to perform halo-interpolations (nvi=5 needed)")

!calculate values at corner special points
!store them in sepparate arrays to not to spoil original f%p

!Halo-procedures along edges
if(lhalo(3)) then
    do k=1,klev
        do j=jsv, jev
            zbufy(j,1:hw,k) = pf(0:1-hw:-1,j,k)
        end do
    end do
    call ecs_ext_halo_1e(zbufy,jsv,jev,js,je,hw,klev,wy,wsy,wey,whw,indy,n)
    do k=1,klev
        do j=1,hw
            pf(1-j,jsv:jev,k) = zbufy(jsv:jev,j,k)
        end do
    end do
end if
if(lhalo(4)) then
    do k=1, klev
        do j=jsv,jev
            zbufy(j,1:hw,k) = pf(n+1:n+hw,j,k)
        end do
    end do
    call ecs_ext_halo_1e(zbufy,jsv,jev,js,je,hw,klev,wy,wsy,wey,whw,indy,n)
    do k=1,klev
        do j=1,hw
            pf(n+j,jsv:jev,k) = zbufy(jsv:jev,j,k)
        end do
    end do
end if

end subroutine ecs_halo_edges_y

subroutine ecs_ext_halo_1e(zf,i1v,i2v,i1,i2,hw,klev,w,ws,we,whw,ind,n)
!input
integer(kind=4), intent(in) :: i1v,i2v,i1,i2,hw,klev
!in-output:
real(kind=8), intent(inout) :: zf(i1v:i2v,hw,klev)!input: source face values, output: interpolated target face values
!input
integer(kind=4), intent(in) :: ws,we,whw
real(kind=8),    intent(in) :: w(-1:2,ws:we,whw)
integer(kind=4), intent(in) :: ind(ws:we,whw)
integer(kind=4), intent(in) :: n
!locals
real(kind=8) zh(i1-hw:i2+hw)!buffer for interpolated values
integer i, j, k, ii
integer ihs(hw), ihe(hw)

!exclude "special" points at first row, where we initially have not enough surrounding values to use cubic-interpolation
ihs(2:) = 1; ihs(1) = 2
ihe(2:) = n; ihe(1) = n-1

do k=1, klev
    do j=1, hw
        do i = max(i1-hw,ihs(j)),min(i2+hw,ihe(j))
            ii = ind(i,j)
            zh(i) = sum(w(:,i,j)*zf(ii-1:ii+2,j,k))
        end do
        zf(max(i1-hw,ihs(j)):min(i2+hw,ihe(j)),j,k) = zh(max(i1-hw,ihs(j)):min(i2+hw,ihe(j)))
    end do
end do

end subroutine ecs_ext_halo_1e

end module ecs_halo_mod

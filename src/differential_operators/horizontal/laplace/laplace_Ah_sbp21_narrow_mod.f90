module laplace_Ah_sbp21_narrow_mod

use abstract_laplace_mod,   only : laplace_operator_t
use grid_field_mod,         only : grid_field_t
use domain_mod,             only : domain_t
use sbp_operator_mod,       only : sbp_operator_t
use exchange_abstract_mod,  only : exchange_t
use halo_mod,               only : halo_t

implicit none

type, extends(laplace_operator_t) :: laplace_Ah_sbp21_narrow_t
    type(grid_field_t)             :: d2f_x, d2f_y, d1f_x, d1f_y, d2f_xy, d2f_yx
    type(sbp_operator_t)           :: sbp_d1_op
    class(halo_t), allocatable     :: edge_sync
    class(exchange_t), allocatable :: exchange
contains
    procedure :: calc_laplace
end type laplace_Ah_sbp21_narrow_t

contains

subroutine calc_laplace(this, f1, f, domain)

    use vec_math_mod,   only : divide_by_J_self

    class(laplace_Ah_sbp21_narrow_t), intent(inout) :: this
    type(domain_t),                   intent(in)    :: domain
    type(grid_field_t),               intent(inout) :: f
    !output:
    type(grid_field_t),               intent(inout) :: f1

    integer(kind=4) :: t

    call this%exchange%do(f, domain%parcomm)
    do t = domain%partition%ts, domain%partition%te

        call calc_tile_d1x_d1y(this%d1f_x%tile(t), this%d1f_y%tile(t), f%tile(t), this%sbp_d1_op, domain%mesh_xy%tile(t))

    end do

    call this%exchange%do(this%d1f_x, domain%parcomm)
    call this%exchange%do(this%d1f_y, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te

        call calc_laplace_tile(f1%tile(t), this%d2f_x%tile(t), this%d2f_y%tile(t), &
                                      this%d2f_xy%tile(t), this%d2f_yx%tile(t), this%d1f_x%tile(t), this%d1f_y%tile(t), &
                                f%tile(t), &
                                this%sbp_d1_op, domain%mesh_xy%tile(t), domain%mesh_xy%scale)

    end do

    call divide_by_J_self(f1, domain%mesh_xy)

    call this%edge_sync%get_halo_scalar(f1,domain,0)

end subroutine calc_laplace

subroutine calc_tile_d1x_d1y(d1f_x, d1f_y, f, sbp_d1_op, mesh_xy)

    use grid_field_mod, only : grid_field_t, tile_field_t
    use mesh_mod,       only : tile_mesh_t
    use tile_mod,       only : tile_t

    !result
    type(tile_field_t), intent(inout) :: d1f_x, d1f_y
    !input
    type(tile_field_t),   intent(in) :: f
    type(sbp_operator_t), intent(in) :: sbp_d1_op
    type(tile_mesh_t),    intent(in) :: mesh_xy

    integer(kind=4) :: i, j, k

    type(tile_t) :: work

    work = tile_t(is = mesh_xy%is, ie = mesh_xy%ie, &
                  js = mesh_xy%js, je = mesh_xy%je, &
                  ks = mesh_xy%ks, ke = mesh_xy%ke)

    call sbp_d1_op%apply(d1f_x, work, mesh_xy%nx, "x", f)
    call sbp_d1_op%apply(d1f_y, work, mesh_xy%ny, "y", f)

    do k = mesh_xy%ks, mesh_xy%ke
        do j = mesh_xy%js, mesh_xy%je
            do i = mesh_xy%is, mesh_xy%ie
                d1f_x%p(i,j,k) = d1f_x%p(i,j,k)*mesh_xy%J(i,j,k)*mesh_xy%Qi(2,i,j,k)
                d1f_y%p(i,j,k) = d1f_y%p(i,j,k)*mesh_xy%J(i,j,k)*mesh_xy%Qi(2,i,j,k)
            end do
        end do
    end do

end subroutine calc_tile_d1x_d1y
subroutine calc_laplace_tile(fout,d2f_x, d2f_y, d2f_xy, d2f_yx, d1f_x, d1f_y, f, sbp_d1_op, mesh_xy, scale)

    use grid_field_mod, only : grid_field_t, tile_field_t
    use mesh_mod,       only : tile_mesh_t
    use tile_mod,       only : tile_t

    !result
    type(tile_field_t), intent(inout) :: fout, d2f_x, d2f_y, d2f_xy, d2f_yx, d1f_x, d1f_y
    !input
    type(tile_field_t),   intent(in) :: f
    type(sbp_operator_t), intent(in) :: sbp_d1_op
    type(tile_mesh_t),    intent(in) :: mesh_xy
    real(kind=8),         intent(in) :: scale

    integer(kind=4) :: i, j, k
    real(kind=8)    :: penalty, bm, bc, bp, cm, cc, cp, cmm, cpp

    type(tile_t) :: work

    work = tile_t(is = mesh_xy%is, ie = mesh_xy%ie, &
                  js = mesh_xy%js, je = mesh_xy%je, &
                  ks = mesh_xy%ks, ke = mesh_xy%ke)

    do k = mesh_xy%ks, mesh_xy%ke

      if (mesh_xy%is==1) then
          do j = mesh_xy%js, mesh_xy%je

              bc  = mesh_xy%Qi(1,1,j,k)*mesh_xy%J(1,j,k)
              bp  = mesh_xy%Qi(1,2,j,k)*mesh_xy%J(2,j,k)

              cc  =  -(bc+bp)
              cp  = + (bc+bp)

              d2f_x%p(1,j,k) = (cc*f%p(1,j,k)+cp*f%p(2,j,k))
          end do
      end if
      do j = mesh_xy%js, mesh_xy%je
          do i = max(mesh_xy%is,2), min(mesh_xy%ie,mesh_xy%nx-1)

              bm = mesh_xy%Qi(1,i-1,j,k)*mesh_xy%J(i-1,j,k)
              bc = mesh_xy%Qi(1,i  ,j,k)*mesh_xy%J(i  ,j,k)
              bp = mesh_xy%Qi(1,i+1,j,k)*mesh_xy%J(i+1,j,k)

              cm =  (bm+bc)/2
              cc = -(bm+2*bc+bp)/2
              cp =  (bp+bc)/2

              d2f_x%p(i,j,k) = (cp*f%p(i+1,j,k)+cc*f%p(i,j,k)+cm*f%p(i-1,j,k))
          end do
      end do
      if (mesh_xy%ie==mesh_xy%nx) then
          do j = mesh_xy%js, mesh_xy%je

              bc = mesh_xy%Qi(1,mesh_xy%nx  ,j,k)*mesh_xy%J(mesh_xy%nx  ,j,k)
              bm = mesh_xy%Qi(1,mesh_xy%nx-1,j,k)*mesh_xy%J(mesh_xy%nx-1,j,k)

              cc  =  -(bc+bm)
              cm  =  +(bc+bm)

              d2f_x%p(mesh_xy%nx,j,k) = cm*f%p(mesh_xy%nx-1,j,k)+cc*f%p(mesh_xy%nx,j,k)
          end do
      end if

      if (mesh_xy%js==1) then
          do i = mesh_xy%is, mesh_xy%ie

              bc  = mesh_xy%Qi(3,i,1,k)*mesh_xy%J(i,1,k)
              bp  = mesh_xy%Qi(3,i,2,k)*mesh_xy%J(i,2,k)

              cc  = -(bc+bp)
              cp  = +(bc+bp)


              d2f_y%p(i,1,k) = cc*f%p(i,1,k)+cp*f%p(i,2,k)
          end do
      end if
      do j = max(mesh_xy%js,2), min(mesh_xy%je,mesh_xy%ny-1)
          do i = mesh_xy%is, mesh_xy%ie

              bm = mesh_xy%Qi(3,i,j-1,k)*mesh_xy%J(i,j-1,k)
              bc = mesh_xy%Qi(3,i,j  ,k)*mesh_xy%J(i,j  ,k)
              bp = mesh_xy%Qi(3,i,j+1,k)*mesh_xy%J(i,j+1,k)

              cm =  (bm+bc)/2
              cc = -(bm+2*bc+bp)/2
              cp =  (bp+bc)/2

              d2f_y%p(i,j,k) = (cp*f%p(i,j+1,k)+cc*f%p(i,j,k)+cm*f%p(i,j-1,k))
          end do
      end do
      if (mesh_xy%je==mesh_xy%ny) then
          do i = mesh_xy%is, mesh_xy%ie

              bc = mesh_xy%Qi(3,i,mesh_xy%ny  ,k)*mesh_xy%J(i,mesh_xy%ny  ,k)
              bm = mesh_xy%Qi(3,i,mesh_xy%ny-1,k)*mesh_xy%J(i,mesh_xy%ny-1,k)

              cc  = -(bc+bm)
              cm  = +(bc+bm)

              d2f_y%p(i,mesh_xy%ny, k) = cm*f%p(i,mesh_xy%ny-1,k)+cc*f%p(i,mesh_xy%ny,k)

          end do
      end if
    end do

    do k = mesh_xy%ks, mesh_xy%ke

      if (mesh_xy%is==1) then
          do j = mesh_xy%js, mesh_xy%je
              penalty = 2.0_8*(-2.0_8*d1f_y%p(1,j,k))
              d2f_x%p(1,j,k) = d2f_x%p(1,j,k)-0.5_8*penalty
          end do
      end if
      if (mesh_xy%ie==mesh_xy%nx) then
          do j = mesh_xy%js, mesh_xy%je
              penalty = 2.0_8*(2*d1f_y%p(mesh_xy%nx,j,k))
              d2f_x%p(mesh_xy%nx,j,k) = d2f_x%p(mesh_xy%nx,j,k) -0.5_8*penalty
          end do
      end if

      if (mesh_xy%js==1) then
          do i = mesh_xy%is, mesh_xy%ie
              penalty = 2.0_8*(-2*d1f_x%p(i,1,k))
              d2f_x%p(i,1,k) = d2f_x%p(i,1,k) -0.5_8*penalty
          end do
      end if
      if (mesh_xy%je==mesh_xy%ny) then
          do i = mesh_xy%is, mesh_xy%ie
              penalty = 2.0_8*(2*d1f_x%p(i,mesh_xy%ny,k))
              d2f_x%p(i,mesh_xy%ny,k) = d2f_x%p(i,mesh_xy%ny,k) -0.5_8*penalty
          end do
      end if

    end do
    call sbp_d1_op%apply(d2f_yx, work, mesh_xy%ny, "y", d1f_x)
    call sbp_d1_op%apply(d2f_xy, work, mesh_xy%nx, "x", d1f_y)

    call fout%assign(1.0_8, d2f_x, &
                     1.0_8, d2f_y, mesh_xy)

    call fout%update(1.0_8, d2f_xy, &
                     1.0_8, d2f_yx, mesh_xy)


     do k = mesh_xy%ks, mesh_xy%ke
         do j = mesh_xy%js, mesh_xy%je
             do i = mesh_xy%is, mesh_xy%ie
                 fout%p(i,j,k) = fout%p(i,j,k)/(mesh_xy%hx*scale)**2
             end do
         end do
     end do

end subroutine calc_laplace_tile
! subroutine calc_laplace_tile(fout,d2f_x, d2f_y, d2f_xy, d2f_yx, d1f_x, d1f_y, f, sbp_d1_op, mesh_xy, scale)
!
!     use grid_field_mod, only : grid_field_t, tile_field_t
!     use mesh_mod,       only : tile_mesh_t
!     use tile_mod,       only : tile_t
!
!     !result
!     type(tile_field_t), intent(inout) :: fout, d2f_x, d2f_y, d2f_xy, d2f_yx, d1f_x, d1f_y
!     !input
!     type(tile_field_t),   intent(in) :: f
!     type(sbp_operator_t), intent(in) :: sbp_d1_op
!     type(tile_mesh_t),    intent(in) :: mesh_xy
!     real(kind=8),         intent(in) :: scale
!
!     integer(kind=4) :: i, j, k
!     real(kind=8)    :: penalty, bm, bc, bp, cm, cc, cp, cmm, cpp
!
!     type(tile_t) :: work
!
!     work = tile_t(is = mesh_xy%is, ie = mesh_xy%ie, &
!                   js = mesh_xy%js, je = mesh_xy%je, &
!                   ks = mesh_xy%ks, ke = mesh_xy%ke)
!
!     do k = mesh_xy%ks, mesh_xy%ke
!
!       if (mesh_xy%is==1) then
!           do j = mesh_xy%js, mesh_xy%je
!
!               bc  = mesh_xy%Qi(1,1,j,k)*mesh_xy%J(1,j,k)
!               bp  = mesh_xy%Qi(1,2,j,k)*mesh_xy%J(2,j,k)
!
!               cc  =  -(bc+bp)
!               cp  = + (bc+bp)
!               cpp =    0
!
!               d2f_x%p(1,j,k) = (cc*f%p(1,j,k)+cp*f%p(2,j,k)+cpp*f%p(3,j,k))
!           end do
!       end if
!       do j = mesh_xy%js, mesh_xy%je
!           do i = max(mesh_xy%is,2), min(mesh_xy%ie,mesh_xy%nx-1)
!
!               bm = mesh_xy%Qi(1,i-1,j,k)*mesh_xy%J(i-1,j,k)
!               bc = mesh_xy%Qi(1,i  ,j,k)*mesh_xy%J(i  ,j,k)
!               bp = mesh_xy%Qi(1,i+1,j,k)*mesh_xy%J(i+1,j,k)
!
!               cm =  (bm+bc)/2
!               cc = -(bm+2*bc+bp)/2
!               cp =  (bp+bc)/2
!
!               d2f_x%p(i,j,k) = (cp*f%p(i+1,j,k)+cc*f%p(i,j,k)+cm*f%p(i-1,j,k))
!           end do
!       end do
!       if (mesh_xy%ie==mesh_xy%nx) then
!           do j = mesh_xy%js, mesh_xy%je
!
!               bc = mesh_xy%Qi(1,mesh_xy%nx  ,j,k)*mesh_xy%J(mesh_xy%nx  ,j,k)
!               bm = mesh_xy%Qi(1,mesh_xy%nx-1,j,k)*mesh_xy%J(mesh_xy%nx-1,j,k)
!
!               cc  =  -(bc+bm)
!               cm  =  +(bc+bm)
!               cmm =    0
!
!               d2f_x%p(mesh_xy%nx,j,k) = (cmm*f%p(mesh_xy%nx-2,j,k)+cm*f%p(mesh_xy%nx-1,j,k)+cc*f%p(mesh_xy%nx,j,k))
!           end do
!       end if
!
!       if (mesh_xy%js==1) then
!           do i = mesh_xy%is, mesh_xy%ie
!
!               bc  = mesh_xy%Qi(3,i,1,k)*mesh_xy%J(i,1,k)
!               bp  = mesh_xy%Qi(3,i,2,k)*mesh_xy%J(i,2,k)
!
!               cc  = -(bc+bp)
!               cp  = +(bc+bp)
!               cpp =    0
!
!               d2f_y%p(i,1,k) = (cc*f%p(i,1,k)+cp*f%p(i,2,k)+cpp*f%p(i,3,k))
!           end do
!       end if
!       do j = max(mesh_xy%js,2), min(mesh_xy%je,mesh_xy%ny-1)
!           do i = mesh_xy%is, mesh_xy%ie
!
!               bm = mesh_xy%Qi(3,i,j-1,k)*mesh_xy%J(i,j-1,k)
!               bc = mesh_xy%Qi(3,i,j  ,k)*mesh_xy%J(i,j  ,k)
!               bp = mesh_xy%Qi(3,i,j+1,k)*mesh_xy%J(i,j+1,k)
!
!               cm =  (bm+bc)/2
!               cc = -(bm+2*bc+bp)/2
!               cp =  (bp+bc)/2
!
!               d2f_y%p(i,j,k) = (cp*f%p(i,j+1,k)+cc*f%p(i,j,k)+cm*f%p(i,j-1,k))
!           end do
!       end do
!       if (mesh_xy%je==mesh_xy%ny) then
!           do i = mesh_xy%is, mesh_xy%ie
!
!               bc = mesh_xy%Qi(3,i,mesh_xy%ny  ,k)*mesh_xy%J(i,mesh_xy%ny  ,k)
!               bm = mesh_xy%Qi(3,i,mesh_xy%ny-1,k)*mesh_xy%J(i,mesh_xy%ny-1,k)
!
!               cc  = -(bc+bm)
!               cm  = +(bc+bm)
!               cmm =    0
!
!               d2f_y%p(i,mesh_xy%ny, k) = (cmm*f%p(i,mesh_xy%ny-2,k)+cm*f%p(i,mesh_xy%ny-1,k)+cc*f%p(i,mesh_xy%ny,k))
!
!           end do
!       end if
!     end do
!
!     do k = mesh_xy%ks, mesh_xy%ke
!
!       if (mesh_xy%is==1) then
!           do j = mesh_xy%js, mesh_xy%je
!
!               bc  = mesh_xy%Qi(1,1,j,k)*mesh_xy%J(1,j,k)
!
!               penalty = (( 3*f%p(0,j,k) -4*f%p(-1,j,k) + f%p(-2,j,k)) &
!                         +( 3*f%p(1,j,k) -4*f%p (2,j,k) + f%p( 3,j,k)))*bc*0 &
!                         + 2.0_8*(-2.0_8*d1f_y%p(1,j,k))
!
!                d2f_x%p(1,j,k) = d2f_x%p(1,j,k)-0.5_8*penalty
!           end do
!       end if
!       if (mesh_xy%ie==mesh_xy%nx) then
!           do j = mesh_xy%js, mesh_xy%je
!
!               bc = mesh_xy%Qi(1,mesh_xy%nx,j,k)*mesh_xy%J(mesh_xy%nx,j,k)
!
!               penalty = ( ( 3*f%p(mesh_xy%nx  ,j,k)-4*f%p(mesh_xy%nx-1,j,k) + f%p(mesh_xy%nx-2,j,k)) &
!                          +( 3*f%p(mesh_xy%nx+1,j,k)-4*f%p(mesh_xy%nx+2,j,k) + f%p(mesh_xy%nx+3,j,k)))*bc*0 &
!                          + 2.0_8*(2*d1f_y%p(mesh_xy%nx,j,k))
!               d2f_x%p(mesh_xy%nx,j,k) = d2f_x%p(mesh_xy%nx,j,k) -0.5_8*penalty
!           end do
!       end if
!
!       if (mesh_xy%js==1) then
!           do i = mesh_xy%is, mesh_xy%ie
!
!               bc  = mesh_xy%Qi(3,i,1,k)*mesh_xy%J(i,1,k)
!               penalty = ( ( 3*f%p(i,0,k) -4*f%p(i,-1,k) + f%p(i,-2,k)) &
!                          +( 3*f%p(i,1,k) -4*f%p(i, 2,k) + f%p(i, 3,k)))*bc*0 &
!                          + 2.0_8*(-2*d1f_x%p(i,1,k))
!
!               d2f_x%p(i,1,k) = d2f_x%p(i,1,k) -0.5_8*penalty
!           end do
!       end if
!       if (mesh_xy%je==mesh_xy%ny) then
!           do i = mesh_xy%is, mesh_xy%ie
!
!               bc = mesh_xy%Qi(3,i,mesh_xy%ny,k)*mesh_xy%J(i,mesh_xy%ny,k)
!               penalty = (( 3*f%p(i,mesh_xy%ny  ,k) -4*f%p(i,mesh_xy%ny-1,k) + f%p(i,mesh_xy%ny-2,k)) &
!                         +( 3*f%p(i,mesh_xy%ny+1,k) -4*f%p(i,mesh_xy%ny+2,k) + f%p(i,mesh_xy%ny+3,k)))*bc*0 &
!                         + 2.0_8*(2*d1f_x%p(i,mesh_xy%ny,k))
!
!               d2f_x%p(i,mesh_xy%ny,k) = d2f_x%p(i,mesh_xy%ny,k) -0.5_8*penalty
!
!           end do
!       end if
!
!     end do
!     call sbp_d1_op%apply(d2f_yx, work, mesh_xy%ny, "y", d1f_x)
!     call sbp_d1_op%apply(d2f_xy, work, mesh_xy%nx, "x", d1f_y)
!
!     call fout%assign(1.0_8, d2f_x, &
!                      1.0_8, d2f_y, mesh_xy)
!
!     call fout%update(1.0_8, d2f_xy, &
!                      1.0_8, d2f_yx, mesh_xy)
!
!
!      do k = mesh_xy%ks, mesh_xy%ke
!          do j = mesh_xy%js, mesh_xy%je
!              do i = mesh_xy%is, mesh_xy%ie
!                  fout%p(i,j,k) = fout%p(i,j,k)/(mesh_xy%hx*scale)**2
!              end do
!          end do
!      end do
!
! end subroutine calc_laplace_tile

end module laplace_Ah_sbp21_narrow_mod

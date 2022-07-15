module laplace_Ah_sbp42_narrow_mod

use abstract_laplace_mod,   only : laplace_operator_t
use grid_field_mod,         only : grid_field_t
use domain_mod,             only : domain_t
use sbp_operator_mod,       only : sbp_operator_t
use exchange_abstract_mod,  only : exchange_t
use halo_mod,               only : halo_t

implicit none

type, extends(laplace_operator_t) :: laplace_Ah_sbp42_narrow_t
    type(grid_field_t)             :: d2f_x, d2f_y, d1f_x, d1f_y, d2f_xy, d2f_yx
    type(sbp_operator_t)           :: sbp_d1_op
    class(halo_t), allocatable     :: edge_sync
    class(exchange_t), allocatable :: exchange
contains
    procedure :: calc_laplace
end type laplace_Ah_sbp42_narrow_t

contains

subroutine calc_laplace(this, f1, f, domain)

    use vec_math_mod,   only : divide_by_J_self

    class(laplace_Ah_sbp42_narrow_t), intent(inout) :: this
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
    real(kind=8)    :: penalty, bmm, bm, bc, bp, bpp, cm, cc, cp, cmm, cpp
    real(kind=8) :: b1, b2, b3, b4, b5, b6, b7, b8
    real(kind=8) :: m11, m12, m13, m14, m15, m16,          &
                    m21, m22, m23, m24, m25, m26,          &
                    m31, m32, m33, m34, m35, m36,          &
                    m41, m42, m43, m44, m45, m46,          &
                    m51, m52, m53, m54, m55, m56, m57,     &
                    m61, m62, m63, m64, m65, m66, m67, m68
    real(kind=8) :: h11, h22, h33, h44, h55, h66
    type(tile_t) :: work

    work = tile_t(is = mesh_xy%is, ie = mesh_xy%ie, &
                  js = mesh_xy%js, je = mesh_xy%je, &
                  ks = mesh_xy%ks, ke = mesh_xy%ke)

    do k = mesh_xy%ks, mesh_xy%ke
      if (mesh_xy%is==1) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,1,j,k)*mesh_xy%J(1,j,k)
              b2  = mesh_xy%Qi(1,2,j,k)*mesh_xy%J(2,j,k)
              b3  = mesh_xy%Qi(1,3,j,k)*mesh_xy%J(3,j,k)
              b4  = mesh_xy%Qi(1,4,j,k)*mesh_xy%J(4,j,k)

              h11 = 17.0_8/48

              m11 = 12.0_8/17.0_8*b1 + 59.0_8/192.0_8*b2 + 27010400129.0_8/345067064608.0_8*b3 + 69462376031.0_8/2070402387648.0_8*b4
              m12 = -59.0_8/68.0_8*b1 -  6025413881.0_8/21126554976.0_8*b3 - 537416663.0_8/7042184992.0_8*b4
              m13 = 2.0_8/17.0_8*b1 - 59.0_8/192.0_8*b2 + 213318005.0_8/16049630912.0_8*b4 + 2083938599.0_8/8024815456.0_8*b3
              m14 = 3.0_8/68.0_8*b1 - 1244724001.0_8/21126554976.0_8*b3 + 752806667.0_8/21126554976.0_8*b4
              m15 = 49579087.0_8/10149031312.0_8*b3 - 49579087.0_8/10149031312.0_8*b4
              m16 = -1.0_8/784.0_8*b4 + 1.0_8/784.0_8*b3

              d2f_x%p(1,j,k) = -1.0_8/h11*(m11*f%p(1,j,k)+m12*f%p(2,j,k)+m13*f%p(3,j,k) + &
                                           m14*f%p(4,j,k)+m15*f%p(5,j,k)+m16*f%p(6,j,k) )
          end do
      end if
      if (mesh_xy%is<=2 .and. mesh_xy%ie>=2) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,1,j,k)*mesh_xy%J(1,j,k)
              b2  = mesh_xy%Qi(1,2,j,k)*mesh_xy%J(2,j,k)
              b3  = mesh_xy%Qi(1,3,j,k)*mesh_xy%J(3,j,k)
              b4  = mesh_xy%Qi(1,4,j,k)*mesh_xy%J(4,j,k)

              h22 = 59.0_8/48

              m21 = -59.0_8/68.0_8*b1 -  6025413881.0_8/21126554976.0_8*b3 - 537416663.0_8/7042184992.0_8*b4
              m22 = 3481.0_8/3264.0_8*b1 + 9258282831623875.0_8/7669235228057664.0_8*b3 + 236024329996203.0_8/1278205871342944.0_8*b4
              m23 = - 59.0_8/408.0_8*b1 - 29294615794607.0_8/29725717938208.0_8*b3 - 2944673881023.0_8/29725717938208.0_8*b4
              m24 = - 59.0_8/1088.0_8*b1 + 260297319232891.0_8/2556411742685888.0_8*b3 - 60834186813841.0_8/1278205871342944.0_8*b4
              m25 = - 1328188692663.0_8/37594290333616.0_8*b3 + 1328188692663.0_8/37594290333616.0_8*b4
              m26 = - 8673.0_8/2904112.0_8*b3 + 8673.0_8/2904112.0_8*b4

              d2f_x%p(2,j,k) = -1.0_8/h22*(m21*f%p(1,j,k)+m22*f%p(2,j,k)+m23*f%p(3,j,k) + &
                                          m24*f%p(4,j,k)+m25*f%p(5,j,k)+m26*f%p(6,j,k) )
          end do
      end if
      if (mesh_xy%is<=3 .and. mesh_xy%ie>=3) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,1,j,k)*mesh_xy%J(1,j,k)
              b2  = mesh_xy%Qi(1,2,j,k)*mesh_xy%J(2,j,k)
              b3  = mesh_xy%Qi(1,3,j,k)*mesh_xy%J(3,j,k)
              b4  = mesh_xy%Qi(1,4,j,k)*mesh_xy%J(4,j,k)
              b5  = mesh_xy%Qi(1,5,j,k)*mesh_xy%J(5,j,k)

              h33 = 43.0_8/48

              m31 = 2.0_8/17.0_8*b1 - 59.0_8/192.0_8*b2 + 213318005.0_8/16049630912.0_8*b4 + 2083938599.0_8/8024815456.0_8*b3
              m32 = - 59.0_8/408.0_8*b1 - 29294615794607.0_8/29725717938208.0_8*b3 - 2944673881023.0_8/29725717938208.0_8*b4
              m33 = 1.0_8/51.0_8*b1 + 59.0_8/192.0_8*b2 + 13777050223300597.0_8/26218083221499456.0_8*b4 + 564461.0_8/13384296.0_8*b5 + 378288882302546512209.0_8/270764341349677687456.0_8*b3
              m34 = 1.0_8/136.0_8*b1 - 125059.0_8/743572.0_8*b5 - 4836340090442187227.0_8/5525802884687299744.0_8*b3 - 17220493277981.0_8/89177153814624.0_8*b4
              m35 = - 10532412077335.0_8/42840005263888.0_8*b4 + 1613976761032884305.0_8/7963657098519931984.0_8*b3 + 564461.0_8/4461432.0_8*b5
              m36 = - 960119.0_8/1280713392.0_8*b4 - 3391.0_8/6692148.0_8*b5 + 33235054191.0_8/26452850508784.0_8*b3

              d2f_x%p(3,j,k) = -1.0_8/h33*(m31*f%p(1,j,k)+m32*f%p(2,j,k)+m33*f%p(3,j,k) + &
                                           m34*f%p(4,j,k)+m35*f%p(5,j,k)+m36*f%p(6,j,k) )
          end do
      end if
      if (mesh_xy%is<=4 .and. mesh_xy%ie>=4) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,1,j,k)*mesh_xy%J(1,j,k)
              b2  = mesh_xy%Qi(1,2,j,k)*mesh_xy%J(2,j,k)
              b3  = mesh_xy%Qi(1,3,j,k)*mesh_xy%J(3,j,k)
              b4  = mesh_xy%Qi(1,4,j,k)*mesh_xy%J(4,j,k)
              b5  = mesh_xy%Qi(1,5,j,k)*mesh_xy%J(5,j,k)
              b6  = mesh_xy%Qi(1,6,j,k)*mesh_xy%J(6,j,k)

              h44 = 49.0_8/48

              m41 = 3.0_8/68.0_8*b1 - 1244724001.0_8/21126554976.0_8*b3 + 752806667.0_8/21126554976.0_8*b4
              m42 = - 59.0_8/1088.0_8*b1 + 260297319232891.0_8/2556411742685888.0_8*b3 - 60834186813841.0_8/1278205871342944.0_8*b4
              m43 = 1.0_8/136.0_8*b1 - 125059.0_8/743572.0_8*b5 - 4836340090442187227.0_8/5525802884687299744.0_8*b3 - 17220493277981.0_8/89177153814624.0_8*b4
              m44 = 3.0_8/1088.0_8*b1 + 507284006600757858213.0_8/475219048083107777984.0_8*b3 + 1869103.0_8/2230716.0_8*b5 + 1.0_8/24.0_8*b6 + 1950062198436997.0_8/3834617614028832.0_8*b4
              m45 = - 4959271814984644613.0_8/20965546238960637264.0_8*b3 - 1.0_8/6.0_8*b6 - 15998714909649.0_8/37594290333616.0_8*b4 - 375177.0_8/743572.0_8*b5
              m46 = - 368395.0_8/2230716.0_8*b5 + 752806667.0_8/539854092016.0_8*b3 + 1063649.0_8/8712336.0_8*b4 + 1.0_8/8.0_8*b6

              d2f_x%p(4,j,k) = -1.0_8/h44*(m41*f%p(1,j,k)+m42*f%p(2,j,k)+m43*f%p(3,j,k) + &
                                          m44*f%p(4,j,k)+m45*f%p(5,j,k)+m46*f%p(6,j,k) )
          end do
      end if
      if (mesh_xy%is<=5 .and. mesh_xy%ie>=5) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,1,j,k)*mesh_xy%J(1,j,k)
              b2  = mesh_xy%Qi(1,2,j,k)*mesh_xy%J(2,j,k)
              b3  = mesh_xy%Qi(1,3,j,k)*mesh_xy%J(3,j,k)
              b4  = mesh_xy%Qi(1,4,j,k)*mesh_xy%J(4,j,k)
              b5  = mesh_xy%Qi(1,5,j,k)*mesh_xy%J(5,j,k)
              b6  = mesh_xy%Qi(1,6,j,k)*mesh_xy%J(6,j,k)
              b7  = mesh_xy%Qi(1,7,j,k)*mesh_xy%J(7,j,k)

              h55 = 1

              m51 = 49579087.0_8/10149031312.0_8*b3 - 49579087.0_8/10149031312.0_8*b4
              m52 = - 1328188692663.0_8/37594290333616.0_8*b3 + 1328188692663.0_8/37594290333616.0_8*b4
              m53 = - 10532412077335.0_8/42840005263888.0_8*b4 + 1613976761032884305.0_8/7963657098519931984.0_8*b3 + 564461.0_8/4461432.0_8*b5
              m54 = - 4959271814984644613.0_8/20965546238960637264.0_8*b3 - 1.0_8/6.0_8*b6 - 15998714909649.0_8/37594290333616.0_8*b4 - 375177.0_8/743572.0_8*b5
              m55 = 8386761355510099813.0_8/128413970713633903242.0_8*b3 + 2224717261773437.0_8/2763180339520776.0_8*b4 + 5.0_8/6.0_8*b6 + 1.0_8/24.0_8*b7 + 280535.0_8/371786.0_8*b5
              m56 = - 35039615.0_8/213452232.0_8*b4 - 1.0_8/6.0_8*b7 - 13091810925.0_8/13226425254392.0_8*b3 - 1118749.0_8/2230716.0_8*b5 - 1.0_8/2.0_8*b6
              m57 = -1.0_8/6.0_8*b6 + 1.0_8/8.0_8*b5 + 1.0_8/8.0_8*b7

              d2f_x%p(5,j,k) = -1.0_8/h55*(m51*f%p(1,j,k)+m52*f%p(2,j,k)+m53*f%p(3,j,k) + &
                                          m54*f%p(4,j,k)+m55*f%p(5,j,k)+m56*f%p(6,j,k) + m57*f%p(7,j,k))
          end do
      end if
      if (mesh_xy%is<=6 .and. mesh_xy%ie>=6) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,1,j,k)*mesh_xy%J(1,j,k)
              b2  = mesh_xy%Qi(1,2,j,k)*mesh_xy%J(2,j,k)
              b3  = mesh_xy%Qi(1,3,j,k)*mesh_xy%J(3,j,k)
              b4  = mesh_xy%Qi(1,4,j,k)*mesh_xy%J(4,j,k)
              b5  = mesh_xy%Qi(1,5,j,k)*mesh_xy%J(5,j,k)
              b6  = mesh_xy%Qi(1,6,j,k)*mesh_xy%J(6,j,k)
              b7  = mesh_xy%Qi(1,7,j,k)*mesh_xy%J(7,j,k)
              b8  = mesh_xy%Qi(1,8,j,k)*mesh_xy%J(8,j,k)

              h66 = 1

              m61 = -1.0_8/784.0_8*b4 + 1.0_8/784.0_8*b3
              m62 = -8673.0_8/2904112.0_8*b3 + 8673.0_8/2904112.0_8*b4
              m63 = -960119.0_8/1280713392.0_8*b4 - 3391.0_8/6692148.0_8*b5 + 33235054191.0_8/26452850508784.0_8*b3
              m64 = -368395.0_8/2230716.0_8*b5 + 752806667.0_8/539854092016.0_8*b3 + 1063649.0_8/8712336.0_8*b4 + 1.0_8/8.0_8*b6
              m65 = -35039615.0_8/213452232.0_8*b4 - 1.0_8/6.0_8*b7 - 13091810925.0_8/13226425254392.0_8*b3 - 1118749.0_8/2230716.0_8*b5 - 1.0_8/2.0_8*b6
              m66 =  3290636.0_8/80044587.0_8*b4 + 5580181.0_8/6692148.0_8*b5 + 5.0_8/6.0_8*b7 + 1.0_8/24.0_8*b8 + 660204843.0_8/13226425254392.0_8*b3 + 3.0_8/4.0_8*b6
              m67 = -1.0_8/6.0_8*b5 - 1.0_8/6.0_8*b8 - 1.0_8/2.0_8*b6 - 1.0_8/2.0_8*b7
              m68 = -1.0_8/6.0_8*b7 + 1.0_8/8.0_8*b6 + 1.0_8/8.0_8*b8

              d2f_x%p(6,j,k) = -1.0_8/h66*(m61*f%p(1,j,k)+m62*f%p(2,j,k)+m63*f%p(3,j,k) + &
                                           m64*f%p(4,j,k)+m65*f%p(5,j,k)+m66*f%p(6,j,k) + &
                                           m67*f%p(7,j,k)+m68*f%p(8,j,k) )
          end do
      end if
      do j = mesh_xy%js, mesh_xy%je
          do i = max(mesh_xy%is,7), min(mesh_xy%ie,mesh_xy%nx-6)

              bmm = mesh_xy%Qi(1,i-2,j,k)*mesh_xy%J(i-2,j,k)
              bm  = mesh_xy%Qi(1,i-1,j,k)*mesh_xy%J(i-1,j,k)
              bc  = mesh_xy%Qi(1,i  ,j,k)*mesh_xy%J(i  ,j,k)
              bp  = mesh_xy%Qi(1,i+1,j,k)*mesh_xy%J(i+1,j,k)
              bpp = mesh_xy%Qi(1,i+2,j,k)*mesh_xy%J(i+2,j,k)


              cmm = 1.0_8/6.0_8*bm - 1.0_8/8.0_8*bmm - 1.0_8/8.0_8*bc
              cm  = 1.0_8/6.0_8*bmm + 1.0_8/6.0_8*bp + 1.0_8/2.0_8*bm + 1.0_8/2.0_8*bc
              cc  = -1.0_8/24.0_8*bmm - 5.0_8/6.0_8*bm - 5.0_8/6.0_8*bp - 1.0_8/24.0_8*bpp - 3.0_8/4.0_8*bc
              cp  = 1.0_8/6.0_8*bm + 1.0_8/6.0_8*bpp + 1.0_8/2.0_8*bc + 1.0_8/2.0_8*bp
              cpp = 1.0_8/6.0_8*bp - 1.0_8/8.0_8*bc - 1.0_8/8.0_8*bpp

              d2f_x%p(i,j,k) = (cmm*f%p(i-2,j,k)+cm*f%p(i-1,j,k)+cc*f%p(i,j,k) + &
                                 cp*f%p(i+1,j,k)+cpp*f%p(i+2,j,k))

          end do
      end do
      if (mesh_xy%ie==mesh_xy%nx) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,mesh_xy%nx  ,j,k)*mesh_xy%J(mesh_xy%nx  ,j,k)
              b2  = mesh_xy%Qi(1,mesh_xy%nx-1,j,k)*mesh_xy%J(mesh_xy%nx-1,j,k)
              b3  = mesh_xy%Qi(1,mesh_xy%nx-2,j,k)*mesh_xy%J(mesh_xy%nx-2,j,k)
              b4  = mesh_xy%Qi(1,mesh_xy%nx-3,j,k)*mesh_xy%J(mesh_xy%nx-3,j,k)

              h11 = 17.0_8/48

              m11 = 12.0_8/17.0_8*b1 + 59.0_8/192.0_8*b2 + 27010400129.0_8/345067064608.0_8*b3 + 69462376031.0_8/2070402387648.0_8*b4
              m12 = -59.0_8/68.0_8*b1 -  6025413881.0_8/21126554976.0_8*b3 - 537416663.0_8/7042184992.0_8*b4
              m13 = 2.0_8/17.0_8*b1 - 59.0_8/192.0_8*b2 + 213318005.0_8/16049630912.0_8*b4 + 2083938599.0_8/8024815456.0_8*b3
              m14 = 3.0_8/68.0_8*b1 - 1244724001.0_8/21126554976.0_8*b3 + 752806667.0_8/21126554976.0_8*b4
              m15 = 49579087.0_8/10149031312.0_8*b3 - 49579087.0_8/10149031312.0_8*b4
              m16 = -1.0_8/784.0_8*b4 + 1.0_8/784.0_8*b3

              d2f_x%p(mesh_xy%nx,j,k) = -1.0_8/h11*(m11*f%p(mesh_xy%nx  ,j,k)+m12*f%p(mesh_xy%nx-1,j,k)+m13*f%p(mesh_xy%nx-2,j,k) + &
                                                    m14*f%p(mesh_xy%nx-3,j,k)+m15*f%p(mesh_xy%nx-4,j,k)+m16*f%p(mesh_xy%nx-5,j,k) )
          end do
      end if
      if (mesh_xy%is<=mesh_xy%nx-1 .and. mesh_xy%ie>=mesh_xy%nx-1) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,mesh_xy%nx  ,j,k)*mesh_xy%J(mesh_xy%nx  ,j,k)
              b2  = mesh_xy%Qi(1,mesh_xy%nx-1,j,k)*mesh_xy%J(mesh_xy%nx-1,j,k)
              b3  = mesh_xy%Qi(1,mesh_xy%nx-2,j,k)*mesh_xy%J(mesh_xy%nx-2,j,k)
              b4  = mesh_xy%Qi(1,mesh_xy%nx-3,j,k)*mesh_xy%J(mesh_xy%nx-3,j,k)

              h22 = 59.0_8/48

              m21 = -59.0_8/68.0_8*b1 -  6025413881.0_8/21126554976.0_8*b3 - 537416663.0_8/7042184992.0_8*b4
              m22 = 3481.0_8/3264.0_8*b1 + 9258282831623875.0_8/7669235228057664.0_8*b3 + 236024329996203.0_8/1278205871342944.0_8*b4
              m23 = - 59.0_8/408.0_8*b1 - 29294615794607.0_8/29725717938208.0_8*b3 - 2944673881023.0_8/29725717938208.0_8*b4
              m24 = - 59.0_8/1088.0_8*b1 + 260297319232891.0_8/2556411742685888.0_8*b3 - 60834186813841.0_8/1278205871342944.0_8*b4
              m25 = - 1328188692663.0_8/37594290333616.0_8*b3 + 1328188692663.0_8/37594290333616.0_8*b4
              m26 = - 8673.0_8/2904112.0_8*b3 + 8673.0_8/2904112.0_8*b4

              d2f_x%p(mesh_xy%nx-1,j,k) = -1.0_8/h22*(m21*f%p(mesh_xy%nx  ,j,k)+m22*f%p(mesh_xy%nx-1,j,k)+m23*f%p(mesh_xy%nx-2,j,k) + &
                                                      m24*f%p(mesh_xy%nx-3,j,k)+m25*f%p(mesh_xy%nx-4,j,k)+m26*f%p(mesh_xy%nx-5,j,k) )
          end do
      end if
      if (mesh_xy%is<=mesh_xy%nx-2 .and. mesh_xy%ie>=mesh_xy%nx-2) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,mesh_xy%nx  ,j,k)*mesh_xy%J(mesh_xy%nx  ,j,k)
              b2  = mesh_xy%Qi(1,mesh_xy%nx-1,j,k)*mesh_xy%J(mesh_xy%nx-1,j,k)
              b3  = mesh_xy%Qi(1,mesh_xy%nx-2,j,k)*mesh_xy%J(mesh_xy%nx-2,j,k)
              b4  = mesh_xy%Qi(1,mesh_xy%nx-3,j,k)*mesh_xy%J(mesh_xy%nx-3,j,k)
              b5  = mesh_xy%Qi(1,mesh_xy%nx-4,j,k)*mesh_xy%J(mesh_xy%nx-4,j,k)

              h33 = 43.0_8/48

              m31 = 2.0_8/17.0_8*b1 - 59.0_8/192.0_8*b2 + 213318005.0_8/16049630912.0_8*b4 + 2083938599.0_8/8024815456.0_8*b3
              m32 = - 59.0_8/408.0_8*b1 - 29294615794607.0_8/29725717938208.0_8*b3 - 2944673881023.0_8/29725717938208.0_8*b4
              m33 = 1.0_8/51.0_8*b1 + 59.0_8/192.0_8*b2 + 13777050223300597.0_8/26218083221499456.0_8*b4 + 564461.0_8/13384296.0_8*b5 + 378288882302546512209.0_8/270764341349677687456.0_8*b3
              m34 = 1.0_8/136.0_8*b1 - 125059.0_8/743572.0_8*b5 - 4836340090442187227.0_8/5525802884687299744.0_8*b3 - 17220493277981.0_8/89177153814624.0_8*b4
              m35 = - 10532412077335.0_8/42840005263888.0_8*b4 + 1613976761032884305.0_8/7963657098519931984.0_8*b3 + 564461.0_8/4461432.0_8*b5
              m36 = - 960119.0_8/1280713392.0_8*b4 - 3391.0_8/6692148.0_8*b5 + 33235054191.0_8/26452850508784.0_8*b3

              d2f_x%p(mesh_xy%nx-2,j,k) = -1.0_8/h33*(m31*f%p(mesh_xy%nx  ,j,k)+m32*f%p(mesh_xy%nx-1,j,k)+m33*f%p(mesh_xy%nx-2,j,k) + &
                                                     m34*f%p(mesh_xy%nx-3,j,k)+m35*f%p(mesh_xy%nx-4,j,k)+m36*f%p(mesh_xy%nx-5,j,k) )
          end do
      end if
      if (mesh_xy%is<=mesh_xy%nx-3 .and. mesh_xy%ie>=mesh_xy%nx-3) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,mesh_xy%nx  ,j,k)*mesh_xy%J(mesh_xy%nx  ,j,k)
              b2  = mesh_xy%Qi(1,mesh_xy%nx-1,j,k)*mesh_xy%J(mesh_xy%nx-1,j,k)
              b3  = mesh_xy%Qi(1,mesh_xy%nx-2,j,k)*mesh_xy%J(mesh_xy%nx-2,j,k)
              b4  = mesh_xy%Qi(1,mesh_xy%nx-3,j,k)*mesh_xy%J(mesh_xy%nx-3,j,k)
              b5  = mesh_xy%Qi(1,mesh_xy%nx-4,j,k)*mesh_xy%J(mesh_xy%nx-4,j,k)
              b6  = mesh_xy%Qi(1,mesh_xy%nx-5,j,k)*mesh_xy%J(mesh_xy%nx-5,j,k)

              h44 = 49.0_8/48

              m41 = 3.0_8/68.0_8*b1 - 1244724001.0_8/21126554976.0_8*b3 + 752806667.0_8/21126554976.0_8*b4
              m42 = - 59.0_8/1088.0_8*b1 + 260297319232891.0_8/2556411742685888.0_8*b3 - 60834186813841.0_8/1278205871342944.0_8*b4
              m43 = 1.0_8/136.0_8*b1 - 125059.0_8/743572.0_8*b5 - 4836340090442187227.0_8/5525802884687299744.0_8*b3 - 17220493277981.0_8/89177153814624.0_8*b4
              m44 = 3.0_8/1088.0_8*b1 + 507284006600757858213.0_8/475219048083107777984.0_8*b3 + 1869103.0_8/2230716.0_8*b5 + 1.0_8/24.0_8*b6 + 1950062198436997.0_8/3834617614028832.0_8*b4
              m45 = - 4959271814984644613.0_8/20965546238960637264.0_8*b3 - 1.0_8/6.0_8*b6 - 15998714909649.0_8/37594290333616.0_8*b4 - 375177.0_8/743572.0_8*b5
              m46 = - 368395.0_8/2230716.0_8*b5 + 752806667.0_8/539854092016.0_8*b3 + 1063649.0_8/8712336.0_8*b4 + 1.0_8/8.0_8*b6

              d2f_x%p(mesh_xy%nx-3,j,k) = -1.0_8/h44*(m41*f%p(mesh_xy%nx  ,j,k)+m42*f%p(mesh_xy%nx-1,j,k)+m43*f%p(mesh_xy%nx-2,j,k) + &
                                                      m44*f%p(mesh_xy%nx-3,j,k)+m45*f%p(mesh_xy%nx-4,j,k)+m46*f%p(mesh_xy%nx-5,j,k) )
          end do
      end if
      if (mesh_xy%is<=mesh_xy%nx-4 .and. mesh_xy%ie>=mesh_xy%nx-4) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,mesh_xy%nx  ,j,k)*mesh_xy%J(mesh_xy%nx  ,j,k)
              b2  = mesh_xy%Qi(1,mesh_xy%nx-1,j,k)*mesh_xy%J(mesh_xy%nx-1,j,k)
              b3  = mesh_xy%Qi(1,mesh_xy%nx-2,j,k)*mesh_xy%J(mesh_xy%nx-2,j,k)
              b4  = mesh_xy%Qi(1,mesh_xy%nx-3,j,k)*mesh_xy%J(mesh_xy%nx-3,j,k)
              b5  = mesh_xy%Qi(1,mesh_xy%nx-4,j,k)*mesh_xy%J(mesh_xy%nx-4,j,k)
              b6  = mesh_xy%Qi(1,mesh_xy%nx-5,j,k)*mesh_xy%J(mesh_xy%nx-5,j,k)
              b7  = mesh_xy%Qi(1,mesh_xy%nx-6,j,k)*mesh_xy%J(mesh_xy%nx-6,j,k)

              h55 = 1

              m51 = 49579087.0_8/10149031312.0_8*b3 - 49579087.0_8/10149031312.0_8*b4
              m52 = - 1328188692663.0_8/37594290333616.0_8*b3 + 1328188692663.0_8/37594290333616.0_8*b4
              m53 = - 10532412077335.0_8/42840005263888.0_8*b4 + 1613976761032884305.0_8/7963657098519931984.0_8*b3 + 564461.0_8/4461432.0_8*b5
              m54 = - 4959271814984644613.0_8/20965546238960637264.0_8*b3 - 1.0_8/6.0_8*b6 - 15998714909649.0_8/37594290333616.0_8*b4 - 375177.0_8/743572.0_8*b5
              m55 = 8386761355510099813.0_8/128413970713633903242.0_8*b3 + 2224717261773437.0_8/2763180339520776.0_8*b4 + 5.0_8/6.0_8*b6 + 1.0_8/24.0_8*b7 + 280535.0_8/371786.0_8*b5
              m56 = - 35039615.0_8/213452232.0_8*b4 - 1.0_8/6.0_8*b7 - 13091810925.0_8/13226425254392.0_8*b3 - 1118749.0_8/2230716.0_8*b5 - 1.0_8/2.0_8*b6
              m57 = -1.0_8/6.0_8*b6 + 1.0_8/8.0_8*b5 + 1.0_8/8.0_8*b7

              d2f_x%p(mesh_xy%nx-4,j,k) = -1.0_8/h55*(m51*f%p(mesh_xy%nx  ,j,k)+m52*f%p(mesh_xy%nx-1,j,k)+m53*f%p(mesh_xy%nx-2,j,k) + &
                                                      m54*f%p(mesh_xy%nx-3,j,k)+m55*f%p(mesh_xy%nx-4,j,k)+m56*f%p(mesh_xy%nx-5,j,k) + m57*f%p(mesh_xy%nx-6,j,k))
          end do
      end if
      if (mesh_xy%is<=mesh_xy%nx-5 .and. mesh_xy%ie>=mesh_xy%nx-5) then
          do j = mesh_xy%js, mesh_xy%je

              b1  = mesh_xy%Qi(1,mesh_xy%nx  ,j,k)*mesh_xy%J(mesh_xy%nx  ,j,k)
              b2  = mesh_xy%Qi(1,mesh_xy%nx-1,j,k)*mesh_xy%J(mesh_xy%nx-1,j,k)
              b3  = mesh_xy%Qi(1,mesh_xy%nx-2,j,k)*mesh_xy%J(mesh_xy%nx-2,j,k)
              b4  = mesh_xy%Qi(1,mesh_xy%nx-3,j,k)*mesh_xy%J(mesh_xy%nx-3,j,k)
              b5  = mesh_xy%Qi(1,mesh_xy%nx-4,j,k)*mesh_xy%J(mesh_xy%nx-4,j,k)
              b6  = mesh_xy%Qi(1,mesh_xy%nx-5,j,k)*mesh_xy%J(mesh_xy%nx-5,j,k)
              b7  = mesh_xy%Qi(1,mesh_xy%nx-6,j,k)*mesh_xy%J(mesh_xy%nx-6,j,k)
              b8  = mesh_xy%Qi(1,mesh_xy%nx-7,j,k)*mesh_xy%J(mesh_xy%nx-7,j,k)

              h66 = 1

              m61 = -1.0_8/784.0_8*b4 + 1.0_8/784.0_8*b3
              m62 = - 8673.0_8/2904112.0_8*b3 + 8673.0_8/2904112.0_8*b4
              m63 = - 960119.0_8/1280713392.0_8*b4 - 3391.0_8/6692148.0_8*b5 + 33235054191.0_8/26452850508784.0_8*b3
              m64 = - 368395.0_8/2230716.0_8*b5 + 752806667.0_8/539854092016.0_8*b3 + 1063649.0_8/8712336.0_8*b4 + 1.0_8/8.0_8*b6
              m65 = - 35039615.0_8/213452232.0_8*b4 - 1.0_8/6.0_8*b7 - 13091810925.0_8/13226425254392.0_8*b3 - 1118749.0_8/2230716.0_8*b5 - 1.0_8/2.0_8*b6
              m66 = 3290636.0_8/80044587.0_8*b4 + 5580181.0_8/6692148.0_8*b5 + 5.0_8/6.0_8*b7 + 1.0_8/24.0_8*b8 + 660204843.0_8/13226425254392.0_8*b3 + 3.0_8/4.0_8*b6
              m67 = -1.0_8/6.0_8*b5 - 1.0_8/6.0_8*b8 - 1.0_8/2.0_8*b6 - 1.0_8/2.0_8*b7
              m68 = -1.0_8/6.0_8*b7 + 1.0_8/8.0_8*b6 + 1.0_8/8.0_8*b8

              d2f_x%p(mesh_xy%nx-5,j,k) = -1.0_8/h66*(m61*f%p(mesh_xy%nx  ,j,k)+m62*f%p(mesh_xy%nx-1,j,k)+m63*f%p(mesh_xy%nx-2,j,k) + &
                                                      m64*f%p(mesh_xy%nx-3,j,k)+m65*f%p(mesh_xy%nx-4,j,k)+m66*f%p(mesh_xy%nx-5,j,k) + &
                                                      m67*f%p(mesh_xy%nx-6,j,k)+m68*f%p(mesh_xy%nx-7,j,k))
          end do
      end if

      if (mesh_xy%js==1) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,1,k)*mesh_xy%J(i,1,k)
              b2  = mesh_xy%Qi(3,i,2,k)*mesh_xy%J(i,2,k)
              b3  = mesh_xy%Qi(3,i,3,k)*mesh_xy%J(i,3,k)
              b4  = mesh_xy%Qi(3,i,4,k)*mesh_xy%J(i,4,k)

              h11 = 17.0_8/48

              m11 = 12.0_8/17.0_8*b1 + 59.0_8/192.0_8*b2 + 27010400129.0_8/345067064608.0_8*b3 + 69462376031.0_8/2070402387648.0_8*b4
              m12 = -59.0_8/68.0_8*b1 -  6025413881.0_8/21126554976.0_8*b3 - 537416663.0_8/7042184992.0_8*b4
              m13 = 2.0_8/17.0_8*b1 - 59.0_8/192.0_8*b2 + 213318005.0_8/16049630912.0_8*b4 + 2083938599.0_8/8024815456.0_8*b3
              m14 = 3.0_8/68.0_8*b1 - 1244724001.0_8/21126554976.0_8*b3 + 752806667.0_8/21126554976.0_8*b4
              m15 = 49579087.0_8/10149031312.0_8*b3 - 49579087.0_8/10149031312.0_8*b4
              m16 = -1.0_8/784.0_8*b4 + 1.0_8/784.0_8*b3

              d2f_y%p(i,1,k) = -1.0_8/h11*(m11*f%p(i,1,k)+m12*f%p(i,2,k)+m13*f%p(i,3,k) + &
                                           m14*f%p(i,4,k)+m15*f%p(i,5,k)+m16*f%p(i,6,k) )
          end do
      end if
      if (mesh_xy%js<=2 .and. mesh_xy%je>=2) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,1,k)*mesh_xy%J(i,1,k)
              b2  = mesh_xy%Qi(3,i,2,k)*mesh_xy%J(i,2,k)
              b3  = mesh_xy%Qi(3,i,3,k)*mesh_xy%J(i,3,k)
              b4  = mesh_xy%Qi(3,i,4,k)*mesh_xy%J(i,4,k)

              h22 = 59.0_8/48

              m21 = -59.0_8/68.0_8*b1 -  6025413881.0_8/21126554976.0_8*b3 - 537416663.0_8/7042184992.0_8*b4
              m22 = 3481.0_8/3264.0_8*b1 + 9258282831623875.0_8/7669235228057664.0_8*b3 + 236024329996203.0_8/1278205871342944.0_8*b4
              m23 = - 59.0_8/408.0_8*b1 - 29294615794607.0_8/29725717938208.0_8*b3 - 2944673881023.0_8/29725717938208.0_8*b4
              m24 = - 59.0_8/1088.0_8*b1 + 260297319232891.0_8/2556411742685888.0_8*b3 - 60834186813841.0_8/1278205871342944.0_8*b4
              m25 = - 1328188692663.0_8/37594290333616.0_8*b3 + 1328188692663.0_8/37594290333616.0_8*b4
              m26 = - 8673.0_8/2904112.0_8*b3 + 8673.0_8/2904112.0_8*b4

              d2f_y%p(i,2,k) = -1.0_8/h22*(m21*f%p(i,1,k)+m22*f%p(i,2,k)+m23*f%p(i,3,k) + &
                                           m24*f%p(i,4,k)+m25*f%p(i,5,k)+m26*f%p(i,6,k) )
          end do
      end if
      if (mesh_xy%js<=3 .and. mesh_xy%je>=3) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,1,k)*mesh_xy%J(i,1,k)
              b2  = mesh_xy%Qi(3,i,2,k)*mesh_xy%J(i,2,k)
              b3  = mesh_xy%Qi(3,i,3,k)*mesh_xy%J(i,3,k)
              b4  = mesh_xy%Qi(3,i,4,k)*mesh_xy%J(i,4,k)
              b5  = mesh_xy%Qi(3,i,5,k)*mesh_xy%J(i,5,k)

              h33 = 43.0_8/48

              m31 = 2.0_8/17.0_8*b1 - 59.0_8/192.0_8*b2 + 213318005.0_8/16049630912.0_8*b4 + 2083938599.0_8/8024815456.0_8*b3
              m32 = - 59.0_8/408.0_8*b1 - 29294615794607.0_8/29725717938208.0_8*b3 - 2944673881023.0_8/29725717938208.0_8*b4
              m33 = 1.0_8/51.0_8*b1 + 59.0_8/192.0_8*b2 + 13777050223300597.0_8/26218083221499456.0_8*b4 + 564461.0_8/13384296.0_8*b5 + 378288882302546512209.0_8/270764341349677687456.0_8*b3
              m34 = 1.0_8/136.0_8*b1 - 125059.0_8/743572.0_8*b5 - 4836340090442187227.0_8/5525802884687299744.0_8*b3 - 17220493277981.0_8/89177153814624.0_8*b4
              m35 = - 10532412077335.0_8/42840005263888.0_8*b4 + 1613976761032884305.0_8/7963657098519931984.0_8*b3 + 564461.0_8/4461432.0_8*b5
              m36 = - 960119.0_8/1280713392.0_8*b4 - 3391.0_8/6692148.0_8*b5 + 33235054191.0_8/26452850508784.0_8*b3

              d2f_y%p(i,3,k) = -1.0_8/h33*(m31*f%p(i,1,k)+m32*f%p(i,2,k)+m33*f%p(i,3,k) + &
                                           m34*f%p(i,4,k)+m35*f%p(i,5,k)+m36*f%p(i,6,k) )
          end do
      end if
      if (mesh_xy%js<=4 .and. mesh_xy%je>=4) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,1,k)*mesh_xy%J(i,1,k)
              b2  = mesh_xy%Qi(3,i,2,k)*mesh_xy%J(i,2,k)
              b3  = mesh_xy%Qi(3,i,3,k)*mesh_xy%J(i,3,k)
              b4  = mesh_xy%Qi(3,i,4,k)*mesh_xy%J(i,4,k)
              b5  = mesh_xy%Qi(3,i,5,k)*mesh_xy%J(i,5,k)
              b6  = mesh_xy%Qi(3,i,6,k)*mesh_xy%J(i,6,k)

              h44 = 49.0_8/48

              m41 = 3.0_8/68.0_8*b1 - 1244724001.0_8/21126554976.0_8*b3 + 752806667.0_8/21126554976.0_8*b4
              m42 = - 59.0_8/1088.0_8*b1 + 260297319232891.0_8/2556411742685888.0_8*b3 - 60834186813841.0_8/1278205871342944.0_8*b4
              m43 = 1.0_8/136.0_8*b1 - 125059.0_8/743572.0_8*b5 - 4836340090442187227.0_8/5525802884687299744.0_8*b3 - 17220493277981.0_8/89177153814624.0_8*b4
              m44 = 3.0_8/1088.0_8*b1 + 507284006600757858213.0_8/475219048083107777984.0_8*b3 + 1869103.0_8/2230716.0_8*b5 + 1.0_8/24.0_8*b6 + 1950062198436997.0_8/3834617614028832.0_8*b4
              m45 = - 4959271814984644613.0_8/20965546238960637264.0_8*b3 - 1.0_8/6.0_8*b6 - 15998714909649.0_8/37594290333616.0_8*b4 - 375177.0_8/743572.0_8*b5
              m46 = - 368395.0_8/2230716.0_8*b5 + 752806667.0_8/539854092016.0_8*b3 + 1063649.0_8/8712336.0_8*b4 + 1.0_8/8.0_8*b6

              d2f_y%p(i,4,k) = -1.0_8/h44*(m41*f%p(i,1,k)+m42*f%p(i,2,k)+m43*f%p(i,3,k) + &
                                           m44*f%p(i,4,k)+m45*f%p(i,5,k)+m46*f%p(i,6,k) )
          end do
      end if
      if (mesh_xy%js<=5 .and. mesh_xy%je>=5) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,1,k)*mesh_xy%J(i,1,k)
              b2  = mesh_xy%Qi(3,i,2,k)*mesh_xy%J(i,2,k)
              b3  = mesh_xy%Qi(3,i,3,k)*mesh_xy%J(i,3,k)
              b4  = mesh_xy%Qi(3,i,4,k)*mesh_xy%J(i,4,k)
              b5  = mesh_xy%Qi(3,i,5,k)*mesh_xy%J(i,5,k)
              b6  = mesh_xy%Qi(3,i,6,k)*mesh_xy%J(i,6,k)
              b7  = mesh_xy%Qi(3,i,7,k)*mesh_xy%J(i,7,k)

              h55 = 1

              m51 = 49579087.0_8/10149031312.0_8*b3 - 49579087.0_8/10149031312.0_8*b4
              m52 = - 1328188692663.0_8/37594290333616.0_8*b3 + 1328188692663.0_8/37594290333616.0_8*b4
              m53 = - 10532412077335.0_8/42840005263888.0_8*b4 + 1613976761032884305.0_8/7963657098519931984.0_8*b3 + 564461.0_8/4461432.0_8*b5
              m54 = - 4959271814984644613.0_8/20965546238960637264.0_8*b3 - 1.0_8/6.0_8*b6 - 15998714909649.0_8/37594290333616.0_8*b4 - 375177.0_8/743572.0_8*b5
              m55 = 8386761355510099813.0_8/128413970713633903242.0_8*b3 + 2224717261773437.0_8/2763180339520776.0_8*b4 + 5.0_8/6.0_8*b6 + 1.0_8/24.0_8*b7 + 280535.0_8/371786.0_8*b5
              m56 = - 35039615.0_8/213452232.0_8*b4 - 1.0_8/6.0_8*b7 - 13091810925.0_8/13226425254392.0_8*b3 - 1118749.0_8/2230716.0_8*b5 - 1.0_8/2.0_8*b6
              m57 = -1.0_8/6.0_8*b6 + 1.0_8/8.0_8*b5 + 1.0_8/8.0_8*b7

              d2f_y%p(i,5,k) = -1.0_8/h55*(m51*f%p(i,1,k)+m52*f%p(i,2,k)+m53*f%p(i,3,k) + &
                                           m54*f%p(i,4,k)+m55*f%p(i,5,k)+m56*f%p(i,6,k) + m57*f%p(i,7,k))
          end do
      end if
      if (mesh_xy%js<=6 .and. mesh_xy%je>=6) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,1,k)*mesh_xy%J(i,1,k)
              b2  = mesh_xy%Qi(3,i,2,k)*mesh_xy%J(i,2,k)
              b3  = mesh_xy%Qi(3,i,3,k)*mesh_xy%J(i,3,k)
              b4  = mesh_xy%Qi(3,i,4,k)*mesh_xy%J(i,4,k)
              b5  = mesh_xy%Qi(3,i,5,k)*mesh_xy%J(i,5,k)
              b6  = mesh_xy%Qi(3,i,6,k)*mesh_xy%J(i,6,k)
              b7  = mesh_xy%Qi(3,i,7,k)*mesh_xy%J(i,7,k)
              b8  = mesh_xy%Qi(3,i,8,k)*mesh_xy%J(i,8,k)

              h66 = 1

              m61 = -1.0_8/784.0_8*b4 + 1.0_8/784.0_8*b3
              m62 = - 8673.0_8/2904112.0_8*b3 + 8673.0_8/2904112.0_8*b4
              m63 = - 960119.0_8/1280713392.0_8*b4 - 3391.0_8/6692148.0_8*b5 + 33235054191.0_8/26452850508784.0_8*b3
              m64 = - 368395.0_8/2230716.0_8*b5 + 752806667.0_8/539854092016.0_8*b3 + 1063649.0_8/8712336.0_8*b4 + 1.0_8/8.0_8*b6
              m65 = - 35039615.0_8/213452232.0_8*b4 - 1.0_8/6.0_8*b7 - 13091810925.0_8/13226425254392.0_8*b3 - 1118749.0_8/2230716.0_8*b5 - 1.0_8/2.0_8*b6
              m66 = 3290636.0_8/80044587.0_8*b4 + 5580181.0_8/6692148.0_8*b5 + 5.0_8/6.0_8*b7 + 1.0_8/24.0_8*b8 + 660204843.0_8/13226425254392.0_8*b3 + 3.0_8/4.0_8*b6
              m67 = -1.0_8/6.0_8*b5 - 1.0_8/6.0_8*b8 - 1.0_8/2.0_8*b6 - 1.0_8/2.0_8*b7
              m68 = -1.0_8/6.0_8*b7 + 1.0_8/8.0_8*b6 + 1.0_8/8.0_8*b8

              d2f_y%p(i,6,k) = -1.0_8/h66*(m61*f%p(i,1,k)+m62*f%p(i,2,k)+m63*f%p(i,3,k) + &
                                           m64*f%p(i,4,k)+m65*f%p(i,5,k)+m66*f%p(i,6,k) + &
                                           m67*f%p(i,7,k)+m68*f%p(i,8,k))
          end do
      end if
      do j = max(mesh_xy%js,7), min(mesh_xy%je,mesh_xy%ny-6)
          do i = mesh_xy%is, mesh_xy%ie

              bmm = mesh_xy%Qi(3,i,j-2,k)*mesh_xy%J(i,j-2,k)
              bm  = mesh_xy%Qi(3,i,j-1,k)*mesh_xy%J(i,j-1,k)
              bc  = mesh_xy%Qi(3,i,j  ,k)*mesh_xy%J(i,j  ,k)
              bp  = mesh_xy%Qi(3,i,j+1,k)*mesh_xy%J(i,j+1,k)
              bpp = mesh_xy%Qi(3,i,j+2,k)*mesh_xy%J(i,j+2,k)


              cmm = 1.0_8/6.0_8*bm - 1.0_8/8.0_8*bmm - 1.0_8/8.0_8*bc
              cm  = 1.0_8/6.0_8*bmm + 1.0_8/6.0_8*bp + 1.0_8/2.0_8*bm + 1.0_8/2.0_8*bc
              cc  = -1.0_8/24.0_8*bmm - 5.0_8/6.0_8*bm - 5.0_8/6.0_8*bp - 1.0_8/24.0_8*bpp - 3.0_8/4.0_8*bc
              cp  = 1.0_8/6.0_8*bm + 1.0_8/6.0_8*bpp + 1.0_8/2.0_8*bc + 1.0_8/2.0_8*bp
              cpp = 1.0_8/6.0_8*bp - 1.0_8/8.0_8*bc - 1.0_8/8.0_8*bpp

              d2f_y%p(i,j,k) = (cmm*f%p(i,j-2,k)+cm*f%p(i,j-1,k)+cc*f%p(i,j,k) + &
                                 cp*f%p(i,j+1,k)+cpp*f%p(i,j+2,k))
          end do
      end do
      if (mesh_xy%je==mesh_xy%ny) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,mesh_xy%ny  ,k)*mesh_xy%J(i,mesh_xy%ny  ,k)
              b2  = mesh_xy%Qi(3,i,mesh_xy%ny-1,k)*mesh_xy%J(i,mesh_xy%ny-1,k)
              b3  = mesh_xy%Qi(3,i,mesh_xy%ny-2,k)*mesh_xy%J(i,mesh_xy%ny-2,k)
              b4  = mesh_xy%Qi(3,i,mesh_xy%ny-3,k)*mesh_xy%J(i,mesh_xy%ny-3,k)

              h11 = 17.0_8/48

              m11 = 12.0_8/17.0_8*b1 + 59.0_8/192.0_8*b2 + 27010400129.0_8/345067064608.0_8*b3 + 69462376031.0_8/2070402387648.0_8*b4
              m12 = -59.0_8/68.0_8*b1 -  6025413881.0_8/21126554976.0_8*b3 - 537416663.0_8/7042184992.0_8*b4
              m13 = 2.0_8/17.0_8*b1 - 59.0_8/192.0_8*b2 + 213318005.0_8/16049630912.0_8*b4 + 2083938599.0_8/8024815456.0_8*b3
              m14 = 3.0_8/68.0_8*b1 - 1244724001.0_8/21126554976.0_8*b3 + 752806667.0_8/21126554976.0_8*b4
              m15 = 49579087.0_8/10149031312.0_8*b3 - 49579087.0_8/10149031312.0_8*b4
              m16 = -1.0_8/784.0_8*b4 + 1.0_8/784.0_8*b3

              d2f_y%p(i,mesh_xy%ny,k) = -1.0_8/h11*(m11*f%p(i,mesh_xy%ny  ,k)+m12*f%p(i,mesh_xy%ny-1,k)+m13*f%p(i,mesh_xy%ny-2,k) + &
                                                    m14*f%p(i,mesh_xy%ny-3,k)+m15*f%p(i,mesh_xy%ny-4,k)+m16*f%p(i,mesh_xy%ny-5,k) )
          end do
      end if
      if (mesh_xy%js<=mesh_xy%ny-1 .and. mesh_xy%je>=mesh_xy%ny-1) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,mesh_xy%ny  ,k)*mesh_xy%J(i,mesh_xy%ny  ,k)
              b2  = mesh_xy%Qi(3,i,mesh_xy%ny-1,k)*mesh_xy%J(i,mesh_xy%ny-1,k)
              b3  = mesh_xy%Qi(3,i,mesh_xy%ny-2,k)*mesh_xy%J(i,mesh_xy%ny-2,k)
              b4  = mesh_xy%Qi(3,i,mesh_xy%ny-3,k)*mesh_xy%J(i,mesh_xy%ny-3,k)

              h22 = 59.0_8/48

              m21 = -59.0_8/68.0_8*b1 -  6025413881.0_8/21126554976.0_8*b3 - 537416663.0_8/7042184992.0_8*b4
              m22 = 3481.0_8/3264.0_8*b1 + 9258282831623875.0_8/7669235228057664.0_8*b3 + 236024329996203.0_8/1278205871342944.0_8*b4
              m23 = - 59.0_8/408.0_8*b1 - 29294615794607.0_8/29725717938208.0_8*b3 - 2944673881023.0_8/29725717938208.0_8*b4
              m24 = - 59.0_8/1088.0_8*b1 + 260297319232891.0_8/2556411742685888.0_8*b3 - 60834186813841.0_8/1278205871342944.0_8*b4
              m25 = - 1328188692663.0_8/37594290333616.0_8*b3 + 1328188692663.0_8/37594290333616.0_8*b4
              m26 = - 8673.0_8/2904112.0_8*b3 + 8673.0_8/2904112.0_8*b4

              d2f_y%p(i,mesh_xy%ny-1,k) = -1.0_8/h22*(m21*f%p(i,mesh_xy%ny  ,k)+m22*f%p(i,mesh_xy%ny-1,k)+m23*f%p(i,mesh_xy%ny-2,k) + &
                                                      m24*f%p(i,mesh_xy%ny-3,k)+m25*f%p(i,mesh_xy%ny-4,k)+m26*f%p(i,mesh_xy%ny-5,k) )
          end do
      end if
      if (mesh_xy%js<=mesh_xy%ny-2 .and. mesh_xy%je>=mesh_xy%ny-2) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,mesh_xy%ny  ,k)*mesh_xy%J(i,mesh_xy%ny  ,k)
              b2  = mesh_xy%Qi(3,i,mesh_xy%ny-1,k)*mesh_xy%J(i,mesh_xy%ny-1,k)
              b3  = mesh_xy%Qi(3,i,mesh_xy%ny-2,k)*mesh_xy%J(i,mesh_xy%ny-2,k)
              b4  = mesh_xy%Qi(3,i,mesh_xy%ny-3,k)*mesh_xy%J(i,mesh_xy%ny-3,k)
              b5  = mesh_xy%Qi(3,i,mesh_xy%ny-4,k)*mesh_xy%J(i,mesh_xy%ny-4,k)

              h33 = 43.0_8/48

              m31 = 2.0_8/17.0_8*b1 - 59.0_8/192.0_8*b2 + 213318005.0_8/16049630912.0_8*b4 + 2083938599.0_8/8024815456.0_8*b3
              m32 = - 59.0_8/408.0_8*b1 - 29294615794607.0_8/29725717938208.0_8*b3 - 2944673881023.0_8/29725717938208.0_8*b4
              m33 = 1.0_8/51.0_8*b1 + 59.0_8/192.0_8*b2 + 13777050223300597.0_8/26218083221499456.0_8*b4 + 564461.0_8/13384296.0_8*b5 + 378288882302546512209.0_8/270764341349677687456.0_8*b3
              m34 = 1.0_8/136.0_8*b1 - 125059.0_8/743572.0_8*b5 - 4836340090442187227.0_8/5525802884687299744.0_8*b3 - 17220493277981.0_8/89177153814624.0_8*b4
              m35 = - 10532412077335.0_8/42840005263888.0_8*b4 + 1613976761032884305.0_8/7963657098519931984.0_8*b3 + 564461.0_8/4461432.0_8*b5
              m36 = - 960119.0_8/1280713392.0_8*b4 - 3391.0_8/6692148.0_8*b5 + 33235054191.0_8/26452850508784.0_8*b3

              d2f_y%p(i,mesh_xy%ny-2,k) = -1.0_8/h33*(m31*f%p(i,mesh_xy%ny  ,k)+m32*f%p(i,mesh_xy%ny-1,k)+m33*f%p(i,mesh_xy%ny-2,k) + &
                                                      m34*f%p(i,mesh_xy%ny-3,k)+m35*f%p(i,mesh_xy%ny-4,k)+m36*f%p(i,mesh_xy%ny-5,k) )
          end do
      end if
      if (mesh_xy%js<=mesh_xy%ny-3 .and. mesh_xy%je>=mesh_xy%ny-3) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,mesh_xy%ny  ,k)*mesh_xy%J(i,mesh_xy%ny  ,k)
              b2  = mesh_xy%Qi(3,i,mesh_xy%ny-1,k)*mesh_xy%J(i,mesh_xy%ny-1,k)
              b3  = mesh_xy%Qi(3,i,mesh_xy%ny-2,k)*mesh_xy%J(i,mesh_xy%ny-2,k)
              b4  = mesh_xy%Qi(3,i,mesh_xy%ny-3,k)*mesh_xy%J(i,mesh_xy%ny-3,k)
              b5  = mesh_xy%Qi(3,i,mesh_xy%ny-4,k)*mesh_xy%J(i,mesh_xy%ny-4,k)
              b6  = mesh_xy%Qi(3,i,mesh_xy%ny-5,k)*mesh_xy%J(i,mesh_xy%ny-5,k)

              h44 = 49.0_8/48

              m41 = 3.0_8/68.0_8*b1 - 1244724001.0_8/21126554976.0_8*b3 + 752806667.0_8/21126554976.0_8*b4
              m42 = - 59.0_8/1088.0_8*b1 + 260297319232891.0_8/2556411742685888.0_8*b3 - 60834186813841.0_8/1278205871342944.0_8*b4
              m43 = 1.0_8/136.0_8*b1 - 125059.0_8/743572.0_8*b5 - 4836340090442187227.0_8/5525802884687299744.0_8*b3 - 17220493277981.0_8/89177153814624.0_8*b4
              m44 = 3.0_8/1088.0_8*b1 + 507284006600757858213.0_8/475219048083107777984.0_8*b3 + 1869103.0_8/2230716.0_8*b5 + 1.0_8/24.0_8*b6 + 1950062198436997.0_8/3834617614028832.0_8*b4
              m45 = - 4959271814984644613.0_8/20965546238960637264.0_8*b3 - 1.0_8/6.0_8*b6 - 15998714909649.0_8/37594290333616.0_8*b4 - 375177.0_8/743572.0_8*b5
              m46 = - 368395.0_8/2230716.0_8*b5 + 752806667.0_8/539854092016.0_8*b3 + 1063649.0_8/8712336.0_8*b4 + 1.0_8/8.0_8*b6

              d2f_y%p(i,mesh_xy%ny-3,k) = -1.0_8/h44*(m41*f%p(i,mesh_xy%ny  ,k)+m42*f%p(i,mesh_xy%ny-1,k)+m43*f%p(i,mesh_xy%ny-2,k) + &
                                                     m44*f%p(i,mesh_xy%ny-3,k)+m45*f%p(i,mesh_xy%ny-4,k)+m46*f%p(i,mesh_xy%ny-5,k) )
          end do
      end if
      if (mesh_xy%js<=mesh_xy%ny-4 .and. mesh_xy%je>=mesh_xy%ny-4) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,mesh_xy%ny  ,k)*mesh_xy%J(i,mesh_xy%ny  ,k)
              b2  = mesh_xy%Qi(3,i,mesh_xy%ny-1,k)*mesh_xy%J(i,mesh_xy%ny-1,k)
              b3  = mesh_xy%Qi(3,i,mesh_xy%ny-2,k)*mesh_xy%J(i,mesh_xy%ny-2,k)
              b4  = mesh_xy%Qi(3,i,mesh_xy%ny-3,k)*mesh_xy%J(i,mesh_xy%ny-3,k)
              b5  = mesh_xy%Qi(3,i,mesh_xy%ny-4,k)*mesh_xy%J(i,mesh_xy%ny-4,k)
              b6  = mesh_xy%Qi(3,i,mesh_xy%ny-5,k)*mesh_xy%J(i,mesh_xy%ny-5,k)
              b7  = mesh_xy%Qi(3,i,mesh_xy%ny-6,k)*mesh_xy%J(i,mesh_xy%ny-6,k)

              h55 = 1

              m51 = 49579087.0_8/10149031312.0_8*b3 - 49579087.0_8/10149031312.0_8*b4
              m52 = - 1328188692663.0_8/37594290333616.0_8*b3 + 1328188692663.0_8/37594290333616.0_8*b4
              m53 = - 10532412077335.0_8/42840005263888.0_8*b4 + 1613976761032884305.0_8/7963657098519931984.0_8*b3 + 564461.0_8/4461432.0_8*b5
              m54 = - 4959271814984644613.0_8/20965546238960637264.0_8*b3 - 1.0_8/6.0_8*b6 - 15998714909649.0_8/37594290333616.0_8*b4 - 375177.0_8/743572.0_8*b5
              m55 = 8386761355510099813.0_8/128413970713633903242.0_8*b3 + 2224717261773437.0_8/2763180339520776.0_8*b4 + 5.0_8/6.0_8*b6 + 1.0_8/24.0_8*b7 + 280535.0_8/371786.0_8*b5
              m56 = - 35039615.0_8/213452232.0_8*b4 - 1.0_8/6.0_8*b7 - 13091810925.0_8/13226425254392.0_8*b3 - 1118749.0_8/2230716.0_8*b5 - 1.0_8/2.0_8*b6
              m57 = -1.0_8/6.0_8*b6 + 1.0_8/8.0_8*b5 + 1.0_8/8.0_8*b7

              d2f_y%p(i,mesh_xy%ny-4,k) = -1.0_8/h55*(m51*f%p(i,mesh_xy%ny  ,k)+m52*f%p(i,mesh_xy%ny-1,k)+m53*f%p(i,mesh_xy%ny-2,k) + &
                                                      m54*f%p(i,mesh_xy%ny-3,k)+m55*f%p(i,mesh_xy%ny-4,k)+m56*f%p(i,mesh_xy%ny-5,k) + m57*f%p(i,mesh_xy%ny-6,k))
          end do
      end if
      if (mesh_xy%js<=mesh_xy%ny-5 .and. mesh_xy%je>=mesh_xy%ny-5) then
          do i = mesh_xy%is, mesh_xy%ie
              b1  = mesh_xy%Qi(3,i,mesh_xy%ny  ,k)*mesh_xy%J(i,mesh_xy%ny  ,k)
              b2  = mesh_xy%Qi(3,i,mesh_xy%ny-1,k)*mesh_xy%J(i,mesh_xy%ny-1,k)
              b3  = mesh_xy%Qi(3,i,mesh_xy%ny-2,k)*mesh_xy%J(i,mesh_xy%ny-2,k)
              b4  = mesh_xy%Qi(3,i,mesh_xy%ny-3,k)*mesh_xy%J(i,mesh_xy%ny-3,k)
              b5  = mesh_xy%Qi(3,i,mesh_xy%ny-4,k)*mesh_xy%J(i,mesh_xy%ny-4,k)
              b6  = mesh_xy%Qi(3,i,mesh_xy%ny-5,k)*mesh_xy%J(i,mesh_xy%ny-5,k)
              b7  = mesh_xy%Qi(3,i,mesh_xy%ny-6,k)*mesh_xy%J(i,mesh_xy%ny-6,k)
              b8  = mesh_xy%Qi(3,i,mesh_xy%ny-7,k)*mesh_xy%J(i,mesh_xy%ny-7,k)

              h66 = 1

              m61 = -1.0_8/784.0_8*b4 + 1.0_8/784.0_8*b3
              m62 = - 8673.0_8/2904112.0_8*b3 + 8673.0_8/2904112.0_8*b4
              m63 = - 960119.0_8/1280713392.0_8*b4 - 3391.0_8/6692148.0_8*b5 + 33235054191.0_8/26452850508784.0_8*b3
              m64 = - 368395.0_8/2230716.0_8*b5 + 752806667.0_8/539854092016.0_8*b3 + 1063649.0_8/8712336.0_8*b4 + 1.0_8/8.0_8*b6
              m65 = - 35039615.0_8/213452232.0_8*b4 - 1.0_8/6.0_8*b7 - 13091810925.0_8/13226425254392.0_8*b3 - 1118749.0_8/2230716.0_8*b5 - 1.0_8/2.0_8*b6
              m66 = 3290636.0_8/80044587.0_8*b4 + 5580181.0_8/6692148.0_8*b5 + 5.0_8/6.0_8*b7 + 1.0_8/24.0_8*b8 + 660204843.0_8/13226425254392.0_8*b3 + 3.0_8/4.0_8*b6
              m67 = -1.0_8/6.0_8*b5 - 1.0_8/6.0_8*b8 - 1.0_8/2.0_8*b6 - 1.0_8/2.0_8*b7
              m68 = -1.0_8/6.0_8*b7 + 1.0_8/8.0_8*b6 + 1.0_8/8.0_8*b8

              d2f_y%p(i,mesh_xy%ny-5,k) = -1.0_8/h66*(m61*f%p(i,mesh_xy%ny  ,k)+m62*f%p(i,mesh_xy%ny-1,k)+m63*f%p(i,mesh_xy%ny-2,k) + &
                                                      m64*f%p(i,mesh_xy%ny-3,k)+m65*f%p(i,mesh_xy%ny-4,k)+m66*f%p(i,mesh_xy%ny-5,k) + &
                                                      m67*f%p(i,mesh_xy%ny-6,k)+m68*f%p(i,mesh_xy%ny-7,k) )
          end do
      end if
    end do

    do k = mesh_xy%ks, mesh_xy%ke
      if (mesh_xy%is==1) then
          do j = mesh_xy%js, mesh_xy%je
              penalty = 48.0_8/17*(-2.0_8*d1f_y%p(1,j,k))
              d2f_x%p(1,j,k) = d2f_x%p(1,j,k)-0.5_8*penalty
          end do
      end if
      if (mesh_xy%ie==mesh_xy%nx) then
          do j = mesh_xy%js, mesh_xy%je
              penalty = 48.0_8/17*(2*d1f_y%p(mesh_xy%nx,j,k))
              d2f_x%p(mesh_xy%nx,j,k) = d2f_x%p(mesh_xy%nx,j,k) -0.5_8*penalty
          end do
      end if

      if (mesh_xy%js==1) then
          do i = mesh_xy%is, mesh_xy%ie
              penalty = 48.0_8/17*(-2*d1f_x%p(i,1,k))
              d2f_x%p(i,1,k) = d2f_x%p(i,1,k) -0.5_8*penalty
          end do
      end if
      if (mesh_xy%je==mesh_xy%ny) then
          do i = mesh_xy%is, mesh_xy%ie
              penalty = 48.0_8/17*(2*d1f_x%p(i,mesh_xy%ny,k))
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

end module laplace_Ah_sbp42_narrow_mod

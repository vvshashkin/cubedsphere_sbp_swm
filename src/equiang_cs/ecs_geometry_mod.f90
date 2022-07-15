!Module for equiangular cubed sphere geometrical parameters
module ecs_geometry_mod
implicit none

contains

function ecs_ab2xyz_proto(alpha,beta) result(xyz)
    real(kind=8)                   xyz(3)
    real(kind=8), intent(in)    :: alpha, beta
    !locals:
    real(kind=8) r
    integer(kind=4) i

    !grid point at proto cube face: (assumes x = tan(alpha), y = tan(beta), z = 1)
    r = sqrt(1._8+tan(alpha)**2+tan(beta)**2)
    !coordinates at prototype spherical face: vec(r)/||r||, r = (x, y, z)
    xyz(1) = tan(alpha)/r; xyz(2) = tan(beta)/r; xyz(3) = 1._8/r

end function ecs_ab2xyz_proto

function ecs_acov_proto(alpha,beta) result(acov)
    !covariant a-vector at prototype face given alpha&beta
    real(kind=8)              :: acov(3)
    real(kind=8), intent(in)  :: alpha, beta
    !local
    real(kind=8) zta, ztb

    zta = tan(alpha);   ztb = tan(beta)
    !d (xyz)^T/ d alpha
    acov(1) = (1d0+zta**2d0)*(1d0+ztb**2d0)/(1d0+zta**2d0+ztb**2d0)**(3d0/2d0)
    acov(2) = -zta*ztb*(1d0+zta**2d0)/(1d0+zta**2d0+ztb**2d0)**(3d0/2d0)
    acov(3) = -zta*(1d0+zta**2d0) / (1d0+zta**2d0+ztb**2d0)**(3d0/2d0)

end function ecs_acov_proto

function ecs_bcov_proto(alpha,beta) result(bcov)
    !covariant b-vector at prototype face given alpha&beta
    real(kind=8)              :: bcov(3)
    real(kind=8), intent(in)  :: alpha, beta
    !local
    real(kind=8) zta, ztb

    zta = tan(alpha);   ztb = tan(beta)
    !d (xyz)^T/ d beta
    bcov(1) = -zta*ztb*(1d0+ztb**2d0)/(1d0+zta**2d0+ztb**2d0)**(3d0/2d0)
    bcov(2) = (1d0+zta**2d0)*(1d0+ztb**2d0)/(1d0+zta**2d0+ztb**2d0)**(3d0/2d0)
    bcov(3) = -ztb*(1d0+ztb**2d0) / (1d0+zta**2d0+ztb**2d0)**(3d0/2d0)

end function ecs_bcov_proto

function ecs_actv_proto(alpha,beta) result(actv)
    !covariant a-vector at prototype face given alpha&beta
    real(kind=8)              :: actv(3)
    real(kind=8), intent(in)  :: alpha, beta
    !local
    real(kind=8) zta, ztb, zsigm

    zta = tan(alpha);   ztb = tan(beta)
    zsigm = sqrt(1._8+zta**2+ztb**2)
    !d (xyz)^T/ d alpha
    actv(1) = (1._8+ztb**2-(zta*ztb)**2/(1._8+zta**2)) / zsigm
    actv(2) = 0._8
    actv(3) =-zta*(1._8+ztb**2/(1+zta**2))/zsigm

end function ecs_actv_proto

function ecs_bctv_proto(alpha,beta) result(bctv)
    !covariant b-vector at prototype face given alpha&beta
    real(kind=8)              :: bctv(3)
    real(kind=8), intent(in)  :: alpha, beta
    !local
    real(kind=8) zta, ztb, zsigm

    zta = tan(alpha);   ztb = tan(beta)
    zsigm = sqrt(1._8+zta**2+ztb**2)
    !d (xyz)^T/ d beta
    bctv(1) = 0._8
    bctv(2) = ((1+zta**2)-(zta*ztb)**2/(1+ztb**2))/zsigm
    bctv(3) =-ztb*(zta**2/(1+ztb**2)+1._8)/zsigm

end function ecs_bctv_proto

function ecs_metric_tensor(alpha,beta) result(Q)
    !equiangular cubsphere metric tensor (only
    !3-elem stored due to symmetricity)
    real(kind=8)              :: Q(3) !Q = |q11 q12| = |Q(1) Q(2)|
                                      !    |q21 q22|   |Q(2) Q(3)|
    real(kind=8), intent(in)  :: alpha, beta
    !local
    real(kind=8) zta, ztb, zsigm

    zta = tan(alpha);   ztb = tan(beta)
    zsigm = sqrt(1._8+zta**2+ztb**2)

    Q(1) = (1._8+zta**2)**2*(1._8+ztb**2)/zsigm**4
    Q(2) =-(1._8+zta**2)*(1._8+ztb**2)*zta*ztb/zsigm**4
    Q(3) = (1._8+zta**2)*(1._8+ztb**2)**2/zsigm**4

end function ecs_metric_tensor

function ecs_invmetric_tensor(alpha,beta) result(Qi)
    !equiangular cubsphere inverse metric tensor (only
    !3-elem stored due to symmetricity)
    real(kind=8)              :: Qi(3) !Qi = |q11 q12| = |Q(1) Q(2)|
                                       !     |q21 q22|   |Q(2) Q(3)|
    real(kind=8), intent(in)  :: alpha, beta
    !local
    real(kind=8) zta, ztb, zsigm

    zta = tan(alpha);   ztb = tan(beta)
    zsigm = sqrt(1._8+zta**2+ztb**2)

    Qi(1) = zsigm**2 / (1._8+zta**2)
    Qi(2) = zsigm**2*zta*ztb / ((1+zta**2)*(1+ztb**2))
    Qi(3) = zsigm**2 / (1._8+ztb**2)

end function ecs_invmetric_tensor

function ecs_G(alpha,beta) result(G)
    !equiangular cubsphere sqrt of metric tensor determinant
    !area of dalpha*dbeta element = G*dalpha*dbeta +O(h**2)
    real(kind=8)              :: G
    real(kind=8), intent(in)  :: alpha, beta
    !local
    real(kind=8) zta, ztb, zsigm

    zta = tan(alpha);   ztb = tan(beta)
    zsigm = sqrt(1._8+zta**2+ztb**2)

    G = (1._8+zta**2)*(1._8+ztb**2) / zsigm**3

end function ecs_G

function ecs_proto2face(xyz,panel_ind) result(xyz1)
    !transform prototype-face xyz to real spherical face
    use topology_mod, only : ex, ey, n
    real(kind=8)                :: xyz1(3)
    real(kind=8), intent(in)    :: xyz(3)
    integer(kind=4), intent(in) :: panel_ind

    xyz1(1) = ex(1,panel_ind)*xyz(1)+ey(1,panel_ind)*xyz(2)-n(1,panel_ind)*xyz(3)
    xyz1(2) = ex(2,panel_ind)*xyz(1)+ey(2,panel_ind)*xyz(2)-n(2,panel_ind)*xyz(3)
    xyz1(3) = ex(3,panel_ind)*xyz(1)+ey(3,panel_ind)*xyz(2)-n(3,panel_ind)*xyz(3)

end function ecs_proto2face

subroutine ecs_xyz2ab(pa,pb,px,py,pz,ifc)
use topology_mod, only : ex, ey, n
!transform x,y,z to a,b of panel number ifc
!potentially dangerous: ifc-a,b are not defined for some xyz
real(kind=8), intent(out)   :: pa, pb
real(kind=8), intent(in)    :: px, py, pz
integer(kind=4), intent(in) :: ifc
!locals
real(kind=8) zx, zy, zz

!transform to prototype-face coordinates
zx = ex(1,ifc)*px+ex(2,ifc)*py+ex(3,ifc)*pz
zy = ey(1,ifc)*px+ey(2,ifc)*py+ey(3,ifc)*pz
zz = -n(1,ifc)*px- n(2,ifc)*py- n(3,ifc)*pz

pa = atan(zx/zz);   pb = atan(zy/zz)

end subroutine ecs_xyz2ab

end module ecs_geometry_mod

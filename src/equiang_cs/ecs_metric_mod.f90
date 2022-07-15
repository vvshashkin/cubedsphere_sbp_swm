!Module for equiangular cubed sphere metric parameters
module ecs_metric_mod

use metric_mod,                only : metric_t
use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
use parcomm_mod,               only : parcomm_global

implicit none

type, extends(metric_t) :: ecs_metric_t
    class(cubed_sphere_topology_t), allocatable :: topology
contains
    !New interface
    procedure :: calculate_r_orog  => calculate_ecs_r_orog
    procedure :: calculate_r_2d    => calculate_ecs_r_2d

    procedure :: calculate_h       => return_zero_h

    procedure :: calculate_a1_orog => calculate_ecs_a1_orog
    procedure :: calculate_a1_2d   => calculate_ecs_a1_2d
    procedure :: calculate_a2_orog => calculate_ecs_a2_orog
    procedure :: calculate_a2_2d   => calculate_ecs_a2_2d
    procedure :: calculate_a3_orog => calculate_ecs_a3_orog
    procedure :: calculate_a3_2d   => calculate_ecs_a3_2d

    procedure :: calculate_b1_orog => calculate_ecs_b1_orog
    procedure :: calculate_b1_2d   => calculate_ecs_b1_2d
    procedure :: calculate_b2_orog => calculate_ecs_b2_orog
    procedure :: calculate_b2_2d   => calculate_ecs_b2_2d
    procedure :: calculate_b3_orog => calculate_ecs_b3_orog
    procedure :: calculate_b3_2d   => calculate_ecs_b3_2d

    procedure :: calculate_Q_orog  => calculate_ecs_Q_orog
    procedure :: calculate_Q_2d    => calculate_ecs_Q_2d
    procedure :: calculate_Qi_orog => calculate_ecs_Qi_orog
    procedure :: calculate_Qi_2d   => calculate_ecs_Qi_2d

    procedure :: calculate_J_orog => calculate_ecs_J_orog
    procedure :: calculate_J_2d   => calculate_ecs_J_2d

    procedure :: calculate_G_orog => calculate_ecs_G_orog
    procedure :: calculate_G_2d   => calculate_ecs_G_2d

    procedure :: transform_cartesian_to_native => transform_cartesian_to_native_ecs

end type ecs_metric_t

contains

pure function calculate_ecs_r_orog(this, panel_ind, alpha, beta, eta, h_surf, h_top) result(r)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, h_top
    real(kind=8)                    :: r(3)

    r = calculate_ecs_r_2d(this, panel_ind, alpha, beta)
end function calculate_ecs_r_orog

pure function calculate_ecs_r_2d(this, panel_ind, alpha, beta) result(r)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: r(3)

    r = ecs_point_r_proto(alpha,beta)
    r = ecs_proto2realface(this%topology,this%rotation_matrix,panel_ind,r)

end function calculate_ecs_r_2d

pure function return_zero_h(this,panel_ind,alpha,beta,eta,h_surf,h_top) result(h)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, h_top
    real(kind=8)                    :: h

    h = 0.0_8
end function return_zero_h

pure function calculate_ecs_a1_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dcov_h_surf, h_top) result(a)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dcov_h_surf, h_top
    real(kind=8)                    :: a(4)

    a = calculate_ecs_a1_2d(this, panel_ind, alpha, beta)

end function calculate_ecs_a1_orog

pure function calculate_ecs_a1_2d(this, panel_ind, alpha, beta) result(a)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: a(4)

    a(1:3) = ecs_a1_proto(alpha,beta)
    a(1:3) = ecs_proto2realface(this%topology,this%rotation_matrix,panel_ind,a(1:3))
    a(4) = 0._8
end function calculate_ecs_a1_2d

pure function calculate_ecs_a2_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dcov_h_surf, h_top) result(a)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dcov_h_surf, h_top
    real(kind=8)                    :: a(4)

    a = calculate_ecs_a2_2d(this, panel_ind, alpha, beta)

end function calculate_ecs_a2_orog

pure function calculate_ecs_a2_2d(this, panel_ind, alpha, beta) result(a)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: a(4)

    a(1:3) = ecs_a2_proto(alpha,beta)
    a(1:3) = ecs_proto2realface(this%topology,this%rotation_matrix,panel_ind,a(1:3))
    a(4) = 0._8
end function calculate_ecs_a2_2d

pure function calculate_ecs_a3_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, h_top) result(a)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, h_top
    real(kind=8)                    :: a(4)

    a = [0.0_8, 0.0_8, 0.0_8, 1.0_8]

end function calculate_ecs_a3_orog

pure function calculate_ecs_a3_2d(this, panel_ind, alpha, beta) result(a)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: a(4)

    a = [0.0_8, 0.0_8, 0.0_8, 1.0_8]

end function calculate_ecs_a3_2d

pure function calculate_ecs_b1_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(b)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: b(3)

    b = calculate_ecs_b1_2d(this, panel_ind, alpha, beta)

end function calculate_ecs_b1_orog

pure function calculate_ecs_b1_2d(this, panel_ind, alpha, beta) result(b)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: b(3)

    b(1:3) = ecs_b1_proto(alpha,beta)
    b(1:3) = ecs_proto2realface(this%topology,this%rotation_matrix,panel_ind,b(1:3))
end function calculate_ecs_b1_2d

pure function calculate_ecs_b2_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(b)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: b(3)

    b = calculate_ecs_b2_2d(this, panel_ind, alpha, beta)

end function calculate_ecs_b2_orog

pure function calculate_ecs_b2_2d(this, panel_ind, alpha, beta) result(b)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: b(3)

    b(1:3) = ecs_b2_proto(alpha,beta)
    b(1:3) = ecs_proto2realface(this%topology,this%rotation_matrix,panel_ind,b(1:3))
end function calculate_ecs_b2_2d

pure function calculate_ecs_b3_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(b)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: b(4)

    b = [0.0_8, 0.0_8, 0.0_8, 1.0_8]

end function calculate_ecs_b3_orog

pure function calculate_ecs_b3_2d(this, panel_ind, alpha, beta) result(b)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: b(4)

    b = [0.0_8, 0.0_8, 0.0_8, 1.0_8]

end function calculate_ecs_b3_2d

pure function calculate_ecs_Q_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(Q)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: Q(6)

    Q = calculate_ecs_Q_2d(this, panel_ind, alpha, beta)

end function calculate_ecs_Q_orog

pure function calculate_ecs_Q_2d(this, panel_ind, alpha, beta) result(Q)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: Q(6)

    Q(1:3) = ecs_Q(this,panel_ind,alpha,beta)
    Q(4:5) = 0.0_8
    Q(6)   = 1.0_8
end function calculate_ecs_Q_2d

pure function calculate_ecs_Qi_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(Qi)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: Qi(6)

    Qi = calculate_ecs_Qi_2d(this, panel_ind, alpha, beta)

end function calculate_ecs_Qi_orog

pure function calculate_ecs_Qi_2d(this, panel_ind, alpha, beta) result(Qi)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: Qi(6)

    Qi(1:3) = ecs_Qi(this,panel_ind,alpha,beta)
    Qi(4:5) = 0.0_8
    Qi(6)   = 1.0_8
end function calculate_ecs_Qi_2d

pure function calculate_ecs_J_orog(this, panel_ind, alpha, beta, eta, h_surf, h_top) result(J)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, h_top
    real(kind=8)                    :: J

    J = calculate_ecs_J_2d(this, panel_ind, alpha, beta)
end function calculate_ecs_J_orog

pure function calculate_ecs_J_2d(this, panel_ind, alpha, beta) result(J)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: J

    J = ecs_J(this,panel_ind,alpha,beta)

end function calculate_ecs_J_2d

pure function calculate_ecs_G_orog(this, panel_ind, alpha, beta, eta, &
                                         h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top) result(G)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta, eta
    real(kind=8),        intent(in) :: h_surf, dh_surf_dalpha, dh_surf_dbeta, h_top
    real(kind=8)                    :: G(3,3,3)

    G = calculate_ecs_G_2d(this, panel_ind, alpha, beta)

end function calculate_ecs_G_orog

pure function calculate_ecs_G_2d(this, panel_ind, alpha, beta) result(G)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: G(3,3,3)

    G(1:3,1:3,1:3) = 0.0_8
    G(1:2,1:2,1:2) = ecs_Christoffel(this,panel_ind,alpha,beta)
end function calculate_ecs_G_2d

pure function ecs_proto2realface(topology,rotation_matrix,panel_ind,r) result(r1)
    !transform prototype-face 3d Cartesian coords to real (spherical / elliptical) face
    class(cubed_sphere_topology_t), intent(in) :: topology
    real(kind=8),                   intent(in) :: rotation_matrix(3,3)
    integer(kind=4),                intent(in) :: panel_ind
    real(kind=8),                   intent(in) :: r(3)
    !output
    real(kind=8) :: r1(3)
    !local
    real(kind=8) :: rtemp(3)


    rtemp(1) = topology%ex(1,panel_ind)*r(1)+topology%ey(1,panel_ind)*r(2)- &
                                                    topology%n(1,panel_ind)*r(3)
    rtemp(2) = topology%ex(2,panel_ind)*r(1)+topology%ey(2,panel_ind)*r(2)- &
                                                    topology%n(2,panel_ind)*r(3)
    rtemp(3) = topology%ex(3,panel_ind)*r(1)+topology%ey(3,panel_ind)*r(2)- &
                                                    topology%n(3,panel_ind)*r(3)
    r1(1) = sum(rotation_matrix(1:3,1)*rtemp(1:3))
    r1(2) = sum(rotation_matrix(1:3,2)*rtemp(1:3))
    r1(3) = sum(rotation_matrix(1:3,3)*rtemp(1:3))
end function ecs_proto2realface

pure function ecs_point_r_proto(alpha,beta) result(r)
    real(kind=8), intent(in)    :: alpha, beta
    real(kind=8)                :: r(3)
    !locals:
    real(kind=8) d

    !grid point at proto cube face: (assumes x = tan(alpha), y = tan(beta), z = 1)
    d = sqrt(1._8+tan(alpha)**2+tan(beta)**2)
    !coordinates at prototype spherical face: vec(r)/||r||, r = (x, y, z)
    r(1) = tan(alpha)/d; r(2) = tan(beta)/d; r(3) = 1._8/d

end function ecs_point_r_proto

pure function ecs_a1_proto(alpha,beta) result(a1)
    !dr / d alpha vector at prototype face for given alpha&beta
    real(kind=8), intent(in)  :: alpha, beta
    real(kind=8)              :: a1(3)
    !local
    real(kind=8) ta, tb

    ta = tan(alpha);   tb = tan(beta)
    !d (xyz)^T/ d alpha
    a1(1) = (1d0+ta**2d0)*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
    a1(2) = -ta*tb*(1d0+ta**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
    a1(3) = -ta*(1d0+ta**2d0) / (1d0+ta**2d0+tb**2d0)**(3d0/2d0)

end function ecs_a1_proto

pure function ecs_a2_proto(alpha,beta) result(a2)
    !dr / d beta vector at prototype face for given alpha&beta
    real(kind=8)              :: a2(3)
    real(kind=8), intent(in)  :: alpha, beta
    !local
    real(kind=8) ta, tb

    ta = tan(alpha);   tb = tan(beta)
    !d (xyz)^T/ d beta
    a2(1) = -ta*tb*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
    a2(2) = (1d0+ta**2d0)*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
    a2(3) = -tb*(1d0+tb**2d0) / (1d0+ta**2d0+tb**2d0)**(3d0/2d0)

end function ecs_a2_proto

pure function ecs_b1_proto(alpha,beta) result(b1)
    !contravariant alpha vector at prototype face for given alpha&beta
    real(kind=8), intent(in)  :: alpha, beta
    real(kind=8)              :: b1(3)
    !local
    real(kind=8) ta, tb, sigm

    ta = tan(alpha);   tb = tan(beta)
    sigm = sqrt(1._8+ta**2+tb**2)
    !d (xyz)^T/ d alpha
    b1(1) = (1._8+tb**2-(ta*tb)**2/(1._8+ta**2)) / sigm
    b1(2) = 0._8
    b1(3) =-ta*(1._8+tb**2/(1+ta**2))/sigm

end function ecs_b1_proto

pure function ecs_b2_proto(alpha,beta) result(b2)
    !contravariant beta vector at prototype face for given alpha&beta
    real(kind=8), intent(in)  :: alpha, beta
    real(kind=8)              :: b2(3)
    !local
    real(kind=8) ta, tb, sigm

    ta = tan(alpha);   tb = tan(beta)
    sigm = sqrt(1._8+ta**2+tb**2)
    !d (xyz)^T/ d beta
    b2(1) = 0._8
    b2(2) = ((1+ta**2)-(ta*tb)**2/(1+tb**2))/sigm
    b2(3) =-tb*(ta**2/(1+tb**2)+1._8)/sigm

end function ecs_b2_proto

pure function ecs_Q(this,panel_ind,alpha,beta) result(Q)
    !Calculate metric tensor
    !Q = |a1*a1  a1*a2| = |Q(1) Q(2)|
    !    |a1*a2  a2*a2|   |Q(2) Q(3)|
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: Q(3)
    !local
    real(kind=8) ta, tb, sigm

    ta = tan(alpha);   tb = tan(beta)
    sigm = sqrt(1._8+ta**2+tb**2)

    Q(1) = (1._8+ta**2)**2*(1._8+tb**2)/sigm**4
    Q(2) =-(1._8+ta**2)*(1._8+tb**2)*ta*tb/sigm**4
    Q(3) = (1._8+ta**2)*(1._8+tb**2)**2/sigm**4

end function ecs_Q

pure function ecs_QI(this,panel_ind, alpha, beta) result(QI)
    !Calculate inverse metric tensor
    !QI = inv |a1*a1  a1*a2| = |QI(1) QI(2)|
    !         |a1*a2  a2*a2|   |QI(2) QI(3)|
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha,beta
    real(kind=8)                    :: QI(3)
    !local
    real(kind=8) ta, tb, sigm

    ta = tan(alpha);   tb = tan(beta)
    sigm = sqrt(1._8+ta**2+tb**2)

    Qi(1) = sigm**2 / (1._8+ta**2)
    Qi(2) = sigm**2*ta*tb / ((1+ta**2)*(1+tb**2))
    Qi(3) = sigm**2 / (1._8+tb**2)

end function ecs_QI

pure function ecs_Christoffel(this,panel_ind, alpha, beta) result(G)
    !Calculate christoffel symbols
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha,beta
    real(kind=8)                    :: G(2,2,2)
    !local
    real(kind=8) X, Y, d2

    X = tan(alpha);   Y = tan(beta)
    d2 = 1._8+X**2+Y**2

    G(1,1,1) = 2.0_8*X*Y**2 / d2
    G(1,2,1) =-Y*(1+Y**2) / d2
    G(2,1,1) =-Y*(1+Y**2) / d2
    G(2,2,1) = 0.0_8

    G(1,1,2) = 0.0_8
    G(1,2,2) =-X*(1+X**2) / d2
    G(2,1,2) =-X*(1+X**2) / d2
    G(2,2,2) = 2.0_8*X**2*Y / d2

end function ecs_Christoffel

pure function ecs_J(this,panel_ind,alpha,beta) result(J)
    !Compute sqrt of metric tensor
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: J
    !local
    real(kind=8) ta,tb,sigm

    ta = tan(alpha);   tb = tan(beta)
    sigm = sqrt(1._8+ta**2+tb**2)

    J = (1._8+ta**2)*(1._8+tb**2) / sigm**3
end function ecs_J

subroutine transform_cartesian_to_native_ecs(this,panel_ind, alpha, beta, r)
    import metric_t
    class(ecs_metric_t), intent(in)  :: this
    integer(kind=4),     intent(out) :: panel_ind
    real(kind=8),        intent(out) :: alpha, beta
    real(kind=8),        intent(in)  :: r(3)

    real(kind=8)    :: r_dot_n, r_dot_n_max, r_loc(3)
    integer(kind=8) :: ipanel

    !find panel with maximum projection of r onto cube outer normal
    !this will be the panel the point r belongs to
    panel_ind = 1
    r_dot_n_max = 0.0_8
    do ipanel=1,6
        !minus sign because of inner normal stored in topology
        r_dot_n =-sum(r(1:3)*this%topology%n(1:3,ipanel))
        if(r_dot_n > r_dot_n_max) then
            panel_ind = ipanel
            r_dot_n_max = r_dot_n
        end if
    end do

    !it is assumed that ||r|| == 1
    alpha = atan( sum(this%topology%ex(1:3,panel_ind)*r(1:3)) / r_dot_n_max)
    beta  = atan( sum(this%topology%ey(1:3,panel_ind)*r(1:3)) / r_dot_n_max)
end subroutine transform_cartesian_to_native_ecs

end module ecs_metric_mod

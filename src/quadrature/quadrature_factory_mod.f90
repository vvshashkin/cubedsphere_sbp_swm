module quadrature_factory_mod

use abstract_quadrature_mod, only : quadrature_t
use default_quadrature_mod,  only : default_tile_quadrature_t
use parcomm_mod,             only : parcomm_global
use mesh_mod,                only : mesh_t

implicit none

contains

subroutine create_quadrature(quadrature, quadrature_name, mesh)

    class(quadrature_t), allocatable, intent(out) :: quadrature
    character(len=*),                 intent(in)  :: quadrature_name
    type(mesh_t),                     intent(in)  :: mesh

    integer(kind=4) :: ts, te

    allocate(quadrature)
    if(quadrature_name == "default_quadrature") then
        ts = mesh%ts
        te = mesh%te
        allocate(default_tile_quadrature_t :: quadrature%tile(ts:te))
    else if(quadrature_name == "SBP_Ah21_quadrature" .or. &
            quadrature_name == "SBP_Ah42_quadrature" .or. &
            quadrature_name == "SBP_Ah63_quadrature") then
        call create_ah_sbp_quadrature(quadrature, quadrature_name, mesh)
    else if(quadrature_name == "SBP_C21_quadrature" .or. &
            quadrature_name == "SBP_C42_quadrature") then
        call create_C_sbp_quadrature(quadrature, quadrature_name, mesh)
    else
        call parcomm_global%abort("unknown quadrature:"// quadrature_name //&
                                  " in create_quadrature, quadrature_factory_mod")
    end if

end subroutine create_quadrature

subroutine create_ah_sbp_quadrature(quadrature, quadrature_name, mesh)

    use sbp_quadrature_mod, only : sbp_tile_quadrature_t
    use sbp_operators_collection_mod, only : Q21_A, Q42_A, Q63_A

    class(quadrature_t), allocatable, intent(inout) :: quadrature
    character(len=*),                 intent(in)    :: quadrature_name
    type(mesh_t),                     intent(in)    :: mesh

    type(sbp_tile_quadrature_t), allocatable :: tile_q(:)
    integer(kind=4) :: t, is, ie, js, je, i, j, n
    real(kind=8) A_edge(8)
    integer(kind=4) :: nedge

    if(mesh%tile(mesh%ts)%points_type /= "xy") &
        call parcomm_global%abort("create_ah_sbp_quadrature is only for xy-point type meshes")

    if(quadrature_name == "SBP_Ah21_quadrature") then
        nedge = size(Q21_A)
        A_edge(1:nedge) =  Q21_A(1:nedge)
    else if(quadrature_name == "SBP_Ah42_quadrature") then
        nedge = size(Q42_A)
        A_edge(1:nedge) =  Q42_A(1:nedge)
    else if(quadrature_name == "SBP_Ah63_quadrature") then
        nedge = size(Q63_A)
        A_edge(1:nedge) =  Q63_A(1:nedge)
    else
        call parcomm_global%abort("unknown quadrature:"// quadrature_name //&
                                  " in create_Ah_sbp_quadrature, quadrature_factory_mod")
    end if

    allocate(tile_q(mesh%ts:mesh%te))

    do t = mesh%ts, mesh%te
       is = mesh%tile(t)%is; ie = mesh%tile(t)%ie
       js = mesh%tile(t)%js; je = mesh%tile(t)%je

       allocate(tile_q(t)%Ax(is:ie))
       allocate(tile_q(t)%Ay(js:je))

       call create_mass_coefficients(tile_q(t)%Ax, is, ie, mesh%tile(t)%nx, &
                                                                 A_edge(1:nedge), nedge)
       call create_mass_coefficients(tile_q(t)%Ay, js, je, mesh%tile(t)%ny, &
                                                                 A_edge(1:nedge), nedge)
   end do

   call move_alloc(tile_q, quadrature%tile)
end

subroutine create_C_sbp_quadrature(quadrature, quadrature_name, mesh)

    use sbp_quadrature_mod,           only : sbp_tile_quadrature_t
    use sbp_operators_collection_mod, only : D42_A_interfaces, D42_A_centers, &
                                             D21_A_interfaces, D21_A_centers

    class(quadrature_t), allocatable, intent(inout) :: quadrature
    character(len=*),                 intent(in)    :: quadrature_name
    type(mesh_t),                     intent(in)    :: mesh

    type(sbp_tile_quadrature_t), allocatable :: tile_q(:)
    integer(kind=4) :: t, is, ie, js, je, i, j, n
    real(kind=8) A_edge_c(8)
    real(kind=8) A_edge_i(8)
    integer(kind=4) :: nedge_c, nedge_i

    if(quadrature_name == "SBP_C21_quadrature") then
        nedge_i = size(D21_A_interfaces)
        A_edge_i(1:nedge_i) =  D21_A_interfaces(1:nedge_i)
        nedge_c = size(D21_A_centers)
        A_edge_c(1:nedge_c) =  D21_A_centers(1:nedge_c)
    else if(quadrature_name == "SBP_C42_quadrature") then
        nedge_i = size(D42_A_interfaces)
        A_edge_i(1:nedge_i) =  D42_A_interfaces(1:nedge_i)
        nedge_c = size(D42_A_centers)
        A_edge_c(1:nedge_c) =  D42_A_centers(1:nedge_c)
    else
        call parcomm_global%abort("unknown quadrature:"// quadrature_name //&
                                  " in create_C_sbp_quadrature, quadrature_factory_mod")
    end if

    allocate(tile_q(mesh%ts:mesh%te))

    do t = mesh%ts, mesh%te
        is = mesh%tile(t)%is; ie = mesh%tile(t)%ie
        js = mesh%tile(t)%js; je = mesh%tile(t)%je

        allocate(tile_q(t)%Ax(is:ie))
        allocate(tile_q(t)%Ay(js:je))

        if(mesh%tile(t)%points_type == "o" .or. mesh%tile(t)%points_type == "c") then
            call create_mass_coefficients(tile_q(t)%Ax, is, ie, mesh%tile(t)%nx, &
                                                        A_edge_c(1:nedge_c), nedge_c)
            call create_mass_coefficients(tile_q(t)%Ay, js, je, mesh%tile(t)%ny, &
                                                        A_edge_c(1:nedge_c), nedge_c)
        else if(mesh%tile(t)%points_type == "x") then
            call create_mass_coefficients(tile_q(t)%Ax, is, ie, mesh%tile(t)%nx, &
                                                        A_edge_i(1:nedge_i), nedge_i)
            call create_mass_coefficients(tile_q(t)%Ay, js, je, mesh%tile(t)%ny, &
                                                        A_edge_c(1:nedge_c), nedge_c)
        else if(mesh%tile(t)%points_type == "y") then
            call create_mass_coefficients(tile_q(t)%Ax, is, ie, mesh%tile(t)%nx, &
                                                        A_edge_c(1:nedge_c), nedge_c)
            call create_mass_coefficients(tile_q(t)%Ay, js, je, mesh%tile(t)%ny, &
                                                        A_edge_i(1:nedge_i), nedge_i)
        else if(mesh%tile(t)%points_type == "xy") then
            call create_mass_coefficients(tile_q(t)%Ax, is, ie, mesh%tile(t)%nx, &
                                                        A_edge_i(1:nedge_i), nedge_i)
            call create_mass_coefficients(tile_q(t)%Ay, js, je, mesh%tile(t)%ny, &
                                                        A_edge_i(1:nedge_i), nedge_i)
        else
            call parcomm_global%abort("unknown mesh points type in create_C_sbp_quadrature, "// &
                                      "quadrature_factory_mod: "// mesh%tile(t)%points_type)
        end if
   end do

   call move_alloc(tile_q, quadrature%tile)
end


subroutine create_mass_coefficients(A, is, ie, n, A_edge, n_edge)

    integer(kind=4), intent(in)    :: is, ie, n, n_edge
    real(kind=8),    intent(in)    :: A_edge(n_edge)

    real(kind=8),    intent(inout) :: A(is:ie)

    integer(kind=8) :: i

    A(is:ie) = 1.0_8

    do i = max(1, is), min(n_edge,ie)
        A(i) = A_edge(i)
    end do

    do i = max(n-n_edge+1, is), min(n,ie)
        A(i) = A_edge(n-i+1)
    end do

end subroutine create_mass_coefficients

end module quadrature_factory_mod

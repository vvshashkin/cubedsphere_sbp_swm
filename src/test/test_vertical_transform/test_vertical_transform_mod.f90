module test_vertical_transform_mod

implicit none

contains

subroutine test_vertical_transform()
    use abstract_vertical_transform_mod, only : vertical_transform_t
    use vertical_transform_factory_mod,  only : create_vertical_transform

    class(vertical_transform_t), allocatable :: vert_transform
    real(kind=8),    parameter :: h_surf = 1e3_8, h_top=30e3_8
    integer(kind=4), parameter :: nz = 100
    real(kind=8),    parameter :: deta = 1.0_8 / nz
    real(kind=8),    parameter :: tol1 = 1e-16_8, tol2 = 1e-5
    real(kind=8) :: dz, e1, e2, e3
    logical :: is_passed
    integer(kind=4) :: k

    vert_transform = create_vertical_transform("vertical_transform_default")

    e1 = abs(vert_transform%calc_z(h_surf, h_top, eta=0.0_8) - h_surf)+ &
         abs(vert_transform%calc_dz_dh_surf(eta=0.0_8) - 1.0_8) +         &
         abs(vert_transform%calc_dz_dh_surf(eta=1.0_8) - 0.0_8)
    e2 = abs(vert_transform%calc_z(h_surf, h_top, eta=1.0_8) - h_top)+ &
         abs(vert_transform%calc_dz_dh_top(eta=0.0_8) - 0.0_8) +         &
         abs(vert_transform%calc_dz_dh_top(eta=1.0_8) - 1.0_8)

    e3 = 0.0_8
    do k=1, nz
        dz = vert_transform%calc_z(h_surf, h_top, eta=k*deta)- &
             vert_transform%calc_z(h_surf, h_top, eta=(k-1)*deta)
        e3 = max(e3, abs(dz - vert_transform%calc_dz_deta(h_surf, h_top, eta=(k-0.5_8)*deta)*deta))
    end do

    print *, "vertical transform test: err_hsurf = ", e1, " err h_top = ", e2, &
             "err dz = ", e3
    if(e1 < tol1 .and. e2 < tol1 .and. e3 < tol2) then
        print *, "vertical transform test passed"
    else
        print *, "vertical_transform test failed"
    end if
end subroutine test_vertical_transform

end module test_vertical_transform_mod

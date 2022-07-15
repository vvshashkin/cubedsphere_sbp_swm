module tile_mod
implicit none

type, public :: tile_t
    integer(kind=4) :: is, ie, js, je, ks, ke
contains
    procedure, public :: init
    procedure, public :: get_points_num
    procedure, public :: check
    procedure, public :: print
    procedure, public :: getind
end type tile_t

private

contains

subroutine init(this, is, ie, js, je, ks, ke)
    class(tile_t), intent(inout) :: this
    integer(kind=4), intent(in) :: is, ie, js, je, ks, ke

    this%is = is; this%ie = ie
    this%js = js; this%je = je
    this%ks = ks; this%ke = ke

    ! this%panel_number = panel_number

    call this%check()

end subroutine init

subroutine check(this)
    class(tile_t), intent(in) :: this
    logical :: is_test_passed

    is_test_passed = .true.
    if (this%ie<this%is) is_test_passed = .false.
    if (this%je<this%js) is_test_passed = .false.
    if (this%ke<this%ks) is_test_passed = .false.

    if (not(is_test_passed)) then
        print*, 'Error in tile_mod!!!'
    end if

end subroutine check

pure function get_points_num(this) result(points_num)

    class(tile_t), intent(in) :: this
    integer(kind=4) :: points_num

    points_num = (this%ke - this%ks + 1)* &
                 (this%je - this%js + 1)* &
                 (this%ie - this%is + 1)

end function get_points_num

subroutine print(this)
    class(tile_t), intent(in) :: this
    character(len=1000) :: is, ie, js, je
    write(is,*) this%is; write(ie,*) this%ie
    write(js,*) this%js; write(je,*) this%je
    print*, ''
    print '(2(1x,a, a), /, 2(1x,a, a))', 'is = ', trim(adjustl(is)), 'ie = ', trim(adjustl(ie)), &
                                         'js = ', trim(adjustl(js)), 'je = ', trim(adjustl(je))

end subroutine print

subroutine getind(this, is,ie,js,je,ks,ke)
    class(tile_t), intent(in) :: this
    integer(kind=4), intent(out), optional :: is, ie, js, je, ks, ke

    if(present(is)) is = this%is
    if(present(ie)) ie = this%ie
    if(present(js)) js = this%js
    if(present(je)) je = this%je
    if(present(ks)) ks = this%ks
    if(present(ke)) ke = this%ke

end subroutine getind


end module

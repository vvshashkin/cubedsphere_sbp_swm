module grid_field_mod

use mesh_mod, only : mesh_t, tile_mesh_t

implicit none

type, public :: grid_field_t
    integer(kind=4)                 :: ts, te ! tile_field_t-array bounds
    type(tile_field_t), allocatable :: tile(:)
contains
    procedure, public :: update_s1 => update_grid_field_s1 !v = v + s1*unity
    procedure, public :: update_s1v1 => update_grid_field_s1v1 !v = v + s1*v1
    procedure, public :: update_s1v1s2v2 => update_grid_field_s1v1s2v2 !v = v + s1*v1+s2*v2
    generic :: update => update_s1, update_s1v1, update_s1v1s2v2

    procedure, public :: assign_s1v1       => assign_grid_field_s1v1 !v = s1*v1
    procedure, public :: assign_s1         => assign_grid_field_s1 !v = s1
    procedure, public :: assign_v1         => assign_grid_field_v1 !v = v1
    procedure, public :: assign_s1v1s2v2   => assign_grid_field_s1v1s2v2 !v = s1*v1+s2*v2
    generic :: assign => assign_s1v1, assign_s1, assign_v1, assign_s1v1s2v2

    procedure, public :: assign_prod  => assign_grid_field_prod_s1v1v2
    procedure, public :: assign_ratio => assign_grid_field_ratio_s1v1v2

    procedure, public :: copy => copy_grid_field
    procedure, public :: create_similar => create_similar_grid_field
    procedure, public :: algebraic_norm2 => compute_grid_field_algebraic_norm2
    procedure, public :: algebraic_dot   => compute_grid_field_algebraic_dot
    procedure, public :: maximum => compute_grid_field_maximum
    procedure, public :: maxabs  => compute_grid_field_maxabs
    procedure, public :: minimum => compute_grid_field_minimum
    ! procedure, public :: init => grid_field_init
end type grid_field_t

type, public :: tile_field_t
    real(kind=8), allocatable  :: p(:,:,:)
!    integer(kind=4)            :: panel_ind ! index of grid panel hosted the tile_field
    integer(kind=4)            :: is, ie, js, je, ks, ke ! p-array bounds
contains
    procedure, public :: init => tile_field_init

    procedure, public :: update_s1 => tile_field_update_s1!v = v + s1*unity
    procedure, public :: update_s1v1 => tile_field_update_s1v1!v = v + s1*v1
    procedure, public :: update_s1v1s2v2 => tile_field_update_s1v1s2v2!v = v + s1*v1+s2*v2
    !update -- generic procedure for v = v + ... operations
    generic :: update => update_s1, update_s1v1, update_s1v1s2v2

    procedure, public :: assign_s1       => tile_field_assign_s1  !v = s1
    procedure, public :: assign_v1       => tile_field_assign_v1  !v = v1
    procedure, public :: assign_s1v1     => tile_field_assign_s1v1!v = s1*v1
    procedure, public :: assign_s1v1s2v2 => tile_field_assign_s1v1s2v2!v = s1*v1+s2*v2
    !assign -- generic procedure for v = ... operations
    generic :: assign => assign_s1v1, assign_s1, assign_v1, assign_s1v1s2v2

    procedure, public :: assign_prod  => tile_field_assign_prod_s1v1v2
    procedure, public :: assign_ratio => tile_field_assign_ratio_s1v1v2

    procedure, public :: algebraic_norm2 => compute_tile_field_algebraic_norm2
    procedure, public :: algebraic_dot => compute_tile_field_algebraic_dot
    procedure, public :: maximum => compute_tile_field_maximum
    procedure, public :: maxabs  => compute_tile_field_maxabs
    procedure, public :: minimum => compute_tile_field_minimum
end type tile_field_t

contains

subroutine tile_field_init(this, is, ie, js, je, ks, ke)

    class(tile_field_t),  intent(out) :: this
    integer(kind=4),      intent(in)  :: is, ie, js, je, ks, ke

    if (is>ie .or. js>je .or. ks>ke) then
        print*, 'Error! Problem with grid function initialization! Abort!'
        stop
    end if

    allocate( this%p(is : ie, js : je, ks : ke) )

    this%is = is; this%ie = ie
    this%js = js; this%je = je
    this%ks = ks; this%ke = ke

end subroutine tile_field_init

function create_similar_grid_field(this) result(grid_field)

    class(grid_field_t), intent(in)  :: this
    !type(mesh_t),        intent(in)  :: mesh
    type(grid_field_t)               :: grid_field


    integer(kind=4) :: t

    allocate(grid_field%tile(this%ts:this%te))

    grid_field%ts = this%ts
    grid_field%te = this%te

    do t = this%ts, this%te

        call grid_field%tile(t)%init(this%tile(t)%is, this%tile(t)%ie, &
                                     this%tile(t)%js, this%tile(t)%je, &
                                     this%tile(t)%ks, this%tile(t)%ke)
    end do

end function create_similar_grid_field

function copy_grid_field(this) result(grid_field)

    class(grid_field_t), intent(in)  :: this
    !type(mesh_t),        intent(in)  :: mesh
    type(grid_field_t)               :: grid_field

    integer(kind=4) :: t, i, j, k

    allocate(grid_field%tile(this%ts:this%te))

    do t = this%ts, this%te

        call grid_field%tile(t)%init(this%tile(t)%is, this%tile(t)%ie, &
                                     this%tile(t)%js, this%tile(t)%je, &
                                     this%tile(t)%ks, this%tile(t)%ke)

        do k = this%tile(t)%ks,this%tile(t)%ke
            do j = this%tile(t)%js,this%tile(t)%je
                do i = this%tile(t)%is,this%tile(t)%ie
                    grid_field%tile(t)%p(i,j,k) = this%tile(t)%p(i,j,k)
                end do
            end do
        end do

    end do

end function copy_grid_field

function compute_grid_field_algebraic_norm2(this, mesh, parcomm) result(norm2)
        use parcomm_mod, only : parcomm_t
        use mpi

        class(grid_field_t), intent(in)  :: this
        type(mesh_t),        intent(in)  :: mesh
        type(parcomm_t),     intent(in)  :: parcomm
        real(kind=8)                     :: norm2

        real(kind=8) :: local_norm2
        integer(kind=4) :: t
        integer(kind=4) :: ierr

        local_norm2 = 0.0

        do t = mesh%ts, mesh%te
            local_norm2 = local_norm2 + this%tile(t)%algebraic_norm2(mesh%tile(t))
        end do

        call mpi_allreduce(local_norm2, norm2, 1, mpi_double, mpi_sum, parcomm%comm_w, ierr)
        norm2 = sqrt(norm2)
end function compute_grid_field_algebraic_norm2

function compute_tile_field_algebraic_norm2(this, mesh) result(norm2)
        use parcomm_mod, only : parcomm_t

        class(tile_field_t), intent(in)  :: this
        type(tile_mesh_t),   intent(in)  :: mesh
        real(kind=8)                     :: norm2

        integer(kind=4) :: i, j, k

        norm2 = 0.0_8

        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    norm2 = norm2 + this%p(i,j,k)**2
                end do
            end do
        end do
end function compute_tile_field_algebraic_norm2

function compute_grid_field_algebraic_dot(this, other, mesh, parcomm) result(dot_product)
        use parcomm_mod, only : parcomm_t
        use mpi

        class(grid_field_t), intent(in)  :: this
        type(grid_field_t),  intent(in)  :: other
        type(mesh_t),        intent(in)  :: mesh
        type(parcomm_t),     intent(in)  :: parcomm
        real(kind=8)                     :: dot_product

        real(kind=8) :: local_dot
        integer(kind=4) :: t
        integer(kind=4) :: ierr

        local_dot = 0.0

        do t = mesh%ts, mesh%te
            local_dot = local_dot + this%tile(t)%algebraic_dot(other%tile(t),mesh%tile(t))
        end do

        call mpi_allreduce(local_dot, dot_product, 1, mpi_double, mpi_sum, parcomm%comm_w, ierr)
end function compute_grid_field_algebraic_dot

function compute_tile_field_algebraic_dot(this, other, mesh) result(dot_product)
        use parcomm_mod, only : parcomm_t

        class(tile_field_t), intent(in)  :: this
        type(tile_field_t),  intent(in)  :: other
        type(tile_mesh_t),   intent(in)  :: mesh
        real(kind=8)                     :: dot_product

        integer(kind=4) :: i, j, k

        dot_product = 0.0_8

        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    dot_product = dot_product + this%p(i,j,k)*other%p(i,j,k)
                end do
            end do
        end do
end function compute_tile_field_algebraic_dot

subroutine update_grid_field_s1(this, scalar1, mesh)

    class(grid_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%update(scalar1, mesh%tile(t))
    end do

end subroutine update_grid_field_s1

subroutine update_grid_field_s1v1(this, scalar1, f1, mesh)

    class(grid_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    type(grid_field_t),  intent(in)    :: f1
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%update(scalar1, f1%tile(t), mesh%tile(t))
    end do

end subroutine update_grid_field_s1v1

subroutine update_grid_field_s1v1s2v2(this, scalar1, f1, scalar2, f2, mesh)

    class(grid_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1, scalar2
    type(grid_field_t),  intent(in)    :: f1, f2
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%update(scalar1, f1%tile(t), scalar2, f2%tile(t), mesh%tile(t))
    end do

end subroutine update_grid_field_s1v1s2v2

subroutine assign_grid_field_s1(this, scalar1, mesh)

    class(grid_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%assign(scalar1, mesh%tile(t))
    end do

end subroutine assign_grid_field_s1

subroutine assign_grid_field_v1(this, v1, mesh)

    class(grid_field_t), intent(inout) :: this
    class(grid_field_t), intent(in)    :: v1
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%assign(v1%tile(t), mesh%tile(t))
    end do

end subroutine assign_grid_field_v1

subroutine assign_grid_field_s1v1(this, scalar1, f1, mesh)

    class(grid_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    type(grid_field_t),  intent(in)    :: f1
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%assign(scalar1, f1%tile(t), mesh%tile(t))
    end do

end subroutine assign_grid_field_s1v1

subroutine assign_grid_field_s1v1s2v2(this, scalar1, f1, scalar2, f2, mesh)

    class(grid_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1, scalar2
    type(grid_field_t),  intent(in)    :: f1, f2
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%assign(scalar1, f1%tile(t), scalar2, f2%tile(t), mesh%tile(t))
    end do

end subroutine assign_grid_field_s1v1s2v2

subroutine assign_grid_field_prod_s1v1v2(this, scalar1, f1, f2, mesh)

    class(grid_field_t), intent(inout) :: this
    type(grid_field_t),  intent(in)    :: f1, f2
    real(kind=8),        intent(in)    :: scalar1
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%assign_prod(scalar1, f1%tile(t), f2%tile(t), mesh%tile(t))
    end do

end subroutine assign_grid_field_prod_s1v1v2

subroutine assign_grid_field_ratio_s1v1v2(this, scalar1, f1, f2, mesh)

    class(grid_field_t), intent(inout) :: this
    type(grid_field_t),  intent(in)    :: f1, f2
    real(kind=8),        intent(in)    :: scalar1
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%assign_ratio(scalar1, f1%tile(t), f2%tile(t), mesh%tile(t))
    end do

end subroutine assign_grid_field_ratio_s1v1v2

subroutine tile_field_assign_prod_s1v1v2(this, scalar1, v1, v2, mesh)

    class(tile_field_t), intent(inout) :: this
    type(tile_field_t),  intent(in)    :: v1, v2
    real(kind=8),        intent(in)    :: scalar1
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = scalar1*v1%p(i,j,k)*v2%p(i,j,k)
            end do
        end do
    end do

end subroutine tile_field_assign_prod_s1v1v2

subroutine tile_field_assign_ratio_s1v1v2(this, scalar1, v1, v2, mesh)

    class(tile_field_t), intent(inout) :: this
    type(tile_field_t),  intent(in)    :: v1, v2
    real(kind=8),        intent(in)    :: scalar1
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = scalar1*v1%p(i,j,k)/v2%p(i,j,k)
            end do
        end do
    end do

end subroutine tile_field_assign_ratio_s1v1v2

subroutine tile_field_update_s1(this, scalar1, mesh)

    class(tile_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = this%p(i,j,k)+scalar1
            end do
        end do
    end do

end subroutine tile_field_update_s1

subroutine tile_field_update_s1v1(this, scalar1, v1, mesh)

    class(tile_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    type(tile_field_t),  intent(in)    :: v1
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = this%p(i,j,k) + scalar1*v1%p(i,j,k)
            end do
        end do
    end do

end subroutine tile_field_update_s1v1

subroutine tile_field_update_s1v1s2v2(this, scalar1, v1, scalar2, v2, mesh)

    class(tile_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1, scalar2
    type(tile_field_t),  intent(in)    :: v1, v2
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = this%p(i,j,k)+scalar1*v1%p(i,j,k)+scalar2*v2%p(i,j,k)
            end do
        end do
    end do

end subroutine tile_field_update_s1v1s2v2

subroutine tile_field_assign_s1v1(this, scalar1, v1, mesh)

    class(tile_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    type(tile_field_t),  intent(in)    :: v1
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = scalar1*v1%p(i,j,k)
            end do
        end do
    end do

end subroutine tile_field_assign_s1v1

subroutine tile_field_assign_s1v1s2v2(this, scalar1, v1, scalar2, v2, mesh)

    class(tile_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1, scalar2
    type(tile_field_t),  intent(in)    :: v1, v2
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = scalar1*v1%p(i,j,k)+scalar2*v2%p(i,j,k)
            end do
        end do
    end do

end subroutine tile_field_assign_s1v1s2v2

subroutine tile_field_assign_s1(this, scalar1, mesh)

    class(tile_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = scalar1
            end do
        end do
    end do

end subroutine tile_field_assign_s1

subroutine tile_field_assign_v1(this, v1, mesh)

    class(tile_field_t), intent(inout) :: this
    class(tile_field_t),  intent(in)   :: v1
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = v1%p(i,j,k)
            end do
        end do
    end do

end subroutine tile_field_assign_v1

function compute_grid_field_maximum(this, mesh, parcomm) result(maximum_value)
    use parcomm_mod, only : parcomm_t
    use mpi

    class(grid_field_t), intent(in)  :: this
    type(mesh_t),        intent(in)  :: mesh
    type(parcomm_t),     intent(in)  :: parcomm
    real(kind=8)                     :: maximum_value

    real(kind=8) :: local_max
    integer(kind=4) :: t
    integer(kind=4) :: ierr

    local_max = this%tile(mesh%ts)%maximum(mesh%tile(mesh%ts))
    do t = mesh%ts+1, mesh%te
        local_max = max(local_max,this%tile(t)%maximum(mesh%tile(t)))
    end do

    call mpi_allreduce(local_max, maximum_value, 1, mpi_double, mpi_max, parcomm%comm_w, ierr)
end function compute_grid_field_maximum

function compute_tile_field_maximum(this, mesh) result(maximum_value)
    use parcomm_mod, only : parcomm_t

    class(tile_field_t), intent(in)  :: this
    type(tile_mesh_t),   intent(in)  :: mesh
    real(kind=8)                     :: maximum_value

    integer(kind=4) :: is,ie,js,je,ks,ke
    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    maximum_value = maxval(this%p(is:ie,js:je,ks:ke))

end function compute_tile_field_maximum

function compute_grid_field_minimum(this, mesh, parcomm) result(minimum_value)
    use parcomm_mod, only : parcomm_t
    use mpi

    class(grid_field_t), intent(in)  :: this
    type(mesh_t),        intent(in)  :: mesh
    type(parcomm_t),     intent(in)  :: parcomm
    real(kind=8)                     :: minimum_value

    real(kind=8) :: local_min
    integer(kind=4) :: t
    integer(kind=4) :: ierr

    local_min = this%tile(mesh%ts)%minimum(mesh%tile(mesh%ts))
    do t = mesh%ts+1, mesh%te
        local_min = min(local_min,this%tile(t)%minimum(mesh%tile(t)))
    end do

    call mpi_allreduce(local_min, minimum_value, 1, mpi_double, mpi_min, parcomm%comm_w, ierr)
end function compute_grid_field_minimum

function compute_tile_field_minimum(this, mesh) result(minimum_value)
    use parcomm_mod, only : parcomm_t

    class(tile_field_t), intent(in)  :: this
    type(tile_mesh_t),   intent(in)  :: mesh
    real(kind=8)                     :: minimum_value

    integer(kind=4) :: is,ie,js,je,ks,ke
    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    minimum_value = minval(this%p(is:ie,js:je,ks:ke))

end function compute_tile_field_minimum

function compute_grid_field_maxabs(this, mesh, parcomm) result(maxabs_value)
    use parcomm_mod, only : parcomm_t
    use mpi

    class(grid_field_t), intent(in)  :: this
    type(mesh_t),        intent(in)  :: mesh
    type(parcomm_t),     intent(in)  :: parcomm
    real(kind=8)                     :: maxabs_value

    real(kind=8) :: local_maxabs
    integer(kind=4) :: t
    integer(kind=4) :: ierr

    local_maxabs = 0.0_8
    do t = mesh%ts, mesh%te
        local_maxabs = max(local_maxabs,this%tile(t)%maxabs(mesh%tile(t)))
    end do

    call mpi_allreduce(local_maxabs, maxabs_value, 1, mpi_double, mpi_max, parcomm%comm_w, ierr)
end function compute_grid_field_maxabs

function compute_tile_field_maxabs(this, mesh) result(maxabs_value)
    use parcomm_mod, only : parcomm_t

    class(tile_field_t), intent(in)  :: this
    type(tile_mesh_t),   intent(in)  :: mesh
    real(kind=8)                     :: maxabs_value

    integer(kind=4) :: is,ie,js,je,ks,ke
    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    maxabs_value = maxval(abs(this%p(is:ie,js:je,ks:ke)))

end function compute_tile_field_maxabs
end module grid_field_mod

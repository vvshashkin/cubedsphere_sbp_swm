module vertical_transform_factory_mod

use abstract_vertical_transform_mod, only : vertical_transform_t
use vertical_transform_mod,          only : vertical_transform_default_t
use parcomm_mod,                     only : parcomm_global

contains

function create_vertical_transform(vert_transform_name) result(vert_transform)

    character(len=*) :: vert_transform_name
    class(vertical_transform_t), allocatable :: vert_transform

    select case(vert_transform_name)
    case("vertical_transform_default")
        vert_transform = vertical_transform_default_t()
    case default
        call parcomm_global%abort("create_vertical_transform, unknown vertical transform:"// &
                                                                        vert_transform_name)
    end select
end function create_vertical_transform

end module vertical_transform_factory_mod

module config_postnh_mod

use config_mod,        only : config_t

implicit none

type, extends(config_t) :: config_postnh_t

    integer(kind=4) :: Nlon, Nlat
    character(:), allocatable :: outputer_name

    contains
    procedure, public :: parse

end type config_postnh_t

contains

subroutine parse(this, config_string)

    class(config_postnh_t),  intent(inout) :: this
    character(len=*),        intent(in)    :: config_string

    integer(kind=4)    :: Nlon=180, Nlat=90
    character(len=256) :: outputer_name

    namelist /nh_postprocessing/ Nlon, Nlat, outputer_name

    read(config_string, nh_postprocessing)

    this%Nlon = Nlon
    this%Nlat = Nlat
    this%outputer_name = trim(outputer_name)

end subroutine parse


end module config_postnh_mod

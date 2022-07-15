module latlon_functions_mod
    implicit none
contains
    real(kind=8) pure function sin_phi(x,y,z) result(s)
        real(kind=8), intent(in) :: x,y,z
        s = z / sqrt(x**2+y**2+z**2)
    end function sin_phi

    real(kind=8) pure function cos_phi(x,y,z) result(c)
        real(kind=8), intent(in) :: x,y,z
        c = sqrt((x**2+y**2)/(x**2+y**2+z**2))
    end function cos_phi

    real(kind=8) pure function cos_lam(x,y,z) result(c)
        real(kind=8), intent(in) :: x,y,z
        c = x / max(sqrt(x**2+y**2),1e-14)
    end function cos_lam

    real(kind=8) pure function sin_lam(x,y,z) result(s)
        real(kind=8), intent(in) :: x,y,z
        !s = cos_lam(x,y,z)
        !s = sqrt(1.0_8-s**2)
        s = y / max(sqrt(x**2+y**2),1e-14)
        if(x**2+y**2 == 0) s = 1.0_8
    end function sin_lam

    subroutine transform_cartesian_hor_vector_to_latlon(u,v,Npts,nlev,vx,vy,vz,lat,lon)
        integer(kind=4), intent(in)  :: Npts, nlev
        real(kind=8),    intent(in)  :: vx(Npts,nlev), vy(Npts,nlev), vz(Npts,nlev)
        real(kind=8),    intent(in)  :: lon(Npts), lat(Npts)
        real(kind=8),    intent(out) :: u(Npts,nlev), v(Npts,nlev)

        integer(kind=4) :: i, k

        do k=1, nlev
            do i=1, Npts
                u(i,k) = -sin(lon(i))*vx(i,k)+cos(lon(i))*vy(i,k)
                v(i,k) = -sin(lat(i))*cos(lon(i))*vx(i,k)- &
                          sin(lat(i))*sin(lon(i))*vy(i,k)+cos(lat(i))*vz(i,k)
            end do
        end do
    end subroutine transform_cartesian_hor_vector_to_latlon

end module latlon_functions_mod

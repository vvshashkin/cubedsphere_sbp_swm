module solid_rotation_fields_factory_mod

use test_fields_3d_mod,            only : scalar_field3d_t, vector_field3d_t
use solid_rotation3d_therm_mod,    only : solid_rotation_theta_Nb_t, &
                                          solid_rotation_PExner_Nb_t,&
                                          solid_rotation_theta_isoT_t, &
                                          solid_rotation_PExner_isoT_t
use solid_rotation_wind_field_mod, only : solid_rotation_wind_field_t
use parcomm_mod,                   only : parcomm_global

implicit none

contains

subroutine create_solid_rotation_field_generators(background_type,&
                                                  u0,omega,sphere_rad,Nb,grav,alpha,&
                                                  theta_gen, Pexner_gen, wind_gen)
    character(len=*) :: background_type
    real(kind=8), intent(in) :: u0,omega,sphere_rad,Nb,grav,alpha
    class(scalar_field3d_t), optional, allocatable, intent(out) :: theta_gen, Pexner_gen
    class(vector_field3d_t), optional, allocatable, intent(out) :: wind_gen

    select case(background_type)
    case("Nb_const")
        if(present(theta_gen)) &
            theta_gen = solid_rotation_theta_Nb_t(u0=u0,omega=omega,sphere_rad=sphere_rad,&
                                                  Nb=Nb,grav=grav,alpha=alpha)
        if(present(PExner_gen)) &
            PExner_gen = solid_rotation_PExner_Nb_t(u0=u0,omega=omega,sphere_rad=sphere_rad,&
                                                   Nb=Nb,grav=grav,alpha=alpha)
        if(present(wind_gen)) &
            wind_gen =  solid_rotation_wind_field_t(u0=u0,alpha=alpha)
    case("isoterm")
        if(present(theta_gen)) &
            theta_gen = solid_rotation_theta_isoT_t(u0=u0,omega=omega,sphere_rad=sphere_rad,&
                                                  T0=300._8,grav=grav,alpha=alpha)
        if(present(PExner_gen)) &
            PExner_gen = solid_rotation_PExner_isoT_t(u0=u0,omega=omega,sphere_rad=sphere_rad,&
                                                   T0=300._8,grav=grav,alpha=alpha)
        if(present(wind_gen)) &
            wind_gen =  solid_rotation_wind_field_t(u0=u0,alpha=alpha)
    case default
        call parcomm_global%abort("create_solid_rotation_field_generators error, "//&
                                  "unknown background_type: "// background_type)
    end select
end subroutine create_solid_rotation_field_generators

end module solid_rotation_fields_factory_mod

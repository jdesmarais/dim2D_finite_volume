      module far_field_perturbation_module

        use parameters_input, only :
     $     obc_perturbation_T0_ac,
     $     obc_perturbation_vx0_ac,
     $     obc_perturbation_vy0_ac,
     $     obc_perturbation_T0_amp,
     $     obc_perturbation_vx0_amp,
     $     obc_perturbation_vy0_amp

        use parameters_kind, only :
     $     ikind,
     $     rkind


        implicit none


        private
        public :: add_far_field_perturbation


        contains


        subroutine add_far_field_perturbation(
     $       temperature,
     $       velocity_x,
     $       velocity_y)

          implicit none

          real(rkind), intent(inout) :: temperature
          real(rkind), intent(inout) :: velocity_x
          real(rkind), intent(inout) :: velocity_y


          if(obc_perturbation_T0_ac) then
             temperature = temperature + obc_perturbation_T0_amp
          end if

          if(obc_perturbation_vx0_ac) then
             velocity_x =  velocity_x + obc_perturbation_vx0_amp
          end if

          if(obc_perturbation_vy0_ac) then
             velocity_y =  velocity_y + obc_perturbation_vy0_amp
          end if

        end subroutine add_far_field_perturbation

      end module far_field_perturbation_module

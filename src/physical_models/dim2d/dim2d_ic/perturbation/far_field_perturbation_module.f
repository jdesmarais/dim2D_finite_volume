      !> @file
      !> add a perturbation to the far-field values
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> add a perturbation to the far-field values
      !
      !> @date
      !> 09_04_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add a perturbation to the far-field values
        !> \f[ T_\infty = 
        !> \begin{cases}
        !> (1 + \epsilon_T) T_0 & \textrm{if temperature perturbation activated} \\\
        !> T_0 & \textrm{otherwise}
        !> \end{cases}
        !> \f]
        !> \f[ u_\infty = 
        !> \begin{cases}
        !> (1 + \epsilon_u) u_0 & \textrm{if x-velocity perturbation activated} \\\
        !> u_0 & \textrm{otherwise}
        !> \end{cases}
        !> \f]
        !> \f[ v_\infty = 
        !> \begin{cases}
        !> (1 + \epsilon_v) v_0 & \textrm{if y-velocity perturbation activated} \\\
        !> v_0 & \textrm{otherwise}
        !> \end{cases}
        !> \f]
        !
        !> @date
        !> 09_04_2015 - initial version - J.L. Desmarais
        !
        !>@param temperature
        !> far-field temperature, \f$ T_\infty \f$
        !
        !>@param velocity_x
        !> far-field x-velocity, \f$ u_\infty \f$
        !
        !>@param velocity_y
        !> far-field y-velocity, \f$ v_\infty \f$
        !---------------------------------------------------------------
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
             velocity_x =  velocity_x  + obc_perturbation_vx0_amp
          end if

          if(obc_perturbation_vy0_ac) then
             velocity_y =  velocity_y  + obc_perturbation_vy0_amp
          end if

        end subroutine add_far_field_perturbation

      end module far_field_perturbation_module

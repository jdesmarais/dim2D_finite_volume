      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for steady state
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for steady state
      !
      !> @date
      !> 11_12_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_class
       
        use dim2d_parameters, only :
     $       cv_r,
     $       T_c

        use dim2d_state_eq_module, only :
     $       get_mass_density_liquid

        use ic_abstract_class, only :
     $       ic_abstract

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none


        private
        public :: ic

        
        !set the initial temperature in the field
        real(rkind), parameter :: T0_degrees = 100
        real(rkind), parameter :: T0         = (T0_degrees+273.15)/T_c


        !> @class ic
        !> class encapsulating operators to set the initial
        !> conditions and the conditions enforced at the edge of the
        !> computational domain for phase separation
        !
        !> @param apply_initial_conditions
        !> set the initial conditions
        !
        !> @param get_mach_ux_infty
        !> get the mach number along the x-direction in the far field
        !
        !> @param get_mach_uy_infty
        !> get the mach number along the y-direction in the far field
        !
        !> @param get_u_in
        !> get the x-component of the velocity at the edge of the
        !> computational domain
        !
        !> @param get_v_in
        !> get the y-component of the velocity at the edge of the
        !> computational domain
        !
        !> @param get_T_in
        !> get the temperature at the edge of the computational
        !> domain
        !
        !> @param get_P_out
        !> get the pressure at the edge of the computational domain
        !---------------------------------------------------------------
        type, extends(ic_abstract) :: ic

          contains

          procedure, nopass :: apply_ic
          procedure, nopass :: get_mach_ux_infty
          procedure, nopass :: get_mach_uy_infty
          procedure, nopass :: get_u_in
          procedure, nopass :: get_v_in
          procedure, nopass :: get_T_in
          procedure, nopass :: get_P_out
          procedure, nopass :: get_far_field

        end type ic


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for a steady state
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param x_map
        !> map of x-coordinates
        !
        !>@param y_map
        !> map of y-coordinates
        !---------------------------------------------------------------
        subroutine apply_ic(nodes,x_map,y_map)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map          

          
          integer(ikind) :: i,j
          real(rkind)    :: d_liq, E_liq
          real(rkind)    :: x_s, y_s

          x_s = x_map(1)
          y_s = y_map(1)

          !< compute the corresponding reduced temperature
          !> reduced liquid mas density and the total energy
          !> at the temperature asked by the user
          d_liq = get_mass_density_liquid(T0)
          E_liq = get_total_energy(d_liq,T0)

          if(rkind.eq.8) then
             do j=1, ny
                do i=1, nx
                   
                   nodes(i,j,1) =  d_liq
                   nodes(i,j,2) =  0.0d0
                   nodes(i,j,3) =  0.0d0
                   nodes(i,j,4) =  E_liq
                   
                end do
             end do
          else
             do j=1, ny
                do i=1, nx
                   
                   nodes(i,j,1) =  d_liq
                   nodes(i,j,2) =  0.0
                   nodes(i,j,3) =  0.0
                   nodes(i,j,4) =  E_liq
                   
                end do
             end do
          end if

        end subroutine apply_ic


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_ux_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical     :: side_s

          side_s = side

          var = 0.0d0

        end function get_mach_ux_infty


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_uy_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical :: side_s

          side_s = side

          var = 0.0d0

        end function get_mach_uy_infty


        !get the x-component of the velocity enforced
        !at the edge of the computational domain
        function get_u_in(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          var = 0.0d0

        end function get_u_in


        !get the y-component of the velocity enforced
        !at the edge of the computational domain
        function get_v_in(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          var = 0.0d0

        end function get_v_in

      
        !get the temperature enforced at the edge of the
        !computational domain
        function get_T_in(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          var = T0

        end function get_T_in


        !get the pressure enforced at the edge of the
        !computational domain
        function get_P_out(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s
          real(rkind) :: mass

          t_s = t
          x_s = x
          y_s = y

          mass = get_mass_density_liquid(T0)

          if(rkind.eq.8) then
             var = 8.0d0*mass*T0/(3.0d0-mass) - 3.0d0*mass**2
          else
             var = 8.0*mass*T0/(3.0-mass) - 3.0*mass**2
          end if

        end function get_P_out


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the governing variables imposed in the far field
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@param t
        !> time
        !
        !>@param x
        !> x-coordinate
        !
        !>@param y
        !> y-coordinate
        !
        !>@return var
        !> governing variables in the far-field
        !--------------------------------------------------------------
        function get_far_field(t,x,y) result(var)

          implicit none

          real(rkind)   , intent(in) :: t
          real(rkind)   , intent(in) :: x
          real(rkind)   , intent(in) :: y
          real(rkind), dimension(ne) :: var


          real(rkind) :: t_s,x_s,y_s
          real(rkind) :: d_liq
          
          t_s = t
          x_s = x
          y_s = y

          d_liq = get_mass_density_liquid(T0)

          if(rkind.eq.8) then

             var(1) = d_liq
             var(2) = 0.0d0
             var(3) = 0.0d0
             var(4) = d_liq*(8.0d0/3.0d0*cv_r*T0-3.0d0*d_liq)

          else

             var(1) = d_liq
             var(2) = 0.0d0
             var(3) = 0.0d0
             var(4) = d_liq*(8.0/3.0*cv_r*T0-3.0*d_liq)

          end if

        end function get_far_field


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the total energy
        !> for a homogeneous liquid given its mass
        !> density and its temperature
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !>@param mass_density
        !> mass density
        !
        !>@param temperature
        !> temperature
        !---------------------------------------------------------------
        function get_total_energy(mass_density, temperature)
     $     result(total_energy)

          implicit none

          real(rkind), intent(in) :: mass_density
          real(rkind), intent(in) :: temperature
          real(rkind)             :: total_energy


          if(rkind.eq.8) then
             total_energy = mass_density*(
     $            8.0d0/3.0d0*cv_r*temperature - 3.0d0*mass_density)
          else
             total_energy = mass_density*(
     $            8.0/3.0*cv_r*temperature - 3.0*mass_density)
          end if

        end function get_total_energy

      end module ic_class

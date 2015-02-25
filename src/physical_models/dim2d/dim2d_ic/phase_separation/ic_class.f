      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for phase separation
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for phase separation
      !
      !> @date
      !> 11_12_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_class
       
        use dim2d_parameters, only :
     $       cv_r,
     $       we

        use dim2d_state_eq_module, only :
     $       get_mass_density_liquid,
     $       get_mass_density_vapor

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
        real(rkind), parameter :: T0 = 0.99d0


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
          procedure,   pass :: get_far_field

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

          
          !local variables for the droplet/bubble
          real(rkind) :: dliq,dvap


          !< local variables for the perturbation
          real(rkind) :: xMin, xMax
          real(rkind) :: Tx, kx, Ax
          
          real(rkind) :: yMin, yMax
          real(rkind) :: Ty, ky, Ay

          !< local variables for the initialization
          integer(ikind) :: i,j
          real(rkind)    :: x,y
          real(rkind)    :: noise


          !<get the mass densities corresponding to the
          !>liquid and vapor phases for the initial
          !>temperature field
          dliq = get_mass_density_liquid(T0)
          dvap = get_mass_density_vapor(T0)

          !<set the perturbation properties
          if(rkind.eq.8) then
             xMin = -0.5d0
             xMax =  0.5d0
             Tx   =  xMax - xMin
             kx   =  2.0d0 * acos(-1.0d0)/Tx
             Ax   =  0.7d0*(dliq-dvap)
          else
             xMin = -0.5
             xMax =  0.5
             Tx   =  xMax - xMin
             kx   =  2.0 * acos(-1.0)/Tx
             Ax   =  0.7*(dliq-dvap)
          end if

          if(rkind.eq.8) then
             yMin = -0.5d0
             yMax =  0.5d0
             Ty   =  yMax - yMin
             ky   =  2.0d0 * acos(-1.0d0)/Ty
             Ay   =  0.7d0*(dliq-dvap)
          else
             yMin = -0.5
             yMax =  0.5
             Ty   =  yMax - yMin
             ky   =  2.0 * acos(-1.0)/Ty
             Ay   =  0.7*(dliq-dvap)
          end if


          !<initialize the fields
          do j=1, ny
             !DEC$ IVDEP
             do i=1, nx

                !<coordinates
                x = x_map(i)
                y = y_map(j)

                !<unstable mass density
                if(rkind.eq.8) then
                   nodes(i,j,1)=(dliq+dvap)/2.0d0
                else
                   nodes(i,j,1)=(dliq+dvap)/2.0
                end if

                !<adding the sinusoidal perturbation to 
                !the initial unstable mass density
                noise = perturbation(x,xMin,xMax,kx,Ax)*
     $               perturbation(y,yMin,yMax,ky,Ay)
                nodes(i,j,1)=nodes(i,j,1)+noise

                !<null velocity field
                nodes(i,j,2)=0.0d0
                nodes(i,j,3)=0.0d0

                !<total energy field corresponding to the
                !<temperature and the mass density fields
                nodes(i,j,4)=energy_phase_separation(
     $               x,y,
     $               nodes(i,j,1),T0,
     $               xMin,xMax,kx,Ax,
     $               yMin,yMax,ky,Ay)

             end do
          end do

        end subroutine apply_ic


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_ux_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical :: side_s

          side_s = side

          if(rkind.eq.8) then
             var = 0.0d0
          else
             var = 0.0
          end if

        end function get_mach_ux_infty


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_uy_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical :: side_s

          side_s = side

          if(rkind.eq.8) then
             var = 0.0d0
          else
             var = 0.0
          end if

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

          if(rkind.eq.8) then
             var = 0.0d0
          else
             var = 0.0
          end if

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

          if(rkind.eq.8) then
             var = 0.0d0
          else
             var = 0.0
          end if

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

          mass = get_mass_density_vapor(T0)

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
        function get_far_field(this,t,x,y) result(var)

            implicit none

            class(ic)     , intent(in) :: this
            real(rkind)   , intent(in) :: t
            real(rkind)   , intent(in) :: x
            real(rkind)   , intent(in) :: y
            real(rkind), dimension(ne) :: var


            real(rkind) :: mass

            mass = get_mass_density_vapor(T0)

            if(rkind.eq.8) then

               var(1) = mass
               var(2) = var(1)*get_u_in(t,x,y)
               var(3) = var(1)*get_v_in(t,x,y)
               var(4) = 0.5d0*(var(2)**2+var(3)**2)/var(1) +
     $                  var(1)*(8.0d0/3.0d0*cv_r*T0-3.0d0*mass)

            else

               var(1) = mass
               var(2) = var(1)*get_u_in(t,x,y)
               var(3) = var(1)*get_v_in(t,x,y)
               var(4) = 0.5*(var(2)**2+var(3)**2)/var(1) +
     $                  var(1)*(8.0/3.0*cv_r*T0-3.0*mass)

            end if

          end function get_far_field


          !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the perturbation
        !
        !> @date
        !> 19_12_2013 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !>@param x_min
        !> lower space border of the perturbation
        !>@param x_max
        !> upper space border of the perturbation
        !>@param kx
        !> wave number for the sinusoidal perturbation
        !>@param Ax
        !> amplitude of the perturbation
        !---------------------------------------------------------------
        function perturbation(x,x_min,x_max,kx,Ax)

          implicit none

          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: x_min
          real(rkind), intent(in) :: x_max
          real(rkind), intent(in) :: kx
          real(rkind), intent(in) :: Ax
          real(rkind)             :: perturbation
          

          if((x.lt.x_min).or.(x.gt.x_max)) then
             if(rkind.eq.8) then
                perturbation = 0.0d0
             else
                perturbation = 0.0
             end if

          else
             if(rkind.eq.8) then
                perturbation = Ax*(
     $               1.0d0+cos(kx*(x-(x_max+x_min)/2.0d0)))
             else
                perturbation = Ax*(
     $               1.0+cos(kx*(x-(x_max+x_min)/2.0)))
             end if
          end if


        end function perturbation

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the perturbation
        !
        !> @date
        !> 19_12_2013 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !>@param x_min
        !> lower space border of the perturbation
        !>@param x_max
        !> upper space border of the perturbation
        !>@param kx
        !> wave number for the sinusoidal perturbation
        !>@param Ax
        !> amplitude of the perturbation
        !---------------------------------------------------------------
        function perturbation_gradient(x,x_min,x_max,kx,Ax)

          implicit none

          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: x_min
          real(rkind), intent(in) :: x_max
          real(rkind), intent(in) :: kx
          real(rkind), intent(in) :: Ax
          real(rkind)             :: perturbation_gradient
          

          if((x.lt.x_min).or.(x.gt.x_max)) then
             if(rkind.eq.8) then
                perturbation_gradient = 0.0d0
             else
                perturbation_gradient = 0.0
             end if

          else
             if(rkind.eq.8) then
                perturbation_gradient = - kx*Ax*
     $               sin(kx*(x-(x_max+x_min)/2.0d0))
             else
                perturbation_gradient = - kx*Ax*
     $               sin(kx*(x-(x_max+x_min)/2.0))
             end if
          end if

        end function perturbation_gradient


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the total energy corresponding
        !> to the initial mass density and temperature fields
        !> leading to phase separation
        !
        !> @date
        !> 19_12_2013 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !>@param y
        !> y-coordinate
        !>@param mass_density
        !> mass density at (x,y)
        !>@param temperature
        !> temperature at (x,y)
        !>@param x_min
        !> lower space border of the perturbation along x
        !>@param x_max
        !> upper space border of the perturbation along x
        !>@param kx
        !> wave number for the sinusoidal perturbation along x
        !>@param Ax
        !> amplitude of the perturbation along x
        !>@param y_min
        !> lower space border of the perturbation along y
        !>@param y_max
        !> upper space border of the perturbation along y
        !>@param ky
        !> wave number for the sinusoidal perturbation along y
        !>@param Ay
        !> amplitude of the perturbation along y
        !---------------------------------------------------------------
        function energy_phase_separation(
     $               x,y,
     $               mass_density,temperature,
     $               x_min,x_max,kx,Ax,
     $               y_min,y_max,ky,Ay)

          implicit none

          real(rkind), intent(in) :: x,y
          real(rkind), intent(in) :: mass_density, temperature
          real(rkind), intent(in) :: x_min, x_max, kx, Ax
          real(rkind), intent(in) :: y_min, y_max, ky, Ay
          real(rkind)             :: energy_phase_separation


          real(rkind) :: d_md_dx, d_md_dy


          !<compute the mass density gradients at (x,y)
          d_md_dx = perturbation_gradient(x,x_min,x_max,kx,Ax)
          d_md_dy = perturbation_gradient(y,y_min,y_max,ky,Ay)

          
          !<compute the total energy assuming no velocity field
          if(rkind.eq.8) then
             energy_phase_separation=
     $            mass_density*(
     $            8.0d0/3.0d0*cv_r*temperature-3.0d0*mass_density)
     $            + 1.0d0/(2.0d0*we)*((d_md_dx)**2+(d_md_dy)**2)
          else
             energy_phase_separation=
     $            mass_density*(
     $            8.0/3.0*cv_r*temperature-3.0*mass_density)
     $            + 1.0/(2.0*we)*((d_md_dx)**2+(d_md_dy)**2)
          end if

        end function energy_phase_separation

      end module ic_class

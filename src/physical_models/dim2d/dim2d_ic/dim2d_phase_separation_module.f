      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> the initial conditions of unstable initial
      !> density with a small perturbation leading
      !> to phase separation
      !
      !> @date
      !> 19_12_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_phase_separation_module

        use dim2d_parameters       , only : cv_r,we
        use dim2d_state_eq_module  , only : get_mass_density_liquid,
     $                                      get_mass_density_vapor
        use parameters_input       , only : nx,ny,ne
        use parameters_kind        , only : ikind, rkind

        implicit none

        private
        public :: apply_phase_separation_ic


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> leading to phase separation
        !
        !> @date
        !> 19_12_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !---------------------------------------------------------------
        subroutine apply_phase_separation_ic(nodes,x_map,y_map)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map


          !< local variables for the unstable state
          real(rkind)    :: T0
          real(rkind)    :: dliq,dvap


          !< local variables for the perturbation
          real(rkind) :: xMin, xMax
          real(rkind) :: Tx, kx, Ax
          
          real(rkind) :: yMin, yMax
          real(rkind) :: Ty, ky, Ay

          !< local variables for the initialization
          integer(ikind) :: i,j
          real(rkind)    :: x,y
          real(rkind)    :: noise


          !<set the initial temperature in the field
          T0 = 0.99d0

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

        end subroutine apply_phase_separation_ic


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

      end module dim2d_phase_separation_module


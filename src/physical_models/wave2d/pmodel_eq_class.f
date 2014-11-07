      !> @file
      !> class encapsulating subroutines to compute
      !> the wave governing equations in 2D
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the wave governing equations in 2D
      !
      !> @date
      !> 28_07_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module pmodel_eq_class
      
        use interface_primary, only :
     $     gradient_x_proc,
     $     gradient_y_proc

        use parameters_bf_layer, only :
     $       bc_interior_pt,
     $       interior_pt

        use parameters_constant, only :
     $       scalar,
     $       vector_x,
     $       vector_y,
     $       oscillatory_forcing

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size,
     $       wave_forcing

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_default_class, only :
     $       pmodel_eq_default

        use sd_operators_class, only :
     $       sd_operators

        use wave2d_ncoords_module, only :
     $       compute_n_gradient_wave2d,
     $       compute_n_eigenvalues_wave2d,
     $       compute_n1_lefteigenvector_wave2d,
     $       compute_n1_righteigenvector_wave2d,
     $       compute_n2_lefteigenvector_wave2d,
     $       compute_n2_righteigenvector_wave2d

        use wave2d_parameters, only :
     $       c,
     $       c_x,
     $       c_y,
     $       epsilon,
     $       x_center,
     $       y_center,
     $       amplitude,
     $       period,
     $       amplitude_force,
     $       period_force,
     $       x_center_force,
     $       y_center_force

        use wave2d_prim_module, only :
     $       position,
     $       velocity_x,
     $       velocity_y
        
        implicit none

        private
        public :: pmodel_eq, peak


        !> @class pmodel_eq
        !> class encapsulating operators to compute
        !> the wave governing equations in 2D
        !>
        !> @param get_model_name
        !> get the name of the physical model
        !>
        !> @param get_var_name
        !> get the name of the main variables
        !> (position, x-velocity, y-velocity)
        !>
        !> @param get_var_longname
        !> get the description of the main variables for the
        !> governing equations of the physical model
        !>
        !> @param get_var_unit
        !> get the units of the main variables
        !>
        !> @param get_var_types
        !> get the type of the main variables
        !> (scalar, vector-x, vector-y)
        !>
        !> @param get_eq_nb
        !> get the number of governing equations: 3
        !>
        !> @param apply_initial_conditions
        !> initialize the main variables
        !>
        !> @param compute_fluxes
        !> compute the fluxes along the x- and y-axis
        !
        !> @param compute_body_forces
        !> compute the forcing term
        !
        !> @param are_openbc_undermined
        !> check whether the open boundary conditions are undermined
        !> at the grid point location
        !---------------------------------------------------------------
        type, extends(pmodel_eq_default) :: pmodel_eq
          
          contains

          procedure, nopass :: get_model_name
          procedure, nopass :: get_var_name
          procedure, nopass :: get_var_longname
          procedure, nopass :: get_var_unit
          procedure, nopass :: get_var_type
          procedure, nopass :: get_eq_nb

          procedure,   pass :: apply_ic

          procedure, nopass :: compute_flux_x
          procedure, nopass :: compute_flux_y
          procedure, nopass :: compute_flux_x_nopt
          procedure, nopass :: compute_flux_y_nopt
          procedure, nopass :: compute_flux_x_oneside
          procedure, nopass :: compute_flux_y_oneside
          procedure, nopass :: compute_body_forces

          procedure, nopass :: get_velocity
          procedure, nopass :: are_openbc_undermined

          procedure, nopass :: compute_x_eigenvalues
          procedure, nopass :: compute_y_eigenvalues
          procedure, nopass :: compute_x_lefteigenvector
          procedure, nopass :: compute_x_righteigenvector
          procedure, nopass :: compute_y_lefteigenvector
          procedure, nopass :: compute_y_righteigenvector

          procedure, nopass :: compute_n1_eigenvalues       => compute_n_eigenvalues_wave2d
          procedure, nopass :: compute_n2_eigenvalues       => compute_n_eigenvalues_wave2d
          procedure, nopass :: compute_n1_lefteigenvector   => compute_n1_lefteigenvector_wave2d
          procedure, nopass :: compute_n1_righteigenvector  => compute_n1_righteigenvector_wave2d
          procedure, nopass :: compute_n2_lefteigenvector   => compute_n2_lefteigenvector_wave2d
          procedure, nopass :: compute_n2_righteigenvector  => compute_n2_righteigenvector_wave2d

          procedure, nopass :: compute_x_gradient
          procedure, nopass :: compute_y_gradient
          procedure, nopass :: compute_n_gradient => compute_n_gradient_wave2d

        end type pmodel_eq


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to get the name of the physical model
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param model_name
        !> character giving the name of the model
        !---------------------------------------------------------------
        function get_model_name() result(model_name)

          implicit none

          character(len=10) :: model_name

          model_name="Wave2D"

        end function get_model_name
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the name of the main variables
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param var_name
        !> characters giving the variable names
        !---------------------------------------------------------------
        function get_var_name() result(var_pties)

          implicit none

          character(len=10), dimension(ne) :: var_pties

          var_pties(1)="position"
          var_pties(2)="x-velocity"
          var_pties(3)="y-velocity"

        end function get_var_name
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the name of the main variables
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param var_name
        !> characters giving the variable names
        !---------------------------------------------------------------
        function get_var_longname() result(var_pties)

          implicit none

          character(len=33), dimension(ne) :: var_pties

          var_pties(1)="position"
          var_pties(2)="velocity along x-axis"
          var_pties(3)="velocity along y-axis"

        end function get_var_longname


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the units of the main variables
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param var_name
        !> characters giving the variable units
        !---------------------------------------------------------------
        function get_var_unit() result(var_pties)

          implicit none

          character(len=23), dimension(ne) :: var_pties

          var_pties(1)= "m/m"
          var_pties(2)= "(m/s)/(m/s)"
          var_pties(3)= "(m/s)/(m/s)"

        end function get_var_unit


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to get the type of the main variables
        !> (scalar, vector_x, vector_y, scalar)
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param var_name
        !> characters giving the variable type
        !---------------------------------------------------------------
        function get_var_type() result(var_type)

          implicit none

          integer, dimension(ne) :: var_type

          var_type(1)=scalar
          var_type(2)=vector_x
          var_type(3)=vector_y

        end function get_var_type
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the number of main variables
        !> in the governing equations: 4
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param eq_nb
        !> number of governing equations
        !---------------------------------------------------------------
        function get_eq_nb() result(eq_nb)
          implicit none
          integer :: eq_nb
          eq_nb=3
        end function get_eq_nb
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the initial conditions to the main
        !> variables of the governing equations
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param x_map
        !> array with the x-ccordinates
        !
        !>@param nodes
        !> array with the y-coordinates
        !---------------------------------------------------------------
        subroutine apply_ic(this,nodes,x_map,y_map)

          implicit none

          class(pmodel_eq)             , intent(in)    :: this
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map


          integer(ikind) :: i,j
          integer :: neq

          neq = this%get_eq_nb()


          if(rkind.eq.8) then

             do j=1, size(y_map,1)
                do i=1, size(x_map,1)

                   nodes(i,j,1) = peak(amplitude,
     $                                 period,
     $                                 x_map(i)-x_center,
     $                                 y_map(j)-y_center)
                   nodes(i,j,2) = 0.0d0
                   nodes(i,j,3) = 0.0d0
                   
                end do
             end do
             
          else
             
             do j=1, size(y_map,1)
                do i=1, size(x_map,1)

                   nodes(i,j,1) = peak(amplitude,
     $                                 period,
     $                                 x_map(i)-x_center,
     $                                 y_map(j)-y_center)
                   nodes(i,j,2) = 0.0
                   nodes(i,j,3) = 0.0
                   
                end do
             end do
             
          end if


c$$$          if(.true.) then
c$$$             nodes(1,1,1) = 0.2
c$$$             nodes(2,1,1) = 2.3
c$$$             nodes(3,1,1) = -6.3
c$$$             nodes(4,1,1) = 7.8
c$$$             nodes(5,1,1) = -4.1
c$$$             nodes(6,1,1) = 6.9
c$$$
c$$$             nodes(1,2,1) = 3.6
c$$$             nodes(2,2,1) = 5.2
c$$$             nodes(3,2,1) = 2.7
c$$$             nodes(4,2,1) = -5.23
c$$$             nodes(5,2,1) = 1.2
c$$$             nodes(6,1,1) = -7.26
c$$$
c$$$             nodes(1,3,1) = 9.26
c$$$             nodes(2,3,1) = -3.2
c$$$             nodes(3,3,1) = 7.89
c$$$             nodes(4,3,1) = -2.3
c$$$             nodes(5,3,1) = 8.62
c$$$             nodes(6,3,1) = 6.23
c$$$
c$$$             nodes(1,4,1) = 0.23
c$$$             nodes(2,4,1) = 9.26
c$$$             nodes(3,4,1) = 7.12
c$$$             nodes(4,4,1) = -5.6
c$$$             nodes(5,4,1) = -6.98
c$$$             nodes(6,4,1) = -4.56
c$$$
c$$$             nodes(1,5,1) = 1.25
c$$$             nodes(2,5,1) =-4.53
c$$$             nodes(3,5,1) = 9.42
c$$$             nodes(4,5,1) = 8.15
c$$$             nodes(5,5,1) = -8.152
c$$$             nodes(6,5,1) = 9.785
c$$$
c$$$             nodes(1,6,1) = 3.12
c$$$             nodes(2,6,1) = 5.146
c$$$             nodes(3,6,1) =-7.584
c$$$             nodes(4,6,1) = 2.36
c$$$             nodes(5,6,1) = 1.20
c$$$             nodes(6,6,1) = -8.12
c$$$
c$$$             do j=1,6
c$$$                do i=1,6
c$$$                   nodes(i,j,2) = nodes(i,j,1)+1.0
c$$$                   nodes(i,j,3) = nodes(i,j,1)*10-5
c$$$                end do
c$$$             end do
c$$$          end if

        end subroutine apply_ic



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the peak needed for the initial conditions
        !> of the wave 2d
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param amplitude
        !> amplitude of the peak
        !
        !>@param period
        !> characteristic size for the compact support of the peak
        !
        !>@param x
        !> x-coordinate identifying where the initial condition is
        !> evaluated
        !
        !>@param y
        !> y-coordinate identifying where the initial condition is
        !> evaluated
        !---------------------------------------------------------------
        function peak(amplitude, period, x, y)

          implicit none

          real(rkind), intent(in) :: amplitude
          real(rkind), intent(in) :: period
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: peak


          real(rkind) :: radius
          real(rkind) :: radius_max
          real(rkind) :: omega
          

          radius = Sqrt(x**2+y**2)

          if(rkind.eq.8) then

             radius_max = period/2.0d0
             omega      = 2.0d0*ACOS(-1.0d0)/period

             if(radius.le.radius_max) then
                peak = amplitude*(1.0d0 + cos(omega*radius))
             else
                peak = 0.0d0                
             end if

          else

             radius_max = period/2.0
             omega      = 2.0*ACOS(-1.0)/period

             if(radius.le.radius_max) then
                peak = amplitude*(1.0 + cos(omega*radius))
             else
                peak = 0.0
             end if
             
          end if

        end function peak
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the fluxes along the x-axis
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param s
        !> space discretization operators
        !
        !>@param flux_x
        !> fluxes along the x-axis
        !---------------------------------------------------------------
        function compute_flux_x(nodes,dx,dy,s) result(flux_x)
        
          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)   :: nodes
          real(rkind)                       , intent(in)   :: dx
          real(rkind)                       , intent(in)   :: dy
          type(sd_operators)                , intent(in)   :: s
          real(rkind), dimension(nx+1,ny,ne)               :: flux_x

          integer     :: i,j
          real(rkind) :: dy_s

          dy_s = dy

          !<fluxes along the x-axis
          do j=bc_size+1, ny-bc_size
             do i=bc_size+1, nx+1-bc_size

                flux_x(i,j,1) = - c**2*s%f(nodes,i,j,velocity_x)
     $                          - epsilon*s%dfdx(nodes,i,j,position,dx)

                flux_x(i,j,2) = - c**2*s%f(nodes,i,j,position)
     $                          - epsilon*s%dfdx(nodes,i,j,velocity_x,dx)

                flux_x(i,j,3) = - c_x*s%f(nodes,i,j,velocity_y)
     $                          - epsilon*s%dfdx(nodes,i,j,velocity_y,dx)

             end do
          end do

        end function compute_flux_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the fluxes along the y-axis
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param s
        !> space discretization operators
        !
        !>@param flux_y
        !> fluxes along the y-axis
        !---------------------------------------------------------------
        function compute_flux_y(nodes,dx,dy,s) result(flux_y)
        
          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)   :: nodes
          real(rkind)                       , intent(in)   :: dx
          real(rkind)                       , intent(in)   :: dy
          type(sd_operators)                , intent(in)   :: s
          real(rkind), dimension(nx,ny+1,ne)               :: flux_y

          integer     :: i,j
          real(rkind) :: dx_s

          dx_s = dx


          !<fluxes along the x-axis
          do j=bc_size+1, ny+1-bc_size
             do i=bc_size+1, nx-bc_size

                flux_y(i,j,1) = - c**2*s%g(nodes,i,j,velocity_y)
     $                          - epsilon*s%dgdy(nodes,i,j,position,dy)

                flux_y(i,j,2) = - c_y*s%g(nodes,i,j,velocity_x)
     $                          - epsilon*s%dgdy(nodes,i,j,velocity_x,dy)

                flux_y(i,j,3) = - c**2*s%g(nodes,i,j,position)
     $                          - epsilon*s%dgdy(nodes,i,j,velocity_y,dy)

             end do
          end do

        end function compute_flux_y


        subroutine compute_flux_x_nopt(
     $     nodes,dx,dy,s,
     $     grdpts_id,
     $     flux_x,
     $     x_borders,
     $     y_borders)
        
          implicit none

          real(rkind)   , dimension(:,:,:), intent(in)    :: nodes
          real(rkind)                     , intent(in)    :: dx
          real(rkind)                     , intent(in)    :: dy
          type(sd_operators)              , intent(in)    :: s
          integer       , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind)   , dimension(:,:,:), intent(inout) :: flux_x
          integer(ikind), dimension(2)    , intent(in)    :: x_borders
          integer(ikind), dimension(2)    , intent(in)    :: y_borders

          integer(ikind) :: i,j
          real(rkind)    :: dy_s

          dy_s = dy

          !<fluxes along the x-axis
          do j=y_borders(1), y_borders(2)
             !DEC$ IVDEP
             do i=x_borders(1), x_borders(2)+1

                if((grdpts_id(i,j).eq.interior_pt).or.
     $               (grdpts_id(i,j).eq.bc_interior_pt)) then

                   flux_x(i,j,1) = - c**2*s%f(nodes,i,j,velocity_x)
     $                             - epsilon*s%dfdx(nodes,i,j,position,dx)

                   flux_x(i,j,2) = - c**2*s%f(nodes,i,j,position)
     $                             - epsilon*s%dfdx(nodes,i,j,velocity_x,dx)

                   flux_x(i,j,3) = - c_x*s%f(nodes,i,j,velocity_y)
     $                             - epsilon*s%dfdx(nodes,i,j,velocity_y,dx)

                end if

             end do
          end do

        end subroutine compute_flux_x_nopt


        subroutine compute_flux_y_nopt(
     $     nodes,dx,dy,s,
     $     grdpts_id,
     $     flux_y,
     $     x_borders,
     $     y_borders)
        
          implicit none

          real(rkind)   , dimension(:,:,:), intent(in)    :: nodes
          real(rkind)                     , intent(in)    :: dx
          real(rkind)                     , intent(in)    :: dy
          type(sd_operators)              , intent(in)    :: s
          integer       , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind)   , dimension(:,:,:), intent(inout) :: flux_y
          integer(ikind), dimension(2)    , intent(in)    :: x_borders
          integer(ikind), dimension(2)    , intent(in)    :: y_borders

          integer(ikind) :: i,j
          real(rkind)    :: dx_s

          dx_s = dx

          !<fluxes along the y-axis
          do j=y_borders(1), y_borders(2)+1
             !DEC$ IVDEP
             do i=x_borders(1), x_borders(2)

                 if((grdpts_id(i,j).eq.interior_pt).or.
     $               (grdpts_id(i,j).eq.bc_interior_pt)) then

                   flux_y(i,j,1) = - c**2*s%g(nodes,i,j,velocity_y)
     $                             - epsilon*s%dgdy(nodes,i,j,position,dy)

                   flux_y(i,j,2) = - c_y*s%g(nodes,i,j,velocity_x)
     $                             - epsilon*s%dgdy(nodes,i,j,velocity_x,dy)

                   flux_y(i,j,3) = - c**2*s%g(nodes,i,j,position)
     $                             - epsilon*s%dgdy(nodes,i,j,velocity_y,dy)

                end if

             end do
          end do

        end subroutine compute_flux_y_nopt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the fluxes along the x-axis
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param i
        !> x-index where the flux_x is computed
        !
        !>@param j
        !> y-index where the flux_x is computed
        !
        !>@param s_oneside
        !> space discretization operators
        !
        !>@param flux_x
        !> fluxes along the x-axis
        !---------------------------------------------------------------
        function compute_flux_x_oneside(nodes,dx,dy,i,j,s_oneside)
     $     result(flux_x)
        
          implicit none

          real(rkind), dimension(:,:,:), intent(in)   :: nodes
          real(rkind)                  , intent(in)   :: dx
          real(rkind)                  , intent(in)   :: dy
          integer(ikind)               , intent(in)   :: i
          integer(ikind)               , intent(in)   :: j
          class(sd_operators)          , intent(in)   :: s_oneside
          real(rkind), dimension(ne)                  :: flux_x

          real(rkind) :: dy_s

          dy_s = dy

          !<fluxes along the x-axis
          flux_x(1) = - c**2*s_oneside%f(nodes,i,j,velocity_x)
     $                - epsilon*s_oneside%dfdx(nodes,i,j,position,dx)

          flux_x(2) = - c**2*s_oneside%f(nodes,i,j,position)
     $                - epsilon*s_oneside%dfdx(nodes,i,j,velocity_x,dx)
          
          flux_x(3) = - c_x*s_oneside%f(nodes,i,j,velocity_y)
     $                - epsilon*s_oneside%dfdx(nodes,i,j,velocity_y,dx)

        end function compute_flux_x_oneside


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the fluxes along the y-axis
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param i
        !> x-index where the flux_x is computed
        !
        !>@param j
        !> y-index where the flux_x is computed
        !
        !>@param s_oneside
        !> space discretization operators
        !
        !>@param flux_y
        !> fluxes along the y-axis
        !---------------------------------------------------------------
        function compute_flux_y_oneside(nodes,dx,dy,i,j,s_oneside)
     $     result(flux_y)
        
          implicit none

          real(rkind), dimension(:,:,:), intent(in)   :: nodes
          real(rkind)                  , intent(in)   :: dx
          real(rkind)                  , intent(in)   :: dy
          integer(ikind)               , intent(in)   :: i
          integer(ikind)               , intent(in)   :: j
          class(sd_operators)          , intent(in)   :: s_oneside
          real(rkind), dimension(ne)                  :: flux_y

          real(rkind) :: dx_s

          dx_s = dx


          !<fluxes along the x-axis
          flux_y(1) = - c**2*s_oneside%g(nodes,i,j,velocity_y)
     $                - epsilon*s_oneside%dgdy(nodes,i,j,position,dy)
          
          flux_y(2) = - c_y*s_oneside%g(nodes,i,j,velocity_x)
     $                - epsilon*s_oneside%dgdy(nodes,i,j,velocity_x,dy)
          
          flux_y(3) = - c**2*s_oneside%g(nodes,i,j,position)
     $                - epsilon*s_oneside%dgdy(nodes,i,j,velocity_y,dy)

        end function compute_flux_y_oneside


        function compute_body_forces(t,x,y,nodes,k) result(body_forces)

          implicit none

          real(rkind)               , intent(in) :: t
          real(rkind)               , intent(in) :: x
          real(rkind)               , intent(in) :: y
          real(rkind), dimension(ne), intent(in) :: nodes
          integer                   , intent(in) :: k
          real(rkind)                            :: body_forces

          real(rkind) :: omega
          real(rkind) :: node_s

          if(wave_forcing.eq.oscillatory_forcing) then

             if(k.eq.1) then
                omega       = 2.0d0*ACOS(-1.0d0)/period_force
                body_forces = peak(
     $               amplitude_force*SIN(omega*t),
     $               period,
     $               x-x_center_force,
     $               y-y_center_force)
             else
                body_forces = 0
             end if

          else
             
             body_forces = 0

          end if

          node_s = nodes(1)

        end function compute_body_forces


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to compute the body forces
        !> acting on the cell
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables at the grid point location
        !
        !>@param velocity
        !> velocity vector at the grid point location
        !--------------------------------------------------------------
        function get_velocity(nodes) result(velocity)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(2)              :: velocity

          velocity(1) = nodes(2)
          velocity(2) = nodes(3)

        end function get_velocity


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the open boundary conditions
        !> are undermined at the grid point location
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param undermined
        !> check if the open boundary conditions are undermined
        !> at the grid point location
        !--------------------------------------------------------------
        function are_openbc_undermined(nodes) result(undermined)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          logical                                :: undermined

          real(rkind) :: node_s

          !real(rkind) :: d_liq, d_vap
          !
          !d_liq = 1.1-0.05*(1.1-0.1)
          !d_vap = 0.1+0.05*(1.1-0.1)
          !
          !if((nodes(1).ge.d_vap).and.(nodes(1).le.d_liq)) then
          !   undermined = .true.
          !else
          !   undermined = .false.
          !end if

          node_s = nodes(1)

          undermined = .false.

        end function are_openbc_undermined


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the x-direction
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_x_eigenvalues(nodes) result(eigenvalues)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne)             :: eigenvalues


          real(rkind) :: node_s

          node_s = nodes(1)

          if(rkind.eq.8) then
             eigenvalues(1) = 0.0d0
          else
             eigenvalues(1) = 0.0
          end if
          eigenvalues(2) = -c**2
          eigenvalues(3) =  c**2

        end function compute_x_eigenvalues


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the y-direction
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_y_eigenvalues(nodes) result(eigenvalues)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne)             :: eigenvalues

          real(rkind) :: node_s

          node_s = nodes(1)

          if(rkind.eq.8) then
             eigenvalues(1) = 0.0d0
          else
             eigenvalues(1) = 0.0
          end if
          eigenvalues(2) = -c**2
          eigenvalues(3) =  c**2

        end function compute_y_eigenvalues


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenvectors for the hyperbolic terms
        !> in the x-direction
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvectors at the location of the grid point
        !--------------------------------------------------------------
        function compute_x_lefteigenvector(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind) :: node_s

          node_s = nodes(1)


          if(rkind.eq.8) then
             eigenvect(1,1) =  0.0d0
             eigenvect(2,1) =  0.0d0
             eigenvect(3,1) =  1.0d0
                               
             eigenvect(1,2) =  0.5d0
             eigenvect(2,2) =  0.5d0
             eigenvect(3,2) =  0.0d0
             
             eigenvect(1,3) = -0.5d0
             eigenvect(2,3) =  0.5d0
             eigenvect(3,3) =  0.0d0
             
          else
             eigenvect(1,1) =  0.0
             eigenvect(2,1) =  0.0
             eigenvect(3,1) =  1.0
                               
             eigenvect(1,2) =  0.5
             eigenvect(2,2) =  0.5
             eigenvect(3,2) =  0.0

             eigenvect(1,3) = -0.5
             eigenvect(2,3) =  0.5
             eigenvect(3,3) =  0.0

          end if

        end function compute_x_lefteigenvector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the right eigenvectors for the hyperbolic terms
        !> in the x-direction
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvectors at the location of the grid point
        !--------------------------------------------------------------
        function compute_x_righteigenvector(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind) :: node_s

          node_s = nodes(1)

          if(rkind.eq.8) then
             eigenvect(1,1) =  0.0d0
             eigenvect(2,1) =  1.0d0
             eigenvect(3,1) = -1.0d0
             
             eigenvect(1,2) =  0.0d0
             eigenvect(2,2) =  1.0d0
             eigenvect(3,2) =  1.0d0
                               
             eigenvect(1,3) =  1.0d0
             eigenvect(2,3) =  0.0d0
             eigenvect(3,3) =  0.0d0
             
          else
             eigenvect(1,1) =  0.0
             eigenvect(2,1) =  1.0
             eigenvect(3,1) = -1.0
             
             eigenvect(1,2) =  0.0
             eigenvect(2,2) =  1.0
             eigenvect(3,2) =  1.0
                               
             eigenvect(1,3) =  1.0
             eigenvect(2,3) =  0.0
             eigenvect(3,3) =  0.0

          end if

        end function compute_x_righteigenvector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenvectors for the hyperbolic terms
        !> in the x-direction
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvectors at the location of the grid point
        !--------------------------------------------------------------
        function compute_y_lefteigenvector(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind) :: node_s

          node_s = nodes(1)


          if(rkind.eq.8) then
             eigenvect(1,1) =  0.0d0
             eigenvect(2,1) =  1.0d0
             eigenvect(3,1) =  0.0d0
                               
             eigenvect(1,2) =  0.5d0
             eigenvect(2,2) =  0.0d0
             eigenvect(3,2) =  0.5d0

             eigenvect(1,3) = -0.5d0
             eigenvect(2,3) =  0.0d0
             eigenvect(3,3) =  0.5d0

          else
             eigenvect(1,1) =  0.0
             eigenvect(2,1) =  1.0
             eigenvect(3,1) =  0.0
                               
             eigenvect(1,2) =  0.5
             eigenvect(2,2) =  0.0
             eigenvect(3,2) =  0.5

             eigenvect(1,3) = -0.5
             eigenvect(2,3) =  0.0
             eigenvect(3,3) =  0.5
          end if

        end function compute_y_lefteigenvector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the right eigenvector for the hyperbolic terms
        !> in the x-direction
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvectors at the location of the grid point
        !--------------------------------------------------------------
        function compute_y_righteigenvector(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind) :: node_s

          node_s = nodes(1)

          if(rkind.eq.8) then
             eigenvect(1,1) =  0.0d0
             eigenvect(2,1) =  1.0d0
             eigenvect(3,1) = -1.0d0

             eigenvect(1,2) =  1.0d0
             eigenvect(2,2) =  0.0d0
             eigenvect(3,2) =  0.0d0
                               
             eigenvect(1,3) =  0.0d0
             eigenvect(2,3) =  1.0d0
             eigenvect(3,3) =  1.0d0

          else
             eigenvect(1,1) =  0.0
             eigenvect(2,1) =  1.0
             eigenvect(3,1) = -1.0

             eigenvect(1,2) =  1.0
             eigenvect(2,2) =  0.0
             eigenvect(3,2) =  0.0
                               
             eigenvect(1,3) =  0.0
             eigenvect(2,3) =  1.0
             eigenvect(3,3) =  1.0
          end if

        end function compute_y_righteigenvector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface for the computation of the gradient of the
        !> governing variables in the x-direction
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> integer identifying the index in the x-direction
        !
        !>@param j
        !> integer identifying the index in the y-direction
        !
        !>@param gradient
        !> procedure used to compute the gradient along the x-axis
        !
        !>@param dx
        !> grid space step along the x-axis
        !
        !>@return grad_var
        !> gradient of the governing variables along the x-axis
        !--------------------------------------------------------------
        function compute_x_gradient(nodes,i,j,gradient,dx) result(grad_var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(gradient_x_proc)                :: gradient
          real(rkind)                  , intent(in) :: dx
          real(rkind), dimension(ne)                :: grad_var


          grad_var(1) = gradient(nodes,i,j,position  ,dx)
          grad_var(2) = gradient(nodes,i,j,velocity_x,dx)
          grad_var(3) = gradient(nodes,i,j,velocity_y,dx)

        end function compute_x_gradient


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface for the computation of the gradient of the
        !> governing variables in the y-direction 
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> integer identifying the index in the x-direction
        !
        !>@param j
        !> integer identifying the index in the y-direction
        !
        !>@param gradient
        !> procedure used to compute the gradient along the y-axis
        !
        !>@param dy
        !> grid space step along the y-axis
        !
        !>@return grad_var
        !> gradient of the governing variables along the x-axis
        !--------------------------------------------------------------
        function compute_y_gradient(nodes,i,j,gradient,dy) result(grad_var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(gradient_y_proc)                :: gradient
          real(rkind)                  , intent(in) :: dy
          real(rkind), dimension(ne)                :: grad_var

          grad_var(1) = gradient(nodes,i,j,position  ,dy)
          grad_var(2) = gradient(nodes,i,j,velocity_x,dy)
          grad_var(3) = gradient(nodes,i,j,velocity_y,dy)

        end function compute_y_gradient        

      end module pmodel_eq_class

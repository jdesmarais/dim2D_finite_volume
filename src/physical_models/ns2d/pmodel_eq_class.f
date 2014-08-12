      !> @file
      !> class encapsulating subroutines to compute
      !> the reduced governing equations of the Navier-Stokes
      !> equations in 2D
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the reduced
      !> governing equations of the Navier-Stokes equations in 2D
      !
      !> @date
      !> 08_08_2014 - initial version               - J.L. Desmarais
      !-----------------------------------------------------------------
      module pmodel_eq_class

        use interface_primary           , only : gradient_x_proc,
     $                                           gradient_y_proc
        use sd_operators_class          , only : sd_operators
        use ns2d_parameters             , only : viscous_r, Re, Pr,
     $                                           gamma, mach_infty,
     $                                           gravity
        use ns2d_vortex_module          , only : apply_vortex_ic
        use ns2d_prim_module            , only : mass_density,
     $                                           momentum_x,
     $                                           momentum_y,
     $                                           total_energy,
     $                                           speed_of_sound
        use ns2d_fluxes_module          , only : flux_x_mass_density,
     $                                           flux_y_mass_density,
     $                                           flux_x_momentum_x,
     $                                           flux_y_momentum_x,
     $                                           flux_x_momentum_y,
     $                                           flux_y_momentum_y,
     $                                           flux_x_total_energy,
     $                                           flux_y_total_energy
c$$$        use ns2d_ncoords_module         , only : compute_n_gradient_ns2d,
c$$$     $                                           compute_n_eigenvalues_ns2d,
c$$$     $                                           compute_n1_lefteigenvector_ns2d,
c$$$     $                                           compute_n1_righteigenvector_ns2d,
c$$$     $                                           compute_n2_lefteigenvector_ns2d,
c$$$     $                                           compute_n2_righteigenvector_ns2d
        use ns2d_steadystate_module     , only : apply_steady_state_ic
        use parameters_bf_layer         , only : interior_pt
        use parameters_constant         , only : scalar,
     $                                           vector_x, vector_y,
     $                                           steady_state,
     $                                           vortex,
     $                                           earth_gravity_choice
        use parameters_input            , only : nx,ny,ne,bc_size,
     $                                           ic_choice,
     $                                           gravity_choice,
     $                                           sigma_P
        use parameters_kind             , only : ikind,rkind
        use pmodel_eq_default_class     , only : pmodel_eq_default


        implicit none

        private
        public :: pmodel_eq


        !> @class pmodel_eq
        !> class encapsulating operators to compute
        !> the governing equations of the Navier-Stokes
        !> equatiosn in 2D
        !
        !> @param get_model_name
        !> get the name of the physcial model
        !
        !> @param get_var_name
        !> get the name of the main variables
        !> (mass, momentum_x, momentum_y, total_energy)
        !
        !> @param get_var_longname
        !> get the description of the main variables for the
        !> governing equations of the physical model
        !
        !> @param get_var_unit
        !> get the units of the main variables
        !> (\f$ kg.m^{-3}, kg.m^{-2}.s^{-1},
        !> kg.m^{-2}.s^{-1}, J.kg.m^{-3}) \f$
        !
        !> @param get_var_types
        !> get the type of the main variables
        !> (scalar, vector_x, vector_y, scalar)
        !
        !> @param get_sim_parameters
        !> get the simulation parameters (ex: Re, Pr, gamma, ...)
        !
        !> @param get_eq_nb
        !> get the number of governing equations: 4
        !
        !> @param initialize
        !> initialize the main parameters of the physical model
        !> (viscosity ratio, Reynolds, Prandtl, Mach numbers and
        !> heat capacity ratio)
        !
        !> @param apply_ic
        !> initialize the main variables of the governing equations
        !> considering the user choices (steady state, vortex ...)
        !
        !> @param compute_flux_x
        !> compute the fluxes along the x-axis with fixed sized arrays
        !
        !> @param compute_flux_y
        !> compute the fluxes along the y-axis with fixed sized arrays    
        !
        !> @param compute_flux_x_nopt
        !> compute the fluxes along the x-axis with non-fixed sized arrays
        !
        !> @param compute_flux_y_nopt
        !> compute the fluxes along the y-axis with non-fixed sized arrays
        !
        !> @param compute_body_forces
        !> compute the body forces at a grid point location
        !
        !> @param get_velocity
        !> compute the velocity at a grid point location
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
          procedure, nopass :: get_sim_parameters
          procedure, nopass :: get_eq_nb
          procedure, nopass :: apply_ic
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
          procedure, nopass :: compute_x_gradient
          procedure, nopass :: compute_y_gradient

          procedure, nopass :: compute_lodi_subsonic_inflow_csttemp
          procedure, nopass :: compute_lodi_subsonic_outflow_cstpressure          

        end type pmodel_eq


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to get the name of the physical model
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
        !
        !>@param model_name
        !> character giving the name of the model
        !---------------------------------------------------------------
        function get_model_name() result(model_name)

          implicit none

          character(len=10) :: model_name

          model_name="NS2D"

        end function get_model_name
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the name of the main variables
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
        !
        !>@param var_name
        !> characters giving the variable names
        !---------------------------------------------------------------
        function get_var_name() result(var_pties)

          implicit none

          character(len=10), dimension(ne) :: var_pties

          var_pties(1)="mass"
          var_pties(2)="momentum_x"
          var_pties(3)="momentum_y"
          var_pties(4)="energy"          

        end function get_var_name
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the name of the main variables
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
        !
        !>@param var_name
        !> characters giving the variable names
        !---------------------------------------------------------------
        function get_var_longname() result(var_pties)

          implicit none

          character(len=33), dimension(ne) :: var_pties

          var_pties(1)="mass density"
          var_pties(2)="momentum density along the x-axis"
          var_pties(3)="momentum density along the y-axis"
          var_pties(4)="total energy density"

        end function get_var_longname


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the units of the main variables
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
        !
        !>@param var_name
        !> characters giving the variable units
        !---------------------------------------------------------------
        function get_var_unit() result(var_pties)

          implicit none

          character(len=23), dimension(ne) :: var_pties

          var_pties(1)= "(kg/m3)/(kg/m3)"
          var_pties(2)= "(kg/(m2.s))/(kg/(m2.s))"
          var_pties(3)= "(kg/(m2.s))/(kg/(m2.s))"
          var_pties(4)= "(J/m3)/(J.m3)"

        end function get_var_unit


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to get the type of the main variables
        !> (scalar, vector_x, vector_y, scalar)
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
          var_type(4)=scalar

        end function get_var_type


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the simulation parameters
        !
        !> @date
        !> 12_08_2014 - initial version - J.L. Desmarais
        !
        !>@param param_name
        !> array with the name of the characteristic parameters
        !> for the simulation
        !
        !>@param param_value
        !> array with the value of the characteristic parameters
        !> for the simulation
        !--------------------------------------------------------------
        subroutine get_sim_parameters(param_name, param_value)

          implicit none

          character(10), dimension(:), allocatable, intent(out) :: param_name
          real(rkind)  , dimension(:), allocatable, intent(out) :: param_value


          allocate(param_name(6))
          allocate(param_value(6))

          param_name(1) = 'viscous_r'
          param_name(2) = 'Re'
          param_name(3) = 'Pr'
          param_name(4) = 'gamma'
          param_name(5) = 'mach_infty'
          param_name(6) = 'gravity'

          param_value(1) = viscous_r
          param_value(2) = Re
          param_value(3) = Pr
          param_value(4) = gamma
          param_value(5) = mach_infty
          param_value(6) = gravity

        end subroutine get_sim_parameters
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the number of main variables
        !> in the governing equations: 4
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
        !
        !>@param eq_nb
        !> number of governing equations
        !---------------------------------------------------------------
        function get_eq_nb() result(eq_nb)
          implicit none
          integer :: eq_nb
          eq_nb=4
        end function get_eq_nb
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the initial conditions to the main
        !> variables of the governing equations
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !---------------------------------------------------------------
        subroutine apply_ic(nodes,x_map,y_map)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map


          real(rkind) :: x_s
          real(rkind) :: y_s
          

          !<initialize the field depending on the user choice
          select case(ic_choice)

            case(steady_state)
               call apply_steady_state_ic(nodes)

            case(vortex)
               call apply_vortex_ic(nodes,x_map,y_map)

            case default
               print '(''pmodel_eq_class'')'
               stop 'ic_choice not recognized'
          end select

          x_s = x_map(1)
          y_s = y_map(1)

        end subroutine apply_ic
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to apply the initial conditions
        !> to the main variables of the governing
        !> equations
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
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

          integer(ikind) :: i,j


          !<fluxes along the x-axis
          do j=1+bc_size, ny-bc_size
             !DEC$ IVDEP
             do i=1+bc_size, nx+1-bc_size

                !DEC$ FORCEINLINE RECURSIVE
                flux_x(i,j,1) =
     $               flux_x_mass_density(
     $               nodes,s,i,j)

                !DEC$ FORCEINLINE RECURSIVE
                flux_x(i,j,2) = flux_x_momentum_x(
     $               nodes,s,i,j,
     $               dx, dy)

                !DEC$ FORCEINLINE RECURSIVE
                flux_x(i,j,3) = flux_x_momentum_y(
     $               nodes,s,i,j,
     $               dx, dy)

                !DEC$ FORCEINLINE RECURSIVE
                flux_x(i,j,4) = flux_x_total_energy(
     $               nodes,s,i,j,
     $               dx, dy)

             end do
          end do

        end function compute_flux_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to apply the initial conditions
        !> to the main variables of the governing
        !> equations
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
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

          integer(ikind) :: i,j


          !<fluxes along the y-axis
          do j=1+bc_size, ny+1-bc_size
             !DEC$ IVDEP
             do i=1+bc_size, nx-bc_size

                !DEC$ FORCEINLINE RECURSIVE
                flux_y(i,j,1) = flux_y_mass_density(
     $               nodes,s,i,j)

                !DEC$ FORCEINLINE RECURSIVE
                flux_y(i,j,2) = flux_y_momentum_x(
     $               nodes,s,i,j,
     $               dx, dy)

                !DEC$ FORCEINLINE RECURSIVE
                flux_y(i,j,3) = flux_y_momentum_y(
     $               nodes,s,i,j,
     $               dx, dy)

                !DEC$ FORCEINLINE RECURSIVE
                flux_y(i,j,4) = flux_y_total_energy(
     $               nodes,s,i,j,
     $               dx, dy)

             end do
          end do

        end function compute_flux_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to apply the initial conditions
        !> to the main variables of the governing
        !> equations
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> physical model
        !>
        !>@param nodes
        !> array with the grid point data
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param s
        !> space discretization operators      
        !
        !>@param grdpts_id
        !> role of the grid points
        !
        !>@param flux_x
        !> fluxes along the x-axis
        !---------------------------------------------------------------
        subroutine compute_flux_x_nopt(nodes,dx,dy,s,grdpts_id,flux_x)
        
          implicit none

          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          real(rkind)                  , intent(in)    :: dx
          real(rkind)                  , intent(in)    :: dy
          type(sd_operators)           , intent(in)    :: s
          integer    , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind), dimension(:,:,:), intent(inout) :: flux_x

          integer(ikind) :: i,j


          !<fluxes along the x-axis
          do j=1+bc_size, size(flux_x,2)-bc_size
             !DEC$ IVDEP
             do i=1+bc_size, size(flux_x,1)-bc_size

                if(grdpts_id(i,j).eq.interior_pt) then

                   !DEC$ FORCEINLINE RECURSIVE
                   flux_x(i,j,1) =
     $                  flux_x_mass_density(
     $                  nodes,s,i,j)
                   
                   !DEC$ FORCEINLINE RECURSIVE
                   flux_x(i,j,2) = flux_x_momentum_x(
     $                  nodes,s,i,j,
     $                  dx, dy)
                   
                   !DEC$ FORCEINLINE RECURSIVE
                   flux_x(i,j,3) = flux_x_momentum_y(
     $                  nodes,s,i,j,
     $                  dx, dy)
                   
                   !DEC$ FORCEINLINE RECURSIVE
                   flux_x(i,j,4) = flux_x_total_energy(
     $                  nodes,s,i,j,
     $                  dx, dy)

                end if

             end do
          end do

        end subroutine compute_flux_x_nopt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to apply the initial conditions
        !> to the main variables of the governing
        !> equations
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> physical model
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> space discretization operators
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param grdpts_id
        !> role of the grid points
        !
        !>@param flux_y
        !> fluxes along the y-axis
        !---------------------------------------------------------------
        subroutine compute_flux_y_nopt(nodes,dx,dy,s,grdpts_id,flux_y)
        
          implicit none

          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          real(rkind)                  , intent(in)    :: dx
          real(rkind)                  , intent(in)    :: dy
          type(sd_operators)           , intent(in)    :: s
          integer    , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind), dimension(:,:,:), intent(inout) :: flux_y

          integer(ikind) :: i,j


          !<fluxes along the y-axis
          do j=1+bc_size, size(flux_y,2)-bc_size
             !DEC$ IVDEP
             do i=1+bc_size, size(flux_y,1)-bc_size

                if(grdpts_id(i,j).eq.interior_pt) then

                   !DEC$ FORCEINLINE RECURSIVE
                   flux_y(i,j,1) = flux_y_mass_density(nodes,s,i,j)
                   
                   !DEC$ FORCEINLINE RECURSIVE
                   flux_y(i,j,2) = flux_y_momentum_x(
     $                  nodes,s,i,j,
     $                  dx,dy)
                   
                   !DEC$ FORCEINLINE RECURSIVE
                   flux_y(i,j,3) = flux_y_momentum_y(
     $                  nodes,s,i,j,
     $                  dx,dy)
                   
                   !DEC$ FORCEINLINE RECURSIVE
                   flux_y(i,j,4) = flux_y_total_energy(
     $                  nodes,s,i,j,
     $                  dx,dy)

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

          real(rkind), dimension(nx,ny,ne), intent(in)   :: nodes
          real(rkind)                     , intent(in)   :: dx
          real(rkind)                     , intent(in)   :: dy
          integer(ikind)                  , intent(in)   :: i
          integer(ikind)                  , intent(in)   :: j
          class(sd_operators)             , intent(in)   :: s_oneside
          real(rkind), dimension(ne)                     :: flux_x

          !<fluxes along the x-axis
          !DEC$ FORCEINLINE RECURSIVE
          flux_x(1) = flux_y_mass_density(
     $         nodes,s_oneside,i,j)
          
          !DEC$ FORCEINLINE RECURSIVE
          flux_x(2) = flux_y_momentum_x(
     $         nodes,s_oneside,i,j,
     $         dx,dy)
          
          !DEC$ FORCEINLINE RECURSIVE
          flux_x(3) = flux_y_momentum_y(
     $         nodes,s_oneside,i,j,
     $         dx,dy)
          
          !DEC$ FORCEINLINE RECURSIVE
          flux_x(4) = flux_x_total_energy(
     $         nodes,s_oneside,i,j,
     $         dx,dy)

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

          real(rkind), dimension(nx,ny,ne), intent(in)   :: nodes
          real(rkind)                     , intent(in)   :: dx
          real(rkind)                     , intent(in)   :: dy
          integer(ikind)                  , intent(in)   :: i
          integer(ikind)                  , intent(in)   :: j
          class(sd_operators)             , intent(in)   :: s_oneside
          real(rkind), dimension(ne)                     :: flux_y


          !DEC$ FORCEINLINE RECURSIVE
          flux_y(1) = flux_y_mass_density(
     $         nodes,s_oneside,i,j)
          
          !DEC$ FORCEINLINE RECURSIVE
          flux_y(2) = flux_y_momentum_x(
     $         nodes,s_oneside,i,j,
     $         dx,dy)
          
          !DEC$ FORCEINLINE RECURSIVE
          flux_y(3) = flux_y_momentum_y(
     $         nodes,s_oneside,i,j,
     $         dx,dy)
          
          !DEC$ FORCEINLINE RECURSIVE
          flux_y(4) = flux_y_total_energy(
     $         nodes,s_oneside,i,j,
     $         dx,dy)

        end function compute_flux_y_oneside


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to compute the body forces
        !> acting on the cell
        !
        !> @date
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !
        !>@param body_forces
        !> body forces evaluated at (i,j)
        !--------------------------------------------------------------
        function compute_body_forces(nodes,k) result(body_forces)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          integer                   , intent(in) :: k
          real(rkind)                            :: body_forces

          select case(k)
            case(1)
               body_forces = 0
            case(2)
               body_forces = 0
            case(3)
               body_forces = -nodes(1)*gravity
            case(4)
               body_forces = -nodes(3)*gravity
            case default
               body_forces = 0
               print '(''ns2d/pmodel_eq_class'')'
               print '(''compute_body_forces'')'
               stop '1=< k =<4 violated'
          end select

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

          velocity(1) = nodes(2)/nodes(1)
          velocity(2) = nodes(3)/nodes(1)

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
        !>@return undermined
        !> check if the open boundary conditions are undermined
        !> at the grid point location
        !--------------------------------------------------------------
        function are_openbc_undermined(nodes) result(undermined)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          logical                                :: undermined

          real(rkind) :: node_s

          undermined = .false.

          node_s=nodes(1)

        end function are_openbc_undermined


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the Jacobian matrix for primitive to
        !> to conservative variables
        !
        !> @date
        !> 11_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return jacPrimCons
        !> jacobian matrix for conservative to primitive
        !> variables \f$ \frac{\partial p}{\partial v} \f$
        !--------------------------------------------------------------
        function compute_jacobian_prim_to_cons(nodes)
     $     result(jacPrimCons)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: jacPrimCons

          real(rkind)                   :: ux
          real(rkind)                   :: uy

          ux = nodes(2)/nodes(1)
          uy = nodes(3)/nodes(1)


          if(rkind.eq.8) then

             jacPrimCons(1,1) = 1.0d0
             jacPrimCons(2,1) = 0.0d0
             jacPrimCons(3,1) = 0.0d0
             jacPrimCons(4,1) = 0.0d0

             jacPrimCons(1,2) = - ux/nodes(1)
             jacPrimCons(2,2) = 1.0d0/nodes(1)
             jacPrimCons(3,2) = 0.0d0
             jacPrimCons(4,2) = 0.0d0

             jacPrimCons(1,3) = - uy/nodes(1)
             jacPrimCons(2,3) = 0.0d0
             jacPrimCons(3,3) = 1.0d0/nodes(1)
             jacPrimCons(4,3) = 0.0d0

             jacPrimCons(1,4) = 0.5d0*(gamma-1.0d0)*(ux**2+uy**2)
             jacPrimCons(2,4) = -(gamma-1.0d0)*ux
             jacPrimCons(3,4) = -(gamma-1.0d0)*uy
             jacPrimCons(4,4) = gamma-1.0d0             

          else

             jacPrimCons(1,1) = 1.0d0
             jacPrimCons(2,1) = 0.0d0
             jacPrimCons(3,1) = 0.0d0
             jacPrimCons(4,1) = 0.0d0

             jacPrimCons(1,2) = - ux/nodes(1)
             jacPrimCons(2,2) = 1.0d0/nodes(1)
             jacPrimCons(3,2) = 0.0d0
             jacPrimCons(4,2) = 0.0d0

             jacPrimCons(1,3) = - uy/nodes(1)
             jacPrimCons(2,3) = 0.0d0
             jacPrimCons(3,3) = 1.0d0/nodes(1)
             jacPrimCons(4,3) = 0.0d0

             jacPrimCons(1,4) = 0.5d0*(gamma-1.0d0)*(ux**2+uy**2)
             jacPrimCons(2,4) = -(gamma-1.0d0)*ux
             jacPrimCons(3,4) = -(gamma-1.0d0)*uy
             jacPrimCons(4,4) = gamma-1.0d0           

          end if

        end function compute_jacobian_prim_to_cons


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the Jacobian matrix for conservative
        !> to primitive variables
        !
        !> @date
        !> 11_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return jacConsPrim
        !> jacobian matrix for primitive to conservative
        !> variables \f$ \frac{\partial v}{\partial p} \f$
        !--------------------------------------------------------------
        function compute_jacobian_cons_to_prim(nodes)
     $     result(jacConsPrim)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: jacConsPrim

          real(rkind) :: ux
          real(rkind) :: uy

          ux = nodes(2)/nodes(1)
          uy = nodes(3)/nodes(1)          


          if(rkind.eq.8) then
             jacConsPrim(1,1) = 1.0d0
             jacConsPrim(2,1) = 0.0d0
             jacConsPrim(3,1) = 0.0d0
             jacConsPrim(4,1) = 0.0d0

             jacConsPrim(1,2) = ux
             jacConsPrim(2,2) = nodes(1)
             jacConsPrim(3,2) = 0.0d0
             jacConsPrim(4,2) = 0.0d0

             jacConsPrim(1,3) = uy
             jacConsPrim(2,3) = 0.0d0
             jacConsPrim(3,3) = nodes(1)
             jacConsPrim(4,3) = 0.0d0

             jacConsPrim(1,4) = 0.5d0*(ux**2+uy**2)
             jacConsPrim(2,4) = nodes(1)*ux
             jacConsPrim(3,4) = nodes(1)*uy
             jacConsPrim(4,4) = 1.0d0/(gamma-1.0d0)

          else
             jacConsPrim(1,1) = 1.0
             jacConsPrim(2,1) = 0.0
             jacConsPrim(3,1) = 0.0
             jacConsPrim(4,1) = 0.0

             jacConsPrim(1,2) = ux
             jacConsPrim(2,2) = nodes(1)
             jacConsPrim(3,2) = 0.0
             jacConsPrim(4,2) = 0.0

             jacConsPrim(1,3) = uy
             jacConsPrim(2,3) = 0.0
             jacConsPrim(3,3) = nodes(1)
             jacConsPrim(4,3) = 0.0

             jacConsPrim(1,4) = 0.5*(ux**2+uy**2)
             jacConsPrim(2,4) = nodes(1)*ux
             jacConsPrim(3,4) = nodes(1)*uy
             jacConsPrim(4,4) = 1.0/(gamma-1.0)

          end if             
             
        end function compute_jacobian_cons_to_prim


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the x-direction
        !
        !> @date
        !> 11_08_2014 - initial version - J.L. Desmarais
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

          real(rkind) :: velocity_x
          real(rkind) :: speedofsound
          
          velocity_x   = nodes(2)/nodes(1)
          speedofsound = speed_of_sound(nodes)

          eigenvalues(1) = velocity_x
          eigenvalues(2) = velocity_x
          eigenvalues(3) = velocity_x - speedofsound
          eigenvalues(4) = velocity_x + speedofsound

        end function compute_x_eigenvalues


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the y-direction
        !
        !> @date
        !> 11_08_2014 - initial version - J.L. Desmarais
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

          real(rkind) :: velocity_y
          real(rkind) :: speedofsound
          
          velocity_y   = nodes(3)/nodes(1)
          speedofsound = speed_of_sound(nodes)

          eigenvalues(1) = velocity_y
          eigenvalues(2) = velocity_y
          eigenvalues(3) = velocity_y - speedofsound
          eigenvalues(4) = velocity_y + speedofsound

        end function compute_y_eigenvalues


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenvectors for the hyperbolic terms
        !> in the x-direction
        !
        !> @date
        !> 11_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvectors at the location of the grid point
        !--------------------------------------------------------------
        function compute_x_lefteigenvector(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne)   , intent(in) :: nodes
          real(rkind), dimension(ne,ne)             :: eigenvect


          real(rkind), dimension(ne,ne) :: jacPrimCons
          real(rkind), dimension(ne,ne) :: leftEigenMPrim
          real(rkind)                   :: c


          !computation of J, the jacobian matrix from primitive
          !to conservative variables
          jacPrimCons = compute_jacobian_prim_to_cons(nodes)


          !left eigenmatrix for the primitive
          !variables, L_p
          c  = speed_of_sound(nodes)
          if(rkind.eq.8) then

             leftEigenMPrim(1,1) =  0.0d0
             leftEigenMPrim(2,1) =  0.0d0
             leftEigenMPrim(3,1) =  1.0d0
             leftEigenMPrim(4,1) =  0.0d0
             
             leftEigenMPrim(1,2) =  1.0d0
             leftEigenMPrim(2,2) =  0.0d0
             leftEigenMPrim(3,2) =  0.0d0
             leftEigenMPrim(4,2) = -1.d0/c**2

             leftEigenMPrim(1,3) =  0.0d0
             leftEigenMPrim(2,3) = -0.5d0*nodes(1)*c
             leftEigenMPrim(3,3) =  0.0d0
             leftEigenMPrim(4,3) =  0.5d0

             leftEigenMPrim(1,4) =  0.0d0
             leftEigenMPrim(2,4) =  0.5d0*nodes(1)*c
             leftEigenMPrim(3,4) =  0.0d0
             leftEigenMPrim(4,4) =  0.5d0

          else

             leftEigenMPrim(1,1) =  0.0
             leftEigenMPrim(2,1) =  0.0
             leftEigenMPrim(3,1) =  1.0
             leftEigenMPrim(4,1) =  0.0
             
             leftEigenMPrim(1,2) =  1.0
             leftEigenMPrim(2,2) =  0.0
             leftEigenMPrim(3,2) =  0.0
             leftEigenMPrim(4,2) = -1./c**2

             leftEigenMPrim(1,3) =  0.0
             leftEigenMPrim(2,3) = -0.5*nodes(1)*c
             leftEigenMPrim(3,3) =  0.0
             leftEigenMPrim(4,3) =  0.5

             leftEigenMPrim(1,4) =  0.0
             leftEigenMPrim(2,4) =  0.5*nodes(1)*c
             leftEigenMPrim(3,4) =  0.0
             leftEigenMPrim(4,4) =  0.5

          end if


          !compute the left eigenmatrix by L = L_p.J
          eigenvect = MATMUL(jacPrimCons,leftEigenMPrim)

        end function compute_x_lefteigenvector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the right eigenvectors for the hyperbolic terms
        !> in the x-direction
        !
        !> @date
        !> 11_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param k
        !> integer identifying the eigenvector
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_x_righteigenvector(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne)   , intent(in) :: nodes
          real(rkind), dimension(ne,ne)             :: eigenvect


          real(rkind), dimension(ne,ne) :: jacConsPrim
          real(rkind), dimension(ne,ne) :: rightEigenMPrim
          real(rkind)                   :: c


          !computation of J, the jacobian matrix from conservative
          !to primitive variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)


          !right eigenmatrix for the primitive
          !variables, R_p
          c  = speed_of_sound(nodes)

          if(rkind.eq.8) then
             rightEigenMPrim(1,1) =  0.0d0
             rightEigenMPrim(2,1) =  1.0d0
             rightEigenMPrim(3,1) =  1.0d0/c**2
             rightEigenMPrim(4,1) =  1.0d0/c**2

             rightEigenMPrim(1,2) =  0.0d0
             rightEigenMPrim(2,2) =  0.0d0
             rightEigenMPrim(3,2) = -1.0d0/(nodes(1)*c)
             rightEigenMPrim(4,2) =  1.0d0/(nodes(1)*c)             

             rightEigenMPrim(1,3) =  1.0d0
             rightEigenMPrim(2,3) =  0.0d0
             rightEigenMPrim(3,3) =  0.0d0
             rightEigenMPrim(4,3) =  0.0d0

             rightEigenMPrim(1,4) =  0.0d0
             rightEigenMPrim(2,4) =  0.0d0
             rightEigenMPrim(3,4) =  1.0d0
             rightEigenMPrim(4,4) =  1.0d0

          else
             rightEigenMPrim(1,1) =  0.0
             rightEigenMPrim(2,1) =  1.0
             rightEigenMPrim(3,1) =  1.0/c**2
             rightEigenMPrim(4,1) =  1.0/c**2

             rightEigenMPrim(1,2) =  0.0
             rightEigenMPrim(2,2) =  0.0
             rightEigenMPrim(3,2) = -1.0/(nodes(1)*c)
             rightEigenMPrim(4,2) =  1.0/(nodes(1)*c)             

             rightEigenMPrim(1,3) =  1.0
             rightEigenMPrim(2,3) =  0.0
             rightEigenMPrim(3,3) =  0.0
             rightEigenMPrim(4,3) =  0.0

             rightEigenMPrim(1,4) =  0.0
             rightEigenMPrim(2,4) =  0.0
             rightEigenMPrim(3,4) =  1.0
             rightEigenMPrim(4,4) =  1.0

          end if

          !right eigenmatrix computed as R = J.R_p
          eigenvect = MATMUL(rightEigenMPrim,jacConsPrim)

        end function compute_x_righteigenvector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenvector for the hyperbolic terms
        !> in the y-direction. By denoting L the left eigenmatrix, the
        !> result of the function is L[k,:]
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param k
        !> integer identifying the eigenvector
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_y_lefteigenvector(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect

          real(rkind), dimension(ne,ne) :: jacPrimCons
          real(rkind), dimension(ne,ne) :: leftEigenMPrim
          real(rkind)                   :: c


          !computation of J, the jacobian matrix from primitive
          !to conservative variables
          jacPrimCons = compute_jacobian_prim_to_cons(nodes)


          !left eigenmatrix for the primitive
          !variables, L_p
          c  = speed_of_sound(nodes)

          if(rkind.eq.8) then             
             leftEigenMPrim(1,1) =  0.0d0
             leftEigenMPrim(2,1) =  1.0d0
             leftEigenMPrim(3,1) =  0.0d0
             leftEigenMPrim(4,1) =  0.0d0

             leftEigenMPrim(1,2) =  1.0d0
             leftEigenMPrim(2,2) =  0.0d0
             leftEigenMPrim(3,2) =  0.0d0
             leftEigenMPrim(4,2) = -1.0d0/c**2

             leftEigenMPrim(1,3) =  0.0d0
             leftEigenMPrim(2,3) =  0.0d0
             leftEigenMPrim(3,3) = -0.5d0*nodes(1)*c
             leftEigenMPrim(4,3) =  0.5d0             

             leftEigenMPrim(1,4) =  0.0d0
             leftEigenMPrim(2,4) =  0.0d0
             leftEigenMPrim(3,4) =  0.5d0*nodes(1)*c
             leftEigenMPrim(4,4) =  0.5d0             

          else
             leftEigenMPrim(1,1) =  0.0
             leftEigenMPrim(2,1) =  1.0
             leftEigenMPrim(3,1) =  0.0
             leftEigenMPrim(4,1) =  0.0

             leftEigenMPrim(1,2) =  1.0
             leftEigenMPrim(2,2) =  0.0
             leftEigenMPrim(3,2) =  0.0
             leftEigenMPrim(4,2) = -1.0/c**2

             leftEigenMPrim(1,3) =  0.0
             leftEigenMPrim(2,3) =  0.0
             leftEigenMPrim(3,3) = -0.5*nodes(1)*c
             leftEigenMPrim(4,3) =  0.5             

             leftEigenMPrim(1,4) =  0.0
             leftEigenMPrim(2,4) =  0.0
             leftEigenMPrim(3,4) =  0.5*nodes(1)*c
             leftEigenMPrim(4,4) =  0.5             

          end if


          !compute the left eigenmatrix by L = L_p.J
          eigenvect = MATMUL(jacPrimCons,leftEigenMPrim)

        end function compute_y_lefteigenvector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenvector for the hyperbolic terms
        !> in the y-direction. By denoting R the right eigenmatrix, the
        !> result of the function is R[k,:]
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
        function compute_y_righteigenvector(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind), dimension(ne,ne) :: jacConsPrim
          real(rkind), dimension(ne,ne) :: rightEigenMPrim
          real(rkind)                   :: c


          !computation of J, the jacobian matrix from conservative
          !to primitive variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)


          !right eigenmatrix for the primitive
          !variables, R_p
          c  = speed_of_sound(nodes)

          if(rkind.eq.8) then
             rightEigenMPrim(1,1) =  0.0d0
             rightEigenMPrim(2,1) =  1.0d0
             rightEigenMPrim(3,1) =  1.0d0/c**2
             rightEigenMPrim(4,1) =  1.0d0/c**2

             rightEigenMPrim(1,2) =  1.0d0
             rightEigenMPrim(2,2) =  0.0d0
             rightEigenMPrim(3,2) =  0.0d0
             rightEigenMPrim(4,2) =  0.0d0

             rightEigenMPrim(1,3) =  0.0d0
             rightEigenMPrim(2,3) =  0.0d0
             rightEigenMPrim(3,3) = -1.0d0/(nodes(1)*c)
             rightEigenMPrim(4,3) =  1.0d0/(nodes(1)*c)

             rightEigenMPrim(1,4) =  0.0d0
             rightEigenMPrim(2,4) =  0.0d0
             rightEigenMPrim(3,4) =  1.0d0
             rightEigenMPrim(4,4) =  1.0d0

          else
             rightEigenMPrim(1,1) =  0.0
             rightEigenMPrim(2,1) =  1.0
             rightEigenMPrim(3,1) =  1.0/c**2
             rightEigenMPrim(4,1) =  1.0/c**2

             rightEigenMPrim(1,2) =  1.0
             rightEigenMPrim(2,2) =  0.0
             rightEigenMPrim(3,2) =  0.0
             rightEigenMPrim(4,2) =  0.0

             rightEigenMPrim(1,3) =  0.0
             rightEigenMPrim(2,3) =  0.0
             rightEigenMPrim(3,3) = -1.0/(nodes(1)*c)
             rightEigenMPrim(4,3) =  1.0/(nodes(1)*c)

             rightEigenMPrim(1,4) =  0.0
             rightEigenMPrim(2,4) =  0.0
             rightEigenMPrim(3,4) =  1.0
             rightEigenMPrim(4,4) =  1.0

          end if

          !right eigenmatrix computed as R = J.R_p
          eigenvect = MATMUL(rightEigenMPrim,jacConsPrim)

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

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          procedure(gradient_x_proc)                   :: gradient
          real(rkind)                     , intent(in) :: dx
          real(rkind), dimension(ne)                   :: grad_var

          grad_var(1) = gradient(nodes,i,j,mass_density,dx)
          grad_var(2) = gradient(nodes,i,j,momentum_x  ,dx)
          grad_var(3) = gradient(nodes,i,j,momentum_y  ,dx)
          grad_var(4) = gradient(nodes,i,j,total_energy,dx)

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

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          procedure(gradient_y_proc)                   :: gradient
          real(rkind)                     , intent(in) :: dy
          real(rkind), dimension(ne)                   :: grad_var

          grad_var(1) = gradient(nodes,i,j,mass_density,dy)
          grad_var(2) = gradient(nodes,i,j,momentum_x  ,dy)
          grad_var(3) = gradient(nodes,i,j,momentum_y  ,dy)
          grad_var(4) = gradient(nodes,i,j,total_energy,dy)

        end function compute_y_gradient


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface for the computation of the LODI amplitudes
        !> for the Poinsot and Lele open b.c.
        !
        !> @date
        !> 12_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables at the location of the open b.c.
        !
        !>@param mass_grad
        !> gradient of the mass density
        !
        !>@param velocity2_grad
        !> gradient of the velocity
        !
        !>@return lodi
        !> vector with the LODI amplitudes
        !--------------------------------------------------------------!
        function compute_lodi_subsonic_inflow_csttemp(
     $   nodes,
     $   eigenvalues,
     $   velocity1_grad, pressure_grad)
     $   result(lodi)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne), intent(in) :: eigenvalues
          real(rkind)               , intent(in) :: velocity1_grad
          real(rkind)               , intent(in) :: pressure_grad
          real(rkind), dimension(ne)             :: lodi

          real(rkind) :: c

          c = speed_of_sound(nodes)

          if(rkind.eq.8) then
             lodi(3) = eigenvalues(3)*(pressure_grad - nodes(1)*c*velocity1_grad)
             lodi(4) = lodi(3)
             lodi(2) = 0.5d0*(gamma-1.0d0)*(lodi(3)+lodi(4))
             lodi(1) = 0.0d0
          else
             lodi(3) = eigenvalues(3)*(pressure_grad - nodes(1)*c*velocity1_grad)
             lodi(4) = lodi(3)
             lodi(2) = 0.5*(gamma-1.0)*(lodi(3)+lodi(4))
             lodi(1) = 0.0
          end if

        end function compute_lodi_subsonic_inflow_csttemp


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface for the computation of the LODI amplitudes
        !> for the Poinsot and Lele open b.c.
        !
        !> @date
        !> 12_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables at the location of the open b.c.
        !
        !>@param mass_grad
        !> gradient of the mass density
        !
        !>@param velocity_grad
        !> gradient of the velocity
        !
        !>@return lodi
        !> vector with the LODI amplitudes
        !--------------------------------------------------------------!
        function compute_lodi_subsonic_outflow_cstpressure(
     $   nodes, eigenvalues,
     $   mass_grad, velocity1_grad, velocity2_grad, pressure_grad)
     $   result(lodi)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne), intent(in) :: eigenvalues
          real(rkind)               , intent(in) :: mass_grad
          real(rkind)               , intent(in) :: velocity1_grad
          real(rkind)               , intent(in) :: velocity2_grad
          real(rkind)               , intent(in) :: pressure_grad
          real(rkind), dimension(ne)             :: lodi

          real(rkind) :: c,p_infty,p,L,M,K

          c = speed_of_sound(nodes)
         
          if(rkind.eq.8) then

             p_infty = 1.0d0
             p       = (gamma-1.0d0)*(nodes(4)-0.5d0/nodes(1)*(
     $            nodes(2)**2+nodes(3)**2))

             L       = 2.0d0
             M       = 0.002d0
             K       = sigma_P*(1.0d0-M**2)*1.0d0/L

          else

             p_infty = 1.0
             p       = (gamma-1.0)*(nodes(4)-0.5/nodes(1)*(
     $            nodes(2)**2+nodes(3)**2))

             L       = 2.0
             M       = 0.002
             K       = sigma_P*(1.0-M**2)*1.0/L

          end if
             
          lodi(1) = eigenvalues(1)*velocity2_grad
          lodi(2) = eigenvalues(2)*(c**2*mass_grad-pressure_grad)
          lodi(3) = K*(p-p_infty)
          lodi(4) = eigenvalues(4)*(pressure_grad+nodes(1)*c*velocity1_grad)
             
        end function compute_lodi_subsonic_outflow_cstpressure

      end module pmodel_eq_class

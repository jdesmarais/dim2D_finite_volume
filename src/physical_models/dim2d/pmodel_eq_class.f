      !> @file
      !> class encapsulating subroutines to compute
      !> the reduced governing equations of the Diffuse
      !> Interface Model in 2D as derivated by J.Desmarais
      !> and J.G.M Kuerten
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the governing equations of the Diffuse
      !> Interface Model in 2D:
      !> “Extension of the 1-D characteristic open boundary
      !> conditions to the diffuse interface model”,
      !> Computational Methods in Multiphase Flow VII,
      !> WIT Press, J. Desmarais and J.G.M. Kuerten
      !
      !> @date
      !> 08_08_2013 - initial version               - J.L. Desmarais
      !> 11_07_2014 - interface for erymanthianboar - J.L. Desmarais
      !> 10_12_2014 - eigensystems                  - J.L. Desmarais
      !-----------------------------------------------------------------
      module pmodel_eq_class

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc,
     $       gradient_n_proc

        use sd_operators_class, only :
     $       sd_operators

        use dim2d_parameters, only :
     $       viscous_r, Re, We, Pr,
     $       cv_r, gravity

        use dim2d_bubble_ascending_module, only :
     $       apply_bubble_ascending_ic

        use dim2d_drop_collision_module , only :
     $       apply_drop_collision_ic 

        use dim2d_phase_separation_module, only :
     $       apply_phase_separation_ic

        !use dim2d_drop_evaporation_module, only:apply_drop_evaporation_ic
        use dim2d_drop_retraction_module, only :
     $       apply_drop_retraction_ic

        use dim2d_prim_module, only :
     $       mass_density,
     $       momentum_x,
     $       momentum_y,
     $       total_energy,
     $       speed_of_sound,
     $       compute_jacobian_cons_to_prim,
     $       compute_jacobian_prim_to_cons

        use dim2d_ncoords_module, only :
     $       compute_n1_eigenvalues_dim2d,
     $       compute_n2_eigenvalues_dim2d,
     $       compute_n1_lefteigenvector_dim2d,
     $       compute_n1_righteigenvector_dim2d,
     $       compute_n2_lefteigenvector_dim2d,
     $       compute_n2_righteigenvector_dim2d,
     $       compute_n1_transM_dim2d,
     $       compute_n2_transM_dim2d

        use dim2d_fluxes_module, only :
     $       flux_x_mass_density,
     $       flux_y_mass_density,
     $       flux_x_momentum_x,
     $       flux_y_momentum_x,
     $       flux_x_momentum_y,
     $       flux_y_momentum_y,
     $       flux_x_total_energy,
     $       flux_y_total_energy

        use dim2d_homogeneous_module, only :
     $       apply_homogeneous_ic

        use dim2d_steadystate_module, only :
     $       apply_steady_state_ic

        use parameters_bf_layer, only :
     $       bc_interior_pt,
     $       interior_pt

        use parameters_constant, only :
     $       scalar,
     $       vector_x, vector_y,
     $       steady_state,
     $       drop_retraction,
     $       bubble_ascending,
     $       homogeneous_liquid,
     $       drop_collision,
     $       phase_separation,
     $       earth_gravity_choice

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       ic_choice,
     $       gravity_choice

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_default_class, only :
     $       pmodel_eq_default


        implicit none

        private
        public :: pmodel_eq


        !> @class pmodel_eq
        !> class encapsulating operators to compute
        !> the governing equations of the Diffuse Interface
        !> Model in 2D
        !>
        !> @param get_model_name
        !> get the name of the physcial model
        !>
        !> @param get_var_name
        !> get the name of the main variables
        !> (mass, momentum_x, momentum_y, total_energy)
        !>
        !> @param get_var_longname
        !> get the description of the main variables for the
        !> governing equations of the physical model
        !>
        !> @param get_var_unit
        !> get the units of the main variables
        !> (\f$ kg.m^{-3}, kg.m^{-2}.s^{-1},
        !> kg.m^{-2}.s^{-1}, J.kg.m^{-3}) \f$
        !>
        !> @param get_var_types
        !> get the type of the main variables
        !> (scalar, vector_x, vector_y, scalar)
        !>
        !> @param get_eq_nb
        !> get the number of governing equations: 4
        !>
        !> @param initialize
        !> initialize the main parameters of the physical model
        !> (viscosity ratio, Reynolds, Prandtl, Weber numbers and
        !> reduced heat capacity)
        !>
        !> @param apply_ic
        !> initialize the main variables of the governing equations
        !> considering the user choices (drop retraction, two drops
        !> collision...)
        !>
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

          !description of the model and its main
          !governing variables
          procedure, nopass :: get_model_name
          procedure, nopass :: get_var_name
          procedure, nopass :: get_var_longname
          procedure, nopass :: get_var_unit
          procedure, nopass :: get_var_type
          procedure, nopass :: get_eq_nb
          procedure, nopass :: get_sim_parameters

          !initial conditions procedures
          procedure,   pass :: apply_ic

          !flux computation
          procedure, nopass :: compute_flux_x
          procedure, nopass :: compute_flux_y
          procedure, nopass :: compute_flux_x_nopt
          procedure, nopass :: compute_flux_y_nopt
          procedure, nopass :: compute_flux_x_oneside
          procedure, nopass :: compute_flux_y_oneside
          procedure, nopass :: compute_body_forces

          !field extension for open b.c.
          procedure, nopass :: get_velocity
          procedure, nopass :: are_openbc_undermined

          !eigenquantities computation
          procedure, nopass :: compute_x_eigenvalues
          procedure, nopass :: compute_y_eigenvalues
          procedure, nopass :: compute_x_lefteigenvector
          procedure, nopass :: compute_x_righteigenvector
          procedure, nopass :: compute_y_lefteigenvector
          procedure, nopass :: compute_y_righteigenvector

          procedure, nopass :: compute_n1_eigenvalues => compute_n1_eigenvalues_dim2d
          procedure, nopass :: compute_n2_eigenvalues => compute_n2_eigenvalues_dim2d
          procedure, nopass :: compute_n1_lefteigenvector  => compute_n1_lefteigenvector_dim2d
          procedure, nopass :: compute_n1_righteigenvector => compute_n1_righteigenvector_dim2d
          procedure, nopass :: compute_n2_lefteigenvector  => compute_n2_lefteigenvector_dim2d
          procedure, nopass :: compute_n2_righteigenvector => compute_n2_righteigenvector_dim2d

          !transverse matrices
          procedure, nopass :: compute_x_transM
          procedure, nopass :: compute_y_transM
          procedure, nopass :: compute_n1_transM => compute_n1_transM_dim2d
          procedure, nopass :: compute_n2_transM => compute_n2_transM_dim2d

          !conservative lodi matrices
          procedure, nopass :: compute_x_leftConsLodiM
          procedure, nopass :: compute_y_leftConsLodiM
          
          !gradient computation
          procedure, nopass :: compute_x_gradient
          procedure, nopass :: compute_y_gradient
          procedure, nopass :: compute_n_gradient
          
        end type pmodel_eq


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to get the name of the physical model
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param model_name
        !> character giving the name of the model
        !---------------------------------------------------------------
        function get_model_name() result(model_name)

          implicit none

          character(len=10) :: model_name

          model_name="DIM2D"

        end function get_model_name
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the name of the main variables
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        !> 08_08_2013 - initial version - J.L. Desmarais
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
          param_name(3) = 'We'
          param_name(4) = 'Pr'
          param_name(5) = 'cv_r'
          param_name(6) = 'gravity'

          param_value(1) = viscous_r
          param_value(2) = Re
          param_value(3) = We
          param_value(4) = Pr
          param_value(5) = cv_r
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
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !---------------------------------------------------------------
        subroutine apply_ic(this,nodes,x_map,y_map)

          implicit none

          class(pmodel_eq)             , intent(in)    :: this
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map

          integer :: neq

          neq = this%get_eq_nb()


          !<initialize the field depending on the user choice
          select case(ic_choice)

            case(steady_state)
               call apply_steady_state_ic(nodes)

            case(drop_retraction)
               call apply_drop_retraction_ic(nodes,x_map,y_map)

            case(bubble_ascending)
               call apply_bubble_ascending_ic(nodes,x_map,y_map)

            case(homogeneous_liquid)
               call apply_homogeneous_ic(nodes)

            case(drop_collision)
               call apply_drop_collision_ic(nodes,x_map,y_map)

            case(phase_separation)
               call apply_phase_separation_ic(nodes,x_map,y_map)

c$$$            case(drop_evaporation)
c$$$               call apply_drop_evaporation_ic(field_used)
            case default
               print '(''pmodel_eq_class'')'
               stop 'ic_choice not recognized'
          end select

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
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        !> 08_08_2013 - initial version - J.L. Desmarais
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


          !<fluxes along the x-axis
          do j=y_borders(1), y_borders(2)
             !DEC$ IVDEP
             do i=x_borders(1), x_borders(2)+1

                if((grdpts_id(i,j).eq.interior_pt).or.
     $               (grdpts_id(i,j).eq.bc_interior_pt))then

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
        !> 08_08_2013 - initial version - J.L. Desmarais
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


          !<fluxes along the y-axis
          do j=y_borders(1), y_borders(2)+1
             !DEC$ IVDEP
             do i=x_borders(1), x_borders(2)

                if((grdpts_id(i,j).eq.interior_pt).or.
     $               (grdpts_id(i,j).eq.bc_interior_pt))then

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

          real(rkind), dimension(:,:,:), intent(in)   :: nodes
          real(rkind)                  , intent(in)   :: dx
          real(rkind)                  , intent(in)   :: dy
          integer(ikind)               , intent(in)   :: i
          integer(ikind)               , intent(in)   :: j
          class(sd_operators)          , intent(in)   :: s_oneside
          real(rkind), dimension(ne)                  :: flux_x

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

          real(rkind), dimension(:,:,:), intent(in)   :: nodes
          real(rkind)                  , intent(in)   :: dx
          real(rkind)                  , intent(in)   :: dy
          integer(ikind)               , intent(in)   :: i
          integer(ikind)               , intent(in)   :: j
          class(sd_operators)          , intent(in)   :: s_oneside
          real(rkind), dimension(ne)                  :: flux_y


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
        function compute_body_forces(t,x,y,nodes,k) result(body_forces)

          implicit none

          real(rkind)               , intent(in) :: t
          real(rkind)               , intent(in) :: x
          real(rkind)               , intent(in) :: y
          real(rkind), dimension(ne), intent(in) :: nodes
          integer                   , intent(in) :: k
          real(rkind)                            :: body_forces

          real(rkind) :: t_s,x_s,y_s

          if(gravity_choice.eq.earth_gravity_choice) then
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
                  print '(''dim2d/pmodel_eq_class'')'
                  print '(''compute_body_forces'')'
                  stop '1=< k =<4 violated'
             end select

          else

             body_forces = 0

          end if

          t_s = t
          x_s = x
          y_s = y

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
        !>@param undermined
        !> check if the open boundary conditions are undermined
        !> at the grid point location
        !--------------------------------------------------------------
        function are_openbc_undermined(nodes) result(undermined)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          logical                                :: undermined

          real(rkind) :: d_liq, d_vap

          d_liq = 1.1-0.05*(1.1-0.1)
          d_vap = 0.1+0.05*(1.1-0.1)

          if((nodes(1).ge.d_vap).and.(nodes(1).le.d_liq)) then
             undermined = .true.
          else
             undermined = .false.
          end if

          !undermined = .true.

        end function are_openbc_undermined


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the x-direction
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
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

          real(rkind) :: ux
          real(rkind) :: c

          ux = nodes(2)/nodes(1)
          c  = speed_of_sound(nodes)

          eigenvalues(1) = ux
          eigenvalues(2) = ux
          eigenvalues(3) = ux-c
          eigenvalues(4) = ux+c

        end function compute_x_eigenvalues


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the y-direction
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
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

          real(rkind) :: uy
          real(rkind) :: c

          uy = nodes(3)/nodes(1)
          c  = speed_of_sound(nodes)

          eigenvalues(1) = uy
          eigenvalues(2) = uy
          eigenvalues(3) = uy-c
          eigenvalues(4) = uy+c

        end function compute_y_eigenvalues


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenmatrix for the hyperbolic
        !> terms in the x-direction
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param k
        !> integer identifying the eigenvector
        !
        !>@return eigenvect
        !> left eigenmatrix
        !--------------------------------------------------------------
        function compute_x_lefteigenvector(nodes) result(eigenvect)

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
             leftEigenMPrim(4,2) = -1.0/c**2

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
        !> computation of the left eigenmatrixr for the hyperbolic
        !> terms in the x-direction
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> right eigenmatrix
        !--------------------------------------------------------------
        function compute_x_righteigenvector(nodes) result(eigenvect)

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
        !> computation of the left eigenmatrix for the hyperbolic
        !> terms in the y-direction
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> left eigenmatrix
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
        !> computation of the right eigenmatrix for the hyperbolic
        !> terms in the y-direction
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> right eigenmatrix
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
        !> computation of the transverse matrix in the x-direction
        !> if the convective part of the governing equations is
        !> written as
        !> \f$ \frac{\partial v}{\partial t} +
        !>    A^x_v \frac{\partial v}{\partial x} + 
        !>    A^y_v \frac{\partial v}{\partial y}
        !> \f$
        !> then the transverse matrix in the x-direction is
        !> \f$ A^y_p = J^v_p \cdot A^y_p \cdot J^p_v\f$
        !> where $J^v_p$ and $J^p_v$ are the Jacobian matrices
        !> and $A^y_p$ is the transverse matrix for the primitive
        !> variables
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> transverse matrix in the x-direction
        !--------------------------------------------------------------
        function compute_x_transM(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind), dimension(ne,ne) :: jacConsPrim
          real(rkind), dimension(ne,ne) :: jacPrimCons
          real(rkind), dimension(ne,ne) :: xTransMPrim

          real(rkind)                   :: uy
          real(rkind)                   :: c


          !computation of J, the jacobian matrix from conservative
          !to primitive variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)
          jacPrimCons = compute_jacobian_prim_to_cons(nodes)


          !transverse matrix for the primitive variables
          uy = nodes(3)/nodes(1)         !velocity_y
          c  = speed_of_sound(nodes)     !speed of sound
          

          if(rkind.eq.8) then

             xTransMPrim(1,1) =  uy
             xTransMPrim(2,1) =  0.0d0
             xTransMPrim(3,1) =  nodes(1)
             xTransMPrim(4,1) =  0.0d0

             xTransMPrim(1,2) =  0.0d0
             xTransMPrim(2,2) =  uy
             xTransMPrim(3,2) =  0.0d0
             xTransMPrim(4,2) =  0.0d0

             xTransMPrim(1,3) =  0.0d0
             xTransMPrim(2,3) =  0.0d0
             xTransMPrim(3,3) =  uy
             xTransMPrim(4,3) =  1.0d0/nodes(1)

             xTransMPrim(1,4) =  0.0d0
             xTransMPrim(2,4) =  0.0d0
             xTransMPrim(3,4) =  c**2*nodes(1)
             xTransMPrim(4,4) =  uy

          else
             
             xTransMPrim(1,1) =  uy
             xTransMPrim(2,1) =  0.0
             xTransMPrim(3,1) =  nodes(1)
             xTransMPrim(4,1) =  0.0

             xTransMPrim(1,2) =  0.0
             xTransMPrim(2,2) =  uy
             xTransMPrim(3,2) =  0.0
             xTransMPrim(4,2) =  0.0

             xTransMPrim(1,3) =  0.0
             xTransMPrim(2,3) =  0.0
             xTransMPrim(3,3) =  uy
             xTransMPrim(4,3) =  1.0/nodes(1)

             xTransMPrim(1,4) =  0.0
             xTransMPrim(2,4) =  0.0
             xTransMPrim(3,4) =  c**2*nodes(1)
             xTransMPrim(4,4) =  uy

          end if

          !transverse matrix computed as A^y_v = J^v_p.A^y_p.J^p_v
          eigenvect = MATMUL(MATMUL(jacPrimCons,xTransMPrim),jacConsPrim)

        end function compute_x_transM


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the transverse matrix in the x-direction
        !> if the convective part of the governing equations is
        !> written as
        !> \f$ \frac{\partial v}{\partial t} +
        !>    A_x \frac{\partial v}{\partial x} + 
        !>    A_y \frac{\partial v}{\partial y}
        !> \f$
        !> then the transverse matrix in the y-direction is
        !> \f$ A^x_p = J^v_p \cdot A^x_p \cdot J^p_v\f$
        !> where $J^v_p$ and $J^p_v$ are the Jacobian matrices
        !> and $A^x_p$ is the transverse matrix for the primitive
        !> variables
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> transverse matrix in the y-direction
        !--------------------------------------------------------------
        function compute_y_transM(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind), dimension(ne,ne) :: jacConsPrim
          real(rkind), dimension(ne,ne) :: jacPrimCons
          real(rkind), dimension(ne,ne) :: yTransMPrim

          real(rkind)                   :: ux
          real(rkind)                   :: c


          !computation of J, the jacobian matrix from conservative
          !to primitive variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)
          jacPrimCons = compute_jacobian_prim_to_cons(nodes)


          !transverse matrix for the primitive variables
          ux = nodes(2)/nodes(1)         !velocity_x
          c  = speed_of_sound(nodes)     !speed of sound
          

          if(rkind.eq.8) then

             yTransMPrim(1,1) =  ux
             yTransMPrim(2,1) =  nodes(1)
             yTransMPrim(3,1) =  0.0d0
             yTransMPrim(4,1) =  0.0d0

             yTransMPrim(1,2) =  0.0d0
             yTransMPrim(2,2) =  ux
             yTransMPrim(3,2) =  0.0d0
             yTransMPrim(4,2) =  1.0d0/nodes(1)

             yTransMPrim(1,3) =  0.0d0
             yTransMPrim(2,3) =  0.0d0
             yTransMPrim(3,3) =  ux
             yTransMPrim(4,3) =  0.0d0

             yTransMPrim(1,4) =  0.0d0
             yTransMPrim(2,4) =  c**2*nodes(1)
             yTransMPrim(3,4) =  0.0d0
             yTransMPrim(4,4) =  ux

          else
             
             yTransMPrim(1,1) =  ux
             yTransMPrim(2,1) =  nodes(1)
             yTransMPrim(3,1) =  0.0
             yTransMPrim(4,1) =  0.0

             yTransMPrim(1,2) =  0.0
             yTransMPrim(2,2) =  ux
             yTransMPrim(3,2) =  0.0
             yTransMPrim(4,2) =  1.0/nodes(1)

             yTransMPrim(1,3) =  0.0
             yTransMPrim(2,3) =  0.0
             yTransMPrim(3,3) =  ux
             yTransMPrim(4,3) =  0.0

             yTransMPrim(1,4) =  0.0
             yTransMPrim(2,4) =  c**2*nodes(1)
             yTransMPrim(3,4) =  0.0
             yTransMPrim(4,4) =  ux

          end if

          !transverse matrix computed as A^x_v = J^v_p.A^x_p.J^p_v
          eigenvect = MATMUL(MATMUL(jacPrimCons,yTransMPrim),jacConsPrim)

        end function compute_y_transM


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

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(gradient_y_proc)                :: gradient
          real(rkind)                  , intent(in) :: dy
          real(rkind), dimension(ne)                :: grad_var

          grad_var(1) = gradient(nodes,i,j,mass_density,dy)
          grad_var(2) = gradient(nodes,i,j,momentum_x  ,dy)
          grad_var(3) = gradient(nodes,i,j,momentum_y  ,dy)
          grad_var(4) = gradient(nodes,i,j,total_energy,dy)

        end function compute_y_gradient


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface for the computation of the gradient of the
        !> governing variables in the n-direction (either (x-y)
        !> or (x+y) direction depending on the procedure passed)
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
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
        !>@param dx
        !> grid space step along the x-axis
        !
        !>@param dy
        !> grid space step along the y-axis
        !
        !>@return grad_var
        !> gradient of the governing variables along the x-axis
        !--------------------------------------------------------------
        function compute_n_gradient(nodes,i,j,gradient,dx,dy) result(grad_var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(gradient_n_proc)                :: gradient
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind), dimension(ne)                :: grad_var

          grad_var(1) = gradient(nodes,i,j,mass_density,dx,dy)
          grad_var(2) = gradient(nodes,i,j,momentum_x  ,dx,dy)
          grad_var(3) = gradient(nodes,i,j,momentum_y  ,dx,dy)
          grad_var(4) = gradient(nodes,i,j,total_energy,dx,dy)

        end function compute_n_gradient

      end module pmodel_eq_class

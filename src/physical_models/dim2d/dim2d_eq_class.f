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
      !> 08_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_eq_class

        use dim2d_steadystate_module, only : apply_steady_state_ic
        use dim2d_dropretract_module, only : apply_drop_retraction_ic
        use dim2d_fluxes_module     , only : flux_x_mass_density,
     $                                       flux_y_mass_density,
     $                                       flux_x_momentum_x,
     $                                       flux_y_momentum_x,
     $                                       flux_x_momentum_y,
     $                                       flux_y_momentum_y,
     $                                       flux_x_total_energy,
     $                                       flux_y_total_energy
                                    
        use field_class             , only : field
        use parameters_constant     , only : scalar, vector_x, vector_y,
     $                                       steady_state,
     $                                       drop_retraction
        use parameters_input        , only : ic_choice
        use parameters_kind         , only : rkind
        use phy_model_eq_class      , only : phy_model_eq
        use cg_operators_class      , only : cg_operators

        implicit none

        private
        public :: dim2d_eq


        !> @class dim2d_eq
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
        !> @param apply_initial_conditions
        !> initialize the main variables of the governing equations
        !> considering the user choices (drop retraction, two drops
        !> collision...)
        !>
        !> @param compute_fluxes
        !> compute the fluxes along the x- and y-axis
        !---------------------------------------------------------------
        type, extends(phy_model_eq) :: dim2d_eq
          
          contains

          procedure, nopass :: get_model_name
          procedure, nopass :: get_var_name
          procedure, nopass :: get_var_longname
          procedure, nopass :: get_var_unit
          procedure, nopass :: get_var_type
          procedure, nopass :: get_eq_nb
          procedure, nopass :: apply_ic
          procedure, nopass :: compute_fluxes

        end type dim2d_eq


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
        subroutine get_var_name(var_pties)

          implicit none

          character(len=10), dimension(:), intent(inout) :: var_pties

          var_pties(1)="mass"
          var_pties(2)="momentum_x"
          var_pties(3)="momentum_y"
          var_pties(4)="energy"          

        end subroutine get_var_name
        
        
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
        subroutine get_var_longname(var_pties)

          implicit none

          character(len=32), dimension(:), intent(inout) :: var_pties

          var_pties(1)="mass density"
          var_pties(2)="momentum density along the x-axis"
          var_pties(3)="momentum density along the y-axis"
          var_pties(4)="total energy density"

        end subroutine get_var_longname


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
        subroutine get_var_unit(var_pties)

          implicit none

          character(len=10), dimension(:), intent(inout) :: var_pties

          var_pties(1)= "(kg/m3)/(kg/m3)"
          var_pties(2)= "(kg/(m2.s))/(kg/(m2.s))"
          var_pties(3)= "(kg/(m2.s))/(kg/(m2.s))"
          var_pties(4)= "(J/m3)/(J.m3)"

        end subroutine get_var_unit


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
        subroutine get_var_type(var_type)

          implicit none

          integer, dimension(:), intent(inout) :: var_type

          var_type(1)=scalar
          var_type(2)=vector_x
          var_type(3)=vector_y
          var_type(4)=scalar

        end subroutine get_var_type
        
        
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
        subroutine apply_ic(field_used)

          implicit none

          class(field), intent(inout) :: field_used


          !<read the input file to know the user choice
          select case(ic_choice)
            case(steady_state)
               call apply_steady_state_ic(field_used)
            case(drop_retraction)
               call apply_drop_retraction_ic(field_used)
            case default
               print '(''dim2d_eq_class'')'
               stop 'ic_choice not recognized'
          end select

          
          !<initialize the field depending on the user choice
          

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
        !>@param this
        !> physical model
        !>
        !>@param field_used
        !> object encapsulating the main variables
        !
        !>@param s
        !> space discretization operators
        !
        !>@param flux_x
        !> fluxes along the x-axis
        !
        !>@param flux_y
        !> fluxes along the y-axis
        !---------------------------------------------------------------
        subroutine compute_fluxes(
     $     field_used,
     $     s,
     $     flux_x,
     $     flux_y)
        
          implicit none

          class(field)                 , intent(in)   :: field_used
          type(cg_operators)           , intent(in)   :: s
          real(rkind), dimension(:,:,:), intent(inout):: flux_x
          real(rkind), dimension(:,:,:), intent(inout):: flux_y

          integer :: i,j
          integer :: bc_size


          !<get the size of the boundary layers
          bc_size = s%get_bc_size()


          !<fluxes along the x-axis
          do j=bc_size, size(flux_x,2)-bc_size
             do i=bc_size, size(flux_x,1)-bc_size

                flux_x(i,j,1) = flux_x_mass_density(field_used,s,i-1,j)
                flux_x(i,j,2) = flux_x_momentum_x(field_used,s,i-1,j)
                flux_x(i,j,3) = flux_x_momentum_y(field_used,s,i-1,j)
                flux_x(i,j,4) = flux_x_total_energy(field_used,s,i-1,j)

             end do
          end do


          !<fluxes along the y-axis
          do j=bc_size, size(flux_y,2)-bc_size
             do i=bc_size, size(flux_y,1)-bc_size

                flux_y(i,j,1) = flux_y_mass_density(field_used,s,i,j-1)
                flux_y(i,j,2) = flux_y_momentum_x(field_used,s,i,j-1)
                flux_y(i,j,3) = flux_y_momentum_y(field_used,s,i,j-1)
                flux_y(i,j,4) = flux_y_total_energy(field_used,s,i,j-1)

             end do
          end do

        end subroutine compute_fluxes

      end module dim2d_eq_class

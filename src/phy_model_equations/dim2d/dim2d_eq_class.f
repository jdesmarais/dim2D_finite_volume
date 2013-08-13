      !> @file
      !> class encapsulating subroutines to compute
      !> the governing equations of the Diffuse
      !> Interface Model in 2D
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the governing equations of the Diffuse
      !> Interface Model in 2D
      !
      !> @date
      !> 08_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_eq_class
      
        use dim2d_prim_module  , only : mass_density, momentum_x, momentum_y
     $                                  velocity_x, velocity_y, pressure,
     $                                  pressure, temperature
        use field_class        , only : field
        use parameters_constant, only : scalar, vector_x, vector_y
        use parameters_kind    , only : rkind
        use sd_operators_class , only : sd_operators

        implicit none

        private
        public :: dim2d_eq


        !> @class dim2d_eq
        !> abstract class encapsulating operators to compute
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
        type, abstract :: dim2d_eq
          
          contains

          procedure, nopass, deferred :: get_model_name
          procedure, nopass, deferred :: get_var_name
          procedure, nopass, deferred :: get_var_longname
          procedure, nopass, deferred :: get_var_unit
          procedure, nopass, deferred :: get_var_type
          procedure, nopass, deferred :: get_eq_nb
          procedure, nopass, deferred :: apply_ic
          procedure,   pass, deferred :: compute_fluxes

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

          model_name="DIM2d"

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
        !>@param sd_operators_used
        !> space discretization operators
        !
        !>@param flux_x
        !> fluxes along the x-axis
        !
        !>@param flux_y
        !> fluxes along the y-axis
        !---------------------------------------------------------------
        subroutine compute_fluxes(
     $     this,
     $     field_used,
     $     sd_operators_used,
     $     flux_x,
     $     flux_y)
        
          import field
          import dim2d_eq
          import rkind
        
          class(dim2d_eq)              , intent(in)   :: this
          class(field)                 , intent(in)   :: field_used
          class(sd_operators)          , intent(in)   :: sd_operators_used
          real(rkind), dimension(:,:,:), intent(inout):: flux_x
          real(rkind), dimension(:,:,:), intent(inout):: flux_y


          !attention au moment du calcul de ces flux
          !car il faut soit faire appel à des fonctions externes
          !pour que l'inline soit correct ou bien tout définir ici
          !mais c'est peu probable vu que l'on va vouloir tout grouper
          !dans de grandes boucles: voir au moment de la compilation...
          field_used%nodes(1,1,1)=1

        end subroutine compute_fluxes

      end module dim2d_eq_class

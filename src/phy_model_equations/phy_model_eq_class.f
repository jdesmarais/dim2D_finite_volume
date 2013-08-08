      !> @file
      !> abstract class encapsulating subroutines to compute
      !> the governing equations of the physical model
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the governing equations of the physical model
      !
      !> @date
      !> 08_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module phy_model_eq_class
      
        use field_class       , only : field
        use parameters_kind   , only : rkind
        use sd_operators_class, only : sd_operators

        implicit none

        private
        public :: phy_model_eq


        !> @class phy_model_eq
        !> abstract class encapsulating operators to compute
        !> the governing equations of the physical model
        !>
        !> @param get_model_name
        !> get the name of the physcial model
        !>
        !> @param get_var_name
        !> get the name of the main variables of the governing
        !> equations of the physical model
        !>
        !> @param get_var_longname
        !> get the description of the main variables for the
        !> governing equations of the physical model
        !>
        !> @param get_var_unit
        !> get the units of the main variables of the governing
        !> equations of the physical model
        !>
        !> @param get_var_types
        !> get the type of the main variables of the governing
        !> equations of the physical model (ex: scalar, vector_x,
        !> vector_y)
        !>
        !> @param get_eq_nb
        !> get the number of governing equations
        !>
        !> @param initialize
        !> initialize the main parameters of the physical model
        !> (ex:reynolds number)
        !>
        !> @param apply_initial_conditions
        !> initialize the main variables of the governing equations
        !> considering the user choices
        !>
        !> @param compute_fluxes
        !> compute the fluxes along the x- and y-axes
        !---------------------------------------------------------------
        type, abstract :: phy_model_equations
          
          contains

          procedure(name_model), nopass, deferred :: get_model_name
          procedure(name_var)  , nopass, deferred :: get_var_name
          procedure(lname_var) , nopass, deferred :: get_var_longname
          procedure(name_var)  , nopass, deferred :: get_var_unit
          procedure(type_var)  , nopass, deferred :: get_var_type
          procedure(gov_eq_nb) , nopass, deferred :: get_eq_nb
          procedure(ini_cond)  , nopass, deferred :: apply_ic
          procedure(fluxes)    ,   pass, deferred :: compute_fluxes

        end type phy_model_equations


        abstract interface

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
          function name_model() result(model_name)
            character(len=10) :: model_name
          end function name_model


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to get the name, description or unit
          !> of the main variables
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
          !
          !>@param var_name
          !> characters giving the variable properties
          !---------------------------------------------------------------
          subroutine name_var(var_pties)
            character(len=10), dimension(:), intent(inout) :: var_pties
          end subroutine name_var

          
           !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to get the description of the
          !> main variables
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
          !
          !>@param var_name
          !> characters giving the variable descriptions
          !---------------------------------------------------------------
          subroutine lname_var(var_pties)
            character(len=32), dimension(:), intent(inout) :: var_pties
          end subroutine lname_var


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to get the type of the main variables
          !> (ex:scalar, vector_x, vector_y)
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
          !
          !>@param var_name
          !> characters giving the variable type
          !---------------------------------------------------------------
          subroutine type_var(var_type)
            integer, dimension(:), intent(inout) :: var_type
          end subroutine type_var


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to get the number of main variables
          !> in the governing equations
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
          !
          !>@param eq_nb
          !> number of governing equations
          !---------------------------------------------------------------
          function gov_eq_nb() result(eq_nb)
            integer :: eq_nb
          end function gov_eq_nb


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
          !>@param field_used
          !> object encapsulating the main variables
          !---------------------------------------------------------------
          subroutine ini_cond(field_used)
            import field
            class(field), intent(inout) :: field_used
          end subroutine ini_cond


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
          subroutine fluxes(
     $       this,
     $       field_used,
     $       sd_operators_used,
     $       flux_x,
     $       flux_y)

            import field
            import phy_model_eq
            import rkind

            class(phy_model_eq)          , intent(in)   :: this
            class(field)                 , intent(in)   :: field_used
            class(sd_operators)          , intent(in)   :: sd_operators_used
            real(rkind), dimension(:,:,:), intent(inout):: flux_x
            real(rkind), dimension(:,:,:), intent(inout):: flux_y

          end subroutine fluxes

        end interface        

      end module phy_model_eq_class

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
        use parameters_input  , only : nx,ny,ne
        use parameters_kind   , only : rkind
        use parameters_input  , only : nx,ny,ne
        use cg_operators_class, only : cg_operators

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
        !> @param compute_flux_x
        !> compute the fluxes along the x-axe
        !>
        !> @param compute_flux_y
        !> compute the fluxes along the y-axe
        !>
        !> @param add_body_forces
        !> add the body forces to the computation of the time derivatives
        !---------------------------------------------------------------
        type, abstract :: phy_model_eq
          
          contains

          procedure(name_model), nopass, deferred :: get_model_name
          procedure(name_var)  , nopass, deferred :: get_var_name
          procedure(lname_var) , nopass, deferred :: get_var_longname
          procedure(name_var)  , nopass, deferred :: get_var_unit
          procedure(type_var)  , nopass, deferred :: get_var_type
          procedure(gov_eq_nb) , nopass, deferred :: get_eq_nb
          procedure(ini_cond)  , nopass, deferred :: apply_ic
          procedure(fluxes_x)  , nopass, deferred :: compute_flux_x
          procedure(fluxes_y)  , nopass, deferred :: compute_flux_y
          procedure(bodyforces), nopass, deferred :: compute_body_forces

        end type phy_model_eq


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
          !--------------------------------------------------------------
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
          !--------------------------------------------------------------
          function name_var() result(var_pties)
            import ne
            character(len=10), dimension(ne) :: var_pties
          end function name_var

          
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
          !--------------------------------------------------------------
          function lname_var() result(var_pties)
            import ne
            character(len=32), dimension(ne) :: var_pties
          end function lname_var


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
          !--------------------------------------------------------------
          function type_var() result(var_type)
            import ne
            integer, dimension(ne) :: var_type
          end function type_var


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
          !--------------------------------------------------------------
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
          !--------------------------------------------------------------
          subroutine ini_cond(field_used)
            import field
            class(field), intent(inout) :: field_used
          end subroutine ini_cond


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the fluxes along the
          !> x-axis
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
          !
          !>@param field_used
          !> object encapsulating the main variables
          !
          !>@param s
          !> space discretization operators
          !
          !>@param flux_x
          !> fluxes along the x-axis
          !--------------------------------------------------------------
          function fluxes_x(field_used, s) result(flux_x)

            import cg_operators
            import field
            import phy_model_eq
            import rkind
            import nx,ny,ne

            class(field)                     , intent(in)   :: field_used
            type(cg_operators)               , intent(in)   :: s
            real(rkind),dimension(nx+1,ny,ne)               :: flux_x

          end function fluxes_x

      
          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the fluxes
          !> along the y-axis
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
          !
          !>@param field_used
          !> object encapsulating the main variables
          !
          !>@param s
          !> space discretization operators
          !
          !>@param flux_y
          !> fluxes along the y-axis
          !--------------------------------------------------------------
          function fluxes_y(field_used, s) result(flux_y)

            import cg_operators
            import field
            import phy_model_eq
            import rkind
            import nx,ny,ne

            class(field)                     , intent(in)   :: field_used
            type(cg_operators)               , intent(in)   :: s
            real(rkind),dimension(nx,ny+1,ne)               :: flux_y

          end function fluxes_y


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
          !> body forces
          !--------------------------------------------------------------
          function bodyforces(field_used,s) result(body_forces)

            import field
            import rkind
            import nx,ny,ne
            import cg_operators

            class(field)                   , intent(in) :: field_used
            type(cg_operators)             , intent(in) :: s
            real(rkind),dimension(nx,ny,ne)             :: body_forces

          end function bodyforces

        end interface        

      end module phy_model_eq_class

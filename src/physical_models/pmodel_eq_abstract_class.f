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
      module pmodel_eq_abstract_class
      
        use parameters_input  , only : nx,ny,ne
        use parameters_kind   , only : rkind
        use parameters_input  , only : nx,ny,ne
        use sd_operators_class, only : sd_operators

        implicit none

        private
        public :: pmodel_eq_abstract


        !> @class pmodel_eq
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
        !> @param apply_ic
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
        type, abstract :: pmodel_eq_abstract
          
          contains

          procedure(name_model), nopass, deferred :: get_model_name
          procedure(name_var)  , nopass, deferred :: get_var_name
          procedure(lname_var) , nopass, deferred :: get_var_longname
          procedure(mname_var) , nopass, deferred :: get_var_unit
          procedure(type_var)  , nopass, deferred :: get_var_type
          procedure(gov_eq_nb) , nopass, deferred :: get_eq_nb
          procedure(ini_cond)  , nopass, deferred :: apply_ic
          procedure(fluxes_x)  , nopass, deferred :: compute_flux_x
          procedure(fluxes_y)  , nopass, deferred :: compute_flux_y
          procedure(fluxes_x_n), nopass, deferred :: compute_flux_x_nopt
          procedure(fluxes_y_n), nopass, deferred :: compute_flux_y_nopt
          procedure(bodyforces), nopass, deferred :: compute_body_forces

        end type pmodel_eq_abstract


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
          function mname_var() result(var_pties)
            import ne
            character(len=23), dimension(ne) :: var_pties
          end function mname_var


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
            character(len=33), dimension(ne) :: var_pties
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
          !>@param nodes
          !> array with the grid point data
          !--------------------------------------------------------------
          subroutine ini_cond(nodes,x_map,y_map)
            import rkind
            real(rkind), dimension(:,:,:), intent(inout) :: nodes
            real(rkind), dimension(:)    , intent(in)    :: x_map
            real(rkind), dimension(:)    , intent(in)    :: y_map
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
          !--------------------------------------------------------------
          function fluxes_x(nodes,dx,dy,s) result(flux_x)

            import rkind
            import nx,ny,ne
            import sd_operators

            real(rkind), dimension(nx,ny,ne)  , intent(in)   :: nodes
            real(rkind)                       , intent(in)   :: dx
            real(rkind)                       , intent(in)   :: dy
            type(sd_operators)                , intent(in)   :: s
            real(rkind), dimension(nx+1,ny,ne)               :: flux_x

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
          !>@param nodes
          !> array with the grid point data
          !
          !>@param s
          !> space discretization operators
          !
          !>@param flux_y
          !> fluxes along the y-axis
          !--------------------------------------------------------------
          function fluxes_y(nodes,dx,dy,s) result(flux_y)

            import rkind
            import nx,ny,ne
            import sd_operators

            real(rkind), dimension(nx,ny,ne) , intent(in)   :: nodes
            real(rkind)                      , intent(in)   :: dx
            real(rkind)                      , intent(in)   :: dy
            type(sd_operators)               , intent(in)   :: s
            real(rkind),dimension(nx,ny+1,ne)               :: flux_y

          end function fluxes_y


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the fluxes
          !> along the x-axis without knowing the exact
          !> dimensions of the input and output arrays
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
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
          !> grid step along the x-axis
          !
          !>@param grdpts_id
          !> role of the grid points
          !
          !>@param flux_x
          !> fluxes along the x-axis
          !--------------------------------------------------------------
          subroutine fluxes_x_n(nodes,dx,dy,s,grdpts_id,flux_x)
        
            import sd_operators
            import rkind

            real(rkind), dimension(:,:,:), intent(in)    :: nodes
            real(rkind)                  , intent(in)    :: dx
            real(rkind)                  , intent(in)    :: dy
            type(sd_operators)           , intent(in)    :: s
            integer    , dimension(:,:)  , intent(in)    :: grdpts_id
            real(rkind), dimension(:,:,:), intent(inout) :: flux_x

          end subroutine fluxes_x_n


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the fluxes
          !> along the y-axis without knowing the exact
          !> dimensions of the input and output arrays
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
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
          !> grid step along the x-axis
          !
          !>@param grdpts_id
          !> role of the grid points
          !
          !>@param flux_y
          !> fluxes along the y-axis
          !--------------------------------------------------------------
          subroutine fluxes_y_n(nodes,dx,dy,s,grdpts_id,flux_y)
        
            import sd_operators
            import rkind

            real(rkind), dimension(:,:,:), intent(in)    :: nodes
            real(rkind)                  , intent(in)    :: dx
            real(rkind)                  , intent(in)    :: dy
            type(sd_operators)           , intent(in)    :: s           
            integer    , dimension(:,:)  , intent(in)    :: grdpts_id
            real(rkind), dimension(:,:,:), intent(inout) :: flux_y

          end subroutine fluxes_y_n



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
          !>@param nodes
          !> array with the grid point data
          !
          !>@param k
          !> governing variables identifier
          !
          !>@param body_forces
          !> body forces
          !--------------------------------------------------------------
          function bodyforces(nodes,k) result(body_forces)

            import rkind
            import ne

            real(rkind), dimension(ne), intent(in) :: nodes
            integer                   , intent(in) :: k
            real(rkind)                            :: body_forces

          end function bodyforces

        end interface        

      end module pmodel_eq_abstract_class

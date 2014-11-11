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
      
        use interface_primary , only : gradient_x_proc,
     $                                 gradient_y_proc,
     $                                 gradient_n_proc
        use parameters_input  , only : nx,ny,ne
        use parameters_kind   , only : ikind, rkind
        use sd_operators_class, only : sd_operators

        implicit none

        private
        public :: pmodel_eq_abstract


        !> @class pmodel_eq
        !> abstract class encapsulating operators to compute
        !> the governing equations of the physical model
        !
        !> @param get_model_name
        !> get the name of the physcial model
        !
        !> @param get_var_name
        !> get the name of the main variables of the governing
        !> equations of the physical model
        !
        !> @param get_var_longname
        !> get the description of the main variables for the
        !> governing equations of the physical model
        !
        !> @param get_var_unit
        !> get the units of the main variables of the governing
        !> equations of the physical model
        !
        !> @param get_var_types
        !> get the type of the main variables of the governing
        !> equations of the physical model (ex: scalar, vector_x,
        !> vector_y)
        !
        !> @param get_sim_parameters
        !> get the simulation parameters (ex: Re, Pr, We, ...)
        !
        !> @param get_eq_nb
        !> get the number of governing equations
        !
        !> @param initialize
        !> initialize the main parameters of the physical model
        !> (ex:reynolds number)
        !
        !> @param apply_ic
        !> initialize the main variables of the governing equations
        !> considering the user choices
        !
        !> @param compute_flux_x
        !> compute the fluxes along the x-axe
        !
        !> @param compute_flux_y
        !> compute the fluxes along the y-axe
        !
        !> @param compute_body_forces
        !> add the body forces to the computation of the time derivatives
        !
        !> @param get_velocity
        !> compute the velocity out of the governing variables at a
        !> gridpoint location
        !
        !> @param are_openbc_undermined
        !> check if the open boundary conditions are undermined at the
        !> grid point location
        !
        !> @param compute_x_eigenvalues
        !> compute the eigenvalues of the hyperbolic terms in the
        !> x-direction
        !
        !> @param compute_y_eigenvalues
        !> compute the eigenvalues of the hyperbolic terms in the
        !> y-direction
        !
        !> @param compute_n1_eigenvalues
        !> compute the eigenvalues of the hyperbolic terms in the
        !> (x-y)-direction
        !
        !> @param compute_n2_eigenvalues
        !> compute the eigenvalues of the hyperbolic terms in the
        !> (x+y)-direction
        !
        !> @param compute_x_lefteigenvector
        !> compute the left eigenvectors of the hyperbolic terms in the
        !> x-direction
        !
        !> @param compute_x_righteigenvector
        !> compute the right eigenvectors of the hyperbolic terms in the
        !> x-direction
        !
        !> @param compute_y_lefteigenvector
        !> compute the left eigenvectors of the hyperbolic terms in the
        !> y-direction
        !
        !> @param compute_y_righteigenvector
        !> compute the right eigenvectors of the hyperbolic terms in the
        !> y-direction
        !
        !> @param compute_n1_lefteigenvector
        !> compute the left eigenvectors of the hyperbolic terms in the
        !> (x-y)-direction
        !
        !> @param compute_n1_righteigenvector
        !> compute the right eigenvectors of the hyperbolic terms in the
        !> (x-y)-direction
        !
        !> @param compute_n2_lefteigenvector
        !> compute the left eigenvectors of the hyperbolic terms in the
        !> (x+y)-direction
        !
        !> @param compute_n2_righteigenvector
        !> compute the right eigenvectors of the hyperbolic terms in the
        !> (x+y)-direction
        !---------------------------------------------------------------
        type, abstract :: pmodel_eq_abstract
          
          contains

          procedure(name_model)      , nopass, deferred :: get_model_name
          procedure(name_var)        , nopass, deferred :: get_var_name
          procedure(lname_var)       , nopass, deferred :: get_var_longname
          procedure(mname_var)       , nopass, deferred :: get_var_unit
          procedure(type_var)        , nopass, deferred :: get_var_type
          procedure(param_sim)       , nopass, deferred :: get_sim_parameters
          procedure(gov_eq_nb)       , nopass, deferred :: get_eq_nb
          procedure(ini_cond)        ,   pass, deferred :: apply_ic
          procedure(fluxes_x)        , nopass, deferred :: compute_flux_x
          procedure(fluxes_y)        , nopass, deferred :: compute_flux_y
          procedure(fluxes_x_n)      , nopass, deferred :: compute_flux_x_nopt
          procedure(fluxes_y_n)      , nopass, deferred :: compute_flux_y_nopt
          procedure(fluxes_x_oneside), nopass, deferred :: compute_flux_x_oneside
          procedure(fluxes_y_oneside), nopass, deferred :: compute_flux_y_oneside
          procedure(fluxes_x_byparts), nopass, deferred :: compute_flux_x_by_parts
          procedure(fluxes_y_byparts), nopass, deferred :: compute_flux_y_by_parts
          procedure(bodyforces)      , nopass, deferred :: compute_body_forces
          procedure(velocity_proc)   , nopass, deferred :: get_velocity
          procedure(v_coeff_proc)    , nopass, deferred :: get_viscous_coeff
          procedure(openbc_proc)     , nopass, deferred :: are_openbc_undermined

          procedure(eigenvalues_proc), nopass, deferred :: compute_x_eigenvalues
          procedure(eigenvalues_proc), nopass, deferred :: compute_y_eigenvalues
          procedure(eigenvalues_proc), nopass, deferred :: compute_n1_eigenvalues
          procedure(eigenvalues_proc), nopass, deferred :: compute_n2_eigenvalues

          procedure(eigenvect_proc)  , nopass, deferred :: compute_x_lefteigenvector
          procedure(eigenvect_proc)  , nopass, deferred :: compute_x_righteigenvector
          procedure(eigenvect_proc)  , nopass, deferred :: compute_y_lefteigenvector
          procedure(eigenvect_proc)  , nopass, deferred :: compute_y_righteigenvector
          procedure(eigenvect_proc)  , nopass, deferred :: compute_n1_lefteigenvector
          procedure(eigenvect_proc)  , nopass, deferred :: compute_n1_righteigenvector
          procedure(eigenvect_proc)  , nopass, deferred :: compute_n2_lefteigenvector
          procedure(eigenvect_proc)  , nopass, deferred :: compute_n2_righteigenvector

          procedure(eigenvect_proc)  , nopass, deferred :: compute_cons_lodi_matrix_x
          procedure(eigenvect_proc)  , nopass, deferred :: compute_cons_lodi_matrix_y

          procedure(x_gradient_proc) , nopass, deferred :: compute_x_gradient
          procedure(y_gradient_proc) , nopass, deferred :: compute_y_gradient
          procedure(n_gradient_proc) , nopass, deferred :: compute_n_gradient

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
          !> interface to get the main parameters of the
          !> simulation
          !
          !> @date
          !> 12_08_2014 - initial version - J.L. Desmarais
          !
          !>@param var_name
          !> characters giving the variable type
          !--------------------------------------------------------------
          subroutine param_sim(param_name, param_value)
            import rkind
            character(10), dimension(:), allocatable, intent(out) :: param_name
            real(rkind)  , dimension(:), allocatable, intent(out) :: param_value
          end subroutine param_sim


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
          subroutine ini_cond(this,nodes,x_map,y_map)
            import pmodel_eq_abstract
            import rkind
            
            class(pmodel_eq_abstract)    , intent(in)    :: this
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
          subroutine fluxes_x_n(
     $      nodes,dx,dy,s,
     $      grdpts_id,
     $      flux_x,
     $      x_borders,y_borders)
        
            import ikind
            import rkind
            import sd_operators

            real(rkind)   , dimension(:,:,:), intent(in)    :: nodes
            real(rkind)                     , intent(in)    :: dx
            real(rkind)                     , intent(in)    :: dy
            type(sd_operators)              , intent(in)    :: s
            integer       , dimension(:,:)  , intent(in)    :: grdpts_id
            real(rkind)   , dimension(:,:,:), intent(inout) :: flux_x
            integer(ikind), dimension(2)    , intent(in)    :: x_borders
            integer(ikind), dimension(2)    , intent(in)    :: y_borders

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
          subroutine fluxes_y_n(
     $      nodes,dx,dy,s,
     $      grdpts_id,
     $      flux_y,
     $      x_borders, y_borders)
        
            import ikind
            import rkind
            import sd_operators

            real(rkind), dimension(:,:,:), intent(in)    :: nodes
            real(rkind)                  , intent(in)    :: dx
            real(rkind)                  , intent(in)    :: dy
            type(sd_operators)           , intent(in)    :: s           
            integer    , dimension(:,:)  , intent(in)    :: grdpts_id
            real(rkind), dimension(:,:,:), intent(inout) :: flux_y
            integer    , dimension(2)    , intent(in)    :: x_borders
            integer    , dimension(2)    , intent(in)    :: y_borders

          end subroutine fluxes_y_n


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
          !>@param i
          !> x-index where the flux_x is computed
          !
          !>@param j
          !> y-index where the flux_x is computed
          !
          !>@param dy
          !> grid size along the y-axis
          !
          !>@param s_oneside
          !> space discretization operators
          !
          !>@param flux_x
          !> fluxes along the x-axis
          !--------------------------------------------------------------
          function fluxes_x_oneside(
     $      nodes,dx,dy,
     $      i,j,
     $      s_oneside)
     $      result(flux_x)

            import ikind
            import ne
            import rkind
            import sd_operators

            real(rkind), dimension(:,:,:), intent(in)   :: nodes
            real(rkind)                  , intent(in)   :: dx
            real(rkind)                  , intent(in)   :: dy
            integer(ikind)               , intent(in)   :: i
            integer(ikind)               , intent(in)   :: j
            class(sd_operators)          , intent(in)   :: s_oneside
            real(rkind), dimension(ne)                  :: flux_x

          end function fluxes_x_oneside


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
          !>@param i
          !> x-index where the flux_x is computed
          !
          !>@param j
          !> y-index where the flux_x is computed
          !
          !>@param dy
          !> grid size along the y-axis
          !
          !>@param s_oneside
          !> space discretization operators
          !
          !>@param flux_y
          !> fluxes along the x-axis
          !--------------------------------------------------------------
          function fluxes_y_oneside(
     $      nodes,dx,dy,
     $      i,j,
     $      s_oneside)
     $      result(flux_y)

            import ikind
            import ne
            import rkind
            import sd_operators

            real(rkind), dimension(:,:,:), intent(in)   :: nodes
            real(rkind)                  , intent(in)   :: dx
            real(rkind)                  , intent(in)   :: dy
            integer(ikind)               , intent(in)   :: i
            integer(ikind)               , intent(in)   :: j
            class(sd_operators)          , intent(in)   :: s_oneside
            real(rkind), dimension(ne)                  :: flux_y

          end function fluxes_y_oneside


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the fluxes along the
          !> x-axis
          !
          !> @date
          !> 10_11_2014 - initial version - J.L. Desmarais
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
          !>@param i
          !> x-index where the flux_x is computed
          !
          !>@param j
          !> y-index where the flux_x is computed
          !
          !>@param dy
          !> grid size along the y-axis
          !
          !>@param s_oneside
          !> space discretization operators
          !
          !>@param flux_x
          !> fluxes along the x-axis
          !--------------------------------------------------------------
          function fluxes_x_byparts(
     $      nodes,dx,dy,i,j,s_oneside,
     $      inviscid_flux, viscid_flux)
     $      result(flux_x)

            import ikind
            import ne
            import rkind
            import sd_operators

            real(rkind), dimension(:,:,:), intent(in)   :: nodes
            real(rkind)                  , intent(in)   :: dx
            real(rkind)                  , intent(in)   :: dy
            integer(ikind)               , intent(in)   :: i
            integer(ikind)               , intent(in)   :: j
            class(sd_operators)          , intent(in)   :: s_oneside
            real(rkind), dimension(ne)   , intent(out)  :: inviscid_flux
            real(rkind), dimension(ne)   , intent(out)  :: viscid_flux
            real(rkind), dimension(ne)                  :: flux_x

          end function fluxes_x_byparts


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the fluxes along the
          !> y-axis
          !
          !> @date
          !> 10_11_2014 - initial version - J.L. Desmarais
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
          !>@param i
          !> x-index where the flux_x is computed
          !
          !>@param j
          !> y-index where the flux_x is computed
          !
          !>@param dy
          !> grid size along the y-axis
          !
          !>@param s_oneside
          !> space discretization operators
          !
          !>@param flux_y
          !> fluxes along the x-axis
          !--------------------------------------------------------------
          function fluxes_y_byparts(
     $      nodes,dx,dy,i,j,s_oneside,
     $      inviscid_flux, viscid_flux)
     $      result(flux_y)

            import ikind
            import ne
            import rkind
            import sd_operators

            real(rkind), dimension(:,:,:), intent(in)   :: nodes
            real(rkind)                  , intent(in)   :: dx
            real(rkind)                  , intent(in)   :: dy
            integer(ikind)               , intent(in)   :: i
            integer(ikind)               , intent(in)   :: j
            class(sd_operators)          , intent(in)   :: s_oneside
            real(rkind), dimension(ne)   , intent(out)  :: inviscid_flux
            real(rkind), dimension(ne)   , intent(out)  :: viscid_flux
            real(rkind), dimension(ne)                  :: flux_y

          end function fluxes_y_byparts



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
          function bodyforces(t,x,y,nodes,k) result(body_forces)

            import rkind
            import ne

            real(rkind)               , intent(in) :: t
            real(rkind)               , intent(in) :: x
            real(rkind)               , intent(in) :: y
            real(rkind), dimension(ne), intent(in) :: nodes
            integer                   , intent(in) :: k
            real(rkind)                            :: body_forces

          end function bodyforces


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
          function velocity_proc(nodes) result(velocity)

            import rkind
            import ne

            real(rkind), dimension(ne), intent(in) :: nodes
            real(rkind), dimension(2)              :: velocity

          end function velocity_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> get the viscous constant
          !
          !> @date
          !> 11_11_2014 - initial version - J.L. Desmarais
          !
          !>@return viscous_coeff
          !> viscous coefficient
          !-------------------------------------------------------------
          function v_coeff_proc() result(viscous_coeff)
          
            import rkind
          
            real(rkind) :: viscous_coeff
          
          end function v_coeff_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface checking whether the open boundary conditions
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
          function openbc_proc(nodes) result(undermined)

            import rkind
            import ne

            real(rkind), dimension(ne), intent(in) :: nodes
            logical                                :: undermined

          end function openbc_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for the local computation of the eigenvalues
          !> for the hyperbolic terms in the x-direction
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
          function eigenvalues_proc(nodes) result(eigenvalues)

            import rkind
            import ne

            real(rkind), dimension(ne), intent(in) :: nodes
            real(rkind), dimension(ne)             :: eigenvalues

          end function eigenvalues_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for the local computation of the eigenvalues
          !> for the hyperbolic terms in the x-direction
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
          function eigenvect_proc(nodes) result(eigenvect)

            import rkind
            import ne

            real(rkind), dimension(ne), intent(in) :: nodes
            real(rkind), dimension(ne,ne)          :: eigenvect

          end function eigenvect_proc


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
          function x_gradient_proc(nodes,i,j,gradient,dx) result(grad_var)

            import gradient_x_proc
            import ikind,rkind
            import ne

            real(rkind), dimension(:,:,:), intent(in) :: nodes
            integer(ikind)               , intent(in) :: i
            integer(ikind)               , intent(in) :: j
            procedure(gradient_x_proc)                :: gradient
            real(rkind)                  , intent(in) :: dx
            real(rkind), dimension(ne)                :: grad_var

          end function x_gradient_proc


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
          function y_gradient_proc(nodes,i,j,gradient,dy) result(grad_var)

            import gradient_y_proc
            import ikind,rkind
            import ne

            real(rkind), dimension(:,:,:), intent(in) :: nodes
            integer(ikind)               , intent(in) :: i
            integer(ikind)               , intent(in) :: j
            procedure(gradient_y_proc)                :: gradient
            real(rkind)                  , intent(in) :: dy
            real(rkind), dimension(ne)                :: grad_var

          end function y_gradient_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for the computation of the gradient of the
          !> governing variables in a diagonal direction
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
          !> procedure used to compute the gradient along the diagonal
          !> direction
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
          function n_gradient_proc(nodes,i,j,gradient,dx,dy) result(grad_var)

            import gradient_n_proc
            import ikind,rkind
            import nx,ny,ne

            real(rkind), dimension(:,:,:), intent(in) :: nodes
            integer(ikind)               , intent(in) :: i
            integer(ikind)               , intent(in) :: j
            procedure(gradient_n_proc)                :: gradient
            real(rkind)                  , intent(in) :: dx
            real(rkind)                  , intent(in) :: dy
            real(rkind), dimension(ne)                :: grad_var            

          end function n_gradient_proc


        end interface        

      end module pmodel_eq_abstract_class

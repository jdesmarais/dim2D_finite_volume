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
      !> 08_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module pmodel_eq_abstract_class
      
        use interface_primary , only :
     $       gradient_proc

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use sd_operators_class, only :
     $       sd_operators

        implicit none

        private
        public :: pmodel_eq_abstract


        !> @class pmodel_eq
        !> abstract class encapsulating operators to compute
        !> the governing equations of the physical model
        !---------------------------------------------------------------
        type, abstract :: pmodel_eq_abstract
          
          contains

          !description of the model
          procedure(name_model)        , nopass, deferred :: get_model_name        !<@brief get the name of the physcial model
          procedure(name_var)          , nopass, deferred :: get_var_name          !<@brief get the name of the main variables of the governing equations of the physical model
          procedure(lname_var)         , nopass, deferred :: get_var_longname      !<@brief get the description of the main variables for the governing equations of the physical model
          procedure(mname_var)         , nopass, deferred :: get_var_unit          !<@brief get the units of the main variables of the governing equations of the physical model
          procedure(type_var)          , nopass, deferred :: get_var_type          !<@brief get the type of the main variables of the governing equations of the physical model (ex: scalar, vector_x, vector_y)
          procedure(param_sim)         , nopass, deferred :: get_sim_parameters    !<@brief get the simulation parameters (ex: Re, Pr, We, ...)
          procedure(gov_eq_nb)         , nopass, deferred :: get_eq_nb             !<@brief get the number of governing equations
              

          !sd operators pattern for the fluxes
          procedure(sd_pattern)        , nopass, deferred :: get_sd_pattern_flux_x !<@brief identify the number of grid points needed from the flux computation along the x-direction
          procedure(sd_pattern)        , nopass, deferred :: get_sd_pattern_flux_y !<@brief identify the number of grid points needed from the flux computation along the y-direction
            

          !initial conditions procedures
          procedure(ini_cond)          ,   pass, deferred :: apply_ic !<@brief initialize the main variables of the governing equations considering the user choices
               

          !flux computation
          procedure(fluxes_x)          , nopass, deferred :: compute_flux_x          !<@brief compute the fluxes along the x-axis
          procedure(fluxes_y)          , nopass, deferred :: compute_flux_y          !<@brief compute the fluxes along the y-axis
          procedure(fluxes_x_n)        , nopass, deferred :: compute_flux_x_nopt     !<@brief compute the fluxes along the x-axis for non-fixed size array
          procedure(fluxes_y_n)        , nopass, deferred :: compute_flux_y_nopt     !<@brief compute the fluxes along the y-axis for non-fixed size array
          procedure(fluxes_x_oneside)  , nopass, deferred :: compute_flux_x_oneside  !<@brief compute the flux along the x-axis at the boundary
          procedure(fluxes_y_oneside)  , nopass, deferred :: compute_flux_y_oneside  !<@brief compute the flux along the y-axis at the boundary
          procedure(fluxes_x_byparts)  , nopass, deferred :: compute_flux_x_by_parts !<@brief compute the flux along the x-axis by distinguishing inviscid, viscid and capillary parts
          procedure(fluxes_y_byparts)  , nopass, deferred :: compute_flux_y_by_parts !<@brief compute the flux along the y-axis by distinguishing inviscid, viscid and capillary parts
          procedure(bodyforces)        , nopass, deferred :: compute_body_forces     !<@brief compute the body forces


          !field extension for openb b.c.
          procedure(velocity_proc)     , nopass, deferred :: get_velocity            !<@brief compute the velocity out of the governing variables at a gridpoint location
          procedure(v_coeff_proc)      , nopass, deferred :: get_viscous_coeff       !<@brief determine the viscous ratio
          procedure(openbc_proc)       , nopass, deferred :: are_openbc_undermined   !<@brief check if the open boundary conditions are undermined at the grid point location
          procedure(farfield_proc)     ,   pass, deferred :: get_far_field           !<@brief determine the governing variables imposed in the far field
          procedure(openbc_prim_proc)  ,   pass, deferred :: get_prim_obc_eigenqties !<@brief determine the primitive variables needed to compute the eigenquantities in the far field
          

          !computations with primitive variables
          procedure(openbc_var_proc)   , nopass, deferred :: compute_prim_var !<@brief compute the primitive variables from the conservative variables
          procedure(openbc_var_proc)   , nopass, deferred :: compute_cons_var !<@brief compute the conservative variables from the primitive variables
          
          procedure(openbc_matrix_proc), nopass, deferred :: compute_jacobian_prim_to_cons !<@brief compute the Jacobian matrix from primitive to conservative variables
          procedure(openbc_matrix_proc), nopass, deferred :: compute_jacobian_cons_to_prim !<@brief compute the Jacobian matrix from conservative to primitive variables

          procedure(openbc_matrix_proc), nopass, deferred :: compute_x_transM_prim !<@brief compute the transverse matrix along the x-direction using primitive variables
          procedure(openbc_matrix_proc), nopass, deferred :: compute_y_transM_prim !<@brief compute the transverse matrix along the y-direction using primitive variables

          procedure(openbc_vector_proc), nopass, deferred :: compute_x_eigenvalues_prim !<@brief compute the eigenvalues of the hyperbolic terms in the x-direction
          procedure(openbc_vector_proc), nopass, deferred :: compute_y_eigenvalues_prim !<@brief compute the eigenvalues of the hyperbolic terms in the y-direction

          procedure(openbc_matrix_proc), nopass, deferred :: compute_x_lefteigenvector_prim  !<@brief compute the left eigenvectors of the hyperbolic terms in the x-direction
          procedure(openbc_matrix_proc), nopass, deferred :: compute_x_righteigenvector_prim !<@brief compute the right eigenvectors of the hyperbolic terms in the x-direction
          procedure(openbc_matrix_proc), nopass, deferred :: compute_y_lefteigenvector_prim  !<@brief compute the left eigenvectors of the hyperbolic terms in the y-direction
          procedure(openbc_matrix_proc), nopass, deferred :: compute_y_righteigenvector_prim !<@brief compute the right eigenvectors of the hyperbolic terms in the y-direction

          procedure(grad_prim_proc)    , nopass, deferred :: compute_gradient_prim !<@brief compute the gradient of the primitive variables


          !variables in the rotated frame
          procedure(xy_to_n_proc)     , nopass, deferred :: compute_xy_to_n_var !<@brief transform a vector from the (x,y)-cartesian frame to the (n1,n2)-frame
          procedure(n_to_xy_proc)     , nopass, deferred :: compute_n_to_xy_var !<@brief transform a vector from the (n1,n2)-frame to the (x,y)-cartesian frame


c$$$          procedure(eigenvalues_proc) , nopass, deferred :: compute_x_eigenvalues
c$$$          procedure(eigenvalues_proc) , nopass, deferred :: compute_y_eigenvalues
c$$$          procedure(eigenvalues_proc) , nopass, deferred :: compute_n1_eigenvalues
c$$$          procedure(eigenvalues_proc) , nopass, deferred :: compute_n2_eigenvalues
c$$$                                      
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_x_lefteigenvector
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_x_righteigenvector
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_y_lefteigenvector
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_y_righteigenvector
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_n1_lefteigenvector
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_n1_righteigenvector
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_n2_lefteigenvector
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_n2_righteigenvector
c$$$                                      
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_x_leftConsLodiM
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_y_leftConslodiM
c$$$          procedure(lodi_td_proc)     , nopass, deferred :: compute_x_timedev_from_LODI_vector
c$$$          procedure(lodi_td_proc)     , nopass, deferred :: compute_y_timedev_from_LODI_vector
c$$$          procedure(lodi_tds_proc)    , nopass, deferred :: compute_timedev_from_LODI_vectors
c$$$                                      
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_x_transM
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_y_transM
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_n1_transM
c$$$          procedure(eigenvect_proc)   , nopass, deferred :: compute_n2_transM
c$$$                                      
c$$$          procedure(farfield_proc)    ,   pass, deferred :: get_far_field
c$$$                                      
c$$$          procedure(x_gradient_proc)  , nopass, deferred :: compute_x_gradient
c$$$          procedure(y_gradient_proc)  , nopass, deferred :: compute_y_gradient
c$$$          procedure(n_gradient_proc)  , nopass, deferred :: compute_n_gradient

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
          !>@return model_name
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
          !>@return var_pties
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
          !>@return var_pties
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
          !> interface to get the long description of the
          !> main variables
          !
          !> @date
          !> 08_08_2013 - initial version - J.L. Desmarais
          !
          !>@return var_pties
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
          !>@return var_type
          !> integer describing the governing variable types
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
          !>@param param_name
          !> name of the main parameters of the simulation
          !
          !>@param param_value
          !> value of the main parameters of the simulation
          !--------------------------------------------------------------
          subroutine param_sim(param_name, param_value)
            import rkind
            character(20), dimension(:), allocatable, intent(out) :: param_name
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
          !>@return eq_nb
          !> number of governing equations
          !--------------------------------------------------------------
          function gov_eq_nb() result(eq_nb)
            integer :: eq_nb
          end function gov_eq_nb


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to get the pattern of gridpoints needed
          !> when computing fluxes
          !
          !> @date
          !> 27_01_2015 - initial version - J.L. Desmarais
          !
          !> @param operator_type
          !> type of operator used
          !
          !> @return pattern
          !> gridpoints needed around the central gridpoint to compute
          !> the fluxes
          !--------------------------------------------------------------
          function sd_pattern(operator_type) result(pattern)

            implicit none

            integer    , intent(in) :: operator_type
            integer, dimension(2,2) :: pattern

          end function sd_pattern


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
          !>@param x_map
          !> array with the x-coordinates
          !
          !>@param y_map
          !> array with the y-coordinates                
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
          !>@return flux_x
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
          !>@param dx
          !> grid size along the x-axis
          !
          !>@param dy
          !> grid size along the y-axis
          !
          !>@param s
          !> space discretization operators
          !
          !>@return flux_y
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
          !>@param dx
          !> grid step along the x-axis
          !
          !>@param dy
          !> grid step along the x-axis
          !
          !>@param s
          !> space discretization operators          
          !
          !>@param grdpts_id
          !> role of the grid points
          !
          !>@param flux_x
          !> fluxes along the x-axis
          !
          !>@param x_borders
          !> integration borders along the x-axis
          !
          !>@param y_borders
          !> integration borders along the y-axis
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
          !>@param dx
          !> grid step along the x-axis
          !
          !>@param dy
          !> grid step along the x-axis
          !
          !>@param s
          !> space discretization operators          
          !
          !>@param grdpts_id
          !> role of the grid points
          !
          !>@param flux_y
          !> fluxes along the x-axis
          !
          !>@param x_borders
          !> integration borders along the x-axis
          !
          !>@param y_borders
          !> integration borders along the y-axis
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
          !>@param s_oneside
          !> space discretization operators
          !
          !>@return flux_x
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
          !>@param s_oneside
          !> space discretization operators
          !
          !>@return flux_y
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
          !>@param s_oneside
          !> space discretization operators
          !
          !>@param inviscid_flux
          !> inviscid part of the flux along the x-axis
          !
          !>@param viscid_flux
          !> viscid part of the flux along the x-axis
          !
          !>@return flux_x
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
          !>@param s_oneside
          !> space discretization operators
          !
          !>@param inviscid_flux
          !> inviscid part of the flux along the y-axis
          !
          !>@param viscid_flux
          !> viscid part of the flux along the y-axis
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
          !>@param t
          !> time
          !
          !>@param x
          !> x-coordinate
          !
          !>@param y
          !> y-coordinate
          !
          !>@param nodes
          !> governing variable vector
          !
          !>@param k
          !> governing variables identifier
          !
          !>@param prim
          !> whether the body forces are computed for primitive
          !> variable governing equations
          !
          !>@return body_forces
          !> body forces
          !--------------------------------------------------------------
          function bodyforces(t,x,y,nodes,k,prim) result(body_forces)

            import rkind
            import ne

            real(rkind)               , intent(in) :: t
            real(rkind)               , intent(in) :: x
            real(rkind)               , intent(in) :: y
            real(rkind), dimension(ne), intent(in) :: nodes
            integer                   , intent(in) :: k
            logical, optional         , intent(in) :: prim
            real(rkind)                            :: body_forces

          end function bodyforces


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the velocity vector
          !> from the conservative variables
          !
          !> @date
          !> 23_09_2013 - initial version - J.L. Desmarais
          !
          !>@param nodes
          !> governing variable vector \f$ (\rho,q_x,q_y,\rho E) \f$
          !
          !>@return velocity
          !> velocity vector \f$ (v_x,v_y) \f$
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
          !> get the viscous constant ratio
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
          !>@param x_map
          !> x-coordinate map
          !
          !>@param y_map
          !> y-coordinate map
          !
          !>@param nodes
          !> array with the grid point data
          !
          !>@return undermined
          !> check if the open boundary conditions are undermined
          !> at the grid point location
          !--------------------------------------------------------------
          function openbc_proc(x_map,y_map,nodes) result(undermined)

            import rkind
            import ne

            real(rkind), dimension(3)     , intent(in) :: x_map
            real(rkind), dimension(3)     , intent(in) :: y_map
            real(rkind), dimension(3,3,ne), intent(in) :: nodes
            logical                                    :: undermined

          end function openbc_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for the computation of the governing variables
          !> in the far field as chosen by the initial conditions
          !
          !> @date
          !> 13_11_2014 - initial version - J.L. Desmarais
          !
          !>@param this
          !> physical model
          !
          !>@param t
          !> time
          !
          !>@param x
          !> x-coordinate
          !
          !>@param y
          !> y-coordinate
          !
          !>@return var
          !> governing variables in the far field
          !--------------------------------------------------------------
          function farfield_proc(this,t,x,y) result(var)

            import ne
            import pmodel_eq_abstract
            import rkind

            class(pmodel_eq_abstract), intent(in) :: this
            real(rkind)              , intent(in) :: t
            real(rkind)              , intent(in) :: x
            real(rkind)              , intent(in) :: y
            real(rkind), dimension(ne)            :: var
            
          end function farfield_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface determining the grid points used to evaluate
          !> the eigenquantities at the edge of the computational
          !> domain
          !
          !> @date
          !> 02_02_2015 - initial version - J.L. Desmarais
          !
          !>@param this
          !> physical model
          !
          !>@param t
          !> time
          !
          !>@param x
          !> x-coordinate of the grid points at the boundary
          !
          !>@param y
          !> y-coordinate of the grid points at the boundary
          !
          !>@param nodes_bc
          !> array with the grid point data at the boundary
          !
          !>@return nodes_prim_extended
          !> extended primitive variable vector to evaluate the
          !> eigenquantities at the boundary \f$ (\rho,v_x,v_y,P,c) \f$
          !--------------------------------------------------------------
          function openbc_prim_proc(this,t,x,y,nodes_bc)
     $      result(nodes_prim_extended)

            import ne
            import pmodel_eq_abstract
            import rkind

            class(pmodel_eq_abstract)   , intent(in) :: this
            real(rkind)                 , intent(in) :: t
            real(rkind)                 , intent(in) :: x
            real(rkind)                 , intent(in) :: y
            real(rkind), dimension(ne)  , intent(in) :: nodes_bc
            real(rkind), dimension(ne+1)             :: nodes_prim_extended

          end function openbc_prim_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for the local computation of
          !> the eigenquantity vectors of the hyperbolic terms
          !
          !> @date
          !> 01_08_2014 - initial version - J.L. Desmarais
          !
          !>@param nodes_in
          !> array with the conservative variables
          !
          !>@return nodes_out
          !> array with eigen-quantities
          !--------------------------------------------------------------
          function openbc_var_proc(nodes_in) result(nodes_out)

            import rkind
            import ne

            real(rkind), dimension(ne), intent(in) :: nodes_in
            real(rkind), dimension(ne)             :: nodes_out

          end function openbc_var_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for the local computation of
          !> eigenquantity matrix of the hyperbolic terms
          !
          !> @date
          !> 01_08_2014 - initial version - J.L. Desmarais
          !
          !>@param nodes_prim_extended
          !> vector with the primitive variables
          !
          !>@return matrix
          !> eigenquantity matrix at the location of the grid point
          !--------------------------------------------------------------
          function openbc_matrix_proc(nodes_prim_extended) result(matrix)

            import rkind
            import ne

            real(rkind), dimension(ne+1) , intent(in) :: nodes_prim_extended
            real(rkind), dimension(ne,ne)             :: matrix

          end function openbc_matrix_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for the local computation of the
          !> eigenquantity vector of the hyperbolic terms
          !
          !> @date
          !> 01_08_2014 - initial version - J.L. Desmarais
          !
          !>@param nodes_prim_extended
          !> vector with the primitive variables
          !
          !>@return vector
          !> eigenquantity vector at the location of the grid point
          !--------------------------------------------------------------
          function openbc_vector_proc(nodes_prim_extended) result(vector)

            import rkind
            import ne

            real(rkind), dimension(ne+1), intent(in) :: nodes_prim_extended
            real(rkind), dimension(ne)               :: vector

          end function openbc_vector_proc


c$$$          !> @author
c$$$          !> Julien L. Desmarais
c$$$          !
c$$$          !> @brief
c$$$          !> interface for the local computation of the contribution
c$$$          !> of the LODI vector in one direction (x or y) to the time
c$$$          !> derivatives
c$$$          !
c$$$          !> @date
c$$$          !> 12_12_2014 - initial version - J.L. Desmarais
c$$$          !
c$$$          !>@param nodes
c$$$          !> array with the grid point data
c$$$          !
c$$$          !>@param lodi
c$$$          !> LODI vector
c$$$          !
c$$$          !>@return timedev
c$$$          !> time derivatives
c$$$          !--------------------------------------------------------------
c$$$          function lodi_td_proc(nodes,lodi) result(timedev)
c$$$
c$$$            import rkind
c$$$            import ne
c$$$
c$$$            real(rkind), dimension(ne), intent(in) :: nodes
c$$$            real(rkind), dimension(ne), intent(in) :: lodi
c$$$            real(rkind), dimension(ne)             :: timedev
c$$$
c$$$          end function lodi_td_proc
c$$$
c$$$
c$$$          !> @author
c$$$          !> Julien L. Desmarais
c$$$          !
c$$$          !> @brief
c$$$          !> interface for the local computation of the contribution
c$$$          !> of the LODI vectors in both directions (x and y) to the time
c$$$          !> derivatives
c$$$          !
c$$$          !> @date
c$$$          !> 12_12_2014 - initial version - J.L. Desmarais
c$$$          !
c$$$          !>@param nodes
c$$$          !> array with the grid point data
c$$$          !
c$$$          !>@param lodi
c$$$          !> LODI vector
c$$$          !
c$$$          !>@return timedev
c$$$          !> time derivatives
c$$$          !--------------------------------------------------------------
c$$$          function lodi_tds_proc(nodes,lodi_x,lodi_y) result(timedev)
c$$$
c$$$            import rkind
c$$$            import ne
c$$$
c$$$            real(rkind), dimension(ne), intent(in) :: nodes
c$$$            real(rkind), dimension(ne), intent(in) :: lodi_x
c$$$            real(rkind), dimension(ne), intent(in) :: lodi_y
c$$$            real(rkind), dimension(ne)             :: timedev
c$$$
c$$$          end function lodi_tds_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for the computation of the gradient of the
          !> governing variables
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
          !> procedure used to compute the gradient along the axis
          !
          !>@param dn
          !> grid space step along the direction for the gradient
          !
          !>@param use_n_dir
          !> choose to use the diagonal direction instead of one
          !> of the cartesian direction (x,y)
          !
          !>@return grad_var
          !> gradient of the governing variables along the x-axis
          !--------------------------------------------------------------
          function grad_prim_proc(nodes,i,j,gradient,dn,use_n_dir)
     $      result(grad_var)

            import gradient_proc
            import ikind,rkind
            import ne

            real(rkind), dimension(:,:,:), intent(in) :: nodes
            integer(ikind)               , intent(in) :: i
            integer(ikind)               , intent(in) :: j
            procedure(gradient_proc)                  :: gradient
            real(rkind)                  , intent(in) :: dn
            logical    , optional        , intent(in) :: use_n_dir
            real(rkind), dimension(ne)                :: grad_var

          end function grad_prim_proc

      
          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for the computation of the variables
          !> in the (n1,n2) frame from the (x,y) frame
          !
          !> @date
          !> 01_08_2014 - initial version - J.L. Desmarais
          !
          !>@param nodes
          !> governing variables vector in the (x,y) frame
          !
          !>@return nodes_n
          !> governing variables vector in the (n1,n2) frame
          !--------------------------------------------------------------
          function xy_to_n_proc(nodes) result(nodes_n)

            import rkind
            import ne
          
            real(rkind), dimension(ne), intent(in) :: nodes
            real(rkind), dimension(ne)             :: nodes_n

          end function xy_to_n_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for the computation of the variables
          !> in the (x,y) frame from the (n1,n2) frame
          !
          !> @date
          !> 01_08_2014 - initial version - J.L. Desmarais
          !
          !>@param nodes_n
          !> governing variables vector in the (n1,n2) frame
          !
          !>@return nodes
          !> governing variables vector in the (x,y) frame
          !--------------------------------------------------------------
          function n_to_xy_proc(nodes_n) result(nodes)

            import rkind
            import ne

            real(rkind), dimension(ne), intent(in) :: nodes_n
            real(rkind), dimension(ne)             :: nodes

          end function n_to_xy_proc

        end interface        

      end module pmodel_eq_abstract_class

      !> @file
      !> class implementing default subroutines for
      !> the physical model to throw exceptions when
      !> used by the program
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class implementing default subroutines for
      !> the physical model to throw exceptions when
      !> used by the program
      !
      !> @date
      !> 08_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module pmodel_eq_default_class
      
        use interface_primary, only :
     $       gradient_proc

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_abstract_class, only :
     $       pmodel_eq_abstract

        use sd_operators_class, only : 
     $       sd_operators

        implicit none

        private
        public :: pmodel_eq_default,
     $            basic


        !> @class pmodel_eq_default
        !> class implementing default subroutines for
        !> the physical model to throw exceptions when
        !> used by the program
        !---------------------------------------------------------------
        type, abstract, extends(pmodel_eq_abstract) :: pmodel_eq_default
          
          contains

          procedure, nopass :: get_sim_parameters => get_sim_parameters_default

          procedure, nopass :: get_sd_pattern_flux_x => get_sd_pattern_flux_default
          procedure, nopass :: get_sd_pattern_flux_y => get_sd_pattern_flux_default

          procedure,   pass :: apply_ic => apply_ic_default

          procedure, nopass :: compute_flux_x_by_parts => compute_flux_x_by_parts_default
          procedure, nopass :: compute_flux_y_by_parts => compute_flux_y_by_parts_default

          procedure, nopass :: get_velocity      => get_velocity_default
          procedure, nopass :: get_viscous_coeff => get_viscous_coeff_default
          
          procedure, nopass :: are_openbc_undermined   => are_openbc_undermined_default
          procedure,   pass :: get_far_field           => get_far_field_default
          procedure,   pass :: get_prim_obc_eigenqties => get_prim_obc_eigenqties_default
          
          !computations with primitive variables
          procedure, nopass :: compute_prim_var => compute_openbc_var_default
          procedure, nopass :: compute_cons_var => compute_openbc_var_default
          
          procedure, nopass :: compute_jacobian_prim_to_cons => compute_openbc_matrix_default
          procedure, nopass :: compute_jacobian_cons_to_prim => compute_openbc_matrix_default

          procedure, nopass :: compute_x_transM_prim => compute_openbc_matrix_default
          procedure, nopass :: compute_y_transM_prim => compute_openbc_matrix_default

          procedure, nopass :: compute_x_eigenvalues_prim => compute_openbc_vector_default
          procedure, nopass :: compute_y_eigenvalues_prim => compute_openbc_vector_default

          procedure, nopass :: compute_x_lefteigenvector_prim  => compute_openbc_matrix_default
          procedure, nopass :: compute_x_righteigenvector_prim => compute_openbc_matrix_default
          procedure, nopass :: compute_y_lefteigenvector_prim  => compute_openbc_matrix_default
          procedure, nopass :: compute_y_righteigenvector_prim => compute_openbc_matrix_default

          procedure, nopass :: compute_gradient_prim => compute_gradient_prim_default

          procedure, nopass :: compute_xy_to_n_var => compute_xy_to_n_var_default
          procedure, nopass :: compute_n_to_xy_var => compute_n_to_xy_var_default

        end type pmodel_eq_default


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default subroutine to get the main
        !> parameters of the simulation
        !
        !> @date
        !> 12_08_2014 - initial version - J.L. Desmarais
        !
        !>@param param_name
        !> names of the main parameters of the simulation
        !
        !>@param param_value
        !> values of the main parameters of the simulation
        !--------------------------------------------------------------
        subroutine get_sim_parameters_default(param_name, param_value)

          implicit none

          character(20), dimension(:), allocatable, intent(out) :: param_name
          real(rkind)  , dimension(:), allocatable, intent(out) :: param_value


          print '(''pmodel_eq_default'')'
          print '(''get_sim_parameters'')'

          if(allocated(param_name)) then
             stop 'get_sim_parameter: param_name allocated'
          end if

          if(allocated(param_value)) then
             stop 'get_sim_parameter: param_value allocated'
          end if

        end subroutine get_sim_parameters_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to get the pattern of
        !> gridpoints needed when computing fluxes
        !
        !> @date
        !> 27_01_2015 - initial version - J.L. Desmarais
        !
        !> @param operator_type
        !> type of operator used
        !
        !> @return
        !> gridpoints needed around the central gridpoint to compute
        !> the fluxes
        !--------------------------------------------------------------
        function get_sd_pattern_flux_default(operator_type) result(pattern)

          implicit none

          integer    , intent(in) :: operator_type
          integer, dimension(2,2) :: pattern

          integer :: operator_type_s

          print '(''pmodel_eq_default'')'
          print '(''get_sd_pattern_flux_default'')'
          stop 'get_sd_pattern not implemented'

          operator_type_s = operator_type
          pattern(1,1)    = 0

        end function get_sd_pattern_flux_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default subroutine to apply the initial conditions
        !> to the computational domain
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
        subroutine apply_ic_default(this,nodes,x_map,y_map)

          implicit none
          
          class(pmodel_eq_default)     , intent(in)    :: this
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map

          integer     :: nb_eq
          real(rkind) :: node_s
          real(rkind) :: x_s
          real(rkind) :: y_s

          print '(''pmodel_eq_default'')'
          print '(''apply_ic_default'')'
          stop 'not implemented'

          nb_eq  = this%get_eq_nb()
          node_s = nodes(1,1,1)
          x_s    = x_map(1)
          y_s    = y_map(1)

        end subroutine apply_ic_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to compute the fluxes along the
        !> x-axis and distinguishing the inviscid,
        !> viscid and capillary contributions
        !
        !> @date
        !> 10_11_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
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
        !>@return
        !> fluxes along the x-axis
        !--------------------------------------------------------------
        function compute_flux_x_by_parts_default(
     $     nodes,dx,dy,i,j,s_oneside,
     $     inviscid_flux, viscid_flux)
     $     result(flux_x)
        
          implicit none
        
          real(rkind), dimension(:,:,:), intent(in)   :: nodes
          real(rkind)                  , intent(in)   :: dx
          real(rkind)                  , intent(in)   :: dy
          integer(ikind)               , intent(in)   :: i
          integer(ikind)               , intent(in)   :: j
          class(sd_operators)          , intent(in)   :: s_oneside
          real(rkind), dimension(ne)   , intent(out)  :: inviscid_flux
          real(rkind), dimension(ne)   , intent(out)  :: viscid_flux
          real(rkind), dimension(ne)                  :: flux_x


          real(rkind)    :: node_s
          real(rkind)    :: ds
          integer(ikind) :: bc_s

          print '(''pmodel_eq_default_class'')'
          print '(''compute_flux_x_by_parts'')'
          stop 'not implemented'

          node_s = nodes(i,j,1)
          ds     = dx+dy
          bc_s   = s_oneside%get_bc_size()
          inviscid_flux(1) = nodes(i,j,1)
          viscid_flux(1)   = nodes(i,j,1)
          flux_x(1)        = nodes(i,j,1)
          
        end function compute_flux_x_by_parts_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to compute the fluxes along the
        !> y-axis and distinguishing the inviscid,
        !> viscid and capillary contributions
        !
        !> @date
        !> 10_11_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@param i
        !> x-index where the flux_y is computed
        !
        !>@param j
        !> y-index where the flux_y is computed
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
        !> fluxes along the y-axis
        !--------------------------------------------------------------
        function compute_flux_y_by_parts_default(
     $     nodes,dx,dy,i,j,s_oneside,
     $     inviscid_flux, viscid_flux)
     $     result(flux_y)
        
          implicit none
        
          real(rkind), dimension(:,:,:), intent(in)   :: nodes
          real(rkind)                  , intent(in)   :: dx
          real(rkind)                  , intent(in)   :: dy
          integer(ikind)               , intent(in)   :: i
          integer(ikind)               , intent(in)   :: j
          class(sd_operators)          , intent(in)   :: s_oneside
          real(rkind), dimension(ne)   , intent(out)  :: inviscid_flux
          real(rkind), dimension(ne)   , intent(out)  :: viscid_flux
          real(rkind), dimension(ne)                  :: flux_y
          
          real(rkind)    :: node_s
          real(rkind)    :: ds
          integer(ikind) :: bc_s

          print '(''pmodel_eq_default_class'')'
          print '(''compute_flux_y_by_parts'')'
          stop 'not implemented'

          node_s = nodes(i,j,1)
          ds     = dx+dy
          bc_s   = s_oneside%get_bc_size()
          inviscid_flux(1) = nodes(i,j,1)
          viscid_flux(1)   = nodes(i,j,1)
          flux_y(1)        = nodes(i,j,1)        
          
        end function compute_flux_y_by_parts_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to compute the velocity vector
        !> from the conservative variables
        !
        !> @date
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variable vector \f$ (\rho,q_x,q_y,\rho E) \f$
        !
        !>@return
        !> velocity vector \f$ (v_x,v_y) \f$
        !--------------------------------------------------------------
        function get_velocity_default(nodes) result(velocity)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(2)              :: velocity

          print '(''pmodel_eq_default_class'')'
          print '(''get_velocity'')'
          stop 'not implemented'   

          velocity(1) = nodes(1)

        end function get_velocity_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to get the viscous constant ratio
        !
        !> @date
        !> 11_11_2014 - initial version - J.L. Desmarais
        !
        !>@return
        !> viscous coefficient
        !-------------------------------------------------------------
        function get_viscous_coeff_default() result(viscous_coeff)

          implicit none

          real(rkind) :: viscous_coeff

          print '(''pmodel_eq_default_class'')'
          print '(''get_viscous_coeff'')'
          stop 'not implemented'

          viscous_coeff = 0.0d0

        end function get_viscous_coeff_default

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to check whether the open boundary
        !> conditions are undermined at the grid point location
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
        !>@return
        !> check if the open boundary conditions are undermined
        !> at the grid point location
        !--------------------------------------------------------------
        function are_openbc_undermined_default(x_map,y_map,nodes) result(undermined)

          implicit none

          real(rkind), dimension(3)     , intent(in) :: x_map
          real(rkind), dimension(3)     , intent(in) :: y_map
          real(rkind), dimension(3,3,ne), intent(in) :: nodes
          logical                                    :: undermined

          real(rkind) :: x_s,y_s,node_s

          print '(''pmodel_eq_default_class'')'
          print '(''are_openbc_undermined_default'')'
          stop 'not implemented'
          
          x_s        = x_map(1)
          y_s        = y_map(1)
          node_s     = nodes(1,1,1)
          undermined = .false.

        end function are_openbc_undermined_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to compute the governing variables
        !> in the far field as imposed by the initial conditions
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
        !>@return
        !> governing variables in the far field
        !> (ex: \f$ (\rho_\infty,{q_x}_\infty,{q_y}_\infty, {\rho E}_\infty)\f$)
        !--------------------------------------------------------------
        function get_far_field_default(this,t,x,y) result(var)

          implicit none

          class(pmodel_eq_default)  , intent(in) :: this
          real(rkind)               , intent(in) :: t
          real(rkind)               , intent(in) :: x
          real(rkind)               , intent(in) :: y
          real(rkind), dimension(ne)             :: var

          real(rkind) :: viscous_coeff
          real(rkind) :: t_s
          real(rkind) :: x_s
          real(rkind) :: y_s

          print '(''pmodel_eq_default_class'')'
          print '(''get_far_field'')'
          stop 'not implemented'

          viscous_coeff = this%get_viscous_coeff()
          t_s = t
          x_s = x
          y_s = y

          var(1) = 0.0

        end function get_far_field_default
        
      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to determine the vector of primitive
        !> variables used to evaluate the eigenquantities
        !> at the edge of the computational domain
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
        !>@return
        !> extended primitive variable vector to evaluate the
        !> eigenquantities at the boundary \f$ (\rho,v_x,v_y,P,c) \f$
        !--------------------------------------------------------------
        function get_prim_obc_eigenqties_default(this,t,x,y,nodes_bc)
     $    result(nodes_prim_extended)

          implicit none

          class(pmodel_eq_default)    , intent(in) :: this
          real(rkind)                 , intent(in) :: t
          real(rkind)                 , intent(in) :: x
          real(rkind)                 , intent(in) :: y
          real(rkind), dimension(ne)  , intent(in) :: nodes_bc
          real(rkind), dimension(ne+1)             :: nodes_prim_extended

          integer     :: nb_eq
          real(rkind) :: t_s,x_s,y_s

          print '(''pmodel_eq_default_class'')'
          print '(''get_prim_obc_eigenqties_default'')'
          stop 'not implemented'

          nb_eq = this%get_eq_nb()
          t_s = t
          x_s = x
          y_s = y
          nodes_prim_extended(1) = nodes_bc(1)          

        end function get_prim_obc_eigenqties_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to compute the eigenquantity
        !> vectors of the hyperbolic terms using
        !> conservative variables
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes_in
        !> array with the conservative variables
        !
        !>@return
        !> eigenvector
        !--------------------------------------------------------------
        function compute_openbc_var_default(nodes_in) result(nodes_out)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes_in
          real(rkind), dimension(ne)             :: nodes_out
          
          print '(''pmodel_eq_default'')'
          print '(''compute_var_default'')'
          stop 'not implemented'

          nodes_out = nodes_in

        end function compute_openbc_var_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to compute the
        !> eigenmatrix of the hyperbolic terms
        !> using primitive variables
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes_prim_extended
        !> vector with the primitive variables
        !
        !>@return
        !> eigenmatrix
        !--------------------------------------------------------------
        function compute_openbc_matrix_default(nodes_prim_extended) result(matrix)

          implicit none

          real(rkind), dimension(ne+1) , intent(in) :: nodes_prim_extended
          real(rkind), dimension(ne,ne)             :: matrix

          print '(''pmodel_eq_default_class'')'
          print '(''compute_openbc_matrix_default'')'
          stop 'not implemented'

          matrix(1,1) = nodes_prim_extended(1)

        end function compute_openbc_matrix_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to compute the
        !> eigenvector of the hyperbolic terms
        !> using primitive variables
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes_prim_extended
        !> vector with the primitive variables
        !
        !>@return
        !> eigenvector
        !--------------------------------------------------------------
        function compute_openbc_vector_default(nodes_prim_extended) result(vector)

          implicit none

          real(rkind), dimension(ne+1), intent(in) :: nodes_prim_extended
          real(rkind), dimension(ne)               :: vector

          print '(''pmodel_eq_default_class'')'
          print '(''compute_openbc_vector_default'')'
          stop 'not implemented'

          vector(1) = nodes_prim_extended(1)

        end function compute_openbc_vector_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to compute the gradient
        !> of the primitive variables
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
        !>@return
        !> gradient of the primitive variables along the x-axis
        !--------------------------------------------------------------
        function compute_gradient_prim_default(nodes,i,j,gradient,dn,use_n_dir)
     $    result(grad_var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(gradient_proc)                  :: gradient
          real(rkind)                  , intent(in) :: dn
          logical    , optional        , intent(in) :: use_n_dir
          real(rkind), dimension(ne)                :: grad_var


          print '(''pmodel_eq_default_class'')'
          print '(''compute_gradient_prim_default'')'
          stop 'not implemented'

          if(present(use_n_dir)) then
             grad_var = gradient(nodes,i,j,basic,dn)
          else
             grad_var(1) = nodes(1,1,1)
          end if

        end function compute_gradient_prim_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default basic function to extract the first component
        !> of the data at a grid point location
        !
        !> @date
        !> 11_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index identifying the grid point location in the x-direction
        !
        !>@param j
        !> index identifying the grid point location in the y-direction
        !
        !>@return
        !> first component of the data at the grid point
        !--------------------------------------------------------------
        function basic(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var
          
          var = nodes(i,j,1)
          
        end function basic


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to compute the variables
        !> in the (n1,n2) frame from the (x,y) frame
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables vector in the (x,y) frame
        !
        !>@return
        !> governing variables vector in the (n1,n2) frame
        !--------------------------------------------------------------
        function compute_xy_to_n_var_default(nodes) result(nodes_n)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne)             :: nodes_n

          print '(''pmodel_eq_default'')'
          print '(''compute_xy_to_n_var'')'
          print '(''not implemented'')'
          stop ''
          
          nodes_n(1) = nodes(1)

        end function compute_xy_to_n_var_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default function to compute the variables
        !> in the (x,y) frame from the (n1,n2) frame
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes_n
        !> governing variables vector in the (n1,n2) frame
        !
        !>@return
        !> governing variables vector in the (x,y) frame
        !--------------------------------------------------------------
        function compute_n_to_xy_var_default(nodes_n) result(nodes)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes_n
          real(rkind), dimension(ne)             :: nodes

          print '(''pmodel_eq_default'')'
          print '(''compute_n_to_xy_var'')'
          print '(''not implemented'')'
          stop ''
          
          nodes(1) = nodes_n(1)

        end function compute_n_to_xy_var_default

      end module pmodel_eq_default_class

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
        !> abstract class encapsulating default operators to
        !> compute the governing equations of the physical model
        !---------------------------------------------------------------
        type, abstract, extends(pmodel_eq_abstract) :: pmodel_eq_default
          
          contains

          procedure, nopass :: get_sim_parameters => get_sim_parameters_default

          procedure, nopass :: get_sd_pattern_flux_x => get_sd_pattern_flux_default
          procedure, nopass :: get_sd_pattern_flux_y => get_sd_pattern_flux_default

          procedure, nopass :: get_viscous_coeff

          procedure, nopass :: compute_flux_x_by_parts
          procedure, nopass :: compute_flux_y_by_parts

          procedure,   pass :: get_nodes_obc_eigenqties

c$$$          procedure, nopass :: compute_x_eigenvalues  => compute_eigenvalues_default
c$$$          procedure, nopass :: compute_y_eigenvalues  => compute_eigenvalues_default
c$$$          procedure, nopass :: compute_n1_eigenvalues => compute_eigenvalues_default
c$$$          procedure, nopass :: compute_n2_eigenvalues => compute_eigenvalues_default
c$$$
c$$$          procedure, nopass :: compute_x_lefteigenvector    => compute_eigenvector_default
c$$$          procedure, nopass :: compute_x_righteigenvector   => compute_eigenvector_default
c$$$          procedure, nopass :: compute_y_lefteigenvector    => compute_eigenvector_default
c$$$          procedure, nopass :: compute_y_righteigenvector   => compute_eigenvector_default
c$$$          procedure, nopass :: compute_n1_lefteigenvector   => compute_eigenvector_default
c$$$          procedure, nopass :: compute_n1_righteigenvector  => compute_eigenvector_default
c$$$          procedure, nopass :: compute_n2_lefteigenvector   => compute_eigenvector_default
c$$$          procedure, nopass :: compute_n2_righteigenvector  => compute_eigenvector_default
c$$$
c$$$          procedure, nopass :: compute_x_transM   => compute_eigenvector_default
c$$$          procedure, nopass :: compute_y_transM   => compute_eigenvector_default
c$$$          procedure, nopass :: compute_n1_transM  => compute_eigenvector_default
c$$$          procedure, nopass :: compute_n2_transM  => compute_eigenvector_default
c$$$
c$$$          procedure, nopass :: compute_x_leftConsLodiM  => compute_eigenvector_default
c$$$          procedure, nopass :: compute_y_leftConsLodiM  => compute_eigenvector_default
c$$$          procedure, nopass :: compute_x_timedev_from_LODI_vector => compute_timedev_from_LODI_vector_default
c$$$          procedure, nopass :: compute_y_timedev_from_LODI_vector => compute_timedev_from_LODI_vector_default
c$$$          procedure, nopass :: compute_timedev_from_LODI_vectors

          procedure,   pass :: get_far_field => get_far_field_default

          procedure, nopass :: compute_xy_to_n_var
          procedure, nopass :: compute_n_to_xy_var

        end type pmodel_eq_default


        contains


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
        subroutine get_sim_parameters_default(param_name, param_value)

          implicit none

          character(20), dimension(:), allocatable, intent(out) :: param_name
          real(rkind)  , dimension(:), allocatable, intent(out) :: param_value


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
        !> interface to get the pattern of gridpoints needed
        !> when computing fluxes
        !
        !> @date
        !> 27_01_2015 - initial version - J.L. Desmarais
        !
        !>@param operator_type
        !> type of operator used
        !
        !> @return pattern
        !> gridpoints needed around the central gridpoint to compute
        !> the fluxes
        !--------------------------------------------------------------
        function get_sd_pattern_flux_default(operator_type) result(pattern)

          implicit none

          integer    , intent(in) :: operator_type
          integer, dimension(2,2) :: pattern

          integer :: operator_type_s

          print '(''get_sd_pattern_flux_default'')'
          stop 'get_sd_pattern not implemented'

          operator_type_s = operator_type
          pattern(1,1)    = 0

        end function get_sd_pattern_flux_default



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
        function get_viscous_coeff() result(viscous_coeff)

          implicit none

          real(rkind) :: viscous_coeff

          print '(''pmodel_eq_default_class'')'
          print '(''get_viscous_coeff'')'
          stop 'not implemented'

          viscous_coeff = 0.0d0

        end function get_viscous_coeff


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the fluxes by parts (get the inviscid and the
        !> viscid parts)
        !
        !> @date
        !> 10_11_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param i
        !> index identifying the nodes along the x-direction
        !
        !>@param j
        !> index identifying the nodes along the y-direction
        !
        !>@param s_oneside
        !> space discretization operator
        !
        !>@param inviscid_flux
        !> inviscid flux at (i+1/2,j)
        !
        !>@param viscid_flux
        !> viscous flux computed at (i+1/2,j)
        !
        !>@return flux_x
        !> flux computed at (i+1/2,j)
        !--------------------------------------------------------------
        function compute_flux_x_by_parts(
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

          print '(''pmodel_eq_abstract_class'')'
          print '(''compute_flux_x_by_parts'')'
          stop 'not implemented'

          node_s = nodes(i,j,1)
          ds     = dx+dy
          bc_s   = s_oneside%get_bc_size()
          inviscid_flux(1) = nodes(i,j,1)
          viscid_flux(1)   = nodes(i,j,1)
          flux_x(1)        = nodes(i,j,1)
          
        end function compute_flux_x_by_parts


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the fluxes by parts (get the inviscid and the
        !> viscid parts)
        !
        !> @date
        !> 10_11_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param i
        !> index identifying the nodes along the x-direction
        !
        !>@param j
        !> index identifying the nodes along the y-direction
        !
        !>@param s_oneside
        !> space discretization operator
        !
        !>@param inviscid_flux
        !> inviscid flux at (i+1/2,j)
        !
        !>@param viscid_flux
        !> viscous flux computed at (i+1/2,j)
        !
        !>@return flux_x
        !> flux computed at (i+1/2,j)
        !--------------------------------------------------------------
        function compute_flux_y_by_parts(
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

          print '(''pmodel_eq_abstract_class'')'
          print '(''compute_flux_y_by_parts'')'
          stop 'not implemented'

          node_s = nodes(i,j,1)
          ds     = dx+dy
          bc_s   = s_oneside%get_bc_size()
          inviscid_flux(1) = nodes(i,j,1)
          viscid_flux(1)   = nodes(i,j,1)
          flux_y(1)        = nodes(i,j,1)        
          
        end function compute_flux_y_by_parts


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the grid points used to evaluate
        !> the eigenquantities at the edge of the
        !> computational domain
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
        !>@param nodes_bc
        !> array with the grid point data at the boundary
        !
        !>@param nodes_bc
        !> array with the grid point data at the boundary
        !
        !>@param nodes_eigenqties
        !> grid points used to evaluate the eigenquantities at the
        !> boundary
        !--------------------------------------------------------------
        function get_nodes_obc_eigenqties(this,t,x,y,nodes_bc) result(nodes_eigenqties)

          implicit none

          class(pmodel_eq_default)  , intent(in) :: this
          real(rkind)               , intent(in) :: t
          real(rkind)               , intent(in) :: x
          real(rkind)               , intent(in) :: y
          real(rkind), dimension(ne), intent(in) :: nodes_bc
          real(rkind), dimension(ne)             :: nodes_eigenqties


          print '(''pmodel_eq_default_class'')'
          print '(''get_nodes_obc_eigenqties'')'
          print '(''not implemented'')'
          stop ''
          
          nodes_eigenqties(1) = nodes_bc(1)+t+x+y+real(this%get_eq_nb())

        end function get_nodes_obc_eigenqties


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues at the location of the
        !> grid point
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
        function compute_eigenvalues_default(nodes) result(eigenvalues)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne)             :: eigenvalues

          
          real(rkind) :: node_s


          stop 'compute_eigenvalues: not implemented'

          node_s = nodes(1)
          eigenvalues(1) = nodes(1)

        end function compute_eigenvalues_default


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvectors at the location of the
        !> grid point
        !
        !> @date
        !> 11_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvect at the location of the grid point
        !--------------------------------------------------------------
        function compute_eigenvector_default(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect

          
          real(rkind) :: node_s


          stop 'compute_eigenvector: not implemented'

          node_s = nodes(1)
          eigenvect(1,1) = nodes(1)

        end function compute_eigenvector_default


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


        function compute_timedev_from_LODI_vector_default(
     $     nodes,lodi) result(timedev)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne), intent(in) :: lodi
          real(rkind), dimension(ne)             :: timedev


          print '(''pmodel_eq_default_class'')'
          print '(''compute_timedev_from_LODI_vector_default'')'
          print '(''not implemented'')'
          stop ''
          
          timedev(1) = nodes(1)+lodi(1)

        end function compute_timedev_from_LODI_vector_default


        function compute_timedev_from_LODI_vectors(
     $     nodes,lodi_x,lodi_y) result(timedev)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne), intent(in) :: lodi_x
          real(rkind), dimension(ne), intent(in) :: lodi_y
          real(rkind), dimension(ne)             :: timedev


          print '(''pmodel_eq_default_class'')'
          print '(''compute_timedev_from_LODI_vectors'')'
          print '(''not implemented'')'
          stop ''
          
          timedev(1) = nodes(1)+lodi_x(1)+lodi_y(1)

        end function compute_timedev_from_LODI_vectors



c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> default computation of the gradient of the
c$$$        !> governing variables in a diagonal direction
c$$$        !
c$$$        !> @date
c$$$        !> 11_08_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param nodes
c$$$        !> array with the grid point data
c$$$        !
c$$$        !>@param i
c$$$        !> integer identifying the index in the x-direction
c$$$        !
c$$$        !>@param j
c$$$        !> integer identifying the index in the y-direction
c$$$        !
c$$$        !>@param gradient
c$$$        !> procedure used to compute the gradient along the diagonal
c$$$        !> direction
c$$$        !
c$$$        !>@param dx
c$$$        !> grid space step along the x-axis
c$$$        !
c$$$        !>@param dy
c$$$        !> grid space step along the y-axis
c$$$        !
c$$$        !>@return grad_var
c$$$        !> gradient of the governing variables along the x-axis
c$$$        !--------------------------------------------------------------
c$$$        function compute_gradient_default(nodes,i,j,gradient,dx,dy) result(grad_var)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind), dimension(:,:,:), intent(in) :: nodes
c$$$          integer(ikind)               , intent(in) :: i
c$$$          integer(ikind)               , intent(in) :: j
c$$$          procedure(gradient_n_proc)                :: gradient
c$$$          real(rkind)                  , intent(in) :: dx
c$$$          real(rkind)                  , intent(in) :: dy
c$$$          real(rkind), dimension(ne)                :: grad_var
c$$$          
c$$$
c$$$          stop '(''compute_n_gradient: not implemented'')'
c$$$          
c$$$          grad_var = gradient(nodes,i,j,basic,dx,dy)
c$$$          
c$$$        end function compute_gradient_default


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
        !>@return var
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



        function compute_xy_to_n_var(nodes) result(nodes_n)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne)             :: nodes_n

          print '(''pmodel_eq_default'')'
          print '(''compute_xy_to_n_var'')'
          print '(''not implemented'')'
          stop ''
          
          nodes_n(1) = nodes(1)

        end function compute_xy_to_n_var


        function compute_n_to_xy_var(nodes_n) result(nodes)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes_n
          real(rkind), dimension(ne)             :: nodes

          print '(''pmodel_eq_default'')'
          print '(''compute_n_to_xy_var'')'
          print '(''not implemented'')'
          stop ''
          
          nodes(1) = nodes_n(1)

        end function compute_n_to_xy_var

      end module pmodel_eq_default_class

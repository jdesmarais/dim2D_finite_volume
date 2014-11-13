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
     $       gradient_n_proc

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
          procedure, nopass :: get_viscous_coeff

          procedure, nopass :: compute_flux_x_by_parts
          procedure, nopass :: compute_flux_y_by_parts

          procedure, nopass :: compute_x_eigenvalues  => compute_eigenvalues_default
          procedure, nopass :: compute_y_eigenvalues  => compute_eigenvalues_default
          procedure, nopass :: compute_n1_eigenvalues => compute_eigenvalues_default
          procedure, nopass :: compute_n2_eigenvalues => compute_eigenvalues_default

          procedure, nopass :: compute_x_lefteigenvector    => compute_eigenvector_default
          procedure, nopass :: compute_x_righteigenvector   => compute_eigenvector_default
          procedure, nopass :: compute_y_lefteigenvector    => compute_eigenvector_default
          procedure, nopass :: compute_y_righteigenvector   => compute_eigenvector_default
          procedure, nopass :: compute_n1_lefteigenvector   => compute_eigenvector_default
          procedure, nopass :: compute_n1_righteigenvector  => compute_eigenvector_default
          procedure, nopass :: compute_n2_lefteigenvector   => compute_eigenvector_default
          procedure, nopass :: compute_n2_righteigenvector  => compute_eigenvector_default

          procedure, nopass :: compute_cons_lodi_matrix_x   => compute_eigenvector_default
          procedure, nopass :: compute_cons_lodi_matrix_y   => compute_eigenvector_default

          procedure, nopass :: compute_x_transM  => compute_eigenvector_default
          procedure, nopass :: compute_y_transM  => compute_eigenvector_default

          procedure, nopass :: compute_n_gradient => compute_gradient_default

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

          character(10), dimension(:), allocatable, intent(out) :: param_name
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default computation of the gradient of the
        !> governing variables in a diagonal direction
        !
        !> @date
        !> 11_08_2014 - initial version - J.L. Desmarais
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
        function compute_gradient_default(nodes,i,j,gradient,dx,dy) result(grad_var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(gradient_n_proc)                :: gradient
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind), dimension(ne)                :: grad_var
          

          stop '(''compute_n_gradient: not implemented'')'
          
          grad_var = gradient(nodes,i,j,basic,dx,dy)
          
        end function compute_gradient_default


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

      end module pmodel_eq_default_class

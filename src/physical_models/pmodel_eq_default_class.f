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
      
        use interface_primary       , only : gradient_n_proc
        use parameters_input        , only : nx,ny,ne
        use parameters_kind         , only : ikind, rkind
        use pmodel_eq_abstract_class, only : pmodel_eq_abstract

        implicit none

        private
        public :: pmodel_eq_default


        !> @class pmodel_eq_default
        !> abstract class encapsulating default operators to
        !> compute the governing equations of the physical model
        !---------------------------------------------------------------
        type, abstract, extends(pmodel_eq_abstract) :: pmodel_eq_default
          
          contains

          procedure, nopass :: get_sim_parameters => get_sim_parameters_default

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

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          procedure(gradient_n_proc)                   :: gradient
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          real(rkind), dimension(ne)                   :: grad_var
          

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

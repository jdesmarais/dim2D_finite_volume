      !> @file
      !> module implemeting subroutines to apply open boundary
      !> conditions at the edge of the computational domain using
      !> hedstrom boundary conditions with special treatment for
      !> the corners
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implemeting subroutines to apply open boundary
      !> conditions at the edge of the computational domain using
      !> hedstrom boundary conditions with special treatment for
      !> the corners
      !
      !> @date
      !> 06_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
       module hedstrom_ncoords_module

        use interface_primary, only :
     $       gradient_n_proc

        use openbc_operators_module, only :
     $       incoming_proc,
     $       incoming_left,
     $       incoming_right,
     $       add_body_forces

        use parameters_constant, only :
     $       x_direction,
     $       y_direction

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       rkind,ikind

        use pmodel_eq_class, only :
     $       pmodel_eq


        implicit none


        private
        public ::
     $       compute_timedev_corner_ncoords,
     $       compute_n_timedev_with_openbc


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for the
        !> SW corner using open boundary conditions
        !
        !> @date
        !> 06_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array of data for the grid points
        !
        !>@param j
        !> index identifying the grid point in the y-direction
        !
        !>@param dx
        !> space step along the x-axis
        !
        !>@param dy
        !> space step along the y-axis
        !
        !>@param p_model
        !> governing equations of the physical model
        !
        !>@param gradient_y
        !> procedure computing the gradient along the y-direction
        !
        !>@param incoming_y
        !> procedure checking the type of characteristic at the edge
        !> in the y-direction
        !
        !>@param timedev
        !> time derivatives modified
        !-------------------------------------------------------------
        subroutine compute_timedev_corner_ncoords(
     $       nodes, i_indices, j, dx, dy, p_model,
     $       gradient1, gradient2,
     $       incoming_wave,
     $       direction,
     $       timedev)

          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          integer(ikind), dimension(2)      , intent(in)    :: i_indices
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          procedure(gradient_n_proc)                        :: gradient1
          procedure(gradient_n_proc)                        :: gradient2
          procedure(incoming_proc)                          :: incoming_wave
          integer                           , intent(in)    :: direction
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          integer(ikind) :: i


          i=i_indices(1)
          timedev(i,j,:) = 
     $         compute_n_timedev_with_openbc(
     $         nodes, i,j, p_model, dx,dy,
     $         gradient1, incoming_wave,direction) +
     $         
     $         add_body_forces(p_model, nodes(i,j,:))

          i=i_indices(2)
          timedev(i,j,:) = 
     $         compute_n_timedev_with_openbc(
     $         nodes, i,j, p_model, dx,dy,
     $         gradient2, incoming_wave,direction) +
     $         
     $         add_body_forces(p_model, nodes(i,j,:))

        end subroutine compute_timedev_corner_ncoords


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the contribution of the hyperbolic terms of the
        !> governing equations in the n-direction to the time
        !> derivatives of the governing variables
        !
        !> @date
        !> 06_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array of grid points
        !
        !>@param i
        !> index identifying the grid point position along the x-axis
        !
        !>@param j
        !> index identifying the grid point position along the y-axis
        !
        !>@param p_model
        !> governing equations of the physical model
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param gradient
        !> procedure computing the gradient along the n-direction
        !
        !>@param incoming_wave
        !> procedure identifying whether the wave is incoming or
        !> outgoing the edge of the computational domain
        !-------------------------------------------------------------
        function compute_n_timedev_with_openbc(
     $     nodes, i, j, p_model, dx, dy,
     $     gradient, incoming_wave, direction)
     $     result(timedev)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          type(pmodel_eq)                 , intent(in) :: p_model
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          procedure(gradient_n_proc)                   :: gradient
          procedure(incoming_proc)                     :: incoming_wave
          integer                         , intent(in) :: direction
          real(rkind), dimension(ne)                   :: timedev


          real(rkind), dimension(ne)    :: eigenvalues
          real(rkind), dimension(ne,ne) :: left_eigenmatrix
          real(rkind), dimension(ne,ne) :: right_eigenmatrix
          integer                       :: k
          real(rkind), dimension(ne)    :: incoming_amp
          real(rkind), dimension(ne)    :: var_gradient


          !determination of the speed of the amplitude waves
          !along the n-direction
          select case(direction)

            case(x_direction)
               eigenvalues       = p_model%compute_n1_eigenvalues(nodes(i,j,:))
               left_eigenmatrix  = p_model%compute_n1_lefteigenvector(nodes(i,j,:))
               right_eigenmatrix = p_model%compute_n1_righteigenvector(nodes(i,j,:))

            case(y_direction)
               eigenvalues       = p_model%compute_n2_eigenvalues(nodes(i,j,:))
               left_eigenmatrix  = p_model%compute_n2_lefteigenvector(nodes(i,j,:))
               right_eigenmatrix = p_model%compute_n2_righteigenvector(nodes(i,j,:))

            case default
               print '(''hedstrom_xy_module'')'
               print '(''compute_n_timedev_with_openbc'')'
               stop 'direction not recognized'

          end select


          !construction of the vector of characteristic amplitudes
          do k=1, ne

             !distinction of the characteristic waves between the
             !incoming and outgoing
             !if the wave is incoming, its amplitude is set to zero
             if(incoming_wave(eigenvalues(k))) then
                if(rkind.eq.8) then
                   incoming_amp(k) = 0.0d0
                else
                   incoming_amp(k) = 0.0
                end if

             !otherwise, the characteristic amplitude is computed using
             !one-side differentiation
             else

                var_gradient     = p_model%compute_n_gradient(
     $               nodes, i,j, gradient, dx,dy)   
             
                incoming_amp(k)  =  - eigenvalues(k)*
     $               DOT_PRODUCT(left_eigenmatrix(:,k), var_gradient)

             end if

          end do


          !determination of the contribution of the hyperbolic terms
          !to the time derivatives
          timedev = MATMUL(incoming_amp, right_eigenmatrix)

        end function compute_n_timedev_with_openbc

      end module hedstrom_ncoords_module

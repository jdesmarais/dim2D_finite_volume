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
       module hedstrom_xy_corners_module

        use interface_primary, only :
     $       gradient_n_proc

        use openbc_operators_module, only :
     $       incoming_proc

        use parameters_constant, only :
     $       n1_direction,
     $       n2_direction

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       bc_size

        use parameters_kind, only :
     $       rkind,
     $       ikind

        use pmodel_eq_class, only :
     $       pmodel_eq


        implicit none


        private
        public ::
     $       compute_n_timedev_with_openbc


        contains


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
     $     gradient_dir1, gradient_dir2,
     $     incoming_wave,
     $     dir1)
     $     result(timedev)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          type(pmodel_eq)                 , intent(in) :: p_model
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          procedure(gradient_n_proc)                   :: gradient_dir1
          procedure(gradient_n_proc)                   :: gradient_dir2
          procedure(incoming_proc)                     :: incoming_wave
          integer                         , intent(in) :: dir1
          real(rkind), dimension(ne)                   :: timedev

          real(rkind), dimension(ne)    :: eigenvalues
          real(rkind), dimension(ne,ne) :: left_eigenmatrix
          real(rkind), dimension(ne,ne) :: right_eigenmatrix
          integer                       :: k
          real(rkind), dimension(ne)    :: incoming_amp
          real(rkind), dimension(ne)    :: var_gradient_dir1

          real(rkind), dimension(ne,ne) :: transM_dir1
          real(rkind), dimension(ne)    :: var_gradient_dir2

          real(rkind), dimension(ne)    :: normal_timedev
          real(rkind), dimension(ne)    :: trans_timedev
                                                        


          !determination of the speed of the amplitude waves
          !along the n-direction
          select case(dir1)

            case(n1_direction)
               eigenvalues       = p_model%compute_n1_eigenvalues(nodes(i,j,:))
               left_eigenmatrix  = p_model%compute_n1_lefteigenvector(nodes(i,j,:))
               right_eigenmatrix = p_model%compute_n1_righteigenvector(nodes(i,j,:))
               transM_dir1       = p_model%compute_n1_transM(nodes(i,j,:))

            case(n2_direction)
               eigenvalues       = p_model%compute_n2_eigenvalues(nodes(i,j,:))
               left_eigenmatrix  = p_model%compute_n2_lefteigenvector(nodes(i,j,:))
               right_eigenmatrix = p_model%compute_n2_righteigenvector(nodes(i,j,:))
               transM_dir1       = p_model%compute_n2_transM(nodes(i,j,:))

            case default
               print '(''hedstrom_xy_corners_module'')'
               print '(''compute_n_timedev_with_openbc'')'
               stop 'direction not recognized'

          end select


          !construction of the vector of characteristic amplitudes
          !in the direction in which the wave are leaving the domain
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

                var_gradient_dir1     = p_model%compute_n_gradient(
     $               nodes, i,j, gradient_dir1, dx,dy)
             
                incoming_amp(k)  =  - eigenvalues(k)*
     $               DOT_PRODUCT(left_eigenmatrix(:,k), var_gradient_dir1)

             end if

          end do


          !in the direction transverse to leaving the domain, the 
          !contribution of the transverse terms to the time derivatives
          !is computed using the transverse matrix and the transverse
          !gradient
          var_gradient_dir2 = p_model%compute_n_gradient(
     $         nodes, i,j, gradient_dir2, dx,dy)


          !determination of the contribution of the hyperbolic terms
          !to the time derivatives
          normal_timedev = MATMUL(incoming_amp, right_eigenmatrix)
          trans_timedev  = MATMUL(var_gradient_dir2, transM_dir1)
          timedev = normal_timedev - trans_timedev

        end function compute_n_timedev_with_openbc

      end module hedstrom_xy_corners_module

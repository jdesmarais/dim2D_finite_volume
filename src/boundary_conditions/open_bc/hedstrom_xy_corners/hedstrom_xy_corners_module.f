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
     $       gradient_proc

        use openbc_operators_module, only :
     $       incoming_proc

        use parameters_constant, only :
     $       n1_direction,
     $       n2_direction

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none


        private
        public ::
     $       compute_timedev_with_openbc


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
        !>@param gradient_x
        !> procedure computing the gradient along the x-direction
        !
        !>@param gradient_y
        !> procedure computing the gradient along the y-direction
        !
        !>@param incoming_wave
        !> procedure identifying whether the wave is incoming or
        !> outgoing the edge of the computational domain
        !
        !>@param dir1
        !> direction in which the vector points to the outside of
        !> the computational domain
        !-------------------------------------------------------------
        function compute_timedev_with_openbc(
     $       t,x,y,
     $       nodes, i,j,
     $       p_model,
     $       direction,
     $       gradient_n, dn,
     $       incoming_wave)
     $       result(timedev_n)

          implicit none

          real(rkind)                  , intent(in) :: t
          real(rkind)                  , intent(in) :: x
          real(rkind)                  , intent(in) :: y
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          type(pmodel_eq)              , intent(in) :: p_model
          integer                      , intent(in) :: direction
          procedure(gradient_proc)                  :: gradient_n
          real(rkind)                  , intent(in) :: dn
          procedure(incoming_proc)                  :: incoming_wave          
          real(rkind), dimension(ne)                :: timedev_n


          real(rkind), dimension(ne+1)  :: obc_prim_var
          real(rkind), dimension(ne+1)  :: obc_n_prim_var
          real(rkind), dimension(ne)    :: eigenvalues
          real(rkind), dimension(ne,ne) :: left_eigenmatrix
          real(rkind), dimension(ne,ne) :: right_eigenmatrix

          real(rkind), dimension(ne)    :: gradient_prim
          integer                       :: k
          real(rkind), dimension(ne)    :: incoming_amp
          real(rkind), dimension(ne)    :: left_eigenvector

          real(rkind), dimension(ne,ne) :: jacConsPrim
          real(rkind), dimension(ne)    :: timedev_n_prim
                                                        
          
          !1) determination of the primitive
          !   variables used to compute the
          !   eigenquantities
          !
          !   obc_prim_var   = [rho,ux,uy,P,c]
          !   obc_n_prim_var = [rho,u_n1,u_n2,P,c]
          !----------------------------------------
          obc_prim_var         = p_model%get_prim_obc_eigenqties(t,x,y,nodes(i,j,:))
          obc_n_prim_var(1:ne) = p_model%compute_xy_to_n_var(obc_prim_var(1:ne))
          obc_n_prim_var(ne+1) = obc_prim_var(ne+1)


          !2) determination of the eigenquantities
          !----------------------------------------
          select case(direction)
            case(n1_direction)
               eigenvalues       = p_model%compute_x_eigenvalues_prim(obc_n_prim_var)
               left_eigenmatrix  = p_model%compute_x_lefteigenvector_prim(obc_n_prim_var)
               right_eigenmatrix = p_model%compute_x_righteigenvector_prim(obc_n_prim_var)
               
            case(n2_direction)
               eigenvalues       = p_model%compute_y_eigenvalues_prim(obc_n_prim_var)
               left_eigenmatrix  = p_model%compute_y_lefteigenvector_prim(obc_n_prim_var)
               right_eigenmatrix = p_model%compute_y_righteigenvector_prim(obc_n_prim_var)

            case default
               print '(''hedstrom_xy_corners_module'')'
               print '(''compute_timedev_with_openbc'')'
               print '(''direction not recognized'')'
               stop ''
          end select


          ! determination of the gradient
          !----------------------------------------
          gradient_prim = p_model%compute_gradient_prim(
     $         nodes,i,j,gradient_n,dn,
     $         use_n_dir=.true.)


          !construction of the vector of
          !characteristic amplitudes
          !----------------------------------------
          do k=1, ne

             !distinction of the characteristic waves
             !between the incoming and outgoing. If
             !the wave is incoming, its amplitude is
             !set to zero
             if(incoming_wave(eigenvalues(k))) then

                if(rkind.eq.8) then
                   incoming_amp(k) = 0.0d0
                else
                   incoming_amp(k) = 0.0
                end if

             !otherwise, the characteristic amplitude
             !is computed using one-side differentiation
             else

                left_eigenvector = left_eigenmatrix(:,k)

                incoming_amp(k)  =  - eigenvalues(k)*DOT_PRODUCT(
     $                                   left_eigenvector,
     $                                   gradient_prim)

             end if

          end do


          !determination of the contribution of the hyperbolic terms
          !to the time derivatives
          timedev_n_prim = MATMUL(incoming_amp, right_eigenmatrix)


          !determination of the time derivatives
          !of the conservative variables
          jacConsPrim = p_model%compute_jacobian_cons_to_prim(obc_n_prim_var)
          timedev_n   = MATMUL(timedev_n_prim, jacConsPrim)

        end function compute_timedev_with_openbc

      end module hedstrom_xy_corners_module

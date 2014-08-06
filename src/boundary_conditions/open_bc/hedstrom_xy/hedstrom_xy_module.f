      !> @file
      !> module implemeting subroutines to apply open boundary
      !> conditions at the edge of the computational domain using
      !> hedstrom boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implemeting subroutines to apply open boundary
      !> conditions at the edge of the computational domain using
      !> hedstrom boundary conditions
      !
      !> @date
      !> 04_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module hedstrom_xy_module

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        use openbc_operators_module, only :
     $       incoming_proc,
     $       incoming_left,
     $       incoming_right,
     $       add_body_forces

        use pmodel_eq_class, only :
     $       pmodel_eq        

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       rkind,ikind

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_L1,
     $       gradient_x_x_oneside_R1,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_L1,
     $       gradient_y_y_oneside_R1,
     $       gradient_y_y_oneside_R0
        
        implicit none

        private
        public ::
     $       compute_timedev_xlayer,
     $       compute_timedev_ylayer,
     $       compute_timedev_corner_W,
     $       compute_timedev_corner_E,
     $       compute_x_timedev_with_openbc,
     $       compute_y_timedev_with_openbc

        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for the
        !> south or north layer using open boundary conditions
        !
        !> @date
        !> 04_08_2014 - initial version - J.L. Desmarais
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
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param gradient_x
        !> procedure computing the gradient along the x-direction
        !
        !>@param incoming_x
        !> procedure checking the type of characteristic at the edge
        !> in the x-direction
        !
        !>@param timedev
        !> time derivatives modified
        !-------------------------------------------------------------
        subroutine compute_timedev_xlayer(
     $     nodes, i,j, dx,dy, p_model, flux_y,
     $     gradient_x, incoming_x,
     $     timedev)

          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          integer(ikind)                    , intent(in)    :: i
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx,ny+1,ne), intent(in)    :: flux_y
          procedure(gradient_x_proc)                        :: gradient_x
          procedure(incoming_proc)                          :: incoming_x
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev


          timedev(i,j,:) =
     $            compute_x_timedev_with_openbc(
     $            nodes, i, j, p_model, dx,
     $            gradient_x, incoming_x) +
     $         
     $            1.0d0/dy*(flux_y(i,j,:) - flux_y(i,j+1,:)) +
     $            
     $            add_body_forces(p_model, nodes(i,j,:))

        end subroutine compute_timedev_xlayer



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for the
        !> south or north layer using open boundary conditions
        !
        !> @date
        !> 04_08_2014 - initial version - J.L. Desmarais
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
        !>@param flux_x
        !> fluxes along the x-direction
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
        subroutine compute_timedev_ylayer(
     $     nodes, j, dx, dy, p_model, flux_x,
     $     gradient_y, incoming_y,
     $     timedev)

          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx+1,ny,ne), intent(in)    :: flux_x
          procedure(gradient_y_proc)                        :: gradient_y
          procedure(incoming_proc)                          :: incoming_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          
          integer(ikind) :: i


          do i=3, nx-bc_size
             timedev(i,j,:) =
     $            1.0d0/dx*(flux_x(i,j,:) - flux_x(i+1,j,:)) +
     $         
     $            compute_y_timedev_with_openbc(
     $            nodes, i, j, p_model, dy,
     $            gradient_y, incoming_y) +
     $            
     $            add_body_forces(p_model, nodes(i,j,:))
          end do

        end subroutine compute_timedev_ylayer



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for the
        !> west corner using open boundary conditions
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
        subroutine compute_timedev_corner_W(
     $     nodes, j, dx, dy, p_model,
     $     gradient_y, incoming_y,
     $     timedev)

          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          procedure(gradient_y_proc)                        :: gradient_y
          procedure(incoming_proc)                          :: incoming_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          
          integer(ikind) :: i
          

          i=1
          timedev(i,j,:) = 
     $         compute_x_timedev_with_openbc(
     $         nodes, i, j, p_model, dx,
     $         gradient_x_x_oneside_L0, incoming_left) + 
     $         
     $         compute_y_timedev_with_openbc(
     $         nodes, i, j, p_model, dy,
     $         gradient_y, incoming_y) +
     $         
     $         add_body_forces(p_model, nodes(i,j,:))

          i=2
          timedev(i,j,:) = 
     $         compute_x_timedev_with_openbc(
     $         nodes, i, j, p_model, dx,
     $         gradient_x_x_oneside_L1, incoming_left) + 
     $         
     $         compute_y_timedev_with_openbc(
     $         nodes, i, j, p_model, dy,
     $         gradient_y, incoming_y) +
     $         
     $         add_body_forces(p_model, nodes(i,j,:))

        end subroutine compute_timedev_corner_W        


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for the
        !> east corner using open boundary conditions
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
        subroutine compute_timedev_corner_E(
     $     nodes, j, dx, dy, p_model,
     $     gradient_y, incoming_y,
     $     timedev)

          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          procedure(gradient_y_proc)                        :: gradient_y
          procedure(incoming_proc)                          :: incoming_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          integer(ikind) :: i

          i=nx-1
          timedev(i,j,:) = 
     $         compute_x_timedev_with_openbc(
     $         nodes, i, j, p_model, dx,
     $         gradient_x_x_oneside_R1, incoming_right) + 
     $         
     $         compute_y_timedev_with_openbc(
     $         nodes, i, j, p_model, dy,
     $         gradient_y, incoming_y) +
     $         
     $         add_body_forces(p_model, nodes(i,j,:))

          i=nx
          timedev(i,j,:) = 
     $         compute_x_timedev_with_openbc(
     $         nodes, i, j, p_model, dx,
     $         gradient_x_x_oneside_R0, incoming_right) + 
     $         
     $         compute_y_timedev_with_openbc(
     $         nodes, i, j, p_model, dy,
     $         gradient_y, incoming_y) +
     $         
     $         add_body_forces(p_model, nodes(i,j,:))

        end subroutine compute_timedev_corner_E



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the contribution of the hyperbolic terms of the
        !> governing equations in the x-direction to the time
        !> derivatives of the governing variables
        !
        !> @date
        !> 04_08_2014 - initial version - J.L. Desmarais
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
        !>@param gradient
        !> procedure computing the gradient along the x-direction
        !
        !>@param incoming_wave
        !> procedure identifying whether the wave is incoming or
        !> outgoing the edge of the computational domain
        !-------------------------------------------------------------
        function compute_x_timedev_with_openbc(
     $     nodes, i, j, p_model, dx,
     $     gradient, incoming_wave)
     $     result(timedev)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          type(pmodel_eq)                 , intent(in) :: p_model
          real(rkind)                     , intent(in) :: dx
          procedure(gradient_x_proc)                   :: gradient
          procedure(incoming_proc)                     :: incoming_wave
          real(rkind), dimension(ne)                   :: timedev


          real(rkind), dimension(ne) :: eigenvalues
          integer                    :: k
          real(rkind), dimension(ne) :: incoming_amp
          real(rkind), dimension(ne) :: left_eigenvector
          real(rkind), dimension(ne) :: var_gradient
          real(rkind), dimension(ne) :: right_eigenvector


          !determination of the speed of the amplitude waves
          !along the x-direction
          eigenvalues = p_model%compute_x_eigenvalues(nodes(i,j,:))


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
                left_eigenvector = p_model%compute_x_lefteigenvector(
     $               nodes(i,j,:),k)
                var_gradient     = p_model%compute_x_gradient(
     $               nodes, i,j, gradient, dx)
                incoming_amp(k)  =  - eigenvalues(k)*
     $               DOT_PRODUCT(left_eigenvector, var_gradient)

             end if

          end do


          !determination of the contribution of the hyperbolic terms
          !to the time derivatives
          do k=1,ne

             right_eigenvector = p_model%compute_x_righteigenvector(
     $               nodes(i,j,:),k)
             timedev(k) = DOT_PRODUCT(right_eigenvector, incoming_amp)

          end do

        end function compute_x_timedev_with_openbc


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the contribution of the hyperbolic terms of the
        !> governing equations in the y-direction to the time
        !> derivatives of the governing variables
        !
        !> @date
        !> 04_08_2014 - initial version - J.L. Desmarais
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
        !>@param dy
        !> space step along the y-direction
        !
        !>@param gradient
        !> procedure computing the gradient along the y-direction
        !
        !>@param incoming_wave
        !> procedure identifying whether the wave is incoming or
        !> outgoing the edge of the computational domain
        !-------------------------------------------------------------
        function compute_y_timedev_with_openbc(
     $     nodes, i, j, p_model, dy,
     $     gradient, incoming_wave)
     $     result(timedev)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          type(pmodel_eq)                 , intent(in) :: p_model
          real(rkind)                     , intent(in) :: dy
          procedure(gradient_y_proc)                   :: gradient
          procedure(incoming_proc)                     :: incoming_wave
          real(rkind), dimension(ne)                   :: timedev


          real(rkind), dimension(ne) :: eigenvalues
          integer                    :: k
          real(rkind), dimension(ne) :: incoming_amp
          real(rkind), dimension(ne) :: left_eigenvector
          real(rkind), dimension(ne) :: var_gradient
          real(rkind), dimension(ne) :: right_eigenvector


          !determination of the speed of the amplitude waves
          !along the x-direction
          eigenvalues = p_model%compute_y_eigenvalues(nodes(i,j,:))


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
                left_eigenvector = p_model%compute_y_lefteigenvector(
     $               nodes(i,j,:),k)
                var_gradient     = p_model%compute_y_gradient(
     $               nodes, i,j, gradient, dy)
                incoming_amp(k)  =  - eigenvalues(k)*
     $               DOT_PRODUCT(left_eigenvector, var_gradient)

             end if

          end do


          !determination of the contribution of the hyperbolic terms
          !to the time derivatives
          do k=1,ne

             right_eigenvector = p_model%compute_y_righteigenvector(
     $               nodes(i,j,:),k)
             timedev(k) = DOT_PRODUCT(right_eigenvector, incoming_amp)

          end do

        end function compute_y_timedev_with_openbc        

      end module hedstrom_xy_module

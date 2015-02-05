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
      !> 04_08_2014 - initial version                  - J.L. Desmarais
      !> 21_10_2014 - local versions for buffer layers - J.L. Desmarais
      !-----------------------------------------------------------------
      module hedstrom_xy_module

        use interface_primary, only :
     $       gradient_proc

        use openbc_operators_module, only :
     $       incoming_proc,
     $       incoming_left,
     $       incoming_right,
     $       add_body_forces

        use pmodel_eq_class, only :
     $       pmodel_eq        

        use parameters_constant, only :
     $       left,
     $       right,
     $       x_direction,
     $       y_direction

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       rkind,ikind

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_R0

        
        implicit none

        private
        public ::
     $       compute_timedev_y_layer_interior,
     $       compute_timedev_x_edge_local,
     $       compute_timedev_y_edge_local,
     $       compute_timedev_corner_local,
     $       compute_timedev_with_openbc

        contains

        subroutine compute_timedev_y_layer_interior(
     $       p_model,
     $       t,nodes,x_map,y_map,
     $       flux_x,
     $       timedev,
     $       j,dx,dy,
     $       gradient_y,
     $       incoming_y)

          implicit none

          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx+1,ny,ne), intent(in)    :: flux_x
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          procedure(gradient_proc)                          :: gradient_y
          procedure(incoming_proc)                          :: incoming_y
          
          integer(ikind) :: i


          !W corner
          do i=1,bc_size

             timedev(i,j,:) = compute_timedev_corner_local(
     $            p_model,
     $            t, x_map, y_map, nodes,
     $            dx,dy, i,j,
     $            incoming_left, incoming_y,
     $            gradient_x_x_oneside_L0, gradient_y)
             
          end do


          !y-edge
          do i=bc_size+1,nx-bc_size

             timedev(i,j,:) = compute_timedev_y_edge_local(
     $            p_model,
     $            t,x_map,y_map,nodes,
     $            dx,dy, i,j,
     $            flux_x,
     $            incoming_y,
     $            gradient_y)

          end do

          !E corner
          do i=nx-bc_size+1,nx

             timedev(i,j,:) = compute_timedev_corner_local(
     $            p_model,
     $            t, x_map, y_map, nodes,
     $            dx,dy, i,j,
     $            incoming_right, incoming_y,
     $            gradient_x_x_oneside_R0, gradient_y)

          end do

        end subroutine compute_timedev_y_layer_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for the
        !> west or east layer using open boundary conditions
        !
        !> @date
        !> 21_10_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> governing equations of the physical model
        !        
        !>@param nodes
        !> array of data for the grid points
        !
        !>@param dx
        !> space step along the x-axis
        !
        !>@param dy
        !> space step along the y-axis
        !
        !>@param i
        !> index identifying the grid point in the x-direction
        !
        !>@param j
        !> index identifying the grid point in the y-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param incoming_x
        !> procedure checking the type of characteristic at the edge
        !> in the x-direction
        !
        !>@param gradient_x
        !> procedure computing the gradient along the x-direction
        !
        !>@param timedev
        !> time derivatives modified
        !-------------------------------------------------------------
        function compute_timedev_x_edge_local(
     $     p_model,
     $     t,
     $     x_map,
     $     y_map,
     $     nodes,
     $     dx,
     $     dy,
     $     i,
     $     j,
     $     flux_y,
     $     incoming_x,
     $     gradient_x)
     $     result(timedev)

          implicit none


          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(:)         , intent(in)    :: x_map
          real(rkind), dimension(:)         , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)     , intent(in)    :: nodes
          integer(ikind)                    , intent(in)    :: i
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(:,:,:)     , intent(in)    :: flux_y
          procedure(incoming_proc)                          :: incoming_x
          procedure(gradient_proc)                          :: gradient_x
          real(rkind), dimension(ne)                        :: timedev

          timedev =
     $         compute_timedev_with_openbc(
     $            t,x_map(i),y_map(j),
     $            nodes,i,j,
     $            p_model,
     $            x_direction,
     $            gradient_x, dx,
     $            incoming_x) +
     $         
     $         1.0d0/dy*(flux_y(i,j,:) - flux_y(i,j+1,:)) +
     $         
     $         add_body_forces(
     $            p_model,
     $            t, x_map(i), y_map(j), nodes(i,j,:))

        end function compute_timedev_x_edge_local


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for the
        !> north or south layer using open boundary conditions
        !
        !> @date
        !> 21_10_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> governing equations of the physical model
        !        
        !>@param nodes
        !> array of data for the grid points
        !
        !>@param dx
        !> space step along the x-axis
        !
        !>@param dy
        !> space step along the y-axis
        !
        !>@param i
        !> index identifying the grid point in the x-direction
        !
        !>@param j
        !> index identifying the grid point in the y-direction
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param incoming_y
        !> procedure checking the type of characteristic at the edge
        !> in the y-direction
        !
        !>@param gradient_y
        !> procedure computing the gradient along the y-direction
        !
        !>@param timedev
        !> time derivatives modified
        !-------------------------------------------------------------
        function compute_timedev_y_edge_local(
     $     p_model,
     $     t,
     $     x_map,
     $     y_map,
     $     nodes,
     $     dx,
     $     dy,
     $     i,
     $     j,
     $     flux_x,
     $     incoming_y,
     $     gradient_y)
     $     result(timedev)

          implicit none

          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(:)         , intent(in)    :: x_map
          real(rkind), dimension(:)         , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)     , intent(in)    :: nodes
          integer(ikind)                    , intent(in)    :: i
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(:,:,:)     , intent(in)    :: flux_x
          procedure(incoming_proc)                          :: incoming_y
          procedure(gradient_proc)                          :: gradient_y
          real(rkind), dimension(ne)                        :: timedev

          timedev =
     $         1.0d0/dx*(flux_x(i,j,:) - flux_x(i+1,j,:)) +
     $         
     $         compute_timedev_with_openbc(
     $            t,x_map(i),y_map(j),
     $            nodes,i,j,
     $            p_model,
     $            y_direction,
     $            gradient_y, dy,
     $            incoming_y) +
     $            
     $         add_body_forces(
     $            p_model,
     $            t, x_map(i), y_map(j), nodes(i,j,:))

        end function compute_timedev_y_edge_local


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for the
        !> NW or SW corner using open boundary conditions
        !
        !> @date
        !> 21_10_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> governing equations of the physical model
        !        
        !>@param nodes
        !> array of data for the grid points
        !
        !>@param dx
        !> space step along the x-axis
        !
        !>@param dy
        !> space step along the y-axis
        !
        !>@param i
        !> index identifying the grid point in the x-direction
        !
        !>@param j
        !> index identifying the grid point in the y-direction
        !
        !>@param incoming_y
        !> procedure checking the type of characteristic at the edge
        !> in the y-direction
        !
        !>@param gradient_x
        !> procedure computing the gradient along the x-direction
        !
        !>@param gradient_y
        !> procedure computing the gradient along the y-direction
        !
        !>@param timedev
        !> time derivatives modified
        !-------------------------------------------------------------
        function compute_timedev_corner_local(
     $     p_model,
     $     t, x_map, y_map, nodes,
     $     dx,dy, i,j,
     $     incoming_x, incoming_y,
     $     gradient_x, gradient_y)
     $     result(timedev)

          implicit none

          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(:)         , intent(in)    :: x_map
          real(rkind), dimension(:)         , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)     , intent(in)    :: nodes
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          integer(ikind)                    , intent(in)    :: i
          integer(ikind)                    , intent(in)    :: j
          procedure(incoming_proc)                          :: incoming_x
          procedure(incoming_proc)                          :: incoming_y
          procedure(gradient_proc)                          :: gradient_x
          procedure(gradient_proc)                          :: gradient_y
          real(rkind), dimension(ne)                        :: timedev

          timedev =
     $         compute_timedev_with_openbc(
     $            t,x_map(i),y_map(j),
     $            nodes,i,j,
     $            p_model,
     $            x_direction,
     $            gradient_x, dx,
     $            incoming_x) + 
     $         
     $         compute_timedev_with_openbc(
     $            t,x_map(i),y_map(j),
     $            nodes,i,j,
     $            p_model,
     $            y_direction,
     $            gradient_y,dy,
     $            incoming_y) +
     $         
     $         add_body_forces(
     $            p_model,
     $            t,x_map(i),y_map(j),
     $            nodes(i,j,:))

        end function compute_timedev_corner_local


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
        !>@param t
        !> time
        !
        !>@param x
        !> x-coordinate
        !
        !>@param y
        !> y_coordinate
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
        function compute_timedev_with_openbc(
     $     t,x,y,
     $     nodes,i,j,
     $     p_model,
     $     direction,
     $     gradient, dn,
     $     incoming_wave)
     $     result(timedev)

          implicit none

          real(rkind)                  , intent(in) :: t
          real(rkind)                  , intent(in) :: x
          real(rkind)                  , intent(in) :: y
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          type(pmodel_eq)              , intent(in) :: p_model
          integer                      , intent(in) :: direction
          procedure(gradient_proc)                  :: gradient
          real(rkind)                  , intent(in) :: dn
          procedure(incoming_proc)                  :: incoming_wave
          real(rkind), dimension(ne)                :: timedev


          real(rkind), dimension(ne+1)  :: obc_prim_var
          
          real(rkind), dimension(ne)    :: eigenvalues
          real(rkind), dimension(ne,ne) :: left_eigenmatrix
          real(rkind), dimension(ne,ne) :: right_eigenmatrix

          real(rkind), dimension(ne)    :: gradient_prim

          integer                       :: k
          real(rkind), dimension(ne)    :: incoming_amp
          real(rkind), dimension(ne)    :: left_eigenvector

          real(rkind), dimension(ne)    :: timedev_prim
          real(rkind), dimension(ne,ne) :: jacConsPrim


          !determine the nodes for the computation
          !of the eigenquantities
          obc_prim_var = p_model%get_prim_obc_eigenqties(t,x,y,nodes(i,j,:))


          !determination of the eigenquantities
          select case(direction)
            case(x_direction)
               eigenvalues       = p_model%compute_x_eigenvalues_prim(obc_prim_var)
               left_eigenmatrix  = p_model%compute_x_lefteigenvector_prim(obc_prim_var)
               right_eigenmatrix = p_model%compute_x_righteigenvector_prim(obc_prim_var)
               
            case(y_direction)
               eigenvalues       = p_model%compute_y_eigenvalues_prim(obc_prim_var)
               left_eigenmatrix  = p_model%compute_y_lefteigenvector_prim(obc_prim_var)
               right_eigenmatrix = p_model%compute_y_righteigenvector_prim(obc_prim_var)

            case default
               print '(''hedstrom_xy_module'')'
               print '(''compute_timedev_openbc'')'
               print '(''direction not recognized'')'
               stop ''
          end select


          !determination of the gradient 
          gradient_prim = p_model%compute_gradient_prim(nodes,i,j,gradient,dn)


          !construction of the vector of
          !characteristic amplitudes
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

                incoming_amp(k)  = -eigenvalues(k)*DOT_PRODUCT(
     $                                 left_eigenvector,
     $                                 gradient_prim)
             end if

          end do


          !determination of the contribution
          !of the hyperbolic terms to the
          !time derivatives of the primitive
          !variables
          timedev_prim = MATMUL(incoming_amp, right_eigenmatrix)


          !determination of the time derivatives
          !of the conservative variables
          jacConsPrim = p_model%compute_jacobian_cons_to_prim(obc_prim_var)
          timedev     = MATMUL(timedev_prim, jacConsPrim)

        end function compute_timedev_with_openbc

      end module hedstrom_xy_module

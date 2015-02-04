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
     $       obc_outgoing_cons,
     $       obc_outgoing_prim

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       obc_outgoing_strategy

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
     $       compute_timedev_xlayer,
     $       compute_timedev_ylayer,
     $       compute_timedev_xlayer_local,
     $       compute_timedev_ylayer_local,
     $       compute_timedev_xlayer_local_hedstrom,
     $       compute_timedev_ylayer_local_hedstrom,
     $       compute_timedev_corner_W,
     $       compute_timedev_corner_E,
     $       compute_timedev_corner_local,
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
     $     t, x_map, y_map, nodes, i,j, dx,dy, p_model, flux_y,
     $     gradient_x, incoming_x,
     $     timedev)

          implicit none

          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          integer(ikind)                    , intent(in)    :: i
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx,ny+1,ne), intent(in)    :: flux_y
          procedure(gradient_proc)                          :: gradient_x
          procedure(incoming_proc)                          :: incoming_x
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev


          timedev(i,j,:) = compute_timedev_xlayer_local(
     $         p_model,
     $         t, x_map, y_map, nodes, dx,dy, i,j,
     $         flux_y,
     $         incoming_x,
     $         gradient_x)

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
     $     t, x_map, y_map, nodes,
     $     j, dx, dy, p_model, flux_x,
     $     gradient_y, incoming_y,
     $     timedev)

          implicit none

          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx+1,ny,ne), intent(in)    :: flux_x
          procedure(gradient_proc)                          :: gradient_y
          procedure(incoming_proc)                          :: incoming_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          
          integer(ikind) :: i


          do i=3, nx-bc_size

             timedev(i,j,:) = compute_timedev_ylayer_local(
     $            p_model,
     $            t, x_map, y_map, nodes, dx,dy, i,j,
     $            flux_x,
     $            incoming_y,
     $            gradient_y)

          end do

        end subroutine compute_timedev_ylayer


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
        function compute_timedev_xlayer_local(
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
     $         compute_x_timedev_with_openbc(
     $            t, x_map(i), y_map(j),
     $            nodes, i, j, p_model, dx,
     $            gradient_x, incoming_x) +
     $         
     $         1.0d0/dy*(flux_y(i,j,:) - flux_y(i,j+1,:)) +
     $         
     $         add_body_forces(
     $            p_model,
     $            t, x_map(i), y_map(j), nodes(i,j,:))

        end function compute_timedev_xlayer_local


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
        function compute_timedev_ylayer_local(
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
     $         compute_y_timedev_with_openbc(
     $            t, x_map(i), y_map(j),
     $            nodes, i, j, p_model, dy,
     $            gradient_y, incoming_y) +
     $            
     $         add_body_forces(
     $            p_model,
     $            t, x_map(i), y_map(j), nodes(i,j,:))

        end function compute_timedev_ylayer_local


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
        function compute_timedev_xlayer_local_hedstrom(
     $     p_model,
     $     t, x_map, y_map, nodes, dx, dy, i,j,
     $     flux_y,
     $     gradient_x,
     $     side_x)
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
          real(rkind), dimension(:,:,:)     , intent(in)    :: flux_y
          procedure(gradient_proc)                          :: gradient_x
          logical                           , intent(in)    :: side_x
          real(rkind), dimension(ne)                        :: timedev

          if(side_x.eqv.left) then

             timedev = compute_timedev_xlayer_local(
     $            p_model,
     $            t,x_map,y_map, nodes, dx, dy, i,j,
     $            flux_y,
     $            incoming_left,
     $            gradient_x)

          else

             timedev = compute_timedev_xlayer_local(
     $            p_model,
     $            t,x_map,y_map, nodes, dx, dy, i,j,
     $            flux_y,
     $            incoming_right,
     $            gradient_x)

          end if

        end function compute_timedev_xlayer_local_hedstrom


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
        function compute_timedev_ylayer_local_hedstrom(
     $     p_model,
     $     t, x_map, y_map, nodes, dx,dy, i,j,
     $     flux_x,
     $     gradient_y,
     $     side_y)
     $     result(timedev)

          implicit none

          type(pmodel_eq)                   , intent(in) :: p_model
          real(rkind)                       , intent(in) :: t
          real(rkind), dimension(:)         , intent(in) :: x_map
          real(rkind), dimension(:)         , intent(in) :: y_map
          real(rkind), dimension(:,:,:)     , intent(in) :: nodes
          real(rkind)                       , intent(in) :: dx
          real(rkind)                       , intent(in) :: dy
          integer(ikind)                    , intent(in) :: i
          integer(ikind)                    , intent(in) :: j
          real(rkind), dimension(:,:,:)     , intent(in) :: flux_x
          procedure(gradient_proc)                       :: gradient_y
          logical                           , intent(in) :: side_y
          real(rkind), dimension(ne)                     :: timedev


          if(side_y.eqv.left) then

             timedev = compute_timedev_ylayer_local(
     $            p_model,
     $            t,x_map,y_map, nodes, dx, dy, i,j,
     $            flux_x,
     $            incoming_left,
     $            gradient_y)

          else

             timedev = compute_timedev_ylayer_local(
     $            p_model,
     $            t,x_map,y_map, nodes, dx, dy, i,j,
     $            flux_x,
     $            incoming_right,
     $            gradient_y)

          end if

        end function compute_timedev_ylayer_local_hedstrom


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
     $     t,x_map,y_map,nodes, j, dx, dy, p_model,
     $     gradient_y, incoming_y,
     $     timedev)

          implicit none

          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          procedure(gradient_proc)                          :: gradient_y
          procedure(incoming_proc)                          :: incoming_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          
          integer(ikind) :: i          

          i=1
          timedev(i,j,:) = 
     $         compute_timedev_corner_local(
     $         p_model,
     $         t,x_map,y_map,
     $         nodes,
     $         dx,dy,
     $         i,j,
     $         incoming_left,
     $         incoming_y,
     $         gradient_x_x_oneside_L0,
     $         gradient_y)

          i=2
          timedev(i,j,:) = 
     $         compute_timedev_corner_local(
     $         p_model,
     $         t,x_map,y_map,
     $         nodes,
     $         dx,dy,
     $         i,j,
     $         incoming_left,
     $         incoming_y,
     $         gradient_x_x_oneside_L0,
     $         gradient_y)

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
     $     t,x_map,y_map,nodes, j, dx, dy, p_model,
     $     gradient_y, incoming_y,
     $     timedev)

          implicit none

          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          integer(ikind)                    , intent(in)    :: j
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          procedure(gradient_proc)                          :: gradient_y
          procedure(incoming_proc)                          :: incoming_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          integer(ikind) :: i

          i=nx-1
          timedev(i,j,:) = 
     $         compute_timedev_corner_local(
     $         p_model,
     $         t,x_map,y_map,
     $         nodes,
     $         dx,dy,
     $         i,j,
     $         incoming_right,
     $         incoming_y,
     $         gradient_x_x_oneside_R0,
     $         gradient_y)

          i=nx
          timedev(i,j,:) = 
     $         compute_timedev_corner_local(
     $         p_model,
     $         t,x_map,y_map,
     $         nodes,
     $         dx,dy,
     $         i,j,
     $         incoming_right,
     $         incoming_y,
     $         gradient_x_x_oneside_R0,
     $         gradient_y)

        end subroutine compute_timedev_corner_E


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

          timedev = compute_x_timedev_with_openbc(
     $         t,x_map(i),y_map(j),
     $         nodes, i, j, p_model, dx,
     $         gradient_x, incoming_x) + 
     $         
     $         compute_y_timedev_with_openbc(
     $         t,x_map(i),y_map(j),
     $         nodes, i, j, p_model, dy,
     $         gradient_y, incoming_y) +
     $         
     $         add_body_forces(
     $         p_model,
     $         t, x_map(i), y_map(j), nodes(i,j,:))

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
        function compute_x_timedev_with_openbc_prim(
     $     t,x,y,
     $     nodes, i,j,
     $     p_model, dx,
     $     gradient, incoming_wave)
     $     result(timedev)

          implicit none

          real(rkind)                  , intent(in) :: t
          real(rkind)                  , intent(in) :: x
          real(rkind)                  , intent(in) :: y
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: dx
          procedure(gradient_proc)                  :: gradient
          procedure(incoming_proc)                  :: incoming_wave
          real(rkind), dimension(ne)                :: timedev


          real(rkind)                   :: md_lin
          real(rkind)                   :: ux_lin
          real(rkind)                   :: uy_lin
          real(rkind)                   :: P_lin
          real(rkind)                   :: c_lin
          
          real(rkind), dimension(ne)    :: eigenvalues
          integer                       :: k
          real(rkind), dimension(ne)    :: incoming_amp
          real(rkind), dimension(ne,ne) :: left_eigenmatrix_prim
          real(rkind), dimension(ne,ne) :: right_eigenmatrix_prim
          real(rkind), dimension(ne)    :: left_eigenvector_prim
          real(rkind), dimension(ne)    :: gradient_prim
          real(rkind), dimension(ne)    :: timedev_prim
          real(rkind), dimension(ne,ne) :: jacConsPrim


          !determine the nodes for the computation of the
          !eigenquantities
          call p_model%get_prim_obc_eigenqties(
     $         t,x,y,nodes(i,j,:),
     $         md_lin,
     $         ux_lin,
     $         uy_lin,
     $         P_lin,
     $         c_lin)


          !determination of the speed of the amplitude waves
          !along the x-direction
          eigenvalues           = p_model%compute_x_eigenvalues_prim(
     $                                ux_lin,
     $                                c_lin)

          left_eigenmatrix_prim = p_model%compute_x_lefteigenvector_prim(
     $                                md_lin,
     $                                c_lin)


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

                left_eigenvector_prim = left_eigenmatrix_prim(:,k)

                gradient_prim         = p_model%compute_x_gradient_prim(
     $                                      nodes,i,j,
     $                                      gradient,dx)

                incoming_amp(k)       = -eigenvalues(k)*DOT_PRODUCT(
     $                                      left_eigenvector_prim,
     $                                      gradient_prim)

             end if

          end do


          !determination of the contribution of the hyperbolic terms
          !to the time derivatives of the primitive variables
          right_eigenmatrix_prim = p_model%compute_x_righteigenvector_prim(
     $                                 md_lin,
     $                                 c_lin)

          timedev_prim = MATMUL(incoming_amp, right_eigenmatrix_prim)


          ! determination of the time  derivatives of the
          ! conservative variables
          jacConsPrim = p_model%compute_jacobian_cons_to_prim(
     $                      md_lin,
     $                      ux_lin,
     $                      uy_lin,
     $                      P_lin)

          timedev     = MATMUL(timedev_prim, jacConsPrim)

        end function compute_x_timedev_with_openbc_prim


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
        function compute_y_timedev_with_openbc_prim(
     $     t,x,y,
     $     nodes, i,j,
     $     p_model, dy,
     $     gradient, incoming_wave)
     $     result(timedev)

          implicit none

          real(rkind)                  , intent(in) :: t
          real(rkind)                  , intent(in) :: x
          real(rkind)                  , intent(in) :: y
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: dy
          procedure(gradient_proc)                  :: gradient
          procedure(incoming_proc)                  :: incoming_wave
          real(rkind), dimension(ne)                :: timedev


          real(rkind)                   :: md_lin
          real(rkind)                   :: ux_lin
          real(rkind)                   :: uy_lin
          real(rkind)                   :: P_lin
          real(rkind)                   :: c_lin
          
          real(rkind), dimension(ne)    :: eigenvalues
          integer                       :: k
          real(rkind), dimension(ne)    :: incoming_amp
          real(rkind), dimension(ne,ne) :: left_eigenmatrix_prim
          real(rkind), dimension(ne,ne) :: right_eigenmatrix_prim
          real(rkind), dimension(ne)    :: left_eigenvector_prim
          real(rkind), dimension(ne)    :: gradient_prim
          real(rkind), dimension(ne)    :: timedev_prim
          real(rkind), dimension(ne,ne) :: jacConsPrim


          !determine the nodes for the computation of the
          !eigenquantities
          call p_model%get_prim_obc_eigenqties(
     $         t,x,y,nodes(i,j,:),
     $         md_lin,
     $         ux_lin,
     $         uy_lin,
     $         P_lin,
     $         c_lin)


          !determination of the speed of the amplitude waves
          !along the y-direction
          eigenvalues           = p_model%compute_y_eigenvalues_prim(
     $                                uy_lin,
     $                                c_lin)

          left_eigenmatrix_prim = p_model%compute_y_lefteigenvector_prim(
     $                                md_lin,
     $                                c_lin)


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

                left_eigenvector_prim = left_eigenmatrix_prim(:,k)

                gradient_prim         = p_model%compute_gradient_prim(
     $                                      nodes,i,j,
     $                                      gradient,dy)

                incoming_amp(k)       = -eigenvalues(k)*DOT_PRODUCT(
     $                                      left_eigenvector_prim,
     $                                      gradient_prim)

             end if

          end do


          !determination of the contribution of the hyperbolic terms
          !to the time derivatives of the primitive variables
          right_eigenmatrix_prim = p_model%compute_y_righteigenvector_prim(
     $                                 md_lin,
     $                                 c_lin)

          timedev_prim = MATMUL(incoming_amp, right_eigenmatrix_prim)


          ! determination of the time  derivatives of the
          ! conservative variables
          jacConsPrim = p_model%compute_jacobian_cons_to_prim(
     $                      md_lin,
     $                      ux_lin,
     $                      uy_lin,
     $                      P_lin)

          timedev     = MATMUL(timedev_prim, jacConsPrim)

        end function compute_y_timedev_with_openbc_prim

      end module hedstrom_xy_module

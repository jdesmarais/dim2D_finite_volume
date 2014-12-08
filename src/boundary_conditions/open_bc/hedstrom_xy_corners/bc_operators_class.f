      !> @file
      !> class encapsulating subroutines to apply open boundary
      !> conditions at the edge of the computational domain using
      !> hedstrom boundary conditions with special treatment for the
      !> corners
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to apply open boundary
      !> conditions at the edge of the computational domain using
      !> hedstrom boundary conditions with special treatment for the
      !> corners
      !
      !> @date
      !> 06_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use bc_operators_openbc_normal_class, only :
     $     bc_operators_openbc_normal

        use hedstrom_xy_module, only :
     $       compute_timedev_xlayer,
     $       compute_timedev_ylayer,
     $       compute_timedev_xlayer_local_hedstrom,
     $       compute_timedev_ylayer_local_hedstrom

        use hedstrom_xy_corners_module, only :
     $       compute_n_timedev_with_openbc

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right

        use pmodel_eq_class, only :
     $       pmodel_eq

        use parameters_constant, only :
     $       bc_timedev_choice,
     $       left, right,
     $       N,S,E,W,
     $       x_direction,
     $       y_direction,
     $       n1_direction,
     $       n2_direction

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       rkind,ikind

        use sd_operators_fd_module, only :
     $       gradient_x_interior,
     $       gradient_y_interior,
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_L1,
     $       gradient_x_x_oneside_R1,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_L1,
     $       gradient_y_y_oneside_R1,
     $       gradient_y_y_oneside_R0

        use sd_operators_x_oneside_L0_class, only :
     $       sd_operators_x_oneside_L0

        use sd_operators_x_oneside_L1_class, only :
     $       sd_operators_x_oneside_L1

        use sd_operators_x_oneside_R1_class, only :
     $       sd_operators_x_oneside_R1

        use sd_operators_x_oneside_R0_class, only :
     $       sd_operators_x_oneside_R0

        use sd_operators_y_oneside_L0_class, only :
     $       sd_operators_y_oneside_L0

        use sd_operators_y_oneside_L1_class, only :
     $       sd_operators_y_oneside_L1

        use sd_operators_y_oneside_R1_class, only :
     $       sd_operators_y_oneside_R1

        use sd_operators_y_oneside_R0_class, only :
     $       sd_operators_y_oneside_R0

        
        implicit none


        private
        public :: bc_operators



        !> @class bc_operators
        !> class encapsulating subroutines to apply
        !> open boundary conditions in the x and
        !> y directions at the edge of the computational
        !> domain
        !
        !> @param ini
        !> initialize the bc_type attribute of the
        !> boundary conditions
        !
        !> @param apply_bc_on_timedev
        !> apply the open boundary conditions for the time derivatives
        !---------------------------------------------------------------
        type, extends(bc_operators_openbc_normal) :: bc_operators

          contains

          procedure, pass :: ini

          !procedure used w/o field extension
          procedure, pass :: apply_bc_on_timedev => apply_bc_on_timedev_2ndorder

          !procedures used w/ field extension
          procedure, pass :: apply_bc_on_timedev_x_edge
          procedure, pass :: apply_bc_on_timedev_y_edge
          procedure, pass :: apply_bc_on_timedev_xy_corner

        end type bc_operators        
      

        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the main attributes
        !> of the boundary conditions
        !
        !> @date
        !> 04_08_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions initialized
        !
        !>@param p_model
        !> physical model to know the type of the main variables
        !--------------------------------------------------------------
        subroutine ini(this,p_model)
        
          implicit none

          class(bc_operators), intent(inout) :: this
          type(pmodel_eq)    , intent(in)    :: p_model
          
          integer :: neq

          neq = p_model%get_eq_nb()

          this%bc_type = [
     $         bc_timedev_choice,
     $         bc_timedev_choice,
     $         bc_timedev_choice,
     $         bc_timedev_choice]

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the time derivatives with space operators
        !> that are 2nd order accurate in the interior
        !> and 1st order accurate at the boundary
        !
        !> @date
        !> 06_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array of grid point data
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
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param timedev
        !> time derivatives
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev_2ndorder(
     $     this,
     $     p_model,
     $     t,nodes,x_map,y_map,
     $     flux_x,flux_y,
     $     timedev)
        
          implicit none
        
          class(bc_operators)               , intent(in)    :: this
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev


          type(sd_operators_x_oneside_L0) :: s_x_L0
          type(sd_operators_x_oneside_L1) :: s_x_L1
          type(sd_operators_x_oneside_R1) :: s_x_R1
          type(sd_operators_x_oneside_R0) :: s_x_R0
          type(sd_operators_y_oneside_L0) :: s_y_L0
          type(sd_operators_y_oneside_L1) :: s_y_L1
          type(sd_operators_y_oneside_R1) :: s_y_R1
          type(sd_operators_y_oneside_R0) :: s_y_R0


          real(rkind)           :: dx,dy
          integer(ikind)        :: i,j
          integer(ikind)        :: i_min,i_max
          integer(ikind)        :: j_min,j_max

          integer, dimension(4) :: bc_s
          real(rkind)           :: t_s

          
          dx = x_map(2)-x_map(1)
          dy = y_map(2)-y_map(1)


          !prevent unsed parameter warnings while being
          !supress by the compiler afterwards
          bc_s = this%bc_type
          t_s  = t


          !compute the fluxes at the edge of the computational
          !domain
          !compute the fluxes at the edge of the
          !computational domain
          !--------------------------------------------
          !S_edge
          i_min = bc_size+1
          i_max = nx-bc_size+1
          j     = 1
          
          call this%compute_fluxes_for_bc_y_edge(
     $         p_model,
     $         nodes,
     $         s_y_L0, s_y_L1,
     $         s_y_R1, s_y_R0,
     $         dx, dy,
     $         i_min, i_max, j,
     $         S,
     $         flux_x)
          
          
          !E+W_edge
          j_min = bc_size+1
          j_max = ny-bc_size+1
          
          call this%compute_fluxes_for_bc_x_edge(
     $         p_model,
     $         nodes,
     $         s_x_L0, s_x_L1,
     $         s_x_R1, s_x_R0,
     $         dx, dy,
     $         j_min, j_max, i,
     $         E+W,
     $         flux_y)
          
          
          !N_edge
          i_min = bc_size+1
          i_max = nx-bc_size+1
          j     = ny-bc_size+1
          
          call this%compute_fluxes_for_bc_y_edge(
     $         p_model,
     $         nodes,
     $         s_y_L0, s_y_L1,
     $         s_y_R1, s_y_R0,
     $         dx, dy,
     $         i_min, i_max, j,
     $         N,
     $         flux_x)


          !apply the boundary conditions on the south
          !layer
          !--------------------------------------------
          j=1

          i=1
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_L0,
     $         incoming_left,
     $         n2_direction)

          i=2
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_L0,
     $         incoming_left,
     $         n2_direction)

          call compute_timedev_ylayer(
     $         t, x_map, y_map, nodes,
     $         j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_L0, incoming_left,
     $         timedev)


          i=nx-1
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_L0,
     $         incoming_right,
     $         n1_direction)

          i=nx
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_L0,
     $         incoming_right,
     $         n1_direction)

          j=2
          i=1
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_interior,
     $         incoming_left,
     $         n2_direction)

          i=2
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_interior,
     $         gradient_y_interior,
     $         incoming_left,
     $         n2_direction)

          call compute_timedev_ylayer(
     $         t, x_map, y_map, nodes,
     $         j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_L1, incoming_left,
     $         timedev)

          i=nx-1
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_interior,
     $         gradient_y_interior,
     $         incoming_right,
     $         n1_direction)

          i=nx
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_interior,
     $         incoming_right,
     $         n1_direction)


          !apply the boundary conditions on the west and east
          !layers
          do j=bc_size+1, ny-bc_size

             i=1
             call compute_timedev_xlayer(
     $            t, x_map, y_map, nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L0, incoming_left,
     $            timedev)

             i=bc_size
             call compute_timedev_xlayer(
     $            t, x_map, y_map, nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L1, incoming_left,
     $            timedev)

             i=nx-1
             call compute_timedev_xlayer(
     $            t, x_map, y_map, nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_R1, incoming_right,
     $            timedev)

             i=nx
             call compute_timedev_xlayer(
     $            t, x_map, y_map, nodes, i,j, dx,dy, p_model,  flux_y,
     $            gradient_x_x_oneside_R0, incoming_right,
     $            timedev)

          end do


          !apply the boundary conditions on the north layer
          j=ny-1
          
          i=1
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_interior,
     $         incoming_left,
     $         n1_direction)


          i=2
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_interior,
     $         gradient_y_interior,
     $         incoming_left,
     $         n1_direction)

          call compute_timedev_ylayer(
     $         t, x_map, y_map, nodes,
     $         j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_R1, incoming_right,
     $         timedev)

          i=nx-1
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_interior,
     $         gradient_y_interior,
     $         incoming_right,
     $         n2_direction)


          i=nx
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_interior,
     $         incoming_right,
     $         n2_direction)


          j=ny

          i=1
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_R0,
     $         incoming_left,
     $         n1_direction)

          i=2
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_R0,
     $         incoming_left,
     $         n1_direction)

          call compute_timedev_ylayer(
     $         t, x_map, y_map, nodes,
     $         j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_R0, incoming_right,
     $         timedev)

          i=nx-1
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_R0,
     $         incoming_right,
     $         n2_direction)


          i=nx
          timedev(i,j,:) = compute_n_timedev_with_openbc(
     $         nodes, i,j,
     $         p_model, dx, dy,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_R0,
     $         incoming_right,
     $         n2_direction)
        
        end subroutine apply_bc_on_timedev_2ndorder


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> and x edge: W_edge or E_edge
        !
        !> @date
        !> 21_10_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param t
        !> simulation time for boundary conditions depending
        !> on time
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param i
        !> grid point index along the x-axis
        !
        !>@param j
        !> grid point index along the y-axis
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param side_x
        !> edge side to determine the boundary normal vector
        !
        !>@param gradient_x
        !> procedure to compute the gradient along the x-direction
        !> at (i,j)
        !
        !>@param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        function apply_bc_on_timedev_x_edge(
     $     this,
     $     p_model, t,
     $     nodes, x_map, y_map, i,j,
     $     flux_y,
     $     side_x,
     $     gradient_x)
     $     result(timedev)

          implicit none

          class(bc_operators)          , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind), dimension(:,:,:), intent(in) :: flux_y
          logical                      , intent(in) :: side_x
          procedure(gradient_x_proc)                :: gradient_x
          real(rkind), dimension(ne)                :: timedev

          integer, dimension(4) :: bc_s
          real(rkind)           :: dx,dy

          bc_s = this%bc_type
          dx   = x_map(2)-x_map(1)
          dy   = y_map(2)-y_map(1)

          timedev = compute_timedev_xlayer_local_hedstrom(
     $         p_model,
     $         t,x_map,y_map, nodes, dx,dy, i,j,
     $         flux_y,
     $         gradient_x,
     $         side_x)

        end function apply_bc_on_timedev_x_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> an y edge: N_edge or S_edge
        !
        !> @date
        !> 21_10_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param t
        !> simulation time for boundary conditions depending
        !> on time
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param i
        !> grid point index along the x-axis
        !
        !>@param j
        !> grid point index along the y-axis
        !
        !>@param flux_x
        !> fluxes along the y-direction
        !
        !>@param side_y
        !> edge side to determine the boundary normal vector
        !
        !>@param gradient_y
        !> procedure to compute the gradient along the y-direction
        !> at (i,j)
        !
        !>@param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        function apply_bc_on_timedev_y_edge(
     $     this, 
     $     p_model, t,
     $     nodes, x_map, y_map, i,j,
     $     flux_x,
     $     side_y,
     $     gradient_y)
     $     result(timedev)

          implicit none

          class(bc_operators)          , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind), dimension(:,:,:), intent(in) :: flux_x
          logical                      , intent(in) :: side_y
          procedure(gradient_y_proc)                :: gradient_y
          real(rkind), dimension(ne)                :: timedev

          integer, dimension(4) :: bc_s
          real(rkind)           :: dx,dy

          bc_s = this%bc_type
          dx   = x_map(2) - x_map(1)
          dy   = y_map(2) - y_map(1)

          timedev = compute_timedev_ylayer_local_hedstrom(
     $         p_model,
     $         t, x_map, y_map, nodes, dx,dy, i,j,
     $         flux_x,
     $         gradient_y,
     $         side_y)

        end function apply_bc_on_timedev_y_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> a corner: SE_corner, SW_corner, NE_corner, NW_corner
        !
        !> @date
        !> 21_10_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param t
        !> simulation time for boundary conditions depending
        !> on time
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param i
        !> grid point index along the x-axis
        !
        !>@param j
        !> grid point index along the y-axis
        !
        !>@param flux_x
        !> fluxes along the y-direction
        !
        !>@param side_y
        !> edge side to determine the boundary normal vector
        !
        !>@param gradient_y
        !> procedure to compute the gradient along the y-direction
        !> at (i,j)
        !
        !>@param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        function apply_bc_on_timedev_xy_corner(
     $     this,
     $     p_model, t,
     $     nodes, x_map, y_map, i,j,
     $     side_x, side_y,
     $     gradient_x, gradient_y)
     $     result(timedev)

          implicit none

          class(bc_operators)          , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: side_x
          logical                      , intent(in) :: side_y
          procedure(gradient_x_proc)                :: gradient_x
          procedure(gradient_y_proc)                :: gradient_y
          real(rkind), dimension(ne)                :: timedev


          integer, dimension(4) :: bc_s
          real(rkind)           :: t_s
          real(rkind)           :: dx,dy

          bc_s = this%bc_type
          t_s  = t
          dx   = x_map(2) - x_map(1)
          dy   = y_map(2) - y_map(1)

          
          if(side_x.eqv.left) then

             if(side_y.eqv.left) then
                
                timedev = compute_n_timedev_with_openbc(
     $               nodes, i,j,
     $               p_model, dx,dy,
     $               gradient_x,
     $               gradient_y,
     $               incoming_left,
     $               n2_direction)

             else

                timedev = compute_n_timedev_with_openbc(
     $               nodes, i,j,
     $               p_model, dx,dy,
     $               gradient_x,
     $               gradient_y,
     $               incoming_left,
     $               n1_direction)

             end if

          else
             
             if(side_y.eqv.left) then
                
                timedev = compute_n_timedev_with_openbc(
     $               nodes, i,j,
     $               p_model, dx,dy,
     $               gradient_x,
     $               gradient_y,
     $               incoming_right,
     $               n1_direction)

             else

                timedev = compute_n_timedev_with_openbc(
     $               nodes, i,j,
     $               p_model, dx,dy,
     $               gradient_x,
     $               gradient_y,
     $               incoming_right,
     $               n2_direction)

             end if

          end if

        end function apply_bc_on_timedev_xy_corner

      end module bc_operators_class

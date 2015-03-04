      !> @file
      !> class encapsulating subroutines to apply open boundary
      !> conditions at the edge of the computational domain
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to apply open boundary
      !> conditions at the edge of the computational domain
      !
      !> @date
      !> 01_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use bc_operators_openbc_normal_class, only :
     $       bc_operators_openbc_normal

        use bf_layer_errors_module, only :
     $       error_bc_section_type

        use hedstrom_xy_module, only :
     $       compute_timedev_y_layer_interior,
     $       compute_timedev_x_edge_local,
     $       compute_timedev_y_edge_local,
     $       compute_timedev_corner_local

        use hedstrom_xy_anti_corner_flux_module, only :
     $       compute_timedev_anti_corner_with_fluxes

        use hedstrom_xy_anti_corner_diag_flux_module, only :
     $       compute_timedev_anti_corner_with_diag_fluxes

        use interface_primary, only :
     $       gradient_proc

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right

        use parameters_bf_layer, only : 
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       NE_corner_type,
     $       NW_corner_type

        use parameters_constant, only :
     $       bc_timedev_choice,
     $       left,right,
     $       N,S,E,W,
     $       obc_edge_xy_corner,
     $       obc_edge_xy_flux,
     $       obc_edge_xy_diag_flux

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       obc_edge_xy_strategy

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
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

          procedure, pass :: compute_timedev_anti_corner

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
        !> 04_08_2014 - initial version - J.L. Desmarais
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


          real(rkind)    :: dx,dy
          integer(ikind) :: i,j
          integer(ikind) :: i_min, i_max
          integer(ikind) :: j_min, j_max


          !1) determine the space step for the
          !   computation of the fluxes: the
          !   grid size is assumed constant
          !----------------------------------------
          dx = x_map(2)-x_map(1)
          dy = y_map(2)-y_map(1)


          !2) compute the fluxes at the edges of
          !   the computational domain
          !----------------------------------------
          ! S_edge
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
     $         [.true.,.true.],
     $         flux_x)
          
          
          ! E+W_edge
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
     $         [.true.,.true.],
     $         flux_y)
          
          
          ! N_edge
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
     $         [.true.,.true.],
     $         flux_x)

          !3) compute the time derivatives at the
          !   edge of the computational domain
          !----------------------------------------
          ! S layer
          do j=1,bc_size
             call compute_timedev_y_layer_interior(
     $            p_model,
     $            t,nodes,x_map,y_map,
     $            flux_x,
     $            timedev,
     $            j,dx,dy,
     $            gradient_y_y_oneside_L0,
     $            incoming_left)
          end do

          ! E+W layers
          do j=bc_size+1,ny-bc_size
             
             ! W layer
             do i=1,bc_size
                timedev(i,j,:) = compute_timedev_x_edge_local(
     $               p_model,
     $               t,x_map,y_map,
     $               nodes,dx,dy,i,j,
     $               flux_y,
     $               incoming_left,
     $               gradient_x_x_oneside_L0)
             end do

             ! E layer
             do i=nx-bc_size+1,nx
                timedev(i,j,:) = compute_timedev_x_edge_local(
     $               p_model,
     $               t,x_map,y_map,
     $               nodes,dx,dy,i,j,
     $               flux_y,
     $               incoming_right,
     $               gradient_x_x_oneside_R0)
             end do

          end do

          ! N layer
          do j=1,bc_size
             call compute_timedev_y_layer_interior(
     $            p_model,
     $            t,nodes,x_map,y_map,
     $            flux_x,
     $            timedev,
     $            j,dx,dy,
     $            gradient_y_y_oneside_R0,
     $            incoming_right)
          end do
        
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
          procedure(gradient_proc)                  :: gradient_x
          real(rkind), dimension(ne)                :: timedev

          integer, dimension(4) :: bc_s
          real(rkind)           :: dx,dy

          bc_s = this%bc_type
          dx   = x_map(2)-x_map(1)
          dy   = y_map(2)-y_map(1)

          ! W edge
          if(side_x.eqv.left) then
               
             timedev = compute_timedev_x_edge_local(
     $            p_model,
     $            t,x_map,y_map,
     $            nodes,dx,dy,i,j,
     $            flux_y,
     $            incoming_left,
     $            gradient_x_x_oneside_L0)

          ! E edge
          else

             timedev = compute_timedev_x_edge_local(
     $            p_model,
     $            t,x_map,y_map,
     $            nodes,dx,dy,i,j,
     $            flux_y,
     $            incoming_right,
     $            gradient_x_x_oneside_R0)
             
          end if
               
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
          procedure(gradient_proc)                  :: gradient_y
          real(rkind), dimension(ne)                :: timedev

          integer, dimension(4) :: bc_s
          real(rkind)           :: dx,dy

          bc_s = this%bc_type
          dx   = x_map(2) - x_map(1)
          dy   = y_map(2) - y_map(1)

          ! S edge
          if(side_y.eqv.left) then
               
             timedev = compute_timedev_y_edge_local(
     $            p_model,
     $            t,x_map,y_map,
     $            nodes,dx,dy,i,j,
     $            flux_x,
     $            incoming_left,
     $            gradient_y_y_oneside_L0)

          ! N edge
          else

             timedev = compute_timedev_y_edge_local(
     $            p_model,
     $            t,x_map,y_map,
     $            nodes,dx,dy,i,j,
     $            flux_x,
     $            incoming_right,
     $            gradient_y_y_oneside_R0)
             
          end if

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
     $     side_x, side_y)
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
          real(rkind), dimension(ne)                :: timedev

          integer, dimension(4) :: bc_s
          real(rkind)           :: dx
          real(rkind)           :: dy          

          bc_s = this%bc_type
          dx   = x_map(2)-x_map(1)
          dy   = y_map(2)-y_map(1)

          if(side_x.eqv.left) then

             ! SW corner
             if(side_y.eqv.left) then

                timedev = compute_timedev_corner_local(
     $               p_model,
     $               t,x_map,y_map, nodes,
     $               dx,dy, i,j,
     $               incoming_left, incoming_left,
     $               gradient_x_x_oneside_L0, gradient_y_y_oneside_L0)

             ! NW corner
             else
                
                timedev = compute_timedev_corner_local(
     $               p_model,
     $               t,x_map,y_map, nodes,
     $               dx,dy, i,j,
     $               incoming_left, incoming_right,
     $               gradient_x_x_oneside_L0, gradient_y_y_oneside_R0)
                
             end if

          else

             ! SE corner
             if(side_y.eqv.left) then

                timedev = compute_timedev_corner_local(
     $               p_model,
     $               t,x_map,y_map, nodes,
     $               dx,dy, i,j,
     $               incoming_right, incoming_left,
     $               gradient_x_x_oneside_R0, gradient_y_y_oneside_L0)

             ! NE corner
             else
                timedev = compute_timedev_corner_local(
     $               p_model,
     $               t,x_map,y_map, nodes,
     $               dx,dy, i,j,
     $               incoming_right, incoming_right,
     $               gradient_x_x_oneside_R0, gradient_y_y_oneside_R0)
             end if

          end if

        end function apply_bc_on_timedev_xy_corner


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> an xy edge: NE_edge, NW_edge, SE_edge, SW_edge
        !
        !> @date
        !> 26_01_2014 - initial version - J.L. Desmarais
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
        !>@param x_map
        !> coordinates along the x-direction
        !
        !>@param y_map
        !> coordinates along the y-direction
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param s_y_L0
        !> space discretization operator with no grid points
        !> on the left side
        !
        !>@param s_y_L1
        !> space discretization operator with one grid point
        !> on the left side
        !
        !>@param s_y_R1
        !> space discretization operator with one grid point
        !> on the right side
        !
        !>@param s_y_R0
        !> space discretization operator with no grid point
        !> on the right side
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param bc_section
        !> type of edge (NE_edge_type, NW_edge_type, SE_edge_type,
        !> SW_edge_type) and localization of the edge
        !
        !>@param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        subroutine compute_timedev_anti_corner(
     $     this,
     $     p_model, t,
     $     interior_nodes,
     $     bf_alignment,
     $     grdpts_id, nodes, x_map, y_map,
     $     flux_x, flux_y,
     $     s_x_L1, s_x_R1,
     $     s_y_L1, s_y_R1,
     $     dx, dy,
     $     bc_section,
     $     timedev)
        
          implicit none
        
          class(bc_operators)                , intent(in)    :: this
          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)                        , intent(in)    :: t
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment 
          integer       , dimension(:,:)     , intent(in)    :: grdpts_id
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: nodes
          real(rkind)   , dimension(:)       , intent(in)    :: x_map
          real(rkind)   , dimension(:)       , intent(in)    :: y_map
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_x
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_y
          type(sd_operators_x_oneside_L1)    , intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1)    , intent(in)    :: s_x_R1
          type(sd_operators_y_oneside_L1)    , intent(in)    :: s_y_L1
          type(sd_operators_y_oneside_R1)    , intent(in)    :: s_y_R1
          real(rkind)                        , intent(in)    :: dx
          real(rkind)                        , intent(in)    :: dy
          integer       , dimension(5)       , intent(in)    :: bc_section
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev
        

          integer, dimension(5) :: bc_section_modified


          select case(obc_edge_xy_strategy)

            !compute the anti-corner as if it was a corner
            case(obc_edge_xy_corner)

               bc_section_modified = bc_section
               
               select case(bc_section(1))

                 case(NW_edge_type)
                    bc_section_modified(1) = NW_corner_type

                 case(NE_edge_type)
                    bc_section_modified(1) = NE_corner_type

                 case(SW_edge_type)
                    bc_section_modified(1) = SW_corner_type

                 case(SE_edge_type)
                    bc_section_modified(1) = SE_corner_type
                    
                 case default
                    call error_bc_section_type(
     $                   'hedstrom_xy/bc_operators_class',
     $                   'compute_timedev_anti_corner',
     $                   bc_section(1))

               end select

               call this%compute_timedev_corner(
     $              p_model,
     $              t,nodes,x_map,y_map,
     $              bc_section_modified,
     $              timedev)

            ! compute the anti-corner using the fluxes
            case(obc_edge_xy_flux)

               call compute_timedev_anti_corner_with_fluxes(
     $              p_model,
     $              t,grdpts_id,nodes,x_map,y_map,
     $              flux_x, flux_y,
     $              timedev,
     $              s_x_L1, s_x_R1,
     $              s_y_L1, s_y_R1,
     $              dx, dy,
     $              bc_section,
     $              interior_nodes,
     $              bf_alignment)

            ! compute the anti-corner using the diagonal fluxes
            case(obc_edge_xy_diag_flux)

               call compute_timedev_anti_corner_with_diag_fluxes(
     $              p_model,
     $              t,grdpts_id,nodes,x_map,y_map,
     $              timedev,
     $              dx,dy,
     $              bc_section,
     $              interior_nodes,
     $              bf_alignment)

            case default
               print '(''bc_operators_class'')'
               print '(''apply_bc_on_timedev_xy_edge'')'
               print '(''obc_edge_xy_strategy not recognized'')'
               print '(''obc_edge_xy_strategy: '',I2)', obc_edge_xy_strategy

               stop 'obc_edge_xy_strategy'

          end select          

        end subroutine compute_timedev_anti_corner

      end module bc_operators_class

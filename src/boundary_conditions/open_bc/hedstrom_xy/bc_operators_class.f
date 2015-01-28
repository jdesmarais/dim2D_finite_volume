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

        use bc_operators_nopt_module, only :
     $       compute_edge_N,
     $       compute_edge_S,
     $       compute_edge_E,
     $       compute_edge_W,
     $       are_grdpts_needed_for_flux_x,
     $       are_grdpts_needed_for_flux_y,
     $       combine_grdpts_to_compute_fluxes

        use bf_layer_bc_procedure_module, only :
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       NE_corner_type,
     $       NW_corner_type

        use bf_layer_bc_sections_class, only :
     $       determine_edge_points_computed

        use bf_layer_errors_module, only :
     $       error_bc_section_type

        use bf_layer_sync_module, only :
     $       get_bf_layer_match_table

        use hedstrom_xy_module, only :
     $       compute_timedev_xlayer,
     $       compute_timedev_ylayer,
     $       compute_timedev_corner_W,
     $       compute_timedev_corner_E,
     $       compute_timedev_xlayer_local_hedstrom,
     $       compute_timedev_ylayer_local_hedstrom,
     $       compute_timedev_corner_local

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right        

        use parameters_constant, only :
     $       bc_timedev_choice,
     $       left,right,
     $       N,S,E,W,
     $       obc_edge_xy_corner,
     $       obc_edge_xy_flux

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       obc_edge_xy_strategy

        use parameters_kind, only :
     $       rkind,ikind

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

          real(rkind) :: t_s

          
          dx = x_map(2)-x_map(1)
          dy = y_map(2)-y_map(1)


          !prevent unsed parameter warnings while being
          !supress by the compiler afterwards
          !--------------------------------------------
          t_s  = t


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
          call compute_timedev_corner_W(
     $         t,x_map,y_map, nodes, j, dx, dy, p_model,
     $         gradient_y_y_oneside_L0, incoming_left,
     $         timedev)

          call compute_timedev_ylayer(
     $         t,x_map,y_map, nodes,
     $         j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_L0, incoming_left,
     $         timedev)

          call compute_timedev_corner_E(
     $         t,x_map,y_map, nodes, j, dx, dy, p_model,
     $         gradient_y_y_oneside_L0, incoming_left,
     $         timedev)

          j=2
          call compute_timedev_corner_W(
     $         t,x_map,y_map, nodes, j, dx, dy, p_model,
     $         gradient_y_y_oneside_L0, incoming_left,
     $         timedev)

          call compute_timedev_ylayer(
     $         t,x_map,y_map, nodes, j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_L0, incoming_left,
     $         timedev)

          call compute_timedev_corner_E(
     $         t,x_map,y_map, nodes, j, dx, dy, p_model,
     $         gradient_y_y_oneside_L0, incoming_left,
     $         timedev)


          !apply the boundary conditions on the west
          !and east layers
          !--------------------------------------------
          do j=bc_size+1, ny-bc_size

             i=1
             call compute_timedev_xlayer(
     $            t,x_map,y_map, nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L0, incoming_left,
     $            timedev)

             i=bc_size
             call compute_timedev_xlayer(
     $            t,x_map,y_map, nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L0, incoming_left,
     $            timedev)

             i=nx-1
             call compute_timedev_xlayer(
     $            t,x_map,y_map, nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_R0, incoming_right,
     $            timedev)

             i=nx
             call compute_timedev_xlayer(
     $            t,x_map,y_map, nodes, i,j, dx,dy, p_model,  flux_y,
     $            gradient_x_x_oneside_R0, incoming_right,
     $            timedev)

          end do


          !apply the boundary conditions on the north
          !layer
          !--------------------------------------------
          j=ny-1
          call compute_timedev_corner_W(
     $         t,x_map,y_map, nodes, j, dx, dy, p_model,
     $         gradient_y_y_oneside_R0, incoming_right,
     $         timedev)

          call compute_timedev_ylayer(
     $         t,x_map,y_map, nodes, j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_R0, incoming_right,
     $         timedev)

          call compute_timedev_corner_E(
     $         t,x_map,y_map, nodes, j, dx, dy, p_model,
     $         gradient_y_y_oneside_R0, incoming_right,
     $         timedev)

          j=ny
          call compute_timedev_corner_W(
     $         t,x_map,y_map, nodes, j, dx, dy, p_model,
     $         gradient_y_y_oneside_R0, incoming_right,
     $         timedev)

          call compute_timedev_ylayer(
     $         t,x_map,y_map, nodes, j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_R0, incoming_right,
     $         timedev)

          call compute_timedev_corner_E(
     $         t,x_map,y_map, nodes, j, dx, dy, p_model,
     $         gradient_y_y_oneside_R0, incoming_right,
     $         timedev)
        
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
     $     nodes, x_map, y_map,
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
          integer       , dimension(4)       , intent(in)    :: bc_section
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev
        

          integer, dimension(4) :: bc_section_modified


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
     $              this,
     $              p_model, t,
     $              interior_nodes,
     $              bf_alignment,
     $              nodes, x_map, y_map, flux_x, flux_y,
     $              s_x_L1, s_x_R1,
     $              s_y_L1, s_y_R1,
     $              dx, dy,
     $              bc_section,
     $              timedev)


            case default
               print '(''bc_operators_class'')'
               print '(''apply_bc_on_timedev_xy_edge'')'
               print '(''obc_edge_xy_strategy not recognized'')'
               print '(''obc_edge_xy_strategy: '',I2)', obc_edge_xy_strategy

               stop 'obc_edge_xy_strategy'

          end select          

        end subroutine compute_timedev_anti_corner
        

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> an xy edge: NE_edge, NW_edge, SE_edge, SW_edge
        !> using the fluxes
        !
        !> @date
        !> 26_01_2014 - initial version - J.L. Desmarais
        !
        !> @param p_model
        !> object encapsulating the physical model
        !
        !> @param t
        !> simulation time for boundary conditions depending
        !> on time
        !
        !> @param interior_nodes
        !> nodes from the interior computational domain
        !
        !> @param bf_alignment
        !> relative position of the buffer layer compared
        !> to the interiro domain
        !
        !> @param nodes
        !> object encapsulating the main variables
        !
        !> @param x_map
        !> coordinates along the x-direction
        !
        !> @param y_map
        !> coordinates along the y-direction
        !
        !> @param flux_x
        !> fluxes along the x-direction
        !
        !> @param flux_y
        !> fluxes along the y-direction
        !
        !> @param s_x_L1
        !> space discretization operator with one grid point
        !> on the left side (x-direction)
        !     
        !> @param s_x_R1
        !> space discretization operator with one grid point
        !> on the right side (x-direction)
        !
        !> @param s_y_L1
        !> space discretization operator with one grid point
        !> on the left side (y-direction)
        !
        !> @param s_y_R1
        !> space discretization operator with one grid point
        !> on the right side (y-direction)
        !
        !> @param dx
        !> space step along the x-direction
        !
        !> @param dy
        !> space step along the y-direction
        !
        !> @param bc_section
        !> type of edge (NE_edge_type, NW_edge_type, SE_edge_type,
        !> SW_edge_type) and localization of the edge
        !
        !> @param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        subroutine compute_timedev_anti_corner_with_fluxes(
     $     this,
     $     p_model, t,
     $     interior_nodes,
     $     bf_alignment,
     $     nodes, x_map, y_map,
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
          integer       , dimension(4)       , intent(in)    :: bc_section
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev

          
          integer(ikind) :: i_min,j_min
          integer(ikind) :: i,j
          logical        :: side_x,side_y
          logical        :: compute_edge

          logical :: compute_point1
          logical :: compute_point2
          logical :: compute_point3
          logical :: compute_point4


          i_min = bc_section(2)
          j_min = bc_section(3)

          
          call determine_edge_points_computed(
     $         bc_section(4),
     $         compute_point1,
     $         compute_point2,
     $         compute_point3,
     $         compute_point4)


          select case(bc_section(1))

            !  ___ ___
            ! |   |CCC|
            ! |___|CCC|
            ! |   |   |  NE_edge
            ! |___|___|
            !------------
            case(NE_edge_type)               
               
               compute_edge =
     $              compute_edge_N(j_min,y_map,bc_timedev_choice).and.
     $              compute_edge_E(i_min,x_map,bc_timedev_choice)
               
               if(compute_edge) then

                  side_x = right
                  side_y = right

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |CCC|   |  NE_edge(1,1): like NE_corner(1,1)
                  ! |CCC|___|
                  !------------
                  if(compute_point1) then

                     i=i_min
                     j=j_min
                     
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_R0)

                  end if


                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |   |CCC|  NE_edge(2,1): like N_edge
                  ! |___|CCC|
                  !------------
                  if(compute_point2) then

                     i=i_min+1
                     j=j_min


                     ! compute fluxes N_edge
                     call compute_flux_x_anti_corner(
     $                    bf_alignment, nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_y_R1,
     $                    p_model,
     $                    i,j,
     $                    flux_x)
                     
                     ! deduce the time derivatives
                     timedev(i,j,:) = 
     $                    this%apply_bc_on_timedev_y_edge(
     $                    p_model,
     $                    t,nodes,
     $                    x_map,y_map,i,j,
     $                    flux_x,
     $                    side_y,
     $                    gradient_y_y_oneside_R0)
                     
                  end if


                  !  ___ ___
                  ! |CCC|   |
                  ! |CCC|___|
                  ! |   |   |  NE_edge(2,1): like E_edge
                  ! |___|___|
                  !------------
                  if(compute_point3) then

                     i=i_min
                     j=j_min+1
                     
                     ! compute fluxes E_edge
                     call compute_flux_y_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_x_R1,
     $                    p_model,
     $                    i,j,
     $                    flux_y)
                     
                     ! deduce the time derivatives
                     timedev(i,j,:) = 
     $                    this%apply_bc_on_timedev_x_edge(
     $                    p_model,
     $                    t,nodes,
     $                    x_map,y_map,i,j,
     $                    flux_y,
     $                    side_x,
     $                    gradient_x_x_oneside_R0)

                  end if


                  !  ___ ___
                  ! |   |CCC|
                  ! |___|CCC|
                  ! |   |   |  NE_edge(2,2): like NE_corner(2,2)
                  ! |___|___|
                  !------------
                  if(compute_point4) then

                     i=i_min+1
                     j=j_min+1

                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_R0)

                  end if
                     
               end if


            !  ___ ___
            ! |CCC|   |
            ! |CCC|___|
            ! |   |   |  NW_edge
            ! |___|___|
            !------------
            case(NW_edge_type)

               compute_edge =
     $              compute_edge_N(j_min,y_map,bc_timedev_choice).and.
     $              compute_edge_W(i_min,x_map,bc_timedev_choice)
               
               if(compute_edge) then

                  side_x = left
                  side_y = right

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |CCC|   |  NW_edge(1,1): like N_edge
                  ! |CCC|___|
                  !------------
                  if(compute_point1) then

                     i=i_min
                     j=j_min

                     ! compute the x-fluxes
                     call compute_flux_x_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_y_R1,
     $                    p_model,
     $                    i,j,
     $                    flux_x)

                     ! deduce the time derivatives
                     timedev(i,j,:) = 
     $                    this%apply_bc_on_timedev_y_edge(
     $                    p_model,
     $                    t,nodes,
     $                    x_map,y_map,i,j,
     $                    flux_x,
     $                    side_y,
     $                    gradient_y_y_oneside_R0)

                  end if


                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |   |CCC|  NW_edge(2,1): like NW_corner(2,1)
                  ! |___|CCC|
                  !------------
                  if(compute_point2) then

                     i=i_min+1
                     j=j_min

                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_R0)

                  end if


                  !  ___ ___
                  ! |CCC|   |
                  ! |CCC|___|
                  ! |   |   |  NW_edge(1,2): like NW_corner(1,2)
                  ! |___|___|
                  !------------
                  if(compute_point3) then

                     i=i_min
                     j=j_min+1

                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_R0)                     

                  end if


                  !  ___ ___
                  ! |   |CCC|
                  ! |___|CCC|
                  ! |   |   |  NW_edge(2,2): like W_edge
                  ! |___|___|
                  !------------
                  if(compute_point4) then

                     i=i_min+1
                     j=j_min+1


                     ! compute the y-fluxes
                     call compute_flux_y_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_x_L1,
     $                    p_model,
     $                    i,j,
     $                    flux_y)

                     ! deduce the time derivatives
                     timedev(i,j,:) = 
     $                    this%apply_bc_on_timedev_x_edge(
     $                    p_model,
     $                    t,nodes,
     $                    x_map,y_map,i,j,
     $                    flux_y,
     $                    side_x,
     $                    gradient_x_x_oneside_L0)

                  end if

               end if


            !  ___ ___
            ! |   |   |
            ! |___|___|
            ! |CCC|   |  SW_edge
            ! |CCC|___|
            !------------
            case(SW_edge_type)

               compute_edge =
     $              compute_edge_S(j_min,y_map,bc_timedev_choice).and.
     $              compute_edge_W(i_min,x_map,bc_timedev_choice)
               
               if(compute_edge) then

                  side_x = left
                  side_y = left

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |CCC|   |  SW_edge(1,1): like SW_corner(1,1)
                  ! |CCC|___|
                  !------------
                  if(compute_point1) then

                     i=i_min
                     j=j_min
                     
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_L0)

                  end if
                     

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |   |CCC|  SW_edge(2,1): like W_edge
                  ! |___|CCC|
                  !------------
                  if(compute_point2) then

                     i=i_min+1
                     j=j_min
                     
                     ! compute the y-fluxes
                     call compute_flux_y_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_x_L1,
     $                    p_model,
     $                    i,j,
     $                    flux_y)

                     ! deduce the time derivatives
                     timedev(i,j,:) = 
     $                    this%apply_bc_on_timedev_x_edge(
     $                    p_model,
     $                    t,nodes,
     $                    x_map,y_map,i,j,
     $                    flux_y,
     $                    side_x,
     $                    gradient_x_x_oneside_L0)

                  end if


                  !  ___ ___
                  ! |CCC|   |
                  ! |CCC|___|
                  ! |   |   |  SW_edge(2,1): like S_edge
                  ! |___|___|
                  !------------
                  if(compute_point3) then

                     i=i_min
                     j=j_min+1

                     ! compute the x-fluxes
                     call compute_flux_x_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_y_L1,
     $                    p_model,
     $                    i,j,
     $                    flux_x)
                     
                     ! deduce the time derivatives
                     timedev(i,j,:) = 
     $                    this%apply_bc_on_timedev_y_edge(
     $                    p_model,
     $                    t,nodes,
     $                    x_map,y_map,i,j,
     $                    flux_x,
     $                    side_y,
     $                    gradient_y_y_oneside_L0)

                  end if


                  !  ___ ___
                  ! |   |CCC|
                  ! |___|CCC|
                  ! |   |   |  SW_edge(2,2): like SW_corner(2,2)
                  ! |___|___|
                  !------------
                  if(compute_point4) then

                     i=i_min+1
                     j=j_min+1
                     
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_L0)

                  end if

               end if


            !  ___ ___
            ! |   |   |
            ! |___|___|
            ! |   |CCC|  SE_edge
            ! |___|CCC|
            !------------
            case(SE_edge_type)

               compute_edge =
     $              compute_edge_S(j_min,y_map,bc_timedev_choice).and.
     $              compute_edge_E(i_min,x_map,bc_timedev_choice)
               
               if(compute_edge) then

                  side_x = right
                  side_y = left

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |CCC|   |  SE_edge(1,1): like E_edge
                  ! |CCC|___|
                  !------------
                  if(compute_point1) then

                     i=i_min
                     j=j_min
                     
                     ! compute the y-fluxes
                     call compute_flux_y_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_x_R1,
     $                    p_model,
     $                    i,j,
     $                    flux_y)

                     ! deduce the time derivatives
                     timedev(i,j,:) = 
     $                    this%apply_bc_on_timedev_x_edge(
     $                    p_model,
     $                    t,nodes,
     $                    x_map,y_map,i,j,
     $                    flux_y,
     $                    side_x,
     $                    gradient_x_x_oneside_R0)

                  end if
                  

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |   |CCC|  SE_edge(2,1): like SE_corner(2,1)
                  ! |___|CCC|
                  !------------
                  if(compute_point2) then

                     i=i_min+1
                     j=j_min
                     
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_L0)

                  end if
                     

                  !  ___ ___
                  ! |CCC|   |
                  ! |CCC|___|
                  ! |   |   |  SE_edge(2,1): like SE_corner(1,2)
                  ! |___|___|
                  !------------
                  if(compute_point3) then
                     
                     i=i_min
                     j=j_min+1
                     
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_L0)

                  end if
                     
                     
                  !  ___ ___
                  ! |   |CCC|
                  ! |___|CCC|
                  ! |   |   |  SE_edge(2,2): like S_edge
                  ! |___|___|
                  !------------
                  if(compute_point4) then

                     i=i_min+1
                     j=j_min+1
                     
                     ! compute the x-fluxes
                     call compute_flux_x_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_y_L1,
     $                    p_model,
     $                    i,j,
     $                    flux_x)
                     
                     ! deduce the time derivatives
                     timedev(i,j,:) = 
     $                    this%apply_bc_on_timedev_y_edge(
     $                    p_model,
     $                    t,nodes,
     $                    x_map,y_map,i,j,
     $                    flux_x,
     $                    side_y,
     $                    gradient_y_y_oneside_L0)

                  end if

               end if

          end select

        end subroutine compute_timedev_anti_corner_with_fluxes


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
          real(rkind)           :: dx
          real(rkind)           :: dy          

          bc_s = this%bc_type
          t_s  = t
          dx   = x_map(2)-x_map(1)
          dy   = y_map(2)-y_map(1)

          if(side_x.eqv.left) then

             if(side_y.eqv.left) then
                timedev = compute_timedev_corner_local(
     $               p_model,
     $               t,x_map,y_map, nodes, dx, dy, i,j,
     $               incoming_left,
     $               incoming_left,
     $               gradient_x,
     $               gradient_y)

             else
                timedev = compute_timedev_corner_local(
     $               p_model,
     $               t,x_map,y_map, nodes, dx, dy, i,j,
     $               incoming_left,
     $               incoming_right,
     $               gradient_x,
     $               gradient_y)
                
             end if

          else

             if(side_y.eqv.left) then
                timedev = compute_timedev_corner_local(
     $               p_model,
     $               t,x_map,y_map, nodes, dx, dy, i,j,
     $               incoming_right,
     $               incoming_left,
     $               gradient_x,
     $               gradient_y)

             else
                timedev = compute_timedev_corner_local(
     $               p_model,
     $               t,x_map,y_map, nodes, dx, dy, i,j,
     $               incoming_right,
     $               incoming_right,
     $               gradient_x,
     $               gradient_y)
             end if

          end if

        end function apply_bc_on_timedev_xy_corner


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the fluxes in the x-direction needed to compute
        !> the time derivative of an anti-corner
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bf_alignment
        !> relative position of the buffer layer compared to
        !> the interior domain
        !
        !> @param bf_nodes
        !> nodes of the buffer layer
        !
        !> @param interior_nodes
        !> nodes of the interior domain
        !
        !> @param dx
        !> space step along the x-direction
        !
        !> @param dy
        !> space step along the y-direction
        !
        !> @param sd_used
        !> space discretization operator
        !
        !>@param p_model
        !> physical model
        !
        !>@param i
        !> x-index of the anti-corner in the buffer layer
        !
        !>@param j
        !> y-index of the anti-corner in the buffer layer
        !
        !>@param flux_x
        !> fluxes along the x-direction for the buffer layer
        !--------------------------------------------------------------
        subroutine compute_flux_x_anti_corner(
     $     bf_alignment,
     $     bf_nodes,
     $     interior_nodes,
     $     dx,dy,
     $     sd_used,
     $     p_model,
     $     i,j,
     $     flux_x)

          implicit none

          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          real(rkind)                        , intent(in)    :: dx
          real(rkind)                        , intent(in)    :: dy
          class(sd_operators)                , intent(in)    :: sd_used
          type(pmodel_eq)                    , intent(in)    :: p_model
          integer(ikind)                     , intent(in)    :: i
          integer(ikind)                     , intent(in)    :: j
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_x


          logical                        :: grdpts_needed
          integer(ikind), dimension(2,2) :: border_coords
          integer(ikind), dimension(2)   :: cpt_coords

          real(rkind), dimension(:,:,:), allocatable :: tmp_nodes

          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2,2) :: gen_coords


          ! determine whether there are enough grid points
          grdpts_needed = are_grdpts_needed_for_flux_x(
     $         p_model,
     $         sd_used%get_operator_type(),
     $         i,j,
     $         size(bf_nodes,1),size(bf_nodes,2),
     $         border_coords,
     $         cpt_coords)


          ! if there are not enough grid points they should be
          ! extracted and the fluxes are computed from these
          ! temporary grid points
          if(grdpts_needed) then
             
             ! allocate space for the temporary gridpoints
             ! extracted
             allocate(tmp_nodes(
     $            border_coords(1,2)-border_coords(1,1)+1,
     $            border_coords(2,2)-border_coords(2,1)+1,
     $            ne))

             ! compute the general coordinates identifying the
             ! the borders of the gridpoints extracted
             match_table = get_bf_layer_match_table(
     $            bf_alignment)
             
             gen_coords(1,1) = border_coords(1,1) + match_table(1)
             gen_coords(1,2) = border_coords(1,2) + match_table(1)
             gen_coords(2,1) = border_coords(2,1) + match_table(2)
             gen_coords(2,2) = border_coords(2,2) + match_table(2)
             

             ! extract the grid points from the current nodes of
             ! the buffer layer and the interior domain
             call combine_grdpts_to_compute_fluxes(
     $            bf_alignment, bf_nodes,
     $            interior_nodes,
     $            gen_coords,
     $            tmp_nodes)


             !compute the x-fluxes
             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            tmp_nodes,dx,dy,
     $            cpt_coords(1),cpt_coords(2),
     $            sd_used)
             
             flux_x(i+1,j,:) = p_model%compute_flux_x_oneside(
     $            tmp_nodes,dx,dy,
     $            cpt_coords(1)+1,cpt_coords(2),
     $            sd_used)

             deallocate(tmp_nodes)


          ! otherwise the fluxes are directly computed from the
          ! existing nodes
          else

             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            bf_nodes,dx,dy,
     $            i,j,
     $            sd_used)
             
             flux_x(i+1,j,:) = p_model%compute_flux_x_oneside(
     $            bf_nodes,dx,dy,
     $            i+1,j,
     $            sd_used)

          end if

        end subroutine compute_flux_x_anti_corner



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the fluxes in the x-direction needed to compute
        !> the time derivative of an anti-corner
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bf_alignment
        !> relative position of the buffer layer compared to
        !> the interior domain
        !
        !> @param bf_nodes
        !> nodes of the buffer layer
        !
        !> @param interior_nodes
        !> nodes of the interior domain
        !
        !> @param dx
        !> space step along the x-direction
        !
        !> @param dy
        !> space step along the y-direction
        !
        !> @param sd_used
        !> space discretization operator
        !
        !>@param p_model
        !> physical model
        !
        !>@param i
        !> x-index of the anti-corner in the buffer layer
        !
        !>@param j
        !> y-index of the anti-corner in the buffer layer
        !
        !>@param flux_y
        !> fluxes along the y-direction for the buffer layer
        !--------------------------------------------------------------
        subroutine compute_flux_y_anti_corner(
     $     bf_alignment,
     $     bf_nodes,
     $     interior_nodes,
     $     dx,dy,
     $     sd_used,
     $     p_model,
     $     i,j,
     $     flux_y)

          implicit none

          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          real(rkind)                        , intent(in)    :: dx
          real(rkind)                        , intent(in)    :: dy
          class(sd_operators)                , intent(in)    :: sd_used
          type(pmodel_eq)                    , intent(in)    :: p_model
          integer(ikind)                     , intent(in)    :: i
          integer(ikind)                     , intent(in)    :: j
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_y


          logical                        :: grdpts_needed
          integer(ikind), dimension(2,2) :: border_coords
          integer(ikind), dimension(2)   :: cpt_coords

          real(rkind), dimension(:,:,:), allocatable :: tmp_nodes

          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2,2) :: gen_coords


          ! determine whether there are enough grid points
          grdpts_needed = are_grdpts_needed_for_flux_y(
     $         p_model,
     $         sd_used%get_operator_type(),
     $         i,j,
     $         size(bf_nodes,1),size(bf_nodes,2),
     $         border_coords,
     $         cpt_coords)


          ! if there are not enough grid points they should be
          ! extracted and the fluxes are computed from these
          ! temporary grid points
          if(grdpts_needed) then
             
             ! allocate space for the temporary gridpoints
             ! extracted
             allocate(tmp_nodes(
     $            border_coords(1,2)-border_coords(1,1)+1,
     $            border_coords(2,2)-border_coords(2,1)+1,
     $            ne))

             ! compute the general coordinates identifying the
             ! the borders of the gridpoints extracted
             match_table = get_bf_layer_match_table(
     $            bf_alignment)
             
             gen_coords(1,1) = border_coords(1,1) + match_table(1)
             gen_coords(1,2) = border_coords(1,2) + match_table(1)
             gen_coords(2,1) = border_coords(2,1) + match_table(2)
             gen_coords(2,2) = border_coords(2,2) + match_table(2)
             

             ! extract the grid points from the current nodes of
             ! the buffer layer and the interior domain
             call combine_grdpts_to_compute_fluxes(
     $            bf_alignment, bf_nodes,
     $            interior_nodes,
     $            gen_coords,
     $            tmp_nodes)


             !compute the y-fluxes
             flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $            tmp_nodes,dx,dy,
     $            cpt_coords(1),cpt_coords(2),
     $            sd_used)
             
             flux_y(i,j+1,:) = p_model%compute_flux_y_oneside(
     $            tmp_nodes,dx,dy,
     $            cpt_coords(1),cpt_coords(2)+1,
     $            sd_used)

             deallocate(tmp_nodes)


          ! otherwise the fluxes are directly computed from the
          ! existing nodes
          else

             flux_y(i,j,:)   = p_model%compute_flux_y_oneside(
     $            bf_nodes,dx,dy,
     $            i,j,
     $            sd_used)
             
             flux_y(i,j+1,:) = p_model%compute_flux_y_oneside(
     $            bf_nodes,dx,dy,
     $            i,j+1,
     $            sd_used)

          end if

        end subroutine compute_flux_y_anti_corner        

      end module bc_operators_class

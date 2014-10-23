      !> @file
      !> class encapsulating subroutines to apply open boundary
      !> conditions at the edge of the computational domain in
      !> the x-direction and reflection boundary conditions in
      !> the y-direction
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to apply open boundary
      !> conditions at the edge of the computational domain in
      !> the x-direction and reflection boundary conditions in
      !> the y-direction
      !
      !> @date
      !> 04_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use bc_operators_openbc_class, only :
     $     bc_operators_openbc

        use hedstrom_xy_module, only :
     $       compute_timedev_xlayer,
     $       compute_timedev_ylayer

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right,
     $       compute_fluxes_at_the_edges_2ndorder

        use pmodel_eq_class, only :
     $       pmodel_eq

        use parameters_constant, only :
     $       bc_nodes_choice,
     $       bc_timedev_choice,
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       rkind,ikind

        use reflection_xy_module, only :
     $       reflection_y_prefactor

        use sd_operators_fd_module, only :
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
        type, extends(bc_operators_openbc) :: bc_operators

          integer, dimension(ne) :: prefactor_y

          contains

          procedure, pass :: ini
          procedure, pass :: apply_bc_on_nodes
          procedure, pass :: apply_bc_on_timedev => apply_bc_on_timedev_2ndorder

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

          this%bc_type(N) = bc_nodes_choice
          this%bc_type(S) = bc_nodes_choice
          this%bc_type(E) = bc_timedev_choice
          this%bc_type(W) = bc_timedev_choice

          this%prefactor_y = reflection_y_prefactor(p_model)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> along the x and y directions at the edge of the
        !> computational domain
        !
        !> @date
        !> 24_09_2013 - initial version - J.L. Desmarais
        !
        !>@param f
        !> object encapsulating the main variables
        !
        !>@param s
        !> space discretization operators
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(this,nodes)

          implicit none

          class(bc_operators)             , intent(in)    :: this
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes


          integer(ikind) :: i,j
          integer        :: k

          !< compute the reflection b.c. in N and S boundary layers
          do k=1, ne
             do j=1, bc_size
                !DEC$ IVDEP
                do i=1, nx
                   
                   nodes(i,j,k) = 
     $                  this%prefactor_y(k)*nodes(i,2*bc_size+1-j,k)
                   nodes(i,ny-bc_size+j,k) = 
     $                  this%prefactor_y(k)*nodes(i,ny-bc_size-j+1,k)
                   
                end do
             end do
          end do

        end subroutine apply_bc_on_nodes


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
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev


          type(sd_operators_x_oneside_L0) :: s_x_L0
          type(sd_operators_x_oneside_L1) :: s_x_L1
          type(sd_operators_x_oneside_R1) :: s_x_R1
          type(sd_operators_x_oneside_R0) :: s_x_R0


          real(rkind)    :: dx,dy,flux
          integer(ikind) :: i,j
          integer(ikind) :: j_min, j_max

          real(rkind)    :: t_s


          dx = x_map(2)-x_map(1)
          dy = y_map(2)-y_map(1)


          !prevent unsed parameter warnings while being
          !supress by the compiler afterwards
          !--------------------------------------------
          t_s  = t
          flux = flux_x(1,1,1)


          !compute the fluxes at the edge of the
          !computational domain
          !--------------------------------------------
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


          !apply the boundary conditions on the west
          !and east layers
          !--------------------------------------------
          do j=bc_size+1, ny-bc_size

             i=1
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L0, incoming_left,
     $            timedev)

             i=bc_size
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L1, incoming_left,
     $            timedev)

             i=nx-1
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_R1, incoming_right,
     $            timedev)

             i=nx
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model,  flux_y,
     $            gradient_x_x_oneside_R0, incoming_right,
     $            timedev)

          end do
        
        end subroutine apply_bc_on_timedev_2ndorder        

      end module bc_operators_class

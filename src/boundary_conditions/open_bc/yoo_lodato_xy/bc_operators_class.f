      module bc_operators_class

        use bc_operators_openbc_class, only :
     $       bc_operators_openbc

        !lodi edge and corner depending on the type
        !of physical model used
        use lodi_corner_inflow_inflow_class, only :
     $       lodi_corner_inflow_inflow

        use lodi_corner_inflow_outflow_class, only :
     $       lodi_corner_inflow_outflow
        
        use lodi_corner_outflow_outflow_class, only :
     $       lodi_corner_outflow_outflow
        
        use lodi_edge_inflow_class, only :
     $       lodi_edge_inflow

        use lodi_edge_outflow_class, only :
     $       lodi_edge_outflow
        
        !parameters
        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind, rkind

        !physical model
        use pmodel_eq_class, only :
     $       pmodel_eq

        !space discretization methods needed for the
        !computation of the fluxes at the edges
        use sd_operators_x_oneside_L0, only :
     $       sd_operators_x_oneside_L0

        use sd_operators_x_oneside_L1, only :
     $       sd_operators_x_oneside_L1

        use sd_operators_x_oneside_R1, only :
     $       sd_operators_x_oneside_R1

        use sd_operators_x_oneside_R0, only :
     $       sd_operators_x_oneside_R0

        use sd_operators_y_oneside_L0, only :
     $       sd_operators_y_oneside_L0

        use sd_operators_y_oneside_L1, only :
     $       sd_operators_y_oneside_L1

        use sd_operators_y_oneside_R1, only :
     $       sd_operators_y_oneside_R1

        use sd_operators_y_oneside_R0, only :
     $       sd_operators_y_oneside_R0


        implicit none


        type, extends(bc_operators_openbc) :: bc_operators

c$$$          !object encapsulating the procedures for the edges
c$$$          type(lodi_edge_inflow)  :: edge_inflow_bc
c$$$          type(lodi_edge_outflow) :: edge_outflow_bc
c$$$
c$$$          !object encapsulating the procedures for the corners
c$$$          type(lodi_corner_inflow_inflow)   :: corner_inflow_inflow_bc
c$$$          type(lodi_corner_inflow_outflow)  :: corner_inflow_outflow_bc
c$$$          type(lodi_corner_outflow_outflow) :: corner_outflow_outflow_bc
          
          !determine the flow procedure chosen by the user
          !let the flow decide whether it is inflow or outflow
          !or let the user impose whether it should be inflow
          !or outflow
          integer :: flow_user_N
          integer :: flow_user_S
          integer :: flow_user_E
          integer :: flow_user_W          

          !for the temporary use of transverse and viscous lodi terms
          !for the application of boundary conditions on bc_sections
          integer(ikind) :: i_edge_offset
          integer(ikind) :: j_edge_offset

          real(rkind), dimension(:,:,:), allocatable :: transverse_lodi
          real(rkind), dimension(:,:,:), allocatable :: viscous_lodi

          contains

          !for the initialization
          procedure, pass :: ini

c$$$          !for the application of the b.c. on bc_sections
c$$$          procedure, pass :: compute_fluxes_for_bc_x_edge
c$$$          procedure, pass :: compute_fluxes_for_bc_y_edge
c$$$
c$$$          procedure, pass :: apply_bc_on_timedev_x_edge
c$$$          procedure, pass :: apply_bc_on_timedev_y_edge
c$$$          procedure, pass :: apply_bc_on_timedev_xy_corner

          !for the application of the b.c. on the interior
          !fixed domain
          procedure, nopass :: compute_edge_fluxes
c$$$          procedure, nopass :: compute_edge_timedev
c$$$          procedure, pass   :: apply_bc_on_timedev

          !for the computation of the intermediate LODI terms
          procedure, nopass :: compute_lodi_terms_x_edge
          procedure, nopass :: compute_lodi_terms_y_edge
c$$$          procedure, nopass :: compute_timedev_x_edge
c$$$          procedure, nopass :: compute_timedev_y_edge
c$$$          procedure, nopass :: compute_timedev_corner

          !for the tests
          procedure, pass :: set_obc_type
        
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
        !> 10_09_2014 - initial version - J.L. Desmarais
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

          this%flow_user_N = obc_type_N
          this%flow_user_S = obc_type_S
          this%flow_user_E = obc_type_E
          this%flow_user_W = obc_type_W

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the main attributes
        !> of the boundary conditions
        !
        !> @date
        !> 10_09_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions initialized
        !
        !>@param obc_type
        !> integer table containing the 
        !--------------------------------------------------------------
        subroutine set_obc_type(this,obc_type)

          implicit none

          class(bc_operators)  , intent(inout) :: this
          integer, dimension(4), intent(in)    :: obc_type

          this%oneside_flow_N = obc_type(N)
          this%oneside_flow_S = obc_type(S)
          this%oneside_flow_E = obc_type(E)
          this%oneside_flow_W = obc_type(W)

        end subroutine set_obc_type


        !compute the fluxes at the edges of the
        !computational domain
        subroutine compute_edge_fluxes(
     $     nodes,dx,dy,
     $     transverse_lodi_N, transverse_lodi_S,
     $     transverse_lodi_E, transverse_lodi_W,
     $     viscous_lodi_N, viscous_lodi_S,
     $     viscous_lodi_E, viscous_lodi_W,
     $     flux_x, flux_y)

          implicit none

          real(rkind), dimension(nx,ny,ne)                 , intent(in)    :: nodes
          real(rkind)                                      , intent(in)    :: dx
          real(rkind)                                      , intent(in)    :: dy
          real(rkind), dimension(nx-2*bc_size+1,bc_size,ne), intent(out)   :: transverse_lodi_N
          real(rkind), dimension(nx-2*bc_size+1,bc_size,ne), intent(out)   :: transverse_lodi_S
          real(rkind), dimension(bc_size,ny-2*bc_size+1,ne), intent(out)   :: transverse_lodi_E
          real(rkind), dimension(bc_size,ny-2*bc_size+1,ne), intent(out)   :: transverse_lodi_W
          real(rkind), dimension(nx-2*bc_size+1,bc_size,ne), intent(out)   :: viscous_lodi_N
          real(rkind), dimension(nx-2*bc_size+1,bc_size,ne), intent(out)   :: viscous_lodi_S
          real(rkind), dimension(bc_size,ny-2*bc_size+1,ne), intent(out)   :: viscous_lodi_E
          real(rkind), dimension(bc_size,ny-2*bc_size+1,ne), intent(out)   :: viscous_lodi_W
          real(rkind), dimension(nx+1,ny,ne)               , intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne)               , intent(inout) :: flux_y


          type(sd_operators_x_oneside_L0) :: s_x_L0
          type(sd_operators_x_oneside_L1) :: s_x_L1
          type(sd_operators_x_oneside_R1) :: s_x_R1
          type(sd_operators_x_oneside_R0) :: s_x_R0
          type(sd_operators_y_oneside_L0) :: s_y_L0
          type(sd_operators_y_oneside_L1) :: s_y_L1
          type(sd_operators_y_oneside_R1) :: s_y_R1
          type(sd_operators_y_oneside_R0) :: s_y_R0


          integer(ikind) :: i,j


          !fluxes computation
          !-------------------------------------------------------------

          !compute the S layers
          j=1
          do i=bc_size+1, nx-bc_size+1
             
             flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $            nodes,i,j,dx,dy,s_y_L0,
     $            transverse_lodi_S(i-bc_size,j,:),
     $            viscous_lodi_S(i-bc_size,j,:))

          end do

          j=2
          do i=bc_size+1, nx-bc_size+1
             
             flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $            nodes,i,j,dx,dy,s_y_L1,
     $            transverse_lodi_S(i-bc_size,j,:),
     $            viscous_lodi_S(i-bc_size,j,:))

          end do


          !compute the E and W layers
          do j=bc_size+1, ny-bc_size+1

             i=1
             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,i,j,dx,dy,s_x_L0,
     $            transverse_lodi_W(i,j-bc_size,:),
     $            viscous_lodi_W(i,j-bc_size,:))

             i=2
             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,i,j,dx,dy,s_x_L1,
     $            transverse_lodi_W(i,j-bc_size,:),
     $            viscous_lodi_W(i,j-bc_size,:))

             i=nx-bc_size+1
             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,i,j,dx,dy,s_x_R1,
     $            transverse_lodi_E(i-(nx-bc_size),j-bc_size,:),
     $            viscous_lodi_E(i-(nx-bc_size),j-bc_size,:))

             i=nx
             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,i,j,dx,dy,s_x_R0,
     $            transverse_lodi_E(i-(nx-bc_size),j-bc_size,:),
     $            viscous_lodi_E(i-(nx-bc_size),j-bc_size,:))

          end do
          

          !compute the N layer
          j=ny-bc_size+1
          do i=bc_size+1, nx-bc_size+1
             
             flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $            nodes,i,j,dx,dy,s_y_R1,
     $            transverse_lodi_N(i-bc_size,j-(ny-bc_size),:),
     $            viscous_lodi_N(i-bc_size,j-(ny-bc_size),:))
             
          end do

          j=ny
          do i=bc_size+1, nx-bc_size+1
             
             flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $            nodes,i,j,dx,dy,s_y_R0,
     $            transverse_lodi_N(i-bc_size,j-(ny-bc_size),:),
     $            viscous_lodi_N(i-bc_size,j-(ny-bc_size),:))
             
          end do


          !LODI computation
          !-------------------------------------------------------------
          !S layer
          call compute_lodi_terms_x_edge(
     $         p_model,
     $         nodes,
     $         dy,
     $         bc_size+1,
     $         1,
     $         transverse_lodi_S,
     $         viscous_lodi_S,
     $         transverse_lodi_S,
     $         viscous_lodi_S)

          !W layer
          call compute_lodi_terms_x_edge(
     $         p_model,
     $         nodes,
     $         dx,
     $         1,
     $         bc_size+1,
     $         transverse_lodi_W,
     $         viscous_lodi_W,
     $         transverse_lodi_W,
     $         viscous_lodi_W)

          !E layer
          call compute_lodi_terms_x_edge(
     $         p_model,
     $         nodes,
     $         dx,
     $         nx-bc_size+1,
     $         bc_size+1,
     $         transverse_lodi_E,
     $         viscous_lodi_E,
     $         transverse_lodi_E,
     $         viscous_lodi_E)

          !N layer
          call compute_lodi_terms_x_edge(
     $         p_model,
     $         nodes,
     $         dy,
     $         bc_size+1,
     $         ny-bc_size+1,
     $         transverse_lodi_N,
     $         viscous_lodi_N,
     $         transverse_lodi_N,
     $         viscous_lodi_N)

        end subroutine compute_edge_fluxes


c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> subroutine computing the fluxes at the edge of the
c$$$        !> computational domain in the y-direction so that
c$$$        !> the time derivatives for an edge in the x-direction
c$$$        !> can be computed
c$$$        !
c$$$        !> @date
c$$$        !> 10_11_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> abstract boundary conditions
c$$$        !
c$$$        !>@param p_model
c$$$        !> object encapsulating the physical model
c$$$        !
c$$$        !>@param nodes
c$$$        !> array containing the grid point data
c$$$        !
c$$$        !>@param dx
c$$$        !> space step along the x-direction
c$$$        !
c$$$        !>@param dy
c$$$        !> space step along the y-direction
c$$$        !
c$$$        !>@param j_min
c$$$        !> index min along the y-direction corresponding
c$$$        !> to the beginning of the edge layer computed
c$$$        !
c$$$        !>@param j_max
c$$$        !> index max along the y-direction corresponding
c$$$        !> to the end of the edge layer computed
c$$$        !
c$$$        !>@param i
c$$$        !> index along the x-direction positioning the
c$$$        !> the edge boundary layer
c$$$        !
c$$$        !>@param edge_card_coord
c$$$        !> cardinal coordinate identifying the type of
c$$$        !> edge boundary layer
c$$$        !
c$$$        !>@param flux_y
c$$$        !> fluxes along the y-direction
c$$$        !-------------------------------------------------------------
c$$$        subroutine compute_fluxes_for_bc_x_edge(
c$$$     $       this,
c$$$     $       p_model,
c$$$     $       nodes,
c$$$     $       s_x_L0, s_x_L1,
c$$$     $       s_x_R1, s_x_R0,
c$$$     $       dx, dy,
c$$$     $       j_min, j_max, i,
c$$$     $       edge_card_coord,
c$$$     $       flux_y)
c$$$        
c$$$          implicit none            
c$$$        
c$$$          class(bc_operators)            , intent(inout) :: this
c$$$          type(pmodel_eq)                , intent(in)    :: p_model
c$$$          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
c$$$          type(sd_operators_x_oneside_L0), intent(in)    :: s_x_L0
c$$$          type(sd_operators_x_oneside_L1), intent(in)    :: s_x_L1
c$$$          type(sd_operators_x_oneside_R1), intent(in)    :: s_x_R1
c$$$          type(sd_operators_x_oneside_R0), intent(in)    :: s_x_R0
c$$$          real(rkind)                    , intent(in)    :: dx
c$$$          real(rkind)                    , intent(in)    :: dy
c$$$          integer(ikind)                 , intent(in)    :: j_min
c$$$          integer(ikind)                 , intent(in)    :: j_max
c$$$          integer(ikind)                 , intent(in)    :: i
c$$$          integer                        , intent(in)    :: edge_card_coord
c$$$          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_y
c$$$
c$$$          integer(ikind)        :: i_f
c$$$          integer(ikind)        :: j
c$$$          integer, dimension(4) :: bc_s
c$$$          
c$$$
c$$$          !allocate space for temporary array: this%edge_inviscid_flux_y
c$$$          !which is needed for the computation of the transverse LODI
c$$$          !terms
c$$$          if(allocated(this%transverse_lodi)) then
c$$$             print '(''bc_operators_yoolodato_class'')'
c$$$             print '(''compute_fluxes_for_bc_x_edge'')'
c$$$             stop 'this%transverse_lodi already allocated'
c$$$          else
c$$$             allocate(this%transverse_lodi(2,j_max-j_min+1,ne))
c$$$          end if
c$$$
c$$$          !allocate space for temporary array: this%edge_viscid_flux_y
c$$$          !which is needed for the computation of the viscous LODI terms
c$$$          if(allocated(this%viscous_lodi)) then
c$$$             print '(''bc_operators_yoolodato_class'')'
c$$$             print '(''compute_fluxes_for_bc_x_edge'')'
c$$$             stop 'this%edge_viscid_flux_y already allocated'
c$$$          else
c$$$             allocate(this%viscous_lodi(2,j_max-j_min+1,ne))
c$$$          end if
c$$$
c$$$          
c$$$          this%i_edge_offset = i
c$$$          this%j_edge_offset = j_min
c$$$
c$$$
c$$$          select case(edge_card_coord)
c$$$            case(W)
c$$$
c$$$               do j=j_min,j_max
c$$$
c$$$                  !compute the inviscid and viscid fluxes
c$$$                  flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
c$$$     $                 nodes,i,j,dx,dy,s_x_L0,
c$$$     $                 this%transverse_lodi(1,j-j_min+1,:),
c$$$     $                 this%viscous_lodi(1,j-j_min+1,:))
c$$$
c$$$
c$$$                  !compute the inviscid and viscid fluxes
c$$$                  flux_y(i+1,j,:) = p_model%compute_flux_y_by_parts(
c$$$     $                 nodes,i+1,j,dx,dy,s_x_L1,
c$$$     $                 this%transverse_lodi(2,j-j_min+1,:),
c$$$     $                 this%viscous_lodi(2,j-j_min+1,:))
c$$$
c$$$               end do
c$$$               
c$$$            case(E)
c$$$
c$$$               do j=j_min,j_max
c$$$
c$$$                  !compute the inviscid and viscid fluxes
c$$$                  flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
c$$$     $                 nodes,i,j,dx,dy,s_x_R1,
c$$$     $                 this%transverse_lodi(1,j-j_min+1,:),
c$$$     $                 this%viscous_lodi(1,j-j_min+1,:))
c$$$
c$$$                  !compute the inviscid and viscid fluxes
c$$$                  flux_y(i+1,j,:) = p_model%compute_flux_y_by_parts(
c$$$     $                 nodes,i+1,j,dx,dy,s_x_R0,
c$$$     $                 this%transverse_lodi(2,j-j_min+1,:),
c$$$     $                 this%viscous_lodi(2,j-j_min+1,:))
c$$$
c$$$               end do
c$$$
c$$$            case default
c$$$               print '(''bc_operators_yoolodato_class'')'
c$$$               print '(''compute_fluxes_for_bc_x_edge'')'
c$$$               stop 'case not recognized'
c$$$          end select
c$$$
c$$$
c$$$          !compute the transverse_lodi and the viscous_lodi
c$$$          !from the inviscid and viscid fluxes previously
c$$$          !saved in 
c$$$          call compute_lodi_terms_x_edge(
c$$$     $         p_model,
c$$$     $         nodes,
c$$$     $         dy,
c$$$     $         this%i_edge_offset,
c$$$     $         this%j_edge_offset,
c$$$     $         this%transverse_lodi,
c$$$     $         this%viscous_lodi,
c$$$     $         this%transverse_lodi,
c$$$     $         this%viscous_lodi)          
c$$$
c$$$        end subroutine compute_fluxes_for_bc_x_edge
c$$$
c$$$
c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> subroutine computing the fluxes at the edge of the
c$$$        !> computational domain in the x-direction so that
c$$$        !> the time derivatives for an edge in the y-direction
c$$$        !> can be computed
c$$$        !
c$$$        !> @date
c$$$        !> 10_11_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> abstract boundary conditions
c$$$        !
c$$$        !>@param p_model
c$$$        !> object encapsulating the physical model
c$$$        !
c$$$        !>@param nodes
c$$$        !> array containing the grid point data
c$$$        !
c$$$        !>@param dx
c$$$        !> space step along the x-direction
c$$$        !
c$$$        !>@param dy
c$$$        !> space step along the y-direction
c$$$        !
c$$$        !>@param i_min
c$$$        !> index min along the x-direction corresponding
c$$$        !> to the beginning of the edge layer computed
c$$$        !
c$$$        !>@param i_max
c$$$        !> index max along the x-direction corresponding
c$$$        !> to the end of the edge layer computed
c$$$        !
c$$$        !>@param j
c$$$        !> index along the y-direction positioning the
c$$$        !> the edge boundary layer
c$$$        !
c$$$        !>@param edge_card_coord
c$$$        !> cardinal coordinate identifying the type of
c$$$        !> edge boundary layer
c$$$        !
c$$$        !>@param flux_x
c$$$        !> fluxes along the x-direction
c$$$        !-------------------------------------------------------------
c$$$        subroutine compute_fluxes_for_bc_y_edge(
c$$$     $     this,
c$$$     $     p_model,
c$$$     $     nodes,
c$$$     $     s_y_L0, s_y_L1,
c$$$     $     s_y_R1, s_y_R0,
c$$$     $     dx, dy,
c$$$     $     i_min, i_max, j,
c$$$     $     edge_card_coord,
c$$$     $     flux_x)
c$$$        
c$$$          implicit none
c$$$        
c$$$          class(bc_operators_yoolodato)  , intent(inout) :: this
c$$$          type(pmodel_eq)                , intent(in)    :: p_model
c$$$          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
c$$$          type(sd_operators_y_oneside_L0), intent(in)    :: s_y_L0
c$$$          type(sd_operators_y_oneside_L1), intent(in)    :: s_y_L1
c$$$          type(sd_operators_y_oneside_R1), intent(in)    :: s_y_R1
c$$$          type(sd_operators_y_oneside_R0), intent(in)    :: s_y_R0
c$$$          real(rkind)                    , intent(in)    :: dx
c$$$          real(rkind)                    , intent(in)    :: dy
c$$$          integer(ikind)                 , intent(in)    :: i_min
c$$$          integer(ikind)                 , intent(in)    :: i_max
c$$$          integer(ikind)                 , intent(in)    :: j
c$$$          integer                        , intent(in)    :: edge_card_coord
c$$$          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_x
c$$$
c$$$          integer(ikind) :: i
c$$$
c$$$          
c$$$          !allocate space for temporary array: this%transverse_lodi
c$$$          !which is needed for the computation of the transverse LODI
c$$$          !terms
c$$$          if(allocated(this%transverse_lodi)) then
c$$$             print '(''bc_operators_yoolodato_class'')'
c$$$             print '(''compute_fluxes_for_bc_y_edge'')'
c$$$             stop 'this%transverse_lodi already allocated'
c$$$          else
c$$$             allocate(this%transverse_lodi(i_max-i_min+1,2,ne))
c$$$          end if
c$$$
c$$$          !allocate space for temporary array: this%viscous_lodi
c$$$          !which is needed for the computation of the viscous LODI terms
c$$$          if(allocated(this%viscous_lodi)) then
c$$$             print '(''bc_operators_yoolodato_class'')'
c$$$             print '(''compute_fluxes_for_bc_y_edge'')'
c$$$             stop 'this%viscous_lodi already allocated'
c$$$          else
c$$$             allocate(this%viscous_lodi(i_max-i_min+1,2,ne))
c$$$          end if
c$$$
c$$$
c$$$          this%i_edge_offset = i_min
c$$$          this%j_edge_offset = j
c$$$
c$$$
c$$$          select case(edge_card_coord)
c$$$            case(S)               
c$$$
c$$$               do i=i_min, i_max
c$$$
c$$$                  flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
c$$$     $                 nodes,i,j,dx,dy,s_y_L0,
c$$$     $                 this%transverse_lodi(i-i_min+1,1,:),
c$$$     $                 this%viscous_lodi(i-i_min+1,1,:))
c$$$
c$$$               end do
c$$$
c$$$               do i=i_min, i_max
c$$$
c$$$                  flux_x(i,j+1,:) = p_model%compute_flux_x_by_parts(
c$$$     $                 nodes,i,j+1,dx,dy,s_y_L1,
c$$$     $                 this%transverse_lodi(i-i_min+1,2,:),
c$$$     $                 this%viscous_lodi(i-i_min+1,2,:))
c$$$
c$$$               end do
c$$$
c$$$            case(N)
c$$$
c$$$               do i=i_min,i_max
c$$$
c$$$                  flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
c$$$     $                 nodes,i,j,dx,dy,s_y_R1,
c$$$     $                 this%transverse_lodi(i-i_min+1,1,:),
c$$$     $                 this%viscous_lodi(i-i_min+1,1,:))
c$$$
c$$$               end do
c$$$          
c$$$               do i=i_min, i_max
c$$$
c$$$                  flux_x(i,j+1,:) = p_model%compute_flux_x_by_parts(
c$$$     $                 nodes,i,j+1,dx,dy,s_y_R0,
c$$$     $                 this%transverse_lodi(i-i_min+1,2,:),
c$$$     $                 this%viscous_lodi(i-i_min+1,2,:))
c$$$
c$$$               end do
c$$$
c$$$            case default
c$$$               print '(''bc_operators_openbc_class'')'
c$$$               print '(''compute_fluxes_for_bc_y_edge'')'
c$$$               stop 'case not recognized'
c$$$          end select
c$$$
c$$$        end subroutine compute_fluxes_for_bc_y_edge


        !compute the transverse and the viscous LODI terms
        !from the inviscid and the viscous fluxes for an x
        !edge
        subroutine compute_lodi_terms_x_edge(
     $     p_model,
     $     nodes,
     $     dy,
     $     i_offset,
     $     j_offset,
     $     inviscid_flux,
     $     viscid_flux,
     $     transverse_lodi,
     $     viscous_lodi)

          implicit none
          
          type(pmodel_eq)              , intent(in)  :: p_model
          real(rkind), dimension(:,:,:), intent(in)  :: nodes
          real(rkind)                  , intent(in)  :: dy
          integer(ikind)               , intent(in)  :: i_offset
          integer(ikind)               , intent(in)  :: j_offset
          real(rkind), dimension(:,:,:), intent(in)  :: inviscid_flux
          real(rkind), dimension(:,:,:), intent(in)  :: viscid_flux
          real(rkind), dimension(:,:,:), intent(out) :: transverse_lodi
          real(rkind), dimension(:,:,:), intent(out) :: viscous_lodi

          integer(ikind)                :: i
          integer(ikind)                :: j
          real(rkind), dimension(ne)    :: dev
          real(rkind), dimension(ne,ne) :: cons_lodi_matrix

          
          do j=1, size(inviscid_flux,2)-1
             do i=1, bc_size
                
                !compute the conservative lodi matrix
                cons_lodi_matrix = p_model%compute_conservative_lodi_matrix(
     $               nodes(i_offset+i-1,j_offset+j-1,:))

                !compute the transverse LODI from the inviscid flux
                do k=1, ne
                   dev(k) = (inviscid_flux(i,j+1,k)-inviscid_flux(i,j,k))/dy
                end do

                transverse_lodi(i,j,:) = -MATMUL(dev,cons_lodi_matrix)


                !compute the viscous LODI from the viscid flux
                do k=1, ne
                   dev(k) = (viscid_flux(i,j+1,k)-viscid_flux(i,j,k))/dy
                end do

                viscous_lodi(i,j,:)   = p_model%get_epsilon()*MATMUL(dev,cons_lodi_matrix)

             end do
          end do

        end subroutine compute_lodi_terms_x_edge


        !compute the transverse and the viscous LODI terms
        !from the inviscid and the viscous fluxes for an y
        !edge
        subroutine compute_lodi_terms_y_edge(
     $     p_model,
     $     nodes,
     $     dx,
     $     i_offset,
     $     j_offset,
     $     inviscid_flux,
     $     viscid_flux,
     $     transverse_lodi,
     $     viscous_lodi)

          implicit none
          
          type(pmodel_eq)              , intent(in)  :: p_model
          real(rkind), dimension(:,:,:), intent(in)  :: nodes
          real(rkind)                  , intent(in)  :: ds
          integer(ikind)               , intent(in)  :: i_offset
          integer(ikind)               , intent(in)  :: j_offset
          real(rkind), dimension(:,:,:), intent(in)  :: inviscid_flux
          real(rkind), dimension(:,:,:), intent(in)  :: viscid_flux
          real(rkind), dimension(:,:,:), intent(out) :: transverse_lodi
          real(rkind), dimension(:,:,:), intent(out) :: viscous_lodi

          integer(ikind)                :: i
          integer(ikind)                :: j
          real(rkind), dimension(ne)    :: dev
          real(rkind), dimension(ne,ne) :: cons_lodi_matrix

          
          do j=1, bc_size
             do i=1, size(inviscid_flux,2)-1
                
                !compute the conservative lodi matrix
                cons_lodi_matrix = p_model%compute_conservative_lodi_matrix(
     $               nodes(i_offset+i-1,j_offset+j-1,:))


                !compute the transverse LODI from the inviscid flux
                do k=1, ne
                   dev(k) = (inviscid_flux(i,j+1,k)-inviscid_flux(i,j,k))/dx
                end do

                transverse_lodi(i,j,:) = -MATMUL(dev,cons_lodi_matrix)


                !compute the viscous LODI from the viscid flux
                do k=1, ne
                   dev(k) = (viscid_flux(i,j+1,k)-viscid_flux(i,j,k))/dy
                end do

                viscous_lodi(i,j,:)   = p_model%get_epsilon()*MATMUL(dev,cons_lodi_matrix)

             end do
          end do

        end subroutine compute_lodi_terms_y_edge


c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> compute the time derivatives at (i,j) resulting
c$$$        !> of the application of the boundary condition on
c$$$        !> and x edge: W_edge or E_edge
c$$$        !
c$$$        !> @date
c$$$        !> 10_11_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param p_model
c$$$        !> object encapsulating the physical model
c$$$        !
c$$$        !>@param t
c$$$        !> simulation time for boundary conditions depending
c$$$        !> on time
c$$$        !
c$$$        !>@param nodes
c$$$        !> object encapsulating the main variables
c$$$        !
c$$$        !>@param dx
c$$$        !> grid size along the x-axis
c$$$        !
c$$$        !>@param dy
c$$$        !> grid size along the y-axis
c$$$        !
c$$$        !>@param i
c$$$        !> grid point index along the x-axis
c$$$        !
c$$$        !>@param j
c$$$        !> grid point index along the y-axis
c$$$        !
c$$$        !>@param flux_y
c$$$        !> fluxes along the y-direction
c$$$        !
c$$$        !>@param side_x
c$$$        !> edge side to determine the boundary normal vector
c$$$        !
c$$$        !>@param gradient_x
c$$$        !> procedure to compute the gradient along the x-direction
c$$$        !> at (i,j)
c$$$        !
c$$$        !>@param timedev
c$$$        !> time derivatives of the grid points
c$$$        !--------------------------------------------------------------
c$$$        function apply_bc_on_timedev_x_edge(
c$$$     $     this,
c$$$     $     p_model, t,
c$$$     $     nodes, x_map, y_map, i,j,
c$$$     $     flux_y,
c$$$     $     side_x,
c$$$     $     gradient_x)
c$$$     $     result(timedev)
c$$$        
c$$$          implicit none
c$$$        
c$$$          class(bc_operators_yoolodato), intent(in) :: this
c$$$          type(pmodel_eq)              , intent(in) :: p_model
c$$$          real(rkind)                  , intent(in) :: t
c$$$          real(rkind), dimension(:,:,:), intent(in) :: nodes
c$$$          real(rkind), dimension(:)    , intent(in) :: x_map
c$$$          real(rkind), dimension(:)    , intent(in) :: y_map
c$$$          integer(ikind)               , intent(in) :: i
c$$$          integer(ikind)               , intent(in) :: j
c$$$          real(rkind), dimension(:,:,:), intent(in) :: flux_y
c$$$          logical                      , intent(in) :: side_x
c$$$          procedure(gradient_x_proc)                :: gradient_x
c$$$          real(rkind), dimension(ne)                :: timedev
c$$$
c$$$
c$$$          !compute the time derivatives from the LODI terms
c$$$          if(side_x.eqv.left) then
c$$$
c$$$             call compute_timedev_x_edge(
c$$$     $            p_model,
c$$$     $            t,nodes,x_map,y_map,i,j,
c$$$     $            flux_y,
c$$$     $            this%transverse_lodi(i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
c$$$     $            this%viscous_lodi(   i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
c$$$     $            side_x,
c$$$     $            gradient_x,
c$$$     $            this%edge_inflow_bc,
c$$$     $            this%edge_outflow_bc,
c$$$     $            this%flow_user_W,
c$$$     $            timedev)
c$$$
c$$$          else
c$$$
c$$$             call compute_timedev_x_edge(
c$$$     $            p_model,
c$$$     $            t,nodes,x_map,y_map,i,j,
c$$$     $            flux_y,
c$$$     $            this%transverse_lodi(i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
c$$$     $            this%viscous_lodi(   i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
c$$$     $            side_x,
c$$$     $            gradient_x,
c$$$     $            this%edge_inflow_bc,
c$$$     $            this%edge_outflow_bc,
c$$$     $            this%flow_user_E,
c$$$     $            timedev)
c$$$
c$$$          end if
c$$$        
c$$$        end function apply_bc_on_timedev_x_edge
c$$$
c$$$
c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> compute the time derivatives at (i,j) resulting
c$$$        !> of the application of the boundary condition on
c$$$        !> an y edge: N_edge or S_edge
c$$$        !
c$$$        !> @date
c$$$        !> 10_11_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param p_model
c$$$        !> object encapsulating the physical model
c$$$        !
c$$$        !>@param t
c$$$        !> simulation time for boundary conditions depending
c$$$        !> on time
c$$$        !
c$$$        !>@param nodes
c$$$        !> object encapsulating the main variables
c$$$        !
c$$$        !>@param dx
c$$$        !> grid size along the x-axis
c$$$        !
c$$$        !>@param dy
c$$$        !> grid size along the y-axis
c$$$        !
c$$$        !>@param i
c$$$        !> grid point index along the x-axis
c$$$        !
c$$$        !>@param j
c$$$        !> grid point index along the y-axis
c$$$        !
c$$$        !>@param flux_x
c$$$        !> fluxes along the y-direction
c$$$        !
c$$$        !>@param side_y
c$$$        !> edge side to determine the boundary normal vector
c$$$        !
c$$$        !>@param gradient_y
c$$$        !> procedure to compute the gradient along the y-direction
c$$$        !> at (i,j)
c$$$        !
c$$$        !>@param timedev
c$$$        !> time derivatives of the grid points
c$$$        !--------------------------------------------------------------
c$$$        function apply_bc_on_timedev_y_edge(
c$$$     $     this,
c$$$     $     p_model, t,
c$$$     $     nodes, x_map, y_map, i,j,
c$$$     $     flux_x, side_y, gradient_y)
c$$$     $     result(timedev)
c$$$        
c$$$          implicit none
c$$$        
c$$$          class(bc_operators_yoolodato), intent(in) :: this
c$$$          type(pmodel_eq)              , intent(in) :: p_model
c$$$          real(rkind)                  , intent(in) :: t
c$$$          real(rkind), dimension(:,:,:), intent(in) :: nodes
c$$$          real(rkind), dimension(:)    , intent(in) :: x_map
c$$$          real(rkind), dimension(:)    , intent(in) :: y_map
c$$$          integer(ikind)               , intent(in) :: i
c$$$          integer(ikind)               , intent(in) :: j
c$$$          real(rkind), dimension(:,:,:), intent(in) :: flux_x
c$$$          logical                      , intent(in) :: side_y
c$$$          procedure(gradient_y_proc)                :: gradient_y
c$$$          real(rkind), dimension(ne)                :: timedev
c$$$
c$$$
c$$$          !compute the time derivatives from the LODI terms
c$$$          if(side_x.eqv.left) then
c$$$
c$$$             call compute_timedev_y_edge(
c$$$     $            p_model,
c$$$     $            t,nodes,x_map,y_map,i,j,
c$$$     $            flux_y,
c$$$     $            this%transverse_lodi(i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
c$$$     $            this%viscous_lodi(   i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
c$$$     $            side_y,
c$$$     $            gradient_y,
c$$$     $            this%edge_inflow_bc,
c$$$     $            this%edge_outflow_bc,
c$$$     $            this%flow_user_S,
c$$$     $            timedev)
c$$$
c$$$          else
c$$$
c$$$             call compute_timedev_y_edge(
c$$$     $            p_model,
c$$$     $            t,nodes,x_map,y_map,i,j,
c$$$     $            flux_y,
c$$$     $            this%transverse_lodi(i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
c$$$     $            this%viscous_lodi(   i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
c$$$     $            side_y,
c$$$     $            gradient_y,
c$$$     $            this%edge_inflow_bc,
c$$$     $            this%edge_outflow_bc,
c$$$     $            this%flow_user_N,
c$$$     $            timedev)
c$$$
c$$$          end if
c$$$        
c$$$        end function apply_bc_on_timedev_y_edge
c$$$        
c$$$        
c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> compute the time derivatives at (i,j) resulting
c$$$        !> of the application of the boundary condition on
c$$$        !> a corner: SE_corner, SW_corner, NE_corner, NW_corner
c$$$        !
c$$$        !> @date
c$$$        !> 10_11_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param p_model
c$$$        !> object encapsulating the physical model
c$$$        !
c$$$        !>@param t
c$$$        !> simulation time for boundary conditions depending
c$$$        !> on time
c$$$        !
c$$$        !>@param nodes
c$$$        !> object encapsulating the main variables
c$$$        !
c$$$        !>@param dx
c$$$        !> grid size along the x-axis
c$$$        !
c$$$        !>@param dy
c$$$        !> grid size along the y-axis
c$$$        !
c$$$        !>@param i
c$$$        !> grid point index along the x-axis
c$$$        !
c$$$        !>@param j
c$$$        !> grid point index along the y-axis
c$$$        !
c$$$        !>@param flux_x
c$$$        !> fluxes along the y-direction
c$$$        !
c$$$        !>@param side_y
c$$$        !> edge side to determine the boundary normal vector
c$$$        !
c$$$        !>@param gradient_y
c$$$        !> procedure to compute the gradient along the y-direction
c$$$        !> at (i,j)
c$$$        !
c$$$        !>@param timedev
c$$$        !> time derivatives of the grid points
c$$$        !--------------------------------------------------------------
c$$$        function apply_bc_on_timedev_xy_corner(
c$$$     $     this,
c$$$     $     p_model, t,
c$$$     $     nodes, x_map, y_map, i,j,
c$$$     $     side_x, side_y,
c$$$     $     gradient_x, gradient_y)
c$$$     $     result(timedev)
c$$$        
c$$$          implicit none
c$$$        
c$$$          class(bc_operators_yoolodato), intent(in) :: this
c$$$          type(pmodel_eq)              , intent(in) :: p_model
c$$$          real(rkind)                  , intent(in) :: t
c$$$          real(rkind), dimension(:,:,:), intent(in) :: nodes
c$$$          real(rkind), dimension(:)    , intent(in) :: x_map
c$$$          real(rkind), dimension(:)    , intent(in) :: y_map
c$$$          integer(ikind)               , intent(in) :: i
c$$$          integer(ikind)               , intent(in) :: j
c$$$          logical                      , intent(in) :: side_x
c$$$          logical                      , intent(in) :: side_y
c$$$          procedure(gradient_x_proc)                :: gradient_x
c$$$          procedure(gradient_y_proc)                :: gradient_y
c$$$          real(rkind), dimension(ne)                :: timedev
c$$$        
c$$$        end function apply_bc_on_timedev_xy_corner

      end module bc_operators_class

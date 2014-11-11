      module bc_operators_class

        use bc_operators_openbc_class, only :
     $       bc_operators_openbc

        !interfaces for the gradient procedures
        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

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

        !computation of the time derivatives using
        !the edge and corner b.c.
        use lodi_timedev_xy_module, only :
     $       compute_timedev_x_edge_local,
     $       compute_timedev_y_edge_local,
     $       compute_timedev_corner_local

        !parameters
        use parameters_constant, only :
     $       N,S,E,W,
     $       bc_timedev_choice,
     $       left, right

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size,
     $       obc_type_N,
     $       obc_type_S,
     $       obc_type_E,
     $       obc_type_W

        use parameters_kind, only :
     $       ikind, rkind

        !physical model
        use pmodel_eq_class, only :
     $       pmodel_eq

        !space discretization methods needed for the
        !computation of the fluxes at the edges
        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_L1,
     $       gradient_x_x_oneside_R1,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_L1,
     $       gradient_y_y_oneside_R1,
     $       gradient_y_y_oneside_R0
        
        !space discretization methods needed for the
        !computation of the fluxes at the edges
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


        type, extends(bc_operators_openbc) :: bc_operators

          !object encapsulating the procedures for the edges
          type(lodi_edge_inflow)  :: edge_inflow_bc
          type(lodi_edge_outflow) :: edge_outflow_bc

          !object encapsulating the procedures for the corners
          type(lodi_corner_inflow_inflow)   :: corner_inflow_inflow_bc
          type(lodi_corner_inflow_outflow)  :: corner_inflow_outflow_bc
          type(lodi_corner_outflow_outflow) :: corner_outflow_outflow_bc
          
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

          !for the application of the b.c. on the interior
          !fixed domain
          procedure, nopass :: compute_edge_fluxes
          procedure, nopass :: compute_edge_timedev
c$$$          procedure, pass   :: apply_bc_on_timedev

c$$$          !for the application of the b.c. on bc_sections
          procedure, pass :: apply_bc_on_timedev_N_edge
          procedure, pass :: apply_bc_on_timedev_S_edge
          procedure, pass :: apply_bc_on_timedev_E_edge
          procedure, pass :: apply_bc_on_timedev_W_edge
          procedure, pass :: apply_bc_on_timedev_xy_corner

c$$$          procedure, pass :: compute_fluxes_for_bc_x_edge
c$$$          procedure, pass :: compute_fluxes_for_bc_y_edge
c$$$
c$$$          procedure, pass :: apply_bc_on_timedev_x_edge
c$$$          procedure, pass :: apply_bc_on_timedev_y_edge
c$$$          procedure, pass :: apply_bc_on_timedev_xy_corner



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

          this%flow_user_N = obc_type(N)
          this%flow_user_S = obc_type(S)
          this%flow_user_E = obc_type(E)
          this%flow_user_W = obc_type(W)

        end subroutine set_obc_type


        !compute the fluxes at the edges of the
        !computational domain
        subroutine compute_edge_fluxes(
     $     p_model,
     $     nodes,dx,dy,
     $     transverse_lodi_N, transverse_lodi_S,
     $     transverse_lodi_E, transverse_lodi_W,
     $     viscous_lodi_N, viscous_lodi_S,
     $     viscous_lodi_E, viscous_lodi_W,
     $     flux_x, flux_y)

          implicit none

          type(pmodel_eq)                                  , intent(in)    :: p_model
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
     $            nodes,dx,dy,i,j,s_y_L0,
     $            transverse_lodi_S(i-bc_size,j,:),
     $            viscous_lodi_S(i-bc_size,j,:))

          end do

          j=2
          do i=bc_size+1, nx-bc_size+1
             
             flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $            nodes,dx,dy,i,j,s_y_L1,
     $            transverse_lodi_S(i-bc_size,j,:),
     $            viscous_lodi_S(i-bc_size,j,:))

          end do


          !compute the E and W layers
          do j=bc_size+1, ny-bc_size+1

             i=1
             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,dx,dy,i,j,s_x_L0,
     $            transverse_lodi_W(i,j-bc_size,:),
     $            viscous_lodi_W(i,j-bc_size,:))

             i=2
             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,dx,dy,i,j,s_x_L1,
     $            transverse_lodi_W(i,j-bc_size,:),
     $            viscous_lodi_W(i,j-bc_size,:))

             i=nx-bc_size+1
             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,dx,dy,i,j,s_x_R1,
     $            transverse_lodi_E(i-(nx-bc_size),j-bc_size,:),
     $            viscous_lodi_E(i-(nx-bc_size),j-bc_size,:))

             i=nx
             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,dx,dy,i,j,s_x_R0,
     $            transverse_lodi_E(i-(nx-bc_size),j-bc_size,:),
     $            viscous_lodi_E(i-(nx-bc_size),j-bc_size,:))

          end do
          

          !compute the N layer
          j=ny-bc_size+1
          do i=bc_size+1, nx-bc_size+1
             
             flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $            nodes,dx,dy,i,j,s_y_R1,
     $            transverse_lodi_N(i-bc_size,j-(ny-bc_size),:),
     $            viscous_lodi_N(i-bc_size,j-(ny-bc_size),:))
             
          end do

          j=ny
          do i=bc_size+1, nx-bc_size+1
             
             flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $            nodes,dx,dy,i,j,s_y_R0,
     $            transverse_lodi_N(i-bc_size,j-(ny-bc_size),:),
     $            viscous_lodi_N(i-bc_size,j-(ny-bc_size),:))
             
          end do


          !LODI computation
          !-------------------------------------------------------------
          !S layer
          call compute_lodi_terms_y_edge(
     $         p_model,
     $         nodes,
     $         dx,
     $         bc_size+1, 1,
     $         transverse_lodi_S,
     $         viscous_lodi_S,
     $         transverse_lodi_S,
     $         viscous_lodi_S)

          !W layer
          call compute_lodi_terms_x_edge(
     $         p_model,
     $         nodes,
     $         dy,
     $         1, bc_size+1,
     $         transverse_lodi_W,
     $         viscous_lodi_W,
     $         transverse_lodi_W,
     $         viscous_lodi_W)

          !E layer
          call compute_lodi_terms_x_edge(
     $         p_model,
     $         nodes,
     $         dy,
     $         nx-bc_size+1, bc_size+1,
     $         transverse_lodi_E,
     $         viscous_lodi_E,
     $         transverse_lodi_E,
     $         viscous_lodi_E)

          !N layer
          call compute_lodi_terms_y_edge(
     $         p_model,
     $         nodes,
     $         dx,
     $         bc_size+1, ny-bc_size+1,
     $         transverse_lodi_N,
     $         viscous_lodi_N,
     $         transverse_lodi_N,
     $         viscous_lodi_N)

        end subroutine compute_edge_fluxes


        !compute the time derivatives at the edges
        



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
          
          type(pmodel_eq)              , intent(in)    :: p_model
          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          real(rkind)                  , intent(in)    :: dy
          integer(ikind)               , intent(in)    :: i_offset
          integer(ikind)               , intent(in)    :: j_offset
          real(rkind), dimension(:,:,:), intent(inout) :: inviscid_flux
          real(rkind), dimension(:,:,:), intent(inout) :: viscid_flux
          real(rkind), dimension(:,:,:), intent(inout) :: transverse_lodi
          real(rkind), dimension(:,:,:), intent(inout) :: viscous_lodi

          integer(ikind)                :: i
          integer(ikind)                :: j
          integer                       :: k
          real(rkind), dimension(ne)    :: flux_diff
          real(rkind), dimension(ne,ne) :: cons_lodi_matrix

          
          do j=1, size(inviscid_flux,2)-1
             do i=1, bc_size
                
                !compute the conservative lodi matrix
                cons_lodi_matrix = p_model%compute_cons_lodi_matrix_y(
     $               nodes(i_offset+i-1,j_offset+j-1,:))

                !compute the transverse LODI from the inviscid flux
                do k=1, ne
                   flux_diff(k) = (
     $                  inviscid_flux(i,j+1,k)-
     $                  inviscid_flux(i,j,k))/dy
                end do

                transverse_lodi(i,j,:) = -MATMUL(flux_diff,cons_lodi_matrix)


                !compute the viscous LODI from the viscid flux
                do k=1, ne
                   flux_diff(k) = (
     $                  viscid_flux(i,j+1,k)-
     $                  viscid_flux(i,j,k))/dy
                end do

                viscous_lodi(i,j,:)   = p_model%get_viscous_coeff()*MATMUL(flux_diff,cons_lodi_matrix)

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
          
          type(pmodel_eq)              , intent(in)    :: p_model
          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          real(rkind)                  , intent(in)    :: dx
          integer(ikind)               , intent(in)    :: i_offset
          integer(ikind)               , intent(in)    :: j_offset
          real(rkind), dimension(:,:,:), intent(inout) :: inviscid_flux
          real(rkind), dimension(:,:,:), intent(inout) :: viscid_flux
          real(rkind), dimension(:,:,:), intent(inout) :: transverse_lodi
          real(rkind), dimension(:,:,:), intent(inout) :: viscous_lodi

          integer(ikind)                :: i
          integer(ikind)                :: j
          integer                       :: k
          real(rkind), dimension(ne)    :: flux_diff
          real(rkind), dimension(ne,ne) :: cons_lodi_matrix

          
          do j=1, bc_size
             do i=1, size(inviscid_flux,2)-1
                
                !compute the conservative lodi matrix
                cons_lodi_matrix = p_model%compute_cons_lodi_matrix_x(
     $               nodes(i_offset+i-1,j_offset+j-1,:))


                !compute the transverse LODI from the inviscid flux
                do k=1, ne
                   flux_diff(k) = (
     $                  inviscid_flux(i+1,j,k)-
     $                  inviscid_flux(i,j,k))/dx
                end do

                transverse_lodi(i,j,:) = -MATMUL(flux_diff,cons_lodi_matrix)


                !compute the viscous LODI from the viscid flux
                do k=1, ne
                   flux_diff(k) = (
     $                  viscid_flux(i+1,j,k)-
     $                  viscid_flux(i,j,k))/dx
                end do

                viscous_lodi(i,j,:)   = 
     $               p_model%get_viscous_coeff()*
     $               MATMUL(flux_diff,cons_lodi_matrix)

             end do
          end do

        end subroutine compute_lodi_terms_y_edge


        subroutine apply_bc_on_timedev_N_edge(
     $       this,
     $       p_model,
     $       t,nodes,
     $       x_map, y_map,
     $       flux_x,
     $       s_y_L0, s_y_L1,
     $       s_y_R1, s_y_R0,
     $       dx,dy,
     $       i_min, i_max, j_min,
     $       timedev)

          implicit none

          class(bc_operators)            , intent(in)    :: this
          type(pmodel_eq)                , intent(in)    :: p_model
          real(rkind)                    , intent(in)    :: t
          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
          real(rkind), dimension(:)      , intent(in)    :: x_map
          real(rkind), dimension(:)      , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_x
          type(sd_operators_y_oneside_L0), intent(in)    :: s_y_L0
          type(sd_operators_y_oneside_L1), intent(in)    :: s_y_L1
          type(sd_operators_y_oneside_R1), intent(in)    :: s_y_R1
          type(sd_operators_y_oneside_R0), intent(in)    :: s_y_R0
          real(rkind)                    , intent(in)    :: dx
          real(rkind)                    , intent(in)    :: dy
          integer(ikind)                 , intent(in)    :: i_min
          integer(ikind)                 , intent(in)    :: i_max
          integer(ikind)                 , intent(in)    :: j_min
          real(rkind), dimension(:,:,:)  , intent(inout) :: timedev

          real(rkind), dimension(:,:,:), allocatable :: transverse_lodi_N
          real(rkind), dimension(:,:,:), allocatable :: viscous_lodi_N
          integer(ikind)                             :: i_edge_offset
          integer(ikind)                             :: j_edge_offset
          integer(ikind)                             :: i,j
          logical                                    :: side_y
          integer(ikind)                             :: bc_s


          bc_s = s_y_L0%get_bc_size() + s_y_L1%get_bc_size()


          !allocate space for the transverse and lodi terms
          allocate(transverse_lodi_N(i_max-i_min+1,2,ne))
          allocate(viscous_lodi_N(i_max-i_min+1,2,ne))


          !identify the offset for matching the LODI terms
          !with the interior nodes
          i_edge_offset = i_min
          j_edge_offset = j_min


          !compute the fluxes
          j=j_min
          !-----------------
          do i=i_min,i_max+1

             flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $            nodes,dx,dy,i,j,s_y_R1,
     $            transverse_lodi_N(i-i_min+1,1,:),
     $            viscous_lodi_N(i-i_min+1,1,:))
             
          end do
          
          j=j_min+1
          !------------------
          do i=i_min, i_max+1
             
             flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $            nodes,dx,dy,i,j,s_y_R0,
     $            transverse_lodi_N(i-i_min+1,2,:),
     $            viscous_lodi_N(i-i_min+1,2,:))
             
          end do

          !compute the extra LODI terms
          call compute_lodi_terms_y_edge(
     $         p_model,
     $         nodes,
     $         dx,
     $         i_edge_offset,
     $         j_edge_offset,
     $         transverse_lodi_N,
     $         viscous_lodi_N,
     $         transverse_lodi_N,
     $         viscous_lodi_N)

          !apply the boundary conditions on the
          !time derivatives
          side_y = right
          j=j_min
          do i=i_min, i_max
             timedev(i,j,:) = compute_timedev_y_edge_local(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_x,
     $            transverse_lodi_N(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            viscous_lodi_N(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            side_y,
     $            gradient_y_y_oneside_R1,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_N)

          end do

          j=j_min+1
          do i=i_min, i_max
             timedev(i,j,:) = compute_timedev_y_edge_local(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_x,
     $            transverse_lodi_N(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            viscous_lodi_N(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            side_y,
     $            gradient_y_y_oneside_R0,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_N)

          end do

          !deallocate spaces allocated fro the temporary tables
          deallocate(transverse_lodi_N)
          deallocate(viscous_lodi_N)
             
        end subroutine apply_bc_on_timedev_N_edge


        subroutine apply_bc_on_timedev_S_edge(
     $       this,
     $       p_model,
     $       t,nodes,
     $       x_map, y_map,
     $       flux_x,
     $       s_y_L0, s_y_L1,
     $       s_y_R1, s_y_R0,
     $       dx,dy,
     $       i_min, i_max, j_min,
     $       timedev)

          implicit none

          class(bc_operators)            , intent(in)    :: this
          type(pmodel_eq)                , intent(in)    :: p_model
          real(rkind)                    , intent(in)    :: t
          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
          real(rkind), dimension(:)      , intent(in)    :: x_map
          real(rkind), dimension(:)      , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_x
          type(sd_operators_y_oneside_L0), intent(in)    :: s_y_L0
          type(sd_operators_y_oneside_L1), intent(in)    :: s_y_L1
          type(sd_operators_y_oneside_R1), intent(in)    :: s_y_R1
          type(sd_operators_y_oneside_R0), intent(in)    :: s_y_R0
          real(rkind)                    , intent(in)    :: dx
          real(rkind)                    , intent(in)    :: dy
          integer(ikind)                 , intent(in)    :: i_min
          integer(ikind)                 , intent(in)    :: i_max
          integer(ikind)                 , intent(in)    :: j_min
          real(rkind), dimension(:,:,:)  , intent(inout) :: timedev


          real(rkind), dimension(:,:,:), allocatable :: transverse_lodi_S
          real(rkind), dimension(:,:,:), allocatable :: viscous_lodi_S
          integer(ikind)                             :: i_edge_offset
          integer(ikind)                             :: j_edge_offset
          integer(ikind)                             :: i,j
          logical                                    :: side_y
          integer(ikind)                             :: bc_s

          bc_s = s_y_R0%get_bc_size() + s_y_R1%get_bc_size()


          !allocate space for the transverse and lodi terms
          allocate(transverse_lodi_S(i_max-i_min+1,2,ne))
          allocate(viscous_lodi_S(i_max-i_min+1,2,ne))


          !identify the offset for matching the LODI terms
          !with the interior nodes
          i_edge_offset = i_min
          j_edge_offset = j_min


          !compute the fluxes
          j=j_min
          !-----------------
          do i=i_min,i_max+1

             flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $            nodes,dx,dy,i,j,s_y_L0,
     $            transverse_lodi_S(i-i_min+1,1,:),
     $            viscous_lodi_S(i-i_min+1,1,:))
             
          end do
          
          j=j_min+1
          !------------------
          do i=i_min, i_max+1
             
             flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $            nodes,dx,dy,i,j,s_y_L1,
     $            transverse_lodi_S(i-i_min+1,2,:),
     $            viscous_lodi_S(i-i_min+1,2,:))
             
          end do

          !compute the extra LODI terms
          call compute_lodi_terms_y_edge(
     $         p_model,
     $         nodes,
     $         dx,
     $         i_edge_offset,
     $         j_edge_offset,
     $         transverse_lodi_S,
     $         viscous_lodi_S,
     $         transverse_lodi_S,
     $         viscous_lodi_S)

          !apply the boundary conditions on the
          !time derivatives
          side_y = left
          j=j_min
          do i=i_min, i_max
             timedev(i,j,:) = compute_timedev_y_edge_local(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_x,
     $            transverse_lodi_S(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            viscous_lodi_S(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            side_y,
     $            gradient_y_y_oneside_L0,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_S)

          end do

          j=j_min+1
          do i=i_min, i_max
             timedev(i,j,:) = compute_timedev_y_edge_local(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_x,
     $            transverse_lodi_S(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            viscous_lodi_S(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            side_y,
     $            gradient_y_y_oneside_L1,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_S)

          end do

          !deallocate spaces allocated fro the temporary tables
          deallocate(transverse_lodi_S)
          deallocate(viscous_lodi_S)

        end subroutine apply_bc_on_timedev_S_edge


        subroutine apply_bc_on_timedev_E_edge(
     $     this,
     $     p_model,
     $     t,nodes,
     $     x_map, y_map,
     $     flux_y,
     $     s_x_L0, s_x_L1,
     $     s_x_R1, s_x_R0,
     $     dx,dy,
     $     j_min, j_max, i_min,
     $     timedev)

          implicit none

          class(bc_operators)            , intent(in)    :: this
          type(pmodel_eq)                , intent(in)    :: p_model
          real(rkind)                    , intent(in)    :: t
          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
          real(rkind), dimension(:)      , intent(in)    :: x_map
          real(rkind), dimension(:)      , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_y
          type(sd_operators_x_oneside_L0), intent(in)    :: s_x_L0
          type(sd_operators_x_oneside_L1), intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1), intent(in)    :: s_x_R1
          type(sd_operators_x_oneside_R0), intent(in)    :: s_x_R0
          real(rkind)                    , intent(in)    :: dx
          real(rkind)                    , intent(in)    :: dy
          integer(ikind)                 , intent(in)    :: j_min
          integer(ikind)                 , intent(in)    :: j_max
          integer(ikind)                 , intent(in)    :: i_min
          real(rkind), dimension(:,:,:)  , intent(inout) :: timedev


          real(rkind), dimension(:,:,:), allocatable :: transverse_lodi_E
          real(rkind), dimension(:,:,:), allocatable :: viscous_lodi_E
          integer(ikind)                             :: i_edge_offset
          integer(ikind)                             :: j_edge_offset
          integer(ikind)                             :: i,j
          logical                                    :: side_x
          integer(ikind)                             :: bc_s

          bc_s = s_x_L0%get_bc_size() + s_x_L1%get_bc_size()


          !allocate space for the transverse and lodi terms
          allocate(transverse_lodi_E(2,j_max-j_min+1,ne))
          allocate(viscous_lodi_E(2,j_max-j_min+1,ne))


          !identify the offset for matching the LODI terms
          !with the interior nodes
          i_edge_offset = i_min
          j_edge_offset = j_min


          !compute the fluxes
          do j=j_min,j_max+1

             i=i_min

             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,dx,dy,i,j,s_x_R1,
     $            transverse_lodi_E(1,j-j_min+1,:),
     $            viscous_lodi_E(1,j-j_min+1,:))


             i=i_min+1

             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,dx,dy,i,j,s_x_R0,
     $            transverse_lodi_E(2,j-j_min+1,:),
     $            viscous_lodi_E(2,j-j_min+1,:))

          end do


          !compute the extra LODI terms
          call compute_lodi_terms_x_edge(
     $         p_model,
     $         nodes,
     $         dy,
     $         i_edge_offset,
     $         j_edge_offset,
     $         transverse_lodi_E,
     $         viscous_lodi_E,
     $         transverse_lodi_E,
     $         viscous_lodi_E)


          !apply the boundary conditions on the
          !time derivatives
          side_x = right
          do j=j_min, j_max

             i=i_min+1
             timedev(i,j,:) = compute_timedev_x_edge_local(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_y,
     $            transverse_lodi_E(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            viscous_lodi_E(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            side_x,
     $            gradient_x_x_oneside_R1,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_E)


             i=i_min+1
             timedev(i,j,:) = compute_timedev_x_edge_local(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_y,
     $            transverse_lodi_E(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            viscous_lodi_E(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            side_x,
     $            gradient_x_x_oneside_R0,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_E)

          end do


          !deallocate spaces allocated fro the temporary tables
          deallocate(transverse_lodi_E)
          deallocate(viscous_lodi_E)


        end subroutine apply_bc_on_timedev_E_edge


        subroutine apply_bc_on_timedev_W_edge(
     $     this,
     $     p_model,
     $     t,nodes,
     $     x_map, y_map,
     $     flux_y,
     $     s_x_L0, s_x_L1,
     $     s_x_R1, s_x_R0,
     $     dx,dy,
     $     j_min, j_max, i_min,
     $     timedev)

          implicit none

          class(bc_operators)            , intent(in)    :: this
          type(pmodel_eq)                , intent(in)    :: p_model
          real(rkind)                    , intent(in)    :: t
          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
          real(rkind), dimension(:)      , intent(in)    :: x_map
          real(rkind), dimension(:)      , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_y
          type(sd_operators_x_oneside_L0), intent(in)    :: s_x_L0
          type(sd_operators_x_oneside_L1), intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1), intent(in)    :: s_x_R1
          type(sd_operators_x_oneside_R0), intent(in)    :: s_x_R0
          real(rkind)                    , intent(in)    :: dx
          real(rkind)                    , intent(in)    :: dy
          integer(ikind)                 , intent(in)    :: j_min
          integer(ikind)                 , intent(in)    :: j_max
          integer(ikind)                 , intent(in)    :: i_min
          real(rkind), dimension(:,:,:)  , intent(inout) :: timedev


          real(rkind), dimension(:,:,:), allocatable :: transverse_lodi_W
          real(rkind), dimension(:,:,:), allocatable :: viscous_lodi_W
          integer(ikind)                             :: i_edge_offset
          integer(ikind)                             :: j_edge_offset
          integer(ikind)                             :: i,j
          logical                                    :: side_x
          integer(ikind)                             :: bc_s

          bc_s = s_x_R0%get_bc_size() + s_x_R1%get_bc_size()


          !allocate space for the transverse and lodi terms
          allocate(transverse_lodi_W(2,j_max-j_min+1,ne))
          allocate(viscous_lodi_W(2,j_max-j_min+1,ne))


          !identify the offset for matching the LODI terms
          !with the interior nodes
          i_edge_offset = i_min
          j_edge_offset = j_min


          !compute the fluxes
          do j=j_min,j_max+1

             i=i_min

             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,dx,dy,i,j,s_x_L0,
     $            transverse_lodi_W(1,j-j_min+1,:),
     $            viscous_lodi_W(1,j-j_min+1,:))


             i=i_min+1

             flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $            nodes,dx,dy,i,j,s_x_L1,
     $            transverse_lodi_W(2,j-j_min+1,:),
     $            viscous_lodi_W(2,j-j_min+1,:))

          end do


          !compute the extra LODI terms
          call compute_lodi_terms_x_edge(
     $         p_model,
     $         nodes,
     $         dy,
     $         i_edge_offset,
     $         j_edge_offset,
     $         transverse_lodi_W,
     $         viscous_lodi_W,
     $         transverse_lodi_W,
     $         viscous_lodi_W)


          !apply the boundary conditions on the
          !time derivatives
          side_x = left
          do j=j_min, j_max

             i=i_min+1            
             timedev(i,j,:) = compute_timedev_x_edge_local(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_y,
     $            transverse_lodi_W(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            viscous_lodi_W(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            side_x,
     $            gradient_x_x_oneside_L0,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_W)


             i=i_min+1
             timedev(i,j,:) = compute_timedev_x_edge_local(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_y,
     $            transverse_lodi_W(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            viscous_lodi_W(i-i_edge_offset+1,j-j_edge_offset+1,:),
     $            side_x,
     $            gradient_x_x_oneside_L1,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_W)

          end do


          !deallocate spaces allocated fro the temporary tables
          deallocate(transverse_lodi_W)
          deallocate(viscous_lodi_W)


        end subroutine apply_bc_on_timedev_W_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> a corner: SE_corner, SW_corner, NE_corner, NW_corner
        !
        !> @date
        !> 10_11_2014 - initial version - J.L. Desmarais
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
           

           if(side_x.eqv.left) then
              
              if(side_y.eqv.left) then
                 timedev = compute_timedev_corner_local(
     $                p_model,
     $                t, nodes, x_map, y_map, i,j,
     $                side_x, side_y,
     $                gradient_x, gradient_y,
     $                this%corner_inflow_inflow_bc,
     $                this%corner_inflow_outflow_bc,
     $                this%corner_outflow_outflow_bc,
     $                this%flow_user_W, this%flow_user_S)
                 
              else
                 timedev = compute_timedev_corner_local(
     $                p_model,
     $                t, nodes, x_map, y_map, i,j,
     $                side_x, side_y,
     $                gradient_x, gradient_y,
     $                this%corner_inflow_inflow_bc,
     $                this%corner_inflow_outflow_bc,
     $                this%corner_outflow_outflow_bc,
     $                this%flow_user_W, this%flow_user_N)
                 
              end if

           else

              if(side_y.eqv.left) then
                 timedev = compute_timedev_corner_local(
     $                p_model,
     $                t, nodes, x_map, y_map, i,j,
     $                side_x, side_y,
     $                gradient_x, gradient_y,
     $                this%corner_inflow_inflow_bc,
     $                this%corner_inflow_outflow_bc,
     $                this%corner_outflow_outflow_bc,
     $                this%flow_user_E, this%flow_user_S)

              else
                 timedev = compute_timedev_corner_local(
     $                p_model,
     $                t, nodes, x_map, y_map, i,j,
     $                side_x, side_y,
     $                gradient_x, gradient_y,
     $                this%corner_inflow_inflow_bc,
     $                this%corner_inflow_outflow_bc,
     $                this%corner_outflow_outflow_bc,
     $                this%flow_user_E, this%flow_user_N)
              end if

           end if

        end function apply_bc_on_timedev_xy_corner

      end module bc_operators_class

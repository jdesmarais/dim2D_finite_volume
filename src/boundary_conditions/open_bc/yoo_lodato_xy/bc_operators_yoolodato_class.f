      module bc_operators_class

        implicit none

        use bc_operators_openbc_class, only : bc_operators_openbc


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
          procedure, nopass :: compute_edge_timedev
          procedure, pass   :: apply_bc_on_timedev          

          !for the computation of the intermediate LODI terms
          procedure, nopass :: compute_lodi_terms_x_edge
          procedure, nopass :: compute_lodi_terms_y_edge
          procedure, nopass :: compute_timedev_x_edge
          procedure, nopass :: compute_timedev_y_edge
          procedure, nopass :: compute_timedev_corner

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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the fluxes at the edge of the
        !> computational domain in the y-direction so that
        !> the time derivatives for an edge in the x-direction
        !> can be computed
        !
        !> @date
        !> 10_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param j_min
        !> index min along the y-direction corresponding
        !> to the beginning of the edge layer computed
        !
        !>@param j_max
        !> index max along the y-direction corresponding
        !> to the end of the edge layer computed
        !
        !>@param i
        !> index along the x-direction positioning the
        !> the edge boundary layer
        !
        !>@param edge_card_coord
        !> cardinal coordinate identifying the type of
        !> edge boundary layer
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !-------------------------------------------------------------
        subroutine compute_fluxes_for_bc_x_edge(
     $       this,
     $       p_model,
     $       nodes,
     $       s_x_L0, s_x_L1,
     $       s_x_R1, s_x_R0,
     $       dx, dy,
     $       j_min, j_max, i,
     $       edge_card_coord,
     $       flux_y)
        
          implicit none            
        
          class(bc_operators_yoolodato)  , intent(inout) :: this
          type(pmodel_eq)                , intent(in)    :: p_model
          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
          type(sd_operators_x_oneside_L0), intent(in)    :: s_x_L0
          type(sd_operators_x_oneside_L1), intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1), intent(in)    :: s_x_R1
          type(sd_operators_x_oneside_R0), intent(in)    :: s_x_R0
          real(rkind)                    , intent(in)    :: dx
          real(rkind)                    , intent(in)    :: dy
          integer(ikind)                 , intent(in)    :: j_min
          integer(ikind)                 , intent(in)    :: j_max
          integer(ikind)                 , intent(in)    :: i
          integer                        , intent(in)    :: edge_card_coord
          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_y

          integer(ikind)        :: i_f
          integer(ikind)        :: j
          integer, dimension(4) :: bc_s
          

          !allocate space for temporary array: this%edge_inviscid_flux_y
          !which is needed for the computation of the transverse LODI
          !terms
          if(allocated(this%transverse_lodi)) then
             print '(''bc_operators_yoolodato_class'')'
             print '(''compute_fluxes_for_bc_x_edge'')'
             stop 'this%transverse_lodi already allocated'
          else
             allocate(this%transverse_lodi(2,j_max-j_min+1,ne))
          end if

          !allocate space for temporary array: this%edge_viscid_flux_y
          !which is needed for the computation of the viscous LODI terms
          if(allocated(this%viscous_lodi)) then
             print '(''bc_operators_yoolodato_class'')'
             print '(''compute_fluxes_for_bc_x_edge'')'
             stop 'this%edge_viscid_flux_y already allocated'
          else
             allocate(this%viscous_lodi(2,j_max-j_min+1,ne))
          end if

          
          this%i_edge_offset = i
          this%j_edge_offset = j_min


          select case(edge_card_coord)
            case(W)

               do j=j_min,j_max

                  !compute the inviscid and viscid fluxes
                  flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $                 nodes,i,j,dx,dy,s_x_L0,
     $                 this%transverse_lodi(1,j-j_min+1,:),
     $                 this%viscous_lodi(1,j-j_min+1,:))


                  !compute the inviscid and viscid fluxes
                  flux_y(i+1,j,:) = p_model%compute_flux_y_by_parts(
     $                 nodes,i+1,j,dx,dy,s_x_L1,
     $                 this%transverse_lodi(2,j-j_min+1,:),
     $                 this%viscous_lodi(2,j-j_min+1,:))

               end do
               
            case(E)

               do j=j_min,j_max

                  !compute the inviscid and viscid fluxes
                  flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $                 nodes,i,j,dx,dy,s_x_R1,
     $                 this%transverse_lodi(1,j-j_min+1,:),
     $                 this%viscous_lodi(1,j-j_min+1,:))

                  !compute the inviscid and viscid fluxes
                  flux_y(i+1,j,:) = p_model%compute_flux_y_by_parts(
     $                 nodes,i+1,j,dx,dy,s_x_R0,
     $                 this%transverse_lodi(2,j-j_min+1,:),
     $                 this%viscous_lodi(2,j-j_min+1,:))

               end do

            case default
               print '(''bc_operators_yoolodato_class'')'
               print '(''compute_fluxes_for_bc_x_edge'')'
               stop 'case not recognized'
          end select


          !compute the transverse_lodi and the viscous_lodi
          !from the inviscid and viscid fluxes previously
          !saved in 
          call compute_lodi_terms_x_edge(
     $         p_model,
     $         nodes,
     $         dy,
     $         this%i_edge_offset,
     $         this%j_edge_offset,
     $         this%transverse_lodi,
     $         this%viscous_lodi,
     $         this%transverse_lodi,
     $         this%viscous_lodi)          

        end subroutine compute_fluxes_for_bc_x_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the fluxes at the edge of the
        !> computational domain in the x-direction so that
        !> the time derivatives for an edge in the y-direction
        !> can be computed
        !
        !> @date
        !> 10_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param i_min
        !> index min along the x-direction corresponding
        !> to the beginning of the edge layer computed
        !
        !>@param i_max
        !> index max along the x-direction corresponding
        !> to the end of the edge layer computed
        !
        !>@param j
        !> index along the y-direction positioning the
        !> the edge boundary layer
        !
        !>@param edge_card_coord
        !> cardinal coordinate identifying the type of
        !> edge boundary layer
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !-------------------------------------------------------------
        subroutine compute_fluxes_for_bc_y_edge(
     $     this,
     $     p_model,
     $     nodes,
     $     s_y_L0, s_y_L1,
     $     s_y_R1, s_y_R0,
     $     dx, dy,
     $     i_min, i_max, j,
     $     edge_card_coord,
     $     flux_x)
        
          implicit none
        
          class(bc_operators_yoolodato)  , intent(inout) :: this
          type(pmodel_eq)                , intent(in)    :: p_model
          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
          type(sd_operators_y_oneside_L0), intent(in)    :: s_y_L0
          type(sd_operators_y_oneside_L1), intent(in)    :: s_y_L1
          type(sd_operators_y_oneside_R1), intent(in)    :: s_y_R1
          type(sd_operators_y_oneside_R0), intent(in)    :: s_y_R0
          real(rkind)                    , intent(in)    :: dx
          real(rkind)                    , intent(in)    :: dy
          integer(ikind)                 , intent(in)    :: i_min
          integer(ikind)                 , intent(in)    :: i_max
          integer(ikind)                 , intent(in)    :: j
          integer                        , intent(in)    :: edge_card_coord
          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_x

          integer(ikind) :: i

          
          !allocate space for temporary array: this%transverse_lodi
          !which is needed for the computation of the transverse LODI
          !terms
          if(allocated(this%transverse_lodi)) then
             print '(''bc_operators_yoolodato_class'')'
             print '(''compute_fluxes_for_bc_y_edge'')'
             stop 'this%transverse_lodi already allocated'
          else
             allocate(this%transverse_lodi(i_max-i_min+1,2,ne))
          end if

          !allocate space for temporary array: this%viscous_lodi
          !which is needed for the computation of the viscous LODI terms
          if(allocated(this%viscous_lodi)) then
             print '(''bc_operators_yoolodato_class'')'
             print '(''compute_fluxes_for_bc_y_edge'')'
             stop 'this%viscous_lodi already allocated'
          else
             allocate(this%viscous_lodi(i_max-i_min+1,2,ne))
          end if


          this%i_edge_offset = i_min
          this%j_edge_offset = j


          select case(edge_card_coord)
            case(S)               

               do i=i_min, i_max

                  flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $                 nodes,i,j,dx,dy,s_y_L0,
     $                 this%transverse_lodi(i-i_min+1,1,:),
     $                 this%viscous_lodi(i-i_min+1,1,:))

               end do

               do i=i_min, i_max

                  flux_x(i,j+1,:) = p_model%compute_flux_x_by_parts(
     $                 nodes,i,j+1,dx,dy,s_y_L1,
     $                 this%transverse_lodi(i-i_min+1,2,:),
     $                 this%viscous_lodi(i-i_min+1,2,:))

               end do

            case(N)

               do i=i_min,i_max

                  flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $                 nodes,i,j,dx,dy,s_y_R1,
     $                 this%transverse_lodi(i-i_min+1,1,:),
     $                 this%viscous_lodi(i-i_min+1,1,:))

               end do
          
               do i=i_min, i_max

                  flux_x(i,j+1,:) = p_model%compute_flux_x_by_parts(
     $                 nodes,i,j+1,dx,dy,s_y_R0,
     $                 this%transverse_lodi(i-i_min+1,2,:),
     $                 this%viscous_lodi(i-i_min+1,2,:))

               end do

            case default
               print '(''bc_operators_openbc_class'')'
               print '(''compute_fluxes_for_bc_y_edge'')'
               stop 'case not recognized'
          end select

        end subroutine compute_fluxes_for_bc_y_edge


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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> and x edge: W_edge or E_edge
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
        
          class(bc_operators_yoolodato), intent(in) :: this
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


          !compute the time derivatives from the LODI terms
          if(side_x.eqv.left) then

             call compute_timedev_x_edge(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_y,
     $            this%transverse_lodi(i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
     $            this%viscous_lodi(   i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
     $            side_x,
     $            gradient_x,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_W,
     $            timedev)

          else

             call compute_timedev_x_edge(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_y,
     $            this%transverse_lodi(i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
     $            this%viscous_lodi(   i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
     $            side_x,
     $            gradient_x,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_E,
     $            timedev)

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
        function apply_bc_on_timedev_y_edge(
     $     this,
     $     p_model, t,
     $     nodes, x_map, y_map, i,j,
     $     flux_x, side_y, gradient_y)
     $     result(timedev)
        
          implicit none
        
          class(bc_operators_yoolodato), intent(in) :: this
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


          !compute the time derivatives from the LODI terms
          if(side_x.eqv.left) then

             call compute_timedev_y_edge(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_y,
     $            this%transverse_lodi(i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
     $            this%viscous_lodi(   i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
     $            side_y,
     $            gradient_y,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_S,
     $            timedev)

          else

             call compute_timedev_y_edge(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_y,
     $            this%transverse_lodi(i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
     $            this%viscous_lodi(   i-this%i_edge_offset+1,j-this%j_edge_offset+1,:),
     $            side_y,
     $            gradient_y,
     $            this%edge_inflow_bc,
     $            this%edge_outflow_bc,
     $            this%flow_user_N,
     $            timedev)

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
        
          class(bc_operators_yoolodato), intent(in) :: this
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
        
        end function apply_bc_on_timedev_xy_corner

      end module bc_operators_class

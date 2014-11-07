      module bc_operators_yoolodato_class

        implicit none

        use bc_operators_openbc_class, only : bc_operators_openbc


        type, abstract, extends(bc_operators_openbc) :: bc_operators_yoolodato

          real(rkind), dimension(:,:,:), allocatable :: edge_inviscid_flux_x
          real(rkind), dimension(:,:,:), allocatable :: edge_inviscid_flux_y
          real(rkind), dimension(:,:,:), allocatable :: edge_viscid_flux_x
          real(rkind), dimension(:,:,:), allocatable :: edge_viscid_flux_y

          contains

          procedure, pass :: compute_fluxes_for_bc_x_edge
          procedure, pass :: compute_fluxes_for_bc_y_edge
        
        end type bc_operators_yoolodato


        contains

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
          if(allocated(this%edge_inviscid_flux_y)) then
             print '(''bc_operators_yoolodato_class'')'
             print '(''compute_fluxes_for_bc_x_edge'')'
             stop 'this%edge_inviscid_flux_y already allocated'
          else
             allocate(this%edge_inviscid_flux_y(2,j_max-j_min+1,ne))
          end if

          !allocate space for temporary array: this%edge_viscid_flux_y
          !which is needed for the computation of the viscous LODI terms
          if(allocated(this%edge_viscid_flux_y)) then
             print '(''bc_operators_yoolodato_class'')'
             print '(''compute_fluxes_for_bc_x_edge'')'
             stop 'this%edge_viscid_flux_y already allocated'
          else
             allocate(this%edge_viscid_flux_y(2,j_max-j_min+1,ne))
          end if


          select case(edge_card_coord)
            case(W)

               do j=j_min,j_max

                  !compute the inviscid and viscid fluxes
                  flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $                 nodes,i,j,dx,dy,s_x_L0,
     $                 this%edge_inviscid_flux_y(1,j-j_min+1,:),
     $                   this%edge_viscid_flux_y(1,j-j_min+1,:))


                  !compute the inviscid and viscid fluxes
                  flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $                 nodes,i+1,j,dx,dy,s_x_L1,
     $                 this%edge_inviscid_flux_y(2,j-j_min+1,:),
     $                   this%edge_viscid_flux_y(2,j-j_min+1,:))

               end do
               
            case(E)

               do j=j_min,j_max

                  !compute the inviscid and viscid fluxes
                  flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $                 nodes,i,j,dx,dy,s_x_R1,
     $                 this%edge_inviscid_flux_y(1,j-j_min+1,:),
     $                   this%edge_viscid_flux_y(1,j-j_min+1,:))

                  !compute the inviscid and viscid fluxes
                  flux_y(i,j,:) = p_model%compute_flux_y_by_parts(
     $                 nodes,i+1,j,dx,dy,s_x_R0,
     $                 this%edge_inviscid_flux_y(2,j-j_min+1,:),
     $                   this%edge_viscid_flux_y(2,j-j_min+1,:))

               end do

            case default
               print '(''bc_operators_yoolodato_class'')'
               print '(''compute_fluxes_for_bc_x_edge'')'
               stop 'case not recognized'
          end select

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

          
          !allocate space for temporary array: this%edge_inviscid_flux_y
          !which is needed for the computation of the transverse LODI
          !terms
          if(allocated(this%edge_inviscid_flux_x)) then
             print '(''bc_operators_yoolodato_class'')'
             print '(''compute_fluxes_for_bc_y_edge'')'
             stop 'this%edge_inviscid_flux_x already allocated'
          else
             allocate(this%edge_inviscid_flux_x(i_max-i_min+1,2,ne))
          end if

          !allocate space for temporary array: this%edge_viscid_flux_y
          !which is needed for the computation of the viscous LODI terms
          if(allocated(this%edge_viscid_flux_x)) then
             print '(''bc_operators_yoolodato_class'')'
             print '(''compute_fluxes_for_bc_y_edge'')'
             stop 'this%edge_viscid_flux_x already allocated'
          else
             allocate(this%edge_viscid_flux_x(i_max-i_min+1,2,ne))
          end if


          select case(edge_card_coord)
            case(S)

               do i=i_min, i_max

                  flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $                 nodes,i,j,dx,dy,s_y_L0,
     $                 this%edge_inviscid_flux_x(i-i_min+1,1,:),
     $                   this%edge_viscid_flux_x(i-i_min+1,1,:))

               end do

               do i=i_min, i_max

                  flux_x(i,j+1,:) = p_model%compute_flux_x_by_parts(
     $                 nodes,i,j+1,dx,dy,s_y_L1,
     $                 this%edge_inviscid_flux_x(i-i_min+1,2,:),
     $                   this%edge_viscid_flux_x(i-i_min+1,2,:))

               end do

            case(N)

               do i=i_min,i_max

                  flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $                 nodes,i,j,dx,dy,s_y_R1,
     $                 this%edge_inviscid_flux_x(i-i_min+1,1,:),
     $                   this%edge_viscid_flux_x(i-i_min+1,1,:))

               end do
          
               do i=i_min, i_max

                  flux_x(i,j,:) = p_model%compute_flux_x_by_parts(
     $                 nodes,i,j+1,dx,dy,s_y_R0,
     $                 this%edge_inviscid_flux_x(i-i_min+1,2,:),
     $                   this%edge_viscid_flux_x(i-i_min+1,2,:))

               end do

            case default
               print '(''bc_operators_openbc_class'')'
               print '(''compute_fluxes_for_bc_y_edge'')'
               stop 'case not recognized'
          end select

        end subroutine compute_fluxes_for_bc_y_edge

      end module bc_operators_yoolodato_class

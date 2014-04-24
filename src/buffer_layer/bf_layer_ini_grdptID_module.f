      !> @file
      !> module encapsulating the subroutines related to initialization
      !> of the grdpt ID
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> subroutines related to initialization of the grdpt ID
      !
      !> @date
      ! 07_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_ini_grdptID_module

        use parameters_bf_layer, only : interior_pt,
     $                                  bc_interior_pt, bc_pt,
     $                                  exchange_pt
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input   , only : bc_size
        use parameters_kind    , only : ikind

        private
        public :: ini_grdpts_id_N,
     $            ini_grdpts_id_S,
     $            ini_grdpts_id_E,
     $            ini_grdpts_id_W,
     $            ini_grdpts_id_NE,
     $            ini_grdpts_id_NW,
     $            ini_grdpts_id_SE,
     $            ini_grdpts_id_SW

        contains

        
        !< initialize the grdpts_id for the northern buffer layer
        subroutine ini_grdpts_id_N(grdpts_id, neighbors)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id
          logical, dimension(4)  , intent(in)  :: neighbors

          call add_exchange_layer_NS(
     $         grdpts_id, 1, bc_size)

          call add_interior_layer_NS(
     $         grdpts_id, bc_size+1, neighbors)

          call add_bc_interior_layer_NS(
     $         grdpts_id, bc_size+2, neighbors)

          call add_bc_layer_NS(
     $         grdpts_id, bc_size+3, size(grdpts_id,2), neighbors)

        end subroutine ini_grdpts_id_N


        !< initialize the grdpts_id for the southern buffer layer
        subroutine ini_grdpts_id_S(grdpts_id, neighbors)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id
          logical, dimension(4)  , intent(in)  :: neighbors

          call add_bc_layer_NS(
     $         grdpts_id, 1, size(grdpts_id,2)-(bc_size+2), neighbors)

          call add_bc_interior_layer_NS(
     $         grdpts_id, size(grdpts_id,2)-(bc_size+2)+1, neighbors)

          call add_interior_layer_NS(
     $         grdpts_id, size(grdpts_id,2)-(bc_size+2)+2, neighbors)

          call add_exchange_layer_NS(
     $         grdpts_id, size(grdpts_id,2)-(bc_size+2)+3, size(grdpts_id,2))

        end subroutine ini_grdpts_id_S

        !< initialize the grdpts_id for the eastern buffer layer
        subroutine ini_grdpts_id_E(grdpts_id, neighbors)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id
          logical, dimension(4)  , intent(in)  :: neighbors

          call add_edge_layer_E(
     $         grdpts_id, 1, bc_size, neighbors, .false.)

          call add_interior_layer_E(
     $         grdpts_id, bc_size+1, size(grdpts_id,2)-bc_size)

          call add_edge_layer_E(
     $         grdpts_id, size(grdpts_id,2)-bc_size+1, size(grdpts_id,2), neighbors, .true.)
          
        end subroutine ini_grdpts_id_E


        !< initialize the grdpts_id for the eastern buffer layer
        subroutine ini_grdpts_id_W(grdpts_id, neighbors)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id
          logical, dimension(4)  , intent(in)  :: neighbors

          call add_edge_layer_W(
     $         grdpts_id,
     $         1, bc_size,
     $         neighbors,
     $         .false.)

          call add_interior_layer_W(
     $         grdpts_id, bc_size+1, size(grdpts_id,2)-bc_size)

          call add_edge_layer_W(
     $         grdpts_id,
     $         size(grdpts_id,2)-bc_size+1, size(grdpts_id,2),
     $         neighbors,
     $         .true.)
          
        end subroutine ini_grdpts_id_W

      
        !< initialize the grdpts_id for the eastern buffer layer
        subroutine ini_grdpts_id_NE(grdpts_id)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id

          integer(ikind) :: i,j

          !exchange
          do j=1, bc_size
             do i=1, 2*bc_size+1
                grdpts_id(i,j) = exchange_pt
             end do
          end do

          !interior
          j=bc_size+1

          do i=1, bc_size
             grdpts_id(i,j) = exchange_pt
          end do

          i=bc_size+1
          grdpts_id(i,j)= interior_pt

          i=bc_size+2
          grdpts_id(i,j)= bc_interior_pt

          do i=bc_size+3, size(grdpts_id,1)
             grdpts_id(i,j)= bc_pt
          end do

          !bc_interior
          j=bc_size+2

          do i=1, bc_size
             grdpts_id(i,j) = exchange_pt
          end do

          do i=bc_size+1, size(grdpts_id,1)-bc_size+1
             grdpts_id(i,j)= bc_interior_pt
          end do

          do i=size(grdpts_id,1)-bc_size+2, size(grdpts_id,1)
             grdpts_id(i,j)= bc_pt
          end do          

          !bc
          do j=bc_size+3, size(grdpts_id,2)
             do i=1, bc_size
                grdpts_id(i,j) = exchange_pt
             end do

             do i=bc_size+1, size(grdpts_id,1)
                grdpts_id(i,j)= bc_pt
             end do             
          end do
             
        end subroutine ini_grdpts_id_NE


        !< initialize the grdpts_id for the eastern buffer layer
        subroutine ini_grdpts_id_NW(grdpts_id)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id

          integer(ikind) :: i,j

          !exchange
          do j=1, bc_size
             do i=1, 2*bc_size+1
                grdpts_id(i,j) = exchange_pt
             end do
          end do

          !interior
          j=bc_size+1
          do i=1, bc_size-1
             grdpts_id(i,j)= bc_pt
          end do

          i=bc_size
          grdpts_id(i,j)= bc_interior_pt

          i=bc_size+1
          grdpts_id(i,j)= interior_pt

          do i=bc_size+2, size(grdpts_id,1)
             grdpts_id(i,j) = exchange_pt
          end do

          !bc_interior
          j=bc_size+2
          do i=1, bc_size-1
             grdpts_id(i,j)= bc_pt
          end do

          do i=bc_size, bc_size+1
             grdpts_id(i,j)= bc_interior_pt
          end do

          do i=bc_size+2, size(grdpts_id,1)
             grdpts_id(i,j) = exchange_pt
          end do

          !bc
          do j=bc_size+3, size(grdpts_id,2)
             do i=1, size(grdpts_id,1)-bc_size
                grdpts_id(i,j)= bc_pt
             end do
             
             do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                grdpts_id(i,j) = exchange_pt
             end do
          end do
             
        end subroutine ini_grdpts_id_NW


        !< initialize the grdpts_id for the eastern buffer layer
        subroutine ini_grdpts_id_SE(grdpts_id)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id

          integer(ikind) :: i,j

          !bc
          do j=1, size(grdpts_id,2)-(bc_size+2)
             do i=1, bc_size
                grdpts_id(i,j) = exchange_pt
             end do

             do i=bc_size+1, size(grdpts_id,1)
                grdpts_id(i,j)= bc_pt
             end do             
          end do

          !bc_interior
          j=size(grdpts_id,2)-(bc_size+2)+1

          do i=1, bc_size
             grdpts_id(i,j) = exchange_pt
          end do

          do i=bc_size+1, size(grdpts_id,1)-bc_size+1
             grdpts_id(i,j)= bc_interior_pt
          end do

          do i=size(grdpts_id,1)-bc_size+2, size(grdpts_id,1)
             grdpts_id(i,j)= bc_pt
          end do 

          !interior
          j=size(grdpts_id,2)-(bc_size+2)+2

          do i=1, bc_size
             grdpts_id(i,j) = exchange_pt
          end do

          i=bc_size+1
          grdpts_id(i,j)= interior_pt

          i=bc_size+2
          grdpts_id(i,j)= bc_interior_pt

          do i=bc_size+3, size(grdpts_id,1)
             grdpts_id(i,j)= bc_pt
          end do

          !exchange
          do j=size(grdpts_id,2)-(bc_size+2)+3, size(grdpts_id,2)
             do i=1, 2*bc_size+1
                grdpts_id(i,j) = exchange_pt
             end do
          end do          
             
        end subroutine ini_grdpts_id_SE


        !< initialize the grdpts_id for the eastern buffer layer
        subroutine ini_grdpts_id_SW(grdpts_id)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id

          integer(ikind) :: i,j

          !bc
          do j=1, size(grdpts_id,2)-(bc_size+2)
             do i=1, size(grdpts_id,1)-bc_size
                grdpts_id(i,j)= bc_pt
             end do
             
             do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                grdpts_id(i,j) = exchange_pt
             end do
          end do

          !bc_interior
          j=size(grdpts_id,2)-(bc_size+2)+1
          do i=1, bc_size-1
             grdpts_id(i,j)= bc_pt
          end do

          do i=bc_size, bc_size+1
             grdpts_id(i,j)= bc_interior_pt
          end do

          do i=bc_size+2, size(grdpts_id,1)
             grdpts_id(i,j) = exchange_pt
          end do

          !interior
          j=size(grdpts_id,2)-(bc_size+2)+2
          do i=1, bc_size-1
             grdpts_id(i,j)= bc_pt
          end do

          i=bc_size
          grdpts_id(i,j)= bc_interior_pt

          i=bc_size+1
          grdpts_id(i,j)= interior_pt

          do i=bc_size+2, size(grdpts_id,1)
             grdpts_id(i,j) = exchange_pt
          end do

          !exchange
          do j=size(grdpts_id,2)-(bc_size+2)+3, size(grdpts_id,2)
             do i=1, 2*bc_size+1
                grdpts_id(i,j) = exchange_pt
             end do
          end do     
             
        end subroutine ini_grdpts_id_SW


        !< add the interior layer for the initialization
        !> of gridpts_id for eastern buffer layers
        subroutine add_interior_layer_E(
     $     grdpts_id, j_min, j_max)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: j_min
          integer(ikind)         , intent(in)    :: j_max

          integer(ikind) :: i,j

          do j=j_min, j_max
             do i=1, bc_size
                grdpts_id(i,j) = exchange_pt
             end do

             do i=bc_size+1, size(grdpts_id,1)-bc_size
                grdpts_id(i,j) = interior_pt
             end do

             i=size(grdpts_id,1)-bc_size+1
             grdpts_id(i,j) = bc_interior_pt

             do i=size(grdpts_id,1)-bc_size+2, size(grdpts_id,1)
                grdpts_id(i,j) = bc_pt
             end do
          end do          

        end subroutine add_interior_layer_E


        !< add the interior layer for the initialization
        !> of gridpts_id for western buffer layers
        subroutine add_interior_layer_W(
     $     grdpts_id, j_min, j_max)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: j_min
          integer(ikind)         , intent(in)    :: j_max

          integer(ikind) :: i,j

          do j=j_min, j_max
             do i=1, bc_size-1
                grdpts_id(i,j) = bc_pt
             end do

             i=bc_size
             grdpts_id(i,j) = bc_interior_pt

             do i=bc_size+1, size(grdpts_id,1)-bc_size+1
                grdpts_id(i,j) = interior_pt
             end do
             
             do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                grdpts_id(i,j) = exchange_pt
             end do
             
          end do          

        end subroutine add_interior_layer_W


        !< add the bc_pt layer for the initialization
        !> of gridpts_id
        !> type=.true.  = upper layer
        !> type=.false. = lower layer
        subroutine add_edge_layer_E(
     $     grdpts_id, j_min, j_max, neighbors, type)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: j_min
          integer(ikind)         , intent(in)    :: j_max
          logical, dimension(4)  , intent(in)    :: neighbors
          logical                , intent(in)    :: type

          integer(ikind) :: i,j

          if(type) then
             if(neighbors(N)) then
                do j=j_min, j_max
                   do i=1, size(grdpts_id,1)
                      grdpts_id(i,j) = exchange_pt
                   end do
                end do
             else
                j=j_min
                do i=1, bc_size
                   grdpts_id(i,j) = exchange_pt
                end do
                
                do i=bc_size+1, size(grdpts_id,1)-bc_size+1
                   grdpts_id(i,j) = bc_interior_pt
                end do
                
                do i=size(grdpts_id,1)-bc_size+2, size(grdpts_id,1)
                   grdpts_id(i,j) = bc_pt
                end do
                
                do j=j_min+1, j_max
                   do i=1, bc_size
                      grdpts_id(i,j) = exchange_pt
                   end do
                   
                   do i=bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_pt
                   end do
                end do
             end if
          else
             if(neighbors(S)) then
                do j=j_min, j_max
                   do i=1, size(grdpts_id,1)
                      grdpts_id(i,j) = exchange_pt
                   end do
                end do
             else                
                do j=j_min, j_min+bc_size-2
                   do i=1, bc_size
                      grdpts_id(i,j) = exchange_pt
                   end do
                   
                   do i=bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_pt
                   end do
                end do

                do j=j_min+bc_size-1, j_max
                   do i=1, bc_size
                      grdpts_id(i,j) = exchange_pt
                   end do
                   
                   do i=bc_size+1, size(grdpts_id,1)-bc_size+1
                      grdpts_id(i,j) = bc_interior_pt
                   end do
                   
                   do i=size(grdpts_id,1)-bc_size+2, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_pt
                   end do
                end do
             end if
          end if

        end subroutine add_edge_layer_E


        !< add the bc_pt layer for the initialization
        !> of gridpts_id
        !> type=.true.  = upper layer
        !> type=.false. = lower layer
        subroutine add_edge_layer_W(
     $     grdpts_id, j_min, j_max, neighbors, type)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: j_min
          integer(ikind)         , intent(in)    :: j_max
          logical, dimension(4)  , intent(in)    :: neighbors
          logical                , intent(in)    :: type

          integer(ikind) :: i,j

          if(type) then
             if(neighbors(N)) then
                do j=j_min, j_max
                   do i=1, size(grdpts_id,1)
                      grdpts_id(i,j) = exchange_pt
                   end do
                end do
             else
                j=j_min
                do i=1, bc_size-1
                   grdpts_id(i,j) = bc_pt
                end do
                
                do i=bc_size, size(grdpts_id,1)-bc_size
                   grdpts_id(i,j) = bc_interior_pt
                end do

                do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                   grdpts_id(i,j) = exchange_pt
                end do

                do j=j_min+1, j_max
                   do i=1, size(grdpts_id,1)-bc_size
                      grdpts_id(i,j) = bc_pt
                   end do
                   
                   do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = exchange_pt
                   end do             
                end do
             end if

          else
             if(neighbors(S)) then
                do j=j_min, j_max
                   do i=1, size(grdpts_id,1)
                      grdpts_id(i,j) = exchange_pt
                   end do
                end do
             else
                do j=j_min, j_min+bc_size-2
                   do i=1, size(grdpts_id,1)-bc_size
                      grdpts_id(i,j) = bc_pt
                   end do
                   
                   do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = exchange_pt
                   end do             
                end do

                do j=j_min+bc_size-1, j_max
                   do i=1, bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=bc_size, size(grdpts_id,1)-bc_size
                      grdpts_id(i,j) = bc_interior_pt
                   end do

                   do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = exchange_pt
                   end do
                end do
             end if
          end if          

        end subroutine add_edge_layer_W


        !< add the bc_pt layer for the initialization
        !> of gridpts_id
        subroutine add_bc_layer_NS(
     $     grdpts_id, j_min, j_max, neighbors)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: j_min
          integer(ikind)         , intent(in)    :: j_max
          logical, dimension(4)  , intent(in)    :: neighbors

          integer(ikind) :: i,j

          do j=j_min, j_max
             if(neighbors(W)) then
                do i=1, bc_size
                   grdpts_id(i,j) = exchange_pt
                end do
             else
                do i=1, bc_size
                   grdpts_id(i,j) = bc_pt
                end do
             end if

             do i=bc_size+1, size(grdpts_id,1)-bc_size
                grdpts_id(i,j) = bc_pt
             end do

             if(neighbors(E)) then
                do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                   grdpts_id(i,j) = exchange_pt
                end do
             else
                do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                   grdpts_id(i,j) = bc_pt
                end do
             end if             
          end do

        end subroutine add_bc_layer_NS

      
        subroutine add_bc_interior_layer_NS(
     $     grdpts_id, j, neighbors)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: j
          logical, dimension(4)  , intent(in)    :: neighbors

          integer(ikind) :: i

          if(neighbors(W)) then
             do i=1,bc_size
                grdpts_id(i,j) = exchange_pt
             end do
          else
             do i=1, bc_size-1
                grdpts_id(i,j) = bc_pt
             end do
             
             i=bc_size
             grdpts_id(i,j) = bc_interior_pt
          end if

          do i=bc_size+1, size(grdpts_id,1)-bc_size
             grdpts_id(i,j) = bc_interior_pt
          end do

          if(neighbors(E)) then
             do i=size(grdpts_id,1)-bc_size+1,size(grdpts_id,1)
                grdpts_id(i,j) = exchange_pt
             end do
          else
             i=size(grdpts_id,1)-bc_size+1
             grdpts_id(i,j) = bc_interior_pt
             
             do i=size(grdpts_id,1)-bc_size+2,size(grdpts_id,1)
                grdpts_id(i,j) = bc_pt
             end do
          end if

        end subroutine add_bc_interior_layer_NS


        !< add the interior layer for the initialization
        !> of the gridpts_id
        subroutine add_interior_layer_NS(
     $     grdpts_id, j, neighbors)
        
          implicit none
          
          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: j
          logical, dimension(4)  , intent(in)    :: neighbors

          integer(ikind) :: i
          
          if(neighbors(W)) then
             do i=1, bc_size
                grdpts_id(i,j) = exchange_pt
             end do
          else
             do i=1, bc_size-1
                grdpts_id(i,j) = bc_pt
             end do
             
             i=bc_size
             grdpts_id(i,j) = bc_interior_pt
          end if

          do i=bc_size+1, size(grdpts_id,1)-bc_size
             grdpts_id(i,j) = interior_pt
          end do

          if(neighbors(E)) then
             do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                grdpts_id(i,j) = exchange_pt
             end do
          else
             i=size(grdpts_id,1)-bc_size+1
             grdpts_id(i,j) = bc_interior_pt
             
             do i=size(grdpts_id,1)-bc_size+2, size(grdpts_id,1)
                grdpts_id(i,j) = bc_pt
             end do
          end if

        end subroutine add_interior_layer_NS


        !< add the exchnage layer for the initialization
        !> of the gridpts_id
        subroutine add_exchange_layer_NS(
     $     grdpts_id, j_min, j_max)
        
          implicit none
          
          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: j_min, j_max

          integer(ikind) :: i,j

          do j=j_min, j_max
             do i=1, size(grdpts_id,1)
                grdpts_id(i,j) = exchange_pt
             end do
          end do

        end subroutine add_exchange_layer_NS
      

c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> subroutine identifying the grid points of the
c$$$        !> buffer layer after the allocation
c$$$        !
c$$$        !> @date
c$$$        !> 07_04_2013 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> bf_layer_abstract object encapsulating the main
c$$$        !> tables and the integer identifying the
c$$$        !> correspondance between the buffer layer and the
c$$$        !> interior grid points
c$$$        !--------------------------------------------------------------
c$$$        subroutine ini_grdptID(this, neighbors)
c$$$        
c$$$          implicit none
c$$$
c$$$          class(bf_layer_abstract), intent(inout) :: this
c$$$          logical, dimension(4)   , intent(in)    :: neighbors
c$$$
c$$$          integer :: i,j
c$$$
c$$$          do j=1, size(this%grdpts_id,2)
c$$$             do i=1, size(this%grdpts_id,1)
c$$$                this%grdpts_id(i,j) = interior_pt
c$$$             end do
c$$$          end do        
c$$$
c$$$          select case(this%localization)
c$$$            case(N_W)
c$$$               call make_west_layer(this%grdpts_id)
c$$$               call make_north_layer(this%grdpts_id)
c$$$               call add_NE_corner(this%grdpts_id)
c$$$
c$$$               call add_east_layer_exchange(this%grdpts_id)
c$$$               call add_south_layer_exchange(this%grdpts_id)
c$$$
c$$$            case(N)
c$$$               call make_west_layer(this%grdpts_id)
c$$$               call make_east_layer(this%grdpts_id)
c$$$               call make_north_layer(this%grdpts_id)
c$$$
c$$$               call add_south_layer_exchange(this%grdpts_id)
c$$$
c$$$               if(neighbors(E)) then
c$$$                  call add_east_layer_exchange(this%grdpts_id)
c$$$               end if
c$$$               if(neighbors(W)) then
c$$$                  call add_west_layer_exchange(this%grdpts_id)
c$$$               end if
c$$$
c$$$            case(N_E)
c$$$               call make_east_layer(this%grdpts_id)
c$$$               call add_NW_corner(this%grdpts_id)
c$$$               call make_north_layer(this%grdpts_id)
c$$$
c$$$               call add_west_layer_exchange(this%grdpts_id)
c$$$               call add_south_layer_exchange(this%grdpts_id)
c$$$
c$$$            case(W)
c$$$               call make_west_layer(this%grdpts_id)
c$$$               call make_south_layer(this%grdpts_id)
c$$$               call make_north_layer(this%grdpts_id)
c$$$               call add_NE_corner(this%grdpts_id)
c$$$               call add_SE_corner(this%grdpts_id)
c$$$
c$$$               call add_east_layer_exchange(this%grdpts_id)
c$$$
c$$$               if(neighbors(N)) then
c$$$                  call add_north_layer_exchange(this%grdpts_id)
c$$$               end if
c$$$               if(neighbors(S)) then
c$$$                  call add_south_layer_exchange(this%grdpts_id)
c$$$               end if
c$$$
c$$$            case(E)
c$$$               call make_east_layer(this%grdpts_id)
c$$$               call make_south_layer(this%grdpts_id)
c$$$               call make_north_layer(this%grdpts_id)
c$$$               call add_NW_corner(this%grdpts_id)
c$$$               call add_SW_corner(this%grdpts_id)
c$$$
c$$$               call add_west_layer_exchange(this%grdpts_id)
c$$$
c$$$               if(neighbors(N)) then
c$$$                  call add_north_layer_exchange(this%grdpts_id)
c$$$               end if
c$$$               if(neighbors(S)) then
c$$$                  call add_south_layer_exchange(this%grdpts_id)
c$$$               end if
c$$$
c$$$            case(S_W)
c$$$               call make_west_layer(this%grdpts_id)
c$$$               call make_south_layer(this%grdpts_id)
c$$$               call add_SE_corner(this%grdpts_id)
c$$$
c$$$               call add_east_layer_exchange(this%grdpts_id)
c$$$               call add_north_layer_exchange(this%grdpts_id)
c$$$
c$$$            case(S)
c$$$               call make_west_layer(this%grdpts_id)
c$$$               call make_east_layer(this%grdpts_id)
c$$$               call make_south_layer(this%grdpts_id)
c$$$
c$$$               call add_north_layer_exchange(this%grdpts_id)
c$$$
c$$$               if(neighbors(E)) then
c$$$                  call add_east_layer_exchange(this%grdpts_id)
c$$$               end if
c$$$               if(neighbors(W)) then
c$$$                  call add_west_layer_exchange(this%grdpts_id)
c$$$               end if
c$$$
c$$$            case(S_E)
c$$$               call make_east_layer(this%grdpts_id)
c$$$               call make_south_layer(this%grdpts_id)
c$$$               call add_SW_corner(this%grdpts_id)
c$$$
c$$$               call add_west_layer_exchange(this%grdpts_id)
c$$$               call add_north_layer_exchange(this%grdpts_id)
c$$$
c$$$            case default
c$$$               print '(''bf_layer_ini_grdptID_module'')'
c$$$               print '(''ini_grdptID'')'
c$$$               print *, 'localization: ', this%localization
c$$$               stop 'localization not recognized'
c$$$
c$$$          end select
c$$$
c$$$        end subroutine ini_grdptID
c$$$
c$$$
c$$$        subroutine make_west_layer(grdpt_id)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          integer :: j
c$$$
c$$$          do j=1,size(grdpt_id,2)
c$$$             grdpt_id(1,j)=bc_pt
c$$$             grdpt_id(2,j)=bc_interior_pt
c$$$          end do
c$$$
c$$$        end subroutine make_west_layer
c$$$
c$$$
c$$$        subroutine make_east_layer(grdpt_id)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          integer :: j
c$$$
c$$$          do j=1,size(grdpt_id,2)
c$$$             grdpt_id(size(grdpt_id,1)-1,j)=bc_interior_pt
c$$$             grdpt_id(size(grdpt_id,1)  ,j)=bc_pt
c$$$          end do
c$$$
c$$$        end subroutine make_east_layer
c$$$
c$$$
c$$$        subroutine make_north_layer(grdpt_id)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          integer :: i
c$$$
c$$$          do i=bc_size,size(grdpt_id,1)-1
c$$$             grdpt_id(i,size(grdpt_id,2)-1) = bc_interior_pt
c$$$             grdpt_id(i,size(grdpt_id,2))   = bc_pt
c$$$          end do
c$$$
c$$$        end subroutine make_north_layer
c$$$
c$$$
c$$$        subroutine make_south_layer(grdpt_id)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          integer :: i
c$$$
c$$$          do i=bc_size,size(grdpt_id,1)-1
c$$$             grdpt_id(i,1)       = bc_pt
c$$$             grdpt_id(i,bc_size) = bc_interior_pt
c$$$          end do
c$$$        
c$$$        end subroutine make_south_layer
c$$$        
c$$$        subroutine add_NE_corner(grdpt_id)
c$$$        
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          grdpt_id(size(grdpt_id,1),size(grdpt_id,2)-1) = bc_interior_pt
c$$$          grdpt_id(size(grdpt_id,1),size(grdpt_id,2))   = bc_pt
c$$$
c$$$        end subroutine add_NE_corner
c$$$
c$$$
c$$$        subroutine add_NW_corner(grdpt_id)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          grdpt_id(1,size(grdpt_id,2)-1) = bc_interior_pt
c$$$          grdpt_id(1,size(grdpt_id,2))   = bc_pt
c$$$
c$$$        end subroutine add_NW_corner
c$$$
c$$$        subroutine add_SW_corner(grdpt_id)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          grdpt_id(1,1)   = bc_pt
c$$$          grdpt_id(1,bc_size) = bc_interior_pt
c$$$
c$$$        end subroutine add_SW_corner
c$$$      
c$$$        subroutine add_SE_corner(grdpt_id)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          grdpt_id(size(grdpt_id,1),1)       = bc_pt
c$$$          grdpt_id(size(grdpt_id,1),bc_size) = bc_interior_pt
c$$$
c$$$        end subroutine add_SE_corner
c$$$
c$$$
c$$$        subroutine add_north_layer_exchange(grdpt_id)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          integer :: i,j
c$$$
c$$$          do j=size(grdpt_id,2)-bc_size+1, size(grdpt_id,2)
c$$$             do i=1, size(grdpt_id,1)
c$$$                grdpt_id(i,j)=exchange_pt
c$$$             end do
c$$$          end do
c$$$
c$$$        end subroutine add_north_layer_exchange
c$$$
c$$$
c$$$        subroutine add_south_layer_exchange(grdpt_id)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          integer :: i,j
c$$$
c$$$          do j=1,bc_size
c$$$             do i=1, size(grdpt_id,1)
c$$$                grdpt_id(i,j)=exchange_pt
c$$$             end do
c$$$          end do
c$$$
c$$$        end subroutine add_south_layer_exchange
c$$$
c$$$
c$$$        subroutine add_east_layer_exchange(grdpt_id)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          integer :: i,j
c$$$
c$$$          do j=1,size(grdpt_id,2)
c$$$             do i=size(grdpt_id,1)-bc_size+1, size(grdpt_id,1)
c$$$                grdpt_id(i,j)=exchange_pt
c$$$             end do
c$$$          end do
c$$$
c$$$        end subroutine add_east_layer_exchange
c$$$
c$$$
c$$$        subroutine add_west_layer_exchange(grdpt_id)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(inout) :: grdpt_id
c$$$
c$$$          integer :: i,j
c$$$
c$$$          do j=1,size(grdpt_id,2)
c$$$             do i=1,bc_size
c$$$                grdpt_id(i,j)=exchange_pt
c$$$             end do
c$$$          end do
c$$$
c$$$        end subroutine add_west_layer_exchange

      end module bf_layer_ini_grdptID_module

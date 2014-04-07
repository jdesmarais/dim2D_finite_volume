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

        use bf_layer_abstract_class, only : bf_layer_abstract

        use parameters_bf_layer, only : no_pt, interior_pt,
     $                                  bc_interior_pt, bc_pt,
     $                                  exchange_pt
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input   , only : bc_size

        private
        public :: ini_grdptID


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine identifying the grid points of the
        !> buffer layer after the allocation
        !
        !> @date
        !> 07_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !--------------------------------------------------------------
        subroutine ini_grdptID(this, neighbors)
        
          implicit none

          class(bf_layer_abstract), intent(inout) :: this
          logical, dimension(4)   , intent(in)    :: neighbors

          integer :: i,j

          do j=1, size(this%grdpts_id,2)
             do i=1, size(this%grdpts_id,1)
                this%grdpts_id(i,j) = interior_pt
             end do
          end do        

          select case(this%localization(1))
            case(N_W)
               call make_west_layer(this%grdpts_id)
               call make_north_layer(this%grdpts_id)
               call add_NE_corner(this%grdpts_id)

               call add_east_layer_exchange(this%grdpts_id)
               call add_south_layer_exchange(this%grdpts_id)

            case(N)
               call make_west_layer(this%grdpts_id)
               call make_east_layer(this%grdpts_id)
               call make_north_layer(this%grdpts_id)

               call add_south_layer_exchange(this%grdpts_id)

               if(neighbors(E)) then
                  call add_east_layer_exchange(this%grdpts_id)
               end if
               if(neighbors(W)) then
                  call add_west_layer_exchange(this%grdpts_id)
               end if

            case(N_E)
               call make_east_layer(this%grdpts_id)
               call add_NW_corner(this%grdpts_id)
               call make_north_layer(this%grdpts_id)

               call add_west_layer_exchange(this%grdpts_id)
               call add_south_layer_exchange(this%grdpts_id)

            case(W)
               call make_west_layer(this%grdpts_id)
               call make_south_layer(this%grdpts_id)
               call make_north_layer(this%grdpts_id)
               call add_NE_corner(this%grdpts_id)
               call add_SE_corner(this%grdpts_id)

               call add_east_layer_exchange(this%grdpts_id)

               if(neighbors(N)) then
                  call add_north_layer_exchange(this%grdpts_id)
               end if
               if(neighbors(S)) then
                  call add_south_layer_exchange(this%grdpts_id)
               end if

            case(E)
               call make_east_layer(this%grdpts_id)
               call make_south_layer(this%grdpts_id)
               call make_north_layer(this%grdpts_id)
               call add_NW_corner(this%grdpts_id)
               call add_SW_corner(this%grdpts_id)

               call add_west_layer_exchange(this%grdpts_id)

               if(neighbors(N)) then
                  call add_north_layer_exchange(this%grdpts_id)
               end if
               if(neighbors(S)) then
                  call add_south_layer_exchange(this%grdpts_id)
               end if

            case(S_W)
               call make_west_layer(this%grdpts_id)
               call make_south_layer(this%grdpts_id)
               call add_SE_corner(this%grdpts_id)

               call add_east_layer_exchange(this%grdpts_id)
               call add_north_layer_exchange(this%grdpts_id)

            case(S)
               call make_west_layer(this%grdpts_id)
               call make_east_layer(this%grdpts_id)
               call make_south_layer(this%grdpts_id)

               call add_north_layer_exchange(this%grdpts_id)

               if(neighbors(E)) then
                  call add_east_layer_exchange(this%grdpts_id)
               end if
               if(neighbors(W)) then
                  call add_west_layer_exchange(this%grdpts_id)
               end if

            case(S_E)
               call make_east_layer(this%grdpts_id)
               call make_south_layer(this%grdpts_id)
               call add_SW_corner(this%grdpts_id)

               call add_west_layer_exchange(this%grdpts_id)
               call add_north_layer_exchange(this%grdpts_id)

            case default
               print '(''bf_layer_ini_grdptID_module'')'
               print '(''ini_grdptID'')'
               print *, 'localization: ', this%localization(1)
               stop 'localization not recognized'

          end select

        end subroutine ini_grdptID


        subroutine make_west_layer(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: j

          do j=1,size(grdpt_id,2)
             grdpt_id(1,j)=bc_pt
             grdpt_id(2,j)=bc_interior_pt
          end do

        end subroutine make_west_layer


        subroutine make_east_layer(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: j

          do j=1,size(grdpt_id,2)
             grdpt_id(size(grdpt_id,1)-1,j)=bc_interior_pt
             grdpt_id(size(grdpt_id,1)  ,j)=bc_pt
          end do

        end subroutine make_east_layer


        subroutine make_north_layer(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i

          do i=bc_size,size(grdpt_id,1)-1
             grdpt_id(i,size(grdpt_id,2)-1) = bc_interior_pt
             grdpt_id(i,size(grdpt_id,2))   = bc_pt
          end do

        end subroutine make_north_layer


        subroutine make_south_layer(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i

          do i=bc_size,size(grdpt_id,1)-1
             grdpt_id(i,1)       = bc_pt
             grdpt_id(i,bc_size) = bc_interior_pt
          end do
        
        end subroutine make_south_layer
        
        subroutine add_NE_corner(grdpt_id)
        
          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          grdpt_id(size(grdpt_id,1),size(grdpt_id,2)-1) = bc_interior_pt
          grdpt_id(size(grdpt_id,1),size(grdpt_id,2))   = bc_pt

        end subroutine add_NE_corner


        subroutine add_NW_corner(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          grdpt_id(1,size(grdpt_id,2)-1) = bc_interior_pt
          grdpt_id(1,size(grdpt_id,2))   = bc_pt

        end subroutine add_NW_corner

        subroutine add_SW_corner(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          grdpt_id(1,1)   = bc_pt
          grdpt_id(1,bc_size) = bc_interior_pt

        end subroutine add_SW_corner
      
        subroutine add_SE_corner(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          grdpt_id(size(grdpt_id,1),1)       = bc_pt
          grdpt_id(size(grdpt_id,1),bc_size) = bc_interior_pt

        end subroutine add_SE_corner


        subroutine add_north_layer_exchange(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i,j

          do j=size(grdpt_id,2)-bc_size+1, size(grdpt_id,2)
             do i=1, size(grdpt_id,1)
                grdpt_id(i,j)=exchange_pt
             end do
          end do

        end subroutine add_north_layer_exchange


        subroutine add_south_layer_exchange(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i,j

          do j=1,bc_size
             do i=1, size(grdpt_id,1)
                grdpt_id(i,j)=exchange_pt
             end do
          end do

        end subroutine add_south_layer_exchange


        subroutine add_east_layer_exchange(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i,j

          do j=1,size(grdpt_id,2)
             do i=size(grdpt_id,1)-bc_size+1, size(grdpt_id,1)
                grdpt_id(i,j)=exchange_pt
             end do
          end do

        end subroutine add_east_layer_exchange


        subroutine add_west_layer_exchange(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i,j

          do j=1,size(grdpt_id,2)
             do i=1,bc_size
                grdpt_id(i,j)=exchange_pt
             end do
          end do

        end subroutine add_west_layer_exchange

      end module bf_layer_ini_grdptID_module

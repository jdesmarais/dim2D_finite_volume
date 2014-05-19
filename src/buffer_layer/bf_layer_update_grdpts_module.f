      !> @file
      !> module encapsulating the subroutines needed to update
      !> the grid points of the buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the subroutines needed to update
      !> the grid points of the buffer layer
      !
      !> @date
      ! 07_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_update_grdpts_module

        use parameters_bf_layer, only : no_pt, interior_pt, bc_interior_pt, bc_pt
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input   , only : nx,ny,ne,bc_size
        use parameters_kind    , only : ikind, rkind

        implicit none

        
        private
        public :: update_grdpts


        contains
        

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine updating the grid points of the buffer
        !> layer: bc_interior_pt were identified by the
        !> detectors to be turned into interior pt. Therefore,
        !> it is required to make sure that the neighboring
        !> points needed by the normal stencil exist otherwise
        !> they are added to the buffer layer. The new grid points
        !> are then computed
        !
        !> @date
        !> 07_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !
        !>@param selected_grdpts
        !> list of bc_interior_pt grid points that should be turned
        !> into interior_pt and therefore need neighboring grid points
        !
        !>@param match_table
        !> (i_match,j_match) needed to identify correctly the grid
        !> points activated by the detectors even if they were
        !> identified before the reallocation
        !--------------------------------------------------------------
        subroutine update_grdpts(
     $       grdpts_id,
     $       nodes,
     $       selected_grdpts,
     $       match_table)

          implicit none

          integer       , dimension(:,:)  , intent(inout) :: grdpts_id
          real(rkind)   , dimension(:,:,:), intent(inout) :: nodes
          integer(ikind), dimension(:,:)  , intent(in)    :: selected_grdpts
          integer(ikind), dimension(2)    , intent(in)    :: match_table


          integer(ikind) :: i,j,k
          integer(ikind) :: i_prev, j_prev, i_center, j_center


          !we have a list of gridpoints that should be turned
          !from bc_inetrior_pt to interior_pt. For a point to be
          !considered an interior point, we need to make sure
          !that the grid points it needs to compute its time
          !derivatives are available
          !
          !previously, the nodes table of the buffer layer was
          !increased to take into account the space needed for
          !new gridpoints at the boundary
          !
          !in this function, we go through the list of gridpoint
          !whose neighbours need to be tested. As it is possible
          !to reduce the number of tests by considering the
          !previous gridpoint tested, there is a distinction
          !between k=1 and k>1 in the list of gridpoints
          !----------------------------------------------------

          !for the first grid point, there is no simplification
          !possible in checking the neighbours so it is separated
          !from the next loop
          !----------------------------------------------------
          k = 1
          i_center = -match_table(1)+selected_grdpts(1,k)
          j_center = -match_table(2)+selected_grdpts(2,k)
          do j=j_center-bc_size, j_center+bc_size
             do i=i_center-bc_size, i_center+bc_size
                call check_gridpoint(grdpts_id,nodes,i,j,i_center,j_center)
             end do
          end do
          grdpts_id(i_center,j_center) = interior_pt


          !from the second gridpoint, we reduce the number of
          !neighbours to be tested
          !----------------------------------------------------
          do k=2, size(selected_grdpts,2)

             !update the position of the gridpoint previously
             !tested
             i_prev   =   i_center
             j_prev   =   j_center
             i_center = - match_table(1) + selected_grdpts(1,k)
             j_center = - match_table(2) + selected_grdpts(2,k)

             !check its neighbours
             call check_neighbors(grdpts_id,nodes,i_prev,j_prev,i_center,j_center)

             !update the status of the gridpoint
             grdpts_id(i_center,j_center) = interior_pt

          end do             

        end subroutine update_grdpts


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> subroutine checking whether the neighboring points
        !> from a bc_interior_pt that should be turned into
        !> an interior point exists: if they exist, nothing
        !> happens. Otherwise, the grid point ID is updated
        !> and the grid point is computed
        !> the number of neighbors to be checked can be reduced
        !> knowing the previous neighbors checked
        !>
        !> @date
        !> 07_04_2013 - initial version - J.L. Desmarais
        !>
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !>
        !>@param i_prev
        !> x-index of the previous bc_interior_pt whose 
        !> neighbors were checked
        !>
        !>@param j_prev
        !> y-index of the previous bc_interior_pt whose 
        !> neighbors were checked
        !>
        !>@param i_center
        !> x-index of the current bc_interior_pt whose 
        !> neighbors were checked
        !>
        !>@param j_center
        !> y-index of the current bc_interior_pt whose 
        !> neighbors were checked
        !---------------------------------------------------------------
        subroutine check_neighbors(grdpts_id,nodes,i_prev,j_prev,i_center,j_center)

          implicit none

          integer       , dimension(:,:)  , intent(inout) :: grdpts_id
          real(rkind)   , dimension(:,:,:), intent(inout) :: nodes
          integer(ikind)                  , intent(in)    :: i_prev
          integer(ikind)                  , intent(in)    :: j_prev
          integer(ikind)                  , intent(in)    :: i_center
          integer(ikind)                  , intent(in)    :: j_center
          

          integer(ikind) :: min_j, max_j
          integer(ikind) :: i,j


          min_j = min(j_center-j_prev,0)
          max_j = max(j_center-j_prev,0)

          !check neighbors for the extra boundary points bc_pt
          do j=j_center-bc_size, j_prev - bc_size + min(j_center-j_prev+2*bc_size,-1)
             do i=i_center-bc_size,i_center+bc_size
                call check_gridpoint(grdpts_id,nodes,i,j,i_center,j_center)
             end do
          end do

          do j=j_center-bc_size-min_j, j_center+bc_size-max_j
             do i=i_center-bc_size, i_prev-bc_size+min(i_center-i_prev+2*bc_size,-1)
                call check_gridpoint(grdpts_id,nodes,i,j,i_center,j_center)
             end do
          end do

          do j=j_center-bc_size-min_j, j_center+bc_size-max_j
             do i=i_prev+bc_size+max(i_center-i_prev-2*bc_size,1),i_center+bc_size
                call check_gridpoint(grdpts_id,nodes,i,j,i_center,j_center)
             end do
          end do

          do j=j_prev+bc_size+max(j_center-j_prev-2*bc_size,1), j_center+bc_size
             do i=i_center-bc_size,i_center+bc_size
                call check_gridpoint(grdpts_id,nodes,i,j,i_center,j_center)
             end do
          end do


          !update the bc_interior_pt overlapping the previous bc_pt
          do j=min(j_center-j_prev+2*bc_size+1,-1), min(j_prev-bc_size,j_center+1) - j_center
             do i=-1,1
                if(grdpts_id(i_center+i,j_center+j).eq.bc_pt) then
                   grdpts_id(i_center+i,j_center+j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(-bc_size-min_j+1,-1), min(bc_size-max_j-1,1)
             do i=min(i_center-i_prev+2*bc_size+1,-1), min(i_prev-bc_size,i_center+1) - i_center
                if(grdpts_id(i_center+i,j_center+j).eq.bc_pt) then
                   grdpts_id(i_center+i,j_center+j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(-bc_size-min_j+1,-1), min(bc_size-max_j-1,1)
             do i=max(i_prev+bc_size,i_center-1)-i_center, max(i_center-i_prev-2*bc_size-1,1)
                if(grdpts_id(i_center+i,j_center+j).eq.bc_pt) then
                   grdpts_id(i_center+i,j_center+j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(j_prev+bc_size,j_center-1)-j_center, max(j_center-j_prev-2*bc_size-1,1)
             do i=-1,1
                if(grdpts_id(i_center+i,j_center+j).eq.bc_pt) then
                   grdpts_id(i_center+i,j_center+j) = bc_interior_pt
                end if
             end do
          end do

        end subroutine check_neighbors      


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> subroutine checking whether the grid point asked
        !> exists or not. If it does not exist in the buffer
        !> layer, it has to be computed and updated in the
        !> gridpoint ID map
        !>
        !> @date
        !> 07_04_2013 - initial version - J.L. Desmarais
        !>
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !>
        !>@param i
        !> x-index of the grid points checked
        !>
        !>@param j
        !> y-index of the grid points checked
        !---------------------------------------------------------------
        subroutine check_gridpoint(
     $     grdpts_id,nodes,i,j,i_center,j_center)

          implicit none

          integer       , dimension(:,:)  , intent(inout) :: grdpts_id
          real(rkind)   , dimension(:,:,:), intent(out)   :: nodes
          integer(ikind)                  , intent(in)    :: i
          integer(ikind)                  , intent(in)    :: j
          integer(ikind)                  , intent(in)    :: i_center
          integer(ikind)                  , intent(in)    :: j_center

          if (grdpts_id(i,j).eq.no_pt) then

             call this%compute_new_grdpt(i,j)

             !if the grid point is next to the new interior point,
             !it is a bc_interior_pt, otherwise, it is a bc_pt
             if(abs(i_center-1).le.1.and.abs(j_center-j).le.1) then
                grdpts_id(i,j) = bc_interior_pt
             else
                grdpts_id(i,j) = bc_pt
             end if

          end if

        end subroutine check_gridpoint

      end module bf_layer_update_grdpts_module

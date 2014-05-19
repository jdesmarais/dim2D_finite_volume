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

        use parameters_bf_layer, only : no_pt, interior_pt,
     $                                  bc_interior_pt, bc_pt
        use parameters_constant, only : N,S,E,W
        use parameters_input   , only : nx,ny,bc_size
        use parameters_kind    , only : ikind

        private
        public :: ini_grdpts_id_N,
     $            ini_grdpts_id_S,
     $            ini_grdpts_id_E,
     $            ini_grdpts_id_W

        contains

        
        !< initialize the grdpts_id for the northern buffer layer
        subroutine ini_grdpts_id_N(grdpts_id, alignment)

          implicit none

          integer       , dimension(:,:), intent(out) :: grdpts_id
          integer(ikind), dimension(2,2), intent(in)  :: alignment

          call add_interior_layer_NS(
     $         grdpts_id, 1, bc_size, alignment)

          call add_bc_interior_layer_NS(
     $         grdpts_id, bc_size+1, alignment)

          call add_bc_layer_NS(
     $         grdpts_id, bc_size+2)

          call add_no_pt_layer_NS(
     $         grdpts_id, bc_size+3, size(grdpts_id,2))

        end subroutine ini_grdpts_id_N


        !< initialize the grdpts_id for the southern buffer layer
        subroutine ini_grdpts_id_S(grdpts_id, alignment)

          implicit none

          integer       , dimension(:,:), intent(out) :: grdpts_id
          integer(ikind), dimension(2,2), intent(in)  :: alignment

          call add_no_pt_layer_NS(
     $         grdpts_id, 1, size(grdpts_id,2)-(bc_size+2))

          call add_bc_layer_NS(
     $         grdpts_id, size(grdpts_id,2)-(bc_size+2)+1)

          call add_bc_interior_layer_NS(
     $         grdpts_id, size(grdpts_id,2)-(bc_size+2)+2, alignment)

          call add_interior_layer_NS(
     $         grdpts_id,
     $         size(grdpts_id,2)-(bc_size+2)+3,
     $         size(grdpts_id,2), alignment)

        end subroutine ini_grdpts_id_S

        !< initialize the grdpts_id for the eastern buffer layer
        subroutine ini_grdpts_id_E(grdpts_id, alignment)

          implicit none

          integer       , dimension(:,:), intent(out) :: grdpts_id
          integer(ikind), dimension(2,2), intent(in)  :: alignment

          call add_edge_layer_E(
     $         grdpts_id,
     $         1, bc_size,
     $         alignment, .false.)

          call add_interior_layer_E(
     $         grdpts_id, bc_size+1, size(grdpts_id,2)-bc_size)

          call add_edge_layer_E(
     $         grdpts_id,
     $         size(grdpts_id,2)-bc_size+1, size(grdpts_id,2),
     $         alignment, .true.)
          
        end subroutine ini_grdpts_id_E


        !< initialize the grdpts_id for the eastern buffer layer
        subroutine ini_grdpts_id_W(grdpts_id, alignment)

          implicit none

          integer       , dimension(:,:), intent(out) :: grdpts_id
          integer(ikind), dimension(2,2), intent(in)  :: alignment

          call add_edge_layer_W(
     $         grdpts_id,
     $         1, bc_size,
     $         alignment, .false.)

          call add_interior_layer_W(
     $         grdpts_id, bc_size+1, size(grdpts_id,2)-bc_size)

          call add_edge_layer_W(
     $         grdpts_id,
     $         size(grdpts_id,2)-bc_size+1, size(grdpts_id,2),
     $         alignment, .true.)
          
        end subroutine ini_grdpts_id_W


        subroutine add_interior_layer_NS(
     $     grdpts_id, j_min, j_max, alignment)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j_min
          integer(ikind)                , intent(in)    :: j_max
          integer(ikind), dimension(2,2), intent(in)    :: alignment


          integer(ikind) :: i,j
          integer(ikind) :: i_min1,i_min2
          integer(ikind) :: i_max1, i_max2, i_max3

          
          if(alignment(1,1).gt.(bc_size+2)) then
             i_min1=0
             i_min2=i_min1+0
          else
             if(alignment(1,1).eq.(bc_size+2)) then
                i_min1=0
                i_min2=1
             else
                i_min1=bc_size-1
                i_min2=i_min1+1
             end if
          end if

          if(alignment(1,2).lt.(nx-bc_size-1)) then
             i_max1 = size(grdpts_id,1)
             i_max2 = 0
             i_max3 = 0
          else
             if(alignment(1,2).eq.(nx-bc_size-1)) then
                i_max1 = size(grdpts_id,1)-1
                i_max2 = size(grdpts_id,1)
                i_max3 = 0
             else
                i_max1 = size(grdpts_id,1)-bc_size
                i_max2 = size(grdpts_id,1)-bc_size+1
                i_max3 = size(grdpts_id,1)
             end if
          end if


          do j=j_min,j_max
             do i=1,i_min1
                grdpts_id(i,j) = bc_pt
             end do

             do i=i_min1+1,i_min2
                grdpts_id(i,j) = bc_interior_pt
             end do

             do i=i_min2+1, i_max1
                grdpts_id(i,j) = interior_pt
             end do

             do i=i_max1+1, i_max2
                grdpts_id(i,j) = bc_interior_pt
             end do

             do i=i_max2+1, i_max3
                grdpts_id(i,j) = bc_pt
             end do
          end do

        end subroutine add_interior_layer_NS


        subroutine add_bc_interior_layer_NS(
     $     grdpts_id, j, alignment)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j
          integer(ikind), dimension(2,2), intent(in)    :: alignment


          integer(ikind) :: i
          integer(ikind) :: i_min1
          integer(ikind) :: i_max1, i_max2

          
          if(alignment(1,1).gt.(bc_size+1)) then
             i_min1=0
          else
             if(alignment(1,1).eq.(bc_size+1)) then
                i_min1=1
             end if
          end if

          if(alignment(1,2).lt.(nx-bc_size)) then
             i_max1 = size(grdpts_id,1)
             i_max2 = 0
          else
             if(alignment(1,2).eq.(nx-bc_size)) then
                i_max1 = size(grdpts_id,1)-1
                i_max2 = size(grdpts_id,1)
             end if
          end if

          do i=1,i_min1
             grdpts_id(i,j) = bc_pt
          end do

          do i=i_min1+1,i_max1
             grdpts_id(i,j) = bc_interior_pt
          end do

          do i=i_max1+1, i_max2
             grdpts_id(i,j) = bc_pt
          end do
        end subroutine add_bc_interior_layer_NS


        subroutine add_bc_layer_NS(
     $     grdpts_id, j)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j

          integer(ikind) :: i

          do i=1,size(grdpts_id,1)
             grdpts_id(i,j) = bc_pt
          end do

        end subroutine add_bc_layer_NS


        subroutine add_no_pt_layer_NS(
     $     grdpts_id, j_min, j_max)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j_min
          integer(ikind)                , intent(in)    :: j_max

          integer(ikind) :: i,j

          do j=j_min, j_max
             do i=1,size(grdpts_id,1)
                grdpts_id(i,j) = no_pt
             end do
          end do

        end subroutine add_no_pt_layer_NS


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
                grdpts_id(i,j) = interior_pt
             end do

             i=bc_size+1
             grdpts_id(i,j) = bc_interior_pt

             do i=bc_size+2, 2*bc_size
                grdpts_id(i,j) = bc_pt
             end do

             do i=2*bc_size+1, size(grdpts_id,1)
                grdpts_id(i,j) = no_pt
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
             do i=1, size(grdpts_id,1)-(2*bc_size)
                grdpts_id(i,j) = no_pt
             end do

             do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                grdpts_id(i,j) = bc_pt
             end do

             i=size(grdpts_id,1)-bc_size
             grdpts_id(i,j) = bc_interior_pt
             
             do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                grdpts_id(i,j) = interior_pt
             end do
             
          end do          

        end subroutine add_interior_layer_W


        !< add the bc_pt layer for the initialization
        !> of gridpts_id
        !> type=.true.  = upper layer
        !> type=.false. = lower layer
        subroutine add_edge_layer_E(
     $     grdpts_id, j_min, j_max,
     $     alignment, type)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j_min
          integer(ikind)                , intent(in)    :: j_max
          integer(ikind), dimension(2,2), intent(in)    :: alignment
          logical                       , intent(in)    :: type

          integer(ikind) :: i,j

          if(type) then
             if(alignment(2,2).lt.(ny-bc_size-1)) then
                do j=j_min, j_max
                   do i=1, bc_size
                      grdpts_id(i,j) = interior_pt
                   end do
                   i=bc_size+1
                   grdpts_id(i,j) = bc_interior_pt

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                end do
             else

                if(alignment(2,2).eq.(ny-bc_size-1)) then
                   j=j_min
                   do i=1, bc_size
                      grdpts_id(i,j) = interior_pt
                   end do
                   i=bc_size+1
                   grdpts_id(i,j) = bc_interior_pt

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   j=j_min+1
                   do i=1,bc_size+1
                      grdpts_id(i,j)=bc_interior_pt
                   end do

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                else

                   j=j_min
                   do i=1, bc_size+1
                      grdpts_id(i,j) = bc_interior_pt
                   end do

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do

                   j=j_min+1
                   do i=1,2*bc_size
                      grdpts_id(i,j)=bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do                   
                   
                end if
             end if
          else
             if(alignment(2,1).gt.(bc_size+2)) then
                do j=j_min, j_max
                   do i=1, bc_size
                      grdpts_id(i,j) = interior_pt
                   end do
                   i=bc_size+1
                   grdpts_id(i,j) = bc_interior_pt

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                end do
             else

                if(alignment(2,1).eq.(bc_size+2)) then
                   j=j_min
                   do i=1,bc_size+1
                      grdpts_id(i,j)=bc_interior_pt
                   end do

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   j=j_min+1
                    do i=1, bc_size
                      grdpts_id(i,j) = interior_pt
                   end do
                   i=bc_size+1
                   grdpts_id(i,j) = bc_interior_pt

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                else

                   j=j_min
                   do i=1,2*bc_size
                      grdpts_id(i,j)=bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do

                   j=j_min+1
                   do i=1, bc_size+1
                      grdpts_id(i,j) = bc_interior_pt
                   end do

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do                   
                   
                end if
             end if
          end if

        end subroutine add_edge_layer_E


        !< add the bc_pt layer for the initialization
        !> of gridpts_id
        !> type=.true.  = upper layer
        !> type=.false. = lower layer
        subroutine add_edge_layer_W(
     $     grdpts_id, j_min, j_max, alignment, type)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j_min
          integer(ikind)                , intent(in)    :: j_max
          integer(ikind), dimension(2,2), intent(in)    :: alignment
          logical                       , intent(in)    :: type

          integer(ikind) :: i,j

          if(type) then

             if(alignment(2,2).lt.(ny-bc_size-1)) then

                do j=j_min, j_max
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   i=size(grdpts_id,1)-bc_size
                   grdpts_id(i,j) = bc_interior_pt
                   
                   do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = interior_pt
                   end do
                end do

             else

                if(alignment(2,2).eq.(ny-bc_size-1)) then
                   j=j_min
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   i=size(grdpts_id,1)-bc_size
                   grdpts_id(i,j) = bc_interior_pt
                   
                   do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = interior_pt
                   end do

                   j=j_min+1
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=size(grdpts_id,1)-bc_size, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_interior_pt
                   end do
                   
                   
                else

                   j=j_min
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=size(grdpts_id,1)-bc_size, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_interior_pt
                   end do

                   j=j_min+1
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_pt
                   end do

                end if
             end if

          else

             if(alignment(2,1).gt.(bc_size+2)) then
                do j=j_min, j_max
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   i=size(grdpts_id,1)-bc_size
                   grdpts_id(i,j) = bc_interior_pt
                   
                   do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = interior_pt
                   end do
                end do
             else

                if(alignment(2,1).eq.(bc_size+2)) then
                   j=j_min
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=size(grdpts_id,1)-bc_size, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_interior_pt
                   end do

                   j=j_min+1
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   i=size(grdpts_id,1)-bc_size
                   grdpts_id(i,j) = bc_interior_pt
                   
                   do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = interior_pt
                   end do
                   
                else

                   j=j_min
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_pt
                   end do                   

                   j=j_min+1
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=size(grdpts_id,1)-bc_size, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_interior_pt
                   end do
                end if
             end if
          end if

        end subroutine add_edge_layer_W

      end module bf_layer_ini_grdptID_module

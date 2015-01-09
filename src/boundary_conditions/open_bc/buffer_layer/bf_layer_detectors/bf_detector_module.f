      !> @file
      !> module gathering the subroutines needed when computing the
      !> number and position of detectors to be inserted between two
      !> detectors to have a continuous path of increasing detectors
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module gathering the subroutines needed when computing the
      !> number and position of detectors to be inserted between two
      !> detectors to have a continuous path of increasing detectors
      !
      !> @date
      ! 24_11_2014 - documentation update - J.L. Desmarais
      ! 07_01_2015 - re-implementing the rcoords determination
      !              of the intermediate detectors to be rotation
      !              invariant - J.L.Desmarais
      !----------------------------------------------------------------
      module bf_detector_module

        use parameters_input, only :
     $       nx,ny

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: 
     $       determine_local_map_coordinates,
     $       get_inter_detector_param,
     $       get_inter_detector_coords,
     $       
     $       get_rot_coords,
     $       get_minimum_block_nb,
     $       get_intermediate_maps,
     $       get_inter_detector_rot_param,
     $       get_inter_detector_rot_coords

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the cartesian coordinates corresponding
        !> to the index-coordinates of the intermediate
        !> detectors
        !
        !> @date
        !> 09_01_2015 - initial version - J.L. Desmarais
        !
        !>@param interior_map
        !> coordinates of the interior point either along the direction
        !> investigated
        !
        !>@param size_interior_map
        !> size of the coordinate maps along the direction investigated
        !> (either nx or ny)
        !
        !>@param size_local_map
        !> size of the local map of coordinates manufactured from the
        !> interior map
        !
        !>@param first_icoord
        !> index coordinate identifying the first point of the local map
        !
        !>@param icr_sign
        !> boolean determining whether the local map should be made in
        !> increasing or decreasing direction
        !
        !>@param local_map
        !> local map of coordinates manufactured from the interior map
        !--------------------------------------------------------------
        subroutine determine_local_map_coordinates(
     $       interior_map,
     $       size_interior_map,
     $       size_local_map,
     $       first_icoord,
     $       icr_sign,
     $       local_map)

          implicit none

          real(rkind)   , dimension(:), intent(in)  :: interior_map
          integer                     , intent(in)  :: size_interior_map
          integer                     , intent(in)  :: size_local_map
          integer(ikind)              , intent(in)  :: first_icoord
          logical                     , intent(in)  :: icr_sign
          real(rkind)   , dimension(:), intent(out) :: local_map


          integer        :: size_outside
          integer        :: size_inside
          integer(ikind) :: first_icoord_outside
          integer(ikind) :: first_icoord_inside
          integer(ikind) :: last_icoord_inside
          integer        :: local_icoord
          real(rkind)    :: ds
          integer        :: k_start
          integer        :: k


          if(icr_sign) then

             !fill the local map with coordinates on
             !the left side outside the interior_map
             if((first_icoord+size_local_map-1).le.0) then
                size_outside = size_local_map
             else
                size_outside = (0-first_icoord)+1
             end if
             !size_outside = 0-min(1,first_icoord)+1

             if(size_outside.gt.0) then
                
                ds = interior_map(2)-interior_map(1)
                
                do k=1, size_outside
                   local_icoord = first_icoord+(k-1)
                   local_map(k) = interior_map(1) + (local_icoord-1)*ds
                end do
                
                k_start = size_outside

             else

                k_start = 0

             end if


             !fill the local map with coordinates from
             !the interior_map
             first_icoord_inside = max(first_icoord,1)
             last_icoord_inside  = min(first_icoord+size_local_map-1,size_interior_map)
             size_inside         = last_icoord_inside - first_icoord_inside+1

             if(size_inside.gt.0) then

                do k=1, size_inside
                   local_icoord         = first_icoord_inside+(k-1)
                   local_map(k_start+k) = interior_map(local_icoord)
                end do
                
                k_start = k_start+size_inside

             end if


             !fill the local map with coordinates on the right
             !side outside the interior_map
             if(first_icoord.ge.(size_interior_map+1)) then
                size_outside = size_local_map
             else
                size_outside = (first_icoord+size_local_map-1)-(size_interior_map+1)+1
             end if

             if(size_outside.gt.0) then

                ds = interior_map(size_interior_map)-
     $               interior_map(size_interior_map-1)

                if(k_start.ne.0) then
                   first_icoord_outside = size_interior_map+1
                else
                   first_icoord_outside = first_icoord
                end if

                do k=1, size_outside
                   local_icoord         = first_icoord_outside+(k-1)
                   local_map(k_start+k) = interior_map(size_interior_map) +
     $                  (local_icoord-size_interior_map)*ds
                end do

             end if

          else

             !fill the local map with coordinates on the right
             !side outside the interior_map
             if((first_icoord-size_local_map+1).ge.(size_interior_map+1)) then
                size_outside = size_local_map
             else
                size_outside = first_icoord-(size_interior_map+1)+1
             end if

             if(size_outside.gt.0) then

                ds = interior_map(size_interior_map)-
     $               interior_map(size_interior_map-1)

                first_icoord_outside = first_icoord

                do k=1, size_outside
                   local_icoord = first_icoord_outside-(k-1)
                   local_map(k) = interior_map(size_interior_map) +
     $                  (local_icoord-size_interior_map)*ds
                end do

                k_start = size_outside

             else

                k_start = 0
                
             end if


             !fill the local map with coordinates from
             !the interior_map
             first_icoord_inside = max(first_icoord-size_local_map+1,1)
             last_icoord_inside  = min(first_icoord,size_interior_map)
             size_inside         = last_icoord_inside - first_icoord_inside+1

             if(size_inside.gt.0) then

                if(k_start.ne.0) then
                   first_icoord_inside = size_interior_map
                else
                   first_icoord_inside = first_icoord
                end if                   

                do k=1, size_inside
                   local_icoord         = first_icoord_inside-(k-1)
                   local_map(k_start+k) = interior_map(local_icoord)
                end do
                
                k_start = k_start+size_inside

             end if
             

             !fill the local map with coordinates on
             !the left side outside the interior_map
             if(first_icoord.le.0) then
                size_outside = size_local_map
             else
                size_outside = 0-(first_icoord-size_local_map+1) + 1
             end if

             !size_outside = -min(1,first_icoord-size_local_map+1)+1

             if(size_outside.gt.0) then
                
                ds = interior_map(2)-interior_map(1)
                
                if(k_start.ne.0) then
                   first_icoord_outside = 0
                else
                   first_icoord_outside = first_icoord
                end if

                do k=1, size_outside
                   local_icoord = first_icoord_outside-(k-1)
                   local_map(k_start+k) = interior_map(1) + (local_icoord-1)*ds
                end do

             end if

          end if

        end subroutine determine_local_map_coordinates



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the (x,y)-index coordinates and the (x,y)-coordinates of
        !> the intermediate detector corresponding to position k on the
        !> straight line linking the two blocks (\f$k \in
        !> [0,\text{inter_nb}+1]\f$).
        !> For \f$k=0\f$, the coordinates match the first block.
        !> For \f$k=inter_nb+1\f$, the coordinates match the last block
        !
        !> @date
        !> 09_01_2015 - initial version - J.L. Desmarais
        !
        !>@param prev_icoords
        !> (x,y)-index coordinates of the first block to be linked
        !
        !>@param next_icoords
        !> (x,y)-index coordinates of the last block to be linked
        !
        !>@param interior_x_map
        !> map of the interior x-coordinates
        !
        !>@param interior_y_map
        !> map of the interior y-coordinates
        !
        !>@param icoords_icr
        !> slope of the (x,y)-index increase b/w two blocks when determining
        !> the (x,y)-index coordinates of the intermediate detectors
        !
        !>@param rot_icoords_r
        !> (x,y)-index coordinates of the rotation point (as it can be located
        !> at the edge b/w two blocks, the coordinates are expressed with double
        !> to allow for the m+0.5 case, m is an integer)
        !
        !>@param inter_nb
        !> number of intermediate detectors b/w the two blocks
        !
        !>@param x_map_icr
        !> local map of the x-coordinates of the intermediate points b/w
        !> the two detectors
        !
        !>@param y_map_icr
        !> local map of the y-coordinates of the intermediate points b/w
        !> the two detectors
        !
        !>@param rot_icoords_r
        !> (x,y)-index coordinates of the middle point b/w the two blocks
        !> to be linked. In order for the path of intermediate detectors
        !> b/w the two blocks to be rotation invariant, this point is used
        !> as the rotation point
        !
        !>@param icoords_icr
        !> slope of the straight line linking the two blocks expressed as
        !> (x,y)-coordinates
        !
        !>@param k
        !> index identifying the position of the intermediate detector on
        !> the straight line linking the two blocks
        !
        !>@param icoords_inter
        !> (x,y)-index coordinates of the intermediate detector
        !
        !>@param rcoords_inter
        !> (x,y)-coordinates of the intermediate detector
        !--------------------------------------------------------------
        subroutine get_inter_detector_rot_coords(
     $     prev_icoord,
     $     rot_icoords_r,
     $     icoords_icr,
     $     inter_nb,
     $     x_map_icr,
     $     y_map_icr,
     $     k,
     $     icoord_inter,
     $     rcoord_inter)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: prev_icoord
          real(rkind)   , dimension(2), intent(in)  :: rot_icoords_r
          real(rkind)   , dimension(2), intent(in)  :: icoords_icr
          integer(ikind)              , intent(in)  :: inter_nb
          real(rkind)   , dimension(:), intent(in)  :: x_map_icr
          real(rkind)   , dimension(:), intent(in)  :: y_map_icr
          integer                     , intent(in)  :: k
          integer(ikind), dimension(2), intent(out) :: icoord_inter
          real(rkind)   , dimension(2), intent(out) :: rcoord_inter


          icoord_inter = get_inter_detector_icoords(
     $         rot_icoords_r,
     $         icoords_icr,
     $         inter_nb,
     $         k)
          
          rcoord_inter(1) = x_map_icr(abs(icoord_inter(1)-prev_icoord(1))+1)
          rcoord_inter(2) = y_map_icr(abs(icoord_inter(2)-prev_icoord(2))+1)

        end subroutine get_inter_detector_rot_coords


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the (x,y)-index coordinates of the intermediate detector
        !> corresponding to position k on the straight line linking the
        !> two blocks (k \f$\in [0,\text{inter_nb}+1]\$).
        !> For \f$k=0\f$, the coordinates match the first block.
        !> For \f$k=inter_nb+1\f$, the coordinates match the last block
        !
        !> @date
        !> 09_01_2015 - initial version - J.L. Desmarais
        !
        !>@param rot_icoords_r
        !> (x,y)-index coordinates of the middle point b/w the two blocks
        !> to be linked. In order for the path of intermediate detectors
        !> b/w the two blocks to be rotation invariant, this point is used
        !> as the rotation point
        !
        !>@param icoords_icr
        !> slope of the straight line linking the two blocks expressed as
        !> (x,y)-coordinates
        !
        !>@param inter_nb
        !> number of intermediate detectors to link the two blocks
        !
        !>@param k
        !> index identifying the position of the intermediate detector on
        !> the straight line linking the two blocks
        !
        !>@return icoords_inter
        !> (x,y)-index coordinates of the intermediate detector
        !--------------------------------------------------------------
         function get_inter_detector_icoords(
     $     rot_icoords_r,
     $     icoords_icr,
     $     inter_nb,
     $     k)
     $     result(icoords_inter)

          implicit none

          real(rkind)   , dimension(2), intent(in) :: rot_icoords_r
          real(rkind)   , dimension(2), intent(in) :: icoords_icr
          integer                     , intent(in) :: inter_nb
          integer                     , intent(in) :: k
          integer(ikind), dimension(2)             :: icoords_inter

          real(rkind) :: k_rot, k_icr
          integer :: i

          k_rot = real(inter_nb+1)/2.0d0
          k_icr = sign(abs(k-k_rot),k-k_rot)

          do i=1,2
             icoords_inter(i) = nint(rot_icoords_r(i)+k_icr*icoords_icr(i))
          end do

        end function get_inter_detector_icoords


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the parameters allowing to determine the coordinates of
        !> the intermediate detectors linking two blocks: the path
        !> should be rotation invariant to preserve symmetry
        !
        !> @date
        !> 09_01_2015 - initial version - J.L. Desmarais
        !
        !>@param prev_icoords
        !> (x,y)-index coordinates of the first block to be linked
        !
        !>@param next_icoords
        !> (x,y)-index coordinates of the last block to be linked
        !
        !>@param interior_x_map
        !> map of the interior x-coordinates
        !
        !>@param interior_y_map
        !> map of the interior y-coordinates
        !
        !>@param rot_icoords_r
        !> (x,y)-index coordinates of the middle point b/w the two blocks
        !> to be linked. In order for the path of intermediate detectors
        !> b/w the two blocks to be rotation invariant, this point is used
        !> as the rotation point
        !
        !>@param icoords_icr
        !> slope of the straight line linking the two blocks expressed as
        !> (x,y)-coordinates
        !
        !>@param inter_nb
        !> number of intermediate detectors to link the two blocks
        !
        !>@param x_map_icr
        !> local map of the x-coordinates of the intermediate points b/w
        !> the two detectors
        !
        !>@param y_map_icr
        !> local map of the y-coordinates of the intermediate points b/w
        !> the two detectors
        !--------------------------------------------------------------
        subroutine get_inter_detector_rot_param(
     $     prev_icoords,
     $     next_icoords,
     $     interior_x_map,
     $     interior_y_map,
     $     rot_icoords_r,
     $     icoords_icr,
     $     inter_nb,
     $     x_map_icr,
     $     y_map_icr)

          implicit none

          integer(ikind), dimension(2)             , intent(in)  :: prev_icoords
          integer(ikind), dimension(2)             , intent(in)  :: next_icoords
          real(rkind)   , dimension(nx)            , intent(in)  :: interior_x_map
          real(rkind)   , dimension(ny)            , intent(in)  :: interior_y_map
          real(rkind)   , dimension(2)             , intent(out) :: rot_icoords_r
          real(rkind)   , dimension(2)             , intent(out) :: icoords_icr
          integer(ikind)                           , intent(out) :: inter_nb
          real(rkind)   , dimension(:), allocatable, intent(out) :: x_map_icr
          real(rkind)   , dimension(:), allocatable, intent(out) :: y_map_icr

          integer(ikind), dimension(2) :: rot_icoords
          integer(ikind)               :: rot_block_nb
          integer(ikind)               :: total_nb_blocks

          integer     :: i
          integer     :: diff
          real(rkind) :: icr

          !get the (x,y)-index coordinates of the rotation block
          call get_rot_coords(
     $         prev_icoords,
     $         next_icoords,
     $         rot_icoords,
     $         rot_icoords_r)
          
          !get the minimum number of blocks needed to join the first
          !block with prev_icoords coordinates with the last block with
          !rot_icoords coordinates
          rot_block_nb = get_minimum_block_nb(prev_icoords,rot_icoords)

          !if the rotation point is at the intersection b/w blocks
          !  _ _ _ _ 
          ! |_|_|_|2|
          ! |_|0+0|_|   : total_nb_blocks = 2*rot_block_nb
          ! |1|_|_|_|
          !  _ _ _ _ 
          ! |_|_|_|2|
          ! |_|_|0|_|   : total_nb_blocks = 2*rot_block_nb
          ! |_|0|_|_|   
          ! |1|_|_|_|
          !
          !if the rotation point is at the nodal point of a block
          !
          !  _ _ _ _ _ 
          ! |_|_|_|_|2|
          ! |_|0|+|0|_| : total_nb_blocks = 2*rot_block_nb+1
          ! |1|_|_|_|_|
          !
          !--------------------------------------------------------
          if((mod(prev_icoords(1)+next_icoords(1),2).eq.1).or.
     $       (mod(prev_icoords(2)+next_icoords(2),2).eq.1)) then
             total_nb_blocks = 2*(rot_block_nb+1)
          else
             total_nb_blocks = 2*rot_block_nb+1
          end if

          !we now need to verify that this number of blocks satisfy
          !the following condition:
          !when creating the path of intermediate detectors and
          !computing the index-coordinates of these detectors, we
          !should have a continuous path of intermediate detectors
          !i.e. the increase of index-coordinate b/w two detectors
          !should be lower or equal to one
          inter_nb = max(0,total_nb_blocks)

          if(inter_nb.gt.0) then

             do i=1,2
                
                diff = abs(next_icoords(i)-prev_icoords(i))
                icr  = diff/(inter_nb+1)

                do while(icr.gt.(1.0d0))
                   
                   inter_nb = inter_nb + 2
                   icr = diff/(inter_nb+1)
                   
                end do
                
             end do
             
             icoords_icr(1) = real(next_icoords(1)-prev_icoords(1))/real(inter_nb+1)
             icoords_icr(2) = real(next_icoords(2)-prev_icoords(2))/real(inter_nb+1)
             
             !determine the local x_map and y_map corresponding 
             !to the intermediate detectors
             call get_intermediate_maps(
     $            prev_icoords,
     $            interior_x_map,
     $            interior_y_map,
     $            icoords_icr,
     $            rot_icoords_r,
     $            inter_nb,
     $            x_map_icr,
     $            y_map_icr)
             
          end if

        end subroutine get_inter_detector_rot_param


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the local maps needed to compute the (x,y)-coordinates
        !> of the intermediate detectors b/w two blocks
        !
        !> @date
        !> 09_01_2015 - initial version - J.L. Desmarais
        !
        !>@param prev_icoords
        !> (x,y)-index coordinates of the first block to be linked
        !
        !>@param next_icoords
        !> (x,y)-index coordinates of the last block to be linked
        !
        !>@param interior_x_map
        !> map of the interior x-coordinates
        !
        !>@param interior_y_map
        !> map of the interior y-coordinates
        !
        !>@param icoords_icr
        !> slope of the (x,y)-index increase b/w two blocks when determining
        !> the (x,y)-index coordinates of the intermediate detectors
        !
        !>@param rot_icoords_r
        !> (x,y)-index coordinates of the rotation point (as it can be located
        !> at the edge b/w two blocks, the coordinates are expressed with double
        !> to allow for the m+0.5 case, m is an integer)
        !
        !>@param inter_nb
        !> number of intermediate detectors b/w the two blocks
        !
        !>@param x_map_icr
        !> local map of the x-coordinates of the intermediate points b/w
        !> the two detectors
        !
        !>@param y_map_icr
        !> local map of the y-coordinates of the intermediate points b/w
        !> the two detectors
        !--------------------------------------------------------------
        subroutine get_intermediate_maps(
     $     prev_icoords,
     $     interior_x_map,
     $     interior_y_map,
     $     icoords_icr,
     $     rot_icoords_r,
     $     inter_nb,
     $     x_map_icr,
     $     y_map_icr)

          implicit none

          integer(ikind), dimension(2)              , intent(in)  :: prev_icoords
          real(rkind)   , dimension(nx)             , intent(in)  :: interior_x_map
          real(rkind)   , dimension(ny)             , intent(in)  :: interior_y_map
          real(rkind)   , dimension(2)              , intent(in)  :: icoords_icr
          real(rkind)   , dimension(2)              , intent(in)  :: rot_icoords_r
          integer                                   , intent(in)  :: inter_nb
          real(rkind)   , dimension(:) , allocatable, intent(out) :: x_map_icr
          real(rkind)   , dimension(:) , allocatable, intent(out) :: y_map_icr

          integer(ikind), dimension(2) :: first_icoord
          integer(ikind), dimension(2) :: last_icoord
          integer(ikind), dimension(2) :: size_local_map
          integer                      :: i

          !allocation of the tables storing the (x,y)-
          !coordinates corresponding of the index
          first_icoord = get_inter_detector_icoords(
     $         rot_icoords_r,
     $         icoords_icr,
     $         inter_nb,
     $         1)
          
          last_icoord = get_inter_detector_icoords(
     $         rot_icoords_r,
     $         icoords_icr,
     $         inter_nb,
     $         inter_nb)

          do i=1,2

             size_local_map(i) = abs(last_icoord(i)-first_icoord(i))+1
             
             if(prev_icoords(i).ne.first_icoord(i)) then
                size_local_map(i) = size_local_map(i)+1
             end if
             
          end do
          
          allocate(x_map_icr(size_local_map(1)))
          allocate(y_map_icr(size_local_map(2)))

          
          !fill the x_map_icr and the y_map_icr with the
          !coorresponding (x,y)-coordinates
          call determine_local_map_coordinates(
     $         interior_x_map,nx,
     $         size_local_map(1),
     $         prev_icoords(1),
     $         (icoords_icr(1).gt.0),
     $         x_map_icr)
          
          call determine_local_map_coordinates(
     $         interior_y_map,ny,
     $         size_local_map(2),
     $         prev_icoords(2),
     $         (icoords_icr(2).gt.0),
     $         y_map_icr)

        end subroutine get_intermediate_maps

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the minimum number of blocks to link two blocks
        !
        !> @date
        !> 09_01_2015 - initial version - J.L. Desmarais
        !
        !>@param first_icoords
        !> (x,y)-index coordinates of the first block to be linked
        !
        !>@param last_icoords
        !> (x,y)-index coordinates of the last block to be linked
        !
        !>@return block_nb
        !> minimal number of block to link the first and the last blocks
        !--------------------------------------------------------------
        function get_minimum_block_nb(first_icoords,last_icoords)
     $     result(block_nb)

          implicit none

          integer(ikind), dimension(2), intent(in) :: first_icoords
          integer(ikind), dimension(2), intent(in) :: last_icoords
          integer                                  :: block_nb

          !the minimum number of blocks between two blocks whose
          !(x,y)-index coordinates are first_icoords(:) for the
          !first block and last_icoords(:) for the last block
          !is the sum of the maximum number of blocks located
          !along a diagonal and the minium number of block
          !located along the horizontal or vertical direction:
          !
          !            left_block_nb
          !                   |
          !               |-------|
          !    _ _ _ _ _ _ _ _ _ _ _
          !   |_|_|_|_|_|_|_|_|_|_|_|
          !   |_|_|_|_|_|_|0|0|0|0|2|  _
          !   |_|_|_|_|_|0|_|_|_|_|_|  |
          !   |_|_|_|_|0|_|_|_|_|_|_|  | 
          !   |_|_|_|0|_|_|_|_|_|_|_|  | diagonal_block_nb
          !   |_|_|0|_|_|_|_|_|_|_|_|  |
          !   |_|0 _|_|_|_|_|_|_|_|_|  |
          !   |1|_|_|_|_|_|_|_|_|_|_| 
          !
          !-----------------------------------------------------
          integer(ikind) :: diagonal_block_nb
          integer(ikind) :: horizontal_block_nb
          integer(ikind) :: vertical_block_nb
          integer(ikind) :: left_block_nb

          !computation of the number of diagonal blocks
          diagonal_block_nb = min(
     $         abs(last_icoords(1)-first_icoords(1))-1,
     $         abs(last_icoords(2)-first_icoords(2))-1)

          !computation of the number of left blocks
          horizontal_block_nb = abs(
     $         last_icoords(1)-(first_icoords(1)+
     $         sign(1,last_icoords(1)-first_icoords(1))*
     $         diagonal_block_nb))-1

          vertical_block_nb = abs(
     $         last_icoords(2)-(first_icoords(2)+
     $         sign(1,last_icoords(2)-first_icoords(2))*
     $         diagonal_block_nb))-1

          left_block_nb = max(horizontal_block_nb,vertical_block_nb)

          !minimum number of blocks
          block_nb = diagonal_block_nb + left_block_nb

        end function get_minimum_block_nb


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the position of the middle point linking the two blocks
        !> the position is given as both the (x,y)-index coordinates in
        !> double to express the position of the point even if at the
        !> edge b/w two blocks and also as the (x,y)-index to identify
        !> the block closest to the middle point in the direction of
        !> the first block
        !
        !> @date
        !> 09_01_2015 - initial version - J.L. Desmarais
        !
        !>@param prev_icoords
        !> (x,y)-index coordinates of the first block to be linked
        !
        !>@param next_icoords
        !> (x,y)-index coordinates of the last block to be linked
        !
        !>@param rot_icoords
        !> (x,y)-coordinates of the middle point
        !
        !>@param rot_rcoords
        !> (x,y)-coordinates of the middle point as double
        !--------------------------------------------------------------
        subroutine get_rot_coords(
     $     prev_icoords,
     $     next_icoords,
     $     rot_icoords,
     $     rot_rcoords)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: prev_icoords
          integer(ikind), dimension(2), intent(in)  :: next_icoords
          integer(ikind), dimension(2), intent(out) :: rot_icoords
          real(rkind)   , dimension(2), intent(out) :: rot_rcoords

          integer :: i

          do i=1,2

             rot_rcoords(i) = 0.5d0*(prev_icoords(i)+next_icoords(i))

             if(prev_icoords(i).le.rot_rcoords(i)) then
                rot_icoords(i) = floor(rot_rcoords(i))
             else
                rot_icoords(i) = ceiling(rot_rcoords(i))
             end if

          end do

        end subroutine get_rot_coords


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !< get the parameters constraining the addition
        !> of intermediate detectors between the previous
        !> detectors and the new detector to be added
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param prev_coords
        !> general coordinates of the last detector added to the list
        !
        !>@param next_coords
        !> general coordinates of the future detector to be added to the list
        !
        !>@param x_change
        !> change in the x-coordinate between the previous and the next
        !> detector
        !
        !>@param y_change
        !> change in the y-coordinate between the previous and the next
        !> detector
        !
        !>@param inter_nb
        !> number of detectors to be added between the two to ensure a
        !> continuous path
        !--------------------------------------------------------------
        subroutine get_inter_detector_param(
     $       prev_icoord,
     $       next_icoord,
     $       interior_x_map,
     $       interior_y_map,
     $       icoord_icr,
     $       inter_nb,
     $       x_map_icr,
     $       y_map_icr)

          implicit none

          integer(ikind), dimension(2)             , intent(in)  :: prev_icoord
          integer(ikind), dimension(2)             , intent(in)  :: next_icoord
          real(rkind)   , dimension(nx)            , intent(in)  :: interior_x_map
          real(rkind)   , dimension(ny)            , intent(in)  :: interior_y_map
          real(rkind)   , dimension(2)             , intent(out) :: icoord_icr
          integer                                  , intent(out) :: inter_nb
          real(rkind)   , dimension(:), allocatable, intent(out) :: x_map_icr
          real(rkind)   , dimension(:), allocatable, intent(out) :: y_map_icr

          integer :: i_change
          integer :: j_change
          integer :: i_inter_nb
          integer :: j_inter_nb

          integer :: size_x_map_icr
          integer :: size_y_map_icr


          !determination of the (x,y)-index differences
          i_change = next_icoord(1) - prev_icoord(1)
          j_change = next_icoord(2) - prev_icoord(2)

          !determine the number of additional detectors
          !to be added b/w the prev and the next detectors
          i_inter_nb = abs(i_change)-1
          j_inter_nb = abs(j_change)-1
          inter_nb   = max(0,i_inter_nb,j_inter_nb)

          icoord_icr(1) = 0
          icoord_icr(2) = 0

          if(inter_nb.gt.0) then

             !determination of the increase of the (x,y)-index
             !coordinates
             if(i_change.ne.0) then
                icoord_icr(1) = (real(i_change)-sign(1,i_change))/real(inter_nb)
                !icoord_icr(1) = real(i_change)/real(i_inter_nb+1)
             else
                icoord_icr(1) = 0
             end if

             if(j_change.ne.0) then
                icoord_icr(2) = (real(j_change)-sign(1,j_change))/real(inter_nb)
                !icoord_icr(2) = real(j_change)/real(j_inter_nb+1)
             else
                icoord_icr(2) = 0
             end if

             !allocation of the tables storing the (x,y)-
             !coordinates corresponding of the index
             size_x_map_icr = max(1,i_inter_nb+1)
             size_y_map_icr = max(1,j_inter_nb+1)
             
             allocate(x_map_icr(size_x_map_icr))
             allocate(y_map_icr(size_y_map_icr))

             
             !fill the x_map_icr and the y_map_icr with the
             !coorresponding (x,y)-coordinates
             call determine_local_map_coordinates(
     $            interior_x_map,nx,
     $            size_x_map_icr,
     $            prev_icoord(1),
     $            (i_change.gt.0),
     $            x_map_icr)
             
             call determine_local_map_coordinates(
     $            interior_y_map,ny,
     $            size_y_map_icr,
     $            prev_icoord(2),
     $            (j_change.gt.0),
     $            y_map_icr)

          end if          

        end subroutine get_inter_detector_param


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> from the parameters constraining the addition of
        !> intermediate detectors, give the coordinate of 
        !> the intermediate detector identified by the index k
        !> varying between 1 and total number of detectors to
        !> be added
        !
        !> @date
        !> 24_11_2014 - initial version - J.L. Desmarais
        !
        !>@param prev_coords
        !> general coordinates of the last detector added to the list
        !
        !>@param x_change
        !> change in the x-coordinate between the previous and the next
        !> detector
        !
        !>@param y_change
        !> change in the y-coordinate between the previous and the next
        !> detector
        !
        !>@param k
        !> index identifying the detector to be added
        !
        !>@return inter_coords
        !> general coordinates identifying the detector to be added
        !> to the list
        !--------------------------------------------------------------
        subroutine get_inter_detector_coords(
     $     prev_icoord,
     $     icoord_icr,
     $     k,
     $     x_map_icr,
     $     y_map_icr,
     $     icoord_inter,
     $     rcoord_inter)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: prev_icoord
          real(rkind)   , dimension(2), intent(in)  :: icoord_icr
          integer                     , intent(in)  :: k
          real(rkind)   , dimension(:), intent(in)  :: x_map_icr
          real(rkind)   , dimension(:), intent(in)  :: y_map_icr
          integer(ikind), dimension(2), intent(out) :: icoord_inter
          real(rkind)   , dimension(2), intent(out) :: rcoord_inter

          integer(ikind), dimension(2) :: icoord_local_icr


          !determine the index increase for the (x,y)-index coordinates
          !of the intermediate detectors
          icoord_local_icr(1) = nint(icoord_icr(1)*k)
          icoord_local_icr(2) = nint(icoord_icr(2)*k)

          !determine the (x,y)-index coordinates of the intermediate detectors
          icoord_inter(1) = prev_icoord(1) + icoord_local_icr(1)
          icoord_inter(2) = prev_icoord(2) + icoord_local_icr(2)

          !(x,y)-coordinates of the intermediate detectors are such that
          !the detectors are located at the nodal points corresponding
          !to the (x,y)-index coordinates
          rcoord_inter(1) = x_map_icr(1+abs(icoord_local_icr(1)))
          rcoord_inter(2) = y_map_icr(1+abs(icoord_local_icr(2)))
          
        end subroutine get_inter_detector_coords

      end module bf_detector_module

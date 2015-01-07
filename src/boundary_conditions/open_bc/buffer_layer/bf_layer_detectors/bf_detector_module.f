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
      !              of the intermediate detectors - J.L.Desmarais
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
     $       get_inter_detector_coords

        contains


        !determine the cartesian coordinates corresponding
        !to the index-coordinates of the intermediate
        !detectors
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
             size_outside = 0-min(1,first_icoord)+1

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
             size_outside = (first_icoord+size_local_map-1)-size_interior_map

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
             size_outside = first_icoord-max(first_icoord-size_local_map+1,size_interior_map+1)+1

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
             size_outside = -min(1,first_icoord-size_local_map+1)+1

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

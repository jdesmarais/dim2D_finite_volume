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
      !----------------------------------------------------------------
      module bf_detector_module

        use parameters_kind, only :
     $        ikind,
     $        rkind

        implicit none

        private
        public :: get_inter_detector_param,
     $            get_inter_detector_coords

        contains


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
     $     prev_icoord,
     $     prev_rcoord,
     $     next_icoord,
     $     next_rcoord,
     $     icoord_icr,
     $     rcoord_icr,
     $     inter_nb)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: prev_icoord
          real(rkind)   , dimension(2), intent(in)  :: prev_rcoord
          integer(ikind), dimension(2), intent(in)  :: next_icoord
          real(rkind)   , dimension(2), intent(in)  :: next_rcoord
          real(rkind)   , dimension(2), intent(out) :: icoord_icr
          real(rkind)   , dimension(2), intent(out) :: rcoord_icr
          integer                     , intent(out) :: inter_nb

       
          integer                   :: i_change
          integer                   :: j_change
          integer                   :: i_inter_nb
          integer                   :: j_inter_nb


          i_change = next_icoord(1) - prev_icoord(1)
          j_change = next_icoord(2) - prev_icoord(2)

          i_inter_nb = abs(i_change)-1
          j_inter_nb = abs(j_change)-1
          inter_nb   = max(0,i_inter_nb,j_inter_nb)

          icoord_icr(1) = 0
          icoord_icr(2) = 0
          rcoord_icr(1) = 0.0
          rcoord_icr(2) = 0.0

          if(inter_nb.gt.0) then

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

             rcoord_icr(1) = (next_rcoord(1)-prev_rcoord(1))/(inter_nb+1)
             rcoord_icr(2) = (next_rcoord(2)-prev_rcoord(2))/(inter_nb+1)

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
     $     prev_rcoord,
     $     icoord_icr,
     $     rcoord_icr,
     $     k,
     $     icoord_inter,
     $     rcoord_inter)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: prev_icoord
          real(rkind)   , dimension(2), intent(in)  :: prev_rcoord
          real(rkind)   , dimension(2), intent(in)  :: icoord_icr
          real(rkind)   , dimension(2), intent(in)  :: rcoord_icr
          integer                     , intent(in)  :: k
          integer(ikind), dimension(2), intent(out) :: icoord_inter
          real(rkind)   , dimension(2), intent(out) :: rcoord_inter


          if(rkind.eq.8) then
             icoord_inter(1) = prev_icoord(1) + idnint(icoord_icr(1)*k)
             icoord_inter(2) = prev_icoord(2) + idnint(icoord_icr(2)*k)
             rcoord_inter(1) = prev_rcoord(1) + rcoord_icr(1)*k
             rcoord_inter(2) = prev_rcoord(2) + rcoord_icr(2)*k
          else
             icoord_inter(1) = prev_icoord(1) + nint(icoord_icr(1)*k)
             icoord_inter(2) = prev_icoord(2) + nint(icoord_icr(2)*k)
             rcoord_inter(1) = prev_rcoord(1) + rcoord_icr(1)*k
             rcoord_inter(2) = prev_rcoord(2) + rcoord_icr(2)*k
          end if
          
        end subroutine get_inter_detector_coords

      end module bf_detector_module

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
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module bf_detector_module

        use parameters_kind    , only : ikind, rkind

        implicit none

        private
        public :: get_inter_dct_param,
     $            get_inter_dct_coords

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
        subroutine get_inter_dct_param(
     $     prev_coords, next_coords,
     $     x_change, y_change, inter_nb)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: prev_coords
          integer(ikind), dimension(2), intent(in)  :: next_coords
          real(rkind)                 , intent(out) :: x_change
          real(rkind)                 , intent(out) :: y_change
          integer                     , intent(out) :: inter_nb

          integer :: i_change, j_change

          i_change = next_coords(1) - prev_coords(1)
          j_change = next_coords(2) - prev_coords(2)
          inter_nb = max(0, abs(i_change)-1, abs(j_change)-1)

          if(inter_nb.gt.0) then
             if(i_change.ne.0) then
                x_change = (real(i_change)-sign(1,i_change))/real(inter_nb)
             else
                x_change = 0
             end if
             if(j_change.ne.0) then
                y_change = (real(j_change)-sign(1,j_change))/real(inter_nb)
             else
                y_change = 0
             end if
          else
             x_change = 1
             y_change = 1
          end if

        end subroutine get_inter_dct_param


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
        !> 27_06_2014 - initial version - J.L. Desmarais
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
        function get_inter_dct_coords(
     $     prev_coords,
     $     x_change, y_change, k)
     $     result(inter_coords)

          implicit none

          integer(ikind), dimension(2), intent(in) :: prev_coords
          real(rkind)                 , intent(in) :: x_change
          real(rkind)                 , intent(in) :: y_change
          integer                     , intent(in) :: k
          integer(ikind), dimension(2)             :: inter_coords

          
          if(rkind.eq.4) then
             inter_coords(1) = prev_coords(1) + nint(x_change*k)
             inter_coords(2) = prev_coords(2) + nint(y_change*k)
          else
             inter_coords(1) = prev_coords(1) + idnint(x_change*k)
             inter_coords(2) = prev_coords(2) + idnint(y_change*k)
          end if

        end function get_inter_dct_coords

      end module bf_detector_module

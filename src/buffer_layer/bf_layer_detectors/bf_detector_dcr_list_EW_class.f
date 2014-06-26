      !> @file
      !> module implementing the temporary object used to reorganize
      !> the increasing detector list when either an east or west
      !> buffer layer is removed
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the temporary object used to reorganize
      !> the increasing detector list when either an east or west
      !> buffer layer is removed
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module bf_detector_dcr_list_EW_class

        use bf_detector_dcr_list_class, only : bf_detector_dcr_list
        use parameters_constant, only : y_direction
        use parameters_kind    , only : ikind

        implicit none

        private
        public :: bf_detector_dcr_list_EW

        !> @class bf_detector_dcr_list_EW
        !> class encapsulating the temporary variables needed
        !> to create a new detector list out of the previous
        !> detector list when an east or west buffer layer is
        !> removed
        !----------------------------------------------------------------
        type, abstract, extends(bf_detector_dcr_list) :: bf_detector_dcr_list_EW

          contains

          procedure, pass :: get_detector_changes
          procedure, pass :: compute_new_list

        end type bf_detector_dcr_list_EW


        contains


        subroutine get_detector_changes(
     $       this, dct_list,
     $       nb_added_detectors, nb_deleted_detectors,
     $       sign_added_detectors)
        
          implicit none
            
          class(bf_detector_dcr_list_EW) , intent(in) :: this
          integer(ikind), dimension(:,:) , intent(in) :: dct_list
          integer                        , intent(out):: nb_added_detectors
          integer                        , intent(out):: nb_deleted_detectors
          integer                        , intent(out):: sign_added_detectors

          call this%get_detector_changes_g(
     $         dct_list, y_direction,
     $         nb_added_detectors, nb_deleted_detectors,
     $         sign_added_detectors)
          
        end subroutine get_detector_changes


        subroutine compute_new_list(
     $     this, detector_list,
     $     first_pt_linked, last_pt_linked)
            
          implicit none
            
          class(bf_detector_dcr_list_EW)              , intent(inout) :: this
          integer(ikind), dimension(:,:), allocatable , intent(inout) :: detector_list
          integer(ikind), dimension(2)                , intent(in)    :: first_pt_linked
          integer(ikind), dimension(2)                , intent(in)    :: last_pt_linked

          call this%compute_new_list_g(
     $         detector_list,
     $         first_pt_linked,
     $         last_pt_linked,
     $         y_direction)

        end subroutine compute_new_list

      end module bf_detector_dcr_list_EW_class

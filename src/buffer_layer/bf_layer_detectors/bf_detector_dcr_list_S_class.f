      !> @file
      !> module implementing the temporary object used to reorganize
      !> the increasing detector list when a south buffer layer is
      !> removed
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the temporary object used to reorganize
      !> the increasing detector list when a south buffer layer is
      !> removed
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module bf_detector_dcr_list_S_class

        use bf_detector_dcr_list_NS_class, only : bf_detector_dcr_list_NS
        use parameters_bf_layer          , only : dct_icr_S_default
        use parameters_input             , only : bc_size
        use parameters_kind              , only : ikind

        implicit none


        !> @class bf_detector_dcr_list_S
        !> class encapsulating the temporary variables needed
        !> to create a new detector list out of the previous
        !> detector list when a south buffer layer is removed
        !----------------------------------------------------------------
        type, extends(bf_detector_dcr_list_NS) :: bf_detector_dcr_list_S

          contains

          procedure, nopass :: should_be_removed
          procedure, nopass :: get_border_detector

        end type bf_detector_dcr_list_S


        contains


        function should_be_removed(bf_align, g_coords) result(remove)
          
          implicit none

          integer(ikind), dimension(2,2), intent(in) :: bf_align
          integer(ikind), dimension(2)  , intent(in) :: g_coords
          logical                                    :: remove
          

          remove = ((g_coords(1).ge.(bf_align(1,1)-bc_size)).and.
     $              (g_coords(1).le.(bf_align(1,2)+bc_size)))
          remove = remove.and.(g_coords(2).lt.dct_icr_S_default)

        end function should_be_removed


        function get_border_detector(g_coords) result(dct_g_coords)

          implicit none

          integer(ikind), dimension(2), intent(in) :: g_coords
          integer(ikind), dimension(2)             :: dct_g_coords

          dct_g_coords(1) = g_coords(1)
          dct_g_coords(2) = dct_icr_S_default

        end function get_border_detector

      end module bf_detector_dcr_list_S_class

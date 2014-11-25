      !> @file
      !> module implementing the temporary object used to reorganize
      !> the increasing detector list when a buffer layer is removed
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the temporary object used to reorganize
      !> the increasing detector list when a buffer layer is removed
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module bf_detector_dcr_list_class

        use bf_detector_module, only : get_inter_dct_param,
     $                                 get_inter_dct_coords
        use parameters_input  , only : bc_size
        use parameters_kind   , only : ikind, rkind

        implicit none

        private
        public :: bf_detector_dcr_list


        !> @class bf_detector_dcr_list
        !> class encapsulating the temporary variables needed
        !> to create a new detector list out of the previous
        !> detector list when a buffer layer is removed
        !> \image html bf_detector_dcr_list_class.png
        !> \image latex bf_detector_dcr_list_class.eps
        !
        !> @param segment
        !> indices identifying the segment of detectors removed
        !> from the detector list
        !
        !> @param segment_i
        !> integer identifying whether the index to be saved
        !> correspond ot the first or last point of the segment
        !
        !> @param first_detector
        !> general coordinates identifying the first detector
        !> of the new detector list
        !
        !> @param last_detector
        !> general coordinates identifying the last detector
        !> of the new detector list
        !
        !> @param ini
        !> initialize the object
        !
        !>@param compute_new_list_param
        !> compute the segment of detectors removed from the detector
        !> list as well as the first and last detector points for the list
        !> the new detector list
        !
        !>@param get_first_detector
        !> get the first_detector attribute
        !
        !>@param get_last_detector
        !> get the last_detector attribute
        !
        !>@param set_first_detector
        !> set the first_detector attribute
        !
        !>@param set_last_detector
        !> set the last_detector attribute
        !
        !>@param add_deleted_detector
        !> add the index of the detector to be removed
        !> in the segment
        !
        !>@param compute_new_list_g
        !> compute the new list of detectors
        !
        !>@param get_detector_changes_g
        !> compute the parameters determining the size of
        !> the new detector list
        !
        !>@param should_be_removed
        !> determine whether the detector from the old list of
        !> detectors should be removed from the list for the new
        !> list
        !
        !>@param get_border_detector
        !> determine the general coordinates of the first or last
        !> detector in the new list from the general coordinates
        !> of the current detector
        !
        !>@param get_detector_changes
        !> determine the size of the new detector list without
        !> considering the eventual extra points needed to link the
        !> first and the last point of the detector list with the
        !> neighboring detector lists (to get a closed path of detectors)
        !
        !>@param compute_new_list
        !> compute the new list of detectors
        !-------------------------------------------------------
        type, abstract :: bf_detector_dcr_list

          integer, dimension(2), private :: segment
          integer              , private :: segment_i

          integer(ikind), dimension(2), private :: first_detector
          integer(ikind), dimension(2), private :: last_detector

          contains

          procedure, pass :: ini
          procedure, pass :: compute_new_list_param
          procedure, pass :: get_first_detector
          procedure, pass :: get_last_detector
          procedure, pass :: set_first_detector
          procedure, pass :: set_last_detector
          
          procedure, pass, private :: add_deleted_detector
          procedure, pass          :: compute_new_list_g
          procedure, pass          :: get_detector_changes_g

          procedure(should_be_removed_proc), nopass, private, deferred :: should_be_removed
          procedure(get_bdetector_proc)    , nopass, private, deferred :: get_border_detector
          procedure(get_dct_changes_proc)  ,   pass, private, deferred :: get_detector_changes
          procedure(compute_new_list_proc) ,   pass         , deferred :: compute_new_list

        end type bf_detector_dcr_list


        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> determine whether the detector from the old list of
          !> detectors should be removed from the list for the new
          !> list
          !
          !> @date
          !> 18_06_2014 - initial version - J.L. Desmarais
          !
          !>@param g_coords
          !> general coordinates identifying the position of the detector
          !
          !>@return remove
          !> logical stating whether the detector should be removed
          !--------------------------------------------------------------
          function should_be_removed_proc(bf_align, g_coords) result(remove)
          
            import ikind
          
            implicit none

            integer(ikind), dimension(2,2), intent(in) :: bf_align
            integer(ikind), dimension(2)  , intent(in) :: g_coords
            logical                                    :: remove

          end function should_be_removed_proc
        end interface

        
        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> determine the general coordinates of the first or last
          !> detector in the new list from the general coordinates
          !> of the current detector
          !
          !> @date
          !> 18_06_2014 - initial version - J.L. Desmarais
          !
          !>@param g_coords
          !> general coordinates identifying the position of the old
          !> detector
          !
          !>@param dct_g_coords
          !> general coordinates identifying the position of the new
          !> detector
          !--------------------------------------------------------------          
          function get_bdetector_proc(g_coords) result(dct_g_coords)

            import ikind

            implicit none

            integer(ikind), dimension(2), intent(in) :: g_coords
            integer(ikind), dimension(2)             :: dct_g_coords            

          end function get_bdetector_proc
        end interface


        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> determine the size of the new detector list without
          !> considering the eventual extra points needed to link the
          !> first and the last point of the detector list with the
          !> neighboring detector lists (to get a closed path of detectors)
          !
          !> @date
          !> 18_06_2014 - initial version - J.L. Desmarais
          !
          !>@param this
          !> bf_detector_dcr_list object encapsulating the temporary
          !> parameters needed to construct a new detector list after
          !> the removal of a buffer layer
          !
          !>@param dct_list
          !> old detector list
          !
          !>@param nb_added_detectors
          !> number of detectors that should be added to create the new
          !> detector list
          !
          !>@param nb_deleted_detectors
          !> number of detectors that should be removed to create the new
          !> detector list
          !
          !>@param sign_added_detectors
          !> sign identifying whether the detectors added to the new list
          !> correspond to the increasing or the decreasing direction
          !--------------------------------------------------------------          
          subroutine get_dct_changes_proc(
     $       this, dct_list,
     $       nb_added_detectors, nb_deleted_detectors,
     $       sign_added_detectors)

            import bf_detector_dcr_list
            import ikind

            implicit none
            
            class(bf_detector_dcr_list)   , intent(in) :: this
            integer(ikind), dimension(:,:), intent(in) :: dct_list
            integer                       , intent(out):: nb_added_detectors
            integer                       , intent(out):: nb_deleted_detectors
            integer                       , intent(out):: sign_added_detectors
            
          end subroutine get_dct_changes_proc
        end interface


        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> compute the new list of detectors
          !
          !> @date
          !> 18_06_2014 - initial version - J.L. Desmarais
          !
          !>@param this
          !> bf_detector_dcr_list object encapsulating the temporary
          !> parameters needed to construct a new detector list after
          !> the removal of a buffer layer
          !
          !>@param detector_list
          !> old detector list
          !
          !>@param first_pt_linked
          !> general coordinates corresponding to the first point to which
          !> the detector list should be linked
          !
          !>@param last_pt_linked
          !> general coordinates corresponding to the last point to which
          !> the detector list should be linked
          !--------------------------------------------------------------
          subroutine compute_new_list_proc(
     $       this, detector_list,
     $       first_pt_linked, last_pt_linked)

            import bf_detector_dcr_list
            import ikind
            
            implicit none
            
            class(bf_detector_dcr_list)                , intent(inout) :: this
            integer(ikind), dimension(:,:), allocatable, intent(inout) :: detector_list
            integer(ikind), dimension(2)               , intent(in)    :: first_pt_linked
            integer(ikind), dimension(2)               , intent(in)    :: last_pt_linked

          end subroutine compute_new_list_proc
        end interface


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the object
        !
        !> @date
        !> 18_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_list object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(bf_detector_dcr_list), intent(inout) :: this

          this%segment(1) = 0
          this%segment(2) = 0
          this%segment_i = 1

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add the index of the detector to be removed
        !> in the segment
        !
        !> @date
        !> 18_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_list object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !
        !>@param index
        !> index identifying the detector removed in the old detector
        !> list
        !--------------------------------------------------------------
        subroutine add_deleted_detector(this, index)

          implicit none

          class(bf_detector_dcr_list), intent(inout) :: this
          integer                    , intent(in)    :: index

          this%segment(this%segment_i) = index

          if(this%segment_i.eq.1) then
             this%segment_i=2
          end if

        end subroutine add_deleted_detector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the parameters determining the size of
        !> the new detector list
        !
        !> @date
        !> 18_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_list object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !
        !>@param dct_list
        !> old detector list
        !
        !>@param dir
        !> direction in which the detectors are removed
        !
        !>@param nb_added_detectors
        !> number of detectors that should be added to create the new
        !> detector list
        !
        !>@param nb_deleted_detectors
        !> number of detectors that should be removed to create the new
        !> detector list
        !
        !>@param sign_added_detectors
        !> sign identifying whether the detectors added to the new list
        !> correspond to the increasing or the decreasing direction
        !--------------------------------------------------------------
        subroutine get_detector_changes_g(
     $     this, dct_list, dir,
     $     nb_added_detectors, nb_deleted_detectors,
     $     sign_added_detectors)

          implicit none
            
          class(bf_detector_dcr_list)   , intent(in) :: this
          integer(ikind), dimension(:,:), intent(in) :: dct_list
          integer                       , intent(in) :: dir
          integer                       , intent(out):: nb_added_detectors
          integer                       , intent(out):: nb_deleted_detectors
          integer                       , intent(out):: sign_added_detectors

          integer(ikind) :: min_border
          integer(ikind) :: max_border


          !if there are no detectors to be removed or added
          if(this%segment(1).eq.0) then
             nb_deleted_detectors = 0
             nb_added_detectors = 0

          !otherwise
          else
             !the number of deleted detectors is determined by the length 
             !of the segment removed
             nb_deleted_detectors = this%segment(2) - this%segment(1)+1
             
             !the number of added detectors is determined by the difference
             !in the direction imposed (dir) between the first point just
             !before the segment removed and the point just after the segment
             !removed
             if(this%segment(1).eq.1) then
                min_border = this%first_detector(dir)
             else
                min_border = dct_list(dir,this%segment(1)-1)
                min_border = min_border+1
             end if
             
             if(this%segment(2).eq.size(dct_list,2)) then
                max_border = this%last_detector(dir)
             else
                max_border = dct_list(dir,this%segment(2)+1)
                max_border = max_border-1
             end if
             
             nb_added_detectors   = abs(max_border-min_border+1)
             sign_added_detectors = sign(1, max_border-min_border+1)

          end if

        end subroutine get_detector_changes_g


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the segment of detectors removed from the detector
        !> list as well as the first and last detector points for the list
        !> the new detector list
        !
        !> @date
        !> 18_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_list object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !
        !>@param bf_align
        !> alignment of the buffer layer removed
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param icoords
        !> (x,y)-indices of the detectors before removal of the detectors
        !> located in the buffer layer removed
        !
        !>@param rcoords
        !> (x,y)-coordinates of the detectors before removal of the
        !>  detectors located in the buffer layer removed
        !--------------------------------------------------------------
        subroutine compute_new_list_param(
     $     this,
     $     bf_align,
     $     interior_x_map,
     $     interior_y_map,
     $     icoords,
     $     rcoords)

          implicit none

          class(bf_detector_dcr_list), intent(inout) :: this
          integer(ikind), dimension(2,2), intent(in) :: bf_align
          real(rkind)   , dimension(nx) , intent(in) :: interior_x_map
          real(rkind)   , dimension(ny) , intent(in) :: interior_y_map
          integer(ikind), dimension(:,:), intent(in) :: icoords
          integer(rkind), dimension(:,:), intent(in) :: rcoords


          integer(ikind) :: i


          !determine the first point in the new detector list:
          !if the first detector in the previous list should
          !be removed, the first detector of the new list is
          !computed depending on the cardinal coordinate
          !of the sublayer removed
          if(this%should_be_removed(bf_align, icoords(:,1))) then
             call add_deleted_detector(this, 1)
             this%first_detector = this%get_border_detector(dct_list(:,1))

          !otherwise the first detector remains
          else
             this%first_icoord = icoords(:,1)
             this%last_rcoord  = rcoords(:,1)

          end if


          !determine the borders of the segment of detectors
          !removed from the old detector list
          do i=2, size(dct_list,2)-1
             
             !if the detector should be removed, its index is
             !added to the segment of removed detectors
             if(this%should_be_removed(bf_align, dct_list(:,i))) then
                call add_deleted_detector(this, i)
             end if

          end do


          !determine the last point in the new detector list:
          !if the last detector in the previous list should
          !be removed, the last detector of the new list is
          !computed depending on the cardinal coordinate
          !of the sublayer removed
          if(this%should_be_removed(bf_align, dct_list(:,size(dct_list,2)))) then
             call add_deleted_detector(this, size(dct_list,2))
             this%last_detector = this%get_border_detector(dct_list(:,size(dct_list,2)))

          !otherwise the last detector remains
          else
             this%last_detector = dct_list(:,size(dct_list,2))

          end if 
          
        end subroutine compute_new_list_param


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new list of detectors
        !
        !> @date
        !> 18_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_list object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !
        !>@param detector_list
        !> old detector list
        !
        !>@param first_pt_linked
        !> general coordinates corresponding to the first point to which
        !> the detector list should be linked
        !
        !>@param last_pt_linked
        !> general coordinates corresponding to the last point to which
        !> the detector list should be linked
        !
        !>@param dir
        !> direction in which the detectors are removed
        !--------------------------------------------------------------
        subroutine compute_new_list_g(
     $     this, detector_list,
     $     first_pt_linked, last_pt_linked,
     $     dir)

          implicit none

          class(bf_detector_dcr_list)                , intent(inout) :: this
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: detector_list
          integer(ikind), dimension(2)               , intent(in)    :: first_pt_linked
          integer(ikind), dimension(2)               , intent(in)    :: last_pt_linked
          integer                                    , intent(in)    :: dir

          integer                                     :: nb_added_detectors
          integer                                     :: nb_deleted_detectors
          integer                                     :: sign_added_detectors
          real(rkind)                                 :: x_change1
          real(rkind)                                 :: x_change2
          real(rkind)                                 :: y_change1
          real(rkind)                                 :: y_change2
          integer                                     :: inter_nb1
          integer                                     :: inter_nb2
          integer(ikind), dimension(2)                :: segment_first_pt
          integer(ikind), dimension(2)                :: new_coords
          integer(ikind)                              :: new_size
          integer(ikind), dimension(:,:), allocatable :: new_detector_list
          integer                                     :: i,j


          !compute the size of the new detector list
          !without the additional points for the links
          call this%get_detector_changes(
     $         detector_list,
     $         nb_added_detectors,
     $         nb_deleted_detectors,
     $         sign_added_detectors)         

          !compute the parameters for linking the new detector
          !list with the first and the last points
          call get_inter_dct_param(
     $         first_pt_linked, this%first_detector,
     $         x_change1, y_change1, inter_nb1)
          call get_inter_dct_param(
     $         this%last_detector, last_pt_linked,
     $         x_change2, y_change2, inter_nb2)


          !reallocation is needed only if new points are needed
          !for the links or if the number of added detectors is
          !different from the number of deleted detectors
          if((nb_deleted_detectors.eq.nb_added_detectors).and.
     $       (inter_nb1.eq.0).and.
     $       (inter_nb2.eq.0)) then


             !the detectors are modified only if some detectors
             !are removed
             if(nb_deleted_detectors.ne.0) then
                
                !reallocation is not needed and the detectors removed
                !are simply replaced by the added detectors
                segment_first_pt = get_segment_first_pt(this, detector_list, dir)

                call replace_removed_detectors(
     $               detector_list,
     $               nb_added_detectors,
     $               sign_added_detectors,
     $               segment_first_pt,
     $               dir,
     $               this%segment(1)-1)

             end if

          else

             !reallocation is needed and the new size of the table is
             !computed
             new_size = size(detector_list,2) +
     $                  nb_added_detectors -
     $                  nb_deleted_detectors +
     $                  inter_nb1 +
     $                  inter_nb2

             !allocate the new detector list
             allocate(new_detector_list(2,new_size))


             !fill the new detector list

             !1) link with the first point
             j=0
             do i=1, inter_nb1
                new_coords = get_inter_dct_coords(
     $               first_pt_linked,
     $               x_change1, y_change1, i)
                new_detector_list(:,j+i) = new_coords
             end do


             !2) add the detectors from the list that
             !   should not be removed on the left of
             !   the segment
             j=inter_nb1
             do i=1, this%segment(1)-1
                new_detector_list(:,j+i) = detector_list(:,i)
             end do


             !3) replace the detectors that are removed
             if((this%segment(1)-1).gt.0) then
                j=j+this%segment(1)-1
             end if

             if(nb_deleted_detectors.ne.0) then

                segment_first_pt = get_segment_first_pt(this, detector_list, dir)

                call replace_removed_detectors(
     $               new_detector_list,
     $               nb_added_detectors,
     $               sign_added_detectors,
     $               segment_first_pt,
     $               dir,
     $               j)

                j=j+nb_added_detectors

             end if


             !4) add the detectors from the list that
             !   should not be removed on the right of
             !   the segment
             do i=1, size(detector_list,2)-this%segment(2)               
                new_detector_list(:,j+i) = detector_list(:,this%segment(2)+i)
             end do
             
             
             !5) link with the last point
             if((size(detector_list,2)-this%segment(2)).gt.0) then
                j=j+size(detector_list,2)-this%segment(2)
             end if

             do i=1, inter_nb2
                new_coords = get_inter_dct_coords(
     $               this%last_detector,
     $               x_change2, y_change2, i)
                new_detector_list(:,j+i) = new_coords
             end do

             
             !set the new detector list
             call MOVE_ALLOC(new_detector_list, detector_list)

          end if

        end subroutine compute_new_list_g


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the first_detector attribute
        !
        !> @date
        !> 18_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_list object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !
        !>@param coords
        !> general coordinates of the first_detector
        !--------------------------------------------------------------
        function get_first_detector(this) result(coords)

          implicit none

          class(bf_detector_dcr_list), intent(in) :: this
          integer(ikind), dimension(2)            :: coords

          coords = this%first_detector

        end function get_first_detector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the last_detector attribute
        !
        !> @date
        !> 18_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_list object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !
        !>@param coords
        !> general coordinates of the last_detector
        !--------------------------------------------------------------
        function get_last_detector(this) result(coords)

          implicit none

          class(bf_detector_dcr_list), intent(in) :: this
          integer(ikind), dimension(2)            :: coords

          coords = this%last_detector

        end function get_last_detector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the first_detector attribute
        !
        !> @date
        !> 18_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_list object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !
        !>@param coords
        !> general coordinates of the first_detector
        !--------------------------------------------------------------
        subroutine set_first_detector(this, coords)

          implicit none

          class(bf_detector_dcr_list) , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: coords

          this%first_detector = coords

        end subroutine set_first_detector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the last_detector attribute
        !
        !> @date
        !> 18_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_list object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !
        !>@param coords
        !> general coordinates of the last_detector
        !--------------------------------------------------------------
        subroutine set_last_detector(this, coords)

          implicit none

          class(bf_detector_dcr_list) , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: coords

          this%last_detector = coords

        end subroutine set_last_detector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the 
        !
        !> @date
        !> 18_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_list object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !
        !>@param detector_list
        !> old detector list
        !
        !>@param dir
        !> direction in which the detectors are removed
        !
        !>@return segment_first_pt
        !> general coordinate of the detector replacing the first detector
        !> removed from the segment
        !--------------------------------------------------------------
        function get_segment_first_pt(this, detector_list, dir)
     $     result(segment_first_pt)

          implicit none

          class(bf_detector_dcr_list)   , intent(in) :: this
          integer(ikind), dimension(:,:), intent(in) :: detector_list
          integer                       , intent(in) :: dir
          integer(ikind), dimension(2)               :: segment_first_pt

          if(this%segment(1).eq.1) then
             segment_first_pt = this%first_detector
          else
             segment_first_pt = detector_list(:, this%segment(1)-1)
             segment_first_pt(dir) = segment_first_pt(dir)+1
          end if

        end function get_segment_first_pt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the 
        !
        !> @date
        !> 18_06_2014 - initial version - J.L. Desmarais
        !
        !>@param detector_list
        !> new detector list
        !
        !>@param nb_added_detectors
        !> number of detectors that should be added to create the new
        !> detector list
        !
        !>@param nb_deleted_detectors
        !> number of detectors that should be removed to create the new
        !> detector list
        !
        !>@param sign_added_detectors
        !> sign identifying whether the detectors added to the new list
        !> correspond to the increasing or the decreasing direction
        !
        !>@param segment_first_pt
        !> general coordinate of the detector replacing the first detector
        !> removed from the segment
        !
        !>@param dir
        !> direction in which the detectors are removed
        !
        !>@param match_i
        !> integer indicating where the new detector list shall the
        !> first new detector be added
        !--------------------------------------------------------------
        subroutine replace_removed_detectors(
     $     detector_list,
     $     nb_added_detectors,
     $     sign_added_detectors,
     $     segment_first_pt,
     $     dir,
     $     match_i)
        
          implicit none

          integer(ikind), dimension(:,:), intent(inout) :: detector_list
          integer                       , intent(in)    :: nb_added_detectors
          integer                       , intent(in)    :: sign_added_detectors
          integer(ikind), dimension(2)  , intent(in)    :: segment_first_pt
          integer                       , intent(in)    :: dir
          integer                       , intent(in)    :: match_i

          integer :: i
          integer(ikind), dimension(2) :: new_coords

          do i=1, nb_added_detectors
                
             new_coords      = segment_first_pt
             new_coords(dir) = new_coords(dir)+sign_added_detectors*(i-1)
             
             detector_list(:,match_i+i) = new_coords
             
          end do

        end subroutine replace_removed_detectors

      end module bf_detector_dcr_list_class

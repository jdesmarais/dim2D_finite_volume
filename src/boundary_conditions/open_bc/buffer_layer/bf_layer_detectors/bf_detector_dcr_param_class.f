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
      ! 25_11_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module bf_detector_dcr_param_class

        use bf_detector_module, only :
     $     get_inter_detector_param,
     $     get_inter_detector_coords

        use parameters_bf_layer, only :
     $       dct_icr_N_default,
     $       dct_icr_S_default,
     $       dct_icr_E_default,
     $       dct_icr_W_default,
     $       align_N,
     $       align_S,
     $       align_E,
     $       align_W

        use parameters_constant, only :
     $       N,S,E,W,
     $       x_direction,
     $       y_direction

        use parameters_input, only :
     $       nx,ny,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind


        private
        public ::
     $       bf_detector_dcr_param,
     $       
     $       compute_new_list_param_g,
     $       finalize_new_list_g,
     $       pinpoint_detector_for_removal,
     $       get_border_detector_N,
     $       get_border_detector_S,
     $       get_border_detector_E,
     $       get_border_detector_W,
     $       get_detector_changes,
     $       get_segment_first_pt,
     $       replace_removed_detectors,
     $       should_be_removed_N,
     $       should_be_removed_S,
     $       should_be_removed_E,
     $       should_be_removed_W,
     $       get_rcoord



        !> @class bf_detector_dcr_param
        !> class encapsulating the temporary variables needed
        !> to create a new detector list out of the previous
        !> detector list when a buffer layer is removed
        !> \image html bf_detector_dcr_param_class.png
        !> \image latex bf_detector_dcr_param_class.eps
        !
        !> @param segment_borders
        !> array identifying the index borders of the detectors
        !> removed from the detector list
        !
        !> @param segment_i
        !> integer identifying whether the index to be saved
        !> correspond ot the first or last point of the segment
        !
        !> @param first_icoord
        !> (x,y)-indices of the first detector of the new detector
        !> list
        !
        !> @param first_rcoord
        !> (x,y)-coordinates of the first detector of the new
        !> detector list
        !
        !> @param last_icoord
        !> (x,y)-indices of the last detector of the new detector
        !> list
        !
        !> @param last_rcoord
        !> (x,y)-coordinates of the last detector of the new
        !> detector list
        !
        !> @param ini
        !> initialize the object
        !
        !>@param compute_new_list_param
        !> compute the segment of detectors removed from the detector
        !> list as well as the first and last detector points for the list
        !> the new detector list
        !
        !>@param finalize_new_list
        !> combine the new list of detectors out of multiple
        !> pieces: left detectors linking another list with
        !> the new one, the previous detector list where all
        !> the detectors to be removed are excluded, right
        !> detectors linked the end of the new one with another
        !> one
        !-------------------------------------------------------
        type bf_detector_dcr_param

          integer, dimension(2), private :: segment_borders
          integer              , private :: segment_i

          integer(ikind), dimension(2), private :: first_icoord
          real(rkind)   , dimension(2), private :: first_rcoord
          integer(ikind), dimension(2), private :: last_icoord
          real(rkind)   , dimension(2), private :: last_rcoord

          contains

          procedure, pass :: compute_new_list_param
          procedure, pass :: get_first_detector
          procedure, pass :: get_last_detector
          procedure, pass :: finalize_new_list

          procedure, pass :: get_param !only for tests


        end type bf_detector_dcr_param


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
          !>@param bf_align
          !> edges of the buffer layer
          !
          !>@param g_coords
          !> general coordinates identifying the position of the detector
          !
          !>@return remove
          !> logical stating whether the detector should be removed
          !--------------------------------------------------------------
          function dct_removal_proc(bf_align, g_coords) result(remove)
          
            import ikind
          
            implicit none

            integer(ikind), dimension(2,2), intent(in) :: bf_align
            integer(ikind), dimension(2)  , intent(in) :: g_coords
            logical                                    :: remove

          end function dct_removal_proc


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
          !>@param icoords
          !> (x,y)-indices identifying the position of the old
          !> detector
          !
          !>@param rcoords
          !> (x,y)-coordinates identifying the position of the old
          !> detector
          !
          !>@param interior_map
          !> coordinates of the interior domain
          !--------------------------------------------------------------
          subroutine get_border_dct_proc(
     $       icoord,
     $       rcoord,
     $       interior_map,
     $       icoord_n,
     $       rcoord_n)

            import ikind
            import rkind

            integer(ikind), dimension(2), intent(in)  :: icoord
            real(rkind)   , dimension(2), intent(in)  :: rcoord
            real(rkind)   , dimension(:), intent(in)  :: interior_map
            integer(ikind), dimension(2), intent(out) :: icoord_n
            real(rkind)   , dimension(2), intent(out) :: rcoord_n

          end subroutine get_border_dct_proc

        end interface

        contains

        
        subroutine compute_new_list_param(
     $       this,
     $       bf_localization,
     $       bf_align,
     $       interior_x_map,
     $       interior_y_map,
     $       icoords,
     $       rcoords)

          implicit none

          class(bf_detector_dcr_param)  , intent(inout) :: this
          integer                       , intent(in)    :: bf_localization
          integer(ikind), dimension(2,2), intent(in)    :: bf_align
          real(rkind)   , dimension(nx) , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny) , intent(in)    :: interior_y_map
          integer(ikind), dimension(:,:), intent(in)    :: icoords
          real(rkind)   , dimension(:,:), intent(in)    :: rcoords
 

          !initialize the segment parameters
          this%segment_borders(1) = 0
          this%segment_borders(2) = 0
          this%segment_i          = 1


          !compute the sement borders as well as the first
          !and last detectors
          select case(bf_localization)

            case(N)
               call compute_new_list_param_g(
     $              this,
     $              bf_align,
     $              interior_y_map,
     $              icoords,
     $              rcoords,
     $              should_be_removed_N,
     $              get_border_detector_N)

            case(S)
               call compute_new_list_param_g(
     $              this,
     $              bf_align,
     $              interior_y_map,
     $              icoords,
     $              rcoords,
     $              should_be_removed_S,
     $              get_border_detector_S)

            case(E)
               call compute_new_list_param_g(
     $              this,
     $              bf_align,
     $              interior_x_map,
     $              icoords,
     $              rcoords,
     $              should_be_removed_E,
     $              get_border_detector_E)

            case(W)
               call compute_new_list_param_g(
     $              this,
     $              bf_align,
     $              interior_x_map,
     $              icoords,
     $              rcoords,
     $              should_be_removed_W,
     $              get_border_detector_W)

            case default
               print '(''bf_detector_dcr_param_class'')'
               print '(''compute_new_list_param'')'
               print '(''bf_localization not recognized: '',I2)',
     $              bf_localization
               stop ''

           end select

        end subroutine compute_new_list_param


        !compute the new list of detectors
        subroutine finalize_new_list(
     $     this,
     $     bf_localization,
     $     interior_x_map,
     $     interior_y_map,
     $     icoords,
     $     rcoords,
     $     left_icoord,
     $     right_icoord,
     $     remove_first_dct,
     $     remove_last_dct)

          implicit none

          class(bf_detector_dcr_param)               , intent(inout) :: this
          integer                                    , intent(in)    :: bf_localization
          real(rkind)   , dimension(nx)              , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)              , intent(in)    :: interior_y_map
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: icoords
          real(rkind)   , dimension(:,:), allocatable, intent(inout) :: rcoords
          integer(ikind), dimension(2)               , intent(in)    :: left_icoord
          integer(ikind), dimension(2)               , intent(in)    :: right_icoord
          logical                                    , intent(in)    :: remove_first_dct
          logical                                    , intent(in)    :: remove_last_dct


          select case(bf_localization)

            case(N,S)

               call finalize_new_list_g(
     $              this,
     $              icoords,
     $              rcoords,
     $              left_icoord,
     $              right_icoord,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_x_map,
     $              x_direction,
     $              remove_first_dct,
     $              remove_last_dct)

            case(E,W)
               
               call finalize_new_list_g(
     $              this,
     $              icoords,
     $              rcoords,
     $              left_icoord,
     $              right_icoord,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_y_map,
     $              y_direction,
     $              remove_first_dct,
     $              remove_last_dct)

            case default
               print '(''bf_detector_dcr_param_class'')'
               print '(''compute_new_list_param'')'
               print '(''bf_localization not recognized: '',I2)',
     $              bf_localization
               stop ''

           end select

        end subroutine finalize_new_list


!======================================================================
!general versions of the procedures for computing the parameters
!determing the new detector list and the determination of the new
!detector list
!======================================================================

        !general subroutine for determining the parameters for
        !the removal of the detectors belonging to the buffer
        !layer removed
        subroutine compute_new_list_param_g(
     $     this,
     $     bf_align,
     $     interior_map,
     $     icoords,
     $     rcoords,
     $     should_be_removed,
     $     get_border_detector)

          implicit none

          class(bf_detector_dcr_param)  , intent(inout) :: this
          integer(ikind), dimension(2,2), intent(in)    :: bf_align
          real(rkind)   , dimension(:)  , intent(in)    :: interior_map
          integer(ikind), dimension(:,:), intent(in)    :: icoords
          real(rkind)   , dimension(:,:), intent(in)    :: rcoords
          procedure(dct_removal_proc)                   :: should_be_removed
          procedure(get_border_dct_proc)                :: get_border_detector
          
          integer :: i
          integer :: nb_dct


          nb_dct = size(icoords,2)


          !1) check whether the first detector should be removed
          !============================================================
          !if the first detector in the previous list should be
          !removed, the first detector of the new list is computed
          !depending on the cardinal coordinate of the sublayer removed
          !============================================================
          !if the first detector is removed
          if(should_be_removed(bf_align, icoords(:,1))) then

             !pinpoint the first detector for removal
             call pinpoint_detector_for_removal(this,1)

             !determine the first detector that will replace
             !the detector to be removed
             call get_border_detector(
     $            icoords(:,1),
     $            rcoords(:,1),
     $            interior_map,
     $            this%first_icoord,
     $            this%first_rcoord)

          !if the first detector remains
          else

             this%first_icoord = icoords(:,1)
             this%first_rcoord = rcoords(:,1)

          end if

          
          !2) check whether the detectors should be removed except the
          !   first and the last detectors
          !============================================================
          do i=2, nb_dct-1
             
             if(should_be_removed(bf_align, icoords(:,i))) then
                call pinpoint_detector_for_removal(this, i)
             end if

          end do

          
          !3) check whether the last detector should be removed
          !============================================================
          !if the last detector is removed
          if(should_be_removed(bf_align, icoords(:,nb_dct))) then

             !pinpoint the last detector for removal
             call pinpoint_detector_for_removal(this, nb_dct)

             !determine the last detector that will replace
             !the detector to be removed
             call get_border_detector(
     $            icoords(:,nb_dct),
     $            rcoords(:,nb_dct),
     $            interior_map,
     $            this%last_icoord,
     $            this%last_rcoord)

          !if the last detector remains
          else

             this%last_icoord = icoords(:,nb_dct)
             this%last_rcoord = rcoords(:,nb_dct)

          end if


        end subroutine compute_new_list_param_g


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add the index of the detector to be removed
        !> in the segment
        !
        !> @date
        !> 25_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_param object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !
        !>@param index
        !> index identifying the detector removed in the old detector
        !> list
        !--------------------------------------------------------------
        subroutine pinpoint_detector_for_removal(this, index)

          implicit none

          class(bf_detector_dcr_param), intent(inout) :: this
          integer                    , intent(in)    :: index

          this%segment_borders(this%segment_i) = index

          if(this%segment_i.eq.1) then
             this%segment_i=2
          end if

        end subroutine pinpoint_detector_for_removal


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
        !> bf_detector_dcr_param object encapsulating the temporary
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
        subroutine get_detector_changes(
     $     this,
     $     icoords,
     $     dir,
     $     nb_added_detectors,
     $     nb_deleted_detectors,
     $     sign_added_detectors)

          implicit none
            
          class(bf_detector_dcr_param)  , intent(in) :: this
          integer(ikind), dimension(:,:), intent(in) :: icoords
          integer                       , intent(in) :: dir
          integer                       , intent(out):: nb_added_detectors
          integer                       , intent(out):: nb_deleted_detectors
          integer                       , intent(out):: sign_added_detectors

          integer(ikind) :: min_border
          integer(ikind) :: max_border


          !if there are no detectors to be removed or added
          if(this%segment_borders(1).eq.0) then
             nb_deleted_detectors = 0
             nb_added_detectors   = 0
             sign_added_detectors = 0

          !otherwise
          else
             !the number of deleted detectors is determined by the length 
             !of the segment removed
             nb_deleted_detectors =
     $            this%segment_borders(2) -
     $            this%segment_borders(1) + 1
             
             !the number of added detectors is determined by the difference
             !in the direction imposed (dir) between the first point just
             !before the segment removed and the point just after the segment
             !removed
             if(this%segment_borders(1).eq.1) then
                min_border = this%first_icoord(dir)
             else
                min_border = icoords(dir,this%segment_borders(1)-1)
                min_border = min_border+1
             end if
             
             if(this%segment_borders(2).eq.size(icoords,2)) then
                max_border = this%last_icoord(dir)
             else
                max_border = icoords(dir,this%segment_borders(2)+1)
                max_border = max_border-1
             end if
             
             nb_added_detectors   = abs(max_border-min_border+1)
             sign_added_detectors = sign(1, max_border-min_border+1)

          end if

        end subroutine get_detector_changes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new list of detectors
        !
        !> @date
        !> 25_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_dcr_param object encapsulating the temporary
        !> parameters needed to construct a new detector list after
        !> the removal of a buffer layer
        !
        !>@param icoords
        !> (x,y)-indices of the old detectors
        !
        !>@param rcoords
        !> (x,y)-coordinates of the old detectors
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
        !>  (x-direction for N and S)
        !>  (y-direction for E and W)
        !--------------------------------------------------------------
        subroutine finalize_new_list_g(
     $     this,
     $     icoords,
     $     rcoords,
     $     left_icoord,
     $     right_icoord,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_map_dir,
     $     dir,
     $     remove_first_dct,
     $     remove_last_dct)

          implicit none

          class(bf_detector_dcr_param)               , intent(inout) :: this
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: icoords
          real(rkind)   , dimension(:,:), allocatable, intent(inout) :: rcoords
          integer(ikind), dimension(2)               , intent(in)    :: left_icoord
          integer(ikind), dimension(2)               , intent(in)    :: right_icoord
          real(rkind)   , dimension(nx)              , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)              , intent(in)    :: interior_y_map
          real(rkind)   , dimension(:)               , intent(in)    :: interior_map_dir
          integer                                    , intent(in)    :: dir
          logical                                    , intent(in)    :: remove_first_dct
          logical                                    , intent(in)    :: remove_last_dct

          integer                                     :: nb_added_detectors
          integer                                     :: nb_deleted_detectors
          integer                                     :: sign_added_detectors
          real(rkind), dimension(2)                   :: left_icoord_icr
          real(rkind), dimension(2)                   :: right_icoord_icr
          integer                                     :: left_inter_nb
          integer                                     :: right_inter_nb
          real(rkind)   , dimension(:)  , allocatable :: left_x_map_icr
          real(rkind)   , dimension(:)  , allocatable :: left_y_map_icr
          real(rkind)   , dimension(:)  , allocatable :: right_x_map_icr
          real(rkind)   , dimension(:)  , allocatable :: right_y_map_icr
          
          integer(ikind), dimension(:,:), allocatable :: icoords_n
          real(rkind)   , dimension(:,:), allocatable :: rcoords_n

          integer(ikind), dimension(2)                :: seg_icoord
          real(rkind)   , dimension(2)                :: seg_rcoord
          integer(ikind), dimension(2)                :: icoord_inter
          real(rkind)   , dimension(2)                :: rcoord_inter
          integer(ikind)                              :: new_size
          integer                                     :: i,j


          !compute the size of the new detector list
          !without the additional points for the links
          call get_detector_changes(
     $         this,
     $         icoords,
     $         dir,
     $         nb_added_detectors,
     $         nb_deleted_detectors,
     $         sign_added_detectors)

          !compute the parameters for linking the new detector
          !list with the first and the last points
          call get_inter_detector_param(
     $         left_icoord,
     $         this%first_icoord,
     $         interior_x_map,
     $         interior_y_map,
     $         left_icoord_icr,
     $         left_inter_nb,
     $         left_x_map_icr,
     $         left_y_map_icr)

          call get_inter_detector_param(
     $         this%last_icoord,
     $         right_icoord,
     $         interior_x_map,
     $         interior_y_map,
     $         right_icoord_icr,
     $         right_inter_nb,
     $         right_x_map_icr,
     $         right_y_map_icr)


          !reallocation is needed only if new points are needed
          !for the links or if the number of added detectors is
          !different from the number of deleted detectors
          if((nb_deleted_detectors.eq.nb_added_detectors).and.
     $       (left_inter_nb.eq.0).and.
     $       (right_inter_nb.eq.0)) then


             !the detectors are modified only if some detectors
             !are removed
             if(nb_deleted_detectors.ne.0) then
                
                !reallocation is not needed and the detectors removed
                !are simply replaced by the added detectors
                call get_segment_first_pt(
     $               this,
     $               icoords,
     $               rcoords,
     $               interior_map_dir,
     $               dir,
     $               seg_icoord,
     $               seg_rcoord)

                call replace_removed_detectors(
     $               icoords,
     $               rcoords,
     $               nb_added_detectors,
     $               sign_added_detectors,
     $               seg_icoord,
     $               seg_rcoord,
     $               interior_map_dir,
     $               dir,
     $               this%segment_borders(1)-1)

             end if

          else

             !reallocation is needed and the new size of the table is
             !computed
             new_size = size(icoords,2) +
     $                  nb_added_detectors -
     $                  nb_deleted_detectors +
     $                  left_inter_nb +
     $                  right_inter_nb

             !allocate the new detector lists
             allocate(icoords_n(2,new_size))
             allocate(rcoords_n(2,new_size))

             !fill the new detector list

             !1) link with the first point
             j=0

             do i=1, left_inter_nb
                   
                call get_inter_detector_coords(
     $               left_icoord,
     $               left_icoord_icr,
     $               i,
     $               left_x_map_icr,
     $               left_y_map_icr,
     $               icoord_inter,
     $               rcoord_inter)
                   
                icoords_n(:,j+i) = icoord_inter
                rcoords_n(:,j+i) = rcoord_inter
                
             end do

             deallocate(left_x_map_icr)
             deallocate(left_y_map_icr)


             !2) add the detectors from the list that
             !   should not be removed on the left of
             !   the segment
             j=left_inter_nb
             do i=1, this%segment_borders(1)-1
                icoords_n(:,j+i) = icoords(:,i)
                rcoords_n(:,j+i) = rcoords(:,i)
             end do


             !3) replace the detectors that are removed
             if((this%segment_borders(1)-1).gt.0) then
                j=j+this%segment_borders(1)-1
             end if

             if(nb_deleted_detectors.ne.0) then

                call get_segment_first_pt(
     $               this,
     $               icoords,
     $               rcoords,
     $               interior_map_dir,
     $               dir,
     $               seg_icoord,
     $               seg_rcoord)

                call replace_removed_detectors(
     $               icoords_n,
     $               rcoords_n,
     $               nb_added_detectors,
     $               sign_added_detectors,
     $               seg_icoord,
     $               seg_rcoord,
     $               interior_map_dir,
     $               dir,
     $               j)

                j=j+nb_added_detectors

             end if


             !4) add the detectors from the list that
             !   should not be removed on the right of
             !   the segment
             do i=1, size(icoords,2)-this%segment_borders(2)
                icoords_n(:,j+i) = icoords(:,this%segment_borders(2)+i)
                rcoords_n(:,j+i) = rcoords(:,this%segment_borders(2)+i)
             end do
             
             
             !5) link with the last point
             if((size(icoords,2)-this%segment_borders(2)).gt.0) then
                j=j+size(icoords,2)-this%segment_borders(2)
             end if

             do i=1, right_inter_nb

                call get_inter_detector_coords(
     $               this%last_icoord,
     $               right_icoord_icr,
     $               i,
     $               right_x_map_icr,
     $               right_y_map_icr,
     $               icoord_inter,
     $               rcoord_inter)

                icoords_n(:,j+i) = icoord_inter
                rcoords_n(:,j+i) = rcoord_inter
                
             end do

             deallocate(right_x_map_icr)
             deallocate(right_y_map_icr)
             
             !set the new detector list
             call MOVE_ALLOC(icoords_n,icoords)
             call MOVE_ALLOC(rcoords_n,rcoords)             

          end if

          if(remove_first_dct.or.remove_last_dct) then

             if(remove_first_dct) then
                
                if(.not.remove_last_dct) then
                   allocate(icoords_n(2,size(icoords,2)-1))
                   allocate(rcoords_n(2,size(rcoords,2)-1))
                   icoords_n = icoords(:,2:size(icoords,2))
                   rcoords_n = rcoords(:,2:size(rcoords,2))

                else
                   allocate(icoords_n(2,size(icoords,2)-2))
                   allocate(rcoords_n(2,size(rcoords,2)-2))
                   icoords_n = icoords(:,2:size(icoords,2)-1)
                   rcoords_n = rcoords(:,2:size(rcoords,2)-1)

                end if

             else

                allocate(icoords_n(2,size(icoords,2)-1))
                allocate(rcoords_n(2,size(rcoords,2)-1))
                icoords_n = icoords(:,1:size(icoords,2)-1)
                rcoords_n = rcoords(:,1:size(rcoords,2)-1)
                
             end if

             call MOVE_ALLOC(icoords_n,icoords)
             call MOVE_ALLOC(rcoords_n,rcoords)

          end if

       end subroutine finalize_new_list_g


!======================================================================
!procedures specific to a buffer layer for the determination
!of the new parameters for the detector list
!======================================================================

        function should_be_removed_N(bf_align, g_coords) result(remove)
          
          implicit none

          integer(ikind), dimension(2,2), intent(in) :: bf_align
          integer(ikind), dimension(2)  , intent(in) :: g_coords
          logical                                    :: remove
          

          remove = ((g_coords(1).ge.(bf_align(1,1)-bc_size)).and.
     $              (g_coords(1).le.(bf_align(1,2)+bc_size)))
          remove = remove.and.(g_coords(2).gt.dct_icr_N_default)

        end function should_be_removed_N


        subroutine get_border_detector_N(
     $     icoord,
     $     rcoord,
     $     interior_map,
     $     icoord_n,
     $     rcoord_n)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: icoord
          real(rkind)   , dimension(2), intent(in)  :: rcoord
          real(rkind)   , dimension(:), intent(in)  :: interior_map
          integer(ikind), dimension(2), intent(out) :: icoord_n
          real(rkind)   , dimension(2), intent(out) :: rcoord_n

          icoord_n(1) = icoord(1)
          icoord_n(2) = dct_icr_N_default
          
          rcoord_n(1) = rcoord(1)
          rcoord_n(2) = interior_map(dct_icr_N_default)

        end subroutine get_border_detector_N


        function should_be_removed_S(bf_align, g_coords) result(remove)
          
          implicit none

          integer(ikind), dimension(2,2), intent(in) :: bf_align
          integer(ikind), dimension(2)  , intent(in) :: g_coords
          logical                                    :: remove
          

          remove = ((g_coords(1).ge.(bf_align(1,1)-bc_size)).and.
     $              (g_coords(1).le.(bf_align(1,2)+bc_size)))
          remove = remove.and.(g_coords(2).lt.dct_icr_S_default)

        end function should_be_removed_S


        subroutine get_border_detector_S(
     $     icoord,
     $     rcoord,
     $     interior_map,
     $     icoord_n,
     $     rcoord_n)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: icoord
          real(rkind)   , dimension(2), intent(in)  :: rcoord
          real(rkind)   , dimension(:), intent(in)  :: interior_map
          integer(ikind), dimension(2), intent(out) :: icoord_n
          real(rkind)   , dimension(2), intent(out) :: rcoord_n

          icoord_n(1) = icoord(1)
          icoord_n(2) = dct_icr_S_default
          
          rcoord_n(1) = rcoord(1)
          rcoord_n(2) = interior_map(dct_icr_S_default)

        end subroutine get_border_detector_S


        function should_be_removed_E(bf_align, g_coords) result(remove)
          
          implicit none

          integer(ikind), dimension(2,2), intent(in) :: bf_align
          integer(ikind), dimension(2)  , intent(in) :: g_coords
          logical                                    :: remove
          
          logical :: lower_border_condition
          logical :: upper_border_condition


          if(bf_align(2,1).eq.(align_S+1)) then
             lower_border_condition = .true.
          else
             lower_border_condition = g_coords(2).ge.(bf_align(2,1)-bc_size)
          end if

          if(bf_align(2,2).eq.(align_N-1)) then
             upper_border_condition = .true.
          else
             upper_border_condition = g_coords(2).le.(bf_align(2,2)+bc_size)
          end if

          remove = lower_border_condition.and.upper_border_condition
          remove = remove.and.(g_coords(1).gt.dct_icr_E_default)

        end function should_be_removed_E


        subroutine get_border_detector_E(
     $     icoord,
     $     rcoord,
     $     interior_map,
     $     icoord_n,
     $     rcoord_n)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: icoord
          real(rkind)   , dimension(2), intent(in)  :: rcoord
          real(rkind)   , dimension(:), intent(in)  :: interior_map
          integer(ikind), dimension(2), intent(out) :: icoord_n
          real(rkind)   , dimension(2), intent(out) :: rcoord_n

          icoord_n(1) = dct_icr_E_default
          icoord_n(2) = icoord(2)

          rcoord_n(1) = interior_map(dct_icr_E_default)
          rcoord_n(2) = rcoord(2)

        end subroutine get_border_detector_E


        function should_be_removed_W(bf_align, g_coords) result(remove)
          
          implicit none

          integer(ikind), dimension(2,2), intent(in) :: bf_align
          integer(ikind), dimension(2)  , intent(in) :: g_coords
          logical                                    :: remove
          
          logical :: lower_border_condition
          logical :: upper_border_condition

          
          if(bf_align(2,1).eq.(align_S+1)) then
             lower_border_condition = .true.
          else
             lower_border_condition = g_coords(2).ge.(bf_align(2,1)-bc_size)
          end if

          if(bf_align(2,2).eq.(align_N-1)) then
             upper_border_condition = .true.
          else
             upper_border_condition = g_coords(2).le.(bf_align(2,2)+bc_size)
          end if

          remove = lower_border_condition.and.upper_border_condition
          remove = remove.and.(g_coords(1).lt.dct_icr_W_default)

        end function should_be_removed_W


        subroutine get_border_detector_W(
     $     icoord,
     $     rcoord,
     $     interior_map,
     $     icoord_n,
     $     rcoord_n)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: icoord
          real(rkind)   , dimension(2), intent(in)  :: rcoord
          real(rkind)   , dimension(:), intent(in)  :: interior_map
          integer(ikind), dimension(2), intent(out) :: icoord_n
          real(rkind)   , dimension(2), intent(out) :: rcoord_n

          icoord_n(1) = dct_icr_W_default
          icoord_n(2) = icoord(2)

          rcoord_n(1) = interior_map(dct_icr_W_default)
          rcoord_n(2) = rcoord(2)

        end subroutine get_border_detector_W

!=========================================================================
!procedure annexes
!=========================================================================

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
        !> bf_detector_dcr_param object encapsulating the temporary
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
        subroutine get_segment_first_pt(
     $     this,
     $     icoords,
     $     rcoords,
     $     interior_map,
     $     dir,
     $     seg_icoord,
     $     seg_rcoord)

          implicit none

          class(bf_detector_dcr_param)  , intent(in)  :: this
          integer(ikind), dimension(:,:), intent(in)  :: icoords
          real(rkind)   , dimension(:,:), intent(in)  :: rcoords
          real(rkind)   , dimension(:)  , intent(in)  :: interior_map
          integer                       , intent(in)  :: dir
          integer(ikind), dimension(2)  , intent(out) :: seg_icoord
          real(rkind)   , dimension(2)  , intent(out) :: seg_rcoord

          if(this%segment_borders(1).eq.1) then

             seg_icoord = this%first_icoord
             seg_rcoord = this%first_rcoord
             
          else

             seg_icoord = icoords(:, this%segment_borders(1)-1)
             seg_rcoord = rcoords(:, this%segment_borders(1)-1)

             seg_icoord(dir) = seg_icoord(dir)+1
             seg_rcoord(dir) = get_rcoord(seg_icoord(dir), interior_map)

          end if

        end subroutine get_segment_first_pt


        !turn general index into a coordinate
        function get_rcoord(icoord,interior_map)

          implicit none

          integer(ikind)              , intent(in) :: icoord
          real(rkind)   , dimension(:), intent(in) :: interior_map
          real(rkind)                              :: get_rcoord

          real(rkind)    :: ds
          integer(ikind) :: ns


          if(icoord.lt.1) then
             ds = interior_map(2)-interior_map(1)
             get_rcoord = interior_map(1) + (icoord-1)*ds
          else
             if(icoord.le.size(interior_map,1)) then
                get_rcoord = interior_map(icoord)
             else
                ns = size(interior_map,1)
                ds = interior_map(ns) - interior_map(ns-1)
                get_rcoord = interior_map(ns) + (icoord-ns)*ds
             end if
          end if

        end function get_rcoord


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
     $     icoords_n,
     $     rcoords_n,
     $     nb_added_detectors,
     $     sign_added_detectors,
     $     seg_icoord,
     $     seg_rcoord,
     $     interior_map,
     $     dir,
     $     match_i)
        
          implicit none

          integer(ikind), dimension(:,:), intent(inout) :: icoords_n
          real(rkind)   , dimension(:,:), intent(inout) :: rcoords_n
          integer                       , intent(in)    :: nb_added_detectors
          integer                       , intent(in)    :: sign_added_detectors
          integer(ikind), dimension(2)  , intent(in)    :: seg_icoord
          real(rkind)   , dimension(2)  , intent(in)    :: seg_rcoord
          real(rkind)   , dimension(:)  , intent(in)    :: interior_map
          integer                       , intent(in)    :: dir
          integer                       , intent(in)    :: match_i

          integer                      :: i
          integer(ikind), dimension(2) :: icoord_inter
          real(rkind)   , dimension(2) :: rcoord_inter
          

          do i=1, nb_added_detectors
                
             icoord_inter      = seg_icoord
             rcoord_inter      = seg_rcoord

             icoord_inter(dir) = icoord_inter(dir)+sign_added_detectors*(i-1)
             rcoord_inter(dir) = get_rcoord(icoord_inter(dir),interior_map)
             
             icoords_n(:,match_i+i) = icoord_inter
             rcoords_n(:,match_i+i) = rcoord_inter
             
          end do

        end subroutine replace_removed_detectors


        ! get the first detector (x,y)-indices and (x,y)-coordinates
        subroutine get_first_detector(this, icoord, rcoord)

          implicit none

          class(bf_detector_dcr_param), intent(in)  :: this
          integer(ikind), dimension(2), intent(out) :: icoord
          real(rkind)   , dimension(2), intent(out) :: rcoord

          icoord = this%first_icoord
          rcoord = this%first_rcoord

        end subroutine get_first_detector

      
        ! get the last detector (x,y)-indices and (x,y)-coordinates
        subroutine get_last_detector(this, icoord, rcoord)

          implicit none

          class(bf_detector_dcr_param), intent(in)  :: this
          integer(ikind), dimension(2), intent(out) :: icoord
          real(rkind)   , dimension(2), intent(out) :: rcoord

          icoord = this%last_icoord
          rcoord = this%last_rcoord

        end subroutine get_last_detector


        !get the attributes of the objects to comapre with the test
        !data
        subroutine get_param(
     $     this,
     $     segment_borders,
     $     first_icoord,
     $     first_rcoord,
     $     last_icoord,
     $     last_rcoord)

          implicit none

          class(bf_detector_dcr_param), intent(in)  :: this
          integer(ikind), dimension(2), intent(out) :: segment_borders
          integer(ikind), dimension(2), intent(out) :: first_icoord
          real(rkind)   , dimension(2), intent(out) :: first_rcoord
          integer(ikind), dimension(2), intent(out) :: last_icoord
          real(rkind)   , dimension(2), intent(out) :: last_rcoord 

          segment_borders = this%segment_borders
          first_icoord    = this%first_icoord    
          first_rcoord    = this%first_rcoord
          last_icoord     = this%last_icoord 
          last_rcoord     = this%last_rcoord 

        end subroutine get_param

      end module bf_detector_dcr_param_class

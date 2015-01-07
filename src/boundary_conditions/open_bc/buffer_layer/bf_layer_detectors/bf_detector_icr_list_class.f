      !> @file
      !> module implementing the temporary object saving the
      !> general coordinates of a list of new increasing
      !> detectors
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the temporary object saving the
      !> general coordinates of a list of new increasing
      !> detectors
      !
      !> @date
      ! 24_11_2014 - initial version      - J.L. Desmarais
      !----------------------------------------------------------------
      module bf_detector_icr_list_class
      
        use bf_detector_module, only :
     $       get_inter_detector_param,
     $       get_inter_detector_coords

        use parameters_kind, only :
     $       ikind,
     $       rkind
      
        implicit none

        private
        public :: bf_detector_icr_list


        integer, parameter :: nb_detectors_add=5


        !>@class bf_detector_icr_list
        !> object encapuslating the data structure for storing
        !> temporary the new increasing detectors
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying to which mainlayer
        !> the detectors belongs
        !
        !>@param nb_detectors
        !> total number of increasing detectors
        !        
        !>@param icoords
        !> table allocated with the number of detectors at the
        !> previous timestep storing the coordinates of the detectors
        !> as (x,y)-indices (integer(ikind))
        !
        !>@param rcoords
        !> table allocated with the number of detectors at the
        !> previous timestep storing the coordinates of the detectors
        !> as (x,y)-coordinates (real(rkind))
        !
        !>@param ini
        !> initialize the object with the main layer id and the
        !> size of the detector list
        !
        !>@param add_new_detector
        !> add the new detector and its intermediates to the mainlayer
        !> lists to ensure a continuous path
        !
        !>@param add_detector_to_mainlayer
        !> add the new detector general coordinates either to the
        !> detectors_list or the detectors_extra_list
        !
        !>@param get_nb_detectors
        !> get the nb_detectors attribute
        !
        !>@param get_head
        !> get the first detector of the list
        !
        !>@param get_tail
        !> get the last detector of the list
        !
        !>@param fill_new_detector_table
        !> fill the new detector table with the detectors saved
        !> in the bf_detector_icr_list object
        !
        !>@param print_on_matrix
        !> print the coordinates of the detectors saved in the object
        !> as points on a matrix
        !
        !>@param destroy
        !> deallocate the content of the object
        !--------------------------------------------------------------
        type :: bf_detector_icr_list

          integer :: mainlayer_id
          integer, private :: nb_detectors

          integer(ikind), dimension(:,:), allocatable :: icoords
          real(rkind)   , dimension(:,:), allocatable :: rcoords

          contains

          procedure, pass :: ini
          procedure, pass :: add_new_detector

          procedure, pass, private :: add_detector_to_mainlayer
          procedure, pass, private :: make_average_detector

          procedure,   pass :: get_nb_detectors
          procedure,   pass :: get_head
          procedure,   pass :: get_tail
          procedure,   pass :: fill_new_detector_table

          procedure, pass :: print_on_matrix

          procedure, pass :: destroy

        end type bf_detector_icr_list

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the object with the main layer id and the
        !> size of the detector list
        !
        !> @date
        !> 24_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_icr_list object encapsulating the data structure
        !> for temporary storing the position of the new increasing
        !> detectors
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying to which mainlayer
        !> the detectors belongs
        !
        !>@param size_detectors_list
        !> original size guessed for saving all teh new increasing
        !> detectors
        !--------------------------------------------------------------
        subroutine ini(this, mainlayer_id, size_detectors_list)

          implicit none

          class(bf_detector_icr_list), intent(inout) :: this
          integer                    , intent(in)    :: mainlayer_id
          integer(ikind)             , intent(in)    :: size_detectors_list

          this%mainlayer_id = mainlayer_id
          this%nb_detectors = 0

          allocate(this%icoords(2,size_detectors_list))
          allocate(this%rcoords(2,size_detectors_list))

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add the new detector and its intermediates to the mainlayer
        !> lists to ensure a continuous path
        !
        !> @date
        !> 24_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_icr_list object encapsulating the data structure
        !> for temporary storing the position of the new increasing
        !> detectors
        !
        !>@param coords
        !> general coordinates of the new increasing detector added
        !> to the list
        !--------------------------------------------------------------
        subroutine add_new_detector(this, icoord, rcoord)

          implicit none

          class(bf_detector_icr_list) , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: icoord
          real(rkind)   , dimension(2), intent(in)    :: rcoord

          integer(ikind), dimension(2) :: prev_icoord
          real(rkind)   , dimension(2) :: prev_rcoord
          real(rkind)   , dimension(2) :: icoord_icr
          real(rkind)   , dimension(2) :: rcoord_icr
          integer                      :: inter_nb

          real(rkind)   , dimension(:), allocatable :: x_map_icr
          real(rkind)   , dimension(:), allocatable :: y_map_icr

          integer(ikind), dimension(2) :: icoord_inter
          real(rkind)   , dimension(2) :: rcoord_inter
          integer                      :: k
          real(rkind)   , dimension(2) :: rcoord_average

          !if other detectors were saved in the list before,
          !the new detector added should not be too far away
          !from the previous one to retain the closed path
          !figure
          if(this%nb_detectors.gt.0) then

             prev_icoord = this%icoords(:,this%nb_detectors)
             prev_rcoord = this%rcoords(:,this%nb_detectors)

             if((prev_icoord(1).ne.icoord(1)).or.
     $          (prev_icoord(2).ne.icoord(2))) then

                !add intermediate detectors between the previous
                !one and the new one to retain a continuous path
                call get_inter_detector_param(
     $               prev_icoord,
     $               icoord,
     $               interior_x_map,
     $               interior_y_map,
     $               icoord_icr,
     $               inter_nb,
     $               x_map_icr,
     $               y_map_icr)
                
                do k=1, inter_nb
                   call get_inter_detector_coords(
     $                  prev_icoord,
     $                  icoord_icr,
     $                  k,
     $                  x_map_icr,
     $                  y_map_icr,
     $                  icoord_inter,
     $                  rcoord_inter)

                   if((prev_icoord(1).ne.icoord_inter(1)).or.
     $                (prev_icoord(2).ne.icoord_inter(2))) then
                   
                      call add_detector_to_mainlayer(
     $                     this,
     $                     icoord_inter,
     $                     rcoord_inter)

                      prev_icoord = icoord_inter
                      prev_rcoord = rcoord_inter
                      
                   else

                      call this%make_average_detector(prev_rcoord,rcoord_inter)

                   end if

                end do

                if((prev_icoord(1).ne.icoord(1)).or.
     $             (prev_icoord(2).ne.icoord(2))) then

                   call add_detector_to_mainlayer(this, icoord, rcoord)

                else

                   call this%make_average_detector(prev_rcoord,rcoord)

                end if

             !if the detector to be added to the list has the same
             !(x,y)-index coordinates as the previous detector stored
             !in the list, then the previous detector is replaced by
             !a new detector which is an average b/w the two
             else

                call this%make_average_detector(prev_rcoord,rcoord)

             end if

          else

             !add the new detector
             call add_detector_to_mainlayer(this, icoord, rcoord)

          end if

        end subroutine add_new_detector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> if two detectors shared the same (x,y)-index coordinates, an
        !> average detector is created combining the new detector with
        !> the previous detector in the list. This function computes the
        !> (x,y)-coordinates of the average detector and store them in the
        !> list
        !
        !> @date
        !> 07_01_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_icr_list object encapsulating the data structure
        !> for temporary storing the position of the new increasing
        !> detectors
        !
        !>@param prev_rcoord
        !> (x,y)-coordinates of the previous detector in the list sharing
        !> the same (x,y)-index coordinates as the current detector
        !
        !>@param current_rcoord
        !> (x,y)-coordinates of the current detector sharing the same (x,y)-
        !> index coordinates as the previous detector stored in the list
        !--------------------------------------------------------------
        subroutine make_average_detector(this, prev_rcoord, current_rcoord)

          implicit none

          class(bf_detector_icr_list), intent(inout) :: this
          real(rkind), dimension(2)  , intent(in)    :: prev_rcoord
          real(rkind), dimension(2)  , intent(in)    :: current_rcoord

          real(rkind), dimension(2) :: rcoord_average

          if(rkind.eq.8) then
             rcoord_average(1) = 0.5d0*(prev_rcoord(1) + current_rcoord(1))
             rcoord_average(2) = 0.5d0*(prev_rcoord(2) + current_rcoord(2))
          else
             rcoord_average(1) = 0.5*(prev_rcoord(1) + current_rcoord(1))
             rcoord_average(2) = 0.5*(prev_rcoord(2) + current_rcoord(2))
          end if
          
          this%rcoords(:,this%nb_detectors) = rcoord_average

        end subroutine make_average_detector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add the new detector general coordinates either to the
        !> detectors_list or the detectors_extra_list
        !
        !> @date
        !> 24_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_icr_list object encapsulating the data structure
        !> for temporary storing the position of the new increasing
        !> detectors
        !
        !>@param coords
        !> general coordinates of the new increasing detector added
        !> to the list
        !--------------------------------------------------------------
        subroutine add_detector_to_mainlayer(this, icoord, rcoord)

          implicit none

          class(bf_detector_icr_list) , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: icoord
          real(rkind)   , dimension(2), intent(in)    :: rcoord

          integer(ikind), dimension(:,:), allocatable :: tmp_icoords
          real(rkind)   , dimension(:,:), allocatable :: tmp_rcoords


          this%nb_detectors = this%nb_detectors+1

          if(this%nb_detectors.gt.size(this%icoords,2)) then

             allocate(tmp_icoords(2,this%nb_detectors+nb_detectors_add))
             allocate(tmp_rcoords(2,this%nb_detectors+nb_detectors_add))

             tmp_icoords(:,1:size(this%icoords,2)) = this%icoords(:,:)
             tmp_rcoords(:,1:size(this%icoords,2)) = this%rcoords(:,:)

             call MOVE_ALLOC(tmp_icoords,this%icoords)
             call MOVE_ALLOC(tmp_rcoords,this%rcoords)

          end if

          this%icoords(:,this%nb_detectors) = icoord
          this%rcoords(:,this%nb_detectors) = rcoord

        end subroutine add_detector_to_mainlayer

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the nb_detectors attribute
        !
        !> @date
        !> 24_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_icr_list object encapsulating the data structure
        !> for temporary storing the position of the new increasing
        !> detectors
        !--------------------------------------------------------------
        function get_nb_detectors(this)

          implicit none

          class(bf_detector_icr_list), intent(in) :: this
          integer                               :: get_nb_detectors

          get_nb_detectors = this%nb_detectors

        end function get_nb_detectors


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the first detector of the list
        !
        !> @date
        !> 24_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_icr_list object encapsulating the data structure
        !> for temporary storing the position of the new increasing
        !> detectors
        !
        !>@param coords
        !> general coordinates of the detector
        !--------------------------------------------------------------
        subroutine get_head(this,icoords,rcoords)

          implicit none

          class(bf_detector_icr_list) , intent(in)  :: this
          integer(ikind), dimension(2), intent(out) :: icoords
          real(rkind)   , dimension(2), intent(out) :: rcoords

          if(allocated(this%icoords)) then
             icoords = this%icoords(:,1)
             rcoords = this%rcoords(:,1)
          else
             print '(''bf_detector_icr_list_class'')'
             print '(''get_head: arrays not allocated'')'
             stop ''
          end if

        end subroutine get_head


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the last detector of the list
        !
        !> @date
        !> 24_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_icr_list object encapsulating the data structure
        !> for temporary storing the position of the new increasing
        !> detectors
        !
        !>@param coords
        !> general coordinates of the detector
        !--------------------------------------------------------------
        subroutine get_tail(this,icoords,rcoords)

          implicit none

          class(bf_detector_icr_list) , intent(in)  :: this
          integer(ikind), dimension(2), intent(out) :: icoords
          real(rkind)   , dimension(2), intent(out) :: rcoords

          if(allocated(this%icoords)) then

             icoords = this%icoords(:,this%nb_detectors)
             rcoords = this%rcoords(:,this%nb_detectors)

          else

             print '(''bf_detector_icr_list_class'')'
             print '(''get_tail: arrays not allocated'')'
             stop ''

          end if

        end subroutine get_tail


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> fill the new detector table with the detectors saved
        !> in the bf_detector_icr_list object
        !
        !> @date
        !> 24_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_icr_list object encapsulating the data structure
        !> for temporary storing the position of the new increasing
        !> detectors
        !
        !>@param s_index
        !> index where the first detector saved in the bf_detector_icr_list
        !> object should be stored in the new_dt_table
        !
        !>@param new_dt_table
        !> new detector table
        !--------------------------------------------------------------
        subroutine fill_new_detector_table(
     $     this,
     $     start_i,
     $     new_icoords,
     $     new_rcoords)

          implicit none

          class(bf_detector_icr_list)   , intent(in)    :: this
          integer                       , intent(in)    :: start_i
          integer(ikind), dimension(:,:), intent(inout) :: new_icoords
          real(rkind)   , dimension(:,:), intent(inout) :: new_rcoords

          integer :: k

          
          do k=1, this%nb_detectors

             new_icoords(:,start_i+k-1) = this%icoords(:,k)
             new_rcoords(:,start_i+k-1) = this%rcoords(:,k)

          end do

        end subroutine fill_new_detector_table


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate the content of the object
        !
        !> @date
        !> 24_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_icr_list object encapsulating the data structure
        !> for temporary storing the position of the new increasing
        !> detectors
        !--------------------------------------------------------------
        subroutine destroy(this)
        
          implicit none

          class(bf_detector_icr_list), intent(inout) :: this

          if(allocated(this%icoords)) then
             deallocate(this%icoords)
             deallocate(this%rcoords)
          end if

        end subroutine destroy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the coordinates of the detectors saved in the object
        !> as points on a matrix
        !
        !> @date
        !> 24_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_detector_icr_list object encapsulating the data structure
        !> for temporary storing the position of the new increasing
        !> detectors
        !
        !>@param matrix
        !> array where the general coordinates of the detectors are saved
        !--------------------------------------------------------------
        subroutine print_on_matrix(this, matrix)

          implicit none

          class(bf_detector_icr_list), intent(in)  :: this
          real(rkind), dimension(:,:), intent(out) :: matrix


          real(rkind)                  :: color_detector_list          
          integer(ikind), dimension(2) :: icoords
          integer                      :: k

          if(rkind.eq.8) then
             color_detector_list = 0.2d0
          else
             color_detector_list = 0.2
          end if
          
          do k=1, this%nb_detectors
             icoords = this%icoords(:,k)
             matrix(icoords(1),icoords(2)) = color_detector_list
          end do

        end subroutine print_on_matrix

      end module bf_detector_icr_list_class
      

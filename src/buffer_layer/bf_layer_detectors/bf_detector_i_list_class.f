      module bf_detector_i_list_class
      
        use bf_interface_class, only : bf_interface
        use dbf_element_class , only : dbf_element
        use dbf_list_class    , only : dbf_list
        use parameters_kind   , only : ikind, rkind
      
        implicit none

        private
        public :: bf_detector_i_list


        type :: bf_detector_i_list

          integer, private :: mainlayer_id
          integer, private :: nb_detectors

          integer(ikind), dimension(:,:), allocatable, private :: detectors_list

          type(dbf_list), pointer, private :: detectors_extra_list
          type(dbf_list), pointer, private :: detectors_neighbor1_list
          type(dbf_list), pointer, private :: detectors_neighbor2_list

          contains

          procedure, pass :: ini
          procedure, pass :: add_new_detector

          procedure, pass, private :: add_new_detector_to_mainlayer
          procedure, pass, private :: add_new_detector_to_neighbor1
          procedure, pass, private :: add_new_detector_to_neighbor2
          procedure, pass, private :: add_detector_to_mainlayer
          procedure, pass, private :: add_detector_to_neighbor1
          procedure, pass, private :: add_detector_to_neighbor2

          procedure, nopass, private :: get_inter_detector_param
          procedure, nopass, private :: get_inter_detector_coords

          procedure, pass :: print_on_matrix

          procedure, pass :: destroy

        end type bf_detector_i_list

        contains


        !< initialize the object with the main layer id and the
        !> size of the detector list
        subroutine ini(this, mainlayer_id, size_detectors_list)

          implicit none

          class(bf_detector_i_list), intent(inout) :: this
          integer                  , intent(in)    :: mainlayer_id
          integer(ikind)           , intent(in)    :: size_detectors_list

          this%mainlayer_id = mainlayer_id
          this%nb_detectors = 0
          allocate(this%detectors_list(2,size_detectors_list))
          nullify(this%detectors_extra_list)
          nullify(this%detectors_neighbor1_list)
          nullify(this%detectors_neighbor2_list)

        end subroutine ini


        !< add the coordinates from a new detector to the detectors lists
        !> handled by the main layer: depending on its coordinates, it will
        !> be saved with the detectors of the main layer or with the detectors
        !> of a neighboring main layer.
        !> the distance between the current detector and the previous detector
        !> added to the list is computed to know whether additional detectors
        !> are required
        subroutine add_new_detector(this, interface_used, coords)

          implicit none

          class(bf_detector_i_list)   , intent(inout) :: this
          class(bf_interface)         , intent(in)    :: interface_used
          integer(ikind), dimension(2), intent(in)    :: coords

          integer :: mainlayer_id
          integer :: neighbor_id


          !1) check to which list does the new detector belongs to
          !   using its coordinates
          mainlayer_id = interface_used%get_mainlayer_id(coords)

          !2) check whether the new detector belongs to the main layer
          !   of the detectors or to the neighboring main layers
          !   and add the detectors (and potential intermediates) to
          !   the list it belongs to
          if(mainlayer_id.eq.this%mainlayer_id) then
             call add_new_detector_to_mainlayer(
     $            this, coords)
          else
             neighbor_id = interface_used%get_neighbor_id(
     $            this%mainlayer_id, mainlayer_id)
             
             if(neighbor_id.eq.1) then
                call add_new_detector_to_neighbor1(this, coords)
             else
                call add_new_detector_to_neighbor2(this, coords)
             end if
          end if
          
        end subroutine add_new_detector


        !< add the new detectors and its intermediates to the mainlayer
        !> lists
        subroutine add_new_detector_to_mainlayer(this, coords)

          implicit none

          class(bf_detector_i_list), intent(inout) :: this
          integer(ikind), dimension(2), intent(in) :: coords

          type(dbf_element), pointer :: current_element
          integer(ikind), dimension(2) :: prev_coords
          integer(ikind), dimension(2) :: inter_coords
          real(rkind) :: x_change, y_change
          integer :: inter_nb, k


          !if other detectors were saved in the list before,
          !the new detector added should not be too far away
          !from teh previous one to retain the closed path
          !figure
          if(this%nb_detectors.gt.0) then
             if(this%nb_detectors.le.size(this%detectors_list,2)) then
                prev_coords = this%detectors_list(:,this%nb_detectors)
             else
                current_element => this%detectors_extra_list%get_tail()
                prev_coords = current_element%get_coords()
             end if

             !add intermediate detectors between the previous
             !one and the new one to retain a continuous path
             call get_inter_detector_param(
     $            prev_coords, coords,
     $            x_change, y_change, inter_nb)

             do k=1, inter_nb
                inter_coords = get_inter_detector_coords(
     $               prev_coords,
     $               x_change, y_change, k)
                call add_detector_to_mainlayer(this, inter_coords)
             end do

          end if

          !add the new detector
          call add_detector_to_mainlayer(this, coords)

        end subroutine add_new_detector_to_mainlayer

      
        !< add the new detectors and its intermediates to the
        !> neighbor1 list
        subroutine add_new_detector_to_neighbor1(this, coords)

          implicit none

          class(bf_detector_i_list), intent(inout) :: this
          integer(ikind), dimension(2), intent(in) :: coords

          type(dbf_element), pointer :: current_element
          integer(ikind), dimension(2) :: prev_coords
          integer(ikind), dimension(2) :: inter_coords
          real(rkind) :: x_change, y_change
          integer :: inter_nb, k


          !if other detectors were saved in the list before,
          !the new detector added should not be too far away
          !from the previous one to retain the closed path
          !figure
          if(.not.associated(this%detectors_neighbor1_list)) then
             allocate(this%detectors_neighbor1_list)
             call this%detectors_neighbor1_list%ini()
          end if

          if(this%detectors_neighbor1_list%get_nb_elements().gt.0) then

             current_element => this%detectors_neighbor1_list%get_tail()
             prev_coords = current_element%get_coords()

             !add intermediate detectors between the previous
             !one and the new one to retain a continuous path
             call get_inter_detector_param(
     $            prev_coords, coords,
     $            x_change, y_change, inter_nb)

             do k=1, inter_nb
                inter_coords = get_inter_detector_coords(
     $               prev_coords,
     $               x_change, y_change, k)
                call add_detector_to_neighbor1(this, inter_coords)
             end do

          end if

          !add the new detector
          call add_detector_to_neighbor1(this, coords)

        end subroutine add_new_detector_to_neighbor1


        !< add the new detectors and its intermediates to the
        !> neighbor1 list
        subroutine add_new_detector_to_neighbor2(this, coords)

          implicit none

          class(bf_detector_i_list), intent(inout) :: this
          integer(ikind), dimension(2), intent(in) :: coords

          type(dbf_element), pointer :: current_element
          integer(ikind), dimension(2) :: prev_coords
          integer(ikind), dimension(2) :: inter_coords
          real(rkind) :: x_change, y_change
          integer :: inter_nb, k


          !if other detectors were saved in the list before,
          !the new detector added should not be too far away
          !from the previous one to retain the closed path
          !figure
          if(.not.associated(this%detectors_neighbor2_list)) then
             allocate(this%detectors_neighbor2_list)
             call this%detectors_neighbor2_list%ini()
          end if

          if(this%detectors_neighbor2_list%get_nb_elements().gt.0) then

             current_element => this%detectors_neighbor2_list%get_tail()
             prev_coords = current_element%get_coords()

             !add intermediate detectors between the previous
             !one and the new one to retain a continuous path
             call get_inter_detector_param(
     $            prev_coords, coords,
     $            x_change, y_change, inter_nb)

             do k=1, inter_nb
                inter_coords = get_inter_detector_coords(
     $               prev_coords,
     $               x_change, y_change, k)
                call add_detector_to_neighbor2(this, inter_coords)
             end do

          end if

          !add the new detector
          call add_detector_to_neighbor2(this, coords)

        end subroutine add_new_detector_to_neighbor2


        !< get the parameters constraining the addition
        !> of intermediate detectors between the previous
        !> detectors and the new detector to be added
        subroutine get_inter_detector_param(
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

        end subroutine get_inter_detector_param


        !> from the parameters constraining the addition of
        !> intermediate detectors, give the coordinate of 
        !> the intermediate detector identified by the index k
        !> varying between 1 and total number of detectors to
        !> be added
        function get_inter_detector_coords(
     $     prev_coords,
     $     x_change, y_change, k)
     $     result(inter_coords)

          implicit none

          integer(ikind), dimension(2), intent(in) :: prev_coords
          real(rkind)                 , intent(in) :: x_change
          real(rkind)                 , intent(in) :: y_change
          integer                     , intent(in) :: k
          integer(ikind), dimension(2)             :: inter_coords

          
          inter_coords(1) = prev_coords(1) + nint(x_change*k)
          inter_coords(2) = prev_coords(2) + nint(y_change*k)

        end function get_inter_detector_coords


        !< add detector coordinates to the list saving the detectors
        !> from this main layer
        subroutine add_detector_to_mainlayer(this, coords)

          implicit none

          class(bf_detector_i_list), intent(inout) :: this
          integer(ikind), dimension(2), intent(in) :: coords

          this%nb_detectors = this%nb_detectors+1

          if(this%nb_detectors.le.size(this%detectors_list,2)) then
             this%detectors_list(1,this%nb_detectors) = coords(1)
             this%detectors_list(2,this%nb_detectors) = coords(2)
          else
             if(.not.associated(this%detectors_extra_list)) then
                allocate(this%detectors_extra_list)
                call this%detectors_extra_list%ini()
             end if
             call this%detectors_extra_list%add_to_list(coords)
          end if

        end subroutine add_detector_to_mainlayer

      
        !< add detector coordinates to the list saving the detectors
        !> from the neighbor1 main layer
        subroutine add_detector_to_neighbor1(this, coords)

          implicit none

          class(bf_detector_i_list), intent(inout) :: this
          integer(ikind), dimension(2), intent(in) :: coords

          if(.not.associated(this%detectors_neighbor1_list)) then
             allocate(this%detectors_neighbor1_list)
             call this%detectors_neighbor1_list%ini()
          end if
          call this%detectors_neighbor1_list%add_to_list(coords)

        end subroutine add_detector_to_neighbor1


        !< add detector coordinates to the list saving the detectors
        !> from the neighbor2 main layer
        subroutine add_detector_to_neighbor2(this, coords)

          implicit none

          class(bf_detector_i_list), intent(inout) :: this
          integer(ikind), dimension(2), intent(in) :: coords

          if(.not.associated(this%detectors_neighbor2_list)) then
             allocate(this%detectors_neighbor2_list)
             call this%detectors_neighbor2_list%ini()
          end if
          call this%detectors_neighbor2_list%add_to_list(coords)

        end subroutine add_detector_to_neighbor2


        !< destroy the object
        subroutine destroy(this)
        
          implicit none

          class(bf_detector_i_list), intent(inout) :: this

          if(allocated(this%detectors_list)) then
             deallocate(this%detectors_list)
          end if
          if(associated(this%detectors_extra_list)) then
             call this%detectors_extra_list%destroy()
             deallocate(this%detectors_extra_list)
             nullify(this%detectors_extra_list)
          end if
          if(associated(this%detectors_neighbor1_list)) then
             call this%detectors_neighbor1_list%destroy()
             deallocate(this%detectors_neighbor1_list)
             nullify(this%detectors_neighbor1_list)
          end if
          if(associated(this%detectors_neighbor2_list)) then
             call this%detectors_neighbor2_list%destroy()
             deallocate(this%detectors_neighbor2_list)
             nullify(this%detectors_neighbor2_list)
          end if

        end subroutine destroy


        !> print the coordinates of the detectors saved in the object
        !> as points on a matrix
        subroutine print_on_matrix(this, matrix)

          implicit none

          class(bf_detector_i_list)  , intent(in)  :: this
          real(rkind), dimension(:,:), intent(out) :: matrix


          real(rkind) :: color_detector_list
          real(rkind) :: color_detector_extra_list
          real(rkind) :: color_detector_neighbor1_list
          real(rkind) :: color_detector_neighbor2_list
          
          type(dbf_element), pointer   :: current_element
          integer(ikind), dimension(2) :: coords
          integer                      :: k

          color_detector_list           = 0.2d0
          color_detector_extra_list     = 0.3d0
          color_detector_neighbor1_list = 0.4d0
          color_detector_neighbor2_list = 0.5d0

          
          !detector list
          do k=1, size(this%detectors_list,2)
             coords = this%detectors_list(:,k)
             matrix(coords(1), coords(2)) = color_detector_list
          end do

          
          !extra detector list
          if(associated(this%detectors_extra_list)) then
             current_element => this%detectors_extra_list%get_head()
             do k=1, this%detectors_extra_list%get_nb_elements()
                coords = current_element%get_coords()
                matrix(coords(1), coords(2)) = color_detector_extra_list
                current_element => current_element%get_next()
             end do
          end if


          !neighbor1 list
          if(associated(this%detectors_neighbor1_list)) then
             current_element => this%detectors_neighbor1_list%get_head()
             do k=1, this%detectors_neighbor1_list%get_nb_elements()
                coords = current_element%get_coords()
                matrix(coords(1), coords(2)) = color_detector_neighbor1_list
                current_element => current_element%get_next()
             end do
          end if


          !neighbor2 list
          if(associated(this%detectors_neighbor2_list)) then
             current_element => this%detectors_neighbor2_list%get_head()
             do k=1, this%detectors_neighbor2_list%get_nb_elements()
                coords = current_element%get_coords()
                matrix(coords(1), coords(2)) = color_detector_neighbor2_list
                current_element => current_element%get_next()
             end do
          end if

        end subroutine print_on_matrix

      end module bf_detector_i_list_class
      

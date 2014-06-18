      module bf_detector_icr_list_class
      
        use bf_detector_module, only : get_inter_dct_param,
     $                                 get_inter_dct_coords
        use bf_interface_class, only : bf_interface
        use dbf_element_class , only : dbf_element
        use dbf_list_class    , only : dbf_list
        use parameters_kind   , only : ikind, rkind
      
        implicit none

        private
        public :: bf_detector_icr_list


        type :: bf_detector_icr_list

          integer, private :: mainlayer_id
          integer, private :: nb_detectors

          integer(ikind), dimension(:,:), allocatable, private :: detectors_list

          type(dbf_list), pointer, private :: detectors_extra_list

          contains

          procedure, pass :: ini
          procedure, pass :: add_new_detector => add_new_detector_to_mainlayer

          procedure, pass, private :: add_new_detector_to_mainlayer
          procedure, pass, private :: add_detector_to_mainlayer

          procedure, nopass :: get_inter_detector_param
          procedure, nopass :: get_inter_detector_coords

          procedure,   pass :: get_nb_detectors
          procedure,   pass :: get_prev
          procedure,   pass :: get_next
          procedure,   pass :: get_head
          procedure,   pass :: get_tail
          procedure,   pass :: fill_new_detector_table

          procedure, pass :: print_on_matrix

          procedure, pass :: destroy

        end type bf_detector_icr_list

        contains


        !< initialize the object with the main layer id and the
        !> size of the detector list
        subroutine ini(this, mainlayer_id, size_detectors_list)

          implicit none

          class(bf_detector_icr_list), intent(inout) :: this
          integer                  , intent(in)    :: mainlayer_id
          integer(ikind)           , intent(in)    :: size_detectors_list

          this%mainlayer_id = mainlayer_id
          this%nb_detectors = 0
          allocate(this%detectors_list(2,size_detectors_list))
          nullify(this%detectors_extra_list)
c$$$          nullify(this%detectors_neighbor1_list)
c$$$          nullify(this%detectors_neighbor2_list)

        end subroutine ini


        !< add the new detectors and its intermediates to the mainlayer
        !> lists
        subroutine add_new_detector_to_mainlayer(this, coords)

          implicit none

          class(bf_detector_icr_list), intent(inout) :: this
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

             if((prev_coords(1).ne.coords(1)).or.
     $          (prev_coords(2).ne.coords(2))) then

                !add intermediate detectors between the previous
                !one and the new one to retain a continuous path
                call this%get_inter_detector_param(
     $               prev_coords, coords,
     $               x_change, y_change, inter_nb)
                
                do k=1, inter_nb
                   inter_coords = this%get_inter_detector_coords(
     $                  prev_coords,
     $                  x_change, y_change, k)
                   call add_detector_to_mainlayer(this, inter_coords)
                end do

                call add_detector_to_mainlayer(this, coords)

             end if

          else

             !add the new detector
             call add_detector_to_mainlayer(this, coords)

          end if

        end subroutine add_new_detector_to_mainlayer


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

       
          call get_inter_dct_param(
     $         prev_coords, next_coords,
     $         x_change, y_change,
     $         inter_nb)

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


          inter_coords = get_inter_dct_coords(
     $         prev_coords,
     $         x_change, y_change, k)         
          
        end function get_inter_detector_coords


        !< add detector coordinates to the list saving the detectors
        !> from this main layer
        subroutine add_detector_to_mainlayer(this, coords)

          implicit none

          class(bf_detector_icr_list), intent(inout) :: this
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

      
        !< get nb_detectors
        function get_nb_detectors(this)

          implicit none

          class(bf_detector_icr_list), intent(in) :: this
          integer                               :: get_nb_detectors

          get_nb_detectors = this%nb_detectors

        end function get_nb_detectors


        !< get the previous detector
        function get_prev(this, nb_dt, prev_ele_ptr) result(coords)
        
          implicit none

          class(bf_detector_icr_list) , intent(in)    :: this
          integer                   , intent(inout) :: nb_dt
          type(dbf_element), pointer, intent(inout) :: prev_ele_ptr
          integer(ikind), dimension(2)              :: coords


          nb_dt = nb_dt-1

          if(nb_dt.gt.size(this%detectors_list,2)) then

             prev_ele_ptr => prev_ele_ptr%get_prev()
             coords = prev_ele_ptr%get_coords()

          else             
             coords = this%detectors_list(:,nb_dt)
          end if             

        end function get_prev


        !< get the next detector in the list
        function get_next(this, nb_dt, prev_ele_ptr) result(coords)
        
          implicit none

          class(bf_detector_icr_list) , intent(in)    :: this
          integer                   , intent(inout) :: nb_dt
          type(dbf_element), pointer, intent(inout) :: prev_ele_ptr
          integer(ikind), dimension(2)              :: coords


          nb_dt = nb_dt+1

          if(nb_dt.gt.size(this%detectors_list,2)) then

             if(nb_dt.eq.(size(this%detectors_list,2)+1)) then
                prev_ele_ptr => this%detectors_extra_list%get_head()
                coords = prev_ele_ptr%get_coords()
                
             else
                prev_ele_ptr => prev_ele_ptr%get_next()
                coords = prev_ele_ptr%get_coords()

             end if

          else             
             coords = this%detectors_list(:,nb_dt)
          end if             

        end function get_next


        !< get the first element
        function get_head(this) result(coords)

          implicit none

          class(bf_detector_icr_list) , intent(in)  :: this
          integer(ikind), dimension(2)            :: coords

          if(allocated(this%detectors_list)) then
             coords = this%detectors_list(:,1)
          else
             print '(''bf_detector_icr_list'')'
             stop 'get_head: detectors_list not allocated'
          end if

        end function get_head


        !< get the last element
        function get_tail(this, ele_ptr_i) result(coords)

          implicit none

          class(bf_detector_icr_list)           , intent(in)  :: this
          type(dbf_element), pointer, optional, intent(out) :: ele_ptr_i
          integer(ikind), dimension(2)            :: coords

          type(dbf_element), pointer :: ele_ptr

          if(associated(this%detectors_extra_list)) then
             if(this%detectors_extra_list%get_nb_elements().gt.0) then
                ele_ptr => this%detectors_extra_list%get_tail()
                coords = ele_ptr%get_coords()
             else
                if(allocated(this%detectors_list)) then
                   if(this%nb_detectors.gt.0) then
                      coords = this%detectors_list(:,this%nb_detectors)
                   else
                      print '(''bf_detector_icr_list'')'
                      stop 'get_head: no element'
                   end if
                else
                   print '(''bf_detector_icr_list'')'
                   stop 'get_head: detectors_list not allocated'
                end if
             end if
          else
             if(allocated(this%detectors_list)) then
                if(this%nb_detectors.gt.0) then
                   coords = this%detectors_list(:,this%nb_detectors)
                else
                   print '(''bf_detector_icr_list'')'
                   stop 'get_head: no element'
                end if
             else
                print '(''bf_detector_icr_list'')'
                stop 'get_head: detectors_list not allocated'
             end if
          end if

          if(present(ele_ptr_i)) then
             ele_ptr_i => ele_ptr
          end if

        end function get_tail




        !< fill the new detector table with the new detectors saved
        !> in the bf_detector-i_lit object in detector_list and
        !> detector_extra_list attributes
        subroutine fill_new_detector_table(this, s_index, new_dt_table)

          implicit none

          class(bf_detector_icr_list)     , intent(in) :: this
          integer                       , intent(in) :: s_index
          integer(ikind), dimension(:,:), intent(out):: new_dt_table

          integer :: k
          type(dbf_element), pointer   :: current_ele

          if(allocated(this%detectors_list)) then
             do k=1, size(this%detectors_list,2)
                new_dt_table(:,s_index-1+k) =
     $               this%detectors_list(:,k)
             end do
          end if

          if(associated(this%detectors_extra_list)) then

             current_ele => this%detectors_extra_list%get_head()

             do k=1, this%detectors_extra_list%get_nb_elements()

                new_dt_table(:,s_index-1+size(this%detectors_list,2)+k) =
     $               current_ele%get_coords()

                current_ele => current_ele%get_next()

             end do
          end if

        end subroutine fill_new_detector_table


        !< destroy the object
        subroutine destroy(this)
        
          implicit none

          class(bf_detector_icr_list), intent(inout) :: this

          if(allocated(this%detectors_list)) then
             deallocate(this%detectors_list)
          end if
          if(associated(this%detectors_extra_list)) then
             call this%detectors_extra_list%destroy()
             deallocate(this%detectors_extra_list)
             nullify(this%detectors_extra_list)
          end if

        end subroutine destroy


        !> print the coordinates of the detectors saved in the object
        !> as points on a matrix
        subroutine print_on_matrix(this, matrix)

          implicit none

          class(bf_detector_icr_list)  , intent(in)  :: this
          real(rkind), dimension(:,:), intent(out) :: matrix


          real(rkind) :: color_detector_list
          real(rkind) :: color_detector_extra_list
          
          type(dbf_element), pointer   :: current_element
          integer(ikind), dimension(2) :: coords
          integer                      :: k

          color_detector_list           = 0.2d0
          color_detector_extra_list     = 0.3d0
          
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

        end subroutine print_on_matrix

      end module bf_detector_icr_list_class
      

      !> @file
      !> module implementing the bf_mainlayer object
      !> by encapsulating a double chained list of bf_sublayer
      !> elements
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the bf_mainlayer object
      !> by encapsulating a double chained list of bf_sublayer
      !> elements
      !
      !> @date
      ! 11_04_2014 - initial version      - J.L. Desmarais
      ! 26_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_mainlayer_class
      
        use bf_interior_bc_sections_module, only :
     $     ini_interior_bc_sections,
     $     determine_interior_bc_sections,
     $     close_last_bc_section,
     $     set_full_interior_bc_section,
     $     minimize_interior_bc_section

        use bc_operators_class, only :
     $     bc_operators

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use bf_sublayer_class, only :
     $       bf_sublayer

        use interface_integration_step, only :
     $       timeInt_step_nopt

        use parameters_constant, only :
     $       N,S,E,W,
     $       x_direction,
     $       y_direction,
     $       min_border,
     $       max_border

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       bc_size,
     $       debug

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use td_operators_class, only :
     $       td_operators

        implicit none

        private
        public :: bf_mainlayer
        
        
        !> @class bf_mainlayer
        !> class storing bf_sublayers corresponding
        !> to the same cardinal point (N,S,E,W,NE,NW,SE,SW)
        !
        !> @param mainlayer_id
        !> cardinal coordinate identifying the position of
        !> the bf_sublayer stored in the main layer
        !
        !> @param nb_sublayers
        !> number of sublayers stored in the main layer
        !
        !> @param head_sublayer
        !> pointer to the first bf_sublayer of the list
        !
        !> @param tail_sublayer
        !> pointer to the last bf_sublayer of the list
        !
        !> @param ini
        !> initialize the buffer mainlayer by initializing the
        !> number of sublayers to 0 and nullifying the head
        !> and tail attributes
        !
        !> @param get_mainlayer_id
        !> get the mainlayer_id attribute
        !
        !> @param get_nb_sublayers
        !> number of sublayers stored in the main layer
        !
        !> @param get_head_sublayer
        !> get the head_sublayer attribute
        !
        !> @param get_tail_sublayer
        !> get the tail_sublayer attribute
        !
        !> @param add_sublayer
        !> allocate space for a new buffer sublayer in the double
        !> chained list and organize the bf_mainlayer using the
        !> alignment of the buffer layers. A pointer to the newly
        !> added buffer sublayer is returned
        !
        !> @param merge_sublayers
        !> combine two sublayers of the main layer
        !
        !> @param remove_sublayer
        !> remove a sublayer from the doubled chained list
        !
        !> @param determine_interior_bc_layers
        !> determine the interior_bc_sections, i.e. the boundary
        !> grid points of the interior domain that are computed
        !> by the interior domain and not exchanged with the
        !> buffer layers
        !
        !> @param sync_nodes_with_interior
        !> synchronise the nodes between the interior domain and
        !> the buffer layers constituing the buffer main layer
        !
        !> @param print_binary
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on seperate binary output files
        !
        !> @param print_netcdf
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on seperate netcdf files
        !
        !> @param allocate_before_timeInt
        !> allocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @param deallocate_after_timeInt
        !> deallocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @param compute_time_dev
        !> compute the time derivatives
        !
        !> @param compute_integration_step
        !> compute the integration step
        !---------------------------------------------------------------
        type :: bf_mainlayer

          integer, private :: mainlayer_id
          integer, private :: nb_sublayers

          type(bf_sublayer), pointer, private :: head_sublayer
          type(bf_sublayer), pointer, private :: tail_sublayer

          contains

          procedure, pass :: ini

          !get attributes
          procedure, pass :: get_mainlayer_id
          procedure, pass :: get_nb_sublayers
          procedure, pass :: get_head_sublayer
          procedure, pass :: get_tail_sublayer

          !add/merge/remove a sublayer
          procedure, pass :: add_sublayer
          procedure, pass :: merge_sublayers
          procedure, pass :: remove_sublayer

          !interactions with the interior domain
          procedure, pass :: determine_interior_bc_layers
          procedure, pass :: sync_nodes_with_interior

          !i/o management
          procedure, pass :: print_binary
          procedure, pass :: print_netcdf

          !time integration
          procedure, pass :: allocate_before_timeInt
          procedure, pass :: deallocate_after_timeInt
          procedure, pass :: compute_time_dev
          procedure, pass :: compute_integration_step

        end type bf_mainlayer


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the buffer mainlayer by initializing the
        !> number of sublayers to 0 and nullifying the head
        !> and tail attributes
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_mainlayer object encapsulating the double chained
        !> list of sublayers, pointers to the head and tail elements
        !> of the list and the total number of elements in the list
        !
        !> @param mainlayer_id
        !> cardinal coordinate identifying the position of the buffer
        !> sublayers stored in the mainlayer
        !--------------------------------------------------------------
        subroutine ini(this, mainlayer_id)

          implicit none

          class(bf_mainlayer), intent(inout) :: this
          integer            , intent(in)    :: mainlayer_id

          this%mainlayer_id = mainlayer_id
          this%nb_sublayers = 0
          nullify(this%head_sublayer)
          nullify(this%tail_sublayer)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the mainlayer_id attribute
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_mainlayer object encapsulating the double chained
        !> list of sublayers, pointers to the head and tail elements
        !> of the list and the total number of elements in the list
        !
        !>@return get_mainlayer_id
        !> cardinal coordinate identifying the position of the buffer
        !> sublayers stored in the mainlayer
        !--------------------------------------------------------------
        function get_mainlayer_id(this)

          implicit none

          class(bf_mainlayer), intent(in) :: this
          integer                         :: get_mainlayer_id

          get_mainlayer_id = this%mainlayer_id

        end function get_mainlayer_id


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the nb_sublayers attribute
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_mainlayer object encapsulating the double chained
        !> list of sublayers, pointers to the head and tail elements
        !> of the list and the total number of elements in the list
        !
        !>@return get_nb_sublayers
        !> number of sublayers in the mainlayer chained list
        !--------------------------------------------------------------
        function get_nb_sublayers(this)

          implicit none

          class(bf_mainlayer), intent(in) :: this
          integer                         :: get_nb_sublayers

          get_nb_sublayers = this%nb_sublayers

        end function get_nb_sublayers


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the head_sublayer attribute
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_mainlayer object encapsulating the double chained
        !> list of sublayers, pointers to the head and tail elements
        !> of the list and the total number of elements in the list
        !
        !>@return get_head_sublayer
        !> pointer to the first sublayer in the chained list
        !--------------------------------------------------------------
        function get_head_sublayer(this)

          implicit none

          class(bf_mainlayer), intent(in) :: this
          type(bf_sublayer)  , pointer    :: get_head_sublayer

          if(associated(this%head_sublayer)) then
             get_head_sublayer => this%head_sublayer
          else
             nullify(get_head_sublayer)
          end if

        end function get_head_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the tail_sublayer attribute
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_mainlayer object encapsulating the double chained
        !> list of sublayers, pointers to the head and tail elements
        !> of the list and the total number of elements in the list
        !
        !>@return get_head_sublayer
        !> pointer to the last sublayer in the chained list
        !--------------------------------------------------------------
        function get_tail_sublayer(this)

          implicit none

          class(bf_mainlayer), intent(in) :: this
          type(bf_sublayer)  , pointer    :: get_tail_sublayer

          if(associated(this%tail_sublayer)) then
             get_tail_sublayer => this%tail_sublayer
          else
             nullify(get_tail_sublayer)
          end if

        end function get_tail_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate space for a new buffer sublayer in the double
        !> chained list and organize the bf_mainlayer using the
        !> alignment of the buffer layers. A pointer to the newly
        !> added buffer sublayer is returned
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_mainlayer object encapsulating the double chained
        !> list of sublayers, pointers to the head and tail elements
        !> of the list and the total number of elements in the list
        !
        !> @param nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !> @param alignment
        !> table(2,2) of integers identifying the position of the buffer
        !> layer compared to the interior nodes
        !
        !> @return added_sublayer_ptr
        !> pointer to the newly added bf_sublayer
        !--------------------------------------------------------------
        function add_sublayer(
     $     this,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     alignment)
     $     result(added_sublayer_ptr)

          implicit none
          
          class(bf_mainlayer)                , intent(inout) :: this
          real(rkind)   , dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2)     , intent(in)    :: alignment
          type(bf_sublayer), pointer                         :: added_sublayer_ptr


          type(bf_sublayer), pointer :: current_sublayer
          type(bf_sublayer), pointer :: current_sublayer_prev

          integer :: direction_tested
          integer :: i

          !allocate space for the sublayer and initialize it
          allocate(added_sublayer_ptr)
          call added_sublayer_ptr%ini(this%mainlayer_id)
          call added_sublayer_ptr%allocate_bf_layer(
     $         interior_x_map,
     $         interior_y_map, 
     $         interior_nodes,
     $         alignment)


          !if there is already an element in the list of sublayers
          !of the mainlayer, go through the list of sublayers and
          !insert the newly sublayer where its alignment fits
          select case(this%nb_sublayers)
            case(0)

               this%head_sublayer => added_sublayer_ptr
               this%tail_sublayer => added_sublayer_ptr

            case(1)
               
               select case(this%mainlayer_id)
                 case(N,S)
                    direction_tested = x_direction
                 case(E,W)
                    direction_tested = y_direction
                 case default
                    call error_mainlayer_id(
     $                   'bf_mainlayer_class.f',
     $                   'add_sublayer',
     $                   this%mainlayer_id)
               end select

               if(this%head_sublayer%get_alignment(direction_tested,min_border).lt.
     $            alignment(direction_tested,min_border)) then

                  call this%head_sublayer%set_next(added_sublayer_ptr)
                  call added_sublayer_ptr%set_prev(this%head_sublayer)
                  this%tail_sublayer => added_sublayer_ptr

               else
                
                  call added_sublayer_ptr%set_next(this%head_sublayer)
                  call this%head_sublayer%set_prev(added_sublayer_ptr)
                  call this%head_sublayer%nullify_next()
                  this%head_sublayer => added_sublayer_ptr

               end if
                
            case default
             
               select case(this%mainlayer_id)
                 case(N,S)
                    direction_tested = x_direction
                 case(E,W)
                    direction_tested = y_direction
                 case default
                    call error_mainlayer_id(
     $                   'bf_mainlayer_class.f',
     $                   'add_sublayer',
     $                   this%mainlayer_id)
               end select

               current_sublayer => this%head_sublayer

               do i=1, this%nb_sublayers-1
                  
                  if(current_sublayer%get_alignment(direction_tested, min_border).gt.
     $               alignment(direction_tested, min_border)) then
                     exit
                  end if

                  current_sublayer => current_sublayer%get_next()

               end do


               if(i.ne.this%nb_sublayers) then
                  if(associated(current_sublayer%get_prev())) then
                     current_sublayer_prev => current_sublayer%get_prev()
                     call insert_sublayer(
     $                    added_sublayer_ptr,
     $                    current_sublayer_prev,
     $                    current_sublayer)
                  else
                     this%head_sublayer => added_sublayer_ptr
                     call added_sublayer_ptr%set_next(current_sublayer)
                     call current_sublayer%set_prev(added_sublayer_ptr)
                  end if

               else
                  if(current_sublayer%get_alignment(direction_tested, min_border).gt.
     $               alignment(direction_tested, min_border)) then

                     current_sublayer_prev => current_sublayer%get_prev()
                     call insert_sublayer(
     $                    added_sublayer_ptr,
     $                    current_sublayer_prev,
     $                    current_sublayer)

                  else
                     call current_sublayer%set_next(added_sublayer_ptr)
                     call added_sublayer_ptr%set_prev(current_sublayer)
                     this%tail_sublayer => added_sublayer_ptr

                  end if
               end if

          end select

          !update the number of sublayers in the mainlayer
          this%nb_sublayers = this%nb_sublayers+1

        end function add_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine inserting a buffer sublayer between two sublayer
        !> elements
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !> @param inserted_sublayer
        !> sublayer element inserted between two existing sublayer elements
        !> of the double chained list
        !
        !> @param prev_sublayer
        !> sublayer element of the list that will be placed before the
        !> inserted element afterwards
        !
        !> @param next_sublayer
        !> sublayer element of the list that will be placed after the
        !> inserted element afterwards
        !--------------------------------------------------------------
        subroutine insert_sublayer(
     $     inserted_sublayer,
     $     prev_sublayer,
     $     next_sublayer)

          implicit none

          type(bf_sublayer), pointer, intent(inout) :: inserted_sublayer
          type(bf_sublayer), pointer, intent(inout) :: prev_sublayer
          type(bf_sublayer), pointer, intent(inout) :: next_sublayer

          call prev_sublayer%set_next(inserted_sublayer)
          call inserted_sublayer%set_prev(prev_sublayer)

          call inserted_sublayer%set_next(next_sublayer)
          call next_sublayer%set_prev(inserted_sublayer)

        end subroutine insert_sublayer
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> combine two sublayers of the main layer
        !
        !> @date
        !> 09_05_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param bf_sublayer1
        !> pointer to the first sublayer to merge
        !
        !> @param bf_sublayer2
        !> pointer to the second sublayer to merge
        !
        !> @param interior_nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !> @param alignment
        !> table identifying the final position of the merged sublayer
        !> compared to the interior domain
        !
        !>@return merged_sublayer
        !> pointer to the bf_sublayer resulting from the merge of the
        !> buffer sublayers
        !--------------------------------------------------------------
        function merge_sublayers(
     $     this,
     $     bf_sublayer1,
     $     bf_sublayer2,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     alignment)
     $     result(merged_sublayer)
        
          implicit none        
        
          class(bf_mainlayer)                     , intent(inout) :: this
          type(bf_sublayer), pointer              , intent(inout) :: bf_sublayer1
          type(bf_sublayer), pointer              , intent(inout) :: bf_sublayer2
          real(rkind)   , dimension(nx)           , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)           , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne)     , intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2), optional, intent(in)    :: alignment
          type(bf_sublayer), pointer                              :: merged_sublayer
        
          
          integer :: direction_tested
          type(bf_sublayer), pointer :: temp_ptr


          !reorganize the position of the elements in the chained
          !list of sublayers
          select case(bf_sublayer1%get_localization())
            case(N,S)
               direction_tested = x_direction
            case(E,W)
               direction_tested = y_direction
            case default
               call error_mainlayer_id(
     $              'bf_mainlayer_class.f',
     $              'merge_sublayers',
     $              bf_sublayer1%get_localization())
          end select

          if(bf_sublayer1%get_alignment(direction_tested,min_border).lt.
     $         bf_sublayer2%get_alignment(direction_tested,min_border)) then
             
             temp_ptr => bf_sublayer2%get_next()

             if(associated(temp_ptr)) then
                call bf_sublayer1%set_next(temp_ptr)
                call temp_ptr%set_prev(bf_sublayer1)

             else
                call bf_sublayer1%nullify_next()
                this%tail_sublayer => bf_sublayer1

             end if
             
          else

             temp_ptr => bf_sublayer2%get_prev()

             if(associated(temp_ptr)) then
                call bf_sublayer1%set_prev(temp_ptr)
                call temp_ptr%set_next(bf_sublayer1)

             else
                call bf_sublayer1%nullify_prev()
                this%head_sublayer => bf_sublayer1
             end if
             
          end if


          !merge the buffer sublayers 1 and 2
          if(present(alignment)) then
             call bf_sublayer1%merge_bf_layer(
     $            bf_sublayer2,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            alignment)
          else
             call bf_sublayer1%merge_bf_layer(
     $            bf_sublayer2,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes)
          end if

        
          !destroy the bf_sublayer2
          call bf_sublayer2%nullify_next()
          call bf_sublayer2%nullify_prev()
          deallocate(bf_sublayer2)


          !update the number of sublayers in the mainlayer
          this%nb_sublayers = this%nb_sublayers-1


          !return the pointer to the merged sublayer
          merged_sublayer => bf_sublayer1

        end function merge_sublayers


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove a sublayer from the doubled chained list
        !
        !> @date
        !> 09_05_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param bf_sublayer_ptr
        !> pointer to the sublayer to be removed
        !--------------------------------------------------------------
        subroutine remove_sublayer(this, sublayer_ptr)

          implicit none

          class(bf_mainlayer)         , intent(inout) :: this
          type(bf_sublayer)  , pointer, intent(inout) :: sublayer_ptr


          !check if the head of the mainlayer should be changed
          if(associated(this%head_sublayer,sublayer_ptr)) then
             if(associated(sublayer_ptr%get_next())) then
                this%head_sublayer => sublayer_ptr%get_next()
             else
                nullify(this%head_sublayer)
             end if
          end if


          !check if the tail of the mainlayer should be changed
          if(associated(this%tail_sublayer, sublayer_ptr)) then
             if(associated(sublayer_ptr%get_prev())) then
                this%tail_sublayer => sublayer_ptr%get_prev()
             else
                nullify(this%tail_sublayer)
             end if
          end if


          !remove the buffer layer
          call sublayer_ptr%remove()
          deallocate(sublayer_ptr)
          nullify(sublayer_ptr)


          !update the number of sublayers in the mainlayer
          this%nb_sublayers = this%nb_sublayers-1

        end subroutine remove_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the extent of the boundary layers computed
        !> by the interior nodes
        !
        !> @date
        !> 28_10_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param interior_bc_sections
        !> extent of the boundary layers computed by the interior
        !> nodes
        !--------------------------------------------------------------
        subroutine determine_interior_bc_layers(this, bc_sections)

          implicit none

          class(bf_mainlayer)                        , intent(in)    :: this
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections

          integer        :: dir
          integer(ikind) :: interior_inf
          integer(ikind) :: interior_sup

          integer :: nb_bc_sections
          logical :: min_initialized
          logical :: max_initialized
          logical :: no_bf_common_with_interior

          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: k

          integer(ikind), dimension(2,2) :: bf_alignments
          integer(ikind), dimension(2)   :: bf_alignment


          !initialize the variables for the determination of the
          !interior boundary layers depending on the cardinal
          !coordinate of the main layer
          select case(this%mainlayer_id)
            case(N,S)
               dir          = 1
               interior_inf = 1
               interior_sup = nx

            case(E,W)
               dir          = 2
               interior_inf = bc_size+1
               interior_sup = ny-bc_size
               
            case default
               call error_mainlayer_id(
     $              'bf_mainlayer_class.f',
     $              'determine_interior_bc_sections',
     $              this%mainlayer_id)
               
          end select


          !initialize the interior_bc_sections
          call ini_interior_bc_sections(
     $         nb_bc_sections,
     $         min_initialized,
     $         max_initialized,
     $         no_bf_common_with_interior)


          !initialize the pointer to the first sublayer
          !investigated
          current_sublayer => this%head_sublayer


          !go through the chained list and update the
          !extents of the interior boundary layers
          if(this%nb_sublayers.gt.0) then

             do k=1, this%nb_sublayers

                !determine the alignment of the sublayer
                bf_alignments = current_sublayer%get_alignment_tab()

                !determine the alignment relevant for the interior
                bf_alignment(1) = bf_alignments(dir,1)
                bf_alignment(2) = bf_alignments(dir,2)

                !update the extent of the interior boundary
                !layers
                call determine_interior_bc_sections(
     $               bf_alignment,
     $               interior_inf,
     $               interior_sup,
     $               nb_bc_sections,
     $               bc_sections,
     $               min_initialized,
     $               max_initialized,
     $               no_bf_common_with_interior)

                !get the next sublayer in the mainlayer
                current_sublayer => current_sublayer%get_next()

             end do

             !finalize the interior_bc_sections
             call close_last_bc_section(
     $            nb_bc_sections,
     $            bc_sections,
     $            interior_sup,
     $            min_initialized,
     $            max_initialized)

          end if

          if(no_bf_common_with_interior) then
             call set_full_interior_bc_section(
     $            nb_bc_sections,
     $            bc_sections,
     $            interior_inf,
     $            interior_sup)
          else

            !minimize the extent of the interior boundary
            !layers
             call minimize_interior_bc_section(
     $            nb_bc_sections,
     $            bc_sections)

          end if

        end subroutine determine_interior_bc_layers


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> exchange grid points between the buffer main layer and
        !> interior domain
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param interior_nodes
        !> grid points from the interior domain
        !--------------------------------------------------------------
        subroutine sync_nodes_with_interior(this, interior_nodes)

          implicit none

          class(bf_mainlayer)             , intent(inout) :: this
          real(rkind), dimension(nx,ny,ne), intent(inout) :: interior_nodes

          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i


          current_sublayer => this%head_sublayer

          do i=1, this%nb_sublayers

             call current_sublayer%sync_nodes_with_interior(
     $            interior_nodes)

             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine sync_nodes_with_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on seperate binary output files
        !
        !> @date
        !> 09_05_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param suffix_nodes
        !> suffix for the name of the output files storing the nodes
        !> of the bf_sublayers
        !
        !> @param suffix_grdid
        !> suffix for the name of the output files storing the grdpts_id
        !> of the bf_sublayers
        !
        !> @param suffix_sizes
        !> suffix for the name of the output files storing the sizes
        !> of the bf_sublayers        
        !--------------------------------------------------------------
        subroutine print_binary(
     $     this,
     $     suffix_x_map,
     $     suffix_y_map,
     $     suffix_nodes,
     $     suffix_grdid,
     $     suffix_sizes)

          implicit none

          class(bf_mainlayer), intent(in) :: this
          character(*)       , intent(in) :: suffix_x_map
          character(*)       , intent(in) :: suffix_y_map
          character(*)       , intent(in) :: suffix_nodes
          character(*)       , intent(in) :: suffix_grdid
          character(*)       , intent(in) :: suffix_sizes


          character(2), dimension(4) :: bf_layer_char
          integer                    :: i
          character(len=14)          :: filename_format
          character(len=30)          :: sizes_filename
          character(len=30)          :: x_map_filename
          character(len=30)          :: y_map_filename
          character(len=30)          :: nodes_filename
          character(len=30)          :: grdid_filename
          type(bf_sublayer), pointer :: current_sublayer

          bf_layer_char = ['N_','S_','E_','W_']


          current_sublayer => this%head_sublayer


          !go through the chained list and write the content
          !of each element
          do i=1, this%nb_sublayers

             !determine the name of the output files
             write(filename_format,
     $            '(''(A2,I1,A1,A'',I2, '')'')')
     $            len(suffix_nodes)

             write(x_map_filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_', suffix_x_map

             write(y_map_filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_', suffix_y_map

             write(nodes_filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_', suffix_nodes


             write(filename_format,
     $            '(''(A2,I1,A1,A'',I2, '')'')')
     $            len(suffix_grdid)
             write(grdid_filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_', suffix_grdid

             write(filename_format,
     $            '(''(A2,I1,A1,A'',I2, '')'')')
     $            len(suffix_sizes)
             write(sizes_filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_', suffix_sizes


             !write the content of the sublayer corresponding
             !to the index i
             call current_sublayer%print_binary(
     $            x_map_filename,
     $            y_map_filename,
     $            nodes_filename,
     $            grdid_filename,
     $            sizes_filename)

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine print_binary


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on a netcdf output file
        !
        !> @date
        !> 11_07_2013 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param timestep_written
        !> integer identifying the timestep written
        !
        !> @param name_var
        !> table with the short name for the governing variables saved
        !> in the netcdf file
        !
        !> @param longname_var
        !> table with the long name for the governing variables saved
        !> in the netcdf file
        !
        !> @param unit_var
        !> table with the units of the governing variables saved
        !> in the netcdf file
        !
        !> @param x_min_interior
        !> x-coordinate of interior grid point next to the left
        !> boundary layer
        !
        !> @param y_min_interior
        !> y-coordinate of interior grid point next to the lower
        !> boundary layer
        !
        !> @param dx
        !> grid size along the x-coordinate
        !
        !> @param dy
        !> grid size along the y-coordinate
        !
        !> @param time
        !> time corresponding to the data for the grdpts and the nodes
        !--------------------------------------------------------------
        subroutine print_netcdf(
     $     this,
     $     timestep_written,
     $     name_var,
     $     longname_var,
     $     unit_var,
     $     time)

          implicit none

          class(bf_mainlayer)        , intent(in) :: this
          integer                    , intent(in) :: timestep_written
          character(*), dimension(ne), intent(in) :: name_var
          character(*), dimension(ne), intent(in) :: longname_var
          character(*), dimension(ne), intent(in) :: unit_var
          real(rkind)                , intent(in) :: time


          character(2), dimension(4) :: bf_layer_char
          integer                    :: size_timestep
          integer                    :: size_bf_order
          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i
          character(len=16)          :: filename_format
          character(len=30)          :: filename


          bf_layer_char = ['N_','S_','E_','W_']


          !determine the size needed to write the timestep
          if (timestep_written.eq.0) then
             size_timestep  = 1
          else
             size_timestep  = floor(
     $            log(float(timestep_written))/log(10.))+1
          end if

          
          !initialize to the pointer to the first sublayer
          !written
          current_sublayer => this%head_sublayer


          !go through the chained list and write the content
          !of each element
          do i=1, this%nb_sublayers

             !determine the size needed to write the bf_order
             size_bf_order = floor(log(float(i))/log(10.))+1

             !determine the name of the netcdf file
             write(filename_format,
     $            '(''(A2,I'',I1,'',A1,I'',I1,'',A3)'')')
     $            size_bf_order, size_timestep

             write(filename, filename_format)
     $            bf_layer_char(this%mainlayer_id), i, '_',
     $            timestep_written, '.nc'

             !write the content of the sublayer corresponding
             !to the index i
             call current_sublayer%print_netcdf(
     $            trim(filename),
     $            name_var,
     $            longname_var,
     $            unit_var,
     $            i,
     $            time)

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine print_netcdf


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate memory space for the intermediate
        !> variables needed to perform the time integration
        !> for each sublayer contained in this main layer
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !--------------------------------------------------------------
        subroutine allocate_before_timeInt(this)

          implicit none

          class(bf_mainlayer), intent(inout) :: this


          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i

          
          !initialize to the pointer to the first sublayer
          !written
          current_sublayer => this%head_sublayer


          !go through the chained list and allocate the 
          !temporary variables for the time integration
          do i=1, this%nb_sublayers


             !allocate the temporary variables for the 
             !time integration
             call current_sublayer%allocate_before_timeInt()

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine allocate_before_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate the memory space for the intermediate
        !> variables needed to perform the time integration
        !> for each sublayer contained in this main layer
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !--------------------------------------------------------------
        subroutine deallocate_after_timeInt(this)

          implicit none

          class(bf_mainlayer), intent(inout) :: this


          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i

          
          !initialize to the pointer to the first sublayer
          !written
          current_sublayer => this%head_sublayer


          !go through the chained list and allocate the 
          !temporary variables for the time integration
          do i=1, this%nb_sublayers


             !allocate the temporary variables for the 
             !time integration
             call current_sublayer%deallocate_after_timeInt()

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine deallocate_after_timeInt



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the sublayers
        !> contained in this main layer
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !--------------------------------------------------------------
        subroutine compute_time_dev(
     $     this,
     $     td_operators_used,
     $     t,s,p_model,bc_used)

          implicit none

          class(bf_mainlayer), intent(inout) :: this
           type(td_operators), intent(in)    :: td_operators_used
          real(rkind)        , intent(in)    :: t
          type(sd_operators) , intent(in)    :: s
          type(pmodel_eq)    , intent(in)    :: p_model
          type(bc_operators) , intent(in)    :: bc_used

          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i

          
          !initialize to the pointer to the first sublayer
          !written
          current_sublayer => this%head_sublayer


          !go through the chained list and allocate the 
          !temporary variables for the time integration
          do i=1, this%nb_sublayers


             !allocate the temporary variables for the 
             !time integration
             call current_sublayer%compute_time_dev(
     $            td_operators_used,
     $            t,s,p_model,bc_used)

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine compute_time_dev



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the sublayers
        !> contained in this main layer
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param dt
        !> integration time step
        !
        !> @param integration_step_nopt
        !> procedure performing the time integration
        !--------------------------------------------------------------
        subroutine compute_integration_step(
     $     this, dt, integration_step_nopt)

          implicit none

          class(bf_mainlayer), intent(inout) :: this
          real(rkind)        , intent(in)    :: dt
          procedure(timeInt_step_nopt) :: integration_step_nopt

          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i

          
          !initialize to the pointer to the first sublayer
          !written
          current_sublayer => this%head_sublayer


          !go through the chained list and allocate the 
          !temporary variables for the time integration
          do i=1, this%nb_sublayers


             !allocate the temporary variables for the 
             !time integration
             call current_sublayer%compute_integration_step(
     $            dt, integration_step_nopt)

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine compute_integration_step        

      end module bf_mainlayer_class

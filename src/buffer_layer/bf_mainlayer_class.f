      !> @file
      !> module implementing the buffer mainlayer object
      !> encapsulating a double chained list of buffer
      !> sublayer elements
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the buffer mainlayer object
      !> encapsulating a double chained list of buffer
      !> sublayer elements
      !
      !> @date
      ! 11_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_mainlayer_class

        use bf_sublayer_class  , only : bf_sublayer
        use parameters_constant, only : N,S,E,W,
     $                                  x_direction, y_direction,
     $                                  min_border, max_border
        use parameters_input   , only : nx, ny, ne, debug
        use parameters_kind    , only : ikind, rkind

        implicit none

        private
        public :: bf_mainlayer
        
        
        !> @class bf_mainlayer
        !> class encapsulating the buffer sublayers corresponding
        !> to the same cardinal point (N,S,E,W,NE,NW,SE,SW)
        !>
        !> @param nb_sublayers
        !> number of sublayers saved in the main layer
        !>
        !> @param head_sublayer
        !> pointer of the head sublayer of the main layer
        !>
        !> @param tail_sublayer
        !> pointer of the tail sublayer of the main layer
        !---------------------------------------------------------------
        type :: bf_mainlayer

          integer, private :: mainlayer_id
          integer, private :: nb_sublayers

          type(bf_sublayer), pointer, private :: head_sublayer
          type(bf_sublayer), pointer, private :: tail_sublayer

          contains

          procedure, pass :: ini

          procedure, pass :: get_mainlayer_id
          procedure, pass :: get_nb_sublayers
          procedure, pass :: get_head_sublayer
          procedure, pass :: get_tail_sublayer

          procedure, pass :: add_sublayer
          procedure, pass :: merge_sublayers

          procedure, pass :: print_binary

          procedure, nopass, private :: insert_sublayer

        end type bf_mainlayer


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the buffer mainlayer
        !> by initializing the number of sublayers to 0 and
        !> nullifying the head and tail pointers
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer object encapsulating the double chained
        !> list of sublayers, pointers to the head and tail elements
        !> of the list and the total number of elements in the list
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


        !< get the localization of the mainlayer
        function get_mainlayer_id(this)

          implicit none

          class(bf_mainlayer), intent(in) :: this
          integer                         :: get_mainlayer_id

          get_mainlayer_id = this%mainlayer_id

        end function get_mainlayer_id


        !< get the number of sublayers in the mainlayer chained list
        function get_nb_sublayers(this)

          implicit none

          class(bf_mainlayer), intent(in) :: this
          integer                         :: get_nb_sublayers

          get_nb_sublayers = this%nb_sublayers

        end function get_nb_sublayers


        !< get the first sublayer in the chained list of the mainlayer
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


        !< get the last sublayer in the chained list of the mainlayer
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
        !> subroutine allocating space for a new buffer sublayer
        !> in the double chained list and positioning the sublayer
        !> in this double chained list according to the mainlayer
        !> location and its alignment to the interior nodes. A pointer
        !> to the newly added buffer sublayer is returned
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_mainlayer object encapsulating the double chained
        !> list of sublayers, pointers to the head and tail elements
        !> of the list and the total number of elements in the list
        !
        !>@param mainlayer_id
        !> integer identifying the main layer localization
        !
        !>@param alignment
        !> table(2,2) of integers identifying the position of the buffer
        !> layer compared to the interior nodes
        !--------------------------------------------------------------
        function add_sublayer(this, nodes, alignment)
     $     result(added_sublayer_ptr)

          implicit none
          
          class(bf_mainlayer)                , intent(inout) :: this
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: nodes
          integer(ikind), dimension(2,2)     , intent(in)    :: alignment
          type(bf_sublayer), pointer                         :: added_sublayer_ptr


          type(bf_sublayer), pointer :: current_sublayer
          type(bf_sublayer), pointer :: current_sublayer_prev

          integer :: direction_tested
          integer :: i

          !allocate space for the sublayer and initialize it
          allocate(added_sublayer_ptr)
          call added_sublayer_ptr%ini(this%mainlayer_id)
          call added_sublayer_ptr%allocate_bf_layer(nodes, alignment)


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
                    print '(''bf_mainlayer_class'')'
                    print '(''merge_sublayers'')'
                    print '(''corner mainlayers cannot host more'')'
                    print '(''than one sublayer'')'
                    print '(''mainlayer_id: '', I2)', this%mainlayer_id
                    stop 'check mainlayer_id'
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
                    print '(''bf_mainlayer_class'')'
                    print '(''merge_sublayers'')'
                    print '(''corner mainlayers cannot host more'')'
                    print '(''than one sublayer'')'
                    print '(''mainlayer_id: '', I2)', this%mainlayer_id
                    stop 'check mainlayer_id'
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
        !>@param inserted_sublayer
        !> sublayer element inserted between two existing sublayer elements
        !> of the double chained list
        !
        !>@param prev_sublayer
        !> sublayer element of the list that will be placed before the
        !> inserted element afterwards
        !
        !>@param next_sublayer
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
        !> subroutine combining two sublayers of the main layer
        !
        !> @date
        !> 09_05_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !>@param sublayer1
        !> pointer to the first sublayer to merge
        !
        !>@param sublayer2
        !> pointer to the second sublayer to merge
        !
        !>@param alignment
        !> table identifying the final position of the merged sublayer
        !> compared to the interior domain
        !--------------------------------------------------------------
        function merge_sublayers(
     $     this,
     $     bf_sublayer1,
     $     bf_sublayer2,
     $     interior_nodes,
     $     alignment)
     $     result(merged_sublayer)
        
          implicit none        
        
          class(bf_mainlayer)                     , intent(inout) :: this
          type(bf_sublayer), pointer              , intent(inout) :: bf_sublayer1
          type(bf_sublayer), pointer              , intent(inout) :: bf_sublayer2
          real(rkind)   , dimension(nx,ny,ne)     , intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2), optional, intent(in)    :: alignment
          type(bf_sublayer), pointer                              :: merged_sublayer
        

          !merge the buffer sublayers 1 and 2
          if(present(alignment)) then
             call bf_sublayer1%merge_sublayers(
     $            bf_sublayer1, bf_sublayer2,
     $            interior_nodes,
     $            alignment)
          else
             call bf_sublayer1%merge_sublayers(
     $            bf_sublayer1, bf_sublayer2,
     $            interior_nodes)
          end if
        
          !update the number of sublayers in the mainlayer
          this%nb_sublayers = this%nb_sublayers-1

          !return the pointer to the merged sublayer
          merged_sublayer => bf_sublayer1

        end function merge_sublayers


        !< print the content of the sublayers constituing the
        !> buffer main layer on seperate binary output files
        subroutine print_binary(
     $     this, suffix_nodes, suffix_grdid, suffix_sizes)

          implicit none

          class(bf_mainlayer), intent(in) :: this
          character(*)       , intent(in) :: suffix_nodes
          character(*)       , intent(in) :: suffix_grdid
          character(*)       , intent(in) :: suffix_sizes


          character(2), dimension(8) :: bf_layer_char
          integer                    :: i
          character(len=14)          :: filename_format
          character(len=19)          :: sizes_filename
          character(len=19)          :: nodes_filename
          character(len=19)          :: grdid_filename
          type(bf_sublayer), pointer :: current_sublayer

          bf_layer_char = ['N_','S_','E_','W_','NE','NW','SE','SW']


          current_sublayer => this%head_sublayer


          !go through the chained list and write the content
          !of each element
          do i=1, this%nb_sublayers

             !determine the name of the output files
             write(filename_format,
     $            '(''(A2,I1,A1,A'',I2, '')'')')
     $            len(suffix_nodes)

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
     $            nodes_filename,
     $            grdid_filename,
     $            sizes_filename)

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine print_binary
      

      end module bf_mainlayer_class

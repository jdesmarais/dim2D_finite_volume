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

        use bf_layer_errors_module, only : error_mainlayer_id
        use bf_sublayer_class     , only : bf_sublayer
        use parameters_constant   , only : N,S,E,W,
     $                                     x_direction, y_direction,
     $                                     min_border, max_border
        use parameters_input      , only : nx, ny, ne, debug
        use parameters_kind       , only : ikind, rkind

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
        !>@param ini
        !> initialize the buffer mainlayer by initializing the
        !> number of sublayers to 0 and nullifying the head
        !> and tail attributes
        !
        !>@param get_mainlayer_id
        !> get the mainlayer_id attribute
        !
        !>@param get_nb_sublayers
        !> number of sublayers stored in the main layer
        !
        !>@param get_head_sublayer
        !> get the head_sublayer attribute
        !
        !>@param get_tail_sublayer
        !> get the tail_sublayer attribute
        !
        !>@param add_sublayer
        !> allocate space for a new buffer sublayer in the double
        !> chained list and organize the bf_mainlayer using the
        !> alignment of the buffer layers. A pointer to the newly
        !> added buffer sublayer is returned
        !
        !>@param merge_sublayers
        !> combine two sublayers of the main layer
        !
        !>@param remove_sublayer
        !> remove a sublayer from the doubled chained list
        !
        !>@param print_binary
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on seperate binary output files
        !
        !>@param print_netcdf
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on seperate netcdf files
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
          procedure, pass :: remove_sublayer

          procedure, pass :: print_binary
          procedure, pass :: print_netcdf

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
        !>@param this
        !> bf_mainlayer object encapsulating the double chained
        !> list of sublayers, pointers to the head and tail elements
        !> of the list and the total number of elements in the list
        !
        !>@param mainlayer_id
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
        !>@param this
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
        !>@param this
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
        !>@param this
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
        !>@param this
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
        !>@param this
        !> bf_mainlayer object encapsulating the double chained
        !> list of sublayers, pointers to the head and tail elements
        !> of the list and the total number of elements in the list
        !
        !>@param nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@param alignment
        !> table(2,2) of integers identifying the position of the buffer
        !> layer compared to the interior nodes
        !
        !> @return added_sublayer_ptr
        !> pointer to the newly added bf_sublayer
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
        !> combine two sublayers of the main layer
        !
        !> @date
        !> 09_05_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !>@param bf_sublayer1
        !> pointer to the first sublayer to merge
        !
        !>@param bf_sublayer2
        !> pointer to the second sublayer to merge
        !
        !>@param interior_nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@param alignment
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
     $            interior_nodes,
     $            alignment)
          else
             call bf_sublayer1%merge_bf_layer(
     $            bf_sublayer2,
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
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !>@param bf_sublayer_ptr
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
        !> print the content of the bf_sublayers constituing the
        !> bf_mainlayer on seperate binary output files
        !
        !> @date
        !> 09_05_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !>@param suffix_nodes
        !> suffix for the name of the output files storing the nodes
        !> of the bf_sublayers
        !
        !>@param suffix_grdid
        !> suffix for the name of the output files storing the grdpts_id
        !> of the bf_sublayers
        !
        !>@param suffix_sizes
        !> suffix for the name of the output files storing the sizes
        !> of the bf_sublayers        
        !--------------------------------------------------------------
        subroutine print_binary(
     $     this, suffix_nodes, suffix_grdid, suffix_sizes)

          implicit none

          class(bf_mainlayer), intent(in) :: this
          character(*)       , intent(in) :: suffix_nodes
          character(*)       , intent(in) :: suffix_grdid
          character(*)       , intent(in) :: suffix_sizes


          character(2), dimension(4) :: bf_layer_char
          integer                    :: i
          character(len=14)          :: filename_format
          character(len=30)          :: sizes_filename
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
        !>@param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !>@param timestep_written
        !> integer identifying the timestep written
        !
        !>@param name_var
        !> table with the short name for the governing variables saved
        !> in the netcdf file
        !
        !>@param longname_var
        !> table with the long name for the governing variables saved
        !> in the netcdf file
        !
        !>@param unit_var
        !> table with the units of the governing variables saved
        !> in the netcdf file
        !
        !>@param x_min_interior
        !> x-coordinate of interior grid point next to the left
        !> boundary layer
        !
        !>@param y_min_interior
        !> y-coordinate of interior grid point next to the lower
        !> boundary layer
        !
        !>@param dx
        !> grid size along the x-coordinate
        !
        !>@param dy
        !> grid size along the y-coordinate
        !
        !>@param time
        !> time corresponding to the data for the grdpts and the nodes
        !--------------------------------------------------------------
        subroutine print_netcdf(
     $     this,
     $     timestep_written,
     $     name_var,
     $     longname_var,
     $     unit_var,
     $     x_min_interior,
     $     y_min_interior,
     $     dx,
     $     dy,
     $     time)

          implicit none

          class(bf_mainlayer)        , intent(in) :: this
          integer                    , intent(in) :: timestep_written
          character(*), dimension(ne), intent(in) :: name_var
          character(*), dimension(ne), intent(in) :: longname_var
          character(*), dimension(ne), intent(in) :: unit_var
          real(rkind)                , intent(in) :: x_min_interior
          real(rkind)                , intent(in) :: y_min_interior
          real(rkind)                , intent(in) :: dx
          real(rkind)                , intent(in) :: dy
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
     $            x_min_interior,
     $            y_min_interior,
     $            dx,dy,
     $            time)

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine print_netcdf      

      end module bf_mainlayer_class

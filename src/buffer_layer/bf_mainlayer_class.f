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

        use bf_mainlayer_abstract_class, only : bf_mainlayer_abstract
        use bf_sublayer_class          , only : bf_sublayer
        use bf_sublayers_merge_module  , only : merge_sublayers_N,
     $                                          merge_sublayers_S,
     $                                          merge_sublayers_E,
     $                                          merge_sublayers_W
        use parameters_constant        , only : N,S,E,W
        use parameters_input           , only : nx, ny, ne, debug
        use parameters_kind            , only : ikind, rkind

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
        type, extends(bf_mainlayer_abstract) :: bf_mainlayer

          contains

          procedure, pass :: ini
          procedure, pass :: add_sublayer
          procedure, pass :: merge_sublayers

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
        function add_sublayer(this,alignment)
     $     result(added_sublayer_ptr)

          implicit none
          
          class(bf_mainlayer)       , intent(inout) :: this
          integer, dimension(2,2)   , intent(in)    :: alignment
          type(bf_sublayer), pointer                :: added_sublayer_ptr

          type(bf_sublayer), pointer :: current_sublayer
          integer :: current_border
          logical :: is_sublayer_further
          logical :: list_ended


          !< allocate space for the sublayer and initialize it
          allocate(added_sublayer_ptr)
          call added_sublayer_ptr%ini()


          !< if there is already an element in the list of sublayers
          !> of the mainlayer, go through the list of sublayers and
          !> insert the newly sublayer where its alignment fits
          list_ended = .true.
          if(associated(this%head_sublayer)) then

             !< going through the sublayer list starts with the head
             !> of the mainlayer
             current_sublayer => this%head_sublayer

             !< if the mainlayer_id is N or S, the position of the
             !> added sublayer is determined by the [i_min,i_max]
             !> alignment of the sublayer
             if((this%mainlayer_id.eq.N).or.(this%mainlayer_id.eq.S)) then

                current_border = current_sublayer%element%alignment(1,2)
                is_sublayer_further = current_border.lt.alignment(1,1)

                if(is_sublayer_further) then
                   
                   do while (is_sublayer_further)

                      if(.not.associated(current_sublayer%next)) then
                         list_ended=.true.
                         exit
                      end if

                      current_sublayer => current_sublayer%next
                      current_border = current_sublayer%element%alignment(1,2)
                      is_sublayer_further = current_border.lt.alignment(1,1)
                      list_ended=.false.
                   
                   end do
                else
                   list_ended=.false.
                end if

             !< if the mainlayer_id is E or W, the position of the
             !> added sublayer is determined by the [j_min,j_max]
             !> alignment of the sublayer
             else
                if((this%mainlayer_id.eq.E).or.(this%mainlayer_id.eq.W)) then

                   current_border = current_sublayer%element%alignment(2,2)
                   is_sublayer_further = current_border.lt.alignment(2,1)
                   
                   if(is_sublayer_further) then

                      do while(is_sublayer_further)
                         
                         if(.not.associated(current_sublayer%next)) then
                            list_ended=.true.
                            exit
                         end if

                         current_sublayer => current_sublayer%next
                         current_border = current_sublayer%element%alignment(2,2)
                         is_sublayer_further = current_border.lt.alignment(2,1)
                         list_ended=.false.

                      end do
                   else
                      list_ended=.false.
                   end if 

             !< if the mainlayer_id is NE,NW,SE or SW, the newly sublayer
             !> is simply added to the main layer
                else
                   do while(associated(current_sublayer%next))
                      current_sublayer => current_sublayer%next
                   end do
                end if
             end if

             !< if the previous list search stopped in the middle of the list
             !> (i.e. associated(current_sublayer%next)==.true., the new
             !> sublayer should be inserted in the list
             if(.not.(list_ended)) then

                !< if the position where it should be inserted is not the first
                !> position (i.e. associated(current_sublayer%prev==.true.)),
                !> the added_sublayer should be inserted between 
                if(associated(current_sublayer%prev)) then
                   call insert_sublayer(
     $                  added_sublayer_ptr,
     $                  current_sublayer%prev,
     $                  current_sublayer)

                !< if the position where it should be inserted is the first
                !> position (i.e. associated(current_sublayer%prev==.false.)),
                !> the added_sublayer should be simply added to the top of the
                !> list
                else
                   call link_next(added_sublayer_ptr,current_sublayer)
                   this%head_sublayer => added_sublayer_ptr
                end if

             !< if the previous list search stopped at the end of the list
             !> (i.e. associated(current_sublayer%next)==.false., the new
             !> sublayer should be added at the end of the list
             else
                call link_next(current_sublayer,added_sublayer_ptr)
                this%tail_sublayer => added_sublayer_ptr
             end if

          !< if there is no element in the list of sublayers
          !> of the mainlayer, simply add the new element to the
          !> top of the mainlayer
          else
             this%head_sublayer => added_sublayer_ptr
             this%tail_sublayer => added_sublayer_ptr

          end if

          this%nb_sublayers = this%nb_sublayers+1

        end function add_sublayer


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
        !>@param nodes
        !> table encapsulating the data of the interior domain
        !
        !>@param alignment
        !> table identifying the final position of the merged sublayer
        !> compared to the interior domain
        !
        !>@param neighbors_i
        !> table identifying whether the final merged sublayer will have
        !> neighboring sublayers or not
        !--------------------------------------------------------------
        function merge_sublayers(
     $     this,
     $     sublayer1,
     $     sublayer2,
     $     nodes,
     $     alignment,
     $     neighbors_i)
     $     result(merged_sublayer)
        
          implicit none
        
        
          class(bf_mainlayer)                          , intent(inout) :: this
          type(bf_sublayer), pointer                   , intent(inout) :: sublayer1
          type(bf_sublayer), pointer                   , intent(inout) :: sublayer2
          real(rkind)   , dimension(nx,ny,ne), optional, intent(in)    :: nodes
          integer(ikind), dimension(2,2)     , optional, intent(in)    :: alignment
          logical       , dimension(4)       , optional, intent(in)    :: neighbors_i
          type(bf_sublayer), pointer                                   :: merged_sublayer
        
        
          integer :: loc1, loc2
          logical, dimension(4) :: neighbors


          !get the mainlayer ID to which the two sublayers belong
          loc1 = sublayer1%element%localization

          !check input arguments
          if(debug) then

             loc2 = sublayer2%element%localization
             
             
             !if the two sublayers do not belong to the same mainlayer,
             !merging the two is not possible
             if(loc1.ne.loc2) then
                print '(''bf_mainlayer_class'')'
                print '(''merge_sublayers'')'
                print '(''the sublayers do not belong to'')'
                print '(''the same mainlayer:'')'
                print '(''sublayer1%localization: '',I2)', loc1
                print '(''sublayer2%localization: '',I2)', loc2
                stop 'check the sublayers merged'
             end if
             
             
             !if the two sublayers do not belong to this mainlayer,
             !merging the sublayers is not possible
             if(loc1.ne.(this%mainlayer_id)) then
                print '(''bf_mainlayer_class'')'
                print '(''merge_sublayers'')'
                print '(''the sublayers do not belong to'')'
                print '(''the correct mainlayer:'')'
                print '(''sublayer ID: '',I2)', loc1
                print '(''mainlayer ID: '',I2)', this%mainlayer_id
                stop 'check the sublayers and the mainlayer'
             end if
             
             
             !if alignment is passed as optional argument, so should
             !be nodes and vice versa
             if((present(alignment).and..not.present(nodes)).or.
     $            (.not.present(alignment).and.present(nodes))) then
                print '(''bf_mainlayer_class'')'
                print '(''merge_sublayers'')'
                print '(''alignment and nodes arguments should be'')'
                print '(''supplied together'')'
                stop 'change arguments'
             end if

          end if


          !initialize the neighbors table
          if(present(neighbors_i)) then
             neighbors = neighbors_i
          else
             neighbors = [.false., .false., .false., .false.]
          end if


          !select the subroutine adapted to merge the sublayers
          select case(loc1)
            case(N)
               if(present(nodes)) then
                  merged_sublayer => merge_sublayers_N(
     $                 this,
     $                 sublayer1,
     $                 sublayer2,
     $                 nodes,
     $                 alignment,
     $                 neighbor_E_i=neighbors(E),
     $                 neighbor_W_i=neighbors(W))
               else
                  merged_sublayer => merge_sublayers_N(
     $                 this,
     $                 sublayer1,
     $                 sublayer2,
     $                 neighbor_E_i=neighbors(E),
     $                 neighbor_W_i=neighbors(W))
               end if

            case(S)
               if(present(nodes)) then
                  merged_sublayer => merge_sublayers_S(
     $                 this,
     $                 sublayer1,
     $                 sublayer2,
     $                 nodes,
     $                 alignment,
     $                 neighbor_E_i=neighbors(E),
     $                 neighbor_W_i=neighbors(W))
               else
                  merged_sublayer => merge_sublayers_S(
     $                 this,
     $                 sublayer1,
     $                 sublayer2,
     $                 neighbor_E_i=neighbors(E),
     $                 neighbor_W_i=neighbors(W))
               end if

            case(E)
               if(present(nodes)) then
                  merged_sublayer => merge_sublayers_E(
     $                 this,
     $                 sublayer1,
     $                 sublayer2,
     $                 nodes,
     $                 alignment,
     $                 neighbor_N_i=neighbors(N),
     $                 neighbor_S_i=neighbors(S))
               else
                  merged_sublayer => merge_sublayers_E(
     $                 this,
     $                 sublayer1,
     $                 sublayer2,
     $                 neighbor_N_i=neighbors(N),
     $                 neighbor_S_i=neighbors(S))
               end if

            case(W)
               if(present(nodes)) then
                  merged_sublayer => merge_sublayers_W(
     $                 this,
     $                 sublayer1,
     $                 sublayer2,
     $                 nodes,
     $                 alignment,
     $                 neighbor_N_i=neighbors(N),
     $                 neighbor_S_i=neighbors(S))
               else
                  merged_sublayer => merge_sublayers_W(
     $                 this,
     $                 sublayer1,
     $                 sublayer2,
     $                 neighbor_N_i=neighbors(N),
     $                 neighbor_S_i=neighbors(S))
               end if

            case default
               print '(''bf_mainlayer_class'')'
               print '(''merge_sublayers'')'
               print '(''it is not possible to merge sublayers'')'
               print '(''from a corner'')'
               stop 'check the sublayers and the mainlayer'

          end select


          !update the number of sublayers in the mainlayer
          this%nb_sublayers = this%nb_sublayers-1

        end function merge_sublayers



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

          type(bf_sublayer)  , pointer, intent(inout) :: inserted_sublayer
          type(bf_sublayer)  , pointer, intent(inout) :: prev_sublayer
          type(bf_sublayer)  , pointer, intent(inout) :: next_sublayer

          call link_next(prev_sublayer, inserted_sublayer)
          call link_prev(next_sublayer, inserted_sublayer)

        end subroutine insert_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine linking two elements of bf_sublayers
        !> by providing the current and the next element
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_sublayer object encapsulating the buffer layer and
        !> pointers to the previous and next sublayer elements
        !
        !>@param next_sublayer
        !> bf_sublayer object that should become the next 
        !> sublayer element
        !--------------------------------------------------------------
        subroutine link_next(this, next_sublayer)

          implicit none

          type(bf_sublayer), pointer, intent(inout) :: this
          type(bf_sublayer), pointer, intent(inout) :: next_sublayer

          this%next          => next_sublayer
          next_sublayer%prev => this

        end subroutine link_next


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine linking two elements of bf_sublayers
        !> by providing the current and the previous element
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_sublayer object encapsulating the buffer layer and
        !> pointers to the previous and next sublayer elements
        !
        !>@param prev_sublayer
        !> bf_sublayer object that should become the previous
        !> sublayer element
        !--------------------------------------------------------------
        subroutine link_prev(this, prev_sublayer)

          implicit none

          type(bf_sublayer), pointer, intent(inout) :: this
          type(bf_sublayer), pointer, intent(inout) :: prev_sublayer

          call link_next(prev_sublayer,this)

        end subroutine link_prev

      end module bf_mainlayer_class

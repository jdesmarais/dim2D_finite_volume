      !> @file
      !> module implementing the object encapsulating 
      !> a sublayer element, an element of a double
      !> chained list
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the sublyare element, an
      !> element of a double chained list saving a buffer
      !> layer
      !
      !> @date
      ! 11_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_sublayer_class

        use bf_layer_class     , only : bf_layer

        use parameters_constant, only : N,S,E,W,
     $                                  x_direction, y_direction,
     $                                  min_border
        use parameters_input   , only : nx,ny,ne
        use parameters_kind    , only : ikind, rkind

        implicit none

        private
        public :: bf_sublayer


        !> @class bf_sublayer
        !> class encapsulating the bf_layer/interior interface
        !> object but only its main attributes are implemented
        !>
        !> @param element
        !> buffer layer saved in the double chained list element
        !>
        !> @param next_bf_sublayer
        !> pointer of the next element of the double chained list
        !>
        !> @param prev_bf_sublayer
        !> pointer of the previous element of the double chained list
        !>
        !> @param ini
        !> nullify the pointers of the object
        !>
        !> @param link_next
        !> links the current sublayer element with the next sublayer
        !> element given
        !>
        !> @param link_prev
        !> links the current sublayer element with the previous sublayer
        !> element given
        !---------------------------------------------------------------
        type, extends(bf_layer) :: bf_sublayer

          type(bf_sublayer), pointer, private :: next
          type(bf_sublayer), pointer, private :: prev

          type(bf_sublayer), pointer, private :: neighbor1
          type(bf_sublayer), pointer, private :: neighbor2

          contains

          procedure, pass   :: ini

          procedure, pass   :: get_prev
          procedure, pass   :: get_next
          procedure, pass   :: set_prev
          procedure, pass   :: set_next
          procedure, pass   :: nullify_prev
          procedure, pass   :: nullify_next

          procedure, pass   :: get_neighbor
          procedure, pass   :: set_neighbor
          procedure, pass   :: nullify_neighbor
          procedure, nopass :: get_neighbor_index
          procedure, nopass :: link_neighbors

          procedure, nopass :: merge_sublayers

          procedure,   pass :: update_grdpts

          procedure, pass, private :: update_from_neighbor1
          procedure, pass, private :: update_from_neighbor2
          procedure, pass, private :: update_to_neighbor1
          procedure, pass, private :: update_to_neighbor2
          procedure, pass          :: update_from_neighbors
          procedure, pass          :: update_to_neighbors
          
          

        end type bf_sublayer

        contains

        !< initialize the sublayer
        subroutine ini(this, localization)

          implicit none

          class(bf_sublayer), intent(inout) :: this
          integer(ikind)    , intent(in)    :: localization

          call this%bf_layer%ini(localization)

          nullify(this%prev)
          nullify(this%next)

          nullify(this%neighbor1)
          nullify(this%neighbor2)

        end subroutine ini


        !< access the previous list element
        function get_prev(this)

          implicit none

          class(bf_sublayer), intent(in) :: this
          type(bf_sublayer) , pointer    :: get_prev

          if(associated(this%prev)) then
             get_prev => this%prev
          else
             nullify(get_prev)
          end if

        end function get_prev


        !< access the next list element
        function get_next(this)

          implicit none

          class(bf_sublayer), intent(in) :: this
          type(bf_sublayer) , pointer    :: get_next

          if(associated(this%next)) then
             get_next => this%next
          else
             nullify(get_next)
          end if

        end function get_next


        !< link with the previous element in the chained list
        subroutine set_prev(this, bf_sublayer_prev)

          implicit none

          class(bf_sublayer)         , intent(inout) :: this
          type(bf_sublayer) , pointer, intent(in)    :: bf_sublayer_prev

          this%prev => bf_sublayer_prev

        end subroutine set_prev


        !< link with the next element in the chained list
        subroutine set_next(this, bf_sublayer_next)

          implicit none

          class(bf_sublayer)         , intent(inout) :: this
          type(bf_sublayer) , pointer, intent(in)    :: bf_sublayer_next

          this%next => bf_sublayer_next

        end subroutine set_next


        !< nullify the previous list element
        subroutine nullify_prev(this)

          implicit none

          class(bf_sublayer), intent(inout) :: this

          nullify(this%prev)

        end subroutine nullify_prev


        !< nullify the next list element
        subroutine nullify_next(this)

          implicit none

          class(bf_sublayer), intent(inout) :: this

          nullify(this%next)

        end subroutine nullify_next


        !< get a pointer to the neighbor1 or neighbor2 attribute
        function get_neighbor(this, index) result(neighbor_ptr)

          implicit none

          class(bf_sublayer), intent(in) :: this
          integer           , intent(in) :: index
          type(bf_sublayer) , pointer    :: neighbor_ptr

          select case(index)
            case(1)
               if(associated(this%neighbor1)) then
                  neighbor_ptr => this%neighbor1
               else
                  nullify(neighbor_ptr)
               end if

            case(2)
               if(associated(this%neighbor2)) then
                  neighbor_ptr => this%neighbor2
               else
                  nullify(neighbor_ptr)
               end if

            case default
               print '(''bf_sublayer_class'')'
               print '(''get_neighbor'')'
               print '(''index not recognized'')'
               print '(''index '',I2)', index
               stop 'change index for 1 or 2'
          end select
          
        end function get_neighbor


        !< set the neighbor1 or neighbor2 attribute
        subroutine set_neighbor(this, index, neighbor_ptr)

          implicit none

          class(bf_sublayer)         , intent(inout) :: this
          integer                    , intent(in)    :: index
          type(bf_sublayer) , pointer, intent(in)    :: neighbor_ptr

          select case(index)
            case(1)
               if(associated(neighbor_ptr)) then
                  this%neighbor1 => neighbor_ptr
               else
                  nullify(this%neighbor1)
               end if

            case(2)
               if(associated(neighbor_ptr)) then
                  this%neighbor2 => neighbor_ptr
               else
                  nullify(this%neighbor2)
               end if

            case default
               print '(''bf_sublayer_class'')'
               print '(''set_neighbor'')'
               print '(''index not recognized'')'
               print '(''index '',I2)', index
               stop 'change index for 1 or 2'
          end select
          
        end subroutine set_neighbor


        !< nullify the neighbor1 or neighbor2 attribute
        subroutine nullify_neighbor(this, index)

          implicit none

          class(bf_sublayer), intent(inout) :: this
          integer           , intent(in)    :: index

          select case(index)
            case(1)
               nullify(this%neighbor1)
            case(2)
               nullify(this%neighbor2)
            case default
               print '(''bf_sublayer_class'')'
               print '(''nullify_neighbor'')'
               print '(''index not recognized'')'
               print '(''index '',I2)', index
               stop 'change index for 1 or 2'
          end select

        end subroutine nullify_neighbor


        !< get the index (1 or 2) identifying the neighbor
        !> attribute (this%neighbor1 or this%neighbor2) that
        !> can be linked to a sublayer identified by its cardinal
        !> coordinate
        function get_neighbor_index(localization)
     $     result(neighbor_index)

          implicit none

          integer, intent(in) :: localization
          integer             :: neighbor_index


          !for a north sublayer(N), neighbor1 -> W
          !                         neighbor2 -> E
          !
          !for a south sublayer(S), neighbor1 -> W
          !                         neighbor2 -> E
          !
          !for a west sublayer (W), neighbor1 -> S
          !                         neighbor2 -> N
          !
          !for a west sublayer (E), neighbor1 -> S
          !                         neighbor2 -> N
          !---------------------------------------
          select case(localization)
            case(N,E)
               neighbor_index = 1
            case(S,W)
               neighbor_index = 2
            case default
               print '(''bf_sublayer_class'')'
               print '(''get_neighbor_index'')'
               print '(''localization not recognized'')'
               print '(''localization :'',I2)', localization
               stop 'change localization'
          end select

        end function get_neighbor_index


        !< set the current buffer layer pointer to a neighboring
        !> buffer layer and link the neighboring buffer layer to
        !> the current sublayer
        subroutine link_neighbors(bf_sublayer1, nbf_sublayer)

          implicit none

          type(bf_sublayer), pointer, intent(inout) :: bf_sublayer1
          type(bf_sublayer), pointer, intent(inout) :: nbf_sublayer

          integer :: neighbor_index

          !get the index (1 or 2) identifying the neighbor
          !attribute (bf_sublayer%neighbor1 or bf_sublayer%neighbor2) that
          !will be linked to nbf_sublayer
          neighbor_index = get_neighbor_index(nbf_sublayer%get_localization())
          call bf_sublayer1%set_neighbor(neighbor_index, nbf_sublayer)
          
          !get the index (1 or 2) identifying the neighbor
          !attribute (nbf_sublayer%neighbor1 or nbf_sublayer%neighbor2)
          !that will be linked to bf_sublayer
          neighbor_index = get_neighbor_index(bf_sublayer1%get_localization())
          call nbf_sublayer%set_neighbor(neighbor_index, bf_sublayer1)

        end subroutine link_neighbors



        !< merge sublayers
        subroutine merge_sublayers(bf_sublayer1, bf_sublayer2, nodes, alignment)

          implicit none

          type(bf_sublayer), pointer                 , intent(inout) :: bf_sublayer1
          type(bf_sublayer), pointer                 , intent(inout) :: bf_sublayer2
          real(rkind)      , dimension(nx,ny,ne)     , intent(in)    :: nodes
          integer(ikind)   , dimension(2,2), optional, intent(in)    :: alignment

          integer :: direction_tested


          !reorganize the position of the elements in the chained
          !list of sublayers
          select case(bf_sublayer1%get_localization())
            case(N,S)
               direction_tested = x_direction
            case(E,W)
               direction_tested = y_direction
            case default
               print '(''bf_sublayer_class'')'
               print '(''merge_sublayers'')'
               print '(''corner sublayers not eligible for merge'')'
               print '(''localization: '',I2)', bf_sublayer1%get_localization()
               print '(''check localization'')'
          end select

          if(bf_sublayer1%get_alignment(direction_tested,min_border).lt.
     $         bf_sublayer2%get_alignment(direction_tested,min_border)) then
             
             if(associated(bf_sublayer2%next)) then
                bf_sublayer1%next => bf_sublayer2%next
                bf_sublayer2%next%prev => bf_sublayer1
             else
                nullify(bf_sublayer1%next)
             end if
             
          else
             
             if(associated(bf_sublayer2%prev)) then
                bf_sublayer1%prev => bf_sublayer2%prev
                bf_sublayer2%prev%next => bf_sublayer1
             else
                nullify(bf_sublayer1%prev)
             end if
             
          end if


          !update the pointers to the neighboring buffer layers
          if(.not.associated(bf_sublayer1%neighbor1)) then
             if(associated(bf_sublayer2%neighbor1)) then
                call link_neighbors(bf_sublayer1, bf_sublayer2%neighbor1)
             end if
          end if
          if(.not.associated(bf_sublayer1%neighbor2)) then
             if(associated(bf_sublayer2%neighbor2)) then
                call link_neighbors(bf_sublayer1, bf_sublayer2%neighbor2)
             end if
          end if          


          !merge the attributes specific to the bf_layer:
          !localization, alignment, nodes, grdpts_id
          if(present(alignment)) then
             call bf_sublayer1%merge_bf_layer(bf_sublayer2, nodes, alignment)
          else
             call bf_sublayer1%merge_bf_layer(bf_sublayer2, nodes)
          end if

          
          !destroy the bf_sublayer2
          nullify(bf_sublayer2%next)
          nullify(bf_sublayer2%prev)
          deallocate(bf_sublayer2)
            
        end subroutine merge_sublayers


        !< update the grdpts of the buffer sublayer:
        !>
        !> update the nodes: compute the new nodes resulting
        !> from a buffer layer increase
        !>
        !> update the grdpts: update the boundary points
        !> in case of buffer layer increase
        !>
        !> update the neighbor grdpts : update the nodes and the grdpts
        !> of the neighboring buffer layer
        subroutine update_grdpts(this, selected_grdpts)

          implicit none

          class(bf_sublayer)            , intent(inout) :: this
          integer(ikind), dimension(:,:), intent(in)    :: selected_grdpts


          !update the nodes and grdpts_id of the current buffer layer
          call this%bf_layer%update_grdpts(selected_grdpts)

          !update the nodes and grdpts of the neighboring buffer layer
          !that have grid points is common with the current buffer layer
          call this%update_to_neighbors()

        end subroutine update_grdpts



        !< update the content of the sublayer tables
        !> by asking data from the neighbor1
        subroutine update_from_neighbor1(this)

          implicit none

          class(bf_sublayer), intent(inout) :: this

          if(associated(this%neighbor1)) then
             call this%copy_from_neighbor1(this%neighbor1)
          end if

        end subroutine update_from_neighbor1


        !< update the content of the sublayer tables
        !> by asking data from the neighbor2
        subroutine update_from_neighbor2(this)

          implicit none

          class(bf_sublayer), intent(inout) :: this

          if(associated(this%neighbor2)) then
             call this%copy_from_neighbor2(this%neighbor2)
          end if

        end subroutine update_from_neighbor2


        !< update the content of neighbor1
        !> by transfering data from the current
        !> buffer layer
        subroutine update_to_neighbor1(this)

          implicit none

          class(bf_sublayer), intent(inout) :: this

          if(associated(this%neighbor1)) then
             call this%copy_to_neighbor1(this%neighbor1)
          end if

        end subroutine update_to_neighbor1


        !< update the content of neighbor2
        !> by transfering data from the current
        !> buffer layer
        subroutine update_to_neighbor2(this)

          implicit none

          class(bf_sublayer), intent(inout) :: this

          if(associated(this%neighbor2)) then
             call this%copy_to_neighbor2(this%neighbor2)
          end if

        end subroutine update_to_neighbor2


        !< update the content of the sublayer tables
        !> by asking data from the neighbors
        subroutine update_from_neighbors(this)

          implicit none

          class(bf_sublayer), intent(inout) :: this

          call this%update_from_neighbor1()
          call this%update_from_neighbor2()

        end subroutine update_from_neighbors


        !< update the content of the neighboring
        !> buffer layers by transfering data from
        !> the current buffer layer
        subroutine update_to_neighbors(this)

          implicit none

          class(bf_sublayer), intent(inout) :: this

          call this%update_to_neighbor1()
          call this%update_to_neighbor2()

        end subroutine update_to_neighbors


      end module bf_sublayer_class

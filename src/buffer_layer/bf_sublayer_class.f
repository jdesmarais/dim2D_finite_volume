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

        use bf_layer_class, only : bf_layer

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
        type bf_sublayer

          type(bf_layer)             :: element
          type(bf_sublayer), pointer :: next
          type(bf_sublayer), pointer :: prev

          contains

          procedure, pass :: ini

        end type bf_sublayer

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the element of the sublayer
        !> by nullifying the pointers
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_sublayer object encapsulating the buffer layer and
        !> pointers to the previous and next sublayer elements
        !--------------------------------------------------------------
        subroutine ini(this)
        
          implicit none

          class(bf_sublayer), intent(inout) :: this

          nullify(this%next)
          nullify(this%prev)

        end subroutine ini


      end module bf_sublayer_class

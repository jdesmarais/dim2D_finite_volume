      !> @file
      !> module implementing the element of a doubled chained list
      !> of bf_sublayer pointers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the element of a doubled chained list
      !> of bf_sublayer pointers
      !
      !> @date
      ! 27_06_2014 - documentation update     - J.L. Desmarais
      ! 30_10_2014 - synchronization of nodes - J.L. Desmarais
      !-----------------------------------------------------------------
      module nbf_element_class

        use bf_layer_class        , only : bf_layer
        use bf_layer_errors_module, only : error_incompatible_neighbor
        use bf_sublayer_class     , only : bf_sublayer
        use parameters_bf_layer   , only : align_N, align_S
        use parameters_constant   , only : N,S,E,W,
     $                                     x_direction, y_direction
        use parameters_input      , only : bc_size
        use parameters_kind       , only : ikind

        implicit none

        private
        public :: nbf_element


        !>@class nbf_element
        !< element of a double chained list of bf_sublayer
        !> pointers
        !
        !>@param prev
        !> pointer to the previous element in the list
        !
        !>@param next
        !> pointer to the next element in the list
        !
        !>@param ptr
        !> pointer to the bf_sublayer object
        !
        !>@param ini
        !> initialize the nbf_element object by nullify its
        !> ptr attribute and its links to previous and next
        !> elements in the list
        !
        !>@param set_prev
        !> set the link to the previous element of the list
        !
        !>@param set_next
        !> set the link to the next element of the list
        !
        !>@param set_ptr
        !> set the ptr attribute
        !
        !>@param get_prev
        !> access the previous element of the list
        !
        !>@param get_next
        !> access the next element of the list
        !
        !>@param get_ptr
        !> get the ptr attribute
        !
        !>@param nullify_prev
        !> remove the link to the previous element of the list
        !
        !>@param nullify_next
        !> remove the link to the next element of the list
        !
        !>@param is_before
        !> check if the current nbf_element object is before 
        !> the other element in the doubled chained list
        !
        !>@param refers_to
        !> check if the ptr attribute refers to the
        !> buffer sublayer passed as argument
        !
        !> @param copy_from_neighbor1_to
        !> copy the common layer to the ptr attribute
        !> from its neighboring buffer layer identified as of type 1
        !
        !> @param copy_from_neighbor2_to
        !> copy the common layer to the ptr attribute
        !> from its neighboring buffer layer identified as of type 2
        !
        !> @param copy_to_neighbor1_from
        !> copy the common layer from the ptr attribute
        !> to its neighboring buffer layer identified as of type 1
        !
        !> @param copy_to_neighbor2_from
        !> copy the common layer from the ptr attribute
        !> to its neighboring buffer layer identified as of type 2
        !
        !> @param shares_grdpts_along_x_dir_with
        !> check if a neighboring buffer layer (positioned such that
        !> it is either a potential neighboring buffer layer of type
        !> 1 or 2) has indeed grid points in common with the ptr attribute
        !> layer by computing the x-size of the layer to be exchanged
        !
        !> @param get_remain_status
        !> get the can_remain attribute of the ptr attribute
        !
        !> @param sync_nodes_with_neighbor1
        !> synchronize the nodes at the interface between the buffer
        !> layer and another buffer layer which is a neighbor of type 1
        !
        !> @param sync_nodes_with_neighbor2
        !> synchronize the nodes at the interface between the buffer
        !> layer and another buffer layer which is a neighbor of type 2
        !---------------------------------------------------------------
        type :: nbf_element

          type(nbf_element), pointer, private :: prev
          type(nbf_element), pointer, private :: next
          type(bf_sublayer), pointer, private :: ptr

          contains

          procedure, pass :: ini
          procedure, pass :: set_prev
          procedure, pass :: set_next
          procedure, pass :: set_ptr
          procedure, pass :: get_prev
          procedure, pass :: get_next
          procedure, pass :: get_ptr
          procedure, pass :: nullify_prev
          procedure, pass :: nullify_next
          procedure, pass :: nullify_ptr
          procedure, pass :: is_before
          procedure, pass :: refers_to
          procedure, pass :: copy_from_neighbor1_to
          procedure, pass :: copy_from_neighbor2_to
          procedure, pass :: copy_to_neighbor1_from
          procedure, pass :: copy_to_neighbor2_from
          procedure, pass :: shares_grdpts_along_x_dir_with
          procedure, pass :: get_remain_status

          procedure, pass :: sync_nodes_with_neighbor1
          procedure, pass :: sync_nodes_with_neighbor2

        end type nbf_element

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the nbf_element object by nullify its
        !> ptr attribute and its links to previous and next
        !> elements in the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param added_bf_sublayer_ptr
        !> pointer to the bf_sublayer stored in the element
        !--------------------------------------------------------------
        subroutine ini(this, added_bf_sublayer_ptr)

          implicit none

          class(nbf_element)         , intent(inout) :: this
          type(bf_sublayer) , pointer, intent(in)    :: added_bf_sublayer_ptr

          if(associated(added_bf_sublayer_ptr)) then
             this%ptr => added_bf_sublayer_ptr
          else
             nullify(this%ptr)
          end if

          nullify(this%prev)
          nullify(this%next)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the link to the previous element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param prev_ptr
        !> pointer to the previous element of the list
        !--------------------------------------------------------------
        subroutine set_prev(this, prev_ptr)

          implicit none

          class(nbf_element)        , intent(inout) :: this
          type(nbf_element), pointer, intent(in)    :: prev_ptr

          this%prev => prev_ptr

        end subroutine set_prev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the link to the next element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param next_ptr
        !> pointer to the next element of the list
        !--------------------------------------------------------------
        subroutine set_next(this, next_ptr)

          implicit none

          class(nbf_element)        , intent(inout) :: this
          type(nbf_element), pointer, intent(in)    :: next_ptr

          this%next => next_ptr

        end subroutine set_next


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the ptr attribute
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param added_bf_sublayer_ptr
        !> pointer to the bf_sublayer stored in the element
        !--------------------------------------------------------------
        subroutine set_ptr(this, added_bf_sublayer_ptr)

          implicit none

          class(nbf_element)        , intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: added_bf_sublayer_ptr

          if(associated(added_bf_sublayer_ptr)) then
             this%ptr => added_bf_sublayer_ptr
          end if

        end subroutine set_ptr
      
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> access the previous element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@return get_prev
        !> pointer to the previous element of the list
        !--------------------------------------------------------------
        function get_prev(this)

          implicit none

          class(nbf_element), intent(in) :: this
          type(nbf_element) , pointer    :: get_prev

          get_prev => this%prev

        end function get_prev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> access the next element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@return get_next
        !> pointer to the next element of the list
        !--------------------------------------------------------------
        function get_next(this)

          implicit none

          class(nbf_element), intent(in) :: this
          type(nbf_element) , pointer    :: get_next

          get_next => this%next

        end function get_next


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the ptr attribute
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@return get_ptr
        !> pointer to the bf_sublayer stored in the element
        !--------------------------------------------------------------
        function get_ptr(this)

          implicit none

          class(nbf_element), intent(inout) :: this
          type(bf_sublayer) , pointer       :: get_ptr

          get_ptr => this%ptr

        end function get_ptr


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the link to the previous element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !--------------------------------------------------------------
        subroutine nullify_prev(this)

          implicit none

          class(nbf_element), intent(inout) :: this

          nullify(this%prev)

        end subroutine nullify_prev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the link to the next element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !--------------------------------------------------------------
        subroutine nullify_next(this)

          implicit none

          class(nbf_element), intent(inout) :: this

          nullify(this%next)

        end subroutine nullify_next
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the link to the bf_sublayer object
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !--------------------------------------------------------------
        subroutine nullify_ptr(this)

          implicit none

          class(nbf_element), intent(inout) :: this

          nullify(this%ptr)

        end subroutine nullify_ptr

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the current nbf_element object is before 
        !> the other element in the doubled chained list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param element2
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@return is_before
        !> logical stating whether the current element should be placed
        !> before the element2 in the doubled chained list
        !--------------------------------------------------------------
        function is_before(this, element2)

          implicit none

          class(nbf_element), intent(in) :: this
          class(nbf_element), intent(in) :: element2
          logical                        :: is_before

          if(associated(this%ptr)) then
             if(associated(element2%ptr)) then
                
                is_before = 
     $               (this%ptr%get_alignment(x_direction,1)).le.
     $               (element2%ptr%get_alignment(x_direction,1))


             else
                print '(''nbf_element_class'')'
                print '(''is_before'')'
                print '(''element2%ptr not associated'')'
                stop 'cannot compare elements'
             end if
          else
             print '(''nbf_element_class'')'
             print '(''is_before'')'
             print '(''this%ptr not associated'')'
             stop 'is element initialized ?'
          end if

        end function is_before


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the ptr attribute refers to the
        !> buffer sublayer passed as argument
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param added_bf_sublayer_ptr
        !> pointer to the bf_sublayer tested
        !
        !>@return refers_to
        !> logical stating whether the bf_sublayer pointer refers also
        !> to the element pointed by added_bf_sublayer_ptr
        !--------------------------------------------------------------
        function refers_to(this, added_bf_sublayer_ptr)

          implicit none

          class(nbf_element)        , intent(in) :: this
          type(bf_sublayer), pointer, intent(in) :: added_bf_sublayer_ptr
          logical                                :: refers_to

          if(associated(this%ptr)) then
             if(associated(added_bf_sublayer_ptr)) then
                
                refers_to = associated(this%ptr, added_bf_sublayer_ptr)

             else
                print '(''nbf_element_class'')'
                print '(''is_same_element'')'
                print '(''added_bf_sublayer_ptr not associated'')'
                stop 'cannot compare elements'
             end if
          else
             print '(''nbf_element_class'')'
             print '(''is_same_element'')'
             print '(''this%ptr not associated'')'
             stop 'is element initialized ?'
          end if

        end function refers_to


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layer to the ptr attribute
        !> from its neighboring buffer layer identified as of type 1
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param bf_exchanged
        !> bf_layer with which the current element is exchanging data
        !--------------------------------------------------------------
        subroutine copy_from_neighbor1_to(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(in)    :: this
          class(bf_layer)   , intent(inout) :: bf_exchanged

          !if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_from_neighbor1(this%ptr)
          !end if

        end subroutine copy_from_neighbor1_to


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layer to the ptr attribute
        !> from its neighboring buffer layer identified as of type 2
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param bf_exchanged
        !> bf_layer with which the current element is exchanging data
        !--------------------------------------------------------------
        subroutine copy_from_neighbor2_to(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(in)    :: this
          class(bf_layer)   , intent(inout) :: bf_exchanged

          !if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_from_neighbor2(this%ptr)
          !end if

        end subroutine copy_from_neighbor2_to


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layer from the ptr attribute
        !> to its neighboring buffer layer identified as of type 1
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param bf_exchanged
        !> bf_layer with which the current element is exchanging data
        !--------------------------------------------------------------
        subroutine copy_to_neighbor1_from(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(inout) :: this
          class(bf_layer)   , intent(in)    :: bf_exchanged

          !if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_to_neighbor1(this%ptr)
          !end if

        end subroutine copy_to_neighbor1_from


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layer from the ptr attribute
        !> to its neighboring buffer layer identified as of type 2
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param bf_exchanged
        !> bf_layer with which the current element is exchanging data
        !--------------------------------------------------------------
        subroutine copy_to_neighbor2_from(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(inout) :: this
          class(bf_layer)   , intent(in)    :: bf_exchanged

          !if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_to_neighbor2(this%ptr)
          !end if

        end subroutine copy_to_neighbor2_from


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if a neighboring buffer layer (positioned such that
        !> it is either a potential neighboring buffer layer of type
        !> 1 or 2) has indeed grid points in common with the ptr attribute
        !> by computing the x-size of the layer to be exchanged
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param neighbor
        !> bf_layer with which the current element is exchanging data
        !
        !>@return share
        !> logical identifying if the ptr attribute and the bf_layer
        !> object are sharing grid points
        !--------------------------------------------------------------
        function shares_grdpts_along_x_dir_with(this, neighbor)
     $     result(share)

          class(nbf_element), intent(in) :: this
          class(bf_layer)   , intent(in) :: neighbor
          logical                        :: share

          share = this%ptr%shares_grdpts_along_x_dir_with(neighbor)
          
        end function shares_grdpts_along_x_dir_with


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the can_remain attribute of the ptr attribute
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@return get_remain_status
        !> can_remain attribute of the object pointed by the ptr
        !> attribute
        !--------------------------------------------------------------
        function get_remain_status(this)

          implicit none

          class(nbf_element), intent(in) :: this
          logical                        :: get_remain_status

          get_remain_status = this%ptr%get_remain_status()

        end function get_remain_status


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> synchronize the nodes at the interface between the buffer
        !> layer and another buffer layer which is a neighbor of type 1
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param this2
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !--------------------------------------------------------------
        subroutine sync_nodes_with_neighbor1(this,this2)

          implicit none

          class(nbf_element), intent(inout) :: this
          class(nbf_element), intent(inout) :: this2

          call this%ptr%sync_nodes_with_neighbor1(this2%ptr)

        end subroutine sync_nodes_with_neighbor1

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> synchronize the nodes at the interface between the buffer
        !> layer and another buffer layer which is a neighbor of type 2
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !
        !>@param this2
        !> nbf_element object encapsulating the references to
        !> bf_sublayer object and the previous and next elements
        !> in the doubled chained list
        !--------------------------------------------------------------
        subroutine sync_nodes_with_neighbor2(this,this2)

          implicit none

          class(nbf_element), intent(inout) :: this
          class(nbf_element), intent(inout) :: this2

          call this%ptr%sync_nodes_with_neighbor2(this2%ptr)

        end subroutine sync_nodes_with_neighbor2

      end module nbf_element_class

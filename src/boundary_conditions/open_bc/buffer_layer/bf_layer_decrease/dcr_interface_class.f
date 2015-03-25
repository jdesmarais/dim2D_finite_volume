      !> @file
      !> object responsible for managing the
      !> removal of buffer layers in the domain
      !> extension
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> object responsible for managing the
      !> removal of buffer layers in the domain
      !> extension
      !
      !> @date
      ! 25_03_2015 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module dcr_interface_class

        use bf_interface_icr_class, only :
     $       bf_interface_icr

        use bf_mainlayer_class, only :
     $       bf_mainlayer

        use bf_sublayer_class, only :
     $       bf_sublayer

        use dcr_list_class, only :
     $       dcr_list

        implicit none

        private
        public :: dcr_interface


        !>@class dcr_interface
        !> object responsible for managing the removal
        !> of buffer layers in the domain extension
        !
        !>@param no_check_list
        !> pointers references to the buffer layers that cannot be
        !> removed b/c one of its neighbor should remain. It is
        !> not necessary for these references to know whether in the
        !> absence of neighbors, they would be removed.
        !
        !>@param double_check_list
        !> pointers references to the buffer layers that could be
        !> removed if its neighbors are not taking into account.
        !> The references are organized in mainlayers (N,S,E,W)
        !
        !>@param ini
        !> initialize the main attributes and especially the maximum
        !> size for the arrays storing the references to the buffer
        !> layers
        !
        !>@param stage
        !> stage a buffer layer for removal and check whether it can be
        !> removed directly or whether its removal is conditioned by
        !> its neighbors
        !
        !>@param finalize_domain_decrease
        !> remove the buffer layer that have been staged and whose
        !> removal is approved by its neighbors
        !
        !>@param not_in_no_check_list
        !> check whether the bf_sublayer_ptr passed as argument belongs
        !> to the list of buffer layers that should not be checked
        !
        !>@param prevent_neighbor_removal
        !> if the current buffer layer cannot be removed, the removal
        !> of its neighbors should be prevented
        !
        !>@param check_if_neighbors_remain
        !> verify whether the neighbors of the current buffer layer
        !> can be removed
        !--------------------------------------------------------------
        type :: dcr_interface

          type(dcr_list), dimension(4) :: no_check_list
          type(dcr_list), dimension(4) :: double_check_list

          contains

          ! main procedures to stage and remove buffer layers
          procedure,   pass :: ini
          procedure,   pass :: stage
          procedure,   pass :: finalize_domain_decrease

          ! annex procedures needed to determine whether the
          ! removal of buffer layers is approved
          procedure,   pass :: not_in_no_check_list
          procedure,   pass :: prevent_neighbor_removal
          procedure, nopass :: check_if_neighbors_remain

        end type dcr_interface


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the main attributes and especially the maximum
        !> size for the arrays storing the references to the buffer
        !> layers
        !
        !> @date
        !> 25_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object responsible for managing the removal
        !> of buffer layers in the domain extension
        !
        !>@param bf_interface_used
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !--------------------------------------------------------------
        subroutine ini(this,bf_interface_used)

          implicit none

          class(dcr_interface)   , intent(inout) :: this
          class(bf_interface_icr), intent(in)    :: bf_interface_used

          
          integer                     :: k
          type(bf_mainlayer), pointer :: bf_mainlayer_ptr
          integer                     :: nb_sublayers


          !for each main layers, determine the number of sublayers
          !stored. This defines for each dcr_list the maximum number
          !of references stored
          do k=1,4

             bf_mainlayer_ptr => bf_interface_used%get_mainlayer_ptr(k)

             if(associated(bf_mainlayer_ptr)) then
                nb_sublayers = bf_mainlayer_ptr%get_nb_sublayers()
                call this%no_check_list(k)%ini(nb_sublayers)
                call this%double_check_list(k)%ini(nb_sublayers)
             else
                call this%no_check_list(k)%ini()
                call this%double_check_list(k)%ini()
             end if             

          end do

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the bf_sublayer_ptr passed as argument
        !> should be check, i.e. if it does not belong to the
        !> no-check_list
        !
        !> @date
        !> 25_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object responsible for managing the removal
        !> of buffer layers in the domain extension
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying the main layer to which the
        !> bf_sublayer_ptr belongs
        !
        !>@param bf_sublayer_ptr
        !> pointer to a bf_sublayer object
        !--------------------------------------------------------------
        function not_in_no_check_list(this,mainlayer_id,bf_sublayer_ptr)

          implicit none

          class(dcr_interface)      , intent(inout) :: this
          integer                   , intent(in)    :: mainlayer_id
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_ptr
          logical                                   :: not_in_no_check_list


          not_in_no_check_list = this%no_check_list(mainlayer_id)%does_not_contain(bf_sublayer_ptr)

        end function not_in_no_check_list


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> prevent the removal of the neighboring buffer layers of the
        !> buffer layer passed as arguments: their remain_status is
        !> forced to be set to .true. and they should not be checked
        !> to know whether their removal is allowed apart from their
        !> neighbors
        !
        !> @date
        !> 25_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object responsible for managing the removal
        !> of buffer layers in the domain extension
        !
        !>@param bf_interface_used
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying the main layer to which the
        !> bf_sublayer_ptr belongs
        !
        !>@param bf_sublayer_ptr
        !> pointer to a bf_sublayer object
        !--------------------------------------------------------------
        subroutine prevent_neighbor_removal(
     $     this,
     $     bf_interface_used,
     $     mainlayer_id,
     $     bf_sublayer_ptr)

          implicit none

          class(dcr_interface)      , intent(inout) :: this
          class(bf_interface_icr)   , intent(in)    :: bf_interface_used
          integer                   , intent(in)    :: mainlayer_id
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_ptr


          type(bf_sublayer), pointer :: nbf_sublayer_ptr
          integer                    :: nbf_mainlayer_id


          !> check whether the current buffer layer is exchanging
          !> grid-points with a neighbor of type 1
          if(bf_sublayer_ptr%can_exchange_with_neighbor1()) then

             !> extract the reference to the neighbor of type 1
             nbf_sublayer_ptr => bf_interface_used%get_neighbor_sublayer(
     $            mainlayer_id,
     $            1)

             !> force the remain_status of the neighbor to .true.
             call nbf_sublayer_ptr%set_remain_status(.true.)

             !> make sure the remain_status will not be overwritten
             !> by determing the removal status without considering
             !> the neighbors. For this purpose, the reference to
             !> nbf_sublayer_ptr is added to the list of sublayers
             !> that will not be checked later on
             call bf_sublayer_ptr%get_neighbor1_id(nbf_mainlayer_id)
             call this%no_check_list(nbf_mainlayer_id)%add_ele(
     $            nbf_sublayer_ptr)

          end if


          !> check whether the current buffer layer is exchanging
          !> grid-points with a neighbor of type 2
          if(bf_sublayer_ptr%can_exchange_with_neighbor2()) then

             !> extract the reference to the neighbor of type 2
             nbf_sublayer_ptr => bf_interface_used%get_neighbor_sublayer(
     $            mainlayer_id,
     $            2)

             !> force the remain_status of the neighbor to .true.
             call nbf_sublayer_ptr%set_remain_status(.true.)

             !> make sure the remain_status will not be overwritten
             !> by determing the removal status without considering
             !> the neighbors. For this purpose, the reference to
             !> nbf_sublayer_ptr is added to the list of sublayers
             !> that will not be checked later on
             call bf_sublayer_ptr%get_neighbor2_id(nbf_mainlayer_id)
             call this%no_check_list(nbf_mainlayer_id)%add_ele(
     $            nbf_sublayer_ptr)

          end if

        end subroutine prevent_neighbor_removal


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> the buffer layer is staged for removal. If this buffer
        !> layer has neighbors, its removal is conditioned by the
        !> removal of its neighbors
        !
        !> @date
        !> 25_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object responsible for managing the removal
        !> of buffer layers in the domain extension
        !
        !>@param bf_interface_used
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying the main layer to which
        !> the buffe rlayer bf_sublayer_ptr is belonging
        !
        !>@param bf_sublayer_ptr
        !> reference to the buffer layer staged for removal
        !--------------------------------------------------------------
        subroutine stage(
     $     this,
     $     bf_interface_used,
     $     mainlayer_id,
     $     bf_sublayer_ptr)

          implicit none

          class(dcr_interface)      , intent(inout) :: this
          class(bf_interface_icr)   , intent(inout) :: bf_interface_used
          integer                   , intent(in)    :: mainlayer_id
          type(bf_sublayer), pointer, intent(inout) :: bf_sublayer_ptr

          
          logical :: neighbor_dependent


          !> determine whether the removal of the buffer layer depends
          !> on its neighbors, i.e. whether the buffer layer has some
          !> neighbors
          neighbor_dependent = 
     $         (bf_sublayer_ptr%can_exchange_with_neighbor1()).or.
     $         (bf_sublayer_ptr%can_exchange_with_neighbor2())


          !> if the buffer layer depends on its neighbors, its removal
          !> is postpone for the double checked (i.e. once the removal
          !> status of its neighbors is determined)
          if(neighbor_dependent) then

             call this%double_check_list(mainlayer_id)%add_ele(
     $            bf_sublayer_ptr)

          !> otherwise, the buffer layer is removed immediately
          else
             
             call bf_interface_used%remove_sublayer(
     $            bf_sublayer_ptr,
     $            mainlayer_id=mainlayer_id)

             nullify(bf_sublayer_ptr)

          end if

        end subroutine stage


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the neighbors of the buffer layer should
        !> be removed, unlocking the removal of the current
        !> buffer layer
        !
        !> @date
        !> 25_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object responsible for managing the removal
        !> of buffer layers in the domain extension
        !
        !>@param bf_interface_used
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying the main layer to which
        !> the buffe rlayer bf_sublayer_ptr is belonging
        !
        !>@param bf_sublayer_ptr
        !> reference to the buffer layer staged for removal
        !--------------------------------------------------------------
        function check_if_neighbors_remain(
     $     bf_interface_used,
     $     mainlayer_id,
     $     bf_sublayer_ptr)
     $     result(neighbors_remain)

          implicit none

          class(bf_interface_icr)   , intent(inout) :: bf_interface_used
          integer                   , intent(in)    :: mainlayer_id
          type(bf_sublayer), pointer, intent(inout) :: bf_sublayer_ptr
          logical                                   :: neighbors_remain


          type(bf_sublayer), pointer :: nbf_sublayer_ptr


          neighbors_remain = .false.


          !> check whether the neighbor of type 1, if it exists,
          !> prevents the removal of the current buffer layer
          if(bf_sublayer_ptr%can_exchange_with_neighbor1()) then

             !> extract the reference to the neighbor of type 1
             nbf_sublayer_ptr => bf_interface_used%get_neighbor_sublayer(
     $            mainlayer_id,
     $            1)

             !> extract the remain_status of the neighbor
             neighbors_remain = nbf_sublayer_ptr%get_remain_status()

          end if


          !> check whether the neighbor of type 2, if it exists,
          !> prevents the removal of the current buffer layer
          if(.not.neighbors_remain) then
             
             if(bf_sublayer_ptr%can_exchange_with_neighbor2()) then

                !> extract the reference to the neighbor of type 2
                nbf_sublayer_ptr => bf_interface_used%get_neighbor_sublayer(
     $               mainlayer_id,
     $               2)

                !> extract the remain_status of the neighbor
                neighbors_remain = nbf_sublayer_ptr%get_remain_status()

             end if

          end if


        end function check_if_neighbors_remain


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> the buffer layers stored in the double_check_list, i.e. the
        !> buffer layers that could be removed if there were no neighbors,
        !> are analyzed to check whether its neighbors can also be removed.
        !> If so, the buffer layer is removed
        !
        !> @date
        !> 25_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object responsible for managing the removal
        !> of buffer layers in the domain extension
        !
        !>@param bf_interface_used
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying the main layer to which
        !> the buffe rlayer bf_sublayer_ptr is belonging
        !
        !>@param bf_sublayer_ptr
        !> reference to the buffer layer staged for removal
        !--------------------------------------------------------------
        subroutine finalize_domain_decrease(
     $     this,
     $     bf_interface_used)

          implicit none

          class(dcr_interface)   , intent(inout) :: this
          class(bf_interface_icr), intent(inout) :: bf_interface_used


          integer                    :: k
          integer                    :: nb_sublayers
          integer                    :: m
          type(bf_sublayer), pointer :: bf_sublayer_ptr
          logical                    :: neighbors_remain


          !> loop over the cardinal coordinates to check the buffer
          !> layer references in the double_check_lists
          do k=1,4

             !> loop over the sublayers to be checked
             nb_sublayers = this%double_check_list(k)%get_nb_ele()

             do m=1, nb_sublayers

                bf_sublayer_ptr => this%double_check_list(k)%get_ele(m)

                !> check whether the neighbors of the buffer layer
                !> should not be removed
                neighbors_remain = check_if_neighbors_remain(
     $               bf_interface_used,
     $               k,
     $               bf_sublayer_ptr)

                !> if the neighbors can be removed, the current buffer
                !> layer can be removed
                if(.not.neighbors_remain) then

                   call bf_interface_used%remove_sublayer(
     $                  bf_sublayer_ptr,
     $                  mainlayer_id=k)

                end if

             end do

             !> the list containing the references to the sublayers
             !> can be removed for the cardinal coordinate
             call this%no_check_list(k)%remove()
             call this%double_check_list(k)%remove()

          end do

        end subroutine finalize_domain_decrease

      end module dcr_interface_class

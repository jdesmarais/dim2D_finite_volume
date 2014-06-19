      module bf_interface_dcr_class

        use bf_interface_icr_class, only : bf_interface_icr
        use bf_mainlayer_class    , only : bf_mainlayer
        use bf_sublayer_class     , only : bf_sublayer
        use parameters_input      , only : nx,ny,ne
        use parameters_kind       , only : rkind
        use sbf_list_class        , only : sbf_list

        implicit none

        private
        public :: bf_interface_dcr


        type, extends(bf_interface_icr) :: bf_interface_dcr

          contains

          procedure, pass :: update_bf_layers_with_detector_dcr

        end type bf_interface_dcr
        

        contains


        subroutine update_bf_layers_with_detector_dcr(
     $       this, interior_nodes)

          implicit none

          class(bf_interface_dcr)         , intent(inout) :: this
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes


          integer                      :: k,l
          type(bf_mainlayer), pointer  :: mainlayer_ptr
          integer                      :: nb_sublayers
          type(bf_sublayer) , pointer  :: sublayer_ptr
          type(bf_sublayer) , pointer  :: sublayer_next
          type(sbf_list), dimension(4) :: not_to_be_tested_sublayers_list
          type(sbf_list), dimension(4) :: to_be_rechecked_sublayers_list
          logical                      :: can_be_tested
          integer                      :: neighbor1_id
          integer                      :: neighbor2_id
          logical                      :: can_remain
          logical                      :: neighbor_dependent
          integer                      :: nb_sublayers_rechecked


          !0) initialization of the lists containining pointers
          !   to buffer sublayers
          do k=1,4
             
             !get the main layer corresponding to the cardinal
             !coordinate k
             mainlayer_ptr => this%get_mainlayer(k)

             !if there are sublayers inside this mainlayer,
             !the different sublayers are investigated to check
             !whether they should be removed
             if(associated(mainlayer_ptr)) then

                nb_sublayers = mainlayer_ptr%get_nb_sublayers()
                call not_to_be_tested_sublayers_list(k)%ini(nb_sublayers)
                call to_be_rechecked_sublayers_list(k)%ini(nb_sublayers)

             else

                call not_to_be_tested_sublayers_list(k)%ini()
                call to_be_rechecked_sublayers_list(k)%ini()

             end if

          end do


          !1) the buffer layers contained by the mainlayers
          !   are analyzed to see whether they can be removed
          !   if the neighboring buffer layers from which they
          !   depend are not taken in account: their logical
          !   attribute is set accordingly
          !     - we set not_to_be_tested_sublayers_list that it
          !       is only used during the loop
          !     - we set to_be_rechecked_sublayers_list that is
          !       used by the second loop
          do k=1,4

             !get the main layer corresponding to the cardinal
             !coordinate k
             mainlayer_ptr => this%get_mainlayer(k)

             !if there are sublayers inside this mainlayer,
             !the different sublayers are investigated to check
             !whether they should be removed
             if(associated(mainlayer_ptr)) then

                nb_sublayers = mainlayer_ptr%get_nb_sublayers()
                sublayer_ptr => mainlayer_ptr%get_head_sublayer()

                !loop over the sublayers contained in the mainlayer
                do l=1, nb_sublayers

                   !determine if the sublayer found should not be
                   !checked because it is a sublayer found to be 
                   !the neighboring buffer layer of one that cannot
                   !be removed
                   can_be_tested = not_to_be_tested_sublayers_list(k)%does_not_contain(sublayer_ptr)
                   
                   
                   sublayer_next => sublayer_ptr%get_next()

                   !if it can be tested
                   if(can_be_tested) then

                      !check if it is possible to remove the buffer layer
                      can_remain = sublayer_ptr%should_remain(interior_nodes)

                      !set the status to the buffer layer
                      call sublayer_ptr%set_remain_status(can_remain)

                      !if the buffer layer should remain, the neighboring
                      !buffer layers from which the buffer layer depends
                      !should also remain. 
                      if(can_remain) then

                         !as the neighbor buffer layers should remain, they
                         !are added to the list of buffer layers that should
                         !not be tested
                         call sublayer_ptr%get_neighbor1_id(neighbor1_id)
                         call sublayer_ptr%get_neighbor2_id(neighbor2_id)

                         call this%get_nbf_layers_sharing_grdpts_with(
     $                        sublayer_ptr,
     $                        not_to_be_tested_sublayers_list(neighbor1_id),
     $                        not_to_be_tested_sublayers_list(neighbor2_id), k)
                         
                      else

                         !if the buffer layer should be removed, we should
                         !check if it has neighbors. If it is dependent on
                         !other buffer layers, then this buffer layer is
                         !designated as a buffer layer which should be 
                         !checked again once all buffer layers are checked.
                         !In this way it will be clear if its neighbors can
                         !also be removed and so this buffer layer. If it
                         !is dependent on other neighboring buffer layers,
                         !the buffer layer can be removed immediately
                         neighbor_dependent = this%bf_layer_depends_on_neighbors(
     $                        sublayer_ptr, k)

                         if(neighbor_dependent) then
                            call to_be_rechecked_sublayers_list(k)%add_ele(
     $                           sublayer_ptr)
                         else
                            call this%remove_sublayer(sublayer_ptr)
                         end if

                      end if                      

                   !if it should not be tested, the buffer layer should
                   !remain
                   else
                      call sublayer_ptr%set_remain_status(.true.)
                   end if

                   sublayer_ptr => sublayer_next
                   
                end do

             end if

          end do          


          !2) the buffer layers whose removal needed to be
          !   confimed are investigated to see if their
          !   neighbors are also removed

          !loop over the cardinal coordinates
          do k=1,4

             !loop over the sublayers to be re-checked
             !for each cardinal coordinate
             nb_sublayers_rechecked = to_be_rechecked_sublayers_list(k)%get_nb_ele()

             do l=1, nb_sublayers_rechecked

                sublayer_ptr => to_be_rechecked_sublayers_list(k)%get_ele(l)

                !if none of the neihgboring buffer layers remain, the
                !current buffer layer should be removed
                if(.not.(this%does_a_neighbor_remains(sublayer_ptr,k))) then

                   call this%remove_sublayer(sublayer_ptr,k)

                end if

             end do

             call not_to_be_tested_sublayers_list(k)%destroy()
             call to_be_rechecked_sublayers_list(k)%destroy()

          end do          

        end subroutine update_bf_layers_with_detector_dcr


      end module bf_interface_dcr_class

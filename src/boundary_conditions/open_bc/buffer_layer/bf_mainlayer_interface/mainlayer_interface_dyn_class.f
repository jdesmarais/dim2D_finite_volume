      !> @file
      !> mainlayer_interface_sync enhanced with procedures enabling
      !> to check the dynamic nature of the mainlayer interfaces
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> mainlayer_interface_sync enhanced with procedures enabling
      !> to check the dynamic nature of the buffer layer interfaces
      !
      !> @date
      ! 09_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mainlayer_interface_dyn_class

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use bf_layer_exchange_module, only : 
     $       update_alignment_for_exchanges

        use bf_sublayer_class, only :
     $       bf_sublayer

        use mainlayer_interface_sync_class, only :
     $       mainlayer_interface_sync

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W,
     $       
     $       NW_interface_type, NE_interface_type,
     $       SW_interface_type, SE_interface_type

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       bc_size

        use parameters_kind, only :
     $       ikind


        implicit none

        private
        public :: mainlayer_interface_dyn


        !>@class mainlayer_interface_dyn
        !> mainlayer_interface_sync enhanced with procedures
        !> enabling to check the dynamic nature of the
        !> mainlayer interfaces
        !
        !>@param update_alignment_and_sync_properties
        !> check whether the alignment should be modified to
        !> set the buffer layer at the interface between mainlayers
        !
        !>@param uniformize_alignment_for_mainlayer_interface
        !> check whether the alignment should be modified to fit
        !> its neighboring buffer layers
        !------------------------------------------------------------
        type, extends(mainlayer_interface_sync) :: mainlayer_interface_dyn

          contains

          procedure, pass :: update_alignment_and_sync_properties
          procedure, pass :: uniformize_alignment_for_mainlayer_interface

        end type mainlayer_interface_dyn


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the alignment and the synchronization properties
        !> of the parameters for the buffer layer
        !
        !> @date
        !> 09_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating links to buffer layers at the edge
        !> between different main layers
        !
        !> @param bf_mainlayer_id
        !> integer identifying the main layer
        !
        !> @param bf_alignment
        !> respective position of the buffer layer compared to
        !> the interior domain
        !
        !> @param can_exchange_with_neighbor1
        !> logical determining whether the buffer layer is at the
        !> interface between main layers (South interface)
        !
        !> @param mainlayer_interface_type1
        !> integer identifying where the buffer layer is sharing
        !> grid points (SE or SW interface)
        !
        !> @param can_exchange_with_neighbor2
        !> logical determining whether the buffer layer is at the
        !> interface between main layers (North interface)
        !
        !> @param mainlayer_interface_type2
        !> integer identifying where the buffer layer is sharing
        !> grid points (NE or NW interface)
        !------------------------------------------------------------
        subroutine update_alignment_and_sync_properties(
     $       this,
     $       bf_mainlayer_id,
     $       bf_alignment,
     $       can_exchange_with_neighbor1,
     $       mainlayer_interface_type1,
     $       can_exchange_with_neighbor2,
     $       mainlayer_interface_type2)

          implicit none

          class(mainlayer_interface_dyn), intent(in)    :: this
          integer                       , intent(in)    :: bf_mainlayer_id
          integer(ikind), dimension(2,2), intent(inout) :: bf_alignment
          logical                       , intent(out)   :: can_exchange_with_neighbor1
          integer                       , intent(out)   :: mainlayer_interface_type1
          logical                       , intent(out)   :: can_exchange_with_neighbor2
          integer                       , intent(out)   :: mainlayer_interface_type2


          type(bf_sublayer), pointer :: neighbor1
          type(bf_sublayer), pointer :: neighbor2


          mainlayer_interface_type1 = 0
          mainlayer_interface_type2 = 0


          !alignment verification
          !------------------------------------------------------------
          !verify that the alignment corresponds to the bf_mainlayer_id
          !------------------------------------------------------------
          select case(bf_mainlayer_id)
            case(N)
               if(bf_alignment(2,1).ne.align_N) then
                  print '(''mainlayer_interface_dyn_class'')'
                  print '(''update_alignment_and_sync_properties'')'
                  stop 'bf_alignment(2,1).ne.align_N for N buffer layer'
               end if

            case(S)
               if(bf_alignment(2,2).ne.align_S) then
                  print '(''mainlayer_interface_dyn_class'')'
                  print '(''update_alignment_and_sync_properties'')'
                  stop 'bf_alignment(2,2).ne.align_S for S buffer layer'
               end if

            case(E)
               if(
     $              (bf_alignment(1,1).ne.align_E).or.
     $              (bf_alignment(2,1).lt.(align_S+1)).or.
     $              (bf_alignment(2,2).gt.(align_N-1)))then
                  print '(''mainlayer_interface_dyn_class'')'
                  print '(''update_alignment_and_sync_properties'')'
                  print '(''bf_alignment(1,1).ne.align_E    : '',L1)', bf_alignment(1,1).ne.align_E
                  print '(''bf_alignment(2,1).lt.(align_S+1): '',L1)', bf_alignment(2,1).lt.(align_S+1)
                  print '(''bf_alignment(2,2).gt.(align_N-1): '',L1)', bf_alignment(2,2).gt.(align_N-1)
                  print '(''for E buffer layer'')'
                  stop ''
               end if

            case(W)
               if(
     $              (bf_alignment(1,2).ne.align_W).or.
     $              (bf_alignment(2,1).lt.(align_S+1)).or.
     $              (bf_alignment(2,2).gt.(align_N-1)))then
                  print '(''mainlayer_interface_dyn_class'')'
                  print '(''update_alignment_and_sync_properties'')'
                  print '(''bf_alignment(1,2).ne.align_W    : '',L1)', bf_alignment(1,2).ne.align_W
                  print '(''bf_alignment(2,1).lt.(align_S+1): '',L1)', bf_alignment(2,1).lt.(align_S+1)
                  print '(''bf_alignment(2,2).gt.(align_N-1): '',L1)', bf_alignment(2,2).gt.(align_N-1)
                  print '(''for W buffer layer'')'
                  stop ''
               end if

            case default
               call error_mainlayer_id(
     $              'mainlayer_interface_dyn_class',
     $              'update_alignment_and_sync_properties',
     $              bf_mainlayer_id)

          end select


          !alignment update
          !------------------------------------------------------------
          !if the buffer layer is almost at the edge of mainlayer
          !interfaces, the alignment is modified to fit the edge
          !------------------------------------------------------------
          call update_alignment_for_exchanges(
     $         bf_mainlayer_id,
     $         bf_alignment)


          !------------------------------------------------------------
          !if the buffer layer is sharing grid points with the edges
          !of the main layers, the mainlayer_interface_type where its
          !pointer should be saved is determined
          !------------------------------------------------------------
          select case(bf_mainlayer_id)

            case(N,S)

               can_exchange_with_neighbor1 =
     $              (bf_alignment(1,1)-bc_size.le.(align_W+bc_size))

               if(can_exchange_with_neighbor1) then
                  if(bf_mainlayer_id.eq.N) then
                     mainlayer_interface_type1 = NW_interface_type
                  else
                     mainlayer_interface_type1 = SW_interface_type
                  end if
               end if

               can_exchange_with_neighbor2 =
     $              (bf_alignment(1,2)+bc_size.ge.(align_E-bc_size))

               if(can_exchange_with_neighbor2) then
                  if(bf_mainlayer_id.eq.N) then
                     mainlayer_interface_type2 = NE_interface_type
                  else
                     mainlayer_interface_type2 = SE_interface_type
                  end if
               end if


            case(E,W)

               can_exchange_with_neighbor1 =
     $              (bf_alignment(2,1)-bc_size.le.(align_S+bc_size))

               if(can_exchange_with_neighbor1) then
                  if(bf_mainlayer_id.eq.E) then
                     mainlayer_interface_type1 = SE_interface_type
                  else
                     mainlayer_interface_type1 = SW_interface_type
                  end if
               end if

               can_exchange_with_neighbor2 =
     $              (bf_alignment(2,2)+bc_size.ge.(align_N-bc_size))

               if(can_exchange_with_neighbor2) then
                  if(bf_mainlayer_id.eq.E) then
                     mainlayer_interface_type2 = NE_interface_type
                  else
                     mainlayer_interface_type2 = NW_interface_type
                  end if
               end if


            case default
               call error_mainlayer_id(
     $              'mainlayer_interface_dyn_class',
     $              'update_alignment_and_sync_properties',
     $              bf_mainlayer_id)


          end select

          
          !------------------------------------------------------------
          !if the buffer layer is sharing grid points with the edges
          !of the main layers, its alignment along the x-direction
          !is updated to fit its neighboring buffer layers
          !------------------------------------------------------------
          select case(bf_mainlayer_id)
            case(N,S)
               if(can_exchange_with_neighbor1) then
                  neighbor1 => this%get_neighbor_sublayer_ptr(
     $                 bf_mainlayer_id,
     $                 1)
                  
                  if(associated(neighbor1)) then
                     bf_alignment(1,1) = min(
     $                    bf_alignment(1,1),
     $                    neighbor1%get_alignment(1,1))
                  end if
                  
               end if

               if(can_exchange_with_neighbor2) then
                  neighbor2 => this%get_neighbor_sublayer_ptr(
     $                 bf_mainlayer_id,
     $                 2)
                  
                  if(associated(neighbor2)) then
                     bf_alignment(1,2) = max(
     $                    bf_alignment(1,2),
     $                    neighbor2%get_alignment(1,2))
                  end if

               end if

            case(W)
               if(can_exchange_with_neighbor1) then
                  neighbor1 => this%get_neighbor_sublayer_ptr(
     $                 bf_mainlayer_id,
     $                 1)
                  
                  if(associated(neighbor1)) then
                     bf_alignment(1,1) = min(
     $                    bf_alignment(1,1),
     $                    neighbor1%get_alignment(1,1))
                  end if
                  
               end if

               if(can_exchange_with_neighbor2) then
                  neighbor2 => this%get_neighbor_sublayer_ptr(
     $                 bf_mainlayer_id,
     $                 2)
                  
                  if(associated(neighbor2)) then
                     bf_alignment(1,1) = min(
     $                    bf_alignment(1,1),
     $                    neighbor2%get_alignment(1,1))
                  end if

               end if

            case(E)
               if(can_exchange_with_neighbor1) then
                  neighbor1 => this%get_neighbor_sublayer_ptr(
     $                 bf_mainlayer_id,
     $                 1)
                  
                  if(associated(neighbor1)) then
                     bf_alignment(1,2) = max(
     $                    bf_alignment(1,2),
     $                    neighbor1%get_alignment(1,2))
                  end if
                  
               end if

               if(can_exchange_with_neighbor2) then
                  neighbor2 => this%get_neighbor_sublayer_ptr(
     $                 bf_mainlayer_id,
     $                 2)
                  
                  if(associated(neighbor2)) then
                     bf_alignment(1,2) = max(
     $                    bf_alignment(1,2),
     $                    neighbor2%get_alignment(1,2))
                  end if

               end if

            case default
               call error_mainlayer_id(
     $              'mainlayer_interface_dyn_class',
     $              'update_alignment_and_sync_properties',
     $              bf_mainlayer_id)

          end select

        end subroutine update_alignment_and_sync_properties


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the alignment and the synchronization properties
        !> of the parameters for the buffer layer
        !
        !> @date
        !> 09_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating links to buffer layers at the edge
        !> between different main layers
        !
        !> @param bf_sublayer_ptr
        !> sublayer whose alignment may be updated to match
        !> the alignment of the mainlayer interfaces
        !
        !> @param bf_mainlayer_id
        !> integer identifying the main layer
        !------------------------------------------------------------
        function uniformize_alignment_for_mainlayer_interface(
     $     this,
     $     bf_sublayer_ptr,
     $     should_be_reallocated,
     $     mainlayer_id)
     $     result(bf_alignment)

          implicit none

          class(mainlayer_interface_dyn), intent(in) :: this
          type(bf_sublayer), pointer    , intent(in) :: bf_sublayer_ptr
          logical                       , intent(out):: should_be_reallocated
          integer, optional             , intent(in) :: mainlayer_id
          integer(ikind), dimension(2,2)             :: bf_alignment

          integer :: bf_mainlayer_id
          integer(ikind), dimension(2,2) :: bf_alignment_ini

          
          if(present(mainlayer_id)) then
             bf_mainlayer_id = mainlayer_id
          else
             bf_mainlayer_id = bf_sublayer_ptr%get_localization()
          end if


          bf_alignment_ini = bf_sublayer_ptr%get_alignment_tab()
          bf_alignment     = bf_alignment_ini


          select case(bf_mainlayer_id)
            case(N)

               if(bf_sublayer_ptr%can_exchange_with_neighbor1()) then

                  bf_alignment(1,1) = min(
     $                 bf_alignment(1,1),
     $                 this%NW_interface_W_ptr%get_alignment(1,1))
                  !
                  !   _______
                  !  |___N___|
                  !   ___   _____
                  !  |   | |     |
                  !  | W | | int |
                  !  |___| |_____|
                  !   _____
                  !  |__S__|
                  !
                  if(associated(this%SW_interface_W_ptr)) then
                     if(associated(
     $                    this%NW_interface_W_ptr,
     $                    this%SW_interface_W_ptr)) then
                        if(associated(this%SW_interface_S_ptr)) then
                           bf_alignment(1,1) = min(
     $                          bf_alignment(1,1),
     $                          this%SW_interface_S_ptr%get_alignment(1,1))
                        end if
                     end if
                  end if

               end if

               if(bf_sublayer_ptr%can_exchange_with_neighbor2()) then

                  bf_alignment(1,2) = max(
     $                 bf_alignment(1,2),
     $                 this%NE_interface_E_ptr%get_alignment(1,2))
                  !
                  !       _______
                  !      |___N___|
                  !   _____   ___  
                  !  |     | |   |
                  !  | int | | E |
                  !  |_____| |___|
                  !         _____
                  !        |__S__|
                  !
                  if(associated(this%SE_interface_E_ptr)) then
                     if(associated(
     $                    this%NE_interface_E_ptr,
     $                    this%SE_interface_E_ptr)) then
                        if(associated(this%SE_interface_S_ptr)) then
                           bf_alignment(1,2) = max(
     $                          bf_alignment(1,2),
     $                          this%SE_interface_S_ptr%get_alignment(1,2))
                        end if
                     end if
                  end if

               end if
               

            case(S)

               if(bf_sublayer_ptr%can_exchange_with_neighbor1()) then

                  bf_alignment(1,1) = min(
     $                 bf_alignment(1,1),
     $                 this%SW_interface_W_ptr%get_alignment(1,1))
                  !
                  !   _______
                  !  |___N___|
                  !   ___   _____
                  !  |   | |     |
                  !  | W | | int |
                  !  |___| |_____|
                  !   _____
                  !  |__S__|
                  !
                  if(associated(this%NW_interface_W_ptr)) then
                     if(associated(
     $                    this%NW_interface_W_ptr,
     $                    this%SW_interface_W_ptr)) then
                        if(associated(this%NW_interface_N_ptr)) then
                           bf_alignment(1,1) = min(
     $                          bf_alignment(1,1),
     $                          this%NW_interface_N_ptr%get_alignment(1,1))
                        end if
                     end if
                  end if

               end if

               if(bf_sublayer_ptr%can_exchange_with_neighbor2()) then

                  bf_alignment(1,2) = max(
     $                 bf_alignment(1,2),
     $                 this%SE_interface_E_ptr%get_alignment(1,2))
                  !
                  !       _______
                  !      |___N___|
                  !   _____   ___  
                  !  |     | |   |
                  !  | int | | E |
                  !  |_____| |___|
                  !         _____
                  !        |__S__|
                  !
                  if(associated(this%NE_interface_E_ptr)) then
                     if(associated(
     $                    this%NE_interface_E_ptr,
     $                    this%SE_interface_E_ptr)) then
                        if(associated(this%NE_interface_N_ptr)) then
                           bf_alignment(1,2) = max(
     $                          bf_alignment(1,2),
     $                          this%NE_interface_N_ptr%get_alignment(1,2))
                        end if
                     end if
                  end if

               end if

            case(E)

               !
               !       _______
               !      |___N___|
               !   _____   ___  
               !  |     | |   |
               !  | int | | E |
               !  |_____| |___|
               !         _____
               !        |__S__|
               !
               if(bf_sublayer_ptr%can_exchange_with_neighbor1()) then
                  bf_alignment(1,2) = max(
     $                 bf_alignment(1,2),
     $                 this%SE_interface_S_ptr%get_alignment(1,2))

               end if

               if(bf_sublayer_ptr%can_exchange_with_neighbor2()) then
                  bf_alignment(1,2) = max(
     $                 bf_alignment(1,2),
     $                 this%NE_interface_N_ptr%get_alignment(1,2))
               end if

            case(W)

               !
               !   _______
               !  |___N___|
               !   ___   _____
               !  |   | |     |
               !  | W | | int |
               !  |___| |_____|
               !   _____
               !  |__S__|
               !
               if(bf_sublayer_ptr%can_exchange_with_neighbor1()) then
                  bf_alignment(1,1) = min(
     $                 bf_alignment(1,1),
     $                 this%SW_interface_S_ptr%get_alignment(1,1))

               end if

               if(bf_sublayer_ptr%can_exchange_with_neighbor2()) then
                  bf_alignment(1,1) = min(
     $                 bf_alignment(1,1),
     $                 this%NW_interface_N_ptr%get_alignment(1,1))
               end if

            case default
               call error_mainlayer_id(
     $              'mainlayer_interface_dyn_class',
     $              'uniformize_alignment_for_mainlayer_interface',
     $              bf_mainlayer_id)

          end select


          should_be_reallocated = 
     $         (bf_alignment_ini(1,1).ne.bf_alignment(1,1)).or.
     $         (bf_alignment_ini(2,1).ne.bf_alignment(2,1)).or.
     $         (bf_alignment_ini(1,2).ne.bf_alignment(1,2)).or.
     $         (bf_alignment_ini(2,2).ne.bf_alignment(2,2))

        end function uniformize_alignment_for_mainlayer_interface

      end module mainlayer_interface_dyn_class

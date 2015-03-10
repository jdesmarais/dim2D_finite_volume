      !> @file
      !> bf_interface_sync augmented with allocation/reallocation/
      !> merge/remove functions for the buffer layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_interface_sync augmented with allocation/reallocation/
      !> merge/remove functions for the buffer layers
      !
      !> @date
      ! 10_03_2015 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interface_dyn_class

        use bf_interface_sync_class, only :
     $       bf_interface_sync

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: bf_interface_dyn
        

        !>@class bf_interface_dyn
        !> bf_interface_sync augmented with allocation/reallocation/
        !> merge/remove functions for the buffer layers
        !
        !>@param allocate_sublayer
        !> allocate a bf_sublayer and insert it in the corresponding
        !> buffer mainlayer
        !
        !>@param reallocate_sublayer
        !> reallocate a buffer sublayer and check whether the neighboring
        !> buffer layer changed and so if the neighboring links should be
        !> updated
        !
        !>@param merge_sublayers
        !> merge the content of two sublayers and update
        !> the links to the border buffer layers
        !
        !>@param remove_sublayer
        !> remove a sublayer from the buffer main layers
        !
        !>@param uniformize_mainlayer_interfaces
        !> verify that the buffer layer at the interfaces between
        !> main layers have the same alignment
        !---------------------------------------------------------------
        type, extends(bf_interface_sync) :: bf_interface_dyn

          contains

          procedure, pass :: allocate_sublayer
          procedure, pass :: reallocate_sublayer
          procedure, pass :: merge_sublayers
          procedure, pass :: remove_sublayer

          procedure, pass :: uniformize_mainlayer_interfaces

        end type bf_interface_dyn

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate a bf_sublayer and insert it in the corresponding
        !> buffer mainlayer
        !
        !> @date
        !> 10_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        !> the data between them
        !
        !>@param mainlayer_id
        !> cardinal coordinate for the new bf_sublayer
        !
        !>@param nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@param alignment
        !> table of integers characterizing the
        !> correspondance between the interior grid points
        !> and the buffer layer elements
        !
        !>@return added_sublayer
        !> reference to the newly allocated bf_sublayer
        !--------------------------------------------------------------
        function allocate_sublayer(
     $     this,
     $     bf_mainlayer_id,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     bf_alignment)
     $     result(added_sublayer)
        
          class(bf_interface_dyn)         , intent(inout) :: this
          integer                         , intent(in)    :: bf_mainlayer_id
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          integer    , dimension(2,2)     , intent(inout) :: bf_alignment
          type(bf_sublayer), pointer                      :: added_sublayer


          logical :: can_exchange_with_neighbor1
          integer :: mainlayer_interface_type1
          logical :: can_exchange_with_neighbor2
          integer :: mainlayer_interface_type2


          !0) check the mainlayer id
          if(.not.(
     $         (bf_mainlayer_id.eq.N).or.
     $         (bf_mainlayer_id.eq.S).or.
     $         (bf_mainlayer_id.eq.E).or.
     $         (bf_mainlayer_id.eq.W))) then
             
             call error_mainlayer_id(
     $            'bf_interface_dyn_class',
     $            'allocate_sublayer',
     $            bf_mainlayer_id)
             
          end if


          !1) check if the new sublayer is at the interface between several
          !   main layers and update the alignment in case it is just one grid
          !   point away from a corner: exchanges are made easier
          call this%mainlayer_interfaces%update_alignment_and_sync_properties(
     $         bf_mainlayer_id,
     $         bf_alignment,
     $         can_exchange_with_neighbor1,
     $         mainlayer_interface_type1,
     $         can_exchange_with_neighbor2,
     $         mainlayer_interface_type2)


          !2) check if the mainlayer corresponding to the cardinal point is
          !   indeed allocated: if the memory space is not allocated, the space
          !   in memory is first allocated, the pointer identifying the mainlayer
          !   is initialized and the main layer itself is initialized
          if(.not.this%mainlayer_pointers(bf_mainlayer_id)%associated_ptr()) then
             call this%mainlayer_pointers(bf_mainlayer_id)%ini_mainlayer(bf_mainlayer_id)
          end if


          !3) now that we are sure that space is allocated for the main layer,
          !   the sublayer can be integrated to the mainlayer and the buffer
          !   layer can be initialized using the interior data and the
          !   alignment
          added_sublayer => this%mainlayer_pointers(bf_mainlayer_id)%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes, 
     $         bf_alignment)


          !4) if the sublayer newly allocated is indeed a buffer layer at the
          !   interface between main layers and of type 1, the neighboring
          !   buffer layers should know that they will exchange data with this
          !   buffer layer.          
          !   The same is true for an interface buffer layer of type 2      
          if(can_exchange_with_neighbor1) then

             call this%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $            mainlayer_interface_type1, added_sublayer)

          end if

          if(can_exchange_with_neighbor2) then

             call this%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $            mainlayer_interface_type2, added_sublayer)

          end if

          
          !5) Morover, the new buffer layer has been allocated without considering
          !   the other buffer layers already allocated which could share grid
          !   points with this buffer layer. These grid points are now copy from
          !   the other buffer layers to this buffer layer
          call this%mainlayer_interfaces%copy_from_mainlayer_neighbors(added_sublayer)


       end function allocate_sublayer


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> reallocate a buffer sublayer and check whether the neighboring
       !> buffer layer changed and so if the neighboring links should be
       !> updated
       !
       !> @date
       !> 10_03_2015 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param bf_sublayer_r
       !> reference to the bf_sublayer to be reallocated
       !
       !>@param interior_x_map
       !> array encapsulating the coordinates along the x-axis for
       !> the interior computational domain
       !
       !>@param interior_y_map
       !> array encapsulating the coordinates along the y-axis for
       !> the interior computational domain
       !
       !>@param interior_nodes
       !> table encapsulating the data of the grid points of the
       !> interior domain
       !
       !>@param bf_alignment
       !> table of integers localizing the relative position of
       !> the buffer layer compared to the interior domain
       !--------------------------------------------------------------
       subroutine reallocate_sublayer(
     $     this,
     $     bf_sublayer_r,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     bf_alignment)

         implicit none

         class(bf_interface_dyn)         , intent(inout) :: this
         type(bf_sublayer), pointer      , intent(inout) :: bf_sublayer_r
         real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
         real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
         real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
         integer    , dimension(2,2)     , intent(inout) :: bf_alignment


         integer :: bf_mainlayer_id
         logical :: can_exchange_with_neighbor1
         integer :: mainlayer_interface_type1
         logical :: can_exchange_with_neighbor2
         integer :: mainlayer_interface_type2


         bf_mainlayer_id = bf_sublayer_r%get_localization()


         !1) check if the reallocated sublayer is at the interface between
         !   several main layers and update the alignment in case it is just
         !   one grid point away from a corner: exchanges are made easier
         call this%mainlayer_interfaces%update_alignment_and_sync_properties(
     $         bf_mainlayer_id,
     $         bf_alignment,
     $         can_exchange_with_neighbor1,
     $         mainlayer_interface_type1,
     $         can_exchange_with_neighbor2,
     $         mainlayer_interface_type2)

         
         !2) reallocate the buffer sublayer
         call bf_sublayer_r%reallocate_bf_layer(
     $        interior_x_map,
     $        interior_y_map,
     $        interior_nodes,
     $        bf_alignment)


         !3) check if the links to the neighbor1 should be updated
         ! if the bf_sublayer_r was exchanging with neighbor1 before
         if(can_exchange_with_neighbor1) then
            
            call this%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $           mainlayer_interface_type1,
     $           bf_sublayer_r)

         end if

         if(can_exchange_with_neighbor2) then
            
            call this%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $           mainlayer_interface_type2,
     $           bf_sublayer_r)

         end if


         !4) update the grid points new allocated using the
         !   neighboring sublayers
         call this%mainlayer_interfaces%copy_from_mainlayer_neighbors(
     $        bf_sublayer_r)

       end subroutine reallocate_sublayer
       

       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> merge the content of two sublayers and update
       !> the links to the border buffer layers
       !
       !> @date
       !> 10_03_2015 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param bf_sublayer1
       !> reference to the first bf_sublayer merged
       !
       !>@param bf_sublayer2
       !> reference to the second bf_sublayer merged
       !
       !>@param nodes
       !> table encapsulating the data of the grid points of the
       !> interior domain
       !
       !>@param alignment
       !> table of integers characterizing the
       !> correspondance between the interior grid points
       !> and the buffer layer
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

          class(bf_interface_dyn)                    , intent(inout) :: this
          type(bf_sublayer), pointer                 , intent(inout) :: bf_sublayer1
          type(bf_sublayer), pointer                 , intent(inout) :: bf_sublayer2
          real(rkind)      , dimension(nx)           , intent(in)    :: interior_x_map
          real(rkind)      , dimension(ny)           , intent(in)    :: interior_y_map
          real(rkind)      , dimension(nx,ny,ne)     , intent(in)    :: interior_nodes
          integer(ikind)   , dimension(2,2), optional, intent(inout) :: alignment
          type(bf_sublayer), pointer                                 :: merged_sublayer

          integer                        :: bf_mainlayer_id
          integer(ikind), dimension(2,2) :: bf_alignment
          logical                        :: can_exchange_with_neighbor1
          integer                        :: mainlayer_interface_type1
          logical                        :: can_exchange_with_neighbor2
          integer                        :: mainlayer_interface_type2


          bf_mainlayer_id = bf_sublayer1%get_localization()

          if(present(alignment)) then
             bf_alignment = alignment
          else
             bf_alignment = bf_sublayer1%get_alignment_tab()
             bf_alignment(1,1) = min(bf_alignment(1,1),bf_sublayer2%get_alignment(1,1))
             bf_alignment(1,2) = max(bf_alignment(1,2),bf_sublayer2%get_alignment(1,2))
             bf_alignment(2,1) = min(bf_alignment(2,1),bf_sublayer2%get_alignment(2,1))
             bf_alignment(2,2) = max(bf_alignment(2,2),bf_sublayer2%get_alignment(2,2))
          end if


          !1) check if the merged sublayer is at the interface between
          !   several main layers and update the alignment in case it is just
          !   one grid point away from a corner: exchanges are made easier
          call this%mainlayer_interfaces%update_alignment_and_sync_properties(
     $          bf_mainlayer_id,
     $          bf_alignment,
     $          can_exchange_with_neighbor1,
     $          mainlayer_interface_type1,
     $          can_exchange_with_neighbor2,
     $          mainlayer_interface_type2)


          !2) merge the content of the sublayers
          merged_sublayer => this%mainlayer_pointers(bf_mainlayer_id)%merge_sublayers(
     $         bf_sublayer1,
     $         bf_sublayer2,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)


          !3) check if the links to the neighbor1 should be updated
          !   if the bf_sublayer_r was exchanging with neighbor1 before
          if(can_exchange_with_neighbor1) then
             
             call this%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $            mainlayer_interface_type1,
     $            merged_sublayer)
          
          end if
          
          if(can_exchange_with_neighbor2) then
             
             call this%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $            mainlayer_interface_type2,
     $            merged_sublayer)
          
          end if


          !4) update the grid points that have been newly allocated
          !   with grid points from the neighboring sublayers
          call this%mainlayer_interfaces%copy_from_mainlayer_neighbors(merged_sublayer)

       end function merge_sublayers


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> remove a sublayer from the buffer main layers
       !
       !> @date
       !> 10_03_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param sublayer_ptr
       !> reference to the bf_sublayer removed
       !
       !>@param bf_mainlayer_id
       !> cardinal coordinate of the buffer layer removed
       !--------------------------------------------------------------
       subroutine remove_sublayer(
     $     this,
     $     sublayer_ptr,
     $     mainlayer_id)

         implicit none

         class(bf_interface_dyn)   , intent(inout) :: this
         type(bf_sublayer), pointer, intent(inout) :: sublayer_ptr
         integer, optional         , intent(in)    :: mainlayer_id

         integer                        :: bf_mainlayer_id
         integer(ikind), dimension(2,2) :: bf_alignment

         logical :: can_exchange_with_neighbor1
         integer :: mainlayer_interface_type1
         logical :: can_exchange_with_neighbor2
         integer :: mainlayer_interface_type2


         !> identify the mainlayer to which the sublayer belongs
         if(present(mainlayer_id)) then
            bf_mainlayer_id = mainlayer_id
         else
            bf_mainlayer_id = sublayer_ptr%get_localization()
         end if

         bf_alignment = sublayer_ptr%get_alignment_tab()

         
         !> determine whether the sublayer removed was sharing
         !> grid points with the main layer interfaces
         call this%mainlayer_interfaces%update_alignment_and_sync_properties(
     $         bf_mainlayer_id,
     $         bf_alignment,
     $         can_exchange_with_neighbor1,
     $         mainlayer_interface_type1,
     $         can_exchange_with_neighbor2,
     $         mainlayer_interface_type2)


         !> remove the sublayer from the table identifying the
         !> neighboring buffer layers
         if(can_exchange_with_neighbor1) then

            call this%mainlayer_interfaces%remove_mainlayer_interface_bf_layer(
     $           mainlayer_interface_type1,
     $           sublayer_ptr)

         end if

         if(can_exchange_with_neighbor2) then

            call this%mainlayer_interfaces%remove_mainlayer_interface_bf_layer(
     $           mainlayer_interface_type2,
     $           sublayer_ptr)

         end if

         !> remove the sublayer from the main layer
         call this%mainlayer_pointers(bf_mainlayer_id)%remove_sublayer(sublayer_ptr)

       end subroutine remove_sublayer


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> uniformize the alignments of the buffer layers
       !> at the interface between the main layers
       !
       !> @date
       !> 10_03_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !--------------------------------------------------------------
       subroutine uniformize_mainlayer_interfaces(
     $     this,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes)

         implicit none

         class(bf_interface_dyn)         , intent(inout) :: this
         real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
         real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
         real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes         


         type(bf_sublayer), pointer     :: bf_sublayer_ptr
         integer(ikind), dimension(2,2) :: bf_alignment
         logical                        :: should_be_reallocated

         integer :: k
         integer :: l


         !loop over the main layers
         do k=1,4
            !loop over the neighbors
            do l=1,2
               
               !get the sublayer at the interface between the main layers
               bf_sublayer_ptr =>
     $              this%mainlayer_interfaces%get_neighbor_sublayer_ptr(k,l)

               !if the sublayer exists, verify that its alignment
               !fit the other main layers
               if(associated(bf_sublayer_ptr)) then

                  bf_alignment = this%mainlayer_interfaces%uniformize_alignment_for_mainlayer_interface(
     $                 bf_sublayer_ptr,
     $                 should_be_reallocated)

               !if its alignment does not fit the other main layers
               !the buffer layer is reallocated
                  if(should_be_reallocated) then
                     call this%reallocate_sublayer(
     $                    bf_sublayer_ptr,
     $                    interior_x_map,
     $                    interior_y_map,
     $                    interior_nodes,
     $                    bf_alignment)

                  end if

               end if

            end do
         end do


       end subroutine uniformize_mainlayer_interfaces

      end module bf_interface_dyn_class

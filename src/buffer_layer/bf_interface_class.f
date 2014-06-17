      module bf_interface_class

        use bf_layer_errors_module    , only : error_mainlayer_id,
     $                                         error_incompatible_neighbor
        use bf_sublayer_class         , only : bf_sublayer
        use bf_mainlayer_class        , only : bf_mainlayer
        use bf_mainlayer_pointer_class, only : bf_mainlayer_pointer
        use nbf_interface_class       , only : nbf_interface

        use parameters_bf_layer       , only : align_N, align_S,
     $                                         align_E, align_W
        use parameters_constant       , only : N,S,E,W,
     $                                         x_direction, y_direction,
     $                                         interior
        use parameters_input          , only : nx,ny,ne,bc_size,debug
        use parameters_kind           , only : ikind, rkind
        use sbf_list_class            , only : sbf_list


        implicit none

        private
        public :: bf_interface
        

        !> @class bf_interface
        !> class encapsulating the bf_layer/interior interface
        !> object
        !
        !> @param mainlayer_pointers
        !> table with reference to the buffer main layers
        !
        !> @param border_interface
        !> object with references to the sublayers at the interface
        !> between the main layers and ways to exchange data between
        !> these boundary sublayers
        !
        !> @param ini
        !> initialize the interface by nullifying all the 
        !> references to the buffer main layers
        !
        !> @param get_mainlayer
        !> get the reference to the main layer corresponding
        !> to the cardinal coordinate passed
        !
        !> @param add_sublayer
        !> add a new sublayer to the main layer corresponding
        !> to the cardinal coordinate passed
        !---------------------------------------------------------------
        type :: bf_interface

          type(bf_mainlayer_pointer), dimension(4), private :: mainlayer_pointers
          type(nbf_interface)                     , private :: border_interface

          contains

          procedure, pass :: ini
          procedure, pass :: get_mainlayer

          procedure, pass :: allocate_sublayer
          procedure, pass :: reallocate_sublayer
          procedure, pass :: merge_sublayers
          procedure, pass :: remove_sublayer
          procedure, pass :: update_grdpts_after_increase

          procedure, nopass :: get_mainlayer_id
          procedure, pass   :: get_sublayer
          procedure, pass   :: get_nodes

          procedure, pass, private :: update_grdpts_from_neighbors
          procedure, pass, private :: update_neighbor_grdpts
          
          procedure, nopass, private :: shares_with_neighbor1
          procedure, nopass, private :: shares_with_neighbor2

          procedure, pass :: get_nbf_layers_sharing_grdpts_with
          procedure, pass :: bf_layer_depends_on_neighbors
          procedure, pass :: does_a_neighbor_remains

          procedure, pass :: print_binary

        end type bf_interface

        contains


        !< nullify all the pointers to the main layers
        subroutine ini(this)

          implicit none

          class(bf_interface), intent(inout) :: this

          integer :: i
          
          do i=1, size(this%mainlayer_pointers,1)
             call this%mainlayer_pointers(i)%ini()             
          end do

          call this%border_interface%ini()

        end subroutine ini


       !< get main layer corresponding to the cardinal point
       function get_mainlayer(this, mainlayer_id)
       
         implicit none
   
         class(bf_interface), intent(in) :: this
         integer            , intent(in) :: mainlayer_id
         type(bf_mainlayer) , pointer    :: get_mainlayer

         if(debug) then
            if((mainlayer_id.lt.1).or.(mainlayer_id.gt.4)) then
               call error_mainlayer_id(
     $              'bf_interface_class',
     $              'get_mainlayer',
     $              mainlayer_id)
            end if
         end if

         if(this%mainlayer_pointers(mainlayer_id)%associated_ptr()) then
            get_mainlayer => this%mainlayer_pointers(mainlayer_id)%get_ptr()
         else
            nullify(get_mainlayer)
         end if

       end function get_mainlayer


       !< check if the alignment of the future sublayer makes it a potential
       !> buffer layer at the interface between main layers and if so update
       !> its alignment if it is shifted by one grid point from the border
       !> as it makes exchanges easier later
       !> the sublayer is here tested as a neighbor buffer layer of type 1
       !> this function echoes bf_layer_class%shares_grdpts_with_neighbor1
       function shares_with_neighbor1(mainlayer_id, bf_final_alignment)
     $     result(share_grdpts)

         implicit none

         integer                       , intent(in)    :: mainlayer_id
         integer(ikind), dimension(2,2), intent(inout) :: bf_final_alignment
         logical                                       :: share_grdpts


         select case(mainlayer_id)
           case(N,S)
              share_grdpts = bf_final_alignment(1,1).le.(align_W+bc_size)
           case(E,W)
              share_grdpts = bf_final_alignment(2,1).le.(align_S+bc_size)
              if(bf_final_alignment(2,1).eq.(align_S+bc_size)) then
                 bf_final_alignment(2,1)=align_S+1
              end if
           case default
              call error_mainlayer_id(
     $             'bf_layer_class.f',
     $             'share_grdpts_with_neighbor1',
     $             mainlayer_id)
         end select

       end function shares_with_neighbor1


       !< check if the alignment of the future sublayer makes it a potential
       !> buffer layer at the interface between main layers and if so update
       !> its alignment if it is shifted by one grid point from the border
       !> as it makes exchanges easier later
       !> the sublayer is here tested as a neighbor buffer layer of type 2
       !> this function echoes bf_layer_class%shares_grdpts_with_neighbor2
       function shares_with_neighbor2(mainlayer_id, bf_final_alignment)
     $     result(share_grdpts)

         implicit none

         integer                       , intent(in)    :: mainlayer_id
         integer(ikind), dimension(2,2), intent(inout) :: bf_final_alignment
         logical                                       :: share_grdpts


         select case(mainlayer_id)
           case(N,S)
              share_grdpts = bf_final_alignment(1,2).ge.(align_E-bc_size)
           case(E,W)
              share_grdpts = bf_final_alignment(2,2).ge.(align_N-bc_size)
              if(bf_final_alignment(2,2).eq.(align_N-bc_size)) then
                 bf_final_alignment(2,2)=align_N-1
              end if
           case default
              call error_mainlayer_id(
     $             'bf_layer_class.f',
     $             'share_grdpts_with_neighbor1',
     $             mainlayer_id)
         end select

       end function shares_with_neighbor2


       !< allocate a buffer sublayer and insert it in the corresponding
       !> buffer mainlayer
       function allocate_sublayer(
     $     this,
     $     mainlayer_id,
     $     nodes,
     $     alignment)
     $     result(added_sublayer)
        
          class(bf_interface)             , intent(inout) :: this
          integer                         , intent(in)    :: mainlayer_id
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer, dimension(2,2)         , intent(inout) :: alignment

          type(bf_sublayer), pointer                      :: added_sublayer


          logical :: share_with_neighbor1
          logical :: share_with_neighbor2


          !0) debug : check the mainlayer id
          if(debug) then
             if((mainlayer_id.lt.1).or.(mainlayer_id.gt.4)) then
                call error_mainlayer_id(
     $              'bf_interface_class',
     $              'allocate_sublayer',
     $              mainlayer_id)
             end if
          end if

          !1) check if the new sublayer is at the interface between several
          !   main layers and update the alignment in case it is just one grid
          !   point away from a corner: exchanges are made easier
          share_with_neighbor1 = shares_with_neighbor1(mainlayer_id, alignment)
          share_with_neighbor2 = shares_with_neighbor2(mainlayer_id, alignment)


          !2) check if the mainlayer corresponding to the cardinal point is
          !   indeed allocated: if the memory space is not allocated, the space
          !   in memory is first allocated, the pointer identifying the mainlayer
          !   is initialized and the main layer itself is initialized
          if(.not.this%mainlayer_pointers(mainlayer_id)%associated_ptr()) then
             call this%mainlayer_pointers(mainlayer_id)%ini_mainlayer(mainlayer_id)
          end if


          !3) now that we are sure that space is allocated for the main layer,
          !   the sublayer can be integrated to the mainlayer and the buffer
          !   layer can be initialized using the nodes, alignment and neighbors
          !   arguments
          added_sublayer => this%mainlayer_pointers(mainlayer_id)%add_sublayer(
     $         nodes, alignment)
          call added_sublayer%set_neighbor1_share(share_with_neighbor1)
          call added_sublayer%set_neighbor2_share(share_with_neighbor2)

          !4) if the sublayer newly allocated is indeed a buffer layer at the
          !   interface between main layers and of type 1, the neighboring
          !   buffer layers should know that they will exchange data with this
          !   buffer layer.          
          !   The same is true for an interface buffer layer of type 2          
          if(share_with_neighbor1) then
             call this%border_interface%link_neighbor1_to_bf_sublayer(
     $            added_sublayer)
          end if
          if(share_with_neighbor2) then
             call this%border_interface%link_neighbor2_to_bf_sublayer(
     $            added_sublayer)
          end if
          
          !5) Morover, the new buffer layer has been allocated without considering
          !   the other buffer layers already allocated which could share grid
          !   points with this buffer layer. These grid points are now copy from
          !   the other buffer layers to this buffer layer
          call this%border_interface%update_grdpts_from_neighbors(added_sublayer)

       end function allocate_sublayer


       !< reallocate a buffer sublayer and check whether the neighboring buffer
       !> layer changed and so if the neighboring links should be udpated
       subroutine reallocate_sublayer(this, bf_sublayer_r, nodes, alignment)

         implicit none

         class(bf_interface)             , intent(inout) :: this
         type(bf_sublayer), pointer      , intent(inout) :: bf_sublayer_r
         real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
         integer    , dimension(2,2)     , intent(inout) :: alignment

         integer :: mainlayer_id
         logical :: share_with_neighbor1
         logical :: share_with_neighbor2

         mainlayer_id = bf_sublayer_r%get_localization()


         !1) check if the reallocated sublayer is at the interface between
         !   several main layers and update the alignment in case it is just
         !   one grid point away from a corner: exchanges are made easier
         share_with_neighbor1 = shares_with_neighbor1(mainlayer_id, alignment)
         share_with_neighbor2 = shares_with_neighbor2(mainlayer_id, alignment)

         
         !2) reallocate the buffer sublayer
         call bf_sublayer_r%reallocate_bf_layer(
     $        nodes, alignment)

         !3) check if the links to the neighbor1 should be updated
         ! if the bf_sublayer_r was exchanging with neighbor1 before
         if(bf_sublayer_r%can_exchange_with_neighbor1()) then
            
            !and the bf_sublayer_r can no longer exchange
            !the previous links should be removed
            if(.not.share_with_neighbor1) then
               call this%border_interface%remove_link_from_neighbor1_to_bf_sublayer(
     $              bf_sublayer_r)
            end if

         ! if the bf_sublayer_r was not exchanging with neighbor1 before
         else

            !and the bf_sublayer_r can now exchange, links should be added
            !to this bf_sublayer_r
            if(share_with_neighbor1) then
               call this%border_interface%link_neighbor1_to_bf_sublayer(
     $              bf_sublayer_r)
            end if

         end if

         !4) check if the links to the neighbor2 should be updated
         ! if the bf_sublayer_r was exchanging with neighbor2 before
         if(bf_sublayer_r%can_exchange_with_neighbor2()) then
            
            !and the bf_sublayer_r can no longer exchange
            !the previous links should be removed
            if(.not.share_with_neighbor2) then
               call this%border_interface%remove_link_from_neighbor2_to_bf_sublayer(
     $              bf_sublayer_r)
            end if

         ! if the bf_sublayer_r was not exchanging with neighbor2 before
         else

            !and the bf_sublayer_r can now exchange, links should be added
            !to this bf_sublayer_r
            if(share_with_neighbor2) then
               call this%border_interface%link_neighbor2_to_bf_sublayer(
     $              bf_sublayer_r)
            end if

         end if

         !5) update the status of the bf_sublayer_r for its neighbors
         call bf_sublayer_r%set_neighbor1_share(share_with_neighbor1)
         call bf_sublayer_r%set_neighbor2_share(share_with_neighbor2)

         !6) update the grid points new allocated using the neighboring sublayers
         call this%border_interface%update_grdpts_from_neighbors(bf_sublayer_r)
         
       end subroutine reallocate_sublayer
       

       !< merge sublayers : merge the content of the sublayers and update
       !> the links to the border buffer layers
       function merge_sublayers(
     $     this,
     $     bf_sublayer1, bf_sublayer2,
     $     nodes, alignment)
     $     result(merged_sublayer)

          implicit none

          class(bf_interface)                        , intent(inout) :: this
          type(bf_sublayer), pointer                 , intent(inout) :: bf_sublayer1
          type(bf_sublayer), pointer                 , intent(inout) :: bf_sublayer2
          real(rkind)      , dimension(nx,ny,ne)     , intent(in)    :: nodes
          integer(ikind)   , dimension(2,2), optional, intent(inout) :: alignment
          type(bf_sublayer), pointer                                 :: merged_sublayer

          integer :: mainlayer_id
          logical :: share_with_neighbor1
          logical :: share_with_neighbor2


          mainlayer_id = bf_sublayer1%get_localization()
   

          if(present(alignment)) then


             !1) check if the sublayer resulting from the merge is at the interface
             !   between several main layers and update the alignment in case it is just
             !   one grid point away from a corner: exchanges are made easier
             share_with_neighbor1 = shares_with_neighbor1(mainlayer_id, alignment)
             share_with_neighbor2 = shares_with_neighbor2(mainlayer_id, alignment)


             !2) check if the links to neighbor1 should be updated
             !if the bf_sublayer1 was already linked to neighbor1 sublayers
             if(bf_sublayer1%can_exchange_with_neighbor1()) then
                
                !if the bf_sublayer2 was also linked to neighbor1 sublayers
                !as it will be merged, these links should be removed
                if(bf_sublayer2%can_exchange_with_neighbor1()) then
                   call this%border_interface%remove_link_from_neighbor1_to_bf_sublayer(
     $                  bf_sublayer2)

                end if

             !if the bf_sublayer1 was not linked to neighbor1 sublayers
             else

                !and the bf_sublayer2 was linked to this neighbor, the links should
                !be updated from bf_sublayer2 to bf_sublayer1 and the status of the
                !neighbor1 for bf_sublayer1 should be updated
                if(bf_sublayer2%can_exchange_with_neighbor1()) then

                   call this%border_interface%update_link_from_neighbor1_to_bf_sublayer(
     $                  bf_sublayer1, bf_sublayer2)

                   call bf_sublayer1%set_neighbor1_share(.true.)


                !and if the bf_sublayer2 was also not linked to neighbor1
                else
                   
                   !and if the final alignment is such that the merged sublayer
                   !will exchange with neighbor1, links should be created to
                   !bf_sublayer1 and the state of the neighbor1 for bf_sublayer1
                   !should be updated
                   if(share_with_neighbor1) then
                      
                      call this%border_interface%link_neighbor1_to_bf_sublayer(
     $                     bf_sublayer1)

                      call bf_sublayer1%set_neighbor1_share(.true.)

                   end if
                end if
             end if


             !3) check if the links to neighbor2 should be updated
             !if the bf_sublayer1 was already linked to neighbor1 sublayers
             if(bf_sublayer1%can_exchange_with_neighbor2()) then
                
                !if the bf_sublayer2 was also linked to neighbor2 sublayers
                !as it will be merged, these links should be removed
                if(bf_sublayer2%can_exchange_with_neighbor2()) then
                   call this%border_interface%remove_link_from_neighbor2_to_bf_sublayer(
     $                  bf_sublayer2)

                end if

             !if the bf_sublayer1 was not linked to neighbor2 sublayers
             else

                !and on the contrary the bf_sublayer2 was linked to this neighbor
                !the links should be updated from bf_sublayer2 to bf_sublayer1
                !and the status of the neighbor2 for bf_sublayer1 should be updated
                if(bf_sublayer2%can_exchange_with_neighbor2()) then

                   call this%border_interface%update_link_from_neighbor2_to_bf_sublayer(
     $                  bf_sublayer1, bf_sublayer2)

                   call bf_sublayer1%set_neighbor2_share(.true.)


                !and if the bf_sublayer2 was also not linked to neighbor2
                else
                   
                   !and if the final alignment is such that the merged sublayer
                   !will exchange with neighbor2, links should be created to
                   !bf_sublayer1 and the state of the neighbor2 for bf_sublayer1
                   !should be updated
                   if(share_with_neighbor2) then
                      
                      call this%border_interface%link_neighbor2_to_bf_sublayer(
     $                     bf_sublayer1)

                      call bf_sublayer1%set_neighbor2_share(.true.)

                   end if
                end if
             end if


             !4) merge the content of the sublayers
             merged_sublayer => this%mainlayer_pointers(mainlayer_id)%merge_sublayers(
     $            bf_sublayer1, bf_sublayer2,
     $            nodes,alignment)

          else

             !2) check if the links to neighbor1 should be updated
             !if the bf_sublayer1 was already linked to neighbor1 sublayers
             if(bf_sublayer1%can_exchange_with_neighbor1()) then
                
                !if the bf_sublayer2 was also linked to neighbor1 sublayers
                !as it will be merged, these links should be removed
                if(bf_sublayer2%can_exchange_with_neighbor1()) then
                   call this%border_interface%remove_link_from_neighbor1_to_bf_sublayer(
     $                  bf_sublayer2)

                end if

             !if the bf_sublayer1 was not linked to neighbor1 sublayers
             else

                !and the bf_sublayer2 was linked to this neighbor, the links should
                !be updated from bf_sublayer2 to bf_sublayer1 and the status of the
                !neighbor1 for bf_sublayer1 should be updated
                if(bf_sublayer2%can_exchange_with_neighbor1()) then

                   call this%border_interface%update_link_from_neighbor1_to_bf_sublayer(
     $                  bf_sublayer1, bf_sublayer2)

                   call bf_sublayer1%set_neighbor1_share(.true.)

                end if
             end if


             !3) check if the links to neighbor2 should be updated
             !if the bf_sublayer1 was already linked to neighbor1 sublayers
             if(bf_sublayer1%can_exchange_with_neighbor2()) then
                
                !if the bf_sublayer2 was also linked to neighbor2 sublayers
                !as it will be merged, these links should be removed
                if(bf_sublayer2%can_exchange_with_neighbor2()) then
                   call this%border_interface%remove_link_from_neighbor2_to_bf_sublayer(
     $                  bf_sublayer2)

                end if

             !if the bf_sublayer1 was not linked to neighbor2 sublayers
             else

                !and on the contrary the bf_sublayer2 was linked to this neighbor
                !the links should be updated from bf_sublayer2 to bf_sublayer1
                !and the status of the neighbor2 for bf_sublayer1 should be updated
                if(bf_sublayer2%can_exchange_with_neighbor2()) then

                   call this%border_interface%update_link_from_neighbor2_to_bf_sublayer(
     $                  bf_sublayer1, bf_sublayer2)

                   call bf_sublayer1%set_neighbor2_share(.true.)

                end if
             end if


             !4) merge the content of the sublayers
             merged_sublayer => this%mainlayer_pointers(mainlayer_id)%merge_sublayers(
     $            bf_sublayer1, bf_sublayer2,
     $            nodes)

          end if


          !5) update the grid points that have been newly allocated
          !   with grid points from the neighboring sublayers
          call this%border_interface%update_grdpts_from_neighbors(bf_sublayer1)


c$$$          !6) if the merge is such that grid points between the two
c$$$          !   sublayers should be updated to prevent a line of inconsistent
c$$$          !   boundary points
c$$$          !    ________  ________        __________________
c$$$          !   |        ||        |      |                  |
c$$$          !   |        ||        |  ->  |                  |
c$$$          !   |        ||        |      |                  |
c$$$          print '(''bf_interface_class'')'
c$$$          print '(''merge_sublayers'')'
c$$$          stop 'not implemented yet'

       end function merge_sublayers


       !> remove a sublayer from its mainlayer
       subroutine remove_sublayer(this, sublayer_ptr, bf_mainlayer_id)

         implicit none

         class(bf_interface)       , intent(inout) :: this
         type(bf_sublayer), pointer, intent(inout) :: sublayer_ptr
         integer         , optional, intent(in)    :: bf_mainlayer_id

         
         integer :: mainlayer_id


         !> identify the mainlayer to which the sublayer belongs
         if(present(bf_mainlayer_id)) then
            mainlayer_id = bf_mainlayer_id
         else
            mainlayer_id = sublayer_ptr%get_localization()
         end if

         !> remove the sublayer from the table identifying the
         !> neighboring buffer layers
         if(sublayer_ptr%can_exchange_with_neighbor1()) then
            call this%border_interface%remove_link_from_neighbor1_to_bf_sublayer(
     $           sublayer_ptr)
         end if
         if(sublayer_ptr%can_exchange_with_neighbor2()) then
            call this%border_interface%remove_link_from_neighbor2_to_bf_sublayer(
     $           sublayer_ptr)
         end if

         !> remove the sublayer from the main layer
         call this%mainlayer_pointers(mainlayer_id)%remove_sublayer(sublayer_ptr)

       end subroutine remove_sublayer


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> subroutine converting general coordinates into
       !> the main layer ID (N,S,E,W)
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param general_coord
       !> integer table giving the general coordinates
       !
       !>@param mainlayer_id
       !> main layer cardinal coordinates
       !--------------------------------------------------------------
       function get_mainlayer_id(general_coord) result(mainlayer_id)

         implicit none

         integer(ikind), dimension(2), intent(in) :: general_coord
         integer                                  :: mainlayer_id

         if(general_coord(2).le.align_S) then
            mainlayer_id = S

         else
            if(general_coord(2).lt.(align_N)) then

               if(general_coord(1).le.align_W) then
                  mainlayer_id = W

               else
                  if(general_coord(1).lt.align_E) then
                     mainlayer_id = interior

                  else
                     mainlayer_id = E

                  end if
               end if
                     
            else
               mainlayer_id = N

            end if
         end if

       end function get_mainlayer_id


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> subroutine updating the interface pointers
       !> to the main layers
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> interface_abstract class encapsulating the pointers
       !> to the buffer main layers
       !
       !>@param general_coord
       !> table giving the general coordinates of the point analyzed
       !
       !>@param local_coord
       !> table giving the local coordinates of the point analyzed
       !> in the corresponding sublayer
       !
       !>@param tolerance_i
       !> integer indicating how far the gridpoint can be from the
       !> closest sublayer to be considered inside
       !
       !>@param sublayer
       !> pointer to the sublayer matching the general coordinates
       !> of the grid point
       !--------------------------------------------------------------
       function get_sublayer(
     $    this,
     $    general_coord,
     $    local_coord,
     $    tolerance_i,
     $    mainlayer_id_i)
     $    result(sublayer)

         implicit none

         class(bf_interface)         , intent(in)  :: this
         integer(ikind), dimension(2), intent(in)  :: general_coord
         integer(ikind), dimension(2), intent(out) :: local_coord
         integer       , optional    , intent(in)  :: tolerance_i
         integer       , optional    , intent(in)  :: mainlayer_id_i
         type(bf_sublayer), pointer                :: sublayer


         integer                     :: direction_tested
         integer                     :: mainlayer_id
         type(bf_mainlayer), pointer :: mainlayer
         integer                     :: tolerance
         logical                     :: grdpt_in_sublayer


         !< identification of the main layer
         if(present(mainlayer_id_i)) then
            mainlayer_id = mainlayer_id_i
         else
            mainlayer_id = get_mainlayer_id(general_coord)
         end if

         !< if the general coordinates match the interior,
         !> no sublayer matches the general coordinates
         !> and the sublayer pointer is nullified
         if(mainlayer_id.eq.interior) then
            nullify(sublayer)

         !< otherwise, the mainlayers are analyzed
         else

            !< check that the main layer exists
            !< if it does not exist, no sublayer can be saved inside
            !< and the pointer to the sublayer is nullified
            if(.not.this%mainlayer_pointers(mainlayer_id)%associated_ptr()) then
               nullify(sublayer)
                 
            !< if the main layer exists, the sublayers saved inside are
            !< checked to decide whether the grid point asked belongs to
            !< one of them or not
            else
               mainlayer => this%mainlayer_pointers(mainlayer_id)%get_ptr()
            
               !< check if sublayers are saved inside the mainlayer
               !> if no sublayers are saved inside the mainlayer,
               !> no existing sublayer can match the general coord
               !> and so the pointer to sublayer is nullified
               if(.not.associated(mainlayer%get_head_sublayer())) then
                  nullify(sublayer)
            
               !< otherwise, the sublayer corresponding to the general
               !> coordinates is searched by going through the different
               !> element of the doubled chained list
               else
                  sublayer => mainlayer%get_head_sublayer()
            
              	  !< processing the tolerance for matching a sublayer
            	  !> if no tolerence is provided, the default option is 0
                  if(.not.present(tolerance_i)) then
                     tolerance=0
                  else
                     tolerance=tolerance_i
                  end if
                  
                  !< if the mainlayer investigated is N,S,E or W, there can
                  !> be sublayers to be investigated
                  select case(mainlayer_id)
                    case(N,S)
                       direction_tested = x_direction
                    case(E,W)
                       direction_tested = y_direction
                  end select 

                  !check if the grid point belongs to the current sublayer
                  grdpt_in_sublayer =
     $                 (  general_coord(direction_tested).ge.
     $                 (sublayer%get_alignment(direction_tested,1)-bc_size-tolerance))
     $                 .and.(
     $                 general_coord(direction_tested).le.
     $                 (sublayer%get_alignment(direction_tested,2)+bc_size+tolerance))
            	
                  !go through the different sublayers
                  do while(.not.grdpt_in_sublayer)
            	        
            	     !if no matching sublayer can be found
            	     !nullify the corresponding pointer
                     if(.not.associated(sublayer%get_next())) then
                        nullify(sublayer)
                        exit
                     end if
                   
                     sublayer => sublayer%get_next()
                     grdpt_in_sublayer =
     $                    (general_coord(direction_tested).ge.
     $                    (sublayer%get_alignment(direction_tested,1)-bc_size-tolerance))
     $                    .and.(
     $                    general_coord(1).le.
     $                    (sublayer%get_alignment(direction_tested,2)+bc_size+tolerance))
            	
                  end do
            
               end if
            end if
         end if

         !< if a sublayer matching the general coordinates was found
         !> compute the local coordinates in this sublayer
         if(associated(sublayer)) then
            local_coord = sublayer%get_local_coord(general_coord)
         end if           
         
       end function get_sublayer


       !< extract the nodes at a general coordinates asked by the user
       function get_nodes(this, g_coords, interior_nodes) result(var)

         implicit none
         
         class(bf_interface)             , intent(in) :: this
         integer(ikind), dimension(2)    , intent(in) :: g_coords
         real(rkind), dimension(nx,ny,ne), intent(in) :: interior_nodes
         real(rkind), dimension(ne)                   :: var

         integer                      :: mainlayer_id
         type(bf_sublayer), pointer   :: sublayer
         integer(ikind), dimension(2) :: l_coords
         
         mainlayer_id = this%get_mainlayer_id(g_coords)
         if((g_coords(1).ge.1).and.
     $        (g_coords(1).le.nx).and.
     $        (g_coords(2).ge.1).and.
     $        (g_coords(2).le.ny)) then
            var = interior_nodes(g_coords(1),g_coords(2),:)
         else
            sublayer => this%get_sublayer(
     $               g_coords, l_coords, mainlayer_id_i=mainlayer_id)
            if(associated(sublayer)) then
               var = sublayer%get_nodes(l_coords)
            else
               print '(''bf_interface_class'')'
               print '(''get_nodes'')'
               print '(''cannot get sublayer'')'
               stop 'check way to get nodes'
            end if               
         end if

       end function get_nodes    


       !< if the buffer sublayer passed as argument has grid points
       !> in common with the buffer layers from other main layers, the
       !> grid points in common are updated from the neighboring buffer
       !> layers
       subroutine update_grdpts_from_neighbors(this, nbf_sublayer)

         implicit none

         class(bf_interface), intent(in)    :: this
         type(bf_sublayer)  , intent(inout) :: nbf_sublayer

         call this%border_interface%update_grdpts_from_neighbors(
     $        nbf_sublayer)

       end subroutine update_grdpts_from_neighbors


       !< if the buffer sublayer passed as argument has grid points
       !> in common with the buffer layers from other main layers, the
       !> grid points in common are updated in the neighboring buffer
       !> layers from the current buffer layer
       subroutine update_neighbor_grdpts(this, nbf_sublayer)

         implicit none

         class(bf_interface), intent(inout) :: this
         type(bf_sublayer)  , intent(inout) :: nbf_sublayer

         call this%border_interface%update_neighbor_grdpts(nbf_sublayer)

       end subroutine update_neighbor_grdpts


       !< update the grid points of the bf_sublayer after increase
       subroutine update_grdpts_after_increase(
     $     this, bf_sublayer_i, selected_grdpts)

         implicit none

         class(bf_interface)           , intent(inout) :: this
         type(bf_sublayer)             , intent(inout) :: bf_sublayer_i
         integer(ikind), dimension(:,:), intent(in)    :: selected_grdpts

         !compute the new grid points after the increase
         call bf_sublayer_i%update_grdpts_after_increase(
     $        selected_grdpts)

         !update the neighboring buffer layers
         call this%update_neighbor_grdpts(bf_sublayer_i)

       end subroutine update_grdpts_after_increase


       !< determine the neighboring buffer layers sharing grid points
       !> with the current buffer layer
       subroutine get_nbf_layers_sharing_grdpts_with(
     $     this, bf_sublayer_i, nbf1_list, nbf2_list, bf_mainlayer_id)

         implicit none
         
         class(bf_interface)       , intent(in)    :: this
         type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_i
         type(sbf_list)            , intent(inout) :: nbf1_list
         type(sbf_list)            , intent(inout) :: nbf2_list
         integer         , optional, intent(in)    :: bf_mainlayer_id


         !determine inside the list of neighboring buffer layer of
         !type 1 which ones share grid points with the current
         !buffer layer
         if(bf_sublayer_i%can_exchange_with_neighbor1()) then
            if(present(bf_mainlayer_id)) then
               call this%border_interface%get_nbf_layers_sharing_grdpts_with(
     $              1, bf_sublayer_i, nbf1_list, bf_mainlayer_id)
            else
               call this%border_interface%get_nbf_layers_sharing_grdpts_with(
     $              1, bf_sublayer_i, nbf1_list)
            end if
         end if


         !determine inside the list of neighboring buffer layer of
         !type 2 which ones share grid points with the current
         !buffer layer
         if(bf_sublayer_i%can_exchange_with_neighbor2()) then
            if(present(bf_mainlayer_id)) then
               call this%border_interface%get_nbf_layers_sharing_grdpts_with(
     $              2, bf_sublayer_i, nbf2_list, bf_mainlayer_id)
            else
               call this%border_interface%get_nbf_layers_sharing_grdpts_with(
     $              2, bf_sublayer_i, nbf2_list)
            end if
         end if         

       end subroutine get_nbf_layers_sharing_grdpts_with


       !< determine whether a buffer layer depends on its neighboring
       !> buffer layers
       function bf_layer_depends_on_neighbors(this, bf_sublayer_i, bf_mainlayer_id)
     $     result(dependent)

         implicit none

         class(bf_interface)       , intent(in) :: this
         type(bf_sublayer), pointer, intent(in) :: bf_sublayer_i
         integer         , optional, intent(in) :: bf_mainlayer_id
         logical                                :: dependent
       
         
         !determine if the buffer layer is sharing grid points
         !with its neighborign buffer layers of type 1
         if(bf_sublayer_i%can_exchange_with_neighbor1()) then
            
            if(present(bf_mainlayer_id)) then
               dependent = this%border_interface%bf_layer_depends_on_neighbors(
     $              1, bf_sublayer_i, bf_mainlayer_id)
            else
               dependent = this%border_interface%bf_layer_depends_on_neighbors(
     $              1, bf_sublayer_i)
            end if
         else
            dependent = .false.
         end if
         

         !determine if the buffer layer is sharing grid points
         !with its neighborign buffer layers of type 2
         if(.not.dependent) then
            
            if(bf_sublayer_i%can_exchange_with_neighbor2()) then
               
               if(present(bf_mainlayer_id)) then
                  dependent = this%border_interface%bf_layer_depends_on_neighbors(
     $                 2, bf_sublayer_i, bf_mainlayer_id)
               else
                  dependent = this%border_interface%bf_layer_depends_on_neighbors(
     $                 2, bf_sublayer_i)
               end if

            else
               dependent = .false.
            end if
            
         end if

       end function bf_layer_depends_on_neighbors


       !< check if one among the neighboring buffer layers cannot be removed
       function does_a_neighbor_remains(this, bf_sublayer_i, bf_mainlayer_id)
     $     result(a_neighbor_remains)

         implicit none

         class(bf_interface)       , intent(in) :: this
         type(bf_sublayer), pointer, intent(in) :: bf_sublayer_i
         integer         , optional, intent(in) :: bf_mainlayer_id
         logical                                :: a_neighbor_remains


         !determine if the buffer layer is sharing grid points
         !with its neighborign buffer layers of type 1
         if(bf_sublayer_i%can_exchange_with_neighbor1()) then
            
            if(present(bf_mainlayer_id)) then
               a_neighbor_remains = this%border_interface%does_a_neighbor_remains(
     $              1, bf_sublayer_i, bf_mainlayer_id)
            else
               a_neighbor_remains = this%border_interface%does_a_neighbor_remains(
     $              1, bf_sublayer_i)
            end if

         else
            a_neighbor_remains = .false.

         end if
         

         !determine if the buffer layer is sharing grid points
         !with its neighborign buffer layers of type 2
         if(.not.a_neighbor_remains) then
            
            if(bf_sublayer_i%can_exchange_with_neighbor2()) then
               
               if(present(bf_mainlayer_id)) then
                  a_neighbor_remains = this%border_interface%does_a_neighbor_remains(
     $                 2, bf_sublayer_i, bf_mainlayer_id)
               else
                  a_neighbor_remains = this%border_interface%does_a_neighbor_remains(
     $                 2, bf_sublayer_i)
               end if

            else
               a_neighbor_remains = .false.
            end if
            
         end if

       end function does_a_neighbor_remains


       !< print the content of the interface on external binary files
       subroutine print_binary(
     $     this,
     $     suffix_nodes, suffix_grdid, suffix_sizes,
     $     suffix_nb_sublayers_max)

         implicit none

         class(bf_interface), intent(in) :: this
         character(*)       , intent(in) :: suffix_nodes
         character(*)       , intent(in) :: suffix_grdid
         character(*)       , intent(in) :: suffix_sizes
         character(*)       , intent(in) :: suffix_nb_sublayers_max
         

         integer           :: i
         integer           :: nb_sublayers_max
         character(len=18) :: filename_format
         character(len=28) :: nb_sublayers_filename
                  

         !go through the buffer main layers and
         !print the content of each buffer layer
         !in seperate binary output files and 
         !determine the maximum number of sublayers
         nb_sublayers_max = 0
         do i=1, size(this%mainlayer_pointers,1)

            if(this%mainlayer_pointers(i)%associated_ptr()) then
               
               call this%mainlayer_pointers(i)%print_binary(suffix_nodes,
     $                                            suffix_grdid,
     $                                            suffix_sizes)

               nb_sublayers_max = max(
     $              nb_sublayers_max,
     $              this%mainlayer_pointers(i)%get_nb_sublayers())
            end if            
         end do


         !print the maximum number of sublayers in an output
         !binary file
         write(filename_format,
     $        '(''(A12,A'',I2,'')'')')
     $        len(suffix_nb_sublayers_max)

         write(nb_sublayers_filename, filename_format)
     $        'sublayers_nb',
     $        suffix_nb_sublayers_max

         call print_nb_sublayers_max(
     $        nb_sublayers_filename, nb_sublayers_max)

        end subroutine print_binary


        subroutine print_nb_sublayers_max(filename, nb_sublayers)

          implicit none

          character(*), intent(in) :: filename
          integer     , intent(in) :: nb_sublayers

          integer :: ios
          
          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) nb_sublayers
              close(unit=1)
           else
              stop 'file opening pb'
           end if

        end subroutine print_nb_sublayers_max

      end module bf_interface_class

      !> @file
      !> module implementing the subroutines needed to merge
      !> sublayers data
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the subroutines needed to merge
      !> sublayers data
      !
      !> @date
      ! 30_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_sublayers_merge_module

        use bf_mainlayer_abstract_class, only : bf_mainlayer_abstract
        use bf_sublayer_class          , only : bf_sublayer
        use parameters_bf_layer        , only : no_pt, exchange_pt, bc_pt,
     $                                          bc_interior_pt,
     $                                          interior_pt
        use parameters_constant        , only : x_direction, y_direction
        use parameters_input           , only : nx,ny,ne,bc_size
        use parameters_kind            , only : ikind, rkind

        implicit none

        private
        public ::
     $       merge_sublayers_N,
     $       merge_sublayers_S,
     $       merge_sublayers_E,
     $       merge_sublayers_W

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two Northern sublayers
        !> \image html  merge_sublayers_N.png
        !> \image latex merge_sublayers_N.eps width=10cm
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param sublayer1
        !> object encapsulating the sublayer data
        !
        !>@param sublayer2
        !> object encapsulating the sublayer data
        !
        !>@param nodes
        !> object encapsulating the data of the interior domain
        !
        !>@param alignment
        !> table specifying the position of the new merged layers
        !> compared to the previous one
        !--------------------------------------------------------------
        function merge_sublayers_N(
     $       mainlayer,
     $       sublayer1,
     $       sublayer2,
     $       alignment,
     $       neighbor_E_i,
     $       neighbor_W_i)
     $       result(merged_sublayer)

          implicit none

          class(bf_mainlayer_abstract)                 , intent(inout) :: mainlayer
          type(bf_sublayer), pointer                   , intent(inout) :: sublayer1
          type(bf_sublayer), pointer                   , intent(inout) :: sublayer2
          integer(ikind), dimension(2,2)     , optional, intent(in)    :: alignment
          logical                            , optional, intent(in)    :: neighbor_E_i
          logical                            , optional, intent(in)    :: neighbor_W_i
          type(bf_sublayer), pointer                                   :: merged_sublayer


          !local variables
          logical                                       :: neighbor_E, neighbor_W
          real(rkind)   , dimension(:,:,:), allocatable :: temp_nodes
          integer       , dimension(:,:)  , allocatable :: temp_grdptid
          integer(ikind), dimension(2,2)                :: temp_alignment

          integer(ikind), dimension(2) :: new_size

          type(bf_sublayer), pointer :: sublayer_first
          type(bf_sublayer), pointer :: sublayer_second
          type(bf_sublayer), pointer :: sublayer_third
          type(bf_sublayer), pointer :: sublayer_fourth

          integer(ikind) :: j_min
          integer(ikind) :: i_min1, i_min2, i_min3, i_min4
          integer(ikind) :: i_max1, i_max3
          integer(ikind) :: interior_i_max1, interior_i_max2, interior_i_max3
          integer(ikind) :: i_min_first_layer, i_max_first_layer
          integer(ikind) :: i_min_second_layer, i_max_second_layer

          ! initialize the neighbors
          if(present(neighbor_E_i)) then
             neighbor_E=neighbor_E_i
          else
             neighbor_E=.false.
          end if

          if(present(neighbor_W_i)) then
             neighbor_W=neighbor_W_i
          else
             neighbor_W=.false.
          end if


          !< new size of the table
          if(present(alignment)) then
             new_size = get_new_size_NS(sublayer1,sublayer2,alignment)
          else
             new_size = get_new_size_NS(sublayer1,sublayer2)
          end if
          
          !< match betwen the tables
          call get_match_sublayers_N(
     $         sublayer1,
     $         sublayer2,
     $         sublayer_first,
     $         sublayer_second,
     $         sublayer_third,
     $         sublayer_fourth)

          !< allocate the new table for the nodes
          allocate(temp_nodes(new_size(1),new_size(2),ne))
          
          !< determine the matching parameters between the tables
          j_min = min(
     $         size(sublayer_first%element%nodes,2),
     $         size(sublayer_second%element%nodes,2))
          
          if(present(alignment)) then
             call get_match_direction(
     $            i_min1, i_min2, i_min3, i_min4,
     $            i_max1, i_max3,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            sublayer_first, sublayer_second,
     $            x_direction,
     $            alignment=alignment,
     $            j_min_first_layer=i_min_first_layer,
     $            j_max_first_layer=i_max_first_layer,
     $            j_min_second_layer=i_min_second_layer,
     $            j_max_second_layer=i_max_second_layer)
          else
             call get_match_direction(
     $            i_min1, i_min2, i_min3, i_min4,
     $            i_max1, i_max3,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            sublayer_first, sublayer_second,
     $            x_direction,
     $            j_min_first_layer=i_min_first_layer,
     $            j_max_first_layer=i_max_first_layer,
     $            j_min_second_layer=i_min_second_layer,
     $            j_max_second_layer=i_max_second_layer)
          end if
        

          ! fill the new nodes
          ! to fill efficiently in memory the new nodes, we divide
          ! the filling in several parts:
          !
          !                i_min1      i_min2 i_min3        i_min4
          !                   |            |   |              |
          !                   |       interior_i_max2         |
          !                   |            | | |              |
          !       interior_i_max1  i_max1  | | |  i_max2  interior_i_max3
          !                 | | _____|____ | | | _____|______ | |
          !                / \|/          \|/ \|/            \|/ \
          !             _  _______________________________________
          !            /  |   |            |   |              |   |
          !           /   |   |  third     |   |   fourth     |   |
          ! no_pt____/    |   |  sublayer  |   |   sublayer   |   |
          !          \    |   |____________|   |______________|   |_j_min
          !           \   |   |            |   |              |   |
          !            \_ |___|            |___|              |___|
          ! bc_pt       - |3 3|            |3 3|              |3 3|
          ! bc_interior - |2 2|            |2 2|              |2 2|
          ! interior    - |1 1|  first     |1 1|   second     |1 1|
          !          __/- |0 0|  sublayer  |0 0|   sublayer   |0 0|
          ! exchange   \- |0_0|____________|0_0|______________|0_0|
          !                 .                .                  . 
          !                /|\              /|\                /|\
          !                 |________________|__________________| 
          !                                  |
          !                              copy from    
          !                           interior domain
          !
          !--------------------------------------------------------
          !                         Fig. 1
          !--------------------------------------------------------
          if(present(alignment)) then
             call merge_nodes_N(
     $            temp_nodes,
     $            sublayer_first,
     $            sublayer_second,
     $            sublayer_third,
     $            sublayer_fourth,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_max1, i_max3,
     $            i_min1, i_min2, i_min3, i_min4,
     $            j_min,
     $            alignment=alignment)

          else
             call merge_nodes_N(
     $            temp_nodes,
     $            sublayer_first,
     $            sublayer_second,
     $            sublayer_third,
     $            sublayer_fourth,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_max1, i_max3,
     $            i_min1, i_min2, i_min3, i_min4,
     $            j_min)
          end if
          

          !< reallocate the nodes for sublayer_first and deallocate
          !> the one for sublayer_second
          call MOVE_ALLOC(temp_nodes,sublayer_first%element%nodes)
          deallocate(sublayer_second%element%nodes)

c$$$          !< allocate the new table for the gridpt_id
c$$$          allocate(temp_grdptid(new_size(1),new_size(2)))
c$$$
c$$$          !< fill the grdpt_id
c$$$          call merge_grdpts_id_N(
c$$$     $         temp_grdptid,
c$$$     $         sublayer_first,
c$$$     $         sublayer_second,
c$$$     $         sublayer_third,
c$$$     $         sublayer_fourth,
c$$$     $         interior_i_max1, interior_i_max2, interior_i_max3,
c$$$     $         i_max1, i_max3,
c$$$     $         i_min_first_layer, i_max_first_layer,
c$$$     $         i_min_second_layer, i_max_second_layer,
c$$$     $         i_min1, i_min2, i_min3, i_min4,
c$$$     $         j_min,
c$$$     $         neighbor_E,
c$$$     $         neighbor_W)
c$$$
c$$$         !< reallocate the grdpts_id for sublayer_first and deallocate
c$$$         !> the one for sublayer_second
c$$$         call MOVE_ALLOC(temp_grdptid,sublayer_first%element%grdpts_id)
c$$$         deallocate(sublayer_second%element%grdpts_id)
c$$$
c$$$         !fix the frontier between the two sublayers merged
c$$$         !for the grdptid
c$$$         if(interior_i_max2.eq.0) then
c$$$            call fix_frontier_after_merge_N(
c$$$     $           sublayer_first%element%grdpts_id,
c$$$     $           i_min1+i_max1,
c$$$     $           j_min)
c$$$         end if
            
         !< determine the alignment for the merged sublayer
         if(present(alignment)) then
            temp_alignment(1,1) = min(sublayer_first%element%alignment(1,1),alignment(1,1))
            temp_alignment(1,2) = max(sublayer_second%element%alignment(1,2),alignment(1,2))
         else
            temp_alignment(1,1) = sublayer_first%element%alignment(1,1)
            temp_alignment(1,2) = sublayer_second%element%alignment(1,2)
         end if
         temp_alignment(2,1) = nx-1
         temp_alignment(2,2) = nx

         !< replace the alignment of the first sublayer
         sublayer_first%element%alignment = temp_alignment

         !< now all the parameters of the merge sublayer%element
         !> have been reinitialized:
         !> localization : remains the same
         !> alignment    : done
         !> nodes        : done
         !> grdpts_id    : done
         !> we need to change the parameters of the sublayer itself
         !> element : done
         !> next    : not done
         !> prev    : remains the same
         if(associated(sublayer_second%next)) then
            sublayer_first%next      => sublayer_second%next
            sublayer_first%next%prev => sublayer_first
         else
            nullify(sublayer_first%next)
            mainlayer%tail_sublayer => sublayer_first
         end if

         nullify(sublayer_second%prev)
         nullify(sublayer_second%next)
         deallocate(sublayer_second)


         !<return the merged sublayer
         merged_sublayer => sublayer_first

        end function merge_sublayers_N


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two Southern sublayers
        !> \image html  merge_sublayers_S.png
        !> \image latex merge_sublayers_S.eps width=10cm
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param sublayer1
        !> object encapsulating the sublayer data
        !
        !>@param sublayer2
        !> object encapsulating the sublayer data
        !
        !>@param nodes
        !> object encapsulating the data of the interior domain
        !
        !>@param alignment
        !> table specifying the position of the new merged layers
        !> compared to the previous one
        !--------------------------------------------------------------
        function merge_sublayers_S(
     $       mainlayer,
     $       sublayer1,
     $       sublayer2,
     $       nodes,
     $       alignment,
     $       neighbor_E_i,
     $       neighbor_W_i)
     $       result(merged_sublayer)

          implicit none

          class(bf_mainlayer_abstract)                 , intent(inout) :: mainlayer
          type(bf_sublayer), pointer                   , intent(inout) :: sublayer1
          type(bf_sublayer), pointer                   , intent(inout) :: sublayer2
          real(rkind)   , dimension(nx,ny,ne), optional, intent(in)    :: nodes
          integer(ikind), dimension(2,2)     , optional, intent(in)    :: alignment
          logical                            , optional, intent(in)    :: neighbor_E_i
          logical                            , optional, intent(in)    :: neighbor_W_i          
          type(bf_sublayer), pointer                                   :: merged_sublayer


          !local variables
          logical                                       :: neighbor_E, neighbor_W
          real(rkind)   , dimension(:,:,:), allocatable :: temp_nodes
          integer       , dimension(:,:)  , allocatable :: temp_grdptid
          integer(ikind), dimension(2,2)                :: temp_alignment

          integer(ikind), dimension(2) :: new_size

          type(bf_sublayer), pointer :: sublayer_first
          type(bf_sublayer), pointer :: sublayer_second
          type(bf_sublayer), pointer :: sublayer_third
          type(bf_sublayer), pointer :: sublayer_fourth

          integer(ikind) :: j_min, j_min1, j_min3, j_max
          integer(ikind) :: i_min1, i_min2, i_min3, i_min4
          integer(ikind) :: i_max1, i_max3
          integer(ikind) :: interior_i_max1, interior_i_max2, interior_i_max3
          integer(ikind) :: i_min_first_layer, i_max_first_layer
          integer(ikind) :: i_min_second_layer, i_max_second_layer

          if((present(alignment).and..not.present(nodes)).or.
     $         (.not.present(alignment).and.present(nodes))) then
             print '(''bf_sublayers_merge_module'')'
             print '(''merge_sublayers_S'')'
             print '(''alignment and nodes arguments should be'')'
             print '(''supplied together'')'
             stop 'change arguments'
          end if

          !< initialize the neighbors
          if(present(neighbor_E_i)) then
             neighbor_E=neighbor_E_i
          else
             neighbor_E=.false.
          end if

          if(present(neighbor_W_i)) then
             neighbor_W=neighbor_W_i
          else
             neighbor_W=.false.
          end if

          !< new size of the table
          if(present(alignment)) then
             new_size = get_new_size_NS(sublayer1,sublayer2,alignment)
          else
             new_size = get_new_size_NS(sublayer1,sublayer2)
          end if
          
          !< match betwen the tables
          call get_match_sublayers_S(
     $         sublayer1,
     $         sublayer2,
     $         sublayer_first,
     $         sublayer_second,
     $         sublayer_third,
     $         sublayer_fourth,
     $         j_min1, j_min3, j_max)

          !< allocate the new table for the nodes
          allocate(temp_nodes(new_size(1),new_size(2),ne))
          
          !< determine the matching parameters between the tables
          j_min = j_max - min(
     $         size(sublayer_third%element%nodes,2),
     $         size(sublayer_fourth%element%nodes,2))
          
          !< determine the borders in the x-direction
          if(present(alignment)) then
             call get_match_direction(
     $            i_min1, i_min2, i_min3, i_min4,
     $            i_max1, i_max3,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            sublayer_third, sublayer_fourth,
     $            x_direction,
     $            alignment=alignment,
     $            j_min_first_layer=i_min_first_layer,
     $            j_max_first_layer=i_max_first_layer,
     $            j_min_second_layer=i_min_second_layer,
     $            j_max_second_layer=i_max_second_layer)

          else
             call get_match_direction(
     $            i_min1, i_min2, i_min3, i_min4,
     $            i_max1, i_max3,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            sublayer_third, sublayer_fourth,
     $            x_direction,
     $            j_min_first_layer=i_min_first_layer,
     $            j_max_first_layer=i_max_first_layer,
     $            j_min_second_layer=i_min_second_layer,
     $            j_max_second_layer=i_max_second_layer)
          end if         

          ! fill the new nodes
          ! to fill efficiently in memory the new nodes, we divide
          ! the filling in several parts:
          !
          !                i_min1      i_min2 i_min3        i_min4
          !                   |            |   |              |
          !                   |       interior_i_max2         |
          !                   |            | | |              |
          !       interior_i_max1  i_max1  | | |  i_max2  interior_i_max3
          !                 | | _____|____ | | | _____|______ | |
          !                / \|/          \|/ \|/            \|/ \
          !             _  ---------------------------------------
          !          __/- |0 0|            |0 0|              |0 0|
          ! exchange   \- |0 0|            |0 0|              |0 0|
          ! interior    - |1 1|  third     |1 1|   fourth     |1 1|
          ! bc_interior - |2 2|  sublayer  |2 2|   sublayer   |2 2|
          ! bc_pt       - |3_3|            |3_3|              |3_3|
          !            /  |   |____________|   |______________|   |
          !           /   |   |            |   |              |   |
          ! no_pt____/    |   |  first     |   |   second     |   |
          !          \    |   |  sublayer  |   |   sublayer   |   |
          !           \   |   |            |   |              |   |
          !            \_ |___|____________|___|______________|___|
          !                 .                .                  . 
          !                /|\              /|\                /|\
          !                 |________________|__________________| 
          !                                  |
          !                              copy from    
          !                           interior domain
          !
          !--------------------------------------------------------
          !                           Fig. 2
          !--------------------------------------------------------
          if(present(alignment)) then
             call merge_nodes_S(
     $            temp_nodes,
     $            sublayer_first,
     $            sublayer_second,
     $            sublayer_third,
     $            sublayer_fourth,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_max1, i_max3,
     $            i_min1, i_min2, i_min3, i_min4,
     $            j_min,
     $            j_min1, j_min3, j_max,
     $            nodes=nodes,
     $            alignment=alignment)
          else
             call merge_nodes_S(
     $            temp_nodes,
     $            sublayer_first,
     $            sublayer_second,
     $            sublayer_third,
     $            sublayer_fourth,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_max1, i_max3,
     $            i_min1, i_min2, i_min3, i_min4,
     $            j_min,
     $            j_min1, j_min3, j_max)
          end if

          !< reallocate the nodes for sublayer_first and deallocate
          !> the one for sublayer_second
          call MOVE_ALLOC(temp_nodes,sublayer_third%element%nodes)
          deallocate(sublayer_fourth%element%nodes)

          !< allocate the new table for the gridpt_id
          allocate(temp_grdptid(new_size(1),new_size(2)))

          !< fill the grdpt_id
          call merge_grdpts_id_S(
     $         temp_grdptid,
     $         sublayer_first,
     $         sublayer_second,
     $         sublayer_third,
     $         sublayer_fourth,
     $         interior_i_max1, interior_i_max2, interior_i_max3,
     $         i_max1, i_max3,
     $         i_min_first_layer, i_max_first_layer,
     $         i_min_second_layer, i_max_second_layer,
     $         i_min1, i_min2, i_min3, i_min4,
     $         j_min, j_min1, j_min3, j_max,
     $         neighbor_E, neighbor_W)
          

         !< reallocate the grdpts_id for sublayer_first and deallocate
         !> the one for sublayer_second
         call MOVE_ALLOC(temp_grdptid,sublayer_third%element%grdpts_id)
         deallocate(sublayer_fourth%element%grdpts_id)

         !fix the frontier between the two sublayers merged
         !for the grdptid
         if(interior_i_max2.eq.0) then
            call fix_frontier_after_merge_S(
     $           sublayer_third%element%grdpts_id,
     $           i_min1+i_max1,
     $           j_min,
     $           j_max)
         end if
            
         !< determine the alignment for the merged sublayer
         if(present(alignment)) then
            temp_alignment(1,1) = min(sublayer_third%element%alignment(1,1),alignment(1,1))
            temp_alignment(1,2) = max(sublayer_fourth%element%alignment(1,2),alignment(1,2))
         else
            temp_alignment(1,1) = sublayer_third%element%alignment(1,1)
            temp_alignment(1,2) = sublayer_fourth%element%alignment(1,2)
         end if
         temp_alignment(2,1) = 1
         temp_alignment(2,2) = bc_size

         !< replace the alignment of the first sublayer
         sublayer_third%element%alignment = temp_alignment

         !< now all the parameters of the merge sublayer%element
         !> have been reinitialized:
         !> localization : remains the same
         !> alignment    : done
         !> nodes        : done
         !> grdpts_id    : done
         !> we need to change the parameters of the sublayer itself
         !> element : done
         !> next    : not done
         !> prev    : remains the same
         if(associated(sublayer_fourth%next)) then
            sublayer_third%next      => sublayer_fourth%next
            sublayer_third%next%prev => sublayer_third
         else
            nullify(sublayer_third%next)
            mainlayer%tail_sublayer => sublayer_third
         end if

         nullify(sublayer_fourth%prev)
         nullify(sublayer_fourth%next)
         deallocate(sublayer_fourth)


         !<return the merged sublayer
         merged_sublayer => sublayer_third

        end function merge_sublayers_S


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two Eastern sublayers
        !> \image html  merge_sublayers_E.png
        !> \image latex merge_sublayers_E.eps width=10cm
        !
        !> @date
        !> 06_05_2014 - initial version - J.L. Desmarais
        !
        !>@param sublayer1
        !> object encapsulating the sublayer data
        !
        !>@param sublayer2
        !> object encapsulating the sublayer data
        !
        !>@param nodes
        !> object encapsulating the data of the interior domain
        !
        !>@param alignment
        !> table specifying the position of the new merged layers
        !> compared to the previous one
        !--------------------------------------------------------------
        function merge_sublayers_E(
     $       mainlayer,
     $       sublayer1,
     $       sublayer2,
     $       nodes,
     $       alignment,
     $       neighbor_N_i,
     $       neighbor_S_i)
     $       result(merged_sublayer)

          implicit none

          class(bf_mainlayer_abstract)                 , intent(inout) :: mainlayer
          type(bf_sublayer), pointer                   , intent(inout) :: sublayer1
          type(bf_sublayer), pointer                   , intent(inout) :: sublayer2
          real(rkind)   , dimension(nx,ny,ne), optional, intent(in)    :: nodes
          integer(ikind), dimension(2,2)     , optional, intent(in)    :: alignment
          logical                            , optional, intent(in)    :: neighbor_N_i
          logical                            , optional, intent(in)    :: neighbor_S_i
          type(bf_sublayer), pointer                                   :: merged_sublayer


          !local variables
          logical                                       :: neighbor_N, neighbor_S
          real(rkind)   , dimension(:,:,:), allocatable :: temp_nodes
          integer       , dimension(:,:)  , allocatable :: temp_grdptid
          integer(ikind), dimension(2,2)                :: temp_alignment

          integer(ikind), dimension(2) :: new_size

          type(bf_sublayer), pointer :: sublayer_first
          type(bf_sublayer), pointer :: sublayer_third

          integer(ikind) :: i_min, i_max
          integer(ikind) :: j_min1, j_min2, j_min3, j_min4
          integer(ikind) :: j_max1, j_max3
          integer(ikind) :: interior_j_max1, interior_j_max2, interior_j_max3


          ! we need to make sure that the nodes and
          ! the alignment arguments are submitted
          ! together
          if((present(alignment).and..not.present(nodes)).or.
     $         (.not.present(alignment).and.present(nodes))) then
             print '(''bf_sublayers_merge_module'')'
             print '(''merge_sublayers_E'')'
             print '(''alignment and nodes arguments should be'')'
             print '(''supplied together'')'
             stop 'change arguments'
          end if


          ! initialize the neighbors
          if(present(neighbor_N_i)) then
             neighbor_N=neighbor_N_i
          else
             neighbor_N=.false.
          end if

          if(present(neighbor_S_i)) then
             neighbor_S=neighbor_S_i
          else
             neighbor_S=.false.
          end if


          ! new size of the table
          if(present(alignment)) then
             new_size = get_new_size_EW(sublayer1,sublayer2,alignment)
          else
             new_size = get_new_size_EW(sublayer1,sublayer2)
          end if

          
          ! match betwen the tables
          call get_match_sublayers_EW(
     $         sublayer1,
     $         sublayer2,
     $         sublayer_first,
     $         sublayer_third)

          ! allocate the new table for the nodes
          allocate(temp_nodes(new_size(1),new_size(2),ne))
          
          ! determine the matching parameters between the tables
          i_min = min(
     $         size(sublayer_first%element%nodes,1),
     $         size(sublayer_third%element%nodes,1))

          i_max = max(
     $         size(sublayer_first%element%nodes,1),
     $         size(sublayer_third%element%nodes,1))
          
          if(present(alignment)) then
             call get_match_direction(
     $            j_min1, j_min2, j_min3, j_min4,
     $            j_max1, j_max3,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            sublayer_first, sublayer_third,
     $            y_direction,
     $            alignment=alignment)
          else
             call get_match_direction(
     $            j_min1, j_min2, j_min3, j_min4,
     $            j_max1, j_max3,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            sublayer_first, sublayer_third,
     $            y_direction)
          end if
        

          ! fill the new nodes
          ! to fill efficiently in memory the new nodes, we divide
          ! the filling in several parts:
          !
          !              exchange_pt
          !                 | interior_pt
          !                 |  | bc_interior_pt
          !                 |  | | bc_pt
          !                 |  | | |           no_pt
          !                 |  | | |  _____|____
          !                / \ | | | /          \
          !               ------------------------
          !            __\|0 0 1 2 3|            |\___interior_j_max1
          !           |  /|0 0 1 2 3|            |/                  
          !           |   ------------------------                   
          !           |   |                      |\                  
          !           |   | third sublayer       | ---j_max1         
          !           |   |                      |/                  
          !  copy     |   ------------------------                   
          !  from     |__\|0 0 1 2 3|            |\___interior_j_max2
          !  interior |  /|0 0 1 2 3|            |/                  
          !  domain   |   ------------------------                   
          !           |   |                      |\                  
          !           |   | first sublayer       | ---j_max2         
          !           |   |                      |/                  
          !           |   ------------------------                   
          !           |__\|0 0 1 2 3|            |\___interior_j_max3
          !              /|0 0 1 2 3|            |/                   
          !               ------------------------
          !                         |            |
          !                       i_min        i_max
          !
          !--------------------------------------------------------
          !                        Fig. 3
          !--------------------------------------------------------
          if(present(alignment)) then
             call merge_nodes_E(
     $            temp_nodes,
     $            sublayer_first,
     $            sublayer_third,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_max1, j_max3,
     $            j_min1, j_min2, j_min3, j_min4,
     $            nodes=nodes,
     $            alignment=alignment)

          else
             call merge_nodes_E(
     $            temp_nodes,
     $            sublayer_first,
     $            sublayer_third,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_max1, j_max3,
     $            j_min1, j_min2, j_min3, j_min4)
          end if
          

          !< reallocate the nodes for sublayer_first and deallocate
          !> the one for sublayer_third
          call MOVE_ALLOC(temp_nodes,sublayer_first%element%nodes)
          deallocate(sublayer_third%element%nodes)

          !< allocate the new table for the gridpt_id
          allocate(temp_grdptid(new_size(1),new_size(2)))

          !< fill the grdpt_id
          call merge_grdpts_id_E(
     $         temp_grdptid,
     $         sublayer_first,
     $         sublayer_third,
     $         interior_j_max1, interior_j_max2, interior_j_max3,
     $         j_max1, j_max3,
     $         j_min1, j_min2, j_min3, j_min4,
     $         i_max,
     $         neighbor_N, neighbor_S)

         !< reallocate the grdpts_id for sublayer_first and deallocate
         !> the one for sublayer_second
         call MOVE_ALLOC(temp_grdptid,sublayer_first%element%grdpts_id)
         deallocate(sublayer_third%element%grdpts_id)

         !fix the frontier between the two sublayers merged
         !for the grdptid
         if(interior_j_max2.eq.0) then
            call fix_frontier_after_merge_E(
     $           sublayer_first%element%grdpts_id,
     $           j_min1+j_max1,
     $           i_min)
         end if
            
         !< determine the alignment for the merged sublayer
         if(present(alignment)) then
            temp_alignment(2,1) = min(sublayer_first%element%alignment(2,1),alignment(2,1))
            temp_alignment(2,2) = max(sublayer_third%element%alignment(2,2),alignment(2,2))
         else
            temp_alignment(2,1) = sublayer_first%element%alignment(2,1)
            temp_alignment(2,2) = sublayer_third%element%alignment(2,2)
         end if
         temp_alignment(1,1) = nx-1
         temp_alignment(1,2) = nx

         !< replace the alignment of the first sublayer
         sublayer_first%element%alignment = temp_alignment

         !< now all the parameters of the merge sublayer%element
         !> have been reinitialized:
         !> localization : remains the same
         !> alignment    : done
         !> nodes        : done
         !> grdpts_id    : done
         !> we need to change the parameters of the sublayer itself
         !> element : done
         !> next    : not done
         !> prev    : remains the same
         if(associated(sublayer_third%next)) then
            sublayer_first%next      => sublayer_third%next
            sublayer_first%next%prev => sublayer_first
         else
            nullify(sublayer_first%next)
            mainlayer%tail_sublayer => sublayer_first
         end if

         nullify(sublayer_third%prev)
         nullify(sublayer_third%next)
         deallocate(sublayer_third)


         !<return the merged sublayer
         merged_sublayer => sublayer_first

        end function merge_sublayers_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two Western sublayers
        !> \image html  merge_sublayers_W.png
        !> \image latex merge_sublayers_W.eps width=10cm
        !
        !> @date
        !> 06_05_2014 - initial version - J.L. Desmarais
        !
        !>@param sublayer1
        !> object encapsulating the sublayer data
        !
        !>@param sublayer2
        !> object encapsulating the sublayer data
        !
        !>@param nodes
        !> object encapsulating the data of the interior domain
        !
        !>@param alignment
        !> table specifying the position of the new merged layers
        !> compared to the previous one
        !--------------------------------------------------------------
        function merge_sublayers_W(
     $       mainlayer,
     $       sublayer1,
     $       sublayer2,
     $       nodes,
     $       alignment,
     $       neighbor_N_i,
     $       neighbor_S_i)
     $       result(merged_sublayer)

          implicit none

          class(bf_mainlayer_abstract)                 , intent(inout) :: mainlayer
          type(bf_sublayer), pointer                   , intent(inout) :: sublayer1
          type(bf_sublayer), pointer                   , intent(inout) :: sublayer2
          real(rkind)   , dimension(nx,ny,ne), optional, intent(in)    :: nodes
          integer(ikind), dimension(2,2)     , optional, intent(in)    :: alignment
          logical                            , optional, intent(in)    :: neighbor_N_i
          logical                            , optional, intent(in)    :: neighbor_S_i
          type(bf_sublayer), pointer                                   :: merged_sublayer


          !local variables
          logical                                       :: neighbor_N, neighbor_S
          real(rkind)   , dimension(:,:,:), allocatable :: temp_nodes
          integer       , dimension(:,:)  , allocatable :: temp_grdptid
          integer(ikind), dimension(2,2)                :: temp_alignment

          integer(ikind), dimension(2) :: new_size

          type(bf_sublayer), pointer :: sublayer_first
          type(bf_sublayer), pointer :: sublayer_third

          integer(ikind) :: i_max, i_min
          integer(ikind) :: j_min1, j_min2, j_min3, j_min4
          integer(ikind) :: j_max1, j_max3
          integer(ikind) :: interior_j_max1, interior_j_max2, interior_j_max3


          ! we need to make sure that the nodes and
          ! the alignment arguments are submitted
          ! together
          if((present(alignment).and..not.present(nodes)).or.
     $         (.not.present(alignment).and.present(nodes))) then
             print '(''bf_sublayers_merge_module'')'
             print '(''merge_sublayers_W'')'
             print '(''alignment and nodes arguments should be'')'
             print '(''supplied together'')'
             stop 'change arguments'
          end if


          ! initialize the neighbors
          if(present(neighbor_N_i)) then
             neighbor_N=neighbor_N_i
          else
             neighbor_N=.false.
          end if

          if(present(neighbor_S_i)) then
             neighbor_S=neighbor_S_i
          else
             neighbor_S=.false.
          end if


          ! new size of the table
          if(present(alignment)) then
             new_size = get_new_size_EW(sublayer1,sublayer2,alignment)
          else
             new_size = get_new_size_EW(sublayer1,sublayer2)
          end if
          

          ! match betwen the tables
          call get_match_sublayers_EW(
     $         sublayer1,
     $         sublayer2,
     $         sublayer_first,
     $         sublayer_third)


          ! allocate the new table for the nodes
          allocate(temp_nodes(new_size(1),new_size(2),ne))
          
          ! determine the matching parameters between the tables
          i_max = max(
     $         size(sublayer_first%element%nodes,1),
     $         size(sublayer_third%element%nodes,1))

          i_min = i_max-min(
     $         size(sublayer_first%element%nodes,1),
     $         size(sublayer_third%element%nodes,1))
          
          if(present(alignment)) then
             call get_match_direction(
     $            j_min1, j_min2, j_min3, j_min4,
     $            j_max1, j_max3,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            sublayer_first, sublayer_third,
     $            y_direction,
     $            alignment=alignment)
          else
             call get_match_direction(
     $            j_min1, j_min2, j_min3, j_min4,
     $            j_max1, j_max3,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            sublayer_first, sublayer_third,
     $            y_direction)
          end if
        

          ! fill the new nodes
          ! to fill efficiently in memory the new nodes, we divide
          ! the filling in several parts:
          !
          !                            exchange_pt
          !                      interior_pt | 
          !                bc_interior_pt |  | 
          !                      bc _pt | |  | 
          !                  no_pt    | | |  | 
          !                ____|___   | | |  | 
          !               /        \  | | | / \ 
          !               ----------------------
          !            __\|         | 3 2 1 0 0|\___interior_j_max1
          !           |  /|         | 3 2 1 0 0|/                  
          !           |   ----------------------                   
          !           |   |                    |\                  
          !           |   | third_sublayer     | ---j_max1         
          !           |   |                    |/                  
          !  copy     |   ----------------------                   
          !  from     |__\|         | 3 2 1 0 0|\___interior_j_max2
          !  interior |  /|         | 3 2 1 0 0|/                  
          !  domain   |   ----------------------                   
          !           |   |                    |\                  
          !           |   | first_sublayer     | ---j_max2         
          !           |   |                    |/                  
          !           |   ----------------------                   
          !           |__\|         | 3 2 1 0 0|\___interior_j_max3
          !              /|         | 3 2 1 0 0|/                  
          !               ----------------------
          !                         |          |
          !                       i_min      i_max
          !
          !--------------------------------------------------------
          !                         Fig. 4
          !--------------------------------------------------------
          if(present(alignment)) then
             call merge_nodes_W(
     $            temp_nodes,
     $            sublayer_first,
     $            sublayer_third,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_max1, j_max3,
     $            j_min1, j_min2, j_min3, j_min4,
     $            i_max,
     $            nodes=nodes,
     $            alignment=alignment)

          else
             call merge_nodes_W(
     $            temp_nodes,
     $            sublayer_first,
     $            sublayer_third,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_max1, j_max3,
     $            j_min1, j_min2, j_min3, j_min4,
     $            i_max)
          end if
          

          !< reallocate the nodes for sublayer_first and deallocate
          !> the one for sublayer_third
          call MOVE_ALLOC(temp_nodes,sublayer_first%element%nodes)
          deallocate(sublayer_third%element%nodes)

          !< allocate the new table for the gridpt_id
          allocate(temp_grdptid(new_size(1),new_size(2)))

          !< fill the grdpt_id
          call merge_grdpts_id_W(
     $         temp_grdptid,
     $         sublayer_first,
     $         sublayer_third,
     $         interior_j_max1, interior_j_max2, interior_j_max3,
     $         j_max1, j_max3,
     $         j_min1, j_min2, j_min3, j_min4,
     $         i_max,
     $         neighbor_N, neighbor_S)

         !< reallocate the grdpts_id for sublayer_first and deallocate
         !> the one for sublayer_second
         call MOVE_ALLOC(temp_grdptid,sublayer_first%element%grdpts_id)
         deallocate(sublayer_third%element%grdpts_id)

         !fix the frontier between the two sublayers merged
         !for the grdptid
         if(interior_j_max2.eq.0) then
            call fix_frontier_after_merge_W(
     $           sublayer_first%element%grdpts_id,
     $           j_min1+j_max1,
     $           i_min,
     $           i_max)
         end if
            
         !< determine the alignment for the merged sublayer
         if(present(alignment)) then
            temp_alignment(2,1) = min(sublayer_first%element%alignment(2,1),alignment(2,1))
            temp_alignment(2,2) = max(sublayer_third%element%alignment(2,2),alignment(2,2))
         else
            temp_alignment(2,1) = sublayer_first%element%alignment(2,1)
            temp_alignment(2,2) = sublayer_third%element%alignment(2,2)
         end if
         temp_alignment(1,1) = 1
         temp_alignment(1,2) = bc_size

         !< replace the alignment of the first sublayer
         sublayer_first%element%alignment = temp_alignment

         !< now all the parameters of the merge sublayer%element
         !> have been reinitialized:
         !> localization : remains the same
         !> alignment    : done
         !> nodes        : done
         !> grdpts_id    : done
         !> we need to change the parameters of the sublayer itself
         !> element : done
         !> next    : not done
         !> prev    : remains the same
         if(associated(sublayer_third%next)) then
            sublayer_first%next      => sublayer_third%next
            sublayer_first%next%prev => sublayer_first
         else
            nullify(sublayer_first%next)
            mainlayer%tail_sublayer => sublayer_first
         end if

         nullify(sublayer_third%prev)
         nullify(sublayer_third%next)
         deallocate(sublayer_third)


         !<return the merged sublayer
         merged_sublayer => sublayer_first

        end function merge_sublayers_W



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine fixing the frontier between two freshly
        !> merged North sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpt_id
        !> table encapsulating the grid points
        !
        !>@param frontier
        !> index indicating the frontier between the two merged
        !> sublayers
        !
        !>@param size_small_layer
        !> size along the y-direction for the small sublayer (in
        !> terms of the y-direction)
        !--------------------------------------------------------------
        subroutine fix_frontier_after_merge_N(
     $     grdpt_id,
     $     frontier,
     $     size_small_layer)
        
          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpt_id
          integer(ikind)                , intent(in)    :: frontier
          integer(ikind)                , intent(in)    :: size_small_layer

          logical        :: same_size
          integer(ikind) :: i,j,j_fix,j_min,j_max,j_end
          logical        :: j_min_initialized

          same_size         = size_small_layer.eq.size(grdpt_id,2)
          j_end             = min(size_small_layer+1, size(grdpt_id,2)-bc_size+1)
          j_min             = bc_size+1
          j_min_initialized = .true.
          do j=bc_size+1, j_end
             
             if(j_min_initialized.and.(
     $            ((grdpt_id(frontier,j).eq.no_pt).or.
     $            (grdpt_id(frontier+1,j).eq.no_pt)).or.(
     $            j.eq.j_end))
     $            ) then

                if((j.eq.j_end).and.same_size) then
                   j_max = j_end
                else
                   j_max = j-bc_size+1
                end if

                j_fix=j_min-1
                if(j_fix.gt.bc_size) then
                   do i=0,1
                      grdpt_id(frontier+i,j_fix)=bc_interior_pt
                   end do
                end if            


                do j_fix=j_min, j_max-1
                   do i=-1,2
                      grdpt_id(frontier+i,j_fix)=interior_pt
                   end do
                end do

                j_fix=j_max
                do i=0,1
                   grdpt_id(frontier+i,j_fix)=bc_interior_pt
                end do

                j_min_initialized=.false.

             else
                if(.not.j_min_initialized.and.(
     $            (grdpt_id(frontier,j).ne.no_pt).and.
     $            (grdpt_id(frontier+1,j).ne.no_pt))) then
                   j_min=j+bc_size
                   j_min_initialized=.true.
                end if
             end if

          end do

        end subroutine fix_frontier_after_merge_N


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine fixing the frontier between two freshly
        !> merged South sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpt_id
        !> table encapsulating the grid points
        !
        !>@param frontier
        !> index indicating the frontier between the two merged
        !> sublayers
        !
        !>@param j_start
        !> index indicating the starting point when looking for
        !> the frontier to be fixed
        !
        !>@param size_layer
        !> size along the y-direction for the small sublayer (in
        !> terms of the y-direction)
        !--------------------------------------------------------------
        subroutine fix_frontier_after_merge_S(
     $     grdpt_id,
     $     frontier,
     $     j_start,
     $     size_layer)
        
          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpt_id
          integer(ikind)                , intent(in)    :: frontier
          integer(ikind)                , intent(in)    :: j_start
          integer(ikind)                , intent(in)    :: size_layer

          integer(ikind) :: i,j,j_fix,j_min,j_max
          logical        :: j_min_initialized

          j_min=1
          j_min_initialized=.false.
          do j=max(j_start,1), size_layer-bc_size
             
             if(j_min_initialized.and.(
     $            ((grdpt_id(frontier,j).eq.no_pt).or.
     $            (grdpt_id(frontier+1,j).eq.no_pt)).or.
     $            (j.eq.size_layer-bc_size))) then

                if(j.eq.(size_layer-bc_size)) then
                   j_max=size_layer-bc_size
                else
                   j_max = j-(bc_size+1)
                end if

                j_fix=j_min-1
                do i=0,1
                   grdpt_id(frontier+i,j_fix)=bc_interior_pt
                end do

                do j_fix=j_min, j_max
                   do i=-1,2
                      grdpt_id(frontier+i,j_fix)=interior_pt
                   end do
                end do

                j_fix=j_max+1
                if(j_fix.lt.(size_layer-1)) then
                   do i=0,1
                      grdpt_id(frontier+i,j_fix)=bc_interior_pt
                   end do
                end if

                j_min_initialized=.false.

             else
                if(.not.j_min_initialized.and.(
     $            (grdpt_id(frontier,j).ne.no_pt).and.
     $            (grdpt_id(frontier+1,j).ne.no_pt))) then
                   j_min=j+bc_size
                   j_min_initialized=.true.
                end if
             end if

          end do

        end subroutine fix_frontier_after_merge_S


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine fixing the grdpt_id frontier between two
        !> recently merged East sublayers
        !
        !> @date
        !> 07_05_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpt_id
        !> table encapsulating the grid points
        !
        !>@param frontier
        !> index indicating the frontier between the two merged
        !> sublayers
        !
        !>@param size_small_layer
        !> size along the y-direction for the small sublayer (in
        !> terms of the y-direction)
        !--------------------------------------------------------------
        subroutine fix_frontier_after_merge_E(
     $     grdpt_id,
     $     frontier,
     $     size_small_layer)
        
          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpt_id
          integer(ikind)                , intent(in)    :: frontier
          integer(ikind)                , intent(in)    :: size_small_layer

          logical        :: same_size
          integer(ikind) :: i,j,i_fix,i_min,i_max,i_end
          logical        :: i_min_initialized

          same_size         = size_small_layer.eq.size(grdpt_id,1)
          i_end             = min(size_small_layer+1, size(grdpt_id,1)-bc_size+1)
          i_min             = bc_size+1
          i_min_initialized = .true.
          do i=bc_size+1, i_end
             
             if(i_min_initialized.and.(
     $            ((grdpt_id(i,frontier).eq.no_pt).or.
     $            (grdpt_id(i,frontier+1).eq.no_pt)).or.(
     $            i.eq.i_end))
     $            ) then

                if((i.eq.i_end).and.same_size) then
                   i_max = i_end
                else
                   i_max = i-bc_size+1
                end if

                i_fix=i_min-1
                if(i_fix.gt.bc_size) then
                   do j=0,1
                      grdpt_id(i_fix,frontier+j)=bc_interior_pt
                   end do
                end if            

                do j=-1,2
                   do i_fix=i_min, i_max-1
                      grdpt_id(i_fix,frontier+j)=interior_pt
                   end do
                end do

                i_fix=i_max
                do j=0,1
                   grdpt_id(i_fix,frontier+j)=bc_interior_pt
                end do

                i_min_initialized=.false.

             else
                if(.not.i_min_initialized.and.(
     $            (grdpt_id(i,frontier).ne.no_pt).and.
     $            (grdpt_id(i,frontier+1).ne.no_pt))) then
                   i_min=i+bc_size
                   i_min_initialized=.true.
                end if
             end if

          end do

        end subroutine fix_frontier_after_merge_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine fixing the grdpt_id frontier between two
        !> recently merged West sublayers
        !
        !> @date
        !> 08_05_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpt_id
        !> table encapsulating the grid points
        !
        !>@param frontier
        !> index indicating the frontier between the two merged
        !> sublayers
        !
        !>@param i_start
        !> index indicating the starting point when looking for
        !> the frontier to be fixed
        !
        !>@param size_layer
        !> size along the x-direction for the small sublayer (in
        !> terms of the x-direction)
        !--------------------------------------------------------------
        subroutine fix_frontier_after_merge_W(
     $     grdpt_id,
     $     frontier,
     $     i_start,
     $     size_layer)
        
          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpt_id
          integer(ikind)                , intent(in)    :: frontier
          integer(ikind)                , intent(in)    :: i_start
          integer(ikind)                , intent(in)    :: size_layer

          integer(ikind) :: i,j,i_fix,i_min,i_max
          logical        :: i_min_initialized

          i_min=1
          i_min_initialized=.false.
          do i=max(i_start,1), size_layer-bc_size
             
             if(i_min_initialized.and.(
     $            ((grdpt_id(i,frontier).eq.no_pt).or.
     $            (grdpt_id(i,frontier+1).eq.no_pt)).or.
     $            (i.eq.size_layer-bc_size))) then

                if(i.eq.(size_layer-bc_size)) then
                   i_max=size_layer-bc_size
                else
                   i_max = i-(bc_size+1)
                end if

                i_fix=i_min-1
                do j=0,1
                   grdpt_id(i_fix,frontier+j)=bc_interior_pt
                end do

                do j=-1,2
                   do i_fix=i_min, i_max
                      grdpt_id(i_fix,frontier+j)=interior_pt
                   end do
                end do

                i_fix=i_max+1
                if(i_fix.lt.(size_layer-1)) then
                   do j=0,1
                      grdpt_id(i_fix,frontier+j)=bc_interior_pt
                   end do
                end if

                i_min_initialized=.false.

             else
                if(.not.i_min_initialized.and.(
     $            (grdpt_id(i,frontier).ne.no_pt).and.
     $            (grdpt_id(i,frontier+1).ne.no_pt))) then
                   i_min=i+bc_size
                   i_min_initialized=.true.
                end if
             end if

          end do

        end subroutine fix_frontier_after_merge_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine investigating which sublayers correspond
        !> to the subdivision sublayer_first, sublayer_second,
        !> sublayer_third, sublayer_fourth required by the
        !> subroutine merging sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param sublayer1
        !> first sublayer merged
        !
        !>@param sublayer2
        !> second sublayer merged
        !
        !>@param sublayer_first
        !> pointer to the first sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_second
        !> pointer to the second sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_third
        !> pointer to the third sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_fourth
        !> pointer to the fourth sublayer in the convention required
        !> by the subroutine merging sublayers
        !--------------------------------------------------------------
        subroutine get_match_sublayers_N(
     $     sublayer1,
     $     sublayer2,
     $     sublayer_first,
     $     sublayer_second,
     $     sublayer_third,
     $     sublayer_fourth)
        
          implicit none

          type(bf_sublayer), pointer  , intent(in)  :: sublayer1
          type(bf_sublayer), pointer  , intent(in)  :: sublayer2
          type(bf_sublayer), pointer  , intent(out) :: sublayer_first
          type(bf_sublayer), pointer  , intent(out) :: sublayer_second
          type(bf_sublayer), pointer  , intent(out) :: sublayer_third
          type(bf_sublayer), pointer  , intent(out) :: sublayer_fourth

          logical :: sublayer1_is_first
          logical :: sublayer1_is_larger


          ! we divide the new table to be filled into 4 pieces:
          !   ________ ________
          !  |        |        |
          !  |   3    |   4    |
          !  |________|________|
          !  |        |        |
          !  |   1    |   2    |
          !  |________|________|
          !
          ! depending on the relative size and position of the
          ! sublayers, we assign the piece to a certain sublayer

          !< check the relative positions of the two tables
          sublayer1_is_first = sublayer1%element%alignment(1,1).le.
     $         sublayer2%element%alignment(1,1)

          sublayer1_is_larger = size(sublayer1%element%grdpts_id,2).ge.
     $         size(sublayer2%element%grdpts_id,2)

          !< identify the first,second, small and large
          !> sublayers
          if(sublayer1_is_larger) then

             if(sublayer1_is_first) then
                sublayer_first  => sublayer1
                sublayer_second => sublayer2
                sublayer_third  => sublayer1
                nullify(sublayer_fourth)
             else
                sublayer_first  => sublayer2
                sublayer_second => sublayer1
                nullify(sublayer_third)
                sublayer_fourth => sublayer1
             end if
             
          else

             if(sublayer1_is_first) then
                sublayer_first  => sublayer1
                sublayer_second => sublayer2
                nullify(sublayer_third)
                sublayer_fourth => sublayer2
             else
                sublayer_first  => sublayer2
                sublayer_second => sublayer1
                sublayer_third  => sublayer2
                nullify(sublayer_fourth)                
             end if

          end if

        end subroutine get_match_sublayers_N


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine investigating which sublayers correspond
        !> to the subdivision sublayer_first, sublayer_second,
        !> sublayer_third, sublayer_fourth required by the
        !> subroutine merging sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param sublayer1
        !> first sublayer merged
        !
        !>@param sublayer2
        !> second sublayer merged
        !
        !>@param sublayer_first
        !> pointer to the first sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_second
        !> pointer to the second sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_third
        !> pointer to the third sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_fourth
        !> pointer to the fourth sublayer in the convention required
        !> by the subroutine merging sublayers
        !--------------------------------------------------------------
        subroutine get_match_sublayers_S(
     $     sublayer1,
     $     sublayer2,
     $     sublayer_first,
     $     sublayer_second,
     $     sublayer_third,
     $     sublayer_fourth,
     $     j_min1, j_min3, j_max)
        
          implicit none

          type(bf_sublayer), pointer  , intent(in)  :: sublayer1
          type(bf_sublayer), pointer  , intent(in)  :: sublayer2
          type(bf_sublayer), pointer  , intent(out) :: sublayer_first
          type(bf_sublayer), pointer  , intent(out) :: sublayer_second
          type(bf_sublayer), pointer  , intent(out) :: sublayer_third
          type(bf_sublayer), pointer  , intent(out) :: sublayer_fourth
          integer(ikind)              , intent(out) :: j_min1
          integer(ikind)              , intent(out) :: j_min3
          integer(ikind)              , intent(out) :: j_max

          logical :: sublayer1_is_first
          logical :: sublayer1_is_larger


          ! we divide the new table to be filled into 4 pieces:
          !   ________ ________
          !  |        |        |
          !  |   3    |   4    |
          !  |________|________|
          !  |        |        |
          !  |   1    |   2    |
          !  |________|________|
          !
          ! depending on the relative size and position of the
          ! sublayers, we assign the piece to a certain sublayer

          !< check the relative positions of the two tables
          sublayer1_is_first = sublayer1%element%alignment(1,1).le.
     $         sublayer2%element%alignment(1,1)

          sublayer1_is_larger = size(sublayer1%element%grdpts_id,2).ge.
     $         size(sublayer2%element%grdpts_id,2)

          !< identify the first,second, small and large
          !> sublayers
          if(sublayer1_is_larger) then
             j_max = size(sublayer1%element%grdpts_id,2)

             if(sublayer1_is_first) then
                sublayer_first  => sublayer1
                nullify(sublayer_second)
                sublayer_third  => sublayer1
                sublayer_fourth => sublayer2

                j_min1=0
                j_min3=-(j_max-size(sublayer_fourth%element%nodes,2))

             else
                nullify(sublayer_first)
                sublayer_second => sublayer1
                sublayer_third  => sublayer2
                sublayer_fourth => sublayer1

                j_min1=-(j_max-size(sublayer_third%element%nodes,2))
                j_min3=0

             end if
             
          else
             j_max = size(sublayer2%element%grdpts_id,2)

             if(sublayer1_is_first) then
                nullify(sublayer_first)
                sublayer_second => sublayer2
                sublayer_third  => sublayer1
                sublayer_fourth => sublayer2

                j_min1=-(j_max-size(sublayer_third%element%nodes,2))
                j_min3=0

             else
                sublayer_first  => sublayer2
                nullify(sublayer_second)
                sublayer_third  => sublayer2
                sublayer_fourth => sublayer1

                j_min1=0
                j_min3=-(j_max-size(sublayer_fourth%element%nodes,2))

             end if

          end if

        end subroutine get_match_sublayers_S


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine investigating which sublayers correspond
        !> to the subdivision sublayer_first, sublayer_third
        !> required by the subroutine merging sublayers
        !
        !> @date
        !> 06_05_2014 - initial version - J.L. Desmarais
        !
        !>@param sublayer1
        !> first sublayer merged
        !
        !>@param sublayer2
        !> second sublayer merged
        !
        !>@param sublayer_first
        !> pointer to the first sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_third
        !> pointer to the third sublayer in the convention required
        !> by the subroutine merging sublayers
        !--------------------------------------------------------------
        subroutine get_match_sublayers_EW(
     $     sublayer1,
     $     sublayer2,
     $     sublayer_first,
     $     sublayer_third)
        
          implicit none

          type(bf_sublayer), pointer  , intent(in)  :: sublayer1
          type(bf_sublayer), pointer  , intent(in)  :: sublayer2
          type(bf_sublayer), pointer  , intent(out) :: sublayer_first
          type(bf_sublayer), pointer  , intent(out) :: sublayer_third

          logical :: sublayer1_is_first
          logical :: sublayer1_is_larger


          !we divide the new table to be filled into 4 pieces:
          !  ________ ________
          ! |        |        |
          ! |   3    |   4    |
          ! |________|________|
          ! |        |        |
          ! |   1    |   2    |
          ! |________|________|
          !
          !depending on the relative size and position of the
          !sublayers, we assign the piece to a certain sublayer

          !< check the relative positions of the two tables
          sublayer1_is_first = sublayer1%element%alignment(2,1).le.
     $         sublayer2%element%alignment(2,1)

          sublayer1_is_larger = size(sublayer1%element%grdpts_id,1).ge.
     $         size(sublayer2%element%grdpts_id,1)

          !< identify the first,second, small and large
          !> sublayers
          if(sublayer1_is_larger) then

             if(sublayer1_is_first) then
                sublayer_first  => sublayer1
                sublayer_third  => sublayer2

             else
                sublayer_first  => sublayer2
                sublayer_third  => sublayer1

             end if
             
          else

             if(sublayer1_is_first) then
                sublayer_first  => sublayer1
                sublayer_third  => sublayer2

             else
                sublayer_first  => sublayer2
                sublayer_third  => sublayer1

             end if

          end if

        end subroutine get_match_sublayers_EW


        !warning: sublayer1%element%alignment(direction,2) <
        !         sublayer2%element%alignment(direction,1)
        subroutine get_match_direction(
     $     j_min1, j_min2, j_min3, j_min4,
     $     j_max1, j_max3,
     $     interior_j_max1, interior_j_max2, interior_j_max3,
     $     sublayer1, sublayer2,
     $     direction,
     $     alignment,
     $     j_min_first_layer , j_max_first_layer,
     $     j_min_second_layer, j_max_second_layer)

          implicit none

          integer(ikind)                          , intent(out) :: j_min1, j_min2, j_min3, j_min4
          integer(ikind)                          , intent(out) :: j_max1, j_max3
          integer(ikind)                          , intent(out) :: interior_j_max1, interior_j_max2, interior_j_max3
          type(bf_sublayer), pointer              , intent(in)  :: sublayer1, sublayer2
          integer                                 , intent(in)  :: direction
          integer(ikind)                , optional, intent(out) :: j_min_first_layer, j_max_first_layer
          integer(ikind)                , optional, intent(out) :: j_min_second_layer, j_max_second_layer
          integer(ikind), dimension(2,2), optional, intent(in)  :: alignment

          j_max1 = size(sublayer1%element%nodes,direction)
          j_max3 = size(sublayer2%element%nodes,direction)

          if(present(alignment)) then
             interior_j_max1 = max(sublayer1%element%alignment(direction,1)-alignment(direction,1),0)
             interior_j_max3 = max(alignment(direction,2)-sublayer2%element%alignment(direction,2),0)
          else
             interior_j_max1 = 0
             interior_j_max3 = 0
          end if
          interior_j_max2 = max(
     $         sublayer2%element%alignment(direction,1)-
     $         sublayer1%element%alignment(direction,2)-2*bc_size-1,0)

          j_min1 = 0      + interior_j_max1
          j_min2 = j_min1 + j_max1
          j_min3 = j_min2 + interior_j_max2
          j_min4 = j_min3 + j_max3

          if(present(j_min_first_layer)) then
             if(interior_j_max2.gt.0) then
                j_min_first_layer  = 3
                j_max_first_layer  = j_max1-2
                j_min_second_layer = 2
                j_max_second_layer = j_max1-1
             else
                j_min_first_layer  = 1
                j_max_first_layer  = j_max1
                j_min_second_layer = 1
                j_max_second_layer = j_max1
             end if
          end if

        end subroutine get_match_direction

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two nodes tables into one for the
        !> Northern sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param temp_nodes
        !> temporary table where the nodes are copied for the merge
        !
        !>@param nodes
        !> table where the nodes from the interior domain are saved
        !
        !>@param sublayer_first
        !> pointer to the first sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_second
        !> pointer to the second sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_third
        !> pointer to the third sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_fourth
        !> pointer to the fourth sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max3
        !> size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 1)
        !
        !>@param i_max1
        !> size of the block copied from the first sublayer merged
        !>  (see Fig. 1)
        !
        !>@param i_max3
        !> size of the block copied from the second sublayer merged
        !>  (see Fig. 1)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !
        !>@param i_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 1)
        !
        !>@param i_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 1)
        !
        !>@param j_min
        !> index indicating the smaller y-size between the two
        !> nodes tables of the sublayers merged (see Fig. 1)
        !
        !>@param alignment
        !> table indicating the position of the final sublayer after
        !> merging the two sublayers
        !--------------------------------------------------------------
        subroutine merge_nodes_N(
     $     temp_nodes,
     $     sublayer_first,
     $     sublayer_second,
     $     sublayer_third,
     $     sublayer_fourth,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_max1, i_max3,
     $     i_min1, i_min2, i_min3, i_min4,
     $     j_min,
     $     alignment)

          implicit none

          real(rkind), dimension(:,:,:)             , intent(out) :: temp_nodes
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_first
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_second
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_third
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_fourth
          integer(ikind)                            , intent(in)  :: interior_i_max1
          integer(ikind)                            , intent(in)  :: interior_i_max2
          integer(ikind)                            , intent(in)  :: interior_i_max3
          integer(ikind)                            , intent(in)  :: i_max1, i_max3
          integer(ikind)                            , intent(in)  :: i_min1
          integer(ikind)                            , intent(in)  :: i_min2
          integer(ikind)                            , intent(in)  :: i_min3
          integer(ikind)                            , intent(in)  :: i_min4
          integer(ikind)                            , intent(in)  :: j_min
          integer(ikind), dimension(2,2)  , optional, intent(in)  :: alignment


          integer(ikind) :: i,j
          integer        :: k          
          

          do k=1, ne
             do j=1, j_min
                
                do i=1, i_max1
                   temp_nodes(i_min1+i,j,k)=
     $                  sublayer_first%element%nodes(i,j,k)
                end do
                
                do i=1, i_max3
                   temp_nodes(i_min3+i,j,k)=
     $                  sublayer_second%element%nodes(i,j,k)
                end do
                
             end do
             
             if(associated(sublayer_third)) then
                do j=j_min, size(sublayer_third%element%nodes,2)
                   do i=1, i_max1
                      temp_nodes(i_min1+i,j,k)=
     $                     sublayer_third%element%nodes(i,j,k)
                   end do
                end do
             end if
             
             if(associated(sublayer_fourth)) then                   
                do j=j_min, size(sublayer_fourth%element%nodes,2)
                   do i=1, i_max3
                      temp_nodes(i_min3+i,j,k)=
     $                     sublayer_fourth%element%nodes(i,j,k)
                   end do
                end do
             end if
             
          end do

        end subroutine merge_nodes_N


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two nodes tables into one for the
        !> Southern sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param temp_nodes
        !> temporary table where the nodes are copied for the merge
        !
        !>@param nodes
        !> table where the nodes from the interior domain are saved
        !
        !>@param sublayer_first
        !> pointer to the first sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_second
        !> pointer to the second sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_third
        !> pointer to the third sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_fourth
        !> pointer to the fourth sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param interior_i_max1
        !> y-size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max2
        !> y-size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max3
        !> y-size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 1)
        !
        !>@param i_max1
        !> y-size of the block copied from the first sublayer merged
        !>  (see Fig. 1)
        !
        !>@param i_max3
        !> y-size of the block copied from the second sublayer merged
        !>  (see Fig. 1)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !
        !>@param i_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 1)
        !
        !>@param i_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 1)
        !
        !>@param j_min
        !> index indicating the smaller y-size between the two
        !> nodes tables of the sublayers merged (see Fig. 1)
        !
        !>@param j_min1
        !> index indicating the match in y-direction between block2
        !> and the new sublayer (see Fig. 2)
        !
        !>@param j_min3
        !> index indicating the match in y-direction between block4
        !> and the new sublayer (see Fig. 2)
        !
        !>@param j_max
        !> size of the largest sublayer in the y-direction
        !
        !>@param alignment
        !> table indicating the position of the final sublayer after
        !> merging the two sublayers
        !--------------------------------------------------------------
        subroutine merge_nodes_S(
     $     temp_nodes,
     $     sublayer_first,
     $     sublayer_second,
     $     sublayer_third,
     $     sublayer_fourth,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_max1, i_max3,
     $     i_min1, i_min2, i_min3, i_min4,
     $     j_min,
     $     j_min1, j_min3, j_max,
     $     nodes,
     $     alignment)

          implicit none

          real(rkind), dimension(:,:,:)             , intent(out) :: temp_nodes
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_first
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_second
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_third
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_fourth
          integer(ikind)                            , intent(in)  :: interior_i_max1
          integer(ikind)                            , intent(in)  :: interior_i_max2
          integer(ikind)                            , intent(in)  :: interior_i_max3
          integer(ikind)                            , intent(in)  :: i_max1, i_max3
          integer(ikind)                            , intent(in)  :: i_min1
          integer(ikind)                            , intent(in)  :: i_min2
          integer(ikind)                            , intent(in)  :: i_min3
          integer(ikind)                            , intent(in)  :: i_min4
          integer(ikind)                            , intent(in)  :: j_min
          integer(ikind)                            , intent(in)  :: j_min1, j_min3, j_max
          real(rkind), dimension(nx,ny,ne), optional, intent(in)  :: nodes
          integer(ikind), dimension(2,2)  , optional, intent(in)  :: alignment


          integer(ikind) :: i,j
          integer        :: k          
          

          if(
     $         (interior_i_max1.ge.1).or.
     $         (interior_i_max2.ge.1).or.
     $         (interior_i_max3.ge.1)) then

             do k=1, ne

                if(associated(sublayer_first)) then
                   do j=1, j_min
                      do i=1, i_max1
                         temp_nodes(i_min1+i,j,k)=
     $                        sublayer_first%element%nodes(i,j,k)
                      end do
                   end do
                end if

                if(associated(sublayer_second)) then
                   do j=1, j_min
                      do i=1, i_max3
                         temp_nodes(i_min3+i,j,k)=
     $                        sublayer_second%element%nodes(i,j,k)
                      end do
                   end do
                end if

                do j=j_min+1, j_max-2*bc_size

                   do i=1, i_max1
                      temp_nodes(i_min1+i,j,k)=
     $                     sublayer_third%element%nodes(i,j_min1+j,k)
                   end do

                   do i=1, i_max3
                      temp_nodes(i_min3+i,j,k)=
     $                     sublayer_fourth%element%nodes(i,j_min3+j,k)
                   end do

                end do

                do j=j_max-2*bc_size+1, j_max

                   do i=1, interior_i_max1
                      temp_nodes(i,j,k) =
     $                     nodes(alignment(1,1)-(bc_size+1)+i,j-(j_max-2*bc_size),k)
                   end do
                   
                   do i=1, i_max1
                      temp_nodes(i_min1+i,j,k)=
     $                     sublayer_third%element%nodes(i,j_min1+j,k)
                   end do

                   do i=1, interior_i_max2
                      temp_nodes(i_min2+i,j,k)=
     $                     nodes(sublayer_third%element%alignment(1,2)+bc_size+i,j-(j_max-2*bc_size),k)
                   end do

                   do i=1, i_max3
                      temp_nodes(i_min3+i,j,k)=
     $                     sublayer_fourth%element%nodes(i,j_min3+j,k)
                   end do

                   do i=1, interior_i_max3
                      temp_nodes(i_min4+i,j,k) =
     $                     nodes(alignment(1,2)+bc_size-interior_i_max3+i,j-(j_max-2*bc_size),k)
                   end do

                end do

             end do

          !< if no nodes is loaded from the interior
          else

             do k=1, ne

                if(associated(sublayer_first)) then
                   do j=1, j_min
                      do i=1, i_max1
                         temp_nodes(i_min1+i,j,k)=
     $                        sublayer_first%element%nodes(i,j,k)
                      end do
                   end do
                end if

                if(associated(sublayer_second)) then
                   do j=1, j_min
                      do i=1, i_max3
                         temp_nodes(i_min3+i,j,k)=
     $                        sublayer_second%element%nodes(i,j,k)
                      end do
                   end do
                end if

                do j=j_min+1, j_max

                   do i=1, i_max1
                      temp_nodes(i_min1+i,j,k)=
     $                     sublayer_third%element%nodes(i,j_min1+j,k)
                   end do

                   do i=1, i_max3
                      temp_nodes(i_min3+i,j,k)=
     $                     sublayer_fourth%element%nodes(i,j_min3+j,k)
                   end do

                end do

             end do

          end if

        end subroutine merge_nodes_S


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two nodes tables into one for the
        !> Easthern sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param temp_nodes
        !> temporary table where the nodes are copied for the merge
        !
        !>@param nodes
        !> table where the nodes from the interior domain are saved
        !
        !>@param sublayer_first
        !> pointer to the first sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_third
        !> pointer to the third sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param interior_j_max1
        !> y-size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 3)
        !
        !>@param interior_j_max2
        !> y-size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 3)
        !
        !>@param interior_j_max3
        !> y-size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 3)
        !
        !>@param j_max1
        !> y-size of the block copied from the first sublayer merged
        !>  (see Fig. 3)
        !
        !>@param j_max3
        !> y-size of the block copied from the second sublayer merged
        !>  (see Fig. 3)
        !
        !>@param j_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 3)
        !
        !>@param j_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 3)
        !
        !>@param j_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 3)
        !
        !>@param j_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 3)
        !
        !>@param alignment
        !> table indicating the position of the final sublayer after
        !> merging the two sublayers
        !--------------------------------------------------------------
        subroutine merge_nodes_E(
     $     temp_nodes,
     $     sublayer_first,
     $     sublayer_third,
     $     interior_j_max1, interior_j_max2, interior_j_max3,
     $     j_max1, j_max3,
     $     j_min1, j_min2, j_min3, j_min4,
     $     nodes,
     $     alignment)

          implicit none

          real(rkind), dimension(:,:,:)             , intent(out) :: temp_nodes
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_first
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_third
          integer(ikind)                            , intent(in)  :: interior_j_max1
          integer(ikind)                            , intent(in)  :: interior_j_max2
          integer(ikind)                            , intent(in)  :: interior_j_max3
          integer(ikind)                            , intent(in)  :: j_max1, j_max3
          integer(ikind)                            , intent(in)  :: j_min1
          integer(ikind)                            , intent(in)  :: j_min2
          integer(ikind)                            , intent(in)  :: j_min3
          integer(ikind)                            , intent(in)  :: j_min4
          real(rkind), dimension(nx,ny,ne), optional, intent(in)  :: nodes
          integer(ikind), dimension(2,2)  , optional, intent(in)  :: alignment


          integer(ikind) :: i,j
          integer        :: k          
          

          if(
     $         (interior_j_max1.ge.1).or.
     $         (interior_j_max2.ge.1).or.
     $         (interior_j_max3.ge.1)) then

             do k=1, ne
                do j=1, interior_j_max1
                   do i=1,2*bc_size
                      temp_nodes(i,j,k) = nodes(
     $                     nx-2*bc_size+i,
     $                     j+alignment(2,1)-(bc_size+1),
     $                     k)
                   end do                   
                end do
                
                do j=1, j_max1
                   do i=1, size(sublayer_first%element%nodes,1)
                      temp_nodes(i,j_min1+j,k) =
     $                     sublayer_first%element%nodes(i,j,k)
                   end do
                end do
                
                do j=1, interior_j_max2
                   do i=1, 2*bc_size
                      temp_nodes(i,j_min2+j,k) = nodes(
     $                     nx-2*bc_size+i,
     $                     sublayer_first%element%alignment(2,2)+bc_size+j,
     $                     k)
                   end do
                end do
                
                do j=1, j_max3
                   do i=1, size(sublayer_third%element%nodes,1)
                      temp_nodes(i,j_min3+j,k) =
     $                     sublayer_third%element%nodes(i,j,k)
                   end do
                end do
                
                do j=1, interior_j_max3
                   do i=1, 2*bc_size
                      temp_nodes(i,j_min4+j,k) = nodes(
     $                     nx-2*bc_size+i,
     $                     alignment(2,2)+bc_size-interior_j_max3+j,
     $                     k)
                   end do
                end do

             end do

          else

             do k=1, ne
                do j=1, j_max1
                   do i=1, size(sublayer_first%element%nodes,1)
                      temp_nodes(i,j_min1+j,k) =
     $                     sublayer_first%element%nodes(i,j,k)
                   end do
                end do
                
                do j=1, j_max3
                   do i=1, size(sublayer_third%element%nodes,1)
                      temp_nodes(i,j_min3+j,k) =
     $                     sublayer_third%element%nodes(i,j,k)
                   end do
                end do
             end do

          end if
             
        end subroutine merge_nodes_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two nodes tables into one for the
        !> Westhern sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param temp_nodes
        !> temporary table where the nodes are copied for the merge
        !
        !>@param nodes
        !> table where the nodes from the interior domain are saved
        !
        !>@param sublayer_first
        !> pointer to the first sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_third
        !> pointer to the third sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param interior_j_max1
        !> y-size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 3)
        !
        !>@param interior_j_max2
        !> y-size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 3)
        !
        !>@param interior_j_max3
        !> y-size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 3)
        !
        !>@param j_max1
        !> y-size of the block copied from the first sublayer merged
        !>  (see Fig. 3)
        !
        !>@param j_max3
        !> y-size of the block copied from the second sublayer merged
        !>  (see Fig. 3)
        !
        !>@param j_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 3)
        !
        !>@param j_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 3)
        !
        !>@param j_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 3)
        !
        !>@param j_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 3)
        !
        !>@param i_max
        !> x-size of the largest sublayer in the x-direction
        !
        !>@param alignment
        !> table indicating the position of the final sublayer after
        !> merging the two sublayers
        !--------------------------------------------------------------
        subroutine merge_nodes_W(
     $     temp_nodes,
     $     sublayer_first,
     $     sublayer_third,
     $     interior_j_max1, interior_j_max2, interior_j_max3,
     $     j_max1, j_max3,
     $     j_min1, j_min2, j_min3, j_min4,
     $     i_max,
     $     nodes,
     $     alignment)

          implicit none

          real(rkind), dimension(:,:,:)             , intent(out) :: temp_nodes
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_first
          type(bf_sublayer), pointer                , intent(in)  :: sublayer_third
          integer(ikind)                            , intent(in)  :: interior_j_max1
          integer(ikind)                            , intent(in)  :: interior_j_max2
          integer(ikind)                            , intent(in)  :: interior_j_max3
          integer(ikind)                            , intent(in)  :: j_max1, j_max3
          integer(ikind)                            , intent(in)  :: j_min1
          integer(ikind)                            , intent(in)  :: j_min2
          integer(ikind)                            , intent(in)  :: j_min3
          integer(ikind)                            , intent(in)  :: j_min4
          integer(ikind)                            , intent(in)  :: i_max
          real(rkind), dimension(nx,ny,ne), optional, intent(in)  :: nodes
          integer(ikind), dimension(2,2)  , optional, intent(in)  :: alignment


          integer(ikind) :: i,j, i_min
          integer        :: k          
          

          if(
     $         (interior_j_max1.ge.1).or.
     $         (interior_j_max2.ge.1).or.
     $         (interior_j_max3.ge.1)) then

             do k=1, ne
                do j=1, interior_j_max1
                   do i=1,2*bc_size
                      temp_nodes(i_max-2*bc_size+i,j,k) = nodes(
     $                     i,
     $                     j+alignment(2,1)-(bc_size+1),
     $                     k)
                   end do
                end do
                
                i_min= i_max-size(sublayer_first%element%nodes,1)
                do j=1, j_max1
                   do i=i_min+1, i_max
                      temp_nodes(i,j_min1+j,k) =
     $                     sublayer_first%element%nodes(i-i_min,j,k)
                   end do
                end do
                
                do j=1, interior_j_max2
                   do i=1, 2*bc_size
                      temp_nodes(i_max-2*bc_size+i,j_min2+j,k) = nodes(
     $                     i,
     $                     sublayer_first%element%alignment(2,2)+bc_size+j,
     $                     k)
                   end do
                end do
                
                i_min= i_max-size(sublayer_third%element%nodes,1)
                do j=1, j_max3
                   do i=i_min+1, i_max
                      temp_nodes(i,j_min3+j,k) =
     $                     sublayer_third%element%nodes(i-i_min,j,k)
                   end do
                end do
                
                do j=1, interior_j_max3
                   do i=1, 2*bc_size
                      temp_nodes(i_max-2*bc_size+i,j_min4+j,k) = nodes(
     $                     i,
     $                     alignment(2,2)+bc_size-interior_j_max3+j,
     $                     k)
                   end do
                end do

             end do

          else

             do k=1, ne
                
                i_min= i_max-size(sublayer_first%element%nodes,1)
                do j=1, j_max1
                   do i=i_min+1, i_max
                      temp_nodes(i,j_min1+j,k) =
     $                     sublayer_first%element%nodes(i-i_min,j,k)
                   end do
                end do
                
                i_min= i_max-size(sublayer_third%element%nodes,1)
                do j=1, j_max3
                   do i=i_min+1, i_max
                      temp_nodes(i,j_min3+j,k) =
     $                     sublayer_third%element%nodes(i-i_min,j,k)
                   end do
                end do
                
             end do

          end if
             
        end subroutine merge_nodes_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two nodes tables into one for the
        !> Northern sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param temp_grdptid
        !> temporary table where the gridpoint ID are copied for the
        !> merge
        !
        !>@param sublayer_first
        !> pointer to the first sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_second
        !> pointer to the second sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_third
        !> pointer to the third sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_fourth
        !> pointer to the fourth sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max3
        !> size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 1)
        !
        !>@param i_max1
        !> size of the block copied from the first sublayer merged
        !>  (see Fig. 1)
        !
        !>@param i_max3
        !> size of the block copied from the second sublayer merged
        !>  (see Fig. 1)
        !
        !>@param i_min_first_layer
        !> index indicating the lower border for the first layer
        !> above the exchange layer
        !>  (see Fig. 5)
        !
        !>@param i_max_first_layer
        !> index indicating the upper border for the first layer
        !> above the exchange layer
        !>  (see Fig. 5)
        !
        !>@param i_min_second_layer
        !> index indicating the lower border for the second layer
        !> above the exchange layer
        !>  (see Fig. 6)
        !
        !>@param i_max_second_layer
        !> index indicating the upper border for the second layer
        !> above the exchange layer
        !>  (see Fig. 6)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !
        !>@param i_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 1)
        !
        !>@param i_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 1)
        !
        !>@param j_min
        !> index indicating the smaller y-size between the two
        !> nodes tables of the sublayers merged (see Fig. 1)
        !--------------------------------------------------------------
        subroutine merge_grdpts_id_N(
     $     temp_grdptid,
     $     sublayer_first,
     $     sublayer_second,
     $     sublayer_third,
     $     sublayer_fourth,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_max1, i_max3,
     $     i_min_first_layer, i_max_first_layer,
     $     i_min_second_layer, i_max_second_layer,
     $     i_min1, i_min2, i_min3, i_min4,
     $     j_min,
     $     neighbor_E, neighbor_W)

          implicit none

          integer, dimension(:,:)   , intent(out) :: temp_grdptid
          type(bf_sublayer), pointer, intent(in)  :: sublayer_first
          type(bf_sublayer), pointer, intent(in)  :: sublayer_second
          type(bf_sublayer), pointer, intent(in)  :: sublayer_third
          type(bf_sublayer), pointer, intent(in)  :: sublayer_fourth
          integer(ikind)            , intent(in)  :: interior_i_max1
          integer(ikind)            , intent(in)  :: interior_i_max2
          integer(ikind)            , intent(in)  :: interior_i_max3
          integer(ikind)            , intent(in)  :: i_max1, i_max3
          integer(ikind)            , intent(in)  :: i_min_first_layer
          integer(ikind)            , intent(in)  :: i_max_first_layer
          integer(ikind)            , intent(in)  :: i_min_second_layer
          integer(ikind)            , intent(in)  :: i_max_second_layer
          integer(ikind)            , intent(in)  :: i_min1
          integer(ikind)            , intent(in)  :: i_min2
          integer(ikind)            , intent(in)  :: i_min3
          integer(ikind)            , intent(in)  :: i_min4
          integer(ikind)            , intent(in)  :: j_min
          logical                   , intent(in)  :: neighbor_E
          logical                   , intent(in)  :: neighbor_W


          integer(ikind) :: i,j
          integer(ikind) :: neighbor_i_min, neighbor_i_max


          if(
     $         (interior_i_max1.ge.1).or.
     $         (interior_i_max2.ge.1).or.
     $         (interior_i_max3.ge.1)) then

             !exchange layer
             call add_exchange_pt_layer_NS(
     $            temp_grdptid,
     $            sublayer_first%element%grdpts_id,
     $            sublayer_second%element%grdpts_id,
     $            1, bc_size,
     $            interior_i_max1,
     $            interior_i_max2,
     $            interior_i_max3,
     $            i_min1, i_max1, 0,
     $            i_min2,
     $            i_min3, i_max3, 0,
     $            i_min4)
             
             !interior layer = first layer
             call add_interior_pt_layer_NS(
     $            temp_grdptid,
     $            sublayer_first%element%grdpts_id,
     $            sublayer_second%element%grdpts_id,
     $            bc_size+1,
     $            interior_i_max1,
     $            interior_i_max2,
     $            interior_i_max3,
     $            i_min_first_layer,
     $            i_max_first_layer,
     $            i_min1, 0,
     $            i_min2,
     $            i_min3, i_max3, 0,
     $            i_min4,
     $            neighbor_E, neighbor_W)

             !bc_interior = second layer
             call add_bc_interior_pt_layer_NS(
     $            temp_grdptid,
     $            sublayer_first%element%grdpts_id,
     $            sublayer_second%element%grdpts_id,
     $            bc_size+2,
     $            interior_i_max1,
     $            interior_i_max2,
     $            interior_i_max3,
     $            i_min_second_layer,
     $            i_max_second_layer,
     $            i_min1, 0,
     $            i_min2,
     $            i_min3, i_max3, 0,
     $            i_min4,
     $            neighbor_E, neighbor_W)

             !bc layer = third layer
             call add_bc_pt_layer_NS(
     $            temp_grdptid,
     $            sublayer_first%element%grdpts_id,
     $            sublayer_second%element%grdpts_id,
     $            bc_size+3,
     $            interior_i_max1,
     $            interior_i_max2,
     $            interior_i_max3,
     $            i_min1, i_max1, 0,
     $            i_min2,
     $            i_min3, i_max3, 0,
     $            i_min4,
     $            neighbor_E, neighbor_W)

             !no_pt = final layer
             call add_no_pt_layer_NS(
     $            temp_grdptid,
     $            sublayer_first%element%grdpts_id,
     $            sublayer_second%element%grdpts_id,
     $            bc_size+4, j_min,
     $            interior_i_max1,
     $            interior_i_max2,
     $            interior_i_max3,
     $            i_min1, i_max1, 0,
     $            i_min2,
     $            i_min3, i_max3, 0,
     $            i_min4)
             
             !copy from the grdpt tables
             if(associated(sublayer_third)) then
                do j=j_min+1, size(sublayer_third%element%grdpts_id,2)
                   do i=1, interior_i_max1
                      temp_grdptid(i,j)=
     $                     no_pt
                   end do
                   
                   do i=1, i_max1
                      temp_grdptid(i_min1+i,j)=
     $                     sublayer_first%element%grdpts_id(i,j)
                   end do

                   do i=1, interior_i_max2
                      temp_grdptid(i_min2+i,j) =
     $                     no_pt
                   end do
                
                   do i=1, i_max3+interior_i_max3
                      temp_grdptid(i_min3+i,j)= no_pt
                   end do
                end do
             end if
             
             if(associated(sublayer_fourth)) then                   
                do j=j_min+1, size(sublayer_fourth%element%grdpts_id,2)
                   do i=1, interior_i_max1+i_max1
                      temp_grdptid(i,j) = no_pt
                   end do

                   do i=1, interior_i_max2
                      temp_grdptid(i_min2+i,j) =
     $                     no_pt
                   end do

                   do i=1, i_max3
                      temp_grdptid(i_min3+i,j)= 
     $                     sublayer_fourth%element%grdpts_id(i,j)
                   end do

                   do i=1, interior_i_max3
                      temp_grdptid(i_min4+i,j) = no_pt
                   end do
                end do
             end if
          
          !< if no gridpoint ID need to be added
          else
             if(neighbor_E.or.neighbor_W) then

                if(neighbor_W) then
                   neighbor_i_max= bc_size
                else
                   neighbor_i_max= 0
                end if

                if(neighbor_E) then
                   neighbor_i_min= i_max3-bc_size
                else
                   neighbor_i_min= i_max3
                end if

                do j=1, 2*bc_size+1
                   do i=1, neighbor_i_max
                      temp_grdptid(i_min1+i,j)= exchange_pt
                   end do
                   
                   do i=neighbor_i_max+1, i_max1
                      temp_grdptid(i_min1+i,j)=
     $                     sublayer_first%element%grdpts_id(i,j)
                   end do

                   do i=1, neighbor_i_min
                      temp_grdptid(i_min2+i,j)=
     $                     sublayer_second%element%grdpts_id(i,j)
                   end do

                   do i=neighbor_i_min+1, i_max3
                      temp_grdptid(i_min2+i,j)= exchange_pt
                   end do

                end do

                do j=2*bc_size+2, j_min
                   do i=1, i_max1
                      temp_grdptid(i_min1+i,j)=
     $                     sublayer_first%element%grdpts_id(i,j)
                   end do
                   
                   do i=1, i_max3
                      temp_grdptid(i_min2+i,j)=
     $                     sublayer_second%element%grdpts_id(i,j)
                   end do
                end do

             else
                do j=1, j_min
                   do i=1, i_max1
                      temp_grdptid(i_min1+i,j)=
     $                     sublayer_first%element%grdpts_id(i,j)
                   end do
                   
                   do i=1, i_max3
                      temp_grdptid(i_min2+i,j)=
     $                     sublayer_second%element%grdpts_id(i,j)
                   end do
                end do
             end if
             
             if(associated(sublayer_third)) then
                do j=j_min+1, size(sublayer_third%element%grdpts_id,2)
                   do i=1, i_max1
                      temp_grdptid(i_min1+i,j)=
     $                     sublayer_first%element%grdpts_id(i,j)
                   end do
                
                   do i=1, i_max3
                      temp_grdptid(i_min2+i,j)= no_pt
                   end do
                end do
             end if
             
             if(associated(sublayer_fourth)) then                   
                do j=j_min+1, size(sublayer_fourth%element%grdpts_id,2)
                   do i=1, i_max1
                      temp_grdptid(i,j) = no_pt
                   end do

                   do i=1, i_max3
                      temp_grdptid(i_min2+i,j)= 
     $                     sublayer_fourth%element%grdpts_id(i,j)
                   end do
                end do
             end if
          end if

        end subroutine merge_grdpts_id_N


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two nodes tables into one for the
        !> Southern sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param temp_grdptid
        !> temporary table where the gridpoint ID are copied for the
        !> merge
        !
        !>@param sublayer_first
        !> pointer to the first sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_second
        !> pointer to the second sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_third
        !> pointer to the third sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_fourth
        !> pointer to the fourth sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max3
        !> size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 1)
        !
        !>@param i_max1
        !> size of the block copied from the first sublayer merged
        !>  (see Fig. 1)
        !
        !>@param i_max3
        !> size of the block copied from the second sublayer merged
        !>  (see Fig. 1)
        !
        !>@param i_min_first_layer
        !> index indicating the lower border for the first layer
        !> above the exchange layer
        !>  (see Fig. 5)
        !
        !>@param i_max_first_layer
        !> index indicating the upper border for the first layer
        !> above the exchange layer
        !>  (see Fig. 5)
        !
        !>@param i_min_second_layer
        !> index indicating the lower border for the second layer
        !> above the exchange layer
        !>  (see Fig. 6)
        !
        !>@param i_max_second_layer
        !> index indicating the upper border for the second layer
        !> above the exchange layer
        !>  (see Fig. 6)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !
        !>@param i_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 1)
        !
        !>@param i_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 1)
        !
        !>@param j_min
        !> index indicating the smaller y-size between the two
        !> nodes tables of the sublayers merged (see Fig. 1)
        !
        !>@param j_min1
        !> index indicating the match in y-direction between block2
        !> and the new sublayer (see Fig. 2)
        !
        !>@param j_min3
        !> index indicating the match in y-direction between block4
        !> and the new sublayer (see Fig. 2)
        !
        !>@param j_max
        !> size of the largest sublayer in the y-direction
        !--------------------------------------------------------------
        subroutine merge_grdpts_id_S(
     $     temp_grdptid,
     $     sublayer_first,
     $     sublayer_second,
     $     sublayer_third,
     $     sublayer_fourth,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_max1, i_max3,
     $     i_min_first_layer, i_max_first_layer,
     $     i_min_second_layer, i_max_second_layer,
     $     i_min1, i_min2, i_min3, i_min4,
     $     j_min, j_min1, j_min3, j_max,
     $     neighbor_E,
     $     neighbor_W)

          implicit none

          integer, dimension(:,:)   , intent(out) :: temp_grdptid
          type(bf_sublayer), pointer, intent(in)  :: sublayer_first
          type(bf_sublayer), pointer, intent(in)  :: sublayer_second
          type(bf_sublayer), pointer, intent(in)  :: sublayer_third
          type(bf_sublayer), pointer, intent(in)  :: sublayer_fourth
          integer(ikind)            , intent(in)  :: interior_i_max1
          integer(ikind)            , intent(in)  :: interior_i_max2
          integer(ikind)            , intent(in)  :: interior_i_max3
          integer(ikind)            , intent(in)  :: i_max1, i_max3
          integer(ikind)            , intent(in)  :: i_min_first_layer
          integer(ikind)            , intent(in)  :: i_max_first_layer
          integer(ikind)            , intent(in)  :: i_min_second_layer
          integer(ikind)            , intent(in)  :: i_max_second_layer
          integer(ikind)            , intent(in)  :: i_min1
          integer(ikind)            , intent(in)  :: i_min2
          integer(ikind)            , intent(in)  :: i_min3
          integer(ikind)            , intent(in)  :: i_min4
          integer(ikind)            , intent(in)  :: j_min
          integer(ikind)            , intent(in)  :: j_min1, j_min3, j_max
          logical                   , intent(in)  :: neighbor_E
          logical                   , intent(in)  :: neighbor_W

          integer(ikind) :: i,j
          integer(ikind) :: neighbor_i_min, neighbor_i_max
          
          if(
     $         (interior_i_max1.ge.1).or.
     $         (interior_i_max2.ge.1).or.
     $         (interior_i_max3.ge.1)) then

             if(associated(sublayer_first)) then
                do j=1, j_min
                   do i=1, interior_i_max1
                      temp_grdptid(i,j)=
     $                     no_pt
                   end do
                   
                   do i=1, i_max1
                      temp_grdptid(i_min1+i,j)=
     $                     sublayer_first%element%grdpts_id(i,j)
                   end do
                   
                   do i=1, interior_i_max2
                      temp_grdptid(i_min2+i,j) =
     $                     no_pt
                   end do

                   do i=1, i_max3+interior_i_max3
                      temp_grdptid(i_min3+i,j)= no_pt
                   end do

                end do
             end if

             if(associated(sublayer_second)) then
                do j=1, j_min
                   do i=1, interior_i_max1+i_max1
                      temp_grdptid(i,j) = no_pt
                   end do
                   
                   do i=1, interior_i_max2
                      temp_grdptid(i_min2+i,j) =
     $                     no_pt
                   end do
                   
                   do i=1, i_max3
                      temp_grdptid(i_min3+i,j)= 
     $                     sublayer_second%element%grdpts_id(i,j)
                   end do
                   
                   do i=1, interior_i_max3
                      temp_grdptid(i_min4+i,j) = no_pt
                   end do
                end do
             end if

             !no_pt layer
             call add_no_pt_layer_NS(
     $            temp_grdptid,
     $            sublayer_third%element%grdpts_id,
     $            sublayer_fourth%element%grdpts_id,
     $            j_min+1, j_max-(2*bc_size+1),
     $            interior_i_max1,
     $            interior_i_max2,
     $            interior_i_max3,
     $            i_min1, i_max1, j_min1,
     $            i_min2,
     $            i_min3, i_max3, j_min3,
     $            i_min4)

             !bc_pt layer
             call add_bc_pt_layer_NS(
     $            temp_grdptid,
     $            sublayer_third%element%grdpts_id,
     $            sublayer_fourth%element%grdpts_id,
     $            j_max-(2*bc_size),
     $            interior_i_max1,
     $            interior_i_max2,
     $            interior_i_max3,
     $            i_min1, i_max1, j_min1,
     $            i_min2,
     $            i_min3, i_max3, j_min3,
     $            i_min4,
     $            neighbor_E,
     $            neighbor_W)
             
             !bc_interior_pt layer
             call add_bc_interior_pt_layer_NS(
     $            temp_grdptid,
     $            sublayer_third%element%grdpts_id,
     $            sublayer_fourth%element%grdpts_id,
     $            j_max-(2*bc_size)+1,
     $            interior_i_max1,
     $            interior_i_max2,
     $            interior_i_max3,
     $            i_min_second_layer,
     $            i_max_second_layer,
     $            i_min1, j_min1,
     $            i_min2,
     $            i_min3, i_max3, j_min3,
     $            i_min4,
     $            neighbor_E,
     $            neighbor_W)

             !interior_pt layer
             call add_interior_pt_layer_NS(
     $            temp_grdptid,
     $            sublayer_third%element%grdpts_id,
     $            sublayer_fourth%element%grdpts_id,
     $            j_max-bc_size,
     $            interior_i_max1,
     $            interior_i_max2,
     $            interior_i_max3,
     $            i_min_first_layer,
     $            i_max_first_layer,
     $            i_min1, j_min1,
     $            i_min2,
     $            i_min3, i_max3, j_min3,
     $            i_min4,
     $            neighbor_E,
     $            neighbor_W)

             !exchange layer
             call add_exchange_pt_layer_NS(
     $            temp_grdptid,
     $            sublayer_third%element%grdpts_id,
     $            sublayer_fourth%element%grdpts_id,
     $            j_max-bc_size+1, j_max,
     $            interior_i_max1,
     $            interior_i_max2,
     $            interior_i_max3,
     $            i_min1, i_max1, j_min1,
     $            i_min2,
     $            i_min3, i_max3, j_min3,
     $            i_min4)

          !< if there are no gridpt ID to be artificially added
          !> due to a larger alignment than the original merged
          !> sublayers
          else

             if(associated(sublayer_first)) then
                do j=1, j_min
                   do i=1, i_max1
                      temp_grdptid(i_min1+i,j)=
     $                     sublayer_first%element%grdpts_id(i,j_min1+j)
                   end do
                
                   do i=1, i_max3
                      temp_grdptid(i_min2+i,j)= no_pt
                   end do
                end do
             end if
             
             if(associated(sublayer_second)) then
                do j=1, j_min
                   do i=1, i_max1
                      temp_grdptid(i,j) = no_pt
                   end do

                   do i=1, i_max3
                      temp_grdptid(i_min2+i,j)= 
     $                     sublayer_second%element%grdpts_id(i,j)
                   end do
                end do
             end if

             if(neighbor_E.or.neighbor_W) then

                if(neighbor_W) then
                   neighbor_i_max= bc_size
                else
                   neighbor_i_max= 0
                end if

                if(neighbor_E) then
                   neighbor_i_min= i_max3-bc_size
                else
                   neighbor_i_min= i_max3
                end if

                do j=j_min+1, j_max-(2*bc_size+1)
             
                   do i=1, i_max1
                      temp_grdptid(i_min1+i,j) =
     $                     sublayer_third%element%grdpts_id(i,j_min1+j)
                   end do
                   
                   do i=1, i_max3
                      temp_grdptid(i_min3+i,j) =
     $                     sublayer_fourth%element%grdpts_id(i,j_min3+j)
                   end do
                   
                end do

                do j=j_max-(2*bc_size+1)+1, j_max
             
                   do i=1, neighbor_i_max
                      temp_grdptid(i_min1+i,j) = exchange_pt
                   end do

                   do i=neighbor_i_max+1, i_max1
                      temp_grdptid(i_min1+i,j) =
     $                     sublayer_third%element%grdpts_id(i,j_min1+j)
                   end do
                   
                   do i=1, neighbor_i_min
                      temp_grdptid(i_min3+i,j) =
     $                     sublayer_fourth%element%grdpts_id(i,j_min3+j)
                   end do

                   do i=neighbor_i_min+1, i_max3
                      temp_grdptid(i_min3+i,j) = exchange_pt
                   end do
                   
                end do
             else

                do j=j_min+1, j_max
             
                   do i=1, i_max1
                      temp_grdptid(i_min1+i,j) =
     $                     sublayer_third%element%grdpts_id(i,j_min1+j)
                   end do
                   
                   do i=1, i_max3
                      temp_grdptid(i_min3+i,j) =
     $                     sublayer_fourth%element%grdpts_id(i,j_min3+j)
                   end do
                   
                end do
             end if
             
          end if

        end subroutine merge_grdpts_id_S


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two nodes tables into one for the
        !> Easthern sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param temp_grdptid
        !> temporary table where the gridpoint ID are copied for the
        !> merge
        !
        !>@param sublayer_first
        !> pointer to the first sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_third
        !> pointer to the third sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param interior_j_max1
        !> y-size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 3)
        !
        !>@param interior_j_max2
        !> y-size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 3)
        !
        !>@param interior_j_max3
        !> y-size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 3)
        !
        !>@param j_max1
        !> y-size of the block copied from the first sublayer merged
        !>  (see Fig. 3)
        !
        !>@param j_max3
        !> y-size of the block copied from the second sublayer merged
        !>  (see Fig. 3)
        !
        !>@param j_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 3)
        !
        !>@param j_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 3)
        !
        !>@param j_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 3)
        !
        !>@param j_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 3)
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine merge_grdpts_id_E(
     $     temp_grdptid,
     $     sublayer_first,
     $     sublayer_third,
     $     interior_j_max1, interior_j_max2, interior_j_max3,
     $     j_max1, j_max3,
     $     j_min1, j_min2, j_min3, j_min4,
     $     i_max,
     $     neighbor_N,
     $     neighbor_S)

          implicit none

          integer, dimension(:,:)    , intent(out) :: temp_grdptid
          type(bf_sublayer), pointer , intent(in)  :: sublayer_first
          type(bf_sublayer), pointer , intent(in)  :: sublayer_third
          integer(ikind)             , intent(in)  :: interior_j_max1
          integer(ikind)             , intent(in)  :: interior_j_max2
          integer(ikind)             , intent(in)  :: interior_j_max3
          integer(ikind)             , intent(in)  :: j_max1, j_max3
          integer(ikind)             , intent(in)  :: j_min1
          integer(ikind)             , intent(in)  :: j_min2
          integer(ikind)             , intent(in)  :: j_min3
          integer(ikind)             , intent(in)  :: j_min4
          integer(ikind)             , intent(in)  :: i_max
          logical                    , intent(in)  :: neighbor_N
          logical                    , intent(in)  :: neighbor_S

          integer(ikind) :: i,j, j_min_block, j_max_block
          
          !block1: nodes copied from the interior domain
          !        frontier between block1 and block2
          if(interior_j_max1.ge.2) then
             
             !> grdpt_id corresponding to nodes copied
             !> from the interior domain: 1st block
             
             !1st block : j=1, bc_size
             if(neighbor_S) then
                do j=1, bc_size
                   do i=1, 2*bc_size+1
                      temp_grdptid(i,j) = exchange_pt
                   end do
                   do i=2*bc_size+2, i_max
                      temp_grdptid(i,j) = no_pt
                   end do
                end do
             else
                call add_layer1_out_E(
     $               1, temp_grdptid, i_max)
                
                call add_layer2_out_E(
     $               2, temp_grdptid, i_max)
             end if
             
             !1st block : j>bc_size
             call add_layer_interior_E(
     $            3, interior_j_max1, temp_grdptid, i_max)

             !> frontier between the 1st and 2nd blocks
             call add_layer1_in_E(
     $            interior_j_max1+1, temp_grdptid,
     $            1, sublayer_first%element%grdpts_id,
     $            i_max)

             call add_layer1_in_E(
     $            interior_j_max1+2, temp_grdptid,
     $            2, sublayer_first%element%grdpts_id,
     $            i_max)

             !> update the lower border for the next block copied
             j_min_block = bc_size+1

          else

             if(interior_j_max1.eq.1) then
                
                !1st block : j=1, bc_size
                if(neighbor_S) then
                   j=1
                   do i=1, 2*bc_size+1
                      temp_grdptid(i,j) = exchange_pt
                   end do
                   do i=2*bc_size+2, i_max
                      temp_grdptid(i,j) = no_pt
                   end do

                   j=2
                   do i=1, 2*bc_size+1
                      temp_grdptid(i,j) = exchange_pt
                   end do

                   do i=2*bc_size+2, size(sublayer_first%element%grdpts_id,1)
                      temp_grdptid(i,j) =
     $                     sublayer_first%element%grdpts_id(i,j-1)
                   end do

                   do i=size(sublayer_first%element%grdpts_id,1)+1, i_max
                      temp_grdptid(i,j) = no_pt
                   end do

                   
                else

                   !> grdpt_id corresponding to nodes copied
                   !> from the interior domain: 1st block
                   call add_layer1_out_E(
     $                  1, temp_grdptid, i_max)
                   
                   !> frontier between the 1st and 2nd blocks
                   call add_layer2_in_E(
     $                  2, temp_grdptid,
     $                  1, sublayer_first%element%grdpts_id,
     $                  i_max)
                end if                
                
                call add_layer1_in_E(
     $               3, temp_grdptid,
     $               2, sublayer_first%element%grdpts_id,
     $               i_max)

                !> update the lower border for the next block copied
                j_min_block = bc_size+1
                
             else
                
                if(neighbor_S) then
                   
                   do j=1, bc_size
                      do i=1, 2*bc_size+1
                         temp_grdptid(i,j) = exchange_pt
                      end do

                      do i=2*bc_size+2, size(sublayer_first%element%grdpts_id,1)
                         temp_grdptid(i,j) =
     $                        sublayer_first%element%grdpts_id(i,j)
                      end do
                      
                      do i=size(sublayer_first%element%grdpts_id,1)+1, i_max
                         temp_grdptid(i,j) = no_pt
                      end do                      
                      
                   end do

                   !> update the lower border for the next block copied
                   j_min_block = bc_size+1
                else

                   !> update the lower border for the next block copied
                   j_min_block = bc_size+1
                end if

             end if
          end if
                

          !block2: grdpt copied from the sublayer_first
          !        frontier between block2 and block3
          !block3: grdpt corresponding to nodes copied
          !        from interior to the buffer layer
          !        frontier between block3 and block4
          if(interior_j_max2.ge.1) then
             
             !> update the upper border for block2 copied
             j_max_block = j_max1-bc_size

             !> copy the block corresponding to sublayer_
             !> first
             do j=j_min_block, j_max_block
                do i=1, size(sublayer_first%element%grdpts_id,1)
                   temp_grdptid(i,j_min1+j) = sublayer_first%element%grdpts_id(i,j)
                end do
                
                do i=size(sublayer_first%element%grdpts_id,1)+1, i_max
                   temp_grdptid(i,j_min1+j) = no_pt
                end do
             end do

             !> frontier between block2 and block3
             call add_layer1_in_E(
     $            j_min1+j_max1-1, temp_grdptid,
     $            j_max_block+1, sublayer_first%element%grdpts_id,
     $            i_max)

             call add_layer1_in_E(
     $            j_min1+j_max1, temp_grdptid,
     $            j_max_block+2, sublayer_first%element%grdpts_id,
     $            i_max)

             call add_layer_interior_E(
     $            j_min2+1, j_min2+interior_j_max2, temp_grdptid, i_max)

             !> frontier between the 3rd and 4th blocks
             call add_layer1_in_E(
     $            j_min3+1, temp_grdptid,
     $            1, sublayer_third%element%grdpts_id,
     $            i_max)

             call add_layer1_in_E(
     $            j_min3+bc_size, temp_grdptid,
     $            2, sublayer_third%element%grdpts_id,
     $            i_max)

             !> update the lower border for the next block copied
             j_min_block = bc_size+1

          else

             !> update the upper border for block2 copied
             j_max_block = j_max1

             !> copy the block corresponding to sublayer_
             !> first
             do j=j_min_block, j_max_block
                do i=1, size(sublayer_first%element%grdpts_id,1)
                   temp_grdptid(i,j_min1+j) = sublayer_first%element%grdpts_id(i,j)
                end do
                
                do i=size(sublayer_first%element%grdpts_id,1)+1, i_max
                   temp_grdptid(i,j_min1+j) = no_pt
                end do
             end do

             !> update the lower border for the next block copied
             j_min_block = 1

          end if

          
          !block4: grdpt copied from the sublayer_third
          !        frontier between block4 and block5
          !block5: grdpt corresponding to nodes copied
          !        from interior to the buffer layer
          !        frontier between block5 and the edge
          if(interior_j_max3.ge.2) then

             !> update the upper border for block2 copied
             j_max_block = j_max3-bc_size

             !> copy the block corresponding to sublayer_
             !> first
             do j=j_min_block, j_max_block
                do i=1, size(sublayer_third%element%grdpts_id,1)
                   temp_grdptid(i,j_min3+j) = sublayer_third%element%grdpts_id(i,j)
                end do
                
                do i=size(sublayer_third%element%grdpts_id,1)+1, i_max
                   temp_grdptid(i,j_min3+j) = no_pt
                end do
             end do

             !> frontier between block4 and block5
             call add_layer1_in_E(
     $            j_min3+j_max3-1, temp_grdptid,
     $            j_max_block+1, sublayer_third%element%grdpts_id,
     $            i_max)

             call add_layer1_in_E(
     $            j_min3+j_max3, temp_grdptid,
     $            j_max_block+2, sublayer_third%element%grdpts_id,
     $            i_max)

             call add_layer_interior_E(
     $            j_min4+1, j_min4+interior_j_max3-bc_size, temp_grdptid, i_max)

             !> frontier at the edge
             if(neighbor_N) then
                do j=j_min4+interior_j_max3-1, j_min4+interior_j_max3
                   do i=1, 2*bc_size+1
                      temp_grdptid(i,j) = exchange_pt
                   end do
                   do i=2*bc_size+2, i_max
                      temp_grdptid(i,j) = no_pt
                   end do
                end do
             else
                call add_layer2_out_E(
     $               j_min4+interior_j_max3-1, temp_grdptid, i_max)
                
                call add_layer1_out_E(
     $               j_min4+interior_j_max3  , temp_grdptid, i_max)
             end if
          else

             if(interior_j_max3.eq.1) then
                
                !> update the upper border for block2 copied
                j_max_block = j_max3-bc_size

                !> copy the block corresponding to sublayer_
                !> first
                do j=j_min_block, j_max_block
                   do i=1, size(sublayer_third%element%grdpts_id,1)
                      temp_grdptid(i,j_min3+j) = sublayer_third%element%grdpts_id(i,j)
                   end do
                   
                   do i=size(sublayer_third%element%grdpts_id,1)+1, i_max
                      temp_grdptid(i,j_min3+j) = no_pt
                   end do
                end do

                !> frontier between block4 and the edge
                call add_layer1_in_E(
     $               j_min3+j_max_block+1, temp_grdptid,
     $               j_max_block+1, sublayer_third%element%grdpts_id,
     $               i_max)

                !block5: last two layers
                if(neighbor_N) then

                   j=j_min3+j_max_block+2

                   do i=1, 2*bc_size+1
                      temp_grdptid(i,j) = exchange_pt
                   end do
                   
                   do i=2*bc_size+2, size(sublayer_third%element%grdpts_id,1)
                      temp_grdptid(i,j) = 
     $                     sublayer_third%element%grdpts_id(i,j_max_block+2)
                   end do
                   
                   do i=size(sublayer_third%element%grdpts_id,1)+1, i_max
                      temp_grdptid(i,j) = no_pt
                   end do

                   j=j_min3+j_max_block+3

                   do i=1, 2*bc_size+1
                      temp_grdptid(i,j) = exchange_pt
                   end do

                   do i=2*bc_size+2, i_max
                      temp_grdptid(i,j) = no_pt
                   end do

                else

                   call add_layer2_in_E(
     $                  j_min3+j_max_block+2, temp_grdptid,
     $                  j_max_block+2, sublayer_third%element%grdpts_id,
     $                  i_max)

                   call add_layer1_out_E(
     $                  j_min3+j_max_block+3, temp_grdptid, i_max)
                end if

             else

                !> update the upper border for block2 copied
                j_max_block = j_max3

                !> copy the block corresponding to sublayer_
                !> first
                if(neighbor_N) then
                   
                   do j=j_min_block, j_max_block-bc_size
                      do i=1, size(sublayer_third%element%grdpts_id,1)
                         temp_grdptid(i,j_min3+j) =
     $                        sublayer_third%element%grdpts_id(i,j)
                      end do
                      
                      do i=size(sublayer_third%element%grdpts_id,1)+1, i_max
                         temp_grdptid(i,j_min3+j) = no_pt
                      end do
                   end do

                   do j=j_max_block-bc_size+1, j_max_block
                      do i=1, 2*bc_size+1
                         temp_grdptid(i,j_min3+j) = exchange_pt
                      end do
                      
                      do i=2*bc_size+2, size(sublayer_third%element%grdpts_id,1)
                         temp_grdptid(i,j_min3+j) =
     $                        sublayer_third%element%grdpts_id(i,j)
                      end do
                      
                      do i=size(sublayer_third%element%grdpts_id,1)+1, i_max
                         temp_grdptid(i,j_min3+j) = no_pt
                      end do
                   end do

                else

                   do j=j_min_block, j_max_block
                      do i=1, size(sublayer_third%element%grdpts_id,1)
                         temp_grdptid(i,j_min3+j) = sublayer_third%element%grdpts_id(i,j)
                      end do
                      
                      do i=size(sublayer_third%element%grdpts_id,1)+1, i_max
                         temp_grdptid(i,j_min3+j) = no_pt
                      end do
                   end do

                end if
                
             end if

          end if
             
        end subroutine merge_grdpts_id_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine merging two nodes tables into one for the
        !> Westhern sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param temp_grdptid
        !> temporary table where the gridpoint ID are copied for the
        !> merge
        !
        !>@param sublayer_first
        !> pointer to the first sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param sublayer_third
        !> pointer to the third sublayer in the convention required
        !> by the subroutine merging sublayers
        !
        !>@param interior_j_max1
        !> y-size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 3)
        !
        !>@param interior_j_max2
        !> y-size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 3)
        !
        !>@param interior_j_max3
        !> y-size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 3)
        !
        !>@param j_max1
        !> y-size of the block copied from the first sublayer merged
        !>  (see Fig. 3)
        !
        !>@param j_max3
        !> y-size of the block copied from the second sublayer merged
        !>  (see Fig. 3)
        !
        !>@param j_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 3)
        !
        !>@param j_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 3)
        !
        !>@param j_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 3)
        !
        !>@param j_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 3)
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine merge_grdpts_id_W(
     $     temp_grdptid,
     $     sublayer_first,
     $     sublayer_third,
     $     interior_j_max1, interior_j_max2, interior_j_max3,
     $     j_max1, j_max3,
     $     j_min1, j_min2, j_min3, j_min4,
     $     i_max,
     $     neighbor_N,
     $     neighbor_S)

          implicit none

          integer, dimension(:,:)    , intent(out) :: temp_grdptid
          type(bf_sublayer), pointer , intent(in)  :: sublayer_first
          type(bf_sublayer), pointer , intent(in)  :: sublayer_third
          integer(ikind)             , intent(in)  :: interior_j_max1
          integer(ikind)             , intent(in)  :: interior_j_max2
          integer(ikind)             , intent(in)  :: interior_j_max3
          integer(ikind)             , intent(in)  :: j_max1, j_max3
          integer(ikind)             , intent(in)  :: j_min1
          integer(ikind)             , intent(in)  :: j_min2
          integer(ikind)             , intent(in)  :: j_min3
          integer(ikind)             , intent(in)  :: j_min4
          integer(ikind)             , intent(in)  :: i_max
          logical                    , intent(in)  :: neighbor_N
          logical                    , intent(in)  :: neighbor_S


          integer(ikind) :: i,j, i_min, j_min_block, j_max_block
          
          !block1: nodes copied from the interior domain
          !        frontier between block1 and block2
          if(interior_j_max1.ge.2) then
             
             !> grdpt_id corresponding to nodes copied
             !> from the interior domain: 1st block
             !1st block : j=1, bc_size
             if(neighbor_S) then
                do j=1, bc_size
                   do i=1, i_max-(2*bc_size+1)
                      temp_grdptid(i,j) = no_pt
                   end do
                   do i=i_max-(2*bc_size+1)+1, i_max
                      temp_grdptid(i,j) = exchange_pt
                   end do
                end do
             else
                call add_layer1_out_W(
     $               1, temp_grdptid, i_max)

                call add_layer2_out_W(
     $               2, temp_grdptid, i_max)
             end if             

             call add_layer_interior_W(
     $            3, interior_j_max1, temp_grdptid, i_max)

             !> frontier between the 1st and 2nd blocks
             call add_layer1_in_W(
     $            interior_j_max1+1, temp_grdptid,
     $            1, sublayer_first%element%grdpts_id,
     $            i_max)

             call add_layer1_in_W(
     $            interior_j_max1+2, temp_grdptid,
     $            2, sublayer_first%element%grdpts_id,
     $            i_max)

             !> update the lower border for the next block copied
             j_min_block = bc_size+1

          else

             if(interior_j_max1.eq.1) then
                
                !1st block : j=1, bc_size
                if(neighbor_S) then
                   j=1
                   do i=1, i_max-(2*bc_size+1)
                      temp_grdptid(i,j) = no_pt
                   end do
                   do i=i_max-(2*bc_size+1)+1, i_max
                      temp_grdptid(i,j) = exchange_pt
                   end do

                   j=2
                   i_min = i_max-size(sublayer_first%element%grdpts_id,1)

                   do i=1, i_min
                      temp_grdptid(i,j) = no_pt
                   end do

                   do i=i_min+1, i_max-(2*bc_size+1)
                      temp_grdptid(i,j) = 
     $                     sublayer_first%element%grdpts_id(i-i_min,j-1)
                   end do

                   do i=i_max-(2*bc_size+1)+1, i_max
                      temp_grdptid(i,j) = exchange_pt
                   end do
                   
                else

                   !> grdpt_id corresponding to nodes copied
                   !> from the interior domain: 1st block
                   call add_layer1_out_W(
     $                  1, temp_grdptid, i_max)
                   
                   !> frontier between the 1st and 2nd blocks
                   call add_layer2_in_W(
     $                  2, temp_grdptid,
     $                  1, sublayer_first%element%grdpts_id,
     $                  i_max)

                end if
                
                call add_layer1_in_W(
     $               3, temp_grdptid,
     $               2, sublayer_first%element%grdpts_id,
     $               i_max)

                !> update the lower border for the next block copied
                j_min_block = bc_size+1
                
             else

                if(neighbor_S) then
                   
                   i_min = i_max-size(sublayer_first%element%grdpts_id,1)

                   do j=1, bc_size
                      do i=1, i_min
                         temp_grdptid(i,j) = no_pt
                      end do

                      do i=i_min+1, i_max-(2*bc_size+1)
                         temp_grdptid(i,j) = 
     $                        sublayer_first%element%grdpts_id(i-i_min,j)
                      end do
                      
                      do i=i_max-(2*bc_size+1)+1, i_max
                         temp_grdptid(i,j) = exchange_pt
                      end do
                   end do

                   !> update the lower border for the next block copied
                   j_min_block = bc_size+1
                else
                
                   !> update the lower border for the next block copied
                   j_min_block = 1
                end if

             end if
          end if
                

          !block2: grdpt copied from the sublayer_first
          !        frontier between block2 and block3
          !block3: grdpt corresponding to nodes copied
          !        from interior to the buffer layer
          !        frontier between block3 and block4
          if(interior_j_max2.ge.1) then
             
             !> update the upper border for block2 copied
             j_max_block = j_max1-bc_size

             !> copy the block corresponding to sublayer_
             !> first
             i_min = i_max-size(sublayer_first%element%grdpts_id,1)
             do j=j_min_block, j_max_block
                do i=1, i_min
                   temp_grdptid(i,j_min1+j) = no_pt
                end do

                do i=i_min+1, i_max
                   temp_grdptid(i,j_min1+j) = sublayer_first%element%grdpts_id(
     $                  i-i_min,j)
                end do
                
             end do

             !> frontier between block2 and block3
             call add_layer1_in_W(
     $            j_min1+j_max1-1, temp_grdptid,
     $            j_max_block+1, sublayer_first%element%grdpts_id,
     $            i_max)

             call add_layer1_in_W(
     $            j_min1+j_max1, temp_grdptid,
     $            j_max_block+2, sublayer_first%element%grdpts_id,
     $            i_max)

             call add_layer_interior_W(
     $            j_min2+1, j_min2+interior_j_max2, temp_grdptid, i_max)

             !> frontier between the 3rd and 4th blocks
             call add_layer1_in_W(
     $            j_min3+1, temp_grdptid,
     $            1, sublayer_third%element%grdpts_id,
     $            i_max)

             call add_layer1_in_W(
     $            j_min3+bc_size, temp_grdptid,
     $            2, sublayer_third%element%grdpts_id,
     $            i_max)

             !> update the lower border for the next block copied
             j_min_block = bc_size+1

          else

             !> update the upper border for block2 copied
             j_max_block = j_max1

             !> copy the block corresponding to sublayer_
             !> first
             i_min = i_max-size(sublayer_first%element%grdpts_id,1)
             do j=j_min_block, j_max_block
                do i=1,i_min
                   temp_grdptid(i,j_min1+j) = no_pt
                end do

                do i=i_min+1, i_max
                   temp_grdptid(i,j_min1+j) = sublayer_first%element%grdpts_id(
     $                  i-i_min,j)
                end do            
             end do

             !> update the lower border for the next block copied
             j_min_block = 1

          end if

          
          !block4: grdpt copied from the sublayer_third
          !        frontier between block4 and block5
          !block5: grdpt corresponding to nodes copied
          !        from interior to the buffer layer
          !        frontier between block5 and the edge
          if(interior_j_max3.ge.2) then

             !> update the upper border for block2 copied
             j_max_block = j_max3-bc_size

             !> copy the block corresponding to sublayer_
             !> first
             i_min = i_max-size(sublayer_third%element%grdpts_id,1)
             do j=j_min_block, j_max_block
                do i=1, i_min
                   temp_grdptid(i,j_min3+j) = no_pt
                end do
                
                do i=i_min+1,i_max
                   temp_grdptid(i,j_min3+j) = sublayer_third%element%grdpts_id(
     $                  i-i_min,j)
                end do
             end do

             !> frontier between block4 and block5
             call add_layer1_in_W(
     $            j_min3+j_max3-1, temp_grdptid,
     $            j_max_block+1, sublayer_third%element%grdpts_id,
     $            i_max)

             call add_layer1_in_W(
     $            j_min3+j_max3, temp_grdptid,
     $            j_max_block+2, sublayer_third%element%grdpts_id,
     $            i_max)

             call add_layer_interior_W(
     $            j_min4+1, j_min4+interior_j_max3-bc_size, temp_grdptid, i_max)

             !> frontier at the edge
             if(neighbor_N) then
                do j=j_min4+interior_j_max3-1, j_min4+interior_j_max3
                   do i=1, i_max-(2*bc_size+1)
                      temp_grdptid(i,j) = no_pt
                   end do
                   do i=i_max-(2*bc_size+1)+1, i_max
                      temp_grdptid(i,j) = exchange_pt
                   end do
                end do
             else
                call add_layer2_out_W(
     $               j_min4+interior_j_max3-1, temp_grdptid, i_max)
                
                call add_layer1_out_W(
     $               j_min4+interior_j_max3  , temp_grdptid, i_max)
             end if

          else

             if(interior_j_max3.eq.1) then
                
                !> update the upper border for block2 copied
                j_max_block = j_max3-bc_size

                !> copy the block corresponding to sublayer_
                !> first
                i_min = i_max-size(sublayer_third%element%grdpts_id,1)
                do j=j_min_block, j_max_block
                   do i=1, i_min
                      temp_grdptid(i,j_min3+j) = no_pt
                   end do

                   do i=i_min+1, i_max
                      temp_grdptid(i,j_min3+j) = sublayer_third%element%grdpts_id(
     $                     i-i_min,j)
                   end do
                end do

                !> frontier between block4 and the edge
                call add_layer1_in_W(
     $               j_min3+j_max_block+1, temp_grdptid,
     $               j_max_block+1, sublayer_third%element%grdpts_id,
     $               i_max)


                !block5: last two layers
                if(neighbor_N) then

                   j=j_min3+j_max_block+2
                   i_min = i_max-size(sublayer_third%element%grdpts_id,1)

                   do i=1, i_min
                      temp_grdptid(i,j) = no_pt
                   end do
                   
                   do i=i_min+1, i_max-(2*bc_size+1)
                      temp_grdptid(i,j) = 
     $                     sublayer_third%element%grdpts_id(i,j_max_block+2)
                   end do
                   
                   do i=i_max-(2*bc_size+1)+1, i_max
                      temp_grdptid(i,j) = exchange_pt
                   end do

                   j=j_min3+j_max_block+3

                   do i=1, i_max-(2*bc_size+1)
                      temp_grdptid(i,j) = no_pt
                   end do

                   do i=i_max-(2*bc_size+1)+1, i_max
                      temp_grdptid(i,j) = exchange_pt
                   end do

                else

                   call add_layer2_in_W(
     $                  j_min3+j_max_block+2, temp_grdptid,
     $                  j_max_block+2, sublayer_third%element%grdpts_id,
     $                  i_max)

                   call add_layer1_out_W(
     $                  j_min3+j_max_block+3, temp_grdptid, i_max)
                end if
                
             else

                !> update the upper border for block2 copied
                j_max_block = j_max3

                i_min = i_max-size(sublayer_third%element%grdpts_id,1)

                if(neighbor_N) then

                   !> copy the block corresponding to sublayer_
                   !> third                   
                   do j=j_min_block, j_max_block-bc_size
                      do i=1, i_min
                         temp_grdptid(i,j_min3+j) = no_pt
                      end do
                      
                      do i=i_min+1,i_max
                         temp_grdptid(i,j_min3+j) = sublayer_third%element%grdpts_id(
     $                        i-i_min,j)
                      end do
                   end do

                   do j=j_max_block-bc_size+1, j_max_block
                      do i=1, i_min
                         temp_grdptid(i,j_min3+j) = no_pt
                      end do
                      
                      do i=i_min+1,i_max-(2*bc_size+1)
                         temp_grdptid(i,j_min3+j) = sublayer_third%element%grdpts_id(
     $                        i-i_min,j)
                      end do

                      do i=i_max-(2*bc_size+1)+1, i_max
                         temp_grdptid(i,j_min3+j) = exchange_pt
                      end do
                   end do

                else                   
                   
                   !> copy the block corresponding to sublayer_
                   !> third                   
                   do j=j_min_block, j_max_block
                      do i=1, i_min
                         temp_grdptid(i,j_min3+j) = no_pt
                      end do
                      
                      do i=i_min+1,i_max
                         temp_grdptid(i,j_min3+j) = sublayer_third%element%grdpts_id(
     $                        i-i_min,j)
                      end do
                   end do

                end if
                
             end if

          end if
             
        end subroutine merge_grdpts_id_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer1_out_E(
     $     j, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max

          integer(ikind) :: i

          ! exchange
          !    |  bc_pt   no_pt
          !    |   _|_   ___|____
          !   / \ /   \ /        \
          !  ---------------------
          !  |0 0 2 2 2          |
          !  ---------------------

          do i=1, bc_size
             temp_grdptid(i,j) = exchange_pt
          end do

          do i=bc_size+1, 2*bc_size+1
             temp_grdptid(i,j) = bc_pt
          end do

          do i=2*bc_size+2, i_max
             temp_grdptid(i,j) = no_pt
          end do             

        end subroutine add_layer1_out_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer2_out_E(
     $     j, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max

          integer(ikind) :: i

          !
          ! exchange
          !    |  bc_interior_pt
          !    |   |        
          !    |   | bc_pt no_pt
          !    |   |  |  ___|____
          !   / \ / \ | /        \
          !  ---------------------
          !  |0 0 1 1 2          |
          !  ---------------------

          do i=1, bc_size
             temp_grdptid(i,j) = exchange_pt
          end do

          do i=bc_size+1, 2*bc_size
             temp_grdptid(i,j) = bc_interior_pt
          end do

          i=2*bc_size+1
          temp_grdptid(i,j) = bc_pt

          do i=2*bc_size+2, i_max
             temp_grdptid(i,j) = no_pt
          end do             

        end subroutine add_layer2_out_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer:
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j_min
        !> y-index identifying the lower border of the layers
        !> initialized
        !
        !>@param j_max
        !> y-index identifying the upper border of the layers
        !> initialized
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer_interior_E(
     $     j_min, j_max, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j_min
          integer(ikind)         , intent(in)  :: j_max
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max


          integer(ikind) :: i,j


          ! exchange
          !    | interior_pt
          !    |  | bc_interior_pt
          !    |  | | bc_pt
          !    |  | | |   no_pt
          !    |  | | |  ___|____
          !   / \ | | | /        \
          !  ---------------------
          !  |0 0 1 2 3          |
          !  ---------------------

          do j=j_min, j_max
             do i=1, bc_size
                temp_grdptid(i,j) = exchange_pt
             end do

             i=bc_size+1
             temp_grdptid(i,j) = interior_pt

             i=bc_size+2
             temp_grdptid(i,j) = bc_interior_pt

             i=2*bc_size+1
             temp_grdptid(i,j) = bc_pt

             do i=2*bc_size+2, i_max
                temp_grdptid(i,j) = no_pt
             end do
          end do

        end subroutine add_layer_interior_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer:
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param j_match
        !> y-index identifying the layer corresponding in the
        !> sublayer copied for merge
        !
        !>@param grdptid
        !> table for the gridpoint ID of the sublayer copied for
        !> merge
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer1_in_E(
     $     j, temp_grdptid, j_match, grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in) :: j
          integer, dimension(:,:), intent(out):: temp_grdptid
          integer(ikind)         , intent(in) :: j_match
          integer, dimension(:,:), intent(in) :: grdptid
          integer(ikind)         , intent(in) :: i_max

          integer(ikind) :: i


          ! exchange
          !    | interior
          !    |  | bc_interior
          !    |  | | 
          !    |  | |     grdptid   no_pt
          !    |  | |  _____|____  ___|___
          !   / \ | | /          \/       \
          !  ------------------------------
          !  |0 0 1 2                     |
          !  ------------------------------


          do i=1, bc_size
             temp_grdptid(i,j) = exchange_pt
          end do

          i=bc_size+1
          temp_grdptid(i,j) = interior_pt
          
          i=bc_size+2
          temp_grdptid(i,j) = bc_interior_pt
          
          do i=bc_size+3, size(grdptid,1)
             temp_grdptid(i,j) = grdptid(i,j_match)
          end do
          
          do i=size(grdptid,1)+1, i_max
             temp_grdptid(i,j) = no_pt
          end do
          
        end subroutine add_layer1_in_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param j_match
        !> y-index identifying the layer corresponding in the
        !> sublayer copied for merge
        !
        !>@param grdptid
        !> table for the gridpoint ID of the sublayer copied for
        !> merge
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer2_in_E(
     $     j, temp_grdptid, j_match, grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in) :: j
          integer, dimension(:,:), intent(out):: temp_grdptid
          integer(ikind)         , intent(in) :: j_match
          integer, dimension(:,:), intent(in) :: grdptid
          integer(ikind)         , intent(in) :: i_max

          integer(ikind) :: i


          ! exchange
          !    | 
          !    | bc_interior
          !    |   | 
          !    |   |      grdptid   no_pt
          !    |   |   _____|____  ___|___
          !   / \ / \ /          \/       \
          !  ------------------------------
          !  |0 0 1 1                     |
          !  ------------------------------

          do i=1, bc_size
             temp_grdptid(i,j) = exchange_pt
          end do

          do i=bc_size+1, bc_size+2
             temp_grdptid(i,j) = bc_interior_pt
          end do
          
          do i=bc_size+3, size(grdptid,1)
             temp_grdptid(i,j) = grdptid(i,j_match)
          end do
          
          do i=size(grdptid,1)+1, i_max
             temp_grdptid(i,j) = no_pt
          end do
          
        end subroutine add_layer2_in_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Western sublayer
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer1_out_W(
     $     j, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max

          integer(ikind) :: i

          !            
          !                 exchange
          !    no_pt    bc_pt  |  
          !   ___|____   _|_   |  
          !  /        \ /   \ / \ 
          ! ----------------------
          ! |         | 2 2 2 0 0|
          ! ----------------------

          do i=1, i_max-(2*bc_size+1)
             temp_grdptid(i,j) = no_pt
          end do

          do i=i_max-(2*bc_size+1)+1, i_max-bc_size
             temp_grdptid(i,j) = bc_pt
          end do

          do i=i_max-bc_size+1, i_max
             temp_grdptid(i,j) = exchange_pt
          end do

        end subroutine add_layer1_out_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer2_out_W(
     $     j, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max

          integer(ikind) :: i


          !
          !        bc_interior_pt
          !              |
          !  no_pt  bc_pt|  exchange
          !  ___|____ |  |    |        
          ! /        \| / \  / \       
          ! ---------------------      
          ! |        |2 1 1 |0 0|     
          ! ---------------------

          do i=1, i_max-(2*bc_size+1)
             temp_grdptid(i,j) = no_pt
          end do

          i=i_max-(2*bc_size+1)+1
          temp_grdptid(i,j) = bc_pt

          do i=i_max-(2*bc_size+1)+2, i_max-bc_size
             temp_grdptid(i,j) = bc_interior_pt
          end do

          do i=i_max-bc_size+1, i_max
             temp_grdptid(i,j) = exchange_pt
          end do

        end subroutine add_layer2_out_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j_min
        !> y-index identifying the lower border of the layers
        !> initialized
        !
        !>@param j_max
        !> y-index identifying the upper border of the layers
        !> initialized
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer_interior_W(
     $     j_min, j_max, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j_min
          integer(ikind)         , intent(in)  :: j_max
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max


          integer(ikind) :: i,j
          
          !                exchange_pt  
          !         interior_pt|        
          !                 |  |        
          !   bc_interior_pt|  |        
          !    no_pt      | |  |        
          !   ___|__ bc_pt  |  |        
          !  /        \ | | | / \       
          ! -----------------------     
          !           | 3 2 1 0 0 |     
          ! -----------------------     

          do j=j_min, j_max
             do i=1, i_max-bc_size-3
                temp_grdptid(i,j) = no_pt
             end do
             
             i=i_max-bc_size-2
             temp_grdptid(i,j) = bc_pt

             i=i_max-bc_size-1
             temp_grdptid(i,j) = bc_interior_pt

             i=i_max-bc_size
             temp_grdptid(i,j) = interior_pt

             do i=i_max-bc_size+1, i_max
                temp_grdptid(i,j) = exchange_pt
             end do
             
          end do

        end subroutine add_layer_interior_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Western sublayer        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param j_match
        !> y-index identifying the layer corresponding in the
        !> sublayer copied for merge
        !
        !>@param grdptid
        !> table for the gridpoint ID of the sublayer copied for
        !> merge
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer1_in_W(
     $     j, temp_grdptid, j_match, grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in) :: j
          integer, dimension(:,:), intent(out):: temp_grdptid
          integer(ikind)         , intent(in) :: j_match
          integer, dimension(:,:), intent(in) :: grdptid
          integer(ikind)         , intent(in) :: i_max

          integer(ikind) :: i, i_min


          !
          !                      exchange
          !                  interior | 
          !            bc_interior |  | 
          !                      | |  | 
          !   no_pt      grdptid | |  | 
          !  ___|___  _____|____ | |  | 
          ! /       \/          \| | / \
          ! ------------------------------
          ! |        |           |2 1 0 0|
          ! ------------------------------

          i_min = i_max-size(grdptid,1)

          do i=1, i_min
             temp_grdptid(i,j) = no_pt
          end do

          do i=i_min+1, i_max-bc_size-2
             temp_grdptid(i,j) = grdptid(i-i_min,j_match)
          end do

          i=i_max-bc_size-1
          temp_grdptid(i,j) = bc_interior_pt

          i=i_max-bc_size
          temp_grdptid(i,j) = interior_pt

          do i=i_max-bc_size+1, i_max
             temp_grdptid(i,j) = exchange_pt
          end do
          
        end subroutine add_layer1_in_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Western sublayer
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param j_match
        !> y-index identifying the layer corresponding in the
        !> sublayer copied for merge
        !
        !>@param grdptid
        !> table for the gridpoint ID of the sublayer copied for
        !> merge
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer2_in_W(
     $     j, temp_grdptid, j_match, grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in) :: j
          integer, dimension(:,:), intent(out):: temp_grdptid
          integer(ikind)         , intent(in) :: j_match
          integer, dimension(:,:), intent(in) :: grdptid
          integer(ikind)         , intent(in) :: i_max

          integer(ikind) :: i, i_min


          !                       exchange
          !                           |
          !              bc_interior  |       
          !                       |   |  
          !   no_pt      grdptid  |   |  
          !  ___|___  _____|____  |   |  
          ! /       \/          \/ \ / \ 
          ! -----------------------------
          ! |                   |1 1 0 0|
          ! -----------------------------


          i_min = i_max-size(grdptid,1)

          do i=1, i_min
             temp_grdptid(i,j) = no_pt
          end do

          do i=i_min+1, i_max-bc_size-2
             temp_grdptid(i,j) = grdptid(i-i_min,j_match)
          end do          

          do i=i_max-bc_size-1, i_max-bc_size
             temp_grdptid(i,j) = bc_interior_pt
          end do

          do i=i_max-bc_size+1, i_max
             temp_grdptid(i,j) = exchange_pt
          end do          
          
        end subroutine add_layer2_in_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine adding the exchange layer
        !> for Northern and Southern sublayers        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpts_id
        !> table for the merged sublayer
        !
        !>@param grdpts_id1
        !> table copied from the first sublayer
        !
        !>@param grdpts_id2
        !> table copied from the second sublayer
        !
        !>@param j_min_exchange
        !> y-lower border for the layer
        !
        !>@param j_max_exchange
        !> y-upper border for the layer
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max3
        !> size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 1)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param i_max1
        !> x-size of the first sublayer copied
        !
        !>@param j_min1
        !> y-index identifying the correspondance between the merged
        !> sublayer and the first sublayer copied
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !
        !>@param i_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 1)
        !
        !>@param i_max3
        !> x-size of the second sublayer copied
        !
        !>@param j_min3
        !> y-index identifying the correspondance between the merged
        !> sublayer and the second sublayer copied
        !
        !>@param i_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 1)
        !--------------------------------------------------------------
        subroutine add_exchange_pt_layer_NS(
     $     grdpts_id,
     $     grdpts_id1,
     $     grdpts_id2,
     $     j_min_exchange, j_max_exchange,
     $     interior_i_max1,
     $     interior_i_max2,
     $     interior_i_max3,
     $     i_min1, i_max1, j_min1,
     $     i_min2,
     $     i_min3, i_max3, j_min3,
     $     i_min4)

          implicit none

          
          integer, dimension(:,:), intent(out) :: grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_id1
          integer, dimension(:,:), intent(in)  :: grdpts_id2
          integer(ikind)         , intent(in)  :: j_min_exchange, j_max_exchange
          integer(ikind)         , intent(in)  :: interior_i_max1
          integer(ikind)         , intent(in)  :: interior_i_max2
          integer(ikind)         , intent(in)  :: interior_i_max3
          integer(ikind)         , intent(in)  :: i_min1, i_max1, j_min1
          integer(ikind)         , intent(in)  :: i_min2
          integer(ikind)         , intent(in)  :: i_min3, i_max3, j_min3
          integer(ikind)         , intent(in)  :: i_min4

          integer(ikind) :: i,j


          !interior_i_max1  interior_i_max2   interior_i_max3
          !  ___|____          ___|____          ___|____     
          ! /        \ i_max1 /        \ i_max3 /        \    
          ! ---------------------------------------------- j_max_exchange
          ! |exchange|grdptid1|exchange|grdptid2|exchange|
          ! ---------------------------------------------- j_min_exchange
          ! |        |        |        |        |
          ! 1       i_min1  i_min2   i_min3   i_min4

          do j=j_min_exchange, j_max_exchange
             do i=1, interior_i_max1
                grdpts_id(i,j) =
     $               exchange_pt
             end do
             
             do i=1, i_max1
                grdpts_id(i_min1+i,j) =
     $               grdpts_id1(i,j_min1+j)
             end do

             do i=1, interior_i_max2
                grdpts_id(i_min2+i,j) =
     $               exchange_pt
             end do
             
             do i=1, i_max3
                grdpts_id(i_min3+i,j) =
     $               grdpts_id2(i,j_min3+j)
             end do
             
             do i=1, interior_i_max3
                grdpts_id(i_min4+i,j) =
     $               exchange_pt
             end do
          end do

        end subroutine add_exchange_pt_layer_NS


      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine adding the interior layer
        !> for Northern and Southern sublayers
        !  
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpts_id
        !> table for the merged sublayer
        !
        !>@param grdpts_id1
        !> table copied from the first sublayer
        !
        !>@param grdpts_id2
        !> table copied from the second sublayer
        !
        !>@param j
        !> y-index identifying the layer modified in grdpts_id
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max3
        !> size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 1)
        !
        !>@param i_min_first_layer
        !> index indicating the lower border for the first layer
        !> above the exchange layer
        !>  (see Fig. 5)
        !
        !>@param i_max_first_layer
        !> index indicating the upper border for the first layer
        !> above the exchange layer
        !>  (see Fig. 5)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param i_max1
        !> x-size of the first sublayer copied
        !
        !>@param j_min1
        !> y-index identifying the correspondance between the merged
        !> sublayer and the first sublayer copied
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !
        !>@param i_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 1)
        !
        !>@param i_max3
        !> x-size of the second sublayer copied
        !
        !>@param j_min3
        !> y-index identifying the correspondance between the merged
        !> sublayer and the second sublayer copied
        !
        !>@param i_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 1)
        !--------------------------------------------------------------
        subroutine add_interior_pt_layer_NS(
     $     grdpts_id,
     $     grdpts_id1,
     $     grdpts_id2,
     $     j,
     $     interior_i_max1,
     $     interior_i_max2,
     $     interior_i_max3,
     $     i_min_first_layer,
     $     i_max_first_layer,
     $     i_min1, j_min1,
     $     i_min2,
     $     i_min3, i_max3, j_min3,
     $     i_min4,
     $     neighbor_E,
     $     neighbor_W)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_id1
          integer, dimension(:,:), intent(in)  :: grdpts_id2
          integer(ikind)         , intent(in)  :: j
          integer(ikind)         , intent(in)  :: interior_i_max1
          integer(ikind)         , intent(in)  :: interior_i_max2
          integer(ikind)         , intent(in)  :: interior_i_max3
          integer(ikind)         , intent(in)  :: i_min_first_layer
          integer(ikind)         , intent(in)  :: i_max_first_layer
          integer(ikind)         , intent(in)  :: i_min1, j_min1
          integer(ikind)         , intent(in)  :: i_min2
          integer(ikind)         , intent(in)  :: i_min3, i_max3, j_min3
          integer(ikind)         , intent(in)  :: i_min4
          logical                , intent(in)  :: neighbor_E
          logical                , intent(in)  :: neighbor_W


          integer(ikind) :: i


          !interior_i_max1  interior_i_max2   interior_i_max3
          !  ___|____          ___|____          ___|____     
          ! /        \ i_max1 /        \ i_max3 /        \    
          ! ---------------------------------------------- 
          ! |interior|grdptid1|interior|grdptid2|interior| j
          ! ---------------------------------------------- 
          ! |        | |      |        |      | |
          ! 1   i_min1 |    i_min2   i_min3   | i_min4
          !            |                      |
          !    i_max_first_layer        i_min_first_layer
          !
          !-------------------------------------------------
          !                     Fig. 5
          !-------------------------------------------------


          if(neighbor_W) then
             grdpts_id(1,j) = exchange_pt
             grdpts_id(2,j) = exchange_pt
          else
             grdpts_id(1,j) = bc_pt
             grdpts_id(2,j) = bc_interior_pt
          end if

          do i=3, interior_i_max1+2
             grdpts_id(i,j) = interior_pt
          end do

          do i=3, i_max_first_layer
             grdpts_id(i_min1+i,j) =
     $            grdpts_id1(i,j_min1+j)
          end do
          
          do i=-1, interior_i_max2+2
             grdpts_id(i_min2+i,j) =
     $            interior_pt
          end do

          do i=i_min_first_layer,i_max3-2
             grdpts_id(i_min3+i,j) =
     $            grdpts_id2(i,j_min3+j)
          end do

          do i=-1, interior_i_max3-2
             grdpts_id(i_min4+i,j) =
     $            interior_pt
          end do

          if(neighbor_E) then
             grdpts_id(i_min4+interior_i_max3-1,j) = exchange_pt
             grdpts_id(i_min4+interior_i_max3  ,j) = exchange_pt
          else
             grdpts_id(i_min4+interior_i_max3-1,j) = bc_interior_pt
             grdpts_id(i_min4+interior_i_max3  ,j) = bc_pt
          end if

        end subroutine add_interior_pt_layer_NS


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine adding the bc_interior layer
        !> for Northern and Southern sublayers        
        !  
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpts_id
        !> table for the merged sublayer
        !
        !>@param grdpts_id1
        !> table copied from the first sublayer
        !
        !>@param grdpts_id2
        !> table copied from the second sublayer
        !
        !>@param j
        !> y-index identifying the layer modified in grdpts_id
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max3
        !> size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 1)
        !
        !>@param i_min_second_layer
        !> index indicating the lower border for the second layer
        !> above the exchange layer
        !>  (see Fig. 1)
        !
        !>@param i_max_second_layer
        !> index indicating the upper border for the second layer
        !> above the exchange layer
        !>  (see Fig. 1)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param i_max1
        !> x-size of the first sublayer copied
        !
        !>@param j_min1
        !> y-index identifying the correspondance between the merged
        !> sublayer and the first sublayer copied
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !
        !>@param i_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 1)
        !
        !>@param i_max3
        !> x-size of the second sublayer copied
        !
        !>@param j_min3
        !> y-index identifying the correspondance between the merged
        !> sublayer and the second sublayer copied
        !
        !>@param i_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 1)
        !--------------------------------------------------------------
        subroutine add_bc_interior_pt_layer_NS(
     $     grdpts_id,
     $     grdpts_id1,
     $     grdpts_id2,
     $     j,
     $     interior_i_max1,
     $     interior_i_max2,
     $     interior_i_max3,
     $     i_min_second_layer,
     $     i_max_second_layer,
     $     i_min1, j_min1,
     $     i_min2,
     $     i_min3, i_max3, j_min3,
     $     i_min4,
     $     neighbor_E,
     $     neighbor_W)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_id1
          integer, dimension(:,:), intent(in)  :: grdpts_id2
          integer(ikind)         , intent(in)  :: j
          integer(ikind)         , intent(in)  :: interior_i_max1
          integer(ikind)         , intent(in)  :: interior_i_max2
          integer(ikind)         , intent(in)  :: interior_i_max3
          integer(ikind)         , intent(in)  :: i_min_second_layer
          integer(ikind)         , intent(in)  :: i_max_second_layer
          integer(ikind)         , intent(in)  :: i_min1, j_min1
          integer(ikind)         , intent(in)  :: i_min2
          integer(ikind)         , intent(in)  :: i_min3, i_max3, j_min3
          integer(ikind)         , intent(in)  :: i_min4
          logical                , intent(in)  :: neighbor_E
          logical                , intent(in)  :: neighbor_W

          integer(ikind) :: i


          ! interior_i_max1      interior_i_max2      interior_i_max3
          !  _____|_____          _____|_____          _____|_____
          ! /           \ i_max1 /           \ i_max3 /           \    
          ! -------------------------------------------------------
          ! |bc_interior|grdptid1|bc_interior|grdptid2|bc_interior| j
          ! -------------------------------------------------------
          ! |           | |      |           |      | |
          ! 1      i_min1 |    i_min2      i_min3   | i_min4
          !               |                         |
          !      i_max_second_layer         i_min_second_layer
          !
          !---------------------------------------------------------
          !                          Fig. 6
          !---------------------------------------------------------

          !block1 + block2
          if(neighbor_W) then

             grdpts_id(1,j) = exchange_pt
             grdpts_id(2,j) = exchange_pt

             do i=3, interior_i_max1+1
                grdpts_id(i,j) =
     $               bc_interior_pt
             end do

             if(interior_i_max1.ge.1) then
                do i=2, i_max_second_layer
                   grdpts_id(i_min1+i,j)=
     $                  grdpts_id1(i,j_min1+j)
                end do
             else
                do i=3, i_max_second_layer
                   grdpts_id(i_min1+i,j)=
     $                  grdpts_id1(i,j_min1+j)
                end do
             end if
                
          else

             grdpts_id(1,j) = bc_pt
             
             do i=2, interior_i_max1+1
                grdpts_id(i,j) =
     $               bc_interior_pt
             end do

             do i=2, i_max_second_layer
                grdpts_id(i_min1+i,j)=
     $               grdpts_id1(i,j_min1+j)
             end do
          end if

          !block3
          do i=0, interior_i_max2+1
             grdpts_id(i_min2+i,j) =
     $            bc_interior_pt
          end do

          !block4 + block5
          if(neighbor_E) then

             if(interior_i_max3.ge.1) then
                do i=i_min_second_layer,i_max3-1
                   grdpts_id(i_min3+i,j)=
     $                  grdpts_id2(i,j_min3+j)
                end do
             else
                do i=i_min_second_layer,i_max3-2
                   grdpts_id(i_min3+i,j)=
     $                  grdpts_id2(i,j_min3+j)
                end do
             end if

             do i=0, interior_i_max3-2
                grdpts_id(i_min4+i,j) =
     $               bc_interior_pt
             end do

             grdpts_id(i_min4+interior_i_max3-1,j) = exchange_pt
             grdpts_id(i_min4+interior_i_max3  ,j) = exchange_pt

          else
             do i=i_min_second_layer,i_max3-1
                grdpts_id(i_min3+i,j)=
     $               grdpts_id2(i,j_min3+j)
             end do

             do i=0, interior_i_max3-1
                grdpts_id(i_min4+i,j) =
     $               bc_interior_pt
             end do

             grdpts_id(i_min4+interior_i_max3,j) = bc_pt
          end if

        end subroutine add_bc_interior_pt_layer_NS


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine adding the bc_interior layer
        !> for Northern and Southern sublayers        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpts_id
        !> table for the merged sublayer
        !
        !>@param grdpts_id1
        !> table copied from the first sublayer
        !
        !>@param grdpts_id2
        !> table copied from the second sublayer
        !
        !>@param j
        !> y-index identifying the layer modified in grdpts_id
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max3
        !> size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 1)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param i_max1
        !> x-size of the first sublayer copied
        !
        !>@param j_min1
        !> y-index identifying the correspondance between the merged
        !> sublayer and the first sublayer copied
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !
        !>@param i_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 1)
        !
        !>@param i_max3
        !> x-size of the second sublayer copied
        !
        !>@param j_min3
        !> y-index identifying the correspondance between the merged
        !> sublayer and the second sublayer copied
        !
        !>@param i_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 1)
        !--------------------------------------------------------------
        subroutine add_bc_pt_layer_NS(
     $     grdpts_id,
     $     grdpts_id1,
     $     grdpts_id2,
     $     j,
     $     interior_i_max1,
     $     interior_i_max2,
     $     interior_i_max3,
     $     i_min1, i_max1, j_min1,
     $     i_min2,
     $     i_min3, i_max3, j_min3,
     $     i_min4,
     $     neighbor_E,
     $     neighbor_W)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_id1
          integer, dimension(:,:), intent(in)  :: grdpts_id2
          integer(ikind)         , intent(in)  :: j
          integer(ikind)         , intent(in)  :: interior_i_max1
          integer(ikind)         , intent(in)  :: interior_i_max2
          integer(ikind)         , intent(in)  :: interior_i_max3
          integer(ikind)         , intent(in)  :: i_min1, j_min1, i_max1
          integer(ikind)         , intent(in)  :: i_min2
          integer(ikind)         , intent(in)  :: i_min3, i_max3, j_min3
          integer(ikind)         , intent(in)  :: i_min4
          logical                , intent(in)  :: neighbor_E
          logical                , intent(in)  :: neighbor_W

          integer(ikind) :: i


          ! interior_i_max1      interior_i_max2      interior_i_max3
          !  _____|_____          _____|_____          _____|_____
          ! /           \ i_max1 /           \ i_max3 /           \    
          ! -------------------------------------------------------
          ! |   bc_pt   |grdptid1|    bc_pt  |grdptid2|    bc_pt  | j
          ! -------------------------------------------------------
          ! |           |        |           |        |
          ! 1      i_min1      i_min2      i_min3     i_min4

          if(neighbor_W) then
             grdpts_id(1,j) = exchange_pt
             grdpts_id(2,j) = exchange_pt

             do i=3, interior_i_max1
                grdpts_id(i,j) = bc_pt
             end do

             if(interior_i_max1.ge.1) then
                if(interior_i_max1.eq.1) then
                   do i=2, i_max1
                      grdpts_id(i_min1+i,j)=
     $                     grdpts_id1(i,j_min1+j)
                   end do
                else
                   do i=1, i_max1
                      grdpts_id(i_min1+i,j)=
     $                     grdpts_id1(i,j_min1+j)
                   end do
                end if
             else
                do i=3, i_max1
                   grdpts_id(i_min1+i,j)=
     $                  grdpts_id1(i,j_min1+j)
                end do
             end if
          else
             do i=1, interior_i_max1
                grdpts_id(i,j) = bc_pt
             end do

             do i=1, i_max1
                grdpts_id(i_min1+i,j)=
     $               grdpts_id1(i,j_min1+j)
             end do
          end if          

          do i=1, interior_i_max2
             grdpts_id(i_min2+i,j) =
     $            bc_pt
          end do

          do i=1, i_max3
             grdpts_id(i_min3+i,j)=
     $            grdpts_id2(i,j_min3+j)
          end do

          if(neighbor_E) then
             do i=1, interior_i_max3-2
                grdpts_id(i_min4+i,j) = bc_pt
             end do

             grdpts_id(i_min4+interior_i_max3-1,j) = exchange_pt
             grdpts_id(i_min4+interior_i_max3,  j) = exchange_pt
          else
             do i=1, interior_i_max3
                grdpts_id(i_min4+i,j) = bc_pt
             end do
          end if

        end subroutine add_bc_pt_layer_NS


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine adding the bc_interior layer
        !> for Northern and Southern sublayers
        
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpts_id
        !> table for the merged sublayer
        !
        !>@param grdpts_id1
        !> table copied from the first sublayer
        !
        !>@param grdpts_id2
        !> table copied from the second sublayer
        !
        !>@param j
        !> y-index identifying the layer modified in grdpts_id
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged (see Fig. 1)
        !
        !>@param interior_i_max3
        !> size of the third block copied from the inside if the
        !> alignment of the new sublayer is larger than the two
        !> sublayers after merged
        !>  (see Fig. 1)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param i_max1
        !> x-size of the first sublayer copied
        !
        !>@param j_min1
        !> y-index identifying the correspondance between the merged
        !> sublayer and the first sublayer copied
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !
        !>@param i_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 1)
        !
        !>@param i_max3
        !> x-size of the second sublayer copied
        !
        !>@param j_min3
        !> y-index identifying the correspondance between the merged
        !> sublayer and the second sublayer copied
        !
        !>@param i_min4
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block5 (see Fig. 1)
        !--------------------------------------------------------------
        subroutine add_no_pt_layer_NS(
     $     grdpts_id,
     $     grdpts_id1,
     $     grdpts_id2,
     $     j_min_no_pt, j_max_no_pt,
     $     interior_i_max1,
     $     interior_i_max2,
     $     interior_i_max3,
     $     i_min1, i_max1, j_min1,
     $     i_min2,
     $     i_min3, i_max3, j_min3,
     $     i_min4)

          implicit none

          integer, dimension(:,:), intent(out) :: grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_id1
          integer, dimension(:,:), intent(in)  :: grdpts_id2
          integer(ikind)         , intent(in)  :: j_min_no_pt, j_max_no_pt
          integer(ikind)         , intent(in)  :: interior_i_max1
          integer(ikind)         , intent(in)  :: interior_i_max2
          integer(ikind)         , intent(in)  :: interior_i_max3
          integer(ikind)         , intent(in)  :: i_min1, j_min1, i_max1
          integer(ikind)         , intent(in)  :: i_min2
          integer(ikind)         , intent(in)  :: i_min3, i_max3, j_min3
          integer(ikind)         , intent(in)  :: i_min4

          integer(ikind) :: i,j


          !
          ! interior_i_max1      interior_i_max2      interior_i_max3
          !  _____|_____          _____|_____          _____|_____
          ! /           \ i_max1 /           \ i_max3 /           \    
          ! -------------------------------------------------------
          ! |   no_pt   |grdptid1|    no_pt  |grdptid2|    no_pt  | j
          ! -------------------------------------------------------
          ! |           |        |           |        |
          ! 1      i_min1      i_min2      i_min3     i_min4

          do j=j_min_no_pt, j_max_no_pt
             do i=1, interior_i_max1
                grdpts_id(i,j)=
     $               no_pt
             end do
             
             do i=1, i_max1
                grdpts_id(i_min1+i,j)=
     $               grdpts_id1(i,j_min1+j)
             end do
             
             do i=1, interior_i_max2
                grdpts_id(i_min2+i,j) =
     $               no_pt
             end do

             do i=1, i_max3
                grdpts_id(i_min3+i,j)=
     $               grdpts_id2(i,j_min3+j)
             end do
             
             do i=1, interior_i_max3
                grdpts_id(i_min4+i,j) = no_pt
             end do
          end do

        end subroutine add_no_pt_layer_NS



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine investigating the final size of the 
        !> sublayer after merging two North or South sublayers
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param sublayer1
        !> first sublayer merged
        !
        !>@param sublayer2
        !> second sublayer merged
        !
        !>@param alignment
        !> table specifying the position of the new merged layers
        !> compared to the previous one
        !
        !>@param new_size
        !> size of the sublayer after merging
        !--------------------------------------------------------------
        function get_new_size_NS(sublayer1,sublayer2,alignment)
     $     result(new_size)

          implicit none

          type(bf_sublayer), pointer              , intent(in) :: sublayer1
          type(bf_sublayer), pointer              , intent(in) :: sublayer2
          integer(ikind), dimension(2,2), optional, intent(in) :: alignment
          integer(ikind), dimension(2)                         :: new_size

          if(present(alignment)) then

             new_size(1) = 
     $            max(
     $            sublayer1%element%alignment(1,2), 
     $            sublayer2%element%alignment(1,2),
     $            alignment(1,2))
     $            -
     $            min(
     $            sublayer1%element%alignment(1,1), 
     $            sublayer2%element%alignment(1,1),
     $            alignment(1,1))
     $            +
     $            2*bc_size+1

             new_size(2)=
     $            max(
     $            size(sublayer1%element%grdpts_id,2),
     $            size(sublayer2%element%grdpts_id,2))
             
          else

             new_size(1) = 
     $            max(
     $            sublayer1%element%alignment(1,2), 
     $            sublayer2%element%alignment(1,2))
     $            -
     $            min(
     $            sublayer1%element%alignment(1,1), 
     $            sublayer2%element%alignment(1,1))
     $            +
     $            2*bc_size+1

             new_size(2)=
     $            max(
     $            size(sublayer1%element%grdpts_id,2),
     $            size(sublayer2%element%grdpts_id,2))

          end if

        end function get_new_size_NS


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine investigating the final size of the 
        !> sublayer after merging two East or West sublayers
        !
        !> @date
        !> 06_05_2014 - initial version - J.L. Desmarais
        !
        !>@param sublayer1
        !> first sublayer merged
        !
        !>@param sublayer2
        !> second sublayer merged
        !
        !>@param alignment
        !> table specifying the position of the new merged layers
        !> compared to the previous one
        !
        !>@param new_size
        !> size of the sublayer after merging
        !--------------------------------------------------------------
        function get_new_size_EW(sublayer1,sublayer2,alignment)
     $     result(new_size)

          implicit none

          type(bf_sublayer), pointer              , intent(in) :: sublayer1
          type(bf_sublayer), pointer              , intent(in) :: sublayer2
          integer(ikind), dimension(2,2), optional, intent(in) :: alignment
          integer(ikind), dimension(2)                         :: new_size

          if(present(alignment)) then

             new_size(1)=
     $            max(
     $            size(sublayer1%element%grdpts_id,1),
     $            size(sublayer2%element%grdpts_id,1))

             new_size(2) = 
     $            max(
     $            sublayer1%element%alignment(2,2), 
     $            sublayer2%element%alignment(2,2),
     $            alignment(2,2))
     $            -
     $            min(
     $            sublayer1%element%alignment(2,1), 
     $            sublayer2%element%alignment(2,1),
     $            alignment(2,1))
     $            +
     $            2*bc_size+1
             
          else

             new_size(2)=
     $            max(
     $            size(sublayer1%element%grdpts_id,1),
     $            size(sublayer2%element%grdpts_id,1))

             new_size(1) = 
     $            max(
     $            sublayer1%element%alignment(2,2),
     $            sublayer2%element%alignment(2,2))
     $            -
     $            min(
     $            sublayer1%element%alignment(2,1), 
     $            sublayer2%element%alignment(2,1))
     $            +
     $            2*bc_size+1             

          end if

        end function get_new_size_EW

      end module bf_sublayers_merge_module

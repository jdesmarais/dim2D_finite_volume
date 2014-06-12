      !> @file
      !> module implementing the subroutines needed to manage
      !> the allocated memory for the buffer layers when a
      !> current path is analyzed
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the subroutines needed to manage
      !> the allocated memory for the buffer layers when a
      !> current path is analyzed
      !
      !> @date
      ! 23_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_update_allocation_module

        use bf_path_icr_abstract_class, only : bf_path_icr_abstract
        use bf_mainlayer_class          , only : bf_mainlayer
        use bf_sublayer_class           , only : bf_sublayer
        use interface_abstract_class    , only : interface_abstract
        use parameters_constant         , only : N,S,E,W
        use parameters_input            , only : nx,ny,ne,bc_size
        use parameters_kind             , only : rkind, ikind

        implicit none

        private
        public :: update_allocation_bf_layers


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the neighbor of the current sublayer in the clockwise
        !> direction
        !
        !> @date
        !> 23_04_2013 - initial version - J.L. Desmarais
        !--------------------------------------------------------------
        function update_allocation_bf_layers(
     $       interface_used,
     $       nodes,
     $       current_path)
     $       result(modified_sublayer)

          implicit none

          class(interface_abstract)       , intent(inout) :: interface_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          class(bf_path_icr_abstract)   , intent(inout) :: current_path
          type(bf_sublayer), pointer                      :: modified_sublayer

         
          type(bf_mainlayer), pointer :: mainlayer
          type(bf_sublayer) , pointer :: neighbor_sublayer
          logical                     :: sublayer_merge
          

          !< is there a buffer layer matching the current path ?
          !> if yes, we need to check whether there are neighboring
          !> buffer layers that need to be merged
          !> otherwise, we can simply allocate the new buffer layer
          !> using the alignment of current_path
          if(associated(current_path%matching_sublayer)) then


             !< does this matching sublayer has a neighbor in the
             !> clockwise direction ?
             neighbor_sublayer => get_clockwise_neighbor(
     $            current_path%matching_sublayer,
     $            current_path%mainlayer)

             !< if there is a neighbor, then we should investigate
             !> whether the two sublayers will merge
             if(associated(neighbor_sublayer)) then

                sublayer_merge = shall_sublayers_be_merged(
     $               neighbor_sublayer,
     $               current_path%alignment)

                !< if the two sublayers should be merged, then
                !> the alignment of the final buffer layer has
                !> to be computed
                if(sublayer_merge) then

                   mainlayer => interface_used%get_mainlayer(current_path%mainlayer)
                   modified_sublayer => mainlayer%merge_sublayers(
     $                  current_path%matching_sublayer,
     $                  neighbor_sublayer,
     $                  nodes=nodes,
     $                  alignment=current_path%alignment,
     $                  neighbors_i=current_path%neighbors)

                !< otherwise, the current sublayer should be simply
                !> reallocated
                else
                   
                   !< reallocate the buffer layer
                   modified_sublayer => reallocate_current_sublayer(
     $                  current_path,
     $                  nodes)
                   
                end if

             !< otherwise, the current sublayer can simply be
             !> reallocated according to the current path
             else

                !< reallocate the buffer layer
                modified_sublayer => reallocate_current_sublayer(
     $               current_path,
     $               nodes)

             end if             

          !< there is no buffer layer matching the current path,
          !> we should simply allocate the new buffer layer using
          !> the alignment of the current path
          else
             
             !< add a new sublayer to the mainlayer corresponding
             !> to the path
             modified_sublayer => interface_used%add_sublayer(
     $            current_path%mainlayer,
     $            current_path%alignment,
     $            nodes,
     $            current_path%neighbors)

          end if


        end function update_allocation_bf_layers


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the neighbor of the current sublayer in the clockwise
        !> direction
        !
        !> @date
        !> 23_04_2013 - initial version - J.L. Desmarais
        !
        !>@param current_sublayer
        !> pointer to the sublayer investigated
        !
        !>@param mainlayer_id
        !> integer identifying to which main layer this sublayer
        !> belongs
        !
        !>@param neighbor_sublayer
        !> pointer to the neighboring sublayer in clockwise direction
        !> of the current sublayer investigated
        !--------------------------------------------------------------
        function get_clockwise_neighbor(
     $     current_sublayer,
     $     mainlayer_id)
     $     result(neighbor_sublayer)

          implicit none

          type(bf_sublayer), pointer, intent(in) :: current_sublayer
          integer                   , intent(in) :: mainlayer_id
          type(bf_sublayer), pointer             :: neighbor_sublayer


          select case(mainlayer_id)
            case(N,W)
               if(associated(current_sublayer%next)) then
                  neighbor_sublayer => current_sublayer%next
               else
                  nullify(neighbor_sublayer)
               end if

            case(S,E)
               if(associated(current_sublayer%prev)) then
                  neighbor_sublayer => current_sublayer%prev
               else
                  nullify(neighbor_sublayer)
               end if

            case default
               print '(''bf_layer_update_allocation_module'')'
               print '(''get_clockwise_neighbor'')'
               print '(''mainlayer_id: '',I2)', mainlayer_id
               stop 'mainlayer_id not recognized'

          end select

        end function get_clockwise_neighbor


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> reallocate the buffer layer according to the alignment
        !
        !> @date
        !> 23_04_2013 - initial version - J.L. Desmarais
        !
        !>@param current_sublayer
        !> pointer to the sublayer investigated
        !
        !>@param nodes
        !> real(rkind) array encapsulating the data for the interior
        !> domain
        !
        !>@param match_table
        !> integer, dimension(2) identifying the grid points after the
        !> reallocation
        !
        !>@param modified_sublayer
        !> pointer to the modified sublayer
        !--------------------------------------------------------------
        function reallocate_current_sublayer(
     $     current_path,
     $     nodes,
     $     match_table)
     $     result(modified_sublayer)

          implicit none

          class(bf_path_icr_abstract)         , intent(inout) :: current_path
          real(rkind), dimension(nx,ny,ne)      , intent(in)    :: nodes
          integer(ikind), dimension(2), optional, intent(out)   :: match_table
          type(bf_sublayer), pointer                            :: modified_sublayer

          integer, dimension(2,2)    :: border_changes

          !< compute the border changes
          border_changes(1,1) = min(
     $         current_path%alignment(1,1)-
     $         current_path%matching_sublayer%element%alignment(1,1)
     $         ,0)

          border_changes(1,2) = max(
     $         current_path%alignment(1,2)-
     $         current_path%matching_sublayer%element%alignment(1,2)
     $         ,0)

          border_changes(2,1) = min(
     $         current_path%alignment(2,1)-
     $         current_path%matching_sublayer%element%alignment(2,1)
     $         ,0)

          border_changes(2,2) = max(
     $         current_path%alignment(2,2)-
     $         current_path%matching_sublayer%element%alignment(2,2)
     $         ,0)
          
          !< reallocate the buffer layer
          !> matching the current path
          if(present(match_table)) then
             call current_path%matching_sublayer%element%reallocate_bf_layer(
     $            border_changes, nodes, match_table)
          else
             call current_path%matching_sublayer%element%reallocate_bf_layer(
     $            border_changes, nodes)
          end if

          !< associate the new buffer layer with
          !> the modified_sublayer pointer
          modified_sublayer => current_path%matching_sublayer

        end function reallocate_current_sublayer      


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> reallocate the buffer layer according to the alignment
        !
        !> @date
        !> 23_04_2013 - initial version - J.L. Desmarais
        !
        !>@param neighbor_sublayer
        !> pointer to the neighbor_sublayer investigated
        !
        !>@param alignment
        !> integer array identifying the extrema positions of the
        !> current path
        !
        !>@param sublayer_merge
        !> logical identifying whether the sublayers should be merged
        !--------------------------------------------------------------
        function shall_sublayers_be_merged(
     $     neighbor_sublayer,
     $     alignment)
     $     result(sublayer_merge)

          implicit none

          type(bf_sublayer), pointer, intent(in) :: neighbor_sublayer
          integer, dimension(2,2)   , intent(in) :: alignment
          logical                                :: sublayer_merge


          !< check if the current path in addition to the
          !> current_sublayer will lead to a merge with
          !> the neighbor_sublayer
          !> the check is done by comparing the alignments
          select case(neighbor_sublayer%element%localization)
            case(N)
               sublayer_merge = 
     $              (alignment(1,2)+bc_size+1).ge.
     $              (neighbor_sublayer%element%alignment(1,1)-bc_size)
            case(S)
               sublayer_merge = 
     $              (alignment(1,1)-bc_size-1).le.
     $              (neighbor_sublayer%element%alignment(1,2)+bc_size)
            case(E)
               sublayer_merge = 
     $              (alignment(2,1)-bc_size-1).le.
     $              (neighbor_sublayer%element%alignment(2,2)+bc_size)
            case(W)
               sublayer_merge = 
     $              (alignment(1,2)+bc_size+1).ge.
     $              (neighbor_sublayer%element%alignment(1,1)-bc_size)
            case default
               print '(''bf_layer_update_allocation_module'')'
               print '(''shall_sublayers_be_merged'')'
               stop ''
          end select


        end function shall_sublayers_be_merged

      end module bf_layer_update_allocation_module

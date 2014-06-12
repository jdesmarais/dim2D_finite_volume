      !> @file
      !> module implementing the subroutines needed to manage
      !> whether additional buffer layers need to be implemented
      !> when a corner buffer layer is allocated
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the subroutines needed to manage
      !> whether additional buffer layers need to be implemented
      !> when a corner buffer layer is allocated
      !
      !> @date
      ! 18_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_corner_check_module
      
        use bf_path_icr_abstract_class, only : bf_path_icr_abstract
        use bf_sublayer_class           , only : bf_sublayer
        use bf_mainlayer_class          , only : bf_mainlayer
        use interface_abstract_class    , only : interface_abstract
        use parameters_bf_layer         , only : interior_pt, bc_interior_pt, bc_pt,
     $                                           exchange_pt, no_pt
        use parameters_constant         , only : N,S,E,W,N_E,N_W,S_E,S_W,
     $                                           x_direction, y_direction
        use parameters_input            , only : nx,ny,ne,bc_size
        use parameters_kind             , only : ikind, rkind

        implicit none

        private
        public :: check_corner_bf_layer_neighbors

        logical, parameter :: debug=.false.
        integer, parameter :: head=0
        integer, parameter :: tail=1


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine checking whether additional sublayers
        !> are needed when a new corner buffer layer is
        !> allocated
        !
        !> @date
        !> 18_04_2013 - initial version - J.L. Desmarais
        !
        !>@param corner_id
        !> cardinal coordinate identifying the corner investigated
        !
        !>@param nodes
        !> table encapsulating the interior domain
        !
        !>@param current_path
        !> optional argument giving the path analyzed that was ended
        !> by the corner
        !--------------------------------------------------------------
        subroutine check_corner_bf_layer_neighbors(
     $       corner_id,
     $       nodes,
     $       interface_used,
     $       current_path)

          implicit none

          integer                                , intent(in)    :: corner_id
          real(rkind), dimension(nx,ny,ne)       , intent(in)    :: nodes
          class(interface_abstract)              , intent(inout) :: interface_used
          class(bf_path_icr_abstract), optional, intent(inout) :: current_path

          integer :: k

          type(bf_sublayer), pointer :: neighbor_sublayer
          integer                    :: neighbor_id
          integer                    :: neighbor_position

          logical                 :: new_neighbor_fixed
          logical                 :: compatible_neighbor
          integer, dimension(2,2) :: new_alignment
          logical, dimension(4)   :: new_neighbors
          integer, dimension(2,2) :: border_changes
          integer, dimension(2)   :: match_table


          !< loop over the neighbors of the corner
          do k=1,2

             !< initialization
             new_neighbor_fixed  = .false.


             !< get the neighboring buffer layer corresponding
             !> to the corner and the neighboring index
             neighbor_sublayer => get_neighboring_sublayer(
     $            corner_id, k, interface_used,
     $            neighbor_id, neighbor_position)


             !< check if the neighboring buffer layer exists
             !> if it exists, we need to check whether it is 
             !> possible/needed to reallocate it to have a
             !> neighboring buffer layer
             if(associated(neighbor_sublayer)) then

                !< we check whether the alignment of the
                !> neighboring buffer layer is compatible
                !> with reallocating/modifying the exchanged
                !> grid points of the neighboring buffer layer
                compatible_neighbor = is_neighbor_alignment_compatible(
     $               neighbor_id, neighbor_position,
     $               new_alignment, new_neighbors,
     $               alignment_i=neighbor_sublayer%get_alignment_tab())

                !< if the neighbor is compatible, the
                !> neighboring buffer layer is reallocated
                !> using the new_alignment and new_neighbors
                !> arguments
                if(compatible_neighbor) then

c$$$                   border_changes(1,1) = min(
c$$$     $                  new_alignment(1,1)-
c$$$     $                  neighbor_sublayer%get_alignment(1,1)
c$$$     $                  ,0)
c$$$                   border_changes(1,2) = max(
c$$$     $                  new_alignment(1,2)-
c$$$     $                  neighbor_sublayer%get_alignment(1,2)
c$$$     $                  ,0)
c$$$                   border_changes(2,1) = min(
c$$$     $                  new_alignment(2,1)-
c$$$     $                  neighbor_sublayer%get_alignment(2,1)
c$$$     $                  ,0)
c$$$                   border_changes(2,2) = max(
c$$$     $                  new_alignment(2,2)-
c$$$     $                  neighbor_sublayer%get_alignment(2,2)
c$$$     $                  ,0)                
                   
                   !< reallocate the existing buffer layer to match the
                   !> borders
                   call neighbor_sublayer%reallocate_bf_layer(new_alignment)

                   !
                   print '(''you need to implement a new copy for'')'
                   print '(''the data from the interior domain   '')'

                   !< update newly allocated grdpts_id
                   call update_allocated_grdpts_id(
     $                  neighbor_id,
     $                  neighbor_sublayer%element%grdpts_id,
     $                  new_neighbors,
     $                  border_changes)

                   print '(''bf_corner_check_module'')'
                   print '(''check_corner_bf_layer_neighbors'')'
                   print '(''you need to implement after the'')'
                   print '(''reallocation to have really new'')'
                   print '(''gridpoints computed'')'

                   new_neighbor_fixed = .true.

                end if

             !< if the neighboring buffer layer does not exist,
             !> we need to know the alignment needed to allocate
             !> a new minimal buffer layer
             else
                compatible_neighbor = is_neighbor_alignment_compatible(
     $               neighbor_id, neighbor_position,
     $               new_alignment, new_neighbors)

             end if


             !< if no candidate as neighboring buffer layer was found,
             !> we need to check whether the current path could lead
             !> to a new buffer layer that could be the neighboring
             !> buffer layer
             if(.not.new_neighbor_fixed) then
                
                !< if there are no neighboring buffer layer
                !> or if the neighbor is not compatible, we check
                !> whether the current path will lead to a possible
                !> neighboring buffer layer, otherwise, a new buffer
                !> layer used as neighbor is allocated
                if(present(current_path)) then
                
                   if(current_path%mainlayer.eq.neighbor_id) then
                      
                      !< check whether the current path is compatible
                      !> as a neighboring buffer layer
                      compatible_neighbor = is_neighbor_alignment_compatible(
     $                     neighbor_id, neighbor_position,
     $                     new_alignment, new_neighbors,
     $                     alignment_i=current_path%alignment,
     $                     neighbors_i=current_path%neighbors)
                   else
                      compatible_neighbor = .false.
                   end if
                   
                   !< if the current path is compatible, the current_path
                   !> data are updated
                   if(compatible_neighbor) then
                      current_path%alignment  = new_alignment
                      current_path%neighbors  = new_neighbors
                      new_neighbor_fixed=.true.

                   end if

                end if
             end if


             !< if no existing layer neither the current path could
             !> lead to a compatible neighbor, a new minimal buffer
             !> layer is allocated 
             if(.not.new_neighbor_fixed) then
                
                 neighbor_sublayer => interface_used%add_sublayer(
     $               neighbor_id,
     $               new_alignment,
     $               nodes,
     $               new_neighbors)

             end if

          end do

        end subroutine check_corner_bf_layer_neighbors
        


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine giving the sublayer corresponding to
        !> the neighboring sublayer of the corner. If there
        !> is no sublayer, the pointer is nullified
        !
        !> @date
        !> 18_04_2013 - initial version - J.L. Desmarais
        !
        !>@param corner_id
        !> cardinal coordinate identifying the corner
        !> investigated
        !
        !>@param neighbor_order
        !> integer identifying the neighbor (left or right)
        !
        !>@param neighbor_id
        !> integer identifying the cardinal coordinate of the
        !> neighboring buffer layer
        !--------------------------------------------------------------
        function get_neighboring_sublayer(
     $     corner_id, neighbor_order, interface_used,
     $     neighbor_id, neighbor_position)
     $     result(neighbor_sublayer)

          implicit none

          integer                  , intent(in)  :: corner_id
          integer                  , intent(in)  :: neighbor_order
          class(interface_abstract), intent(in) :: interface_used
          integer                  , intent(out) :: neighbor_id
          integer                  , intent(out) :: neighbor_position
          type(bf_sublayer), pointer             :: neighbor_sublayer

          integer, dimension(2)       :: neighbor_id_table
          integer, dimension(2)       :: neighbor_position_table
          type(bf_mainlayer), pointer :: neighbor_mainlayer
          
          !> determine the corner
          select case(corner_id)

            case(N_W)
               neighbor_id_table       = [W,N]
               neighbor_position_table = [tail,head]

            case(N_E)
               neighbor_id_table       = [N,E]
               neighbor_position_table = [tail,tail]

            case(S_E)
               neighbor_id_table       = [E,S]
               neighbor_position_table = [head,tail]

            case(S_W)
               neighbor_id_table       = [S,W]
               neighbor_position_table = [head,head]

            case default
               print '(''bf_layer_corner_check_module'')'
               print '(''get_neighboring_sublayer'')'
               print '(''corner_id: '', I2)', corner_id
               stop 'corner_id not recognized'

          end select


          !< check whether the neighbor_order input
          !> is correct
          if(debug) then
             if((neighbor_order.ne.1).and.(neighbor_order.ne.2)) then
               print '(''bf_layer_corner_check_module'')'
               print '(''get_neighboring_sublayer'')'
               stop 'wrong neighbor_order'
             end if
          end if


          !< get the neighbor_id and neighbor_position
          neighbor_id       = neighbor_id_table(neighbor_order)
          neighbor_position = neighbor_position_table(neighbor_order)


          !< try to find the neighbor sublayer
          if(associated(interface_used%mainlayer_pointers(neighbor_id)%ptr)) then

             neighbor_mainlayer => interface_used%mainlayer_pointers(neighbor_id)%ptr


             if(neighbor_position.eq.head) then
                if(associated(neighbor_mainlayer%head_sublayer)) then
                   neighbor_sublayer => neighbor_mainlayer%head_sublayer
                else
                   nullify(neighbor_sublayer)
                end if

             else
                if(associated(neighbor_mainlayer%tail_sublayer)) then
                   neighbor_sublayer => neighbor_mainlayer%tail_sublayer
                else
                   nullify(neighbor_sublayer)
                end if
             end if

          else
             nullify(neighbor_sublayer)
          end if

        end function get_neighboring_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function checking whether the alignment is compatible
        !> with the corner investigated and so whether the 
        !> buffer layer corresponding to this alignment could be
        !> turned into a neighboring buffer layer
        !
        !> @date
        !> 18_04_2013 - initial version - J.L. Desmarais
        !
        !>@param neighbor_id
        !> cardinal coordinate identifying the neighbor
        !> investigated
        !
        !>@param neighbor_position
        !> integer identifying the neighbor position in the
        !> mainlayer
        !
        !>@param new_alignment
        !> integer table identifying the position of the buffer
        !> layer (either to be reallocated or allocated) compared
        !> to the interior domain
        !
        !>@param new_neighbor
        !> logical table identifying the neighbor points for the
        !> buffer layer (either to be reallocated or allocated)
        !
        !>@param alignment_i
        !> optinal integer table identifying the position of the
        !> buffer layer investigated compared to the interior domain
        !--------------------------------------------------------------
        function is_neighbor_alignment_compatible(
     $               neighbor_id, neighbor_position,
     $               new_alignment, new_neighbors,
     $               alignment_i, neighbors_i)
     $     result(compatible)

          implicit none

          integer                          , intent(in)  :: neighbor_id
          integer                          , intent(in)  :: neighbor_position
          integer, dimension(2,2)          , intent(out) :: new_alignment
          logical, dimension(4)            , intent(out) :: new_neighbors
          integer, dimension(2,2), optional, intent(in)  :: alignment_i
          logical, dimension(4)  , optional, intent(in)  :: neighbors_i
          logical                                        :: compatible


          integer :: direction
          integer :: n_direction


          !< we investigate whether the grid point needed by the new
          !> corner buffer layer can be included in an existing
          !> buffer layer or not
          !> if the buffer layer investigated is N,S, the condition
          !> relies on i_min or i_max
          !> if the buffer layer investigated is E,W, the condition
          !> relies on j_min or j_max
          select case(neighbor_id)

            case(N,S)
               direction   = x_direction
               n_direction = nx

            case(E,W)
               direction   = y_direction
               n_direction = ny

            case default
               print '(''bf_layer_corner_check_module'')'
               print '(''is_neighbor_alignment_compatible'')'
               print '(''neighbor_id: '',I2)', neighbor_id
               stop '(''wrong neighbor_id'')'

          end select
          

          !if the neighboring buffer layer position
          !is head, then we check the min value
          !for i or j depending on the neighbor id
          if(neighbor_position.eq.head) then

             !checking if compatible
             compatible=.false.
             if(present(alignment_i)) then

                compatible=abs(alignment_i(direction,1)-(bc_size+1)).le.(2*bc_size+1)

             end if

             !modification of existing buffer layer
             !or path
             if(compatible) then
                new_alignment = alignment_i
                new_alignment(direction,1) = bc_size+1

                if(present(neighbors_i)) then
                   new_neighbors = neighbors_i
                else
                   new_neighbors = [.false.,.false.,.false.,.false.]
                end if
             
             !initialization for a new buffer layer
             else
                new_neighbors = [.false.,.false.,.false.,.false.]
                new_alignment(direction,1) = bc_size+1
                new_alignment(direction,2) = bc_size+1
             end if

             select case(neighbor_id)
               case(N,S)
                  new_neighbors(W)=.true.
               case(E,W)
                  new_neighbors(S)=.true.
             end select

          !if the neighboring buffer layer position
          !is tail, then we check the max value
          !for i or j depending on the neighbor id
          else

             !checking if compatible
             compatible=.false.
             if(present(alignment_i)) then

                compatible=abs(alignment_i(direction,2)-(n_direction-bc_size)).le.(2*bc_size+1)
                
             end if

             !modification of existing buffer layer
             !or path
             if(compatible) then
                new_alignment = alignment_i
                new_alignment(direction,2) = n_direction-bc_size

                if(present(neighbors_i)) then
                   new_neighbors = neighbors_i
                else
                   new_neighbors = [.false.,.false.,.false.,.false.]
                end if
             
             !initialization for a new buffer layer
             else
                new_neighbors = [.false.,.false.,.false.,.false.]
                new_alignment(direction,1) = n_direction-bc_size
                new_alignment(direction,2) = n_direction-bc_size
             end if

             select case(neighbor_id)
               case(N,S)
                  new_neighbors(E)=.true.
               case(E,W)
                  new_neighbors(N)=.true.
             end select

          end if

        end function is_neighbor_alignment_compatible


        subroutine update_allocated_grdpts_id(
     $     neighbor_id,
     $     grdpts_id,
     $     new_neighbors,
     $     borders_changes)

          implicit none

          integer                , intent(in)    :: neighbor_id
          integer, dimension(:,:), intent(inout) :: grdpts_id
          logical, dimension(4)  , intent(in)    :: new_neighbors
          integer, dimension(2,2), intent(in)    :: borders_changes

          integer(ikind) :: i,j, i_min,i_max, j_min, j_max


          if(new_neighbors(N)) then

             select case(neighbor_id)
               case(E)

                  j_min = size(grdpts_id,2)-bc_size-borders_changes(2,2)+1
                  j_max = size(grdpts_id,2)-bc_size

                  !changes due to connection with corner
                  j=j_min
                  grdpts_id(bc_size+1,j) = interior_pt
                  grdpts_id(bc_size+2,j) = bc_interior_pt

                  j=j_min+1
                  grdpts_id(bc_size+1,j) = interior_pt
                  grdpts_id(bc_size+2,j) = bc_interior_pt


                  !changes due to reallocation
                  do j=j_min+2, j_max

                     do i=1, bc_size
                        grdpts_id(i,j) = exchange_pt
                     end do

                     grdpts_id(bc_size+1,j) = interior_pt
                     grdpts_id(bc_size+2,j) = bc_interior_pt
                     grdpts_id(bc_size+3,j) = bc_pt
                     
                  end do


                  !N exchange with corner N_E
                  do j=j_max+1,j_max+bc_size

                     do i=bc_size+1,2*bc_size+1
                        grdpts_id(i,j)=exchange_pt
                     end do
                     
                  end do

               case(W)
                  
                  j_min = size(grdpts_id,2)-bc_size-borders_changes(2,2)+1
                  j_max = size(grdpts_id,2)-bc_size


                  !changes due to connection with corner
                  j=j_min
                  grdpts_id(size(grdpts_id,1)-bc_size-1,j) = bc_interior_pt
                  grdpts_id(size(grdpts_id,1)-bc_size  ,j) = interior_pt


                  j=j_min+1
                  grdpts_id(size(grdpts_id,1)-bc_size-1,j) = bc_interior_pt
                  grdpts_id(size(grdpts_id,1)-bc_size  ,j) = interior_pt


                  !changes due to reallocation
                  do j=j_min+2, j_max
                     
                     grdpts_id(size(grdpts_id,1)-bc_size-2,j) = bc_pt
                     grdpts_id(size(grdpts_id,1)-bc_size-1,j) = bc_interior_pt
                     grdpts_id(size(grdpts_id,1)-bc_size  ,j) = interior_pt

                     do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                        grdpts_id(i,j) = exchange_pt
                     end do

                  end do


                  !N exchange with corner N_W
                  do j=j_max+1,j_max+bc_size

                     do i=size(grdpts_id,1)-2*bc_size,size(grdpts_id,1)-bc_size
                        grdpts_id(i,j)=exchange_pt
                     end do

                  end do

               case default
                  print '(''bf_layer_corner_check_module'')'
                  print '(''update_allocated_grdpts_id'')'
                  stop 'choice not recognized'

             end select

          end if


          if(new_neighbors(S)) then

             select case(neighbor_id)
               case(E)

                  j_min = bc_size+1
                  j_max = bc_size-borders_changes(2,1)-2

                  !S exchange with corner S_E
                  do j=1,bc_size
                     do i=bc_size+1,2*bc_size+1
                        grdpts_id(i,j)=exchange_pt
                     end do
                  end do

                  !changes due to reallocation
                  do j=j_min, j_max
                     do i=1, bc_size
                        grdpts_id(i,j) = exchange_pt
                     end do

                     grdpts_id(bc_size+1,j) = interior_pt
                     grdpts_id(bc_size+2,j) = bc_interior_pt
                     grdpts_id(bc_size+3,j) = bc_pt
                     
                     do i=bc_size+4, size(grdpts_id,1)
                        grdpts_id(i,j) = no_pt
                     end do

                  end do

                  !changes due to connection with corner
                  if(borders_changes(2,1).le.-2) then
                     j=j_max+1
                     grdpts_id(bc_size+1,j) = interior_pt
                     grdpts_id(bc_size+2,j) = bc_interior_pt
                     
                     j=j_max+2
                     grdpts_id(bc_size+1,j) = interior_pt
                     grdpts_id(bc_size+2,j) = bc_interior_pt
                  end if

               case(W)
                  
                  j_min = bc_size+1
                  j_max = bc_size-borders_changes(2,1)-2

                  !S exchange with corner S_W
                  do j=1,bc_size

                     do i=size(grdpts_id,1)-2*bc_size,size(grdpts_id,1)-bc_size
                        grdpts_id(i,j)=exchange_pt
                     end do

                  end do

                  !changes due to reallocation
                  do j=j_min, j_max
                     
                     do i=1, size(grdpts_id,1)-bc_size-3
                        grdpts_id(i,j) = no_pt
                     end do

                     grdpts_id(size(grdpts_id,1)-bc_size-2,j) = bc_pt
                     grdpts_id(size(grdpts_id,1)-bc_size-1,j) = bc_interior_pt
                     grdpts_id(size(grdpts_id,1)-bc_size  ,j) = interior_pt

                     do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                        grdpts_id(i,j) = exchange_pt
                     end do

                  end do

                  !changes due to connection with corner
                  if(borders_changes(2,1).le.-2) then
                     j=j_max+1
                     grdpts_id(size(grdpts_id,1)-bc_size-1,j) = bc_interior_pt
                     grdpts_id(size(grdpts_id,1)-bc_size  ,j) = interior_pt


                     j=j_max+2
                     grdpts_id(size(grdpts_id,1)-bc_size-1,j) = bc_interior_pt
                     grdpts_id(size(grdpts_id,1)-bc_size  ,j) = interior_pt
                  end if


               case default
                  print '(''bf_layer_corner_check_module'')'
                  print '(''update_allocated_grdpts_id'')'
                  stop 'choice not recognized'

             end select

          end if


          if(new_neighbors(E)) then

             select case(neighbor_id)
               case(N)

                  i_min = size(grdpts_id,1)-bc_size-borders_changes(1,2)+2
                  i_max = size(grdpts_id,1)-bc_size
                  
                  !first line above exchange with interior
                  j = bc_size+1
                  do i=i_min-2, i_max
                     grdpts_id(i,j) = interior_pt
                  end do

                  do i=i_max+1,size(grdpts_id,1)
                     grdpts_id(i,j) = exchange_pt
                  end do
                  
                  
                  !second line above exchange with interior
                  j = bc_size+2
                  do i=i_min-1, i_max
                     grdpts_id(i,j) = bc_interior_pt
                  end do

                  do i=i_max+1,size(grdpts_id,1)
                     grdpts_id(i,j) = exchange_pt
                  end do

                  
                  !third line above exchange with interior
                  j = bc_size+3
                  do i=i_min, i_max
                     grdpts_id(i,j) = bc_pt
                  end do

                  do i=i_max+1,size(grdpts_id,1)
                     grdpts_id(i,j) = exchange_pt
                  end do
                  
               case(S)

                  i_min = size(grdpts_id,1)-bc_size-borders_changes(1,2)+2
                  i_max = size(grdpts_id,1)-bc_size

                  !third line above exchange with interior
                  j = size(grdpts_id,2)-(bc_size+2)
                  do i=i_min, i_max
                     grdpts_id(i,j) = bc_pt
                  end do

                  do i=i_max+1,size(grdpts_id,1)
                     grdpts_id(i,j) = exchange_pt
                  end do

                  !second line above exchange with interior
                  j = size(grdpts_id,2)-(bc_size+1)
                  do i=i_min-1, i_max
                     grdpts_id(i,j) = bc_interior_pt
                  end do

                  do i=i_max+1,size(grdpts_id,1)
                     grdpts_id(i,j) = exchange_pt
                  end do

                  !first line above exchange with interior
                  j = size(grdpts_id,2)-(bc_size)
                  do i=i_min-2, i_max
                     grdpts_id(i,j) = interior_pt
                  end do

                  do i=i_max+1,size(grdpts_id,1)
                     grdpts_id(i,j) = exchange_pt
                  end do

               case default
                  print '(''bf_layer_corner_check_module'')'
                  print '(''update_allocated_grdpts_id'')'
                  stop 'choice not recognized'

             end select

          end if


          if(new_neighbors(W)) then

             select case(neighbor_id)
               case(N)

                  i_min = bc_size+1
                  i_max = bc_size-borders_changes(1,1)

                  !first line above exchange with interior
                  j = bc_size+1
                  do i=1,i_min-1
                     grdpts_id(i,j) = exchange_pt
                  end do

                  do i=i_min, i_max+1
                     grdpts_id(i,j) = interior_pt
                  end do
                  
                  !second line above exchange with interior
                  j=bc_size+2
                  do i=1,i_min-1
                     grdpts_id(i,j) = exchange_pt
                  end do

                  do i=i_min, i_max
                     grdpts_id(i,j) = bc_interior_pt
                  end do
                  
                  !third line above exchange with interior
                  j = bc_size+3
                  do i=1,i_min-1
                     grdpts_id(i,j) = exchange_pt
                  end do

                  do i=i_min, i_max-1
                     grdpts_id(i,j) = bc_pt
                  end do
                  

               case(S)

                  i_min = bc_size+1
                  i_max = bc_size-borders_changes(1,1)

                  !third line above exchange with interior
                  j = size(grdpts_id,2)-(bc_size+2)
                  do i=1,i_min-1
                     grdpts_id(i,j) = exchange_pt
                  end do

                  do i=i_min, i_max-2
                     grdpts_id(i,j) = bc_pt
                  end do

                  !second line above exchange with interior
                  j = size(grdpts_id,2)-(bc_size+1)
                  do i=1,i_min-1
                     grdpts_id(i,j) = exchange_pt
                  end do

                  do i=i_min, i_max-1
                     grdpts_id(i,j) = bc_interior_pt
                  end do

                  !first line above exchange with interior
                  j = size(grdpts_id,2)-bc_size
                  do i=1,i_min-1
                     grdpts_id(i,j) = exchange_pt
                  end do

                  do i=i_min, i_max
                     grdpts_id(i,j) = interior_pt
                  end do

               case default
                  print '(''bf_layer_corner_check_module'')'
                  print '(''update_allocated_grdpts_id'')'
                  stop 'choice not recognized'

             end select

          end if

        end subroutine update_allocated_grdpts_id


      end module bf_layer_corner_check_module

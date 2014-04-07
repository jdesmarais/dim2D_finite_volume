      !> @file
      !> module encapsulating the subroutines related to allocation
      !> and reallocation of the buffer layer object
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> subroutines related to allocation and reallocation of the
      !> buffer layer object
      !
      !> @date
      ! 04_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_allocate_module

        use bf_layer_abstract_class      , only : bf_layer_abstract
        use bf_layer_exchange_module     , only : first_exchange_with_interior
        use bf_layer_compute_module      , only : compute_new_grdpts
        use bf_layer_ini_grdptID_module  , only : ini_grdptID

        use parameters_bf_layer, only : no_pt, interior_pt, bc_pt
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input   , only : ne, bc_size
        use parameters_kind    , only : ikind, rkind

        implicit none

        private
        public :: allocate_bf_layer,
     $            reallocate_bf_layer


        logical, parameter :: debug = .false.

        contains

                
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine allocating the main tables of the
        !> buffer layer for the first time and initializing
        !> these tables (nodes, grdptid) using the internal
        !> data and the boundary conditions applied at the
        !> edges
        !
        !> @date
        !> 07_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !
        !>@param alignment
        !> table of integers characterizing the
        !> correspondance between the interior grid points
        !> and the buffer layer elements
        !
        !>@param nodes
        !> table encapsulating the data of the internal
        !> grid points
        !
        !>@param neighbors
        !> table encapsulating whether the buffer layer allocated
        !> has neighboring buffer layers with which it should
        !> exchange data grid points
        !--------------------------------------------------------------
        subroutine allocate_bf_layer(this, alignment, nodes, neighbors)
        
          implicit none

          class(bf_layer_abstract)        , intent(inout) :: this
          integer(ikind), dimension(2,2)  , intent(in)    :: alignment
          real(rkind)   , dimension(:,:,:), intent(in)    :: nodes
          logical       , dimension(4)    , intent(in)    :: neighbors

          integer(ikind), dimension(:,:), allocatable :: list_new_grdpts


          !< set the alignment between the interior table
          !> and the buffer layer allowing to correspond
          !> the interior grid points and the ones from
          !> the buffer layer
          this%alignment = alignment


          !< determine the total size needed for the buffer
          !> layer and allocate the main tables
          call allocate_main_tables(this)


          !> allocate the number of gridpoints needed for
          !< the buffer layer, copy gridpoints from the
          !> interior table, and create the list of new
          !> gridpoints that should be computed
          call first_exchange_with_interior(
     $         this,
     $         nodes,
     $         list_new_grdpts)


          !> identify the grid points of the buffer layer
          call ini_grdptID(this, neighbors)


          !> compute new gridpoints
          call compute_new_grdpts(
     $         this,
     $         list_new_grdpts)

        end subroutine allocate_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine reallocating the main tables of the
        !> buffer layer and copying the previous data in
        !> the new table
        !
        !> @date
        !> 07_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !
        !>@param border_change
        !> table giving the border changes of the previous
        !> allocated table
        !>     - border_change(1,1) : change of i_min
        !>     - border_change(1,2) : change of i_max
        !>     - border_change(2,1) : change of j_min
        !>     - border_change(2,2) : change of j_max
        !>
        !> e.g. : i_min = +1 : the first column of the table is deleted
        !>        i_max = +1 : a new column is added to the table
        !>        j_min = -1 : a new line is added to the table
        !>        j_max = -1 : the last line of the table is deleted
        !
        !>@param match_table
        !> indices allowing to identify correctly the
        !> previously allocated nodes and the new nodes
        !>     - match_table(1) : i_match
        !>     - match_table(2) : j_match
        !>
        !> e.g. : new_table(i,j) = prev_table(i_match+i,j_match+j)
        !--------------------------------------------------------------
        subroutine reallocate_bf_layer(
     $     this,
     $     border_changes,
     $     match_table)
        
          implicit none

          class(bf_layer_abstract), intent(inout) :: this
          integer, dimension(2,2) , intent(in)    :: border_changes
          integer, dimension(2)   , intent(out)   :: match_table


          integer(ikind) :: i,j,k
          integer(ikind) :: i_min, i_max, j_min, j_max
          integer(ikind) :: new_size_x, new_size_y
          real(rkind), dimension(:,:,:), allocatable :: new_nodes
          integer    , dimension(:,:)  , allocatable :: new_grdptid


          !< determine the new alignment between the interior and 
          !> buffer layer tables
          select case(this%localization(1))
            case(N,S,E,W)
               do j=1,2
                  do i=1,2
                     this%alignment(i,j) = this%alignment(i,j) + border_changes(i,j)
                  end do
               end do
          end select


          !< determine the borders when filling the new table with the
          !> old data
          i_min = 1 + max(0,border_changes(1,1)) - border_changes(1,1)
          i_max = size(this%nodes,1) + min(0,border_changes(1,2)) - border_changes(1,1)
          j_min = 1 + max(0,border_changes(2,1)) - border_changes(2,1)
          j_max = size(this%nodes,2) + min(0,border_changes(2,2)) - border_changes(2,1)

          match_table(1) = border_changes(1,1)
          match_table(2) = border_changes(2,1)


          !< determine the new size of the nodes and grdptid tables
          new_size_x = size(this%nodes,1) - border_changes(1,1) + border_changes(1,2)
          new_size_y = size(this%nodes,2) - border_changes(2,1) + border_changes(2,2)
          

          !< allocate the new tables
          allocate(new_nodes(new_size_x,new_size_y,ne))
          allocate(new_grdptid(new_size_x,new_size_y))


          !debugging step to check the inputs
          if(debug) then
             
             select case(this%localization(1))
               case(N)
                  if(border_changes(2,1).ne.0) then
                     stop 'N: border_changes(2,1).ne.0: this is wrong'
                  end if

               case(S)
                  if(border_changes(2,2).ne.0) then
                     stop 'S: border_changes(2,2).ne.0: this is wrong'
                  end if
                  
               case(E)
                  if(border_changes(1,1).ne.0) then
                     stop 'E: border_changes(1,1).ne.0: this is wrong'
                  end if
                  
               case(W)
                  if(border_changes(1,2).ne.0) then
                     stop 'W: border_changes(1,2).ne.0: this is wrong'
                  end if
                  
               case(N_E)
                  if((border_changes(1,1).ne.0).and.(border_changes(2,1).ne.0)) then
                     stop 'NE: change.ne.0: this is wrong'
                  end if
                  
               case(N_W)
                  if((border_changes(1,2).ne.0).and.(border_changes(2,1).ne.0)) then
                     stop 'NW: change.ne.0: this is wrong'
                  end if
                  
               case(S_E)
                  if((border_changes(1,2).ne.0).and.(border_changes(2,2).ne.0)) then
                     stop 'SE: change.ne.0: this is wrong'
                  end if

               case(S_W)
                  if((border_changes(1,1).ne.0).and.(border_changes(2,2).ne.0)) then
                     stop 'SW: change.ne.0: this is wrong'
                  end if
                  
               case default
                  print '(''bf_layer_abstract_class'')'
                  print '(''reallocate_bf_layer'')'
                  stop 'localization not recognized'
             end select
          end if


          !fill the new nodes table with the previous data
          !and transfer the new table in the buffer layer
          do k=1, ne
             do j=j_min, j_max
                do i=i_min, i_max
                   new_nodes(i,j,k) = this%nodes(
     $                  match_table(1)+i,
     $                  match_table(2)+j,
     $                  k)
                end do
             end do
          end do
          call MOVE_ALLOC(new_nodes,this%nodes)


          !fill the new grdptid with the previous data
          !and transfer the new table in the buffer layer
          do j=1, j_min-1
             do i=1, size(this%nodes,1)
                new_grdptid(i,j) = no_pt
                if(debug) this%nodes(i,j,1) = 1.0
             end do
          end do
          
          do j=j_min, j_max

             do i=1,i_min-1
                new_grdptid(i,j) = no_pt
                if(debug) this%nodes(i,j,1) = 1.0
             end do

             do i=i_min, i_max
                new_grdptid(i,j) = this%grdpts_id(
     $               match_table(1)+i,
     $               match_table(2)+j)
             end do

             do i=i_max+1, size(this%nodes,1)
                new_grdptid(i,j) = no_pt
                if(debug) this%nodes(i,j,1) = 1.0
             end do

          end do          

          do j=j_max+1,size(this%nodes,2)
             do i=1, size(this%nodes,1)
                new_grdptid(i,j) = no_pt
                if(debug) this%nodes(i,j,1) = 1.0
             end do
          end do
             
          call MOVE_ALLOC(new_grdptid,this%grdpts_id)
        
        end subroutine reallocate_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine allocating the main tables of the
        !> buffer layer for the first time
        !
        !> @date
        !> 04_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !-------------------------------------------------
        subroutine allocate_main_tables(this)
        
          implicit none

          class(bf_layer_abstract), intent(inout) :: this


          integer(ikind), dimension(3) :: size_nodes


          !< determine the size of the main tables
          !> depending on the localization of the
          !> buffer layer and its alignment with the
          !> interior
          select case(this%localization(1))
          
            case(N,S)
               size_nodes(1) = this%alignment(1,2) - this%alignment(1,1) + 2*bc_size
               size_nodes(2) = 2*bc_size+1

            case(E,W)
               size_nodes(1) = 2*bc_size+1
               size_nodes(2) = this%alignment(2,2) - this%alignment(2,1) + 2*bc_size

            case(N_E,N_W,S_E,S_W)
               size_nodes(1) = 2*bc_size+1
               size_nodes(2) = 2*bc_size+1

            case default
               print *, 'bf_layer_allocate_module'
               print *, 'allocate_main_tables'
               stop 'localization not recognized'

          end select

          size_nodes(3) = ne


          !< allocate the nodes table
          allocate(this%nodes(
     $         size_nodes(1),
     $         size_nodes(2),
     $         size_nodes(3)))


          !< allocate the grdptid table
          allocate(this%grdpts_id(
     $         size_nodes(1),
     $         size_nodes(2)))

        end subroutine allocate_main_tables        

      end module bf_layer_allocate_module

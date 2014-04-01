      module bf_layer_abstract_class

        use parameters_bf_layer, only : no_pt, interior_pt, bc_pt
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input   , only : nx,ny,ne,bc_size
        use parameters_kind    , only : ikind, rkind

        implicit none

        private
        public :: bf_layer_abstract

        logical, parameter :: debug = .true.

        !localization   : N,S,E,W,N_E,N_W,S_E,S_W
        !alignment      : position of the buffer layer compared
        !                 to the interior nodes table
        !                   - alignment(1,1) : i_min
        !                   - alignment(1,2) : i_max
        !                   - alignment(2,1) : j_min
        !                   - alignment(2,2) : j_max
        !
        !initialize     : initialize the localization of the buffer
        !                 layer compared to the interior table
        !allocate_nodes : allocate the nodes for the buffer layer
        type :: bf_layer_abstract

          integer                        :: localization
          integer(ikind), dimension(2,2) :: alignment
          
          real(rkind)   , dimension(:,:,:), allocatable :: nodes
          integer       , dimension(:,:)  , allocatable :: grdpts_id

          contains

          procedure, pass :: ini
          procedure, pass :: allocate_bf_layer
          procedure, pass :: reallocate_bf_layer
          procedure, pass :: identify_and_compute_new_gridpoints

          procedure, pass :: print_nodes
          procedure, pass :: print_grdpts_id
          procedure, pass :: print_sizes

        end type bf_layer_abstract


        contains


        subroutine ini(this,localization)

          implicit none

          class(bf_layer_abstract), intent(inout) :: this
          integer(ikind)          , intent(in)    :: localization
          
          this%localization = localization

        end subroutine ini
        

        subroutine allocate_bf_layer(this, alignment, nodes)

          implicit none
          
          class(bf_layer_abstract)        , intent(inout) :: this
          integer(ikind), dimension(2,2)  , intent(in)    :: alignment
          real(rkind)   , dimension(:,:,:), intent(in)    :: nodes

          integer(ikind), dimension(:,:)  , allocatable   :: list_new_grdpts

          !set the alignment between the interior table
          !and the buffer layer
          this%alignment = alignment

          !allocate the number of gridpoints needed for
          !the buffer layer,
          !copy gridpoints from the interior table,
          !and create the list of new gridpoints
          !that should be computed
          call initial_grdpt_copy_from_interior(
     $         this,
     $         nodes,
     $         list_new_grdpts)

          !compute new gridpoints
          call compute_new_gridpoints(this, list_new_grdpts)

        end subroutine allocate_bf_layer


        subroutine reallocate_bf_layer(this, border_changes)
        
          implicit none

          class(bf_layer_abstract), intent(inout) :: this
          integer, dimension(2,2) , intent(in)    :: border_changes


          integer(ikind) :: i,j,k
          integer(ikind) :: i_min, i_max, j_min, j_max
          integer(ikind) :: i_match, j_match
          integer(ikind) :: i_min_change, i_max_change
          integer(ikind) :: j_min_change, j_max_change
          integer(ikind) :: new_size_x, new_size_y
          real(rkind), dimension(:,:,:), allocatable :: new_nodes
          integer    , dimension(:,:)  , allocatable :: new_grdptid

        
          i_min_change = border_changes(1,1)
          i_max_change = border_changes(1,2)
          j_min_change = border_changes(2,1)
          j_max_change = border_changes(2,2)


          !determine the new alignment between the interior and 
          !buffer layer tables
          select case(this%localization)
            case(N,S,E,W)
               this%alignment(1,1) = this%alignment(1,1) + i_min_change
               this%alignment(1,2) = this%alignment(1,2) + i_max_change
               this%alignment(2,1) = this%alignment(2,1) + j_min_change
               this%alignment(2,2) = this%alignment(2,2) + j_max_change
          end select


          !determine the borders when filling the new table with the
          !old data
          i_min = 1 + max(0,i_min_change) - i_min_change
          i_max = size(this%nodes,1) + min(0,i_max_change) - i_min_change
          j_min = 1 + max(0,j_min_change) - j_min_change
          j_max = size(this%nodes,2) + min(0,j_max_change) - j_min_change

          i_match = i_min_change
          j_match = j_min_change


          !determine the new size of the nodes and grdptid tables
          new_size_x = size(this%nodes,1) - i_min_change + i_max_change
          new_size_y = size(this%nodes,2) - j_min_change + j_max_change
          

          !allocate the new tables
          allocate(new_nodes(new_size_x,new_size_y,ne))
          allocate(new_grdptid(new_size_x,new_size_y))


          !debugging step to check the inputs
          if(debug) then
             
             select case(this%localization)
               case(N)
                  if(j_min_change.ne.0) then
                     stop 'N: j_min_change.ne.0: this is wrong'
                  end if

               case(S)
                  if(j_max_change.ne.0) then
                     stop 'S: j_max_change.ne.0: this is wrong'
                  end if
                  
               case(E)
                  if(i_min_change.ne.0) then
                     stop 'E: i_min_change.ne.0: this is wrong'
                  end if
                  
               case(W)
                  if(i_max_change.ne.0) then
                     stop 'W: i_max_change.ne.0: this is wrong'
                  end if
                  
               case(N_E)
                  if((i_min_change.ne.0).and.(j_min_change.ne.0)) then
                     stop 'NE: change.ne.0: this is wrong'
                  end if
                  
               case(N_W)
                  if((i_max_change.ne.0).and.(j_min_change.ne.0)) then
                     stop 'NW: change.ne.0: this is wrong'
                  end if
                  
               case(S_E)
                  if((i_max_change.ne.0).and.(j_max_change.ne.0)) then
                     stop 'SE: change.ne.0: this is wrong'
                  end if

               case(S_W)
                  if((i_min_change.ne.0).and.(j_max_change.ne.0)) then
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
                   new_nodes(i,j,k) = this%nodes(i_match+i,j_match+j,k)
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
                new_grdptid(i,j) = this%grdpts_id(i_match+i,j_match+j)
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


        subroutine initial_grdpt_copy_from_interior(
     $     this,
     $     nodes,
     $     list_new_grdpts)

          implicit none

          class(bf_layer_abstract)                     , intent(inout) :: this
          real(rkind)   , dimension(:,:,:)             , intent(in)    :: nodes
          integer(ikind), dimension(:,:)  , allocatable, intent(out)   :: list_new_grdpts
          
          integer(ikind), dimension(3) :: size_nodes
          integer(ikind) :: i,j,k
          integer(ikind) :: i_match, j_match
          integer(ikind) :: nb_new_grdpts


          !allocate nodes for buffer layer
          select case(this%localization)
          
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
               print *, 'bf_layer_abstract_class'
               print *, 'initial_grdpt_copy_from_interior'
               stop 'localization not recognized'

          end select

          size_nodes(3) = ne
          allocate(this%nodes(size_nodes(1),size_nodes(2),size_nodes(3)))

          select case(this%localization)

            case(N)

               !copy of grid points from the interior
               i_match = this%alignment(1,1)-bc_size-1
               j_match = ny-2*bc_size
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, size(this%nodes,1)
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,1)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do i=1, size(this%nodes,1)
                  list_new_grdpts(i,1) = i
                  list_new_grdpts(i,2) = 2*bc_size+1
               end do

            case(S)               

               !copy of grid points from the interior
               i_match = this%alignment(1,1)-bc_size-1
               j_match = 0
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, size(this%nodes,1)
                        this%nodes(i,j+1,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,1)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do i=1, size(this%nodes,1)
                  list_new_grdpts(i,1) = i
                  list_new_grdpts(i,2) = 1
               end do
               
            case(E)

               !copy of grid points from the interior
               i_match = nx-2*bc_size
               j_match = this%alignment(2,1)-bc_size-1
               do k=1, ne
                  do j=1, size(this%nodes,2)
                     do i=1, 2*bc_size
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,2)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do j=1, size(this%nodes,2)
                  list_new_grdpts(j,1) = 2*bc_size+1
                  list_new_grdpts(j,2) = j
               end do

            case(W)

               !copy of grid points from the interior
               i_match = 0
               j_match = this%alignment(2,1)-bc_size-1
               do k=1, ne
                  do j=1, size(this%nodes,2)
                     do i=1, 2*bc_size
                        this%nodes(i+1,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do
               
               !determination of the list of new gridpoints
               nb_new_grdpts = size(this%nodes,2)
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do j=1, size(this%nodes,2)
                  list_new_grdpts(j,1) = 1
                  list_new_grdpts(j,2) = j
               end do

            case(N_E)

               !copy of grid points from the interior
               i_match = nx-2*bc_size
               j_match = ny-2*bc_size
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = 4*bc_size + 1
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do j=1, 2*bc_size
                  list_new_grdpts(j,1) = 2*bc_size+1
                  list_new_grdpts(j,2) = j
               end do
               do i=1, 2*bc_size+1
                  list_new_grdpts(2*bc_size+i,1) = i
                  list_new_grdpts(2*bc_size+i,2) = 2*bc_size+1
               end do
               

            case(N_W)

               !copy of grid points from the interior
               i_match = 0
               j_match = ny-2*bc_size
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(1+i,j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = 4*bc_size + 1
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do j=1, 2*bc_size
                  list_new_grdpts(j,1) = 1
                  list_new_grdpts(j,2) = j
               end do
               do i=1, 2*bc_size+1
                  list_new_grdpts(2*bc_size+i,1) = i
                  list_new_grdpts(2*bc_size+i,2) = 2*bc_size+1
               end do

            case(S_E)

               !copy of grid points from the interior
               i_match = nx-2*bc_size
               j_match = 0
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(i,1+j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = 4*bc_size+1
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do i=1, 2*bc_size+1
                  list_new_grdpts(i,1) = i
                  list_new_grdpts(i,2) = 1
               end do
               do j=2, 2*bc_size+1
                  list_new_grdpts(2*bc_size+j,1) = 2*bc_size+1
                  list_new_grdpts(2*bc_size+j,2) = j
               end do
               

            case(S_W)

               !copy of grid points from the interior
               i_match = 0
               j_match = 0
               do k=1, ne
                  do j=1, 2*bc_size
                     do i=1, 2*bc_size
                        this%nodes(1+i,1+j,k) = nodes(i_match+i,j_match+j,k)
                     end do
                  end do
               end do

               !determination of the list of new gridpoints
               nb_new_grdpts = 4*bc_size+1
               allocate(list_new_grdpts(nb_new_grdpts,2))
               do i=1, 2*bc_size+1
                  list_new_grdpts(i,1) = i
                  list_new_grdpts(i,2) = 1
               end do
               do j=2, 2*bc_size+1
                  list_new_grdpts(2*bc_size+j,1) = 1
                  list_new_grdpts(2*bc_size+j,2) = j
               end do

          end select


          !identification of the buffer layer gridpoints
          allocate(this%grdpts_id(size_nodes(1),size_nodes(2)))
          do j=1, size_nodes(2)
             do i=1, size_nodes(1)
                this%grdpts_id(i,j) = bc_pt
             end do
          end do

          select case(this%localization)

            case(N,S)
               j = bc_size+1
               do i=bc_size+1, size_nodes(1)-bc_size
                  this%grdpts_id(i,j) = interior_pt
               end do

            case(E,W)
               i = bc_size+1
               do j=bc_size+1, size_nodes(2)-bc_size
                  this%grdpts_id(i,j) = interior_pt
               end do

            case(N_E,N_W,S_E,S_W)
               i = bc_size+1
               j = bc_size+1
               this%grdpts_id(i,j) = interior_pt

          end select

        end subroutine initial_grdpt_copy_from_interior

      
        subroutine compute_new_gridpoints(this, list_new_gridpts)

          implicit none

          class(bf_layer_abstract)      , intent(inout) :: this
          integer(ikind), dimension(:,:), intent(in)    :: list_new_gridpts
          
          integer(ikind) :: i,j,k,l

          do l=1, size(list_new_gridpts,1)
             
             i = list_new_gridpts(l,1)
             j = list_new_gridpts(l,2)
          
             do k=1, ne
                this%nodes(i,j,k) = 3
             end do

          end do

        end subroutine compute_new_gridpoints


        subroutine compute_new_gridpoint(this, i, j)
          implicit none
          
          class(bf_layer_abstract), intent(inout) :: this
          integer(ikind)          , intent(in)    :: i
          integer(ikind)          , intent(in)    :: j

          integer :: k

          do k=1, ne
             this%nodes(i,j,k) = 3
          end do

        end subroutine compute_new_gridpoint


        subroutine identify_and_compute_new_gridpoints(
     $     this,
     $     selected_grdpts)

          implicit none

          class(bf_layer_abstract)      , intent(inout) :: this
          integer(ikind), dimension(:,:), intent(in)    :: selected_grdpts


          integer(ikind) :: i,j,k
          integer(ikind) :: i_prev, j_prev, i_center, j_center


          !we have a list of gridpoints that should be turned
          !from interior to boundary points. For a point to be
          !considered an interior point, we need to make sure
          !that the grid points I need to compute its time
          !derivatives are available
          !
          !previously, the nodes table of the buffer layer was
          !increased to take into account the space needed for
          !new gridpoints at the boundary
          !
          !in this function, we go through the list of gridpoint
          !whose neighbours need to be tested. As it is possible
          !to reduce the number of tests by considering the
          !previous gridpoint tested, there is a distinction
          !between k=1 and k>1 in the list of gridpoints
          !----------------------------------------------------

          !for the first grid point, there is no simplification
          !possible in checking the neighbours so it is taken
          !outside the next loop
          !----------------------------------------------------
          k = 1
          i_center = selected_grdpts(k,1)
          j_center = selected_grdpts(k,2)
          this%grdpts_id(i_center,j_center) = interior_pt
          do j=-bc_size, bc_size
             do i=-bc_size, bc_size
                call check_gridpoint(this,i,j)
             end do
          end do


          !from the second gridpoint, we reduce the number of
          !neighbours to be tested
          !----------------------------------------------------
          do k=2, size(selected_grdpts,1)

             !update the position of the gridpoint previously
             !tested
             i_prev   = i_center
             j_prev   = j_center
             i_center = selected_grdpts(k,1)
             j_center = selected_grdpts(k,2)

             !update the status of the gridpoint
             this%grdpts_id(i_center,j_center) = interior_pt

             !check its neighbours
             call check_neighbors(this,i_prev,j_prev,i_center,j_center)

          end do             

        end subroutine identify_and_compute_new_gridpoints


        !algorithm checking the neighbours of the
        !gridpoint and considering the previous
        !gridpoints tested not to test twice the
        !same
        !i_prev:   x-index of the previous gridpoint changed
        !j_prev:   y-index of the previous gridpoint changed
        !i_center: x-index of the current gridpoint changed
        !j_center: y-index of the current gridpoint changed
        subroutine check_neighbors(this,i_prev,j_prev,i_center,j_center)

          implicit none

          class(bf_layer_abstract), intent(inout) :: this
          integer(ikind)          , intent(in)    :: i_prev
          integer(ikind)          , intent(in)    :: j_prev
          integer(ikind)          , intent(in)    :: i_center
          integer(ikind)          , intent(in)    :: j_center
          

          integer(ikind) :: min_j, max_j
          integer(ikind) :: i,j


          min_j = min(j_center-j_prev,0)
          max_j = max(j_center-j_prev,0)

          do j=j_center-j_prev, -1
             do i=-bc_size,bc_size
                call check_gridpoint(this, i_center+i, j_prev-bc_size+j)
             end do
          end do

          do j=-bc_size-min_j, bc_size-max_j
             do i=i_center-i_prev,-1
                call check_gridpoint(this, i_prev-bc_size+i, j_center+j)
             end do
          end do

          do j=-bc_size-min_j, bc_size-max_j
             do i=1,i_center-i_prev
                call check_gridpoint(this, i_prev+bc_size+i, j_center+j)
             end do
          end do

          do j=1, j_center-j_prev
             do i=-bc_size,bc_size
                call check_gridpoint(this, i_center+i, j_prev+bc_size+j)
             end do
          end do

        end subroutine check_neighbors
      


        !check if the gridpoint asked exists or not
        !if it does not exist in the buffer layer
        !it has to be computed and updated in the
        !gridpoint ID map
        subroutine check_gridpoint(this,i,j)

          implicit none

          class(bf_layer_abstract), intent(inout) :: this
          integer(ikind)          , intent(in)    :: i,j

          if (this%grdpts_id(i,j).eq.no_pt) then
             call compute_new_gridpoint(this,i,j)
             this%grdpts_id(i,j) = bc_pt
          end if

        end subroutine check_gridpoint


      
        subroutine print_nodes(this,filename)

          implicit none

          class(bf_layer_abstract), intent(in) :: this
          character(*)            , intent(in) :: filename

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
              write(unit=1, iostat=ios) this%nodes
              close(unit=1)
           else
              stop 'file opening pb'
           end if

        end subroutine print_nodes


        subroutine print_grdpts_id(this,filename)

          implicit none
          
          class(bf_layer_abstract), intent(in) :: this
          character(*)            , intent(in) :: filename

          integer :: ios

          open(unit=2,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

          if(ios.eq.0) then
             write(unit=2, iostat=ios) this%grdpts_id
             close(unit=2)
          else
             stop 'file opening pb'
          end if

        end subroutine print_grdpts_id
      

        subroutine print_sizes(this, filename)

          implicit none
          
          class(bf_layer_abstract), intent(in) :: this
          character(*)            , intent(in) :: filename

          integer :: ios

          open(unit=2,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

          if(ios.eq.0) then
             write(unit=2, iostat=ios)
     $            size(this%nodes,1),
     $            size(this%nodes,2),
     $            size(this%nodes,3),
     $            this%alignment(1,1),
     $            this%alignment(1,2),
     $            this%alignment(2,1),
     $            this%alignment(2,2)
             close(unit=2)
          else
             stop 'file opening pb'
          end if

        end subroutine print_sizes

      end module bf_layer_abstract_class

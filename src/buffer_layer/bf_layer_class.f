      !> @file
      !> module encapsulating the buffer layer object where its
      !> internal procedures are initialized
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the buffer layer object where its
      !> internal procedures are initialized
      !
      !> @date
      ! 07_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_class

        use bf_layer_ini_grdptID_module, only : ini_grdpts_id_N,
     $                                          ini_grdpts_id_S,
     $                                          ini_grdpts_id_E,
     $                                          ini_grdpts_id_W,
     $                                          ini_grdpts_id_NE,
     $                                          ini_grdpts_id_NW,
     $                                          ini_grdpts_id_SE,
     $                                          ini_grdpts_id_SW

        use bf_layer_merge_module      , only : merge_bf_layers_N,
     $                                          merge_bf_layers_S,
     $                                          merge_bf_layers_E,
     $                                          merge_bf_layers_W

        use parameters_bf_layer        , only : exchange_pt,
     $                                          bc_pt, bc_interior_pt,
     $                                          interior_pt, no_pt
        use parameters_constant        , only : N,S,E,W,N_W,N_E,S_E,S_W
        use parameters_input           , only : nx,ny,ne,bc_size
        use parameters_kind            , only : ikind, rkind

        private
        public :: bf_layer

        logical, parameter :: debug = .true.


        !> @class bf_layer
        !> class encapsulating the buffer layer which extends the
        !> interior nodes in a particular direction
        !>
        !> @param allocate_bf_layer
        !> allocate the main tables of the buffer layer (nodes,
        !> grdptid)
        !>
        !> @param reallocate_bf_layer
        !> reallocate the main tables of the buffer layer to increase
        !> the buffer layer or remove parts of it
        !>
        !> @param print_nodes
        !> procedure used to print the nodes on a binary file
        !>
        !> @param print_grdpts_id
        !> procedure used to print the grdpt_id on a binary file
        !> 
        !> @param print_sizes
        !> procedure used to print the sizes of the main tables of the
        !> buffer layer
        !---------------------------------------------------------------
        type :: bf_layer

          integer                       , private :: localization
          integer(ikind), dimension(2,2), private :: alignment
          
          real(rkind)   , dimension(:,:,:), allocatable, private :: nodes
          integer       , dimension(:,:)  , allocatable, private :: grdpts_id

          contains

          procedure, pass :: ini

          procedure, pass :: get_localization
          procedure, pass :: get_sizes
          procedure, pass :: set_nodes

          procedure, pass :: get_local_coord
          procedure, pass :: compute_new_grdpts
          procedure, pass :: allocate_bf_layer
          procedure, pass :: reallocate_bf_layer
          procedure, pass :: merge_bf_layer

          procedure, pass, private :: ini_alignment
          procedure, pass, private :: allocate_nodes_and_grdpts_id
          procedure, pass, private :: ini_nodes
          procedure, pass, private :: ini_grdpts_id

          procedure, pass :: print_binary

        end type bf_layer


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine print the nodes table in a binary
        !> file
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
        !>@param localization
        !> localization of the buffer layer
        !>     - localization(1) : main layer (N,S,E,W,N_E,N_W,S_E,S_W)
        !>     - localization(2) : sub-layer
        !--------------------------------------------------------------        
        subroutine ini(this,localization)

          implicit none

          class(bf_layer), intent(inout) :: this
          integer(ikind) , intent(in)    :: localization
          
          this%localization = localization

        end subroutine ini


        !> get the localization attribute
        function get_localization(this)
        
          implicit none

          class(bf_layer), intent(in) :: this
          integer                     :: get_localization
          
          get_localization = this%localization
          
        end function get_localization

      
        !< get the sizes of the nodes and grdpts_id
        !> attributes
        function get_sizes(this)

          class(bf_layer)             , intent(in) :: this
          integer(ikind), dimension(2)             :: get_sizes

          get_sizes(1) = size(this%grdpts_id,1)
          get_sizes(2) = size(this%grdpts_id,2)

        end function get_sizes

      
        !> set the nodes attribute
        subroutine set_nodes(this, nodes)
        
          implicit none

          class(bf_layer)                           , intent(inout) :: this
          real(rkind), dimension(:,:,:), allocatable                :: nodes

          call MOVE_ALLOC(nodes, this%nodes)

        end subroutine set_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine giving the local coordinates in the
        !> buffer layer considering the general coordinates
        !> compared to the interior nodes
        !
        !> @date
        !> 11_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !
        !>@param general_coord
        !> table of integers encapsulating the general
        !> coordinates
        !--------------------------------------------------------------
        function get_local_coord(this, general_coord) result(local_coord)

          implicit none

          class(bf_layer)             , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: general_coord
          integer(ikind), dimension(2)             :: local_coord

          select case(this%localization)
            case(N)
               local_coord(1) = general_coord(1) - this%alignment(1,1) + bc_size + 1
               local_coord(2) = general_coord(2) - ny + bc_size
            case(S)
               local_coord(1) = general_coord(1) - this%alignment(1,1) + bc_size + 1
               local_coord(2) = general_coord(2) + size(this%nodes,2)  - bc_size
            case(E)
               local_coord(1) = general_coord(1) - nx + bc_size
               local_coord(2) = general_coord(2) - this%alignment(2,1) + bc_size + 1
            case(W)
               local_coord(1) = general_coord(1) + size(this%nodes,1) - bc_size
               local_coord(2) = general_coord(2) - this%alignment(2,1) + bc_size + 1
            case(N_E)
               local_coord(1) = general_coord(1) - nx + bc_size
               local_coord(2) = general_coord(2) - ny + bc_size
            case(N_W)
               local_coord(1) = general_coord(1) + size(this%nodes,1) - bc_size
               local_coord(2) = general_coord(2) - ny + bc_size
            case(S_E)
               local_coord(1) = general_coord(1) - nx + bc_size
               local_coord(2) = general_coord(2) + size(this%nodes,2) - bc_size
            case(S_W)
               local_coord(1) = general_coord(1) + size(this%nodes,1) - bc_size
               local_coord(2) = general_coord(2) + size(this%nodes,2) - bc_size
            case default
               print '(''bf_layer_class'')'
               print '(''get_local_coord'')'
               print '(''localization not recognized'',I2)', this%localization
               stop '(was the buffer layer initialized ?)'
          end select               

        end function get_local_coord
        

        !< compute the new grid points corresponding to
        !> the list_new_grdpts asked
        subroutine compute_new_grdpts(this, list_new_grdpts)

          implicit none

          class(bf_layer)               , intent(inout) :: this
          integer(ikind), dimension(:,:), intent(in) :: list_new_grdpts

          
          integer(ikind) :: i,j,l
          integer        :: k

          do l=1, size(list_new_grdpts,1)
             
             i = list_new_grdpts(l,1)
             j = list_new_grdpts(l,2)
          
             do k=1, ne
                this%nodes(i,j,k) = 3
             end do

          end do          

        end subroutine compute_new_grdpts


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
        !> bf_layer object encapsulating the main
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
        subroutine allocate_bf_layer(this, nodes, alignment, neighbors)

          implicit none
          
          class(bf_layer)                 , intent(inout) :: this
          real(rkind)   , dimension(:,:,:), intent(in)    :: nodes
          integer(ikind), dimension(2,2)  , intent(in)    :: alignment
          logical       , dimension(4)    , intent(in)    :: neighbors

          !initialize the alignment
          call this%ini_alignment(alignment)

          !allocate the nodes and the grdpts_id attributes
          call this%allocate_nodes_and_grdpts_id()

          !initialize the nodes
          call this%ini_nodes(nodes)

          !initialize the grdpts_id
          call this%ini_grdpts_id(neighbors)

        end subroutine allocate_bf_layer


        !< reallocate the buffer layer using the
        !> border changes asked
        subroutine reallocate_bf_layer(this, alignment)

          implicit none

          class(bf_layer)        , intent(inout) :: this
          integer, dimension(2,2), intent(in)    :: alignment

          integer(ikind), dimension(2,2) :: border_changes
          integer(ikind) :: i,j
          integer        :: k
          integer(ikind) :: i_min, i_max, j_min, j_max
          integer(ikind) :: new_size_x, new_size_y
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdptid
          integer(ikind), dimension(2) :: match_table


          !< determine the border changes corresponding to the
          !> alignment asked
          border_changes(1,1) = alignment(1,1) - this%alignment(1,1)
          border_changes(2,1) = alignment(2,1) - this%alignment(2,1)
          border_changes(1,2) = alignment(1,2) - this%alignment(1,2)
          border_changes(2,2) = alignment(2,2) - this%alignment(2,2)


          !< determine the new alignment between the interior and 
          !> buffer layer tables
          select case(this%localization)
            case(N)
               border_changes(1,1) = alignment(1,1) - this%alignment(1,1)
               border_changes(2,1) = 0 !alignment(2,1) - this%alignment(2,1)
               border_changes(1,2) = alignment(1,2) - this%alignment(1,2)
               border_changes(2,2) = alignment(2,2) - this%alignment(2,2)

               this%alignment(1,1) = alignment(1,1)
               !this%alignment(2,1) = alignment(2,1)
               this%alignment(1,2) = alignment(1,2)
               this%alignment(2,2) = alignment(2,2)
           
            case(S)
               border_changes(1,1) = alignment(1,1) - this%alignment(1,1)
               border_changes(2,1) = alignment(2,1) - this%alignment(2,1)
               border_changes(1,2) = alignment(1,2) - this%alignment(1,2)
               border_changes(2,2) = 0 !alignment(2,2) - this%alignment(2,2)

               this%alignment(1,1) = alignment(1,1)
               this%alignment(2,1) = alignment(2,1)
               this%alignment(1,2) = alignment(1,2)
               !this%alignment(2,2) = alignment(2,2)

            case(E)
               border_changes(1,1) = 0 !alignment(1,1) - this%alignment(1,1)
               border_changes(2,1) = alignment(2,1) - this%alignment(2,1)
               border_changes(1,2) = alignment(1,2) - this%alignment(1,2)
               border_changes(2,2) = alignment(2,2) - this%alignment(2,2)

               !this%alignment(1,1) = alignment(1,1)
               this%alignment(2,1) = alignment(2,1)
               this%alignment(1,2) = alignment(1,2)
               this%alignment(2,2) = alignment(2,2)

            case(W)
               border_changes(1,1) = alignment(1,1) - this%alignment(1,1)
               border_changes(2,1) = alignment(2,1) - this%alignment(2,1)
               border_changes(1,2) = 0 !alignment(1,2) - this%alignment(1,2)
               border_changes(2,2) = alignment(2,2) - this%alignment(2,2)

               this%alignment(1,1) = alignment(1,1)
               this%alignment(2,1) = alignment(2,1)
               !this%alignment(1,2) = alignment(1,2)
               this%alignment(2,2) = alignment(2,2)

            case(N_E)
               border_changes(1,1) = 0 !alignment(1,1) - this%alignment(1,1)
               border_changes(2,1) = 0 !alignment(2,1) - this%alignment(2,1)
               border_changes(1,2) = alignment(1,2) - this%alignment(1,2)
               border_changes(2,2) = alignment(2,2) - this%alignment(2,2)

               !this%alignment(1,1) = alignment(1,1)
               !this%alignment(2,1) = alignment(2,1)
               this%alignment(1,2) = alignment(1,2)
               this%alignment(2,2) = alignment(2,2)

            case(N_W)
               border_changes(1,1) = alignment(1,1) - this%alignment(1,1)
               border_changes(2,1) = 0 !alignment(2,1) - this%alignment(2,1)
               border_changes(1,2) = 0 !alignment(1,2) - this%alignment(1,2)
               border_changes(2,2) = alignment(2,2) - this%alignment(2,2)

               this%alignment(1,1) = alignment(1,1)
               !this%alignment(2,1) = alignment(2,1)
               !this%alignment(1,2) = alignment(1,2)
               this%alignment(2,2) = alignment(2,2)

            case(S_E)
               border_changes(1,1) = 0 !alignment(1,1) - this%alignment(1,1)
               border_changes(2,1) = alignment(2,1) - this%alignment(2,1)
               border_changes(1,2) = alignment(1,2) - this%alignment(1,2)
               border_changes(2,2) = 0 !alignment(2,2) - this%alignment(2,2)

               !this%alignment(1,1) = alignment(1,1)
               this%alignment(2,1) = alignment(2,1)
               this%alignment(1,2) = alignment(1,2)
               !this%alignment(2,2) = alignment(2,2)

            case(S_W)
               border_changes(1,1) = alignment(1,1) - this%alignment(1,1)
               border_changes(2,1) = alignment(2,1) - this%alignment(2,1)
               border_changes(1,2) = 0 !alignment(1,2) - this%alignment(1,2)
               border_changes(2,2) = 0 !alignment(2,2) - this%alignment(2,2)

               this%alignment(1,1) = alignment(1,1)
               this%alignment(2,1) = alignment(2,1)
               !this%alignment(1,2) = alignment(1,2)
               !this%alignment(2,2) = alignment(2,2)

            case default
               print '(''bf_layer_class'')'
               print '(''reallocate_bf_layer'')'
               print '(''localization not recognized: '')', this%localization
               stop 'wrong initialization'

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
             end do
          end do
          
          do j=j_min, j_max

             do i=1,i_min-1
                new_grdptid(i,j) = no_pt
             end do

             do i=i_min, i_max
                new_grdptid(i,j) = this%grdpts_id(
     $               match_table(1)+i,
     $               match_table(2)+j)
             end do

             do i=i_max+1, size(this%nodes,1)
                new_grdptid(i,j) = no_pt
             end do

          end do          

          do j=j_max+1,size(this%nodes,2)
             do i=1, size(this%nodes,1)
                new_grdptid(i,j) = no_pt
             end do
          end do
             
          call MOVE_ALLOC(new_grdptid,this%grdpts_id)

        end subroutine reallocate_bf_layer


        !< reallocate the buffer layer using the
        !> border changes asked
        subroutine merge_bf_layer(this, bf_layer2, alignment)

          implicit none

          class(bf_layer)        , intent(inout) :: this
          class(bf_layer)        , intent(inout) :: bf_layer2
          integer, dimension(2,2), intent(in)    :: alignment

          !check if the two buffer layers have the same localization
          if(debug) then
             if(this%localization.ne.bf_layer2%localization) then
                print '(''bf_layer_class'')'
                print '(''merge_bf_layer'')'
                print '(''localization do not match'')'
                print '(''loc1: '',I2)', this%localization
                print '(''loc2: '',I2)', bf_layer2%localization
                stop 'check localizations'
             end if
          end if

          !merge sublayers
          select case(this%localization)
            case(N)
               call merge_bf_layers_N(
     $              this%nodes    , bf_layer2%nodes,
     $              this%grdpts_id, bf_layer2%grdpts_id,
     $              this%alignment, bf_layer2%alignment, alignment)
            case(S)
               call merge_bf_layers_S(
     $              this%nodes    , bf_layer2%nodes,
     $              this%grdpts_id, bf_layer2%grdpts_id,
     $              this%alignment, bf_layer2%alignment, alignment)
            case(E)
               call merge_bf_layers_E(
     $              this%nodes    , bf_layer2%nodes,
     $              this%grdpts_id, bf_layer2%grdpts_id,
     $              this%alignment, bf_layer2%alignment, alignment)

            case(W)
               call merge_bf_layers_W(
     $              this%nodes    , bf_layer2%nodes,
     $              this%grdpts_id, bf_layer2%grdpts_id,
     $              this%alignment, bf_layer2%alignment, alignment)

            case default
               print '(''bf_layer_class'')'
               print '(''merge_bf_layer'')'
               print '(''corner buffer layers cannot be merged'')'
               print '(''loc: '',I2)', this%localization
               stop 'check localizations'
          end select

        end subroutine merge_bf_layer



        !< set the alignment between the interior table
        !> and the buffer layer allowing to correspond
        !> the interior grid points and the ones from
        !> the buffer layer
        subroutine ini_alignment(this, alignment)

          implicit none

          class(bf_layer)               , intent(inout) :: this
          integer(ikind), dimension(2,2), intent(in)    :: alignment
          
          select case(this%localization)
            case(N)
               this%alignment(1,1) = alignment(1,1)
               this%alignment(2,1) = ny+1
               this%alignment(1,2) = alignment(1,2)
               this%alignment(2,2) = ny+1
            case(S)
               this%alignment(1,1) = alignment(1,1)
               this%alignment(2,1) = 0
               this%alignment(1,2) = alignment(1,2)
               this%alignment(2,2) = 0
            case(E)
               this%alignment(1,1) = nx+1
               this%alignment(2,1) = alignment(2,1)
               this%alignment(1,2) = nx+1
               this%alignment(2,2) = alignment(2,2)
            case(W)
               this%alignment(1,1) = 0
               this%alignment(2,1) = alignment(2,1)
               this%alignment(1,2) = 0
               this%alignment(2,2) = alignment(2,2)
            case(N_E)
               this%alignment(1,1) = nx+1
               this%alignment(2,1) = ny+1
               this%alignment(1,2) = nx+1
               this%alignment(2,2) = ny+1
            case(N_W)
               this%alignment(1,1) = 0
               this%alignment(2,1) = ny+1
               this%alignment(1,2) = 0
               this%alignment(2,2) = ny+1
            case(S_E)
               this%alignment(1,1) = nx+1
               this%alignment(2,1) = 0
               this%alignment(1,2) = nx+1
               this%alignment(2,2) = 0
            case(S_W)
               this%alignment(1,1) = 0
               this%alignment(2,1) = 0
               this%alignment(1,2) = 0
               this%alignment(2,2) = 0
          end select

        end subroutine ini_alignment


        !< allocate the nodes and the grdpts_id attributes
        !> of the bf_layer object
        subroutine allocate_nodes_and_grdpts_id(this)

          implicit none

          class(bf_layer), intent(inout) :: this


          integer(ikind), dimension(3) :: size_nodes


          !determine the size of the main tables depending on the
          !localization of the buffer layer and its alignment with
          !the interior
          select case(this%localization)
          
            case(N,S)
               size_nodes(1) = this%alignment(1,2) - this%alignment(1,1) + 1 + 2*bc_size
               size_nodes(2) = 2*bc_size+1

            case(E,W)
               size_nodes(1) = 2*bc_size+1
               size_nodes(2) = this%alignment(2,2) - this%alignment(2,1) + 1 + 2*bc_size

            case(N_E,N_W,S_E,S_W)
               size_nodes(1) = 2*bc_size+1
               size_nodes(2) = 2*bc_size+1

            case default
               print '(''bf_layer_allocate_module'')'
               print '(''allocate_main_tables'')'
               stop 'localization not recognized'

          end select

          size_nodes(3) = ne


          !allocate the nodes table
          allocate(this%nodes(
     $         size_nodes(1),
     $         size_nodes(2),
     $         size_nodes(3)))


          !allocate the grdptid table
          allocate(this%grdpts_id(
     $         size_nodes(1),
     $         size_nodes(2)))

        end subroutine allocate_nodes_and_grdpts_id

      
        !< initialize the nodes of the object from the
        !> interior nodes and determine the new gridpoints
        subroutine ini_nodes(
     $     this,
     $     nodes)

          implicit none

          class(bf_layer)                              , intent(inout) :: this
          real(rkind)   , dimension(nx,ny,ne)          , intent(in)    :: nodes

          integer(ikind) :: i,j,k
          integer(ikind) :: i_match, j_match
          integer(ikind) :: nb_new_grdpts
          integer(ikind), dimension(:,:), allocatable :: list_new_grdpts

          !copy of the interior grid points
          !and identification of the new grid
          !points
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

             case default
                print '(''bf_layer_class'')'
                print '(''ini_nodes'')'
                print '(''localization not recognized'')'
                print '(''localization: '', I2)', this%localization
                stop 'change localization'

          end select

          !compute the new grid points
          call this%compute_new_grdpts(list_new_grdpts)

        end subroutine ini_nodes


        !< initialization of the grdpts_id
        subroutine ini_grdpts_id(this, neighbors)

          implicit none

          class(bf_layer)      , intent(inout) :: this
          logical, dimension(4), intent(in)    :: neighbors

          select case(this%localization)

            case(N)
               call ini_grdpts_id_N(this%grdpts_id,neighbors)

            case(S)
               call ini_grdpts_id_S(this%grdpts_id,neighbors)

            case(E)
               call ini_grdpts_id_E(this%grdpts_id,neighbors)

            case(W)
               call ini_grdpts_id_W(this%grdpts_id,neighbors)

            case(N_E)
               call ini_grdpts_id_NE(this%grdpts_id)
               
            case(N_W)
               call ini_grdpts_id_NW(this%grdpts_id)
               
            case(S_E)
               call ini_grdpts_id_SE(this%grdpts_id)
               
            case(S_W)
               call ini_grdpts_id_SW(this%grdpts_id)              
            
            case default
               print '(''bf_layer_class'')'
               print '(''ini_grdpts_id'')'
               print '(''localization not recognized'')'
               print '(''localization: '', I2)', this%localization
               stop 'change localization'
          
          end select

        end subroutine ini_grdpts_id


        subroutine print_binary(
     $     this,
     $     filename_nodes,
     $     filename_grdpts_id,
     $     filename_sizes)

          implicit none

          class(bf_layer), intent(in) :: this
          character(*)   , intent(in) :: filename_nodes
          character(*)   , intent(in) :: filename_grdpts_id
          character(*)   , intent(in) :: filename_sizes

          integer :: ios
          
          !nodes
          open(unit=3,
     $          file=filename_nodes,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=3, iostat=ios) this%nodes
              close(unit=3)
           else
              stop 'file opening pb'
           end if

           !grdpts_id
           open(unit=2,
     $          file=filename_grdpts_id,
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

          !sizes
          open(unit=2,
     $          file=filename_sizes,
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

        end subroutine print_binary


c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> subroutine reallocating the main tables of the
c$$$        !> buffer layer and copying the previous data in
c$$$        !> the new table
c$$$        !
c$$$        !> @date
c$$$        !> 07_04_2013 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> bf_layer_abstract object encapsulating the main
c$$$        !> tables and the integer identifying the
c$$$        !> correspondance between the buffer layer and the
c$$$        !> interior grid points
c$$$        !
c$$$        !>@param border_change
c$$$        !> table giving the border changes of the previous
c$$$        !> allocated table
c$$$        !>     - border_change(1,1) : change of i_min
c$$$        !>     - border_change(1,2) : change of i_max
c$$$        !>     - border_change(2,1) : change of j_min
c$$$        !>     - border_change(2,2) : change of j_max
c$$$        !>
c$$$        !> e.g. : i_min = +1 : the first column of the table is deleted
c$$$        !>        i_max = +1 : a new column is added to the table
c$$$        !>        j_min = -1 : a new line is added to the table
c$$$        !>        j_max = -1 : the last line of the table is deleted
c$$$        !
c$$$        !>@param match_table
c$$$        !> indices allowing to identify correctly the
c$$$        !> previously allocated nodes and the new nodes
c$$$        !>     - match_table(1) : i_match
c$$$        !>     - match_table(2) : j_match
c$$$        !>
c$$$        !> e.g. : new_table(i,j) = prev_table(i_match+i,j_match+j)
c$$$        !--------------------------------------------------------------
c$$$        subroutine reallocate_bf_layer_proc(
c$$$     $     this,
c$$$     $     border_changes,
c$$$     $     nodes,
c$$$     $     match_table)
c$$$        
c$$$          implicit none
c$$$
c$$$          class(bf_layer)                        , intent(inout) :: this
c$$$          integer       , dimension(2,2)         , intent(in)    :: border_changes
c$$$          real(rkind)   , dimension(nx,ny,ne)    , intent(in)    :: nodes
c$$$          integer(ikind), dimension(2), optional , intent(out)   :: match_table
c$$$
c$$$          if(present(match_table)) then
c$$$             call reallocate_bf_layer(
c$$$     $            this,
c$$$     $            border_changes,
c$$$     $            nodes,
c$$$     $            match_table)
c$$$          else
c$$$             call reallocate_bf_layer(
c$$$     $            this,
c$$$     $            border_changes,
c$$$     $            nodes)
c$$$          end if
c$$$
c$$$        end subroutine reallocate_bf_layer_proc
c$$$
c$$$
c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> subroutine print the nodes table in a binary
c$$$        !> file
c$$$        !
c$$$        !> @date
c$$$        !> 07_04_2013 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> bf_layer_abstract object encapsulating the main
c$$$        !> tables and the integer identifying the
c$$$        !> correspondance between the buffer layer and the
c$$$        !> interior grid points
c$$$        !
c$$$        !>@param filename
c$$$        !> name of the binary file where the nodes are
c$$$        !> written
c$$$        !--------------------------------------------------------------
c$$$        subroutine print_nodes_proc(this,filename)
c$$$
c$$$          implicit none
c$$$
c$$$          class(bf_layer), intent(in) :: this
c$$$          character(*)   , intent(in) :: filename
c$$$
c$$$          call print_nodes(this,filename)
c$$$
c$$$        end subroutine print_nodes_proc
c$$$
c$$$
c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> subroutine print the grdpt_id table in a binary
c$$$        !> file
c$$$        !
c$$$        !> @date
c$$$        !> 07_04_2013 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> bf_layer_abstract object encapsulating the main
c$$$        !> tables and the integer identifying the
c$$$        !> correspondance between the buffer layer and the
c$$$        !> interior grid points
c$$$        !
c$$$        !>@param filename
c$$$        !> name of the binary file where the grdpt_id are
c$$$        !> written
c$$$        !--------------------------------------------------------------
c$$$        subroutine print_grdpts_id_proc(this,filename)
c$$$
c$$$          implicit none
c$$$          
c$$$          class(bf_layer), intent(in) :: this
c$$$          character(*)   , intent(in) :: filename
c$$$
c$$$          call print_grdpts_id(this, filename)
c$$$
c$$$        end subroutine print_grdpts_id_proc
c$$$
c$$$
c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> subroutine print the sizes of the main table in
c$$$        !> a binary file
c$$$        !
c$$$        !> @date
c$$$        !> 07_04_2013 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param this
c$$$        !> bf_layer_abstract object encapsulating the main
c$$$        !> tables and the integer identifying the
c$$$        !> correspondance between the buffer layer and the
c$$$        !> interior grid points
c$$$        !
c$$$        !>@param filename
c$$$        !> name of the binary file where the sizes are
c$$$        !> written
c$$$        !--------------------------------------------------------------
c$$$        subroutine print_sizes_proc(this, filename)
c$$$
c$$$          implicit none
c$$$          
c$$$          class(bf_layer), intent(in) :: this
c$$$          character(*)   , intent(in) :: filename
c$$$
c$$$          call print_sizes(this, filename)
c$$$
c$$$        end subroutine print_sizes_proc

      end module bf_layer_class

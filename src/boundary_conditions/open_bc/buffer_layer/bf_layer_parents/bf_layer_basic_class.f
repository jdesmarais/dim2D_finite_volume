      !> @file
      !> module encapsulating the bf_layer object. It encapsulates the
      !> nodes needed when extending the interior domain
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the bf_layer object. It encapsulates the
      !> nodes needed when extending the interior domain
      !
      !> @date
      ! 07_04_2014 - initial version       - J.L. Desmarais
      ! 26_06_2014 - documentation update  - J.L. Desmarais
      ! 16_07_2014 - time integration fcts - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_basic_class

         use bf_layer_extract_module, only :
     $       get_bf_layer_match_table

        use parameters_constant, only :
     $       min_border_type,
     $       max_border_type,
     $       x_direction,
     $       y_direction

        use parameters_input, only :
     $       bc_size,
     $       debug,
     $       nx,
     $       ny,
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        private
        public :: bf_layer_basic


        !> @class bf_layer
        !> class encapsulating the buffer layer which extends the
        !> interior nodes in a definite direction
        !
        !> @param localization
        !> cardinal coordinate identifying the position of the buffer
        !> layer : N,S,E, or W
        !
        !> @param alignment
        !> integer identifying the position of the buffer layer
        !> compared to the interior domain. The coordinates of the
        !> four border points are stored as general coordinates
        !>\image html bf_layer_alignment.png "Buffer layer alignment"
        !>\image latex bf_layer_alignment.eps "Buffer layer alignment"
        !
        !> @param x_map
        !> array containing the coordinates in the x-direction
        !
        !> @param y_map
        !> array containing the coordinates in the y-direction
        !
        !> @param nodes
        !> array where the governing variables are saved at each grid
        !> point
        !
        !> @param grdpts_id
        !> array where the role of each grid point is stored (no_pt,
        !> bc_pt, bc_interior_pt, interior_pt)
        !
        !> @param ini
        !> initialize the buffer layer by setting its cardinal
        !> coordinate
        !
        !> @param get_localization
        !> get the localization attribute
        !
        !> @param get_sizes
        !> get the sizes of the nodes and grdpts_id
        !> attributes
        !
        !> @param set_nodes
        !> set the nodes attribute (only for tests)
        !
        !> @param set_nodes_pt
        !> set a specific nodes_pt (only for tests)
        !
        !> @param set_grdpts_id
        !> set the grdpts_id attribute (only for tests)
        !
        !> @param set_grdpts_id_pt
        !> set a specific grid point role (only for tests)
        !
        !> @param get_alignment
        !> get the position of the buffer layer compared
        !> to the interior domain for one of its border
        !
        !> @param get_alignment_tab
        !> get the alignment attribute
        !
        !> @param get_local_coord
        !> get the indices for the tables in the buffer layer
        !> buffer layer considering the general indices identifying
        !> both interior and buffer layer grid points
        !
        !> @param get_general_to_local_coord_tab
        !> get the table matching the general to the local coordinates
        !> local_coord(i) = general_coord(i) - match_table(i)
        !
        !> @param get_nodes
        !> get the governing variables for a specific grid point
        !> knowing its local coordinates
        !
        !> @param get_nodes_nonlocal
        !> get the governing variables for a specific grid point
        !> and its nearest neighbors knowing its local coordinates
        !
        !> @param get_grdpts_id
        !> get the grdpts_id attribute (only for tests)
        !-------------------------------------------------------------
        type :: bf_layer_basic

          integer                        :: localization
          integer(ikind), dimension(2,2) :: alignment

          real(rkind), dimension(:)    , allocatable :: x_map
          real(rkind), dimension(:)    , allocatable :: y_map
          real(rkind), dimension(:,:,:), allocatable :: nodes
          integer    , dimension(:,:)  , allocatable :: grdpts_id
          contains

          procedure, pass :: ini
                       
          !procedures to modify the main attributes:
          !localization, alignment, grdpts_id, nodes
          procedure, pass :: get_localization
          procedure, pass :: set_localization
          procedure, pass :: get_sizes
          procedure, pass :: set_nodes
          procedure, pass :: set_nodes_pt
          procedure, pass :: set_grdpts_id
          procedure, pass :: set_grdpts_id_pt
          procedure, pass :: get_alignment
          procedure, pass :: get_alignment_tab
          procedure, pass :: set_alignment_tab

          procedure, pass :: copy_grdpts_id_to_temp

          !procedures to localize the buffer layer:
          !coordinates and coordinate maps
          procedure, pass :: get_local_coord
          procedure, pass :: get_general_to_local_coord_tab
          procedure, pass :: get_x_map
          procedure, pass :: get_y_map
          procedure, pass :: set_x_map
          procedure, pass :: set_y_map
          procedure, pass :: get_nodes_array
          procedure, pass :: get_nodes
          procedure, pass :: get_nodes_nonlocal
          procedure, pass :: get_grdpts_id

          !procedures to remove the buffer layer
          procedure, pass :: remove

        end type bf_layer_basic


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the buffer layer by setting its cardinal
        !> coordinate
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param localization
        !> localization of the buffer layer (N,S,E, or W)
        !--------------------------------------------------------------
        subroutine ini(this,localization)

          implicit none

          class(bf_layer_basic), intent(inout) :: this
          integer(ikind)       , intent(in)    :: localization
          
          this%localization = localization

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the localization attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@return
        !> cardinal coordinate identifying the buffer layer position
        !--------------------------------------------------------------
        function get_localization(this)
        
          implicit none

          class(bf_layer_basic), intent(in) :: this
          integer                           :: get_localization
          
          get_localization = this%localization
          
        end function get_localization


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the localization attribute
        !
        !> @date
        !> 20_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param localization
        !> cardinal coordinate identifying the buffer layer position
        !--------------------------------------------------------------
        subroutine set_localization(this, localization)
        
          implicit none

          class(bf_layer_basic), intent(inout) :: this
          integer                              :: localization
          
          this%localization = localization
          
        end subroutine set_localization

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the sizes of the nodes and grdpts_id
        !> attributes
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@return
        !> profile of the grdpts_id attribute
        !--------------------------------------------------------------
        function get_sizes(this)

          class(bf_layer_basic)       , intent(in) :: this
          integer(ikind), dimension(2)             :: get_sizes

          if(allocated(this%grdpts_id)) then
             get_sizes(1) = size(this%grdpts_id,1)
             get_sizes(2) = size(this%grdpts_id,2)
          else
             get_sizes = [0,0]
          end if

        end function get_sizes

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the nodes attribute
        !
        !> @warning
        !> as the nodes is passed by reference, this is a very
        !> unefficient way to change the nodes if no reallocation
        !> is required: this subroutine has only been implemented
        !> for test purposes
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param nodes
        !> array containing the governing variables at each
        !> grid point
        !--------------------------------------------------------------
        subroutine set_nodes(this, nodes)
        
          implicit none

          class(bf_layer_basic)                     , intent(inout) :: this
          real(rkind), dimension(:,:,:), allocatable, intent(inout) :: nodes

          call MOVE_ALLOC(nodes, this%nodes)

        end subroutine set_nodes

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set a specific nodes_pt
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i
        !> x-index of the grid point set
        !
        !>@param j
        !> y-index of the grid point set
        !
        !>@param k
        !> index identifying the type of governing variable set
        !
        !>@param nodes_pt
        !> array containing the governing variables at the specific
        !> grid point
        !--------------------------------------------------------------
        subroutine set_nodes_pt(this,i,j,k,nodes_pt)

          implicit none

          class(bf_layer_basic), intent(inout) :: this
          integer(ikind)       , intent(in)    :: i
          integer(ikind)       , intent(in)    :: j
          integer              , intent(in)    :: k
          real(rkind)          , intent(in)    :: nodes_pt

          this%nodes(i,j,k) = nodes_pt

        end subroutine set_nodes_pt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the grdpts_id attribute
        !
        !> @warning
        !> as the grdpts_id is passed by reference, this is a very
        !> unefficient way to change the grdpts_id if no reallocation
        !> is required: this subroutine has only been implemented for
        !> test purposes
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param grdpts_id
        !> array containing the role of the grid points
        !--------------------------------------------------------------
        subroutine set_grdpts_id(this, grdpts_id)
        
          implicit none

          class(bf_layer_basic)               , intent(inout) :: this
          integer, dimension(:,:), allocatable, intent(inout) :: grdpts_id

          call MOVE_ALLOC(grdpts_id, this%grdpts_id)

        end subroutine set_grdpts_id


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set a specific grid point role
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i
        !> x-index of the grid point set
        !
        !>@param j
        !> y-index of the grid point set
        !
        !>@param var
        !> array containing the governing variables at the specific
        !> grid point
        !--------------------------------------------------------------
        subroutine set_grdpts_id_pt(this,i,j,var)
        
          implicit none

          class(bf_layer_basic), intent(inout) :: this
          integer              , intent(in)    :: i,j
          integer              , intent(in)    :: var

          this%grdpts_id(i,j) = var

        end subroutine set_grdpts_id_pt



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the position of the buffer layer compared
        !> to the interior domain
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param direction
        !> direction in which the alignment is asked (x_direction
        !> or y_direction)
        !
        !>@param border_type
        !> either min or max border
        !
        !>@return
        !> coordinate of the border point
        !--------------------------------------------------------------
        function get_alignment(this, direction, border_type)

          implicit none

          class(bf_layer_basic), intent(in) :: this
          integer              , intent(in) :: direction
          integer              , intent(in) :: border_type
          integer(ikind)                    :: get_alignment

          if(debug) then
             if((direction.ne.x_direction).and.
     $          (direction.ne.y_direction)) then
                print '(''bf_layer_class'')'
                print '(''get_alignment'')'
                print '(''direction not recognized'')'
                print '(''direction: '', I2)', direction
                stop 'modify direction'
             end if

             if((border_type.ne.min_border_type).and.
     $          (border_type.ne.max_border_type)) then
                print '(''bf_layer_class'')'
                print '(''get_alignment'')'
                print '(''border type not recognized'')'
                print '(''border_type: '', I2)', border_type
                stop 'modify border_type'
             end if
          end if

          get_alignment = this%alignment(direction,border_type)

        end function get_alignment


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the alignment attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@return alignment_tab
        !> alignment attribute
        !--------------------------------------------------------------
        function get_alignment_tab(this)

          implicit none

          class(bf_layer_basic)         , intent(in) :: this
          integer(ikind), dimension(2,2)             :: get_alignment_tab

          get_alignment_tab = this%alignment

        end function get_alignment_tab

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the alignment attribute
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param alignment
        !> alignment fror the buffer layer
        !--------------------------------------------------------------
        subroutine set_alignment_tab(this, alignment)

          implicit none

          class(bf_layer_basic)         , intent(inout) :: this
          integer(ikind), dimension(2,2), intent(in)    :: alignment

          this%alignment = alignment

        end subroutine set_alignment_tab


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> create a truncated copy of the the grdpts_id. The copy is
        !> a 3x3 array whose center (2,2) is identified by its general
        !> coordinates cpt_coords
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param cpt_coords
        !> integer identifying the indices of the central point
        !
        !>@param nbc_template
        !> truncated 3x3 array copied of the grdpts_id attribute
        !--------------------------------------------------------------
        subroutine copy_grdpts_id_to_temp(this,
     $     cpt_coords,
     $     nbc_template)

          implicit none

          class(bf_layer_basic)         , intent(in)  :: this
          integer(ikind), dimension(2)  , intent(in)  :: cpt_coords
          integer       , dimension(3,3), intent(out) :: nbc_template

          integer(ikind) :: min_border, max_border
          integer(ikind) :: bf_copy_size_x, bf_copy_size_y
          integer(ikind) :: bf_i_min, tbf_i_min
          integer(ikind) :: bf_j_min, tbf_j_min

          integer(ikind) :: i,j
          
          min_border = max(this%alignment(1,1)-bc_size, cpt_coords(1)-1)
          max_border = min(this%alignment(1,2)+bc_size, cpt_coords(1)+1)

          bf_copy_size_x = max_border-min_border+1
          bf_i_min  = min_border - (this%alignment(1,1) - (bc_size+1))
          tbf_i_min = min_border - (cpt_coords(1) - 2)


          min_border = max(this%alignment(2,1)-bc_size, cpt_coords(2)-1)
          max_border = min(this%alignment(2,2)+bc_size, cpt_coords(2)+1)

          bf_copy_size_y = max_border-min_border+1
          bf_j_min  = min_border - (this%alignment(2,1) - (bc_size+1))
          tbf_j_min = min_border - (cpt_coords(2) - 2)

          
          do j=bf_j_min, bf_j_min+bf_copy_size_y-1
             do i=bf_i_min, bf_i_min+bf_copy_size_x-1
                nbc_template(
     $               i-bf_i_min+tbf_i_min,
     $               j-bf_j_min+tbf_j_min) = this%grdpts_id(i,j)
             end do
          end do

        end subroutine copy_grdpts_id_to_temp


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the indices for the tables in the buffer layer
        !> buffer layer considering the general indices identifying
        !> both interior and buffer layer grid points
        !
        !> @date
        !> 11_04_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param general_coord
        !> table of integers encapsulating the general
        !> coordinates
        !
        !>@return local_coord
        !> table of integers encapsulating the local
        !> coordinates
        !--------------------------------------------------------------
        function get_local_coord(this, general_coord) result(local_coord)

          implicit none

          class(bf_layer_basic)       , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: general_coord
          integer(ikind), dimension(2)             :: local_coord

          integer(ikind), dimension(2) :: match_table

          match_table = this%get_general_to_local_coord_tab()
          
          local_coord(1) = general_coord(1) - match_table(1)
          local_coord(2) = general_coord(2) - match_table(2)

        end function get_local_coord


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the table matching the general to the local coordinates
        !> local_coord(i) = general_coord(i) - match_table(i)
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@return match_table
        !> table of integers giving the conversion between local and
        !> general coordinates
        !--------------------------------------------------------------
        function get_general_to_local_coord_tab(this) result(match_table)

          implicit none

          class(bf_layer_basic), intent(in) :: this
          integer(ikind), dimension(2)      :: match_table
          
          
          match_table = get_bf_layer_match_table(
     $         this%alignment)

        end function get_general_to_local_coord_tab


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the x_map attribute
        !
        !> @date
        !> 03_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param x_map
        !> map with the coordinates along the x-direction
        !--------------------------------------------------------------
        subroutine get_x_map(this, x_map)

          implicit none

          class(bf_layer_basic)                 , intent(in) :: this
          real(rkind), dimension(:), allocatable, intent(out):: x_map
          
          if(allocated(x_map)) then
             deallocate(x_map)
          end if

          allocate(x_map(size(this%x_map,1)))
          x_map(:) = this%x_map(:)

        end subroutine get_x_map


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the y_map attribute
        !
        !> @date
        !> 03_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param y_map
        !> map with the coordinates along the y-direction
        !--------------------------------------------------------------
        subroutine get_y_map(this, y_map)

          implicit none

          class(bf_layer_basic)                 , intent(in) :: this
          real(rkind), dimension(:), allocatable, intent(out):: y_map
          
          if(allocated(y_map)) then
             deallocate(y_map)
          end if

          allocate(y_map(size(this%y_map,1)))
          y_map(:) = this%y_map(:)

        end subroutine get_y_map


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the x_map attribute
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param x_map
        !> map with the coordinates along the x-direction
        !--------------------------------------------------------------
        subroutine set_x_map(this, x_map)

          implicit none

          class(bf_layer_basic)                 , intent(inout) :: this
          real(rkind), dimension(:), allocatable, intent(inout) :: x_map
          
          call MOVE_ALLOC(x_map,this%x_map)

        end subroutine set_x_map


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the y_map attribute
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param y_map
        !> map with the coordinates along the y-direction
        !--------------------------------------------------------------
        subroutine set_y_map(this, y_map)

          implicit none

          class(bf_layer_basic)                 , intent(inout) :: this
          real(rkind), dimension(:), allocatable, intent(inout) :: y_map
          
          call MOVE_ALLOC(y_map,this%y_map)

        end subroutine set_y_map


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the nodes attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param l_coords
        !> local coordinate giving the indices of the grid point
        !> in the buffer layer nodes attribute
        !
        !>@return var
        !> governig variables at a specific grid point
        !--------------------------------------------------------------
        subroutine get_nodes_array(this, nodes)

          implicit none

          class(bf_layer_basic)                     , intent(in)    :: this
          real(rkind), dimension(:,:,:), allocatable, intent(inout) :: nodes

          if(allocated(nodes)) then
             deallocate(nodes)
          end if

          allocate(nodes(size(this%nodes,1),size(this%nodes,2),size(this%nodes,3)))
          nodes(:,:,:) = this%nodes(:,:,:)

        end subroutine get_nodes_array


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the governing variables for a specific grid point
        !> knowing its local coordinates
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param l_coords
        !> local coordinate giving the indices of the grid point
        !> in the buffer layer nodes attribute
        !
        !>@return var
        !> governig variables at a specific grid point
        !--------------------------------------------------------------
        function get_nodes(this, l_coords) result(var)

          implicit none

          class(bf_layer_basic)       , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: l_coords
          real(rkind)   , dimension(ne)            :: var
          
          var = this%nodes(l_coords(1), l_coords(2),:)

        end function get_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the governing variables for a specific grid point
        !> knowing its local coordinates
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param l_coords
        !> local coordinate giving the indices of the grid point
        !> in the buffer layer nodes attribute
        !
        !>@return var
        !> governig variables at a specific grid point
        !--------------------------------------------------------------
        subroutine get_nodes_nonlocal(
     $     this,
     $     l_coords,
     $     x_map_local,
     $     y_map_local,
     $     nodes_local)

          implicit none

          class(bf_layer_basic)            , intent(in)  :: this
          integer(ikind), dimension(2)     , intent(in)  :: l_coords
          real(rkind)   , dimension(5)     , intent(out) :: x_map_local
          real(rkind)   , dimension(5)     , intent(out) :: y_map_local
          real(rkind)   , dimension(5,5,ne), intent(out) :: nodes_local

          x_map_local = this%x_map(l_coords(1)-2:l_coords(1)+2)
          y_map_local = this%y_map(l_coords(2)-2:l_coords(2)+2)
          nodes_local = this%nodes(l_coords(1)-2:l_coords(1)+2,
     $                             l_coords(2)-2:l_coords(2)+2,
     $                             :)

        end subroutine get_nodes_nonlocal


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the grdpts_id attribute
        !
        !> @warning
        !> as the grdpts_id is passed by reference, this is a very
        !> unefficient way to change the grdpts_id if no reallocation
        !> is required: this subroutine has only been implemented for
        !> test purposes
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param grdpts_id
        !> array containing the role of the grid points
        !--------------------------------------------------------------
        subroutine get_grdpts_id(this, grdpts_id)

          implicit none

          class(bf_layer_basic)               , intent(in)  :: this
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id

          integer(ikind), dimension(2) :: sizes
          integer(ikind)               :: i,j

          sizes(1) = size(this%grdpts_id,1)
          sizes(2) = size(this%grdpts_id,2)

          allocate(grdpts_id(sizes(1),sizes(2)))

          do j=1, sizes(2)
             do i=1, sizes(1)
                grdpts_id(i,j) = this%grdpts_id(i,j)
             end do
          end do

        end subroutine get_grdpts_id


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the buffer layer by deallocating the main tables
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine remove(this)

          implicit none

          class(bf_layer_basic), intent(inout) :: this

          if(allocated(this%x_map)) deallocate(this%x_map)
          if(allocated(this%y_map)) deallocate(this%y_map)
          if(allocated(this%nodes)) deallocate(this%nodes)
          if(allocated(this%grdpts_id)) deallocate(this%grdpts_id)

        end subroutine remove

      end module bf_layer_basic_class


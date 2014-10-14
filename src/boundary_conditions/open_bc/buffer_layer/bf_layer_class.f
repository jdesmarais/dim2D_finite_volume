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
      module bf_layer_class

        use bf_compute_class, only :
     $       bf_compute
        
        use bf_layer_errors_module, only :
     $       error_mainlayer_id,
     $       error_diff_mainlayer_id

        use bf_layer_allocate_module, only :
     $       allocate_bf_layer_N,
     $       allocate_bf_layer_S,
     $       allocate_bf_layer_E,
     $       allocate_bf_layer_W

        use bf_layer_reallocate_module, only :
     $       reallocate_bf_layer_N,
     $       reallocate_bf_layer_S,
     $       reallocate_bf_layer_E,
     $       reallocate_bf_layer_W
                                        
        use bf_layer_merge_module, only :
     $       merge_bf_layers_N,
     $       merge_bf_layers_S,
     $       merge_bf_layers_E,
     $       merge_bf_layers_W
        
        use bf_layer_exchange_module, only :
     $       do_grdpts_overlap_along_x_dir,
     $       get_match_indices_for_exchange_with_neighbor1,
     $       get_match_indices_for_exchange_with_neighbor2,
     $       copy_from_bf1_to_bf2

        use bf_layer_nf90_operators_module, only :
     $       print_bf_layer_on_netcdf

        use bf_remove_module, only :
     $       check_if_bf_layer_remains

        use interface_integration_step, only :
     $       timeInt_step_nopt

        use parameters_bf_layer, only :
     $       align_E,
     $       align_N,
     $       align_S,
     $       align_W,
     $       bc_interior_pt,
     $       bf_neighbors,
     $       bf_neighbors_id,
     $       bc_pt,
     $       interior_pt,
     $       no_pt

        use parameters_constant, only :
     $       N,S,E,W,
     $       x_direction, y_direction,
     $       min_border, max_border

        use parameters_input, only :
     $       bc_size,
     $       debug,
     $       nx,
     $       ny,
     $       ne
        
        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq


        private
        public :: bf_layer


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
        !>
        !> @param nodes
        !> array where the governing variables are saved at each grid
        !> point
        !
        !> @param grdpts_id
        !> array where the role of each grid point is stored (no_pt,
        !> bc_pt, bc_interior_pt, interior_pt)
        !
        !> @param shares_grdpts_with_neighbor1
        !> logical identifying whether the buffer layer can exchange
        !> grid points with its neighboring buffer layer of type 1
        !
        !> @param shares_grdpts_with_neighbor2
        !> logical identifying whether the buffer layer can exchange
        !> grid points with its neighboring buffer layer of type 2
        !
        !> @param can_remain
        !> logical identifying whether the buffer layer based on its
        !> grid points at the edge with the interior domain are such
        !> that the buffer layer can be removed
        !
        !> @param bf_compute_used
        !> object containing the intermediate variables needed to compute
        !> the time integration steps
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
        !> @param get_grdpts_id
        !> get the grdpts_id attribute (only for tests)
        !
        !> @param compute_new_grdpts
        !> compute the new grid points corresponding to
        !> the list_new_grdpts asked
        !
        !> @param compute_new_grdpt
        !> compute the new grid point corresponding to
        !> indices asked
        !
        !> @param allocate_bf_layer
        !> allocate the main tables of the buffer layer for
        !> the first time and initialize these tables (nodes, grdptid)
        !> using the data of the interior domain and the boundary
        !> conditions applied at the edges
        !
        !> @param reallocate_bf_layer
        !> reallocate the main tables of the buffer layer
        !> and initialize the new grid points (nodes, grdptid)
        !> using the data of the interior domain and the boundary
        !> conditions applied at the edges
        !
        !> @param merge_bf_layer
        !> merge the main tables of two buffers layer and initialize
        !> the new grid points (nodes, grdptid) using the data of the
        !> interior domain and the boundary conditions applied at the
        !> edges
        !
        !> @param set_neighbor1_share
        !> set whether the buffer layer is sharing grid points
        !> with its neighbor2
        !
        !> @param set_neighbor2_share
        !> set whether the buffer layer is sharing grid points
        !> with its neighbor2
        !
        !> @param can_exchange_with_neighbor1
        !> get the share_grdpts_with_neighbor1 attribute
        !
        !> @param can_exchange_with_neighbor2
        !> get the share_grdpts_with_neighbor2 attribute
        !
        !> @param get_neighbor1_id
        !> get the cardinal coordinate coresponding to the neighboring
        !> buffer layers of type 1
        !
        !> @param get_neighbor2_id
        !> get the cardinal coordinate coresponding to the neighboring
        !> buffer layers of type 2
        !
        !> @param shares_grdpts_along_x_dir_with
        !> check if a neighboring buffer layer (positioned such that
        !> it is either a potential neighboring buffer layer of type
        !> 1 or 2) has indeed grid points in common with another buffer
        !> layer by computing the x-size of the layer to be exchanged
        !
        !> @param copy_from_neighbor1
        !> copy the common layer to the current buffer layer
        !> from its neighboring buffer layer identified as of type 1
        !
        !> @param copy_from_neighbor2
        !> copy the common layer to the current buffer layer
        !> from its neighboring buffer layer identified as of type 2
        !
        !> @param copy_to_neighbor1
        !> copy the common layer from the current buffer layer
        !> to its neighboring buffer layer identified as of type 1
        !
        !> @param copy_to_neighbor2
        !> copy the common layer from the current buffer layer
        !> to its neighboring buffer layer identified as of type 2
        !
        !> @param copy_grdpts_id_to_temp
        !> create a truncated copy of the the grdpts_id. The copy is
        !> a 3x3 array whose center (2,2) is identified by its general
        !> coordinates cpt_coords
        !
        !> @param check_neighboring_bc_interior_pts
        !> check if the grid points neighboring a point identified by
        !> its general coordinates (cpt_coords) are bc_interior_pt,
        !> if so, the points are added to a list of bc_interior_pt
        !
        !> @param update_grdpts_after_increase
        !> turn the grdpts_id identified by general coordinates
        !> from bc_interior_pt to interior_pt and reallocate the
        !> buffer layer such that the neighboring points around it
        !> are allocated. Then compute these new grid points
        !
        !> @param set_remain_status
        !> set the can_remain attribute
        !
        !> @param get_remain_status
        !> get the can_remain attribute
        !
        !> @param should_remain
        !> check the grid points at the edge between the buffer layer
        !> and the interior domain and compute whether they undermine
        !> the open boundary conditions. This determines whether the
        !> buffer layer should be removed or not
        !
        !> @param remove
        !> remove the buffe rlayer by deallocating the main tables
        !
        !> @param print_binary
        !> print the nodes and the grdpts_id attributes
        !> as well as the size of the previous tables in
        !> output binary files
        !
        !> @param print_netcdf
        !> print the nodes and the grdpts_id attributes
        !> on a netcdf file
        !
        !> @param ini_for_comput
        !> initialize the grid size for the computations
        !
        !> @param allocate_before_timeInt
        !> allocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @param deallocate_after_timeInt
        !> deallocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @param compute_time_dev
        !> compute the time derivatives
        !
        !> @param compute_integration_step
        !> compute the integration step
        !
        !> @param get_time_dev
        !> get the time derivatives
        !-------------------------------------------------------------
        type :: bf_layer

          integer                       , private :: localization
          integer(ikind), dimension(2,2), private :: alignment

          real(rkind), dimension(:,:,:), allocatable, private :: nodes
          integer, dimension(:,:)      , allocatable, private :: grdpts_id

          logical, private :: shares_grdpts_with_neighbor1
          logical, private :: shares_grdpts_with_neighbor2

          logical, private :: can_remain

          type(bf_compute) :: bf_compute_used

          contains

          procedure,   pass :: ini
                       
          procedure,   pass :: get_localization
          procedure,   pass :: get_sizes
          procedure,   pass :: set_nodes
          procedure,   pass :: set_nodes_pt
          procedure,   pass :: set_grdpts_id
          procedure,   pass :: set_grdpts_id_pt
          procedure,   pass :: get_alignment
          procedure,   pass :: get_alignment_tab
                       
          procedure,   pass :: get_local_coord
          procedure,   pass :: get_general_to_local_coord_tab
          procedure,   pass :: get_nodes
          procedure,   pass :: get_grdpts_id
          
          procedure,   pass :: compute_new_grdpts
          procedure, nopass :: compute_new_grdpt

          procedure,   pass :: allocate_bf_layer
          procedure,   pass :: reallocate_bf_layer
          procedure,   pass :: merge_bf_layer

          procedure,   pass :: set_neighbor1_share
          procedure,   pass :: set_neighbor2_share
          procedure,   pass :: can_exchange_with_neighbor1
          procedure,   pass :: can_exchange_with_neighbor2
          procedure,   pass :: get_neighbor1_id
          procedure,   pass :: get_neighbor2_id
          procedure,   pass :: shares_grdpts_along_x_dir_with

          procedure,   pass :: copy_from_neighbor1
          procedure,   pass :: copy_from_neighbor2
          procedure,   pass :: copy_to_neighbor1
          procedure,   pass :: copy_to_neighbor2

          procedure,   pass :: copy_grdpts_id_to_temp
          procedure,   pass :: check_neighboring_bc_interior_pts

          procedure,   pass :: update_grdpts_after_increase
          
          procedure,   pass :: set_remain_status
          procedure,   pass :: get_remain_status
          procedure,   pass :: should_remain

          procedure,   pass :: remove
          
          procedure,   pass :: print_binary
          procedure,   pass :: print_netcdf

          
          !for time integration
          procedure,   pass, private :: ini_for_comput
          procedure,   pass          :: allocate_before_timeInt
          procedure,   pass          :: deallocate_after_timeInt
          procedure,   pass          :: compute_time_dev
          procedure,   pass          :: compute_integration_step
          procedure,   pass          :: get_time_dev !only for tests

        end type bf_layer


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
        subroutine ini(this,localization,dx,dy)

          implicit none

          class(bf_layer), intent(inout) :: this
          integer(ikind) , intent(in)    :: localization
          real(rkind)    , intent(in)    :: dx
          real(rkind)    , intent(in)    :: dy
          
          this%localization = localization

          call this%ini_for_comput(dx,dy)

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

          class(bf_layer), intent(in) :: this
          integer                     :: get_localization
          
          get_localization = this%localization
          
        end function get_localization

      
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

          class(bf_layer)             , intent(in) :: this
          integer(ikind), dimension(2)             :: get_sizes

          get_sizes(1) = size(this%grdpts_id,1)
          get_sizes(2) = size(this%grdpts_id,2)

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

          class(bf_layer)                           , intent(inout) :: this
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

          class(bf_layer), intent(inout) :: this
          integer(ikind) , intent(in)    :: i
          integer(ikind) , intent(in)    :: j
          integer        , intent(in)    :: k
          real(rkind)    , intent(in)    :: nodes_pt

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

          class(bf_layer)                     , intent(inout) :: this
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

          class(bf_layer), intent(inout) :: this
          integer        , intent(in)    :: i,j
          integer        , intent(in)    :: var

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

          class(bf_layer), intent(in) :: this
          integer        , intent(in) :: direction
          integer        , intent(in) :: border_type
          integer(ikind)              :: get_alignment

          if(debug) then
             if((direction.ne.x_direction).and.
     $          (direction.ne.y_direction)) then
                print '(''bf_layer_class'')'
                print '(''get_alignment'')'
                print '(''direction not recognized'')'
                print '(''direction: '', I2)', direction
                stop 'modify direction'
             end if

             if((border_type.ne.min_border).and.
     $          (border_type.ne.max_border)) then
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

          class(bf_layer)                , intent(in) :: this
          integer(ikind) , dimension(2,2)             :: get_alignment_tab

          get_alignment_tab = this%alignment

        end function get_alignment_tab


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

          class(bf_layer)             , intent(in) :: this
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

          class(bf_layer), intent(in)  :: this
          integer(ikind), dimension(2) :: match_table
          
          select case(this%localization)
            case(N)
               match_table(1) = this%alignment(1,1) - bc_size - 1
               match_table(2) = ny - 2*bc_size
            case(S)
               match_table(1) = this%alignment(1,1) - bc_size - 1
               !match_table(2) = -size(this%nodes,2) + bc_size
               match_table(2) = this%alignment(2,1) - bc_size - 1
            case(E)
               match_table(1) = nx - 2*bc_size
               match_table(2) = this%alignment(2,1) - bc_size - 1
            case(W)
               !match_table(1) = -size(this%nodes,1) + bc_size
               match_table(1) = this%alignment(1,1) - bc_size - 1
               match_table(2) = this%alignment(2,1) - bc_size - 1
           case default
              call error_mainlayer_id(
     $             'bf_layer_class.f',
     $             'get_general_to_local_coord_tab',
     $             this%localization)
          end select

        end function get_general_to_local_coord_tab


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

          class(bf_layer)             , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: l_coords
          real(rkind)   , dimension(ne)            :: var
          
          var = this%nodes(l_coords(1), l_coords(2),:)

        end function get_nodes


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

          class(bf_layer)                     , intent(in)  :: this
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
        !> compute the new grid points corresponding to
        !> the list_new_grdpts asked
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param list_new_grdpts
        !> array containing the indices of the new grid points
        !> to be computed
        !--------------------------------------------------------------
        subroutine compute_new_grdpts(this, list_new_grdpts)

          implicit none

          class(bf_layer)               , intent(inout) :: this
          integer(ikind), dimension(:,:), intent(in) :: list_new_grdpts

          
          integer(ikind) :: i,j,l
          integer        :: k

          do l=1, size(list_new_grdpts,1)
             
             i = list_new_grdpts(1,l)
             j = list_new_grdpts(2,l)
          
             do k=1, ne
                this%nodes(i,j,k) = 3
             end do

          end do          

        end subroutine compute_new_grdpts


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid point corresponding to
        !> indices asked
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i
        !> x-index identifying the grid point in the nodes array
        !
        !>@param j
        !> y-index identifying the grid point in the nodes array
        !--------------------------------------------------------------
        subroutine compute_new_grdpt(nodes,i,j)

          implicit none

          real(rkind), dimension(:,:,:), intent(out) :: nodes
          integer(ikind)               , intent(in)  :: i
          integer(ikind)               , intent(in)  :: j
          
          integer :: k

          do k=1, ne
             nodes(i,j,k) = 3
          end do

        end subroutine compute_new_grdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate the main tables of the buffer layer for
        !> the first time and initialize these tables (nodes, grdptid)
        !> using the data of the interior domain and the boundary
        !> conditions applied at the edges
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@param alignment
        !> table of integers characterizing the
        !> correspondance between the interior grid points
        !> and the buffer layer elements
        !--------------------------------------------------------------
        subroutine allocate_bf_layer(this, nodes, alignment)

          implicit none
          
          class(bf_layer)                 , intent(inout) :: this
          real(rkind)   , dimension(:,:,:), intent(in)    :: nodes
          integer(ikind), dimension(2,2)  , intent(in)    :: alignment


          select case(this%localization)
            case(N)
               call allocate_bf_layer_N(
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)
            case(S)
               call allocate_bf_layer_S(
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)
            case(E)
               call allocate_bf_layer_E(
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)
            case(W)
               call allocate_bf_layer_W(
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)
            case default
               print '(''bf_layer_class'')'
               print '(''allocate_bf_layer'')'
               print '(''localization not recognized'')'
               print '(''localization:'',I2)', this%localization
          end select

        end subroutine allocate_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> reallocate the main tables of the buffer layer
        !> and initialize the new grid points (nodes, grdptid)
        !> using the data of the interior domain and the boundary
        !> conditions applied at the edges
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@param alignment
        !> table of integers characterizing the
        !> correspondance between the interior grid points
        !> and the buffer layer elements
        !--------------------------------------------------------------
        subroutine reallocate_bf_layer(this, nodes, alignment)

          implicit none

          class(bf_layer)                 , intent(inout) :: this
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer    , dimension(2,2)     , intent(in)    :: alignment

          select case(this%localization)

            case(N)
               call reallocate_bf_layer_N(
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(S)
               call reallocate_bf_layer_S(
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(E)
               call reallocate_bf_layer_E(
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(W)
               call reallocate_bf_layer_W(
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

          end select

        end subroutine reallocate_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> merge the main tables of two buffers layer and initialize
        !> the new grid points (nodes, grdptid) using the data of the
        !> interior domain and the boundary conditions applied at the
        !> edges
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param bf_layer2
        !> second bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@param alignment
        !> table of integers characterizing the
        !> correspondance between the interior grid points
        !> and the buffer layer elements
        !--------------------------------------------------------------
        subroutine merge_bf_layer(this, bf_layer2, nodes, alignment)

          implicit none

          class(bf_layer)                               , intent(inout) :: this
          class(bf_layer)                               , intent(inout) :: bf_layer2
          real(rkind)    , dimension(nx,ny,ne)          , intent(in)    :: nodes
          integer(ikind) , dimension(2,2)     , optional, intent(in)    :: alignment

          !check if the two buffer layers have the same localization
          if(debug) then
             if(this%localization.ne.bf_layer2%localization) then
                call error_diff_mainlayer_id(
     $               'bf_layer_class.f',
     $               'merge_bf_layer',
     $               this%localization,
     $               bf_layer2%localization)
             end if
          end if

          !merge sublayers
          select case(this%localization)
            case(N)
               if(present(alignment)) then
                  call merge_bf_layers_N(
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment, alignment)
               else
                  call merge_bf_layers_N(
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment)
               end if

            case(S)
               if(present(alignment)) then
                  call merge_bf_layers_S(
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment, alignment)
               else
                  call merge_bf_layers_S(
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment)
               end if

            case(E)
               if(present(alignment)) then
                  call merge_bf_layers_E(
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment, alignment)
               else
                  call merge_bf_layers_E(
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment)
               end if

            case(W)
               if(present(alignment)) then
                  call merge_bf_layers_W(
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment, alignment)
               else
                  call merge_bf_layers_W(
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment)
               end if

            case default
               call error_mainlayer_id(
     $              'bf_layer_class.f',
     $              'merge_bf_layer',
     $              this%localization)
          end select

        end subroutine merge_bf_layer        


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the shares_grdpts_with_neighbor1 attribute
        !> or let the program compute whether the buffer layer
        !> may be exchanging grid points with the neighboring
        !> buffer layers of type 1
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param neighbor1_share
        !> value set for the shares_grdpts_with_neighbor1 attribute
        !--------------------------------------------------------------
        subroutine set_neighbor1_share(this, neighbor1_share)
          
          implicit none
          
          class(bf_layer)          , intent(inout) :: this
          logical        , optional, intent(in)    :: neighbor1_share

          if(present(neighbor1_share)) then
             this%shares_grdpts_with_neighbor1 = neighbor1_share
          else
             select case(this%localization)
               case(N,S)
                  this%shares_grdpts_with_neighbor1 = this%alignment(1,1).le.(align_W+bc_size)
               case(E,W)
                  this%shares_grdpts_with_neighbor1 = this%alignment(2,1).le.(align_S+bc_size)
               case default
                  call error_mainlayer_id(
     $                 'nbf_interface_class.f',
     $                 'share_grdpts_with_neighbor1',
     $                 this%localization)
             end select
          end if

        end subroutine set_neighbor1_share


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the shares_grdpts_with_neighbor2 attribute
        !> or let the program compute whether the buffer layer
        !> may be exchanging grid points with the neighboring
        !> buffer layers of type 2
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param neighbor2_share
        !> value set for the shares_grdpts_with_neighbor1 attribute
        !--------------------------------------------------------------
        subroutine set_neighbor2_share(this, neighbor2_share)
          
          implicit none
          
          class(bf_layer)          , intent(inout) :: this
          logical        , optional, intent(in)    :: neighbor2_share

          if(present(neighbor2_share)) then
             this%shares_grdpts_with_neighbor2 = neighbor2_share
          else
             select case(this%localization)
               case(N,S)
                  this%shares_grdpts_with_neighbor2 = this%alignment(1,2).ge.(align_E-bc_size)
               case(E,W)
                  this%shares_grdpts_with_neighbor2 = this%alignment(2,2).ge.(align_N-bc_size)
               case default
                  call error_mainlayer_id(
     $                 'nbf_interface_class.f',
     $                 'share_grdpts_with_neighbor2',
     $                 this%localization)
             end select
          end if

        end subroutine set_neighbor2_share


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the share_grdpts_with_neighbor1 attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param can_exchange
        !> share_grdpts_with_neighbor1 attribute
        !--------------------------------------------------------------
        function can_exchange_with_neighbor1(this) result(can_exchange)

          implicit none

          class(bf_layer), intent(in) :: this
          logical                     :: can_exchange

          can_exchange = this%shares_grdpts_with_neighbor1

        end function can_exchange_with_neighbor1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the share_grdpts_with_neighbor1 attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param can_exchange
        !> share_grdpts_with_neighbor1 attribute
        !--------------------------------------------------------------
        function can_exchange_with_neighbor2(this) result(can_exchange)

          implicit none

          class(bf_layer), intent(in) :: this
          logical                     :: can_exchange

          can_exchange = this%shares_grdpts_with_neighbor2

        end function can_exchange_with_neighbor2


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the cardinal coordinate coresponding to the neighboring
        !> buffer layers of type 1
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param neighbor1_id
        !> cardinal coordinate corresponding to neighboring buffer layers
        !> of type 1
        !
        !>@param neighbor_index
        !> index identifying if the original buffer layer is a neighbor
        !> of type 1 or 2 for the neighboring buffer layer
        !--------------------------------------------------------------
        subroutine get_neighbor1_id(this, neighbor1_id, neighbor_index)

          implicit none

          class(bf_layer)          , intent(in) :: this
          integer                  , intent(out):: neighbor1_id
          integer        , optional, intent(out):: neighbor_index


          neighbor1_id    = bf_neighbors(this%localization,1)
          if(present(neighbor_index)) then
             neighbor_index  = bf_neighbors_id(this%localization)
          end if

        end subroutine get_neighbor1_id


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the cardinal coordinate coresponding to the neighboring
        !> buffer layers of type 2
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param neighbor2_id
        !> cardinal coordinate corresponding to neighboring buffer layers
        !> of type 2
        !
        !>@param neighbor_index
        !> index identifying if the original buffer layer is a neighbor
        !> of type 1 or 2 for the neighboring buffer layer
        !--------------------------------------------------------------
        subroutine get_neighbor2_id(this, neighbor2_id, neighbor_index)

          implicit none

          class(bf_layer)          , intent(in) :: this
          integer                  , intent(out):: neighbor2_id
          integer        , optional, intent(out):: neighbor_index

          neighbor2_id    = bf_neighbors(this%localization,2)

          if(present(neighbor_index)) then
             neighbor_index  = bf_neighbors_id(this%localization)
          end if

        end subroutine get_neighbor2_id


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if a neighboring buffer layer
        !> (positioned such that it is either a potential
        !> neighboring buffer layer of type 1 or 2)
        !> has indeed grid points in common with another
        !> buffer layer by computing the x-size of the
        !> layer to be exchanged
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param neighbor
        !> second bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@return share
        !> logical identifying if there are grid points in common
        !> along the x-direction
        !--------------------------------------------------------------
        function shares_grdpts_along_x_dir_with(this, neighbor)
     $     result(share)

          implicit none

          class(bf_layer), intent(in) :: this
          class(bf_layer), intent(in) :: neighbor
          logical                     :: share

          share = do_grdpts_overlap_along_x_dir(
     $         this%alignment, neighbor%alignment)
          
        end function shares_grdpts_along_x_dir_with


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layer to the current buffer layer
        !> from its neighboring buffer layer identified as of type 1
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param neighbor1
        !> second bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine copy_from_neighbor1(this,neighbor1)

          implicit none

          class(bf_layer), intent(inout) :: this
          class(bf_layer), intent(in)    :: neighbor1

          integer(ikind) :: bf_i_min, bf_j_min
          integer(ikind) :: nbf_i_min, nbf_j_min
          integer(ikind) :: bf_copy_size_x, bf_copy_size_y
          
          !get the indices for the match between the tables
          !of the current buffer layer and the neighbor1
          call get_match_indices_for_exchange_with_neighbor1(
     $         this%localization,
     $         this%alignment, size(this%nodes,2),
     $         neighbor1%alignment, size(neighbor1%nodes,2),
     $         bf_i_min, bf_j_min, nbf_i_min, nbf_j_min,
     $         bf_copy_size_x, bf_copy_size_y)
          
          !copy from neighbor1 to the current buffer layer
          call copy_from_bf1_to_bf2(
     $         nbf_i_min, nbf_j_min, bf_i_min, bf_j_min,
     $         bf_copy_size_x, bf_copy_size_y,
     $         neighbor1%nodes, neighbor1%grdpts_id,
     $         this%nodes, this%grdpts_id)

        end subroutine copy_from_neighbor1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layer from the current buffer layer
        !> to its neighboring buffer layer identified as of type 1
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param neighbor1
        !> second bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine copy_to_neighbor1(this, neighbor1)

          implicit none

          class(bf_layer), intent(in)    :: this
          class(bf_layer), intent(inout) :: neighbor1

          integer(ikind) :: bf_i_min, bf_j_min
          integer(ikind) :: nbf_i_min, nbf_j_min
          integer(ikind) :: bf_copy_size_x, bf_copy_size_y
          
          !get the indices for the match between the tables
          !of the current buffer layer and the neighbor1
          call get_match_indices_for_exchange_with_neighbor1(
     $         this%localization,
     $         this%alignment, size(this%nodes,2),
     $         neighbor1%alignment, size(neighbor1%nodes,2),
     $         bf_i_min, bf_j_min,
     $         nbf_i_min, nbf_j_min,
     $         bf_copy_size_x, bf_copy_size_y)
          
          !copy from the current buffer layer to neighbor1
          call copy_from_bf1_to_bf2(
     $         bf_i_min, bf_j_min, nbf_i_min, nbf_j_min,
     $         bf_copy_size_x, bf_copy_size_y,
     $         this%nodes, this%grdpts_id,
     $         neighbor1%nodes, neighbor1%grdpts_id)

        end subroutine copy_to_neighbor1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layer to the current buffer layer
        !> from its neighboring buffer layer identified as of type 2
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param neighbor2
        !> second bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine copy_from_neighbor2(this, neighbor2)

          implicit none

          class(bf_layer), intent(inout) :: this
          class(bf_layer), intent(in)    :: neighbor2

          integer(ikind) :: bf_i_min, bf_j_min
          integer(ikind) :: nbf_i_min, nbf_j_min
          integer(ikind) :: bf_copy_size_x, bf_copy_size_y
          
          !get the indices for the match between the tables
          !of the current buffer layer and the neighbor1
          call get_match_indices_for_exchange_with_neighbor2(
     $         this%localization,
     $         this%alignment, size(this%nodes,2),
     $         neighbor2%alignment, size(neighbor2%nodes,2),
     $         bf_i_min, bf_j_min, nbf_i_min, nbf_j_min,
     $         bf_copy_size_x, bf_copy_size_y)
          
          !copy from neighbor1 to the current buffer layer
          call copy_from_bf1_to_bf2(
     $         nbf_i_min, nbf_j_min, bf_i_min, bf_j_min,
     $         bf_copy_size_x, bf_copy_size_y,
     $         neighbor2%nodes, neighbor2%grdpts_id,
     $         this%nodes, this%grdpts_id)

        end subroutine copy_from_neighbor2


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layer from the current buffer layer
        !> to its neighboring buffer layer identified as of type 2
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param neighbor2
        !> second bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine copy_to_neighbor2(this, neighbor2)

          implicit none

          class(bf_layer), intent(in)    :: this
          class(bf_layer), intent(inout) :: neighbor2

          integer(ikind) :: bf_i_min, bf_j_min
          integer(ikind) :: nbf_i_min, nbf_j_min
          integer(ikind) :: bf_copy_size_x, bf_copy_size_y
          
          if(this%shares_grdpts_with_neighbor2) then
             
             !get the indices for the match between the tables
             !of the current buffer layer and the neighbor1
             call get_match_indices_for_exchange_with_neighbor2(
     $            this%localization,
     $            this%alignment, size(this%nodes,2),
     $            neighbor2%alignment, size(neighbor2%nodes,2),
     $            bf_i_min, bf_j_min, nbf_i_min, nbf_j_min,
     $            bf_copy_size_x, bf_copy_size_y)
            
             !copy from the current buffer layer to neighbor1
             call copy_from_bf1_to_bf2(
     $            bf_i_min, bf_j_min, nbf_i_min, nbf_j_min,
     $            bf_copy_size_x, bf_copy_size_y,
     $            this%nodes, this%grdpts_id,
     $            neighbor2%nodes, neighbor2%grdpts_id)
          end if

        end subroutine copy_to_neighbor2


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

          class(bf_layer)               , intent(in)  :: this
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
        !> check if the grid points neighboring a point identified by its
        !> general coordinates (cpt_coords) are bc_interior_pt, if so,
        !> the points are added to a list of bc_interior_pt
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i_prev
        !> integer identifying the x-index of the previous central point
        !> investigated
        !
        !>@param j_prev
        !> integer identifying the y-index of the previous central point
        !> investigated
        !
        !>@param i_center
        !> integer identifying the x-index of the central point
        !> investigated
        !
        !>@param j_central
        !> integer identifying the y-index of the central point
        !> investigated
        !
        !>@param nb_mgrdpts
        !> number of grid points in the list
        !>
        !>@param mgrdpts
        !> list of indices identifying the bc_interior_pt surrounding
        !> the central grid point
        !--------------------------------------------------------------
        subroutine check_neighboring_bc_interior_pts(
     $     this,
     $     i_prev, j_prev,
     $     i_center, j_center,
     $     nb_mgrdpts,
     $     mgrdpts)

          implicit none

          class(bf_layer)                , intent(in)    :: this
          integer(ikind)                 , intent(in)    :: i_prev
          integer(ikind)                 , intent(in)    :: j_prev
          integer(ikind)                 , intent(in)    :: i_center
          integer(ikind)                 , intent(in)    :: j_center
          integer                        , intent(inout) :: nb_mgrdpts
          integer(ikind) , dimension(:,:), intent(out)   :: mgrdpts


          !radius for the search of bc_interior_pt around the
          !central point identified by (i_center, j_center)
          integer, parameter :: search_r = 1

          integer(ikind), dimension(2) :: match_table
          integer(ikind) :: min_j, max_j
          integer(ikind) :: size_x, size_y
          integer(ikind) :: i,j

          !get the match table converting the general coords
          !into local coords
          match_table = get_general_to_local_coord_tab(this)

          !get the borders of the loops
          min_j = min(j_center-j_prev,0)
          max_j = max(j_center-j_prev,0)

          size_x = size(this%grdpts_id,1)
          size_y = size(this%grdpts_id,2)


          !1.
          do j=max(1,j_center-search_r-match_table(2)),
     $         min(size_y, j_center+search_r-match_table(2), j_prev-search_r-1-match_table(2))

             do i=max(1,i_center-search_r-match_table(1)),
     $            min(size_x, i_center+search_r-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               this%grdpts_id,
     $               nb_mgrdpts,
     $               mgrdpts)
                
             end do
          end do

          
          !2.
          do j=max(1,j_center-search_r-min_j-match_table(2)),
     $         min(size_y, j_center+search_r-max_j-match_table(2))

             do i=max(1,i_center-search_r-match_table(1)),
     $            min(size_x,i_center+bc_size-match_table(1), i_prev-search_r-1-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               this%grdpts_id,
     $               nb_mgrdpts,
     $               mgrdpts)

             end do
          end do


          !3.
          do j=max(1,j_center-search_r-min_j-match_table(2)),
     $         min(size_y,j_center+search_r-max_j-match_table(2))

             do i=max(1,i_center-search_r-match_table(1),i_prev+search_r+1-match_table(1)),
     $            min(size_x,i_center+search_r-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               this%grdpts_id,
     $               nb_mgrdpts,
     $               mgrdpts)
                
             end do
          end do


          !4.
          do j=max(1,j_center-search_r-match_table(2),j_prev+search_r+1-match_table(2)),
     $         min(size_y,j_center+search_r-match_table(2))

             do i=max(1,i_center-search_r-match_table(1)),
     $            min(size_x,i_center+search_r-match_table(1))
                
                call check_bc_interior_pt(
     $               i,j,
     $               match_table,
     $               this%grdpts_id,
     $               nb_mgrdpts,
     $               mgrdpts)
                
             end do
          end do

        end subroutine check_neighboring_bc_interior_pts


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the grid point tested is a bc_interior_pt
        !> and if so save the general coordinates of the grid point
        !> in mgrdpts
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param i
        !> integer identifying the x-index of the central point
        !> investigated
        !
        !>@param j
        !> integer identifying the y-index of the central point
        !> investigated
        !
        !>@param match_table
        !> table converting local coordinate into general coordinates
        !
        !>@param grdpts_id
        !> table where the role of the grid points is stored
        !
        !>@param nb_mgrdpts
        !> number of grid points in the list
        !>
        !>@param mgrdpts
        !> list of indices identifying the bc_interior_pt surrounding
        !> the central grid point
        !--------------------------------------------------------------
        subroutine check_bc_interior_pt(
     $     i,j,
     $     match_table,
     $     grdpts_id,
     $     nb_mgrdpts,
     $     mgrdpts)

          implicit none

          integer(ikind)                , intent(in)    :: i,j
          integer(ikind), dimension(2)  , intent(in)    :: match_table
          integer       , dimension(:,:), intent(in)    :: grdpts_id
          integer                       , intent(inout) :: nb_mgrdpts
          integer(ikind), dimension(:,:), intent(out)   :: mgrdpts

          if(grdpts_id(i,j).eq.bc_interior_pt) then

             nb_mgrdpts = nb_mgrdpts+1
             mgrdpts(1,nb_mgrdpts) = i+match_table(1)
             mgrdpts(2,nb_mgrdpts) = j+match_table(2)
             
          end if

        end subroutine check_bc_interior_pt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> turn the grdpts_id identified by general coordinates
        !> from bc_interior_pt to interior_pt and reallocate the buffer
        !> layer such that the neighboring points around it are allocated.
        !> Then compute these new grid points.
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param selected_grdpts
        !> list containing the general coordinates of the grid points to be
        !> turned from bc_interior_pt to interior_pt
        !--------------------------------------------------------------
        subroutine update_grdpts_after_increase(this, selected_grdpts)

          implicit none

          class(bf_layer)               , intent(inout) :: this
          integer(ikind), dimension(:,:), intent(in)    :: selected_grdpts

          integer(ikind), dimension(2) :: match_table

          match_table = this%get_general_to_local_coord_tab()
          
          call update_bc_interior_pt_to_interior_pt(
     $         this%grdpts_id,
     $         this%nodes,
     $         selected_grdpts,
     $         match_table)

        end subroutine update_grdpts_after_increase


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the grid points of the buffer
        !> layer: bc_interior_pt were identified by the
        !> detectors to be turned into interior pt. Therefore,
        !> it is required to make sure that the neighboring
        !> points needed by the normal stencil exist otherwise
        !> they are added to the buffer layer. The new grid points
        !> are then computed
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@param selected_grdpts
        !> list of bc_interior_pt grid points that should be turned
        !> into interior_pt and therefore need neighboring grid points
        !
        !>@param match_table
        !> (i_match,j_match) needed to identify correctly the grid
        !> points activated by the detectors even if they were
        !> identified before the reallocation
        !--------------------------------------------------------------
        subroutine update_bc_interior_pt_to_interior_pt(
     $       grdpts_id,
     $       nodes,
     $       selected_grdpts,
     $       match_table)

          implicit none

          integer       , dimension(:,:)  , intent(inout) :: grdpts_id
          real(rkind)   , dimension(:,:,:), intent(inout) :: nodes
          integer(ikind), dimension(:,:)  , intent(in)    :: selected_grdpts
          integer(ikind), dimension(2)    , intent(in)    :: match_table


          integer(ikind) :: k
          integer(ikind) :: i_prev, j_prev, i_center, j_center


          !we have a list of gridpoints that should be turned
          !from bc_interior_pt to interior_pt. For a point to be
          !considered an interior point, we need to make sure
          !that the grid points it needs to compute its time
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
          !possible in checking the neighbours so the centers
          !are initialized to do as if the previous central point
          !checked was far away
          !----------------------------------------------------
          i_center = -match_table(1)+nx/2
          j_center = -match_table(2)+ny/2

          do k=1, size(selected_grdpts,2)

             !update the position of the gridpoint previously
             !tested
             i_prev   =   i_center
             j_prev   =   j_center
             i_center = - match_table(1) + selected_grdpts(1,k)
             j_center = - match_table(2) + selected_grdpts(2,k)

             !check its neighbours
             call check_neighbors(grdpts_id,nodes,i_prev,j_prev,i_center,j_center)

             !update the status of the gridpoint
             grdpts_id(i_center,j_center) = interior_pt

          end do  

        end subroutine update_bc_interior_pt_to_interior_pt


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> check whether the neighboring points
        !> from a bc_interior_pt that should be turned into
        !> an interior point exists: if they exist, nothing
        !> happens. Otherwise, the grid point ID is updated
        !> and the grid point is computed
        !> the number of neighbors to be checked can be reduced
        !> knowing the previous neighbors checked
        !>
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !>
        !>@param grdpts_id
        !> table where the role of the grid points is stored
        !
        !>@param nodes
        !> table encapsulating the data of the grid points of the
        !> buffer layer
        !
        !>@param i_prev
        !> x-index of the previous bc_interior_pt whose 
        !> neighbors were checked
        !>
        !>@param j_prev
        !> y-index of the previous bc_interior_pt whose 
        !> neighbors were checked
        !>
        !>@param i_center
        !> x-index of the current bc_interior_pt whose 
        !> neighbors were checked
        !>
        !>@param j_center
        !> y-index of the current bc_interior_pt whose 
        !> neighbors were checked
        !---------------------------------------------------------------
        subroutine check_neighbors(
     $     grdpts_id, nodes,
     $     i_prev, j_prev,
     $     i_center, j_center)

          implicit none

          integer       , dimension(:,:)  , intent(inout) :: grdpts_id
          real(rkind)   , dimension(:,:,:), intent(inout) :: nodes
          integer(ikind)                  , intent(in)    :: i_prev
          integer(ikind)                  , intent(in)    :: j_prev
          integer(ikind)                  , intent(in)    :: i_center
          integer(ikind)                  , intent(in)    :: j_center
          

          integer(ikind) :: min_j, max_j
          integer(ikind) :: i,j


          min_j = min(j_center-j_prev,0)
          max_j = max(j_center-j_prev,0)

          !check neighbors for the extra boundary points bc_pt
          do j=j_center-bc_size, j_prev - bc_size + min(j_center-j_prev+2*bc_size,-1)
             do i=i_center-bc_size,i_center+bc_size
                call check_gridpoint(grdpts_id,nodes,i,j,i_center,j_center)
             end do
          end do

          do j=j_center-bc_size-min_j, j_center+bc_size-max_j
             do i=i_center-bc_size, i_prev-bc_size+min(i_center-i_prev+2*bc_size,-1)
                call check_gridpoint(grdpts_id,nodes,i,j,i_center,j_center)
             end do
          end do

          do j=j_center-bc_size-min_j, j_center+bc_size-max_j
             do i=i_prev+bc_size+max(i_center-i_prev-2*bc_size,1),i_center+bc_size
                call check_gridpoint(grdpts_id,nodes,i,j,i_center,j_center)
             end do
          end do

          do j=j_prev+bc_size+max(j_center-j_prev-2*bc_size,1), j_center+bc_size
             do i=i_center-bc_size,i_center+bc_size
                call check_gridpoint(grdpts_id,nodes,i,j,i_center,j_center)
             end do
          end do


          !update the bc_interior_pt overlapping the previous bc_pt
          do j=max(j_prev-bc_size, j_center-1), min(j_prev-bc_size, j_center+1)
             do i=max(i_prev-bc_size, i_center-1), min(i_prev+bc_size, i_center+1)
                if(grdpts_id(i,j).eq.bc_pt) then
                   grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(j_prev-1, j_center-1), min(j_prev+1, j_center+1)
             do i=max(i_prev-bc_size, i_center-1), min(i_prev-bc_size, i_center+1)
                if(grdpts_id(i,j).eq.bc_pt) then
                   grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(j_prev-1, j_center-1), min(j_prev+1, j_center+1)
             do i=max(i_prev+bc_size, i_center-1), min(i_prev+bc_size, i_center+1)
                if(grdpts_id(i,j).eq.bc_pt) then
                   grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(j_prev+bc_size, j_center-1), min(j_prev+bc_size, j_center+1)
             do i=max(i_center-1, i_prev-bc_size), min(i_prev+bc_size, i_center+1)
                if(grdpts_id(i,j).eq.bc_pt) then
                   grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

        end subroutine check_neighbors


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> check whether the grid point asked exists or not.
        !> If it does not exist in the buffer layer, it has
        !> to be computed and updated in the gridpoint ID map
        !>
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !>
        !>@param grdpts_id
        !> table where the role of the grid points is stored
        !
        !>@param nodes
        !> table encapsulating the data of the grid points of the
        !> buffer layer
        !
        !>@param i_center
        !> x-index of the central point
        !>
        !>@param j_center
        !> y-index of the central point
        !
        !>@param i
        !> x-index of the grid points checked
        !>
        !>@param j
        !> y-index of the grid points checked
        !---------------------------------------------------------------
        subroutine check_gridpoint(
     $     grdpts_id,nodes,i,j,i_center,j_center)

          implicit none

          integer       , dimension(:,:)  , intent(inout) :: grdpts_id
          real(rkind)   , dimension(:,:,:), intent(out)   :: nodes
          integer(ikind)                  , intent(in)    :: i
          integer(ikind)                  , intent(in)    :: j
          integer(ikind)                  , intent(in)    :: i_center
          integer(ikind)                  , intent(in)    :: j_center


          if((grdpts_id(i,j).eq.no_pt).or.(grdpts_id(i,j).eq.bc_pt)) then

             if (grdpts_id(i,j).eq.no_pt) then
                call compute_new_grdpt(nodes,i,j)
             end if

             !if the grid point is next to the new interior point,
             !it is a bc_interior_pt, otherwise, it is a bc_pt
             if((abs(i_center-i).le.1).and.(abs(j_center-j).le.1)) then
                grdpts_id(i,j) = bc_interior_pt
             else
                grdpts_id(i,j) = bc_pt
             end if

          end if

        end subroutine check_gridpoint


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the can_remain attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param remain_state
        !> logical identifying whether the buffer layer should be
        !> removed or not
        !--------------------------------------------------------------
        subroutine set_remain_status(this, remain_state)

          implicit none

          class(bf_layer), intent(inout) :: this
          logical        , intent(in)    :: remain_state

          this%can_remain = remain_state
          
        end subroutine set_remain_status


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the can_remain attribute
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@return get_remain_status
        !> logical identifying whether the buffer layer should be
        !> removed or not
        !--------------------------------------------------------------
        function get_remain_status(this)

          implicit none

          class(bf_layer), intent(in) :: this
          logical                     :: get_remain_status

          get_remain_status = this%can_remain

        end function get_remain_status 


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check the grid points at the edge between the buffer layer
        !> and the interior domain and compute whether they undermine
        !> the open boundary conditions. This determines the buffer layer
        !> should be removed or not
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param interior_nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@return get_remain_status
        !> logical identifying whether the buffer layer should be
        !> removed or not
        !--------------------------------------------------------------
        function should_remain(this, interior_nodes, p_model)

          implicit none

          class(bf_layer)                 , intent(in) :: this
          real(rkind), dimension(nx,ny,ne), intent(in) :: interior_nodes
          type(pmodel_eq)                 , intent(in) :: p_model
          logical                                      :: should_remain
          
          integer(ikind), dimension(2) :: bf_match_table

          bf_match_table = get_general_to_local_coord_tab(this)

          should_remain = check_if_bf_layer_remains(
     $         this%localization,
     $         this%alignment,
     $         bf_match_table,
     $         this%grdpts_id,
     $         this%nodes,
     $         interior_nodes,
     $         p_model)

        end function should_remain


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the buffe rlayer by deallocating the main tables
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

          class(bf_layer), intent(inout) :: this

          deallocate(this%nodes)
          deallocate(this%grdpts_id)

        end subroutine remove      


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the nodes and the grdpts_id attributes
        !> as well as the size of the previous tables in
        !> output binary files
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param filename_nodes
        !> name of the output file for the data of the nodes attribute
        !
        !>@param filename_grdpts_id
        !> name of the output file for the data of the grdpts_id attribute
        !
        !>@param filename_sizes
        !> name of the output file for the sizes of the nodes and
        !> grpdts_id attributes
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the nodes and the grdpts_id attributes
        !> on a netcdf file
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param filename
        !> name of the output file for the grdpts_id of the nodes
        !> attributes
        !
        !>@param name_var
        !> physical model used to know the name of the governing
        !> variables when writing the netcdf file
        !
        !>@param bf_order
        !> identification of the buffer layer in the main layer to
        !> determine the title of the netcdf file
        !
        !>@param x_min_interior
        !> x-coordinate corresponding to the grid point next to the
        !> left boundary layer in the interior domain
        !
        !>@param y_min_interior
        !> y-coordinate corresponding to the grid point next to the
        !> lower boundary layer in the interior domain
        !
        !>@param dx
        !> buffer layer grid size along the x-direction
        !
        !>@param dy
        !> buffer layer grid size along the y-direction
        !
        !>@param time
        !> time corresponding to the buffer layer data
        !--------------------------------------------------------------
        subroutine print_netcdf(
     $     this,
     $     filename,
     $     name_var,
     $     longname_var,
     $     unit_var,
     $     bf_order,
     $     x_min_interior,
     $     y_min_interior,
     $     dx,dy,
     $     time)

          implicit none

          class(bf_layer)            , intent(inout) :: this
          character(*)               , intent(in)    :: filename
          character(*), dimension(ne), intent(in)    :: name_var
          character(*), dimension(ne), intent(in)    :: longname_var
          character(*), dimension(ne), intent(in)    :: unit_var
          integer                    , intent(in)    :: bf_order
          real(rkind)                , intent(in)    :: x_min_interior
          real(rkind)                , intent(in)    :: y_min_interior
          real(rkind)                , intent(in)    :: dx
          real(rkind)                , intent(in)    :: dy
          real(rkind)                , intent(in)    :: time

          real(rkind) :: x_start
          real(rkind) :: y_start

          x_start = (this%alignment(1,1)-2*bc_size-1)*dx +
     $              x_min_interior

          y_start = (this%alignment(2,1)-2*bc_size-1)*dy +
     $              y_min_interior

          call print_bf_layer_on_netcdf(
     $         filename,
     $         name_var, longname_var, unit_var,
     $         this%localization, bf_order,
     $         x_start, y_start, dx, dy,
     $         this%grdpts_id, this%nodes, time)

        end subroutine print_netcdf


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the grid size of the computations
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine ini_for_comput(this,dx,dy)

          implicit none

          class(bf_layer), intent(inout) :: this
          real(rkind)    , intent(in)    :: dx
          real(rkind)    , intent(in)    :: dy

          call this%bf_compute_used%ini(dx,dy)

        end subroutine ini_for_comput


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine allocate_before_timeInt(this)

          implicit none

          class(bf_layer), intent(inout) :: this

          call this%bf_compute_used%allocate_tables(
     $         size(this%nodes,1),
     $         size(this%nodes,2))

        end subroutine allocate_before_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate the memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine deallocate_after_timeInt(this)

          implicit none

          class(bf_layer), intent(inout) :: this

          call this%bf_compute_used%deallocate_tables()

        end subroutine deallocate_after_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine compute_time_dev(this)

          implicit none

          class(bf_layer), intent(inout) :: this

          call this%bf_compute_used%compute_time_dev(
     $         this%nodes,
     $         this%grdpts_id)

        end subroutine compute_time_dev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the integration step
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param dt
        !> integration time step
        !
        !>@param integration_step_nopt
        !> procedure performing the time integration
        !--------------------------------------------------------------
        subroutine compute_integration_step(
     $     this, dt, integration_step_nopt)

          implicit none

          class(bf_layer)              , intent(inout) :: this
          real(rkind)                  , intent(in)    :: dt
          procedure(timeInt_step_nopt) :: integration_step_nopt

          call this%bf_compute_used%compute_integration_step(
     $         this%grdpts_id,
     $         this%nodes,
     $         dt,
     $         integration_step_nopt)

        end subroutine compute_integration_step


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the time derivatives
        !
        !> @date
        !> 16_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine get_time_dev(this, time_dev)

          implicit none

          class(bf_layer), intent(in) :: this
          real(rkind), dimension(:,:,:), allocatable, intent(out)::time_dev

          call this%bf_compute_used%get_time_dev(time_dev)

        end subroutine get_time_dev

      end module bf_layer_class

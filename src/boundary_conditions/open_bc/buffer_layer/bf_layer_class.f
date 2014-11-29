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

        use bc_operators_class, only :
     $       bc_operators

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
     $       copy_from_bf1_to_bf2,
     $       get_sync_indices_with_interior,
     $       sync_nodes,
     $       get_sync_indices_with_neighbor1,
     $       get_sync_indices_with_neighbor2

        use bf_layer_nf90_operators_module, only :
     $       print_bf_layer_on_netcdf

        use bf_remove_module, only :
     $       check_if_bf_layer_remains

        use bf_suspicious_bc_interior_pt_module, only :
     $       verify_if_all_grdpts_exist

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
     $       no_pt,
     $       BF_SUCCESS

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

        use sd_operators_class, only :
     $       sd_operators

        use td_operators_class, only :
     $       td_operators


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
        !> @param x_borders
        !> interval identifying the extent of the integration domain
        !> in the x-direction
        !
        !> @param y_borders
        !> interval identifying the extent of the integration domain
        !> in the y-direction
        !
        !> @param N_bc_sections
        !> North boundary section identifying the grid points that
        !> should be included when defining the bc_procedures for
        !> the buffer layer
        !
        !> @param S_bc_sections
        !> South boundary section identifying the grid points that
        !> should be included when defining the bc_procedures for
        !> the buffer layer
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
        !> @param sync_nodes_with_interior
        !> synchronize the nodes at the interface with the interior
        !
        !> @param sync_nodes_with_neighbor1
        !> synchronize the nodes at the interface with the neighbor1
        !
        !> @param sync_nodes_with_neighbor2
        !> synchronize the nodes at the interface with the neighbor2
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
        !
        !> @param set_x_borders
        !> set the index borders for time integration in the
        !> x-direction
        !
        !> @param set_y_borders
        !> set the index borders for time integration in the
        !> y-direction
        !
        !> @param get_x_borders
        !> get the index borders for time integration in the
        !> x-direction
        !
        !> @param get_y_borders
        !> get the index borders for time integration in the
        !> y-direction
        !
        !> @param set_N_bc_sections
        !> set the borders of the north boundary layer computed
        !> by time integration
        !
        !> @param set_S_bc_sections
        !> set the borders of the south boundary layer computed
        !> by time integration
        !
        !> @param get_N_bc_sections
        !> get the borders of the north boundary layer computed
        !> by time integration
        !
        !> @param get_S_bc_sections
        !> get the borders of the south boundary layer computed
        !> by time integration
        !
        !> @param remove_N_bc_sections
        !> remove the borders of the north boundary layer computed
        !> by time integration
        !
        !> @param remove_S_bc_sections
        !> remove the borders of the south boundary layer computed
        !> by time integration
        !-------------------------------------------------------------
        type :: bf_layer

          integer                       , private :: localization
          integer(ikind), dimension(2,2), private :: alignment

          real(rkind), dimension(:)    , allocatable, private :: x_map
          real(rkind), dimension(:)    , allocatable, private :: y_map
          real(rkind), dimension(:,:,:), allocatable, private :: nodes
          integer    , dimension(:,:)  , allocatable, private :: grdpts_id

          logical, private :: shares_grdpts_with_neighbor1
          logical, private :: shares_grdpts_with_neighbor2

          logical, private :: can_remain

          type(bf_compute)                            :: bf_compute_used
          integer(ikind), dimension(2)                :: x_borders
          integer(ikind), dimension(2)                :: y_borders
          integer(ikind), dimension(:,:), allocatable :: N_bc_sections
          integer(ikind), dimension(:,:), allocatable :: S_bc_sections

          contains

          procedure,   pass :: ini
                       
          !procedures to modify the main attributes:
          !localization, alignment, grdpts_id, nodes
          procedure,   pass :: get_localization
          procedure,   pass :: set_localization
          procedure,   pass :: get_sizes
          procedure,   pass :: set_nodes
          procedure,   pass :: set_nodes_pt
          procedure,   pass :: set_grdpts_id
          procedure,   pass :: set_grdpts_id_pt
          procedure,   pass :: get_alignment
          procedure,   pass :: get_alignment_tab
          procedure,   pass :: set_alignment_tab
                       
          !procedure to localize the buffer layer:
          !coordinates and coordinate-maps
          procedure,   pass :: get_local_coord
          procedure,   pass :: get_general_to_local_coord_tab
          procedure,   pass :: get_x_map
          procedure,   pass :: get_y_map
          procedure,   pass :: set_x_map
          procedure,   pass :: set_y_map
          procedure,   pass :: get_nodes_array
          procedure,   pass :: get_nodes
          procedure,   pass :: get_grdpts_id
          
          !procedure to modify the structure of the buffer
          !layer
          procedure,   pass :: allocate_bf_layer
          procedure,   pass :: reallocate_bf_layer
          procedure,   pass :: merge_bf_layer

          !procedures to set the exchange with the neighboring
          !buffer layers and the interior domain
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

          procedure,   pass :: sync_nodes_with_interior
          procedure,   pass :: sync_nodes_with_neighbor1
          procedure,   pass :: sync_nodes_with_neighbor2

          !procedures for the computation of new grid points
          procedure,   pass :: get_data_for_newgrdpt
          procedure,   pass :: compute_newgrdpt
          procedure,   pass :: does_previous_timestep_exist

          procedure,   pass :: set_new_interior_pt
          procedure,   pass :: finalize_neighboring_grdpts_update
          procedure,   pass :: update_neighboring_grdpt
          procedure,   pass :: get_grdpts_id_part
          procedure,   pass :: verify_if_all_bf_grdpts_exist
          procedure,   pass :: has_a_bc_pt_neighbor
          procedure,   pass :: is_bc_interior_pt

          !procedures for removing the buffer layer
          procedure,   pass :: set_remain_status
          procedure,   pass :: get_remain_status
          procedure,   pass :: should_remain
          procedure,   pass :: remove
          
          !i/o procedures
          procedure,   pass :: print_binary
          procedure,   pass :: print_netcdf

          
          !for time integration: interior + boundaries
          procedure,   pass :: allocate_before_timeInt
          procedure,   pass :: deallocate_after_timeInt
          procedure,   pass :: compute_time_dev
          procedure,   pass :: compute_integration_step
          procedure,   pass :: get_time_dev  !only for tests

          procedure,   pass :: set_x_borders
          procedure,   pass :: set_y_borders
          procedure,   pass :: get_x_borders
          procedure,   pass :: get_y_borders

          procedure,   pass :: set_N_bc_sections
          procedure,   pass :: set_S_bc_sections
          procedure,   pass :: remove_N_bc_sections
          procedure,   pass :: remove_S_bc_sections
          procedure,   pass :: get_N_bc_sections
          procedure,   pass :: get_S_bc_sections

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
        subroutine ini(this,localization)

          implicit none

          class(bf_layer), intent(inout) :: this
          integer(ikind) , intent(in)    :: localization
          
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

          class(bf_layer), intent(in) :: this
          integer                     :: get_localization
          
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

          class(bf_layer), intent(inout) :: this
          integer                        :: localization
          
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

          class(bf_layer)               , intent(inout) :: this
          integer(ikind), dimension(2,2), intent(in)    :: alignment

          this%alignment = alignment

        end subroutine set_alignment_tab


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

          class(bf_layer)                          , intent(in) :: this
          real(rkind)   , dimension(:), allocatable, intent(out):: x_map
          
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

          class(bf_layer)                          , intent(in) :: this
          real(rkind)   , dimension(:), allocatable, intent(out):: y_map
          
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

          class(bf_layer)                          , intent(inout) :: this
          real(rkind)   , dimension(:), allocatable, intent(inout) :: x_map
          
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

          class(bf_layer)                          , intent(inout) :: this
          real(rkind)   , dimension(:), allocatable, intent(inout) :: y_map
          
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

          class(bf_layer)                           , intent(in)    :: this
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
        subroutine allocate_bf_layer(this, x_map, y_map, nodes, alignment)

          implicit none
          
          class(bf_layer)                 , intent(inout) :: this
          real(rkind)   , dimension(:)    , intent(in)    :: x_map
          real(rkind)   , dimension(:)    , intent(in)    :: y_map
          real(rkind)   , dimension(:,:,:), intent(in)    :: nodes
          integer(ikind), dimension(2,2)  , intent(in)    :: alignment


          select case(this%localization)
            case(N)
               call allocate_bf_layer_N(
     $              this%x_map, x_map,
     $              this%y_map, y_map,
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)
            case(S)
               call allocate_bf_layer_S(
     $              this%x_map, x_map,
     $              this%y_map, y_map,
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)
            case(E)
               call allocate_bf_layer_E(
     $              this%x_map, x_map,
     $              this%y_map, y_map,
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)
            case(W)
               call allocate_bf_layer_W(
     $              this%x_map, x_map,
     $              this%y_map, y_map,
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)
            case default
               call error_mainlayer_id(
     $              'bf_layer_class.f',
     $              'allocate_bf_layer',
     $              this%localization)
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
        subroutine reallocate_bf_layer(this, x_map, y_map, nodes, alignment)

          implicit none

          class(bf_layer)                 , intent(inout) :: this
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer    , dimension(2,2)     , intent(in)    :: alignment

          select case(this%localization)

            case(N)
               call reallocate_bf_layer_N(
     $              this%x_map, x_map,
     $              this%y_map, y_map,
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(S)
               call reallocate_bf_layer_S(
     $              this%x_map, x_map,
     $              this%y_map, y_map,
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(E)
               call reallocate_bf_layer_E(
     $              this%x_map, x_map,
     $              this%y_map, y_map,
     $              this%nodes, nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(W)
               call reallocate_bf_layer_W(
     $              this%x_map, x_map,
     $              this%y_map, y_map,
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
        subroutine merge_bf_layer(this, bf_layer2, x_map, y_map, nodes, alignment)

          implicit none

          class(bf_layer)                               , intent(inout) :: this
          class(bf_layer)                               , intent(inout) :: bf_layer2
          real(rkind)    , dimension(nx)                , intent(in)    :: x_map
          real(rkind)    , dimension(ny)                , intent(in)    :: y_map
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
     $                 this%x_map    , bf_layer2%x_map, x_map,
     $                 this%y_map    , bf_layer2%y_map, y_map,
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment, alignment)
               else
                  call merge_bf_layers_N(
     $                 this%x_map    , bf_layer2%x_map, x_map,
     $                 this%y_map    , bf_layer2%y_map, y_map,
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment)
               end if

            case(S)
               if(present(alignment)) then
                  call merge_bf_layers_S(
     $                 this%x_map    , bf_layer2%x_map, x_map,
     $                 this%y_map    , bf_layer2%y_map, y_map,
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment, alignment)
               else
                  call merge_bf_layers_S(
     $                 this%x_map    , bf_layer2%x_map, x_map,
     $                 this%y_map    , bf_layer2%y_map, y_map,
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment)
               end if

            case(E)
               if(present(alignment)) then
                  call merge_bf_layers_E(
     $                 this%x_map    , bf_layer2%x_map, x_map,
     $                 this%y_map    , bf_layer2%y_map, y_map,
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment, alignment)
               else
                  call merge_bf_layers_E(
     $                 this%x_map    , bf_layer2%x_map, x_map,
     $                 this%y_map    , bf_layer2%y_map, y_map,
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment)
               end if

            case(W)
               if(present(alignment)) then
                  call merge_bf_layers_W(
     $                 this%x_map    , bf_layer2%x_map, x_map,
     $                 this%y_map    , bf_layer2%y_map, y_map,
     $                 this%nodes    , bf_layer2%nodes, nodes,
     $                 this%grdpts_id, bf_layer2%grdpts_id,
     $                 this%alignment, bf_layer2%alignment, alignment)
               else
                  call merge_bf_layers_W(
     $                 this%x_map    , bf_layer2%x_map, x_map,
     $                 this%y_map    , bf_layer2%y_map, y_map,
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
                  this%shares_grdpts_with_neighbor1 = this%alignment(1,1).le.(align_W+bc_size+1)
               case(E,W)
                  this%shares_grdpts_with_neighbor1 = this%alignment(2,1).le.(align_S+bc_size+1)
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
                  this%shares_grdpts_with_neighbor2 = this%alignment(1,2).ge.(align_E-bc_size-1)
               case(E,W)
                  this%shares_grdpts_with_neighbor2 = this%alignment(2,2).ge.(align_N-bc_size-1)
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
        !> update the common layers between the buffer
        !> layer and the interior domain
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param interior_nodes
        !> grid points for the interior domain
        !--------------------------------------------------------------
        subroutine sync_nodes_with_interior(this, interior_nodes)

          implicit none

          class(bf_layer)                 , intent(inout) :: this
          real(rkind), dimension(nx,ny,ne), intent(inout) :: interior_nodes

          integer(ikind), dimension(2) :: in_send
          integer(ikind), dimension(2) :: in_recv
          integer(ikind), dimension(2) :: bf_send
          integer(ikind), dimension(2) :: bf_recv
          integer(ikind), dimension(2) :: ex_size

          !get the indices identifying which arrays are exchanged
          call get_sync_indices_with_interior(
     $         this%localization,
     $         this%alignment,
     $         [size(this%nodes,1),size(this%nodes,2)],
     $         in_send,
     $         in_recv,
     $         bf_send,
     $         bf_recv,
     $         ex_size)

          !exchange the arrays between the buffer layer
          !and the interior domain
          call sync_nodes(
     $         interior_nodes,
     $         in_send,
     $         in_recv,
     $         this%nodes,
     $         bf_send,
     $         bf_recv,
     $         ex_size)

        end subroutine sync_nodes_with_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the common layers between the buffer
        !> layer and another buffer layer which must be
        !> of type neighbor1 for the first buffer layer
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param neighbor1
        !> buffer layer which is of neighbor1 type
        !--------------------------------------------------------------
        subroutine sync_nodes_with_neighbor1(this, neighbor1)

          implicit none

          class(bf_layer), intent(inout) :: this
          class(bf_layer), intent(inout) :: neighbor1

          integer(ikind), dimension(2) :: bf_send
          integer(ikind), dimension(2) :: bf_recv
          integer(ikind), dimension(2) :: nbf_send
          integer(ikind), dimension(2) :: nbf_recv
          integer(ikind), dimension(2) :: ex_size

          !get the indices identifying which arrays are exchanged
          call get_sync_indices_with_neighbor1(
     $         this%localization,
     $         this%alignment,
     $         [size(this%nodes,1),size(this%nodes,2)],
     $         bf_send,
     $         bf_recv,
     $         neighbor1%alignment,
     $         [size(neighbor1%nodes,1),size(neighbor1%nodes,2)],
     $         nbf_send,
     $         nbf_recv,
     $         ex_size)

          !exchange the arrays between the buffer layer
          !and the interior domain
          call sync_nodes(
     $         this%nodes,
     $         bf_send,
     $         bf_recv,
     $         neighbor1%nodes,
     $         nbf_send,
     $         nbf_recv,
     $         ex_size)

        end subroutine sync_nodes_with_neighbor1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the common layers between the buffer
        !> layer and another buffer layer which must be
        !> of type neighbor2 for the first buffer layer
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param neighbor1
        !> buffer layer which is of neighbor1 type
        !--------------------------------------------------------------
        subroutine sync_nodes_with_neighbor2(this, neighbor2)

          implicit none

          class(bf_layer), intent(inout) :: this
          class(bf_layer), intent(inout) :: neighbor2

          integer(ikind), dimension(2) :: bf_send
          integer(ikind), dimension(2) :: bf_recv
          integer(ikind), dimension(2) :: nbf_send
          integer(ikind), dimension(2) :: nbf_recv
          integer(ikind), dimension(2) :: ex_size

          !get the indices identifying which arrays are exchanged
          call get_sync_indices_with_neighbor2(
     $         this%localization,
     $         this%alignment,
     $         [size(this%nodes,1),size(this%nodes,2)],
     $         bf_send,
     $         bf_recv,
     $         neighbor2%alignment,
     $         [size(neighbor2%nodes,1),size(neighbor2%nodes,2)],
     $         nbf_send,
     $         nbf_recv,
     $         ex_size)

          !exchange the arrays between the buffer layer
          !and the interior domain
          call sync_nodes(
     $         this%nodes,
     $         bf_send,
     $         bf_recv,
     $         neighbor2%nodes,
     $         nbf_send,
     $         nbf_recv,
     $         ex_size)

        end subroutine sync_nodes_with_neighbor2


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the grdpts_id, the coordinate maps and the
        !> nodes at t-dt and t corresponding to the general
        !> coordinates gen_coords
        !    ___________________
        !   |                  _|_________
        !   |    buffer layer |/|         |
        !   |                 |/|  tmp    |
        !   !                 !/!         !
        !                   overlapping which is copied
        !                     from buffer layer to tmp
        !
        !> @date
        !> 18_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param tmp_grdpts_id0
        !> array with the grdpts_id data
        !
        !>@param tmp_nodes0
        !> array with the grid points data at t-dt
        !
        !>@param tmp_nodes1
        !> array with the grid points data at t
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> tmp arrays computed
        !--------------------------------------------------------------
        subroutine get_data_for_newgrdpt(
     $     this,
     $     tmp_grdpts_id0,
     $     tmp_nodes0,
     $     tmp_nodes1,
     $     gen_coords)

          implicit none

          class(bf_layer)                                       , intent(in)    :: this
          integer        , dimension(2*bc_size+1,2*bc_size+1)   , intent(inout) :: tmp_grdpts_id0
          real(rkind)    , dimension(2*bc_size+1,2*bc_size+1,ne), intent(inout) :: tmp_nodes0
          real(rkind)    , dimension(2*bc_size+1,2*bc_size+1,ne), intent(inout) :: tmp_nodes1
          integer(ikind) , dimension(2,2)                       , intent(in)    :: gen_coords


          integer(ikind) :: size_x,size_y
          integer(ikind) :: i_recv,i_send,j_recv,j_send
          integer(ikind) :: i,j
          integer        :: k


          !synchronize the overlapping at t=t
          call this%bf_compute_used%get_sync_indices_for_newgrdpt_data(
     $         this%alignment,
     $         gen_coords,
     $         size_x, size_y,
     $         i_recv, j_recv,
     $         i_send, j_send)        


          do k=1,ne
             do j=1, size_y
                do i=1, size_x

                   tmp_nodes1(i_recv+i-1,j_recv+j-1,k) =
     $                  this%nodes(i_send+i-1,j_send+j-1,k)

                end do
             end do
          end do


          !synchronize the overlapping at t-dt
          call this%bf_compute_used%get_data_for_newgrdpt(
     $         tmp_grdpts_id0,
     $         tmp_nodes0,
     $         gen_coords)

        end subroutine get_data_for_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid point
        !
        !> @date
        !> 18_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i1
        !> index identifying the x-coordinate of the new grid point
        !
        !>@param j1
        !> index identifying the y-coordinate of the new grid point
        !
        !>@param t
        !> time
        !
        !>@param dt
        !> time step
        !--------------------------------------------------------------
        subroutine compute_newgrdpt(this,p_model,t,dt,i1,j1)

          implicit none

          class(bf_layer), intent(inout) :: this
          type(pmodel_eq), intent(in)    :: p_model
          real(rkind)    , intent(in)    :: t
          real(rkind)    , intent(in)    :: dt
          integer(ikind) , intent(in)    :: i1
          integer(ikind) , intent(in)    :: j1
          
          this%nodes(i1,j1,:) = this%bf_compute_used%compute_newgrdpt(
     $         p_model, t, dt,
     $         this%alignment,
     $         this%x_map,
     $         this%y_map,
     $         this%nodes,
     $         i1,j1)

        end subroutine compute_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the previous time step is stored in the
        !> buffer layer
        !
        !> @date
        !> 20_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param exist
        !> says whether the previous time step is stored in the
        !> buffer layer
        !--------------------------------------------------------------
        function does_previous_timestep_exist(this) result(exist)

          implicit none

          class(bf_layer), intent(in) :: this
          logical                     :: exist

          exist = this%bf_compute_used%does_previous_timestep_exist()

        end function does_previous_timestep_exist


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
        !>
        !> @brief
        !> set the grid point as interior_pt
        !>
        !> @date
        !> 21_11_2014 - initial version - J.L. Desmarais
        !>
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i
        !> integer identifying the x-index of the
        !> updated grid point
        !
        !>@param j
        !> integer identifying the y-index of the
        !> updated grid point
        !--------------------------------------------------------------
        subroutine set_new_interior_pt(this,i,j)

          implicit none

          class(bf_layer), intent(inout) :: this
          integer(ikind) , intent(in)    :: i
          integer(ikind) , intent(in)    :: j

          this%grdpts_id(i,j) = interior_pt

        end subroutine set_new_interior_pt


        !> @author
        !> Julien L. Desmarais
        !>
        !> @brief
        !> as the grid point is set as new interior_pt, the neighboring
        !> gridpoints should be updated: here only the identity of the
        !> neighboring grid points is updated (bc_pt-> bc_interior_pt)
        !
        !> @date
        !> 21_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i_prev
        !> integer identifying the x-index of the
        !> previous central grid point that was
        !> updated to interior_pt
        !
        !>@param j_prev
        !> integer identifying the y-index of the
        !> previous central grid point that was
        !> updated to interior_pt
        !
        !>@param i_center
        !> integer identifying the x-index of the
        !> central grid point which is updated to
        !> interior_pt
        !
        !>@param j_center
        !> integer identifying the y-index of the
        !> central grid point which is updated to
        !> interior_pt
        !--------------------------------------------------------------
        subroutine finalize_neighboring_grdpts_update(
     $     this,
     $     i_prev, j_prev,
     $     i_center, j_center)

          implicit none

          class(bf_layer), intent(inout) :: this
          integer(ikind) , intent(in)    :: i_prev
          integer(ikind) , intent(in)    :: j_prev
          integer(ikind) , intent(in)    :: i_center
          integer(ikind) , intent(in)    :: j_center
          
          integer(ikind) :: i,j

          do j=max(j_prev-bc_size, j_center-1), min(j_prev-bc_size, j_center+1)
             do i=max(i_prev-bc_size, i_center-1), min(i_prev+bc_size, i_center+1)
                if(this%grdpts_id(i,j).eq.bc_pt) then
                   this%grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(j_prev-1, j_center-1), min(j_prev+1, j_center+1)
             do i=max(i_prev-bc_size, i_center-1), min(i_prev-bc_size, i_center+1)
                if(this%grdpts_id(i,j).eq.bc_pt) then
                   this%grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(j_prev-1, j_center-1), min(j_prev+1, j_center+1)
             do i=max(i_prev+bc_size, i_center-1), min(i_prev+bc_size, i_center+1)
                if(this%grdpts_id(i,j).eq.bc_pt) then
                   this%grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

          do j=max(j_prev+bc_size, j_center-1), min(j_prev+bc_size, j_center+1)
             do i=max(i_center-1, i_prev-bc_size), min(i_prev+bc_size, i_center+1)
                if(this%grdpts_id(i,j).eq.bc_pt) then
                   this%grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

        end subroutine finalize_neighboring_grdpts_update


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the grid point asked exists or not.
        !> If it does not exist in the buffer layer, it has
        !> to be computed and updated in the gridpoint ID map
        !
        !> @date
        !> 21_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param i
        !> x-index of the grid points checked
        !
        !>@param j
        !> y-index of the grid points checked
        !
        !>@param i_center
        !> x-index of the central point
        !
        !>@param j_center
        !> y-index of the central point 
        !
        !> @return compute_newgrdpt
        !> logical identifying whether the grid point should be
        !> computed
        !---------------------------------------------------------------
        function update_neighboring_grdpt(this,i,j,i_center,j_center)
     $     result(compute_newgrdpt)

          implicit none

          class(bf_layer), intent(inout) :: this
          integer(ikind) , intent(in)    :: i
          integer(ikind) , intent(in)    :: j
          integer(ikind) , intent(in)    :: i_center
          integer(ikind) , intent(in)    :: j_center
          logical                        :: compute_newgrdpt

          
          compute_newgrdpt = .false.

          if(
     $         (this%grdpts_id(i,j).eq.no_pt).or.
     $         (this%grdpts_id(i,j).eq.bc_pt)) then

             if (this%grdpts_id(i,j).eq.no_pt) then
                compute_newgrdpt = .true.
             end if

             !if the grid point is next to the new interior point,
             !it is a bc_interior_pt, otherwise, it is a bc_pt
             if((abs(i_center-i).le.1).and.(abs(j_center-j).le.1)) then
                this%grdpts_id(i,j) = bc_interior_pt
             else
                this%grdpts_id(i,j) = bc_pt
             end if

          end if        

        end function update_neighboring_grdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the grdpts_id needed to evaluate whether the
        !> bc_interior_pt investigated should be turned into
        !> an interior_pt
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param 
        !> logical identifying whether the buffer layer should be
        !> removed or not
        !--------------------------------------------------------------
        subroutine get_grdpts_id_part(
     $     this,
     $     tmp_grdptsid,
     $     gen_coords)

          implicit none

          class(bf_layer)                            , intent(in)    :: this
          integer, dimension(2*bc_size+1,2*bc_size+1), intent(inout) :: tmp_grdptsid
          integer(ikind), dimension(2,2)             , intent(in)    :: gen_coords

          
          integer(ikind) :: size_x,size_y
          integer(ikind) :: i_recv,i_send,j_recv,j_send
          integer(ikind) :: i,j


          !get the synchronization indices
          call this%bf_compute_used%get_sync_indices_for_newgrdpt_data(
     $         this%alignment,
     $         gen_coords,
     $         size_x, size_y,
     $         i_recv, j_recv,
     $         i_send, j_send)


          !fill the grid points asked
          do j=1, size_y
             do i=1, size_x
                
                tmp_grdptsid(i_recv+i-1,j_recv+j-1) =
     $               this%grdpts_id(i_send+i-1,j_send+j-1)

             end do
          end do

        end subroutine get_grdpts_id_part


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> verify if all the grid points around a central grid point
        !> exist in order to be turned into an interior_pt
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param 
        !> logical identifying whether all grid points around a
        !> central grid point exist in order to be turned into
        !> an interior_pt
        !--------------------------------------------------------------
        function verify_if_all_bf_grdpts_exist(
     $     this,i_c,j_c,ierror)
     $     result(all_grdpts_exist)

          implicit none

          class(bf_layer), intent(in)  :: this
          integer(ikind) , intent(in)  :: i_c
          integer(ikind) , intent(in)  :: j_c
          logical        , intent(out) :: ierror
          logical                      :: all_grdpts_exist

          if(
     $         (i_c.ge.(bc_size+1)).and.
     $         (i_c.le.size(this%grdpts_id,1)-bc_size).and.
     $         (j_c.ge.(bc_size+1)).and.
     $         (j_c.le.size(this%grdpts_id,2)-bc_size)) then

             ierror = BF_SUCCESS

             all_grdpts_exist = verify_if_all_grdpts_exist(
     $            i_c, j_c, this%grdpts_id)
             
          else

             ierror = .not.BF_SUCCESS
             all_grdpts_exist = .false.

          end if

        end function verify_if_all_bf_grdpts_exist


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if there is a bc_pt in the neighboring grid
        !> points
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param local_coords 
        !> local coordinates of the grid point whose neighbors are
        !> tested
        !
        !> @return bc_pt_neighbor
        !> logical stating whether there is a neighboring grid point 
        !> which is a bc_pt
        !--------------------------------------------------------------
        function has_a_bc_pt_neighbor(
     $     this, local_coords, ierror)
     $     result(bc_pt_neighbor)

          implicit none

          class(bf_layer)              , intent(in) :: this
          integer(ikind) , dimension(2), intent(in) :: local_coords
          logical                      , intent(out):: ierror
          logical                                   :: bc_pt_neighbor

          integer :: i,j


          if(
     $         (local_coords(1).gt.1).and.
     $         (local_coords(1).lt.size(this%grdpts_id,1)).and.
     $         (local_coords(2).gt.1).and.
     $         (local_coords(2).lt.size(this%grdpts_id,2))) then

             ierror = BF_SUCCESS
             bc_pt_neighbor = .false.
             
             do j=local_coords(2)-1,local_coords(2)+1
                do i=local_coords(1)-1,local_coords(1)+1
                   
                   if(this%grdpts_id(i,j).eq.bc_pt) then
                      bc_pt_neighbor=.true.
                      exit
                   end if
                   
                end do
             end do

          else

             ierror = .not.BF_SUCCESS
             bc_pt_neighbor = .false.

          end if

        end function has_a_bc_pt_neighbor


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if there is a bc_pt in the neighboring grid
        !> points
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param local_coords 
        !> local coordinates of the grid point whose neighbors are
        !> tested
        !
        !> @return bc_pt_neighbor
        !> logical stating whether there is a neighboring grid point 
        !> which is a bc_pt
        !--------------------------------------------------------------
        function is_bc_interior_pt(
     $     this, i,j)

          implicit none

          class(bf_layer), intent(in) :: this
          integer(ikind) , intent(in) :: i
          integer(ikind) , intent(in) :: j
          logical                     :: is_bc_interior_pt

          is_bc_interior_pt = this%grdpts_id(i,j).eq.bc_interior_pt

        end function is_bc_interior_pt


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

          deallocate(this%x_map)
          deallocate(this%y_map)
          deallocate(this%nodes)
          deallocate(this%grdpts_id)

          call remove_N_bc_sections(this)
          call remove_S_bc_sections(this)

          if(does_previous_timestep_exist(this)) then
             call this%bf_compute_used%deallocate_tables()
          end if

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
     $     filename_x_map,
     $     filename_y_map,
     $     filename_nodes,
     $     filename_grdpts_id,
     $     filename_sizes,
     $     timedev)

          implicit none

          class(bf_layer)  , intent(in) :: this
          character(*)     , intent(in) :: filename_x_map
          character(*)     , intent(in) :: filename_y_map
          character(*)     , intent(in) :: filename_nodes
          character(*)     , intent(in) :: filename_grdpts_id
          character(*)     , intent(in) :: filename_sizes
          logical, optional, intent(in) :: timedev

          integer :: ios
          logical :: timedev_op

          real(rkind), dimension(:,:,:), allocatable :: time_dev

          if(present(timedev)) then
             timedev_op = timedev
          else
             timedev_op = .false.
          end if
          
          !x_map
          open(unit=2,
     $         file=filename_x_map,
     $         action="write", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)

          if(ios.eq.0) then
             write(unit=2, iostat=ios) this%x_map
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          !y_map
          open(unit=2,
     $         file=filename_y_map,
     $         action="write", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)

          if(ios.eq.0) then
             write(unit=2, iostat=ios) this%y_map
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          !nodes
          open(unit=3,
     $          file=filename_nodes,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

          if(timedev_op) then
             call this%bf_compute_used%get_time_dev(time_dev)

             if(ios.eq.0) then
                write(unit=3, iostat=ios) time_dev
                close(unit=3)
             else
                stop 'file opening pb'
             end if

             deallocate(time_dev)

          else
             if(ios.eq.0) then
                write(unit=3, iostat=ios) this%nodes
                close(unit=3)
             else
                stop 'file opening pb'
             end if
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
     $     time)

          implicit none

          class(bf_layer)            , intent(inout) :: this
          character(*)               , intent(in)    :: filename
          character(*), dimension(ne), intent(in)    :: name_var
          character(*), dimension(ne), intent(in)    :: longname_var
          character(*), dimension(ne), intent(in)    :: unit_var
          integer                    , intent(in)    :: bf_order
          real(rkind)                , intent(in)    :: time

          call print_bf_layer_on_netcdf(
     $         filename,
     $         name_var, longname_var, unit_var,
     $         this%localization, bf_order,
     $         this%grdpts_id,
     $         this%x_map,
     $         this%y_map,
     $         this%nodes,
     $         time)

        end subroutine print_netcdf


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
     $         size(this%nodes,2),
     $         this%alignment,
     $         this%x_map,
     $         this%y_map,
     $         this%grdpts_id)

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
        subroutine compute_time_dev(
     $     this,
     $     td_operators_used,
     $     t,s,p_model,bc_used)

          implicit none

          class(bf_layer)                , intent(inout) :: this
          type(td_operators)             , intent(in)    :: td_operators_used
          real(rkind)                    , intent(in)    :: t
          type(sd_operators)             , intent(in)    :: s
          type(pmodel_eq)                , intent(in)    :: p_model
          type(bc_operators)             , intent(in)    :: bc_used

          call this%bf_compute_used%compute_time_dev(
     $         td_operators_used,
     $         t, this%nodes, this%x_map, this%y_map,
     $         s,p_model,bc_used,
     $         this%grdpts_id,
     $         this%x_borders, this%y_borders,
     $         this%N_bc_sections, this%S_bc_sections)

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
     $         this%x_borders,
     $         this%y_borders,
     $         integration_step_nopt,
     $         this%N_bc_sections,
     $         this%S_bc_sections)

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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the x-borders for the integration 
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param x_borders
        !> integration borders along the x-direction
        !--------------------------------------------------------------
        subroutine set_x_borders(this,x_borders)

          implicit none

          class(bf_layer)             , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: x_borders

          this%x_borders = x_borders

        end subroutine set_x_borders


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the y-borders for the integration 
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param y_borders
        !> integration borders along the y-direction
        !--------------------------------------------------------------
        subroutine set_y_borders(this,y_borders)

          implicit none

          class(bf_layer)             , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: y_borders

          this%y_borders = y_borders

        end subroutine set_y_borders


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the x-borders for the integration 
        !
        !> @date
        !> 31_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@return x_borders
        !> integration borders along the x-direction
        !--------------------------------------------------------------
        function get_x_borders(this) result(x_borders)

          implicit none

          class(bf_layer)             , intent(inout) :: this
          integer(ikind), dimension(2)                :: x_borders

          x_borders = this%x_borders

        end function get_x_borders


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the y-borders for the integration 
        !
        !> @date
        !> 31_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@return y_borders
        !> integration borders along the y-direction
        !--------------------------------------------------------------
        function get_y_borders(this) result(y_borders)

          implicit none

          class(bf_layer)             , intent(inout) :: this
          integer(ikind), dimension(2)                :: y_borders

          y_borders = this%y_borders

        end function get_y_borders

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the N_bc_section attribute
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param N_bc_section
        !> north boundary section
        !--------------------------------------------------------------
        subroutine set_N_bc_sections(this,N_bc_sections)

          implicit none

          class(bf_layer)                            , intent(inout) :: this
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: N_bc_sections

          call MOVE_ALLOC(N_bc_sections,this%N_bc_sections)

        end subroutine set_N_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the S_bc_section attribute
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param S_bc_section
        !> south boundary section
        !--------------------------------------------------------------
        subroutine set_S_bc_sections(this,S_bc_sections)

          implicit none

          class(bf_layer)                            , intent(inout) :: this
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: S_bc_sections

          call MOVE_ALLOC(S_bc_sections,this%S_bc_sections)

        end subroutine set_S_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the N_bc_section attribute
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine remove_N_bc_sections(this)

          implicit none

          class(bf_layer), intent(inout) :: this

          if(allocated(this%N_bc_sections)) then
             deallocate(this%N_bc_sections)
          end if

        end subroutine remove_N_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the S_bc_section attribute
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !--------------------------------------------------------------
        subroutine remove_S_bc_sections(this)

          implicit none

          class(bf_layer), intent(inout) :: this

          if(allocated(this%S_bc_sections)) then
             deallocate(this%S_bc_sections)
          end if

        end subroutine remove_S_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the N_bc_section attribute
        !
        !> @date
        !> 31_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param N_bc_section
        !> south boundary section
        !--------------------------------------------------------------
        subroutine get_N_bc_sections(this,N_bc_sections)

          implicit none

          class(bf_layer)                            , intent(inout) :: this
          integer(ikind), dimension(:,:), allocatable, intent(out)   :: N_bc_sections

          if(allocated(this%N_bc_sections)) then
             allocate(N_bc_sections(
     $            size(this%N_bc_sections,1),
     $            size(this%N_bc_sections,2)))
             N_bc_sections(:,:) = this%N_bc_sections(:,:)
          end if

        end subroutine get_N_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the S_bc_section attribute
        !
        !> @date
        !> 31_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param S_bc_section
        !> south boundary section
        !--------------------------------------------------------------
        subroutine get_S_bc_sections(this,S_bc_sections)

          implicit none

          class(bf_layer)                            , intent(inout) :: this
          integer(ikind), dimension(:,:), allocatable, intent(out)   :: S_bc_sections

          if(allocated(this%S_bc_sections)) then
             allocate(S_bc_sections(
     $            size(this%S_bc_sections,1),
     $            size(this%S_bc_sections,2)))
             S_bc_sections(:,:) = this%S_bc_sections(:,:)
          end if

        end subroutine get_S_bc_sections

      end module bf_layer_class


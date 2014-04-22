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

        use bf_layer_abstract_class , only : bf_layer_abstract

        use bf_layer_allocate_module, only : allocate_bf_layer,
     $                                       reallocate_bf_layer
        use bf_layer_print_module   , only : print_nodes,
     $                                       print_grdpts_id,
     $                                       print_sizes

        use parameters_constant     , only : N,S,E,W,N_W,N_E,S_E,S_W
        use parameters_input        , only : nx,ny,ne,bc_size
        use parameters_kind         , only : ikind, rkind

        private
        public :: bf_layer


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
        type, extends(bf_layer_abstract) :: bf_layer

          contains

          procedure, pass :: get_local_coord
          procedure, pass :: allocate_bf_layer   => allocate_bf_layer_proc
          procedure, pass :: reallocate_bf_layer => reallocate_bf_layer_proc

          procedure, pass :: print_nodes         => print_nodes_proc
          procedure, pass :: print_grdpts_id     => print_grdpts_id_proc
          procedure, pass :: print_sizes         => print_sizes_proc

        end type bf_layer


        contains


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
        subroutine allocate_bf_layer_proc(this, alignment, nodes, neighbors)

          implicit none
          
          class(bf_layer)                 , intent(inout) :: this
          integer(ikind), dimension(2,2)  , intent(in)    :: alignment
          real(rkind)   , dimension(:,:,:), intent(in)    :: nodes
          logical       , dimension(4)    , intent(in)    :: neighbors

          call allocate_bf_layer(this, alignment, nodes, neighbors)

        end subroutine allocate_bf_layer_proc


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
        subroutine reallocate_bf_layer_proc(
     $     this,
     $     border_changes,
     $     nodes,
     $     match_table)
        
          implicit none

          class(bf_layer)                 , intent(inout) :: this
          integer, dimension(2,2)         , intent(in)    :: border_changes
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer, dimension(2)           , intent(out)   :: match_table

          call reallocate_bf_layer(
     $         this,
     $         border_changes,
     $         nodes,
     $         match_table)

        end subroutine reallocate_bf_layer_proc


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
        !>@param filename
        !> name of the binary file where the nodes are
        !> written
        !--------------------------------------------------------------
        subroutine print_nodes_proc(this,filename)

          implicit none

          class(bf_layer), intent(in) :: this
          character(*)   , intent(in) :: filename

          call print_nodes(this,filename)

        end subroutine print_nodes_proc


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine print the grdpt_id table in a binary
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
        !>@param filename
        !> name of the binary file where the grdpt_id are
        !> written
        !--------------------------------------------------------------
        subroutine print_grdpts_id_proc(this,filename)

          implicit none
          
          class(bf_layer), intent(in) :: this
          character(*)   , intent(in) :: filename

          call print_grdpts_id(this, filename)

        end subroutine print_grdpts_id_proc


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine print the sizes of the main table in
        !> a binary file
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
        !>@param filename
        !> name of the binary file where the sizes are
        !> written
        !--------------------------------------------------------------
        subroutine print_sizes_proc(this, filename)

          implicit none
          
          class(bf_layer), intent(in) :: this
          character(*)   , intent(in) :: filename

          call print_sizes(this, filename)

        end subroutine print_sizes_proc

      end module bf_layer_class

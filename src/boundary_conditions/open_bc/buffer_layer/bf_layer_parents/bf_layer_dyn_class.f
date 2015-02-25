      !> @file
      !> module encapsulating the bf_layer_dyn object.
      !> bf_layer_sync enhanced with reallocation procedures
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_layer_sync enhanced with reallocation procedures
      !
      !> @date
      !  23_02_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_dyn_class

        use bf_layer_errors_module, only :
     $       error_mainlayer_id,
     $       error_diff_mainlayer_id

        use bf_layer_sync_class, only :
     $       bf_layer_sync

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

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,
     $       debug

        use parameters_kind, only :
     $       ikind,
     $       rkind

        !> @class bf_layer_dyn
        !> class encapsulating the reallocation intrinsic
        !> procedures for the buffer layer object
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
        !-------------------------------------------------------------
        type, extends(bf_layer_sync) :: bf_layer_dyn

          contains

          procedure, pass :: allocate_bf_layer
          procedure, pass :: reallocate_bf_layer
          procedure, pass :: merge_bf_layer

        end type bf_layer_dyn

        contains


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
        subroutine allocate_bf_layer(
     $       this,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       alignment)

          implicit none
          
          class(bf_layer_dyn)                , intent(inout) :: this
          real(rkind)   , dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2)     , intent(in)    :: alignment


          select case(this%localization)

            case(N)
               call allocate_bf_layer_N(
     $              this%x_map, interior_x_map,
     $              this%y_map, interior_y_map,
     $              this%nodes, interior_nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(S)
               call allocate_bf_layer_S(
     $              this%x_map, interior_x_map,
     $              this%y_map, interior_y_map,
     $              this%nodes, interior_nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(E)
               call allocate_bf_layer_E(
     $              this%x_map, interior_x_map,
     $              this%y_map, interior_y_map,
     $              this%nodes, interior_nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(W)
               call allocate_bf_layer_W(
     $              this%x_map, interior_x_map,
     $              this%y_map, interior_y_map,
     $              this%nodes, interior_nodes,
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
        subroutine reallocate_bf_layer(
     $     this,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     alignment)

          implicit none

          class(bf_layer_dyn)             , intent(inout) :: this
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          integer    , dimension(2,2)     , intent(in)    :: alignment

          select case(this%localization)

            case(N)
               call reallocate_bf_layer_N(
     $              this%x_map, interior_x_map,
     $              this%y_map, interior_y_map,
     $              this%nodes, interior_nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(S)
               call reallocate_bf_layer_S(
     $              this%x_map, interior_x_map,
     $              this%y_map, interior_y_map,
     $              this%nodes, interior_nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(E)
               call reallocate_bf_layer_E(
     $              this%x_map, interior_x_map,
     $              this%y_map, interior_y_map,
     $              this%nodes, interior_nodes,
     $              this%grdpts_id,
     $              this%alignment, alignment)

            case(W)
               call reallocate_bf_layer_W(
     $              this%x_map, interior_x_map,
     $              this%y_map, interior_y_map,
     $              this%nodes, interior_nodes,
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

          class(bf_layer_dyn)                                , intent(inout) :: this
          class(bf_layer_dyn)                                , intent(inout) :: bf_layer2
          real(rkind)         , dimension(nx)                , intent(in)    :: x_map
          real(rkind)         , dimension(ny)                , intent(in)    :: y_map
          real(rkind)         , dimension(nx,ny,ne)          , intent(in)    :: nodes
          integer(ikind)      , dimension(2,2)     , optional, intent(in)    :: alignment

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

      end module bf_layer_dyn_class

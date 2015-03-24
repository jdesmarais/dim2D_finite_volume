      !> @file
      !> module encapsulating the bf_layer object.
      !> bf_layer_time enhanced with removal procedures
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_layer_time enhanced with removal procedures
      !
      !> @date
      ! 23_02_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_class

        use bf_layer_icr_class, only :
     $       bf_layer_icr

        use bf_remove_module, only :
     $       check_if_bf_layer_remains

        use parameters_input, only :
     $       nx,ny,ne
        
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
        !> @param can_remain
        !> logical identifying whether the buffer layer based on its
        !> grid points at the edge with the interior domain are such
        !> that the buffer layer can be removed
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
        !> remove the buffer layer by deallocating the main tables
        !-------------------------------------------------------------
        type, extends(bf_layer_icr) :: bf_layer

          logical, private :: can_remain

          contains

          procedure, pass :: set_remain_status
          procedure, pass :: get_remain_status
          procedure, pass :: should_remain

        end type bf_layer

        contains


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
        function should_remain(
     $     this,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     p_model)

          implicit none

          class(bf_layer)                 , intent(in) :: this
          real(rkind), dimension(nx)      , intent(in) :: interior_x_map
          real(rkind), dimension(ny)      , intent(in) :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in) :: interior_nodes
          type(pmodel_eq)                 , intent(in) :: p_model
          logical                                      :: should_remain
          
          integer(ikind), dimension(2) :: bf_match_table

          bf_match_table = this%get_general_to_local_coord_tab()

          should_remain = check_if_bf_layer_remains(
     $         this%localization,
     $         this%alignment,
     $         bf_match_table,
     $         this%grdpts_id,
     $         this%x_map,
     $         this%y_map,
     $         this%nodes,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)

        end function should_remain

      end module bf_layer_class


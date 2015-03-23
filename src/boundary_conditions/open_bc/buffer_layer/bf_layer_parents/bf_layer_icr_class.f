      !> @file
      !> module encapsulating the bf_layer_icr object.
      !> bf_layer_grdpts_id_update enhanced with functions
      !> to determine whether the domain should be extended
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_layer_grdpts_id_update enhanced with functions
      !> to determine whether the domain should be extended
      !
      !> @date
      ! 24_02_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_icr_class

        use bf_layer_grdpts_id_update_class, only :
     $       bf_layer_grdpts_id_update

        use bf_newgrdpt_verification_module, only :
     $       are_grdpts_available

        use parameters_bf_layer, only :
     $       BF_SUCCESS

        use parameters_kind, only :
     $       ikind

        use pmodel_eq_class, only :
     $       pmodel_eq


        implicit none


        private
        public :: bf_layer_icr


        !> @class bf_layer_icr_class
        !> bf_layer_grdpts_id_update enhanced with functions
        !> to determine whether the domain should be extended
        !
        !> @param is_node_activated
        !> determine whether a grid point is activated
        !
        !> @param get_bc_sections
        !> extract the boundary section identifying the position
        !> of the bc_interior_pt grid-points, potentially activated
        !-------------------------------------------------------------
        type, extends(bf_layer_grdpts_id_update) :: bf_layer_icr

           contains

           procedure, pass :: is_node_activated
           procedure, pass :: get_bc_sections 

        end type bf_layer_icr


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the grid-point is activated
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param loc_coords
        !> local coordinates of the grid-point checked, i.e.
        !> coordinates of the grid-point in the current buffer layer
        !
        !>@param p_model
        !> physical model
        !
        !>@param ierror
        !> determine whether the determination of the activation was
        !> possible or not
        !
        !>@return node_activated
        !> logical identifying whether the grid-point is activated or not
        !--------------------------------------------------------------
        function is_node_activated(
     $     this,
     $     loc_coords,
     $     p_model,
     $     ierror)
     $     result(node_activated)

          implicit none

          class(bf_layer_icr)         , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: loc_coords
          type(pmodel_eq)             , intent(in) :: p_model
          logical                     , intent(out):: ierror
          logical                                  :: node_activated

          logical :: grdpts_available
          

          node_activated = .false.


          !> check whether there are indeed grid-points to
          !> determine whether the grid-point is activated
          !> or not
          if(.not.allocated(this%nodes)) then
             
             ierror = .not.BF_SUCCESS
             
          else

             !> check whether there are enough grid-points around
             !> the central grid-point checked to determine whether
             !> the grid-point is activated or not
             if((loc_coords(1).lt.2).or.
     $          (loc_coords(1).gt.(size(this%nodes,1)-1)).or.
     $          (loc_coords(2).lt.2).or.
     $          (loc_coords(2).gt.(size(this%nodes,2)-1))) then

                ierror = .not.BF_SUCCESS

             else

                ierror = BF_SUCCESS

                !> check whether all grdpts are available to evaluate
                !> whether the node is activated or not
                grdpts_available = are_grdpts_available(
     $               this%grdpts_id,
     $               reshape((/
     $                  loc_coords(1)-1,
     $                  loc_coords(2)-1,
     $                  loc_coords(1)+1,
     $                  loc_coords(2)+1/),
     $               (/2,2/)))

                !> if all grid-points are available, we check whether
                !> the grid-point is indeed activated
                if(grdpts_available) then

                   !> check whether the node is activated
                   node_activated = p_model%are_openbc_undermined(
     $                  this%x_map(loc_coords(1)-1:loc_coords(1)+1),
     $                  this%y_map(loc_coords(2)-1:loc_coords(2)+1),
     $                  this%nodes(loc_coords(1)-1:loc_coords(1)+1,
     $                             loc_coords(2)-1:loc_coords(2)+1,
     $                             :))

                !> otherwise, the grid-point is set to desactivated
                else
                   
                   node_activated = .false.

                end if

             end if

          end if

        end function is_node_activated


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> extract the boundary section identifying the position
        !> of the bc_interior_pt grid-points, potentially activated
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_grdpts_id_update enhanced with functions
        !> to determine whether the domain should be extended
        !
        !>@param bc_sections
        !> array with the positions and type of boundary sections
        !--------------------------------------------------------------
        subroutine get_bc_sections(this,bc_sections)

          implicit none

          class(bf_layer_icr)                        , intent(inout) :: this
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections


          if(allocated(this%bc_sections)) then

             call MOVE_ALLOC(this%bc_sections,bc_sections)

          end if

        end subroutine get_bc_sections        

      end module bf_layer_icr_class

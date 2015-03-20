      !> @file
      !> bf_interface_grdpts_id_update augmented with procedures
      !> identifying the sublayer corresponding to grid points
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_interface_grdpts_id_update augmented with procedures
      !> identifying the sublayer corresponding to grid points
      !
      !> @date
      ! 19_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interface_coords_class

        use bf_increase_coords_module, only :
     $       get_mainlayer_coord

        use bf_interface_grdpts_id_update_class, only :
     $       bf_interface_grdpts_id_update

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_constant, only :
     $       interior

        use parameters_kind, only :
     $       ikind


        implicit none

        private
        public :: bf_interface_coords


        !>@class bf_interface_coords
        !> bf_interface_grdpts_id_update augmented with procedures
        !> identifying the sublayer corresponding to grid points
        !
        !>@param get_bf_layer_from_gen_coords
        !> identify the buffer layer corresponding to the grid point
        !> general coordinates
        !---------------------------------------------------------------
        type, extends(bf_interface_grdpts_id_update) :: bf_interface_coords

          contains

          procedure, pass :: get_bf_layer_from_gen_coords

        end type bf_interface_coords

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> identify the buffer layer corresponding to the grid point
        !> general coordinates
        !
        !> @date
        !> 19_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> identifying the sublayer corresponding to grid points
        !--------------------------------------------------------------
        function get_bf_layer_from_gen_coords(
     $       this,
     $       gen_coords,
     $       tolerance,
     $       mainlayer_id)
     $       result(bf_sublayer_ptr)

          implicit none

          class(bf_interface_coords)  , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: gen_coords
          integer(ikind), optional    , intent(in) :: tolerance
          integer       , optional    , intent(in) :: mainlayer_id
          type(bf_sublayer), pointer               :: bf_sublayer_ptr

          integer :: tolerance_op
          integer :: mainlayer_id_op


          if(present(tolerance)) then
             tolerance_op = tolerance
          else
             tolerance_op = 0
          end if


          !1) identify the mainlayer to which the grid-point
          !   as well as the neighboring grid points in a radius
          !   of bc_size should belong
          if(present(mainlayer_id)) then
             mainlayer_id_op = mainlayer_id
          else
             mainlayer_id_op = get_mainlayer_coord(gen_coords)
          end if

          if(mainlayer_id_op.eq.interior) then
             print '(''bf_interface_coords'')'
             print '(''get_bf_layer_from_gen_coords'')'
             print '(''this grid-point belongs to the interior'')'
             print '(''gen_coords: '',2I2)', gen_coords
             stop ''
          end if


          !2) extract the sublayer to which the grid point can
          !   belong
          bf_sublayer_ptr => this%mainlayer_pointers(mainlayer_id_op)%get_sublayer_from_gen_coords(
     $         gen_coords,
     $         tolerance=tolerance_op,
     $         no_check_ID=.true.)

        end function get_bf_layer_from_gen_coords

      end module bf_interface_coords_class

      !> @file
      !> bf_interface_time augmented with procedures
      !> updating the configuration of the grid points
      !> around a new interior_pt
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_interface_time augmented with procedures
      !> updating the configuration of the grid points
      !> around a new interior_pt
      !
      !> @date
      ! 19_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interface_grdpts_id_update_class

        use bf_interface_time_class, only :
     $       bf_interface_time

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: bf_interface_grdpts_id_update


        !> @class bf_interface_grdpts_id_update
        !> bf_interface_time augmented with procedures
        !> updating the configuration of the grid points
        !> around a new interior_pt
        !
        !> @param update_grdpts_id_in_bf_layer
        !> update the configuration of the grid points
        !> around the new interior point
        !------------------------------------------------------------
        type, extends(bf_interface_time) :: bf_interface_grdpts_id_update

          contains

          procedure, pass :: update_grdpts_id_in_bf_layer

        end type bf_interface_grdpts_id_update


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> turn the grdpts_id identified by general coordinates
        !> from bc_interior_pt to interior_pt and reallocate
        !> the buffer layer such that the neighboring points
        !> around it are allocated. Then compute these new
        !> grid points.
        !
        !> @date
        !> 19_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> nbf_interface_newgrdpt encapsulating the 
        !> links to neighboring buffer layers and the
        !> functions for the computation of the new grid
        !> points in the buffer layer
        !
        !> @param bf_sublayer_updated
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !> @param p_model
        !> physical model
        !
        !> @param t
        !> time 
        !
        !> @param dt
        !> time step
        !
        !> @param interior_x_map
        !> array with x-coordinates of the interior domain
        !
        !> @param interior_y_map
        !> array with y-coordinates of the interior domain
        !
        !> @param interior_nodes0
        !> array with the grid point data at t=t-dt
        !
        !> @param interior_nodes1
        !> array with the grid point data at t=t
        !
        !>@param selected_grdpts
        !> list containing the general coordinates of the grid points to be
        !> turned from bc_interior_pt to interior_pt
        !--------------------------------------------------------------
        subroutine update_grdpts_id_in_bf_layer(
     $       this,
     $       p_model,
     $       t,
     $       dt,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes0,
     $       interior_nodes1,
     $       bf_sublayer_ptr,
     $       selected_grdpts)

          implicit none

          class(bf_interface_grdpts_id_update), intent(inout) :: this
          type(pmodel_eq)                     , intent(in)    :: p_model
          real(rkind)                         , intent(in)    :: t
          real(rkind)                         , intent(in)    :: dt
          real(rkind)   , dimension(nx)       , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)       , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne) , intent(in)    :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne) , intent(in)    :: interior_nodes1
          type(bf_sublayer)                   , intent(inout) :: bf_sublayer_ptr
          integer(ikind), dimension(:,:)      , intent(in)    :: selected_grdpts


          call this%mainlayer_interfaces%update_grdpts_id_in_bf_layer(
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         bf_sublayer_ptr,
     $         selected_grdpts)

        end subroutine update_grdpts_id_in_bf_layer

      end module bf_interface_grdpts_id_update_class

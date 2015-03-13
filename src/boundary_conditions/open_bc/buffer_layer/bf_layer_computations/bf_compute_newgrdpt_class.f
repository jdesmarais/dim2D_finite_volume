      !> @file
      !> bf_compute object augmented with procedures to
      !> compute the new grid points belonging to the
      !> buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_compute object augmented with procedures to
      !> compute the new grid points belonging to the
      !> buffer layer
      !
      !> @date
      !> 13_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_compute_newgrdpt_class

        use bf_compute_time_class, only :
     $       bf_compute_time

        use bf_newgrdpt_dispatch_module, only :
     $       ask_bf_layer_to_compute_newgrdpt

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: bf_compute_newgrdpt


        !> @class bf_compute_newgrdpt
        !> class encapsulating the intrinsic procedures to compute the 
        !> new grid points from the data at the previous steps
        !
        !> @param compute_newgrdpt
        !> attempt to compute the new grid point
        !---------------------------------------------------------------
        type, extends(bf_compute_time) :: bf_compute_newgrdpt

          contains

          procedure, pass :: compute_newgrdpt

        end type bf_compute_newgrdpt

        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> ask the buffer layer to compute the new grid point,
        !> if the buffer layer does not have enough grid points
        !> to determine the procedure to apply or enough grid
        !> points to compute the new grid point itself, it sends
        !> an error explaining the limitations encountered
        !
        !> @date
        !> 13_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute_time augmented with procedures for the
        !> computation of the new grid points
        !
        !>@param p_model
        !> physical model
        !
        !>@param t
        !> time
        !
        !>@param dt
        !> time step
        !
        !>@param bf_localization
        !> cardinal coordinate identifying the position of the
        !> buffer layer
        !
        !>@param bf_can_exchange_with_neighbor1
        !> logical identifying whether the buffer layer can
        !> exchange grid points with another buffer layer of
        !> type neighbor1
        !
        !>@param bf_can_exchange_with_neighbor2
        !> logical identifying whether the buffer layer can
        !> exchange grid points with another buffer layer of
        !> type neighbor2
        !
        !>@param bf_alignment1
        !> relative position of the buffer layer compared to
        !> the interior domain at t
        !
        !>@param bf_x_map1
        !> x-coordinates for the buffer layer at t
        !
        !>@param bf_y_map1
        !> y-coordinates for the buffer layer at t
        !
        !>@param bf_nodes1
        !> nodes for the buffer layer at t
        !
        !>@param bf_newgrdpt_coords1
        !> coordinates of the new grid-point in the reference
        !> frame of the buffer layer at t
        !
        !>@param nb_procedures
        !> number of procedures needed to compute the new
        !> grid-point
        !
        !>@param procedure_type
        !> integer identifying the procedure used to compute
        !> the new grid point
        !
        !>@param gradient_type
        !> integer identifying the gradient used to compute
        !> the new grid point
        !
        !>@param grdpts_available
        !> logical determining whether there are enough grid points
        !> to compute the new grid point
        !
        !>@param data_needed_bounds0
        !> bounds determining the extent of the data needed around
        !> the position of the new grid point for its computation
        !
        !>@return ierror
        !> logical determining whether the buffer layer managed to
        !> compute the new grid point or whether a temporary array
        !> collecting the data is required
        !--------------------------------------------------------------
        function compute_newgrdpt(
     $       this,
     $       p_model,
     $       t,
     $       dt,
     $       bf_localization,
     $       bf_can_exchange_with_neighbor1,
     $       bf_can_exchange_with_neighbor2,
     $       bf_alignment1,
     $       bf_x_map1,
     $       bf_y_map1,
     $       bf_nodes1,
     $       bf_newgrdpt_coords1,
     $       nb_procedures,
     $       procedure_type,
     $       gradient_type,
     $       grdpts_available,
     $       data_needed_bounds0)
     $       result(ierror)

          implicit none

          class(bf_compute_newgrdpt)         , intent(in)    :: this
          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)                        , intent(in)    :: t
          real(rkind)                        , intent(in)    :: dt
          integer                            , intent(in)    :: bf_localization
          logical                            , intent(in)    :: bf_can_exchange_with_neighbor1
          logical                            , intent(in)    :: bf_can_exchange_with_neighbor2
          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment1
          real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map1
          real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map1
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: bf_nodes1
          integer(ikind), dimension(2)       , intent(in)    :: bf_newgrdpt_coords1
          integer                            , intent(out)   :: nb_procedures
          integer       , dimension(4)       , intent(out)   :: procedure_type
          integer       , dimension(4)       , intent(out)   :: gradient_type
          logical       , dimension(4)       , intent(out)   :: grdpts_available
          integer(ikind), dimension(2,2)     , intent(out)   :: data_needed_bounds0
          integer                                            :: ierror


          ierror = ask_bf_layer_to_compute_newgrdpt(
     $         p_model,
     $         t,
     $         dt,
     $         bf_localization,
     $         bf_can_exchange_with_neighbor1,
     $         bf_can_exchange_with_neighbor2,
     $         this%alignment_tmp,
     $         this%grdpts_id_tmp,
     $         this%x_map_tmp,
     $         this%y_map_tmp,
     $         this%nodes_tmp,
     $         bf_alignment1,
     $         bf_x_map1,
     $         bf_y_map1,
     $         bf_nodes1,
     $         bf_newgrdpt_coords1,
     $         nb_procedures,
     $         procedure_type,
     $         gradient_type,
     $         grdpts_available,
     $         data_needed_bounds0)

        end function compute_newgrdpt

      end module bf_compute_newgrdpt_class

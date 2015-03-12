      !determine whether it is possible to evaluate
      !the new grdpt using the data from the buffer
      !layer and compute it if it is possible, otherwise
      !ask to compute the new grid-point from the main
      !structure (interior+buffer layer) to gather data
      module bf_newgrdpt_dispatch_module

        use 


        private
        public :: 


        contains


                !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid points
        !
        !> @date
        !> 18_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param grdpts_id
        !> identification of the role of the grid points
        !
        !>@param nodes
        !> grid points of the buffer layer whose time derivatives
        !> are computed
        !
        !>@param dt
        !> time
        !
        !>@param x_borders
        !> border indices along the x-direction for the time 
        !> integration
        !
        !>@param y_borders
        !> border indices along the y-direction for the time 
        !> integration
        !
        !>@param integration_step_nopt
        !> function determining the procedure for the time integration
        !--------------------------------------------------------------
        function compute_newgrdpt_with_buffer_layer(
     $       p_model, t, dt,
     $       
     $       newgrdpt_coords1,
     $       
     $       bf_localization,
     $       bf_alignment1,
     $       bf_nodes1,
     $       can_exchange_with_neighbor1,
     $       can_exchange_with_neighbor2,
     $       
     $       nb_procedures,
     $       procedure_type,
     $       gradient_type,
     $       grdpts_available,
     $       tmp_array_needed,
     $       data_needed_bounds0)
     $     result(ierror)

          implicit none

          type(pmodel_eq)                    , intent(in)  :: p_model
          real(rkind)                        , intent(in)  :: t
          real(rkind)                        , intent(in)  :: dt
          integer(ikind), dimension(2)       , intent(in)  :: newgrdpt_coords1
          integer                            , intent(in)  :: bf_localization
          integer(ikind), dimension(2,2)     , intent(in)  :: bf_alignment1
          real(rkind)   , dimension(:,:,:)   , intent(in)  :: bf_nodes1
          logical                            , intent(in)  :: can_exchange_with_neighbor1
          logical                            , intent(in)  :: can_exchange_with_neighbor2
          integer                            , intent(out) :: nb_procedures
          integer       , dimension(4)       , intent(out) :: procedure_type
          integer       , dimension(4)       , intent(out) :: gradient_type
          logical       , dimension(4)       , intent(out) :: grdpts_available
          logical       , dimension(4)       , intent(out) :: tmp_array_needed
          integer(ikind), dimension(2,2)     , intent(out) :: data_needed_bounds0
          integer                                          :: ierror

          integer               :: k
          integer, dimension(2) :: newgrdpt_coords0


          !> get the coordinates of the new grdpt in the
          !> reference frame at t-dt such that we can
          !> estimate the procedure needed for the 
          !> computation of the new grid point with the
          !> grdpts_id configuration at t-dt
          do k=1,2
             newgrdpt_coords0(k) = newgrdpt_coords1(1) + bf_alignment1(k,1) - this%alignment_tmp(k,1)
          end do


          !> determine which procedure should be used
          !> to compute the new grid point if this is
          !> possible to estimate the procedure
          ierror_proc = get_newgrdpt_procedure_from_grdpts_id0(
     $         bf_localization,
     $         bf_grdpts_id0,
     $         newgrdpt_coords0,
     $         can_exchange_with_neighbor1,
     $         can_exchange_with_neighbor2,
     $         nb_procedures,
     $         procedure_type,
     $         gradient_type)

          
          !> if there were not enough grid points available
          !> to determine the procedure to compute the new
          !> grid point, the error message is initialized
          if(ierror_proc.neqv.BF_SUCCESS) then

             ierror = NEWGRDPT_PROC_NBGRDPTS_ERROR

          !> otherwise we try to compute the data for the
          !> new grid point using the available data in
          !> the buffer layer and the interior domain
          else
             
             ierror_proc = compute_newgrdpt_from_bf_layer(
     $            p_model,
     $            t,
     $            dt,
     $            bf_alignment0,
     $            bf_grdpts_id0,
     $            bf_x_map0,
     $            bf_y_map0,
     $            bf_nodes0,
     $            bf_alignment1,
     $            bf_x_map1,
     $            bf_y_map1,
     $            bf_nodes1,
     $            newgrdpt_coords0,
     $            newgrdpt_coords1,
     $            nb_procedures,
     $            procedure_type,
     $            gradient_type,
     $            grdpts_available,
     $            tmp_array_needed,
     $            data_needed_bounds0)

             !> if there were not enough grid points to
             !> compute the new gridpoint, the error
             !> message is initialized
             if(ierror_proc.neqv.BF_SUCCESS) then

                ierror = NEWGRDPT_DATA_NBGRDPTS_ERROR

             !> otherwise, the error is initialized
             !> to success
             else

                ierror = NEWGRDPT_SUCCESS

             end if

          end if

        end function compute_newgrdpt









      end module bf_newgrdpt_dispatch_module

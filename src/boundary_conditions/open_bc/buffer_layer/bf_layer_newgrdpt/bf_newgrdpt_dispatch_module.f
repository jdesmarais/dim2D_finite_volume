      !determine whether it is possible to evaluate
      !the new grdpt using the data from the buffer
      !layer and compute it if it is possible, otherwise
      !ask to compute the new grid-point from the main
      !structure (interior+buffer layer) to gather data
      module bf_newgrdpt_dispatch_module

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use bf_layer_extract_module, only :
     $       get_bf_layer_match_table,
     $       get_grdpts_id_from_interior,
     $       get_grdpts_id_from_bf_layer

        use bf_newgrdpt_class, only :
     $       bf_newgrdpt

        use bf_newgrdpt_procedure_module, only :
     $       get_newgrdpt_procedure

        use bf_newgrdpt_verification_module, only :
     $       are_grdpts_available,
     $       get_newgrdpt_verification_bounds
        
        use parameters_bf_layer, only :
     $       BF_SUCCESS,     
     $       no_gradient_type,
     $       
     $       NEWGRDPT_NO_ERROR,
     $       NEWGRDPT_PROC_NBGRDPTS_ERROR,
     $       NEWGRDPT_DATA_NBGRDPTS_ERROR

        use parameters_constant, only :
     $       N,S,E,W,
     $       no_bc_procedure_type

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq


        private
        public ::
     $       ask_bf_layer_to_compute_newgrdpt,
     $       get_newgrdpt_procedure_from_bf_layer,
     $       compute_newgrdpt_from_bf_layer,
     $       are_grdpts_available_to_get_newgrdpt_data,
     $       are_grdpts_available_to_get_newgrdpt_proc


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
        !>@param bf_alignment0
        !> relative position of the buffer layer compared to
        !> the interior domain at t-dt
        !
        !>@param bf_grdpts_id0
        !> integer array identifying the configuration of the
        !> grid-points at t-dt
        !
        !>@param bf_x_map0
        !> x-coordinates for the buffer layer at t-dt
        !
        !>@param bf_y_map0
        !> y-coordinates for the buffer layer at t-dt
        !
        !>@param bf_nodes0
        !> nodes for the buffer layer at t-dt
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
        function ask_bf_layer_to_compute_newgrdpt(
     $       p_model,
     $       t,
     $       dt,
     $       bf_localization,
     $       bf_can_exchange_with_neighbor1,
     $       bf_can_exchange_with_neighbor2,
     $       bf_alignment0,
     $       bf_grdpts_id0,
     $       bf_x_map0,
     $       bf_y_map0,
     $       bf_nodes0,
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

          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)                        , intent(in)    :: t
          real(rkind)                        , intent(in)    :: dt
          integer                            , intent(in)    :: bf_localization
          logical                            , intent(in)    :: bf_can_exchange_with_neighbor1
          logical                            , intent(in)    :: bf_can_exchange_with_neighbor2
          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment0
          integer       , dimension(:,:)     , intent(in)    :: bf_grdpts_id0
          real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map0
          real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map0
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes0
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

          integer               :: k
          integer, dimension(2) :: bf_newgrdpt_coords0
          logical               :: ierror_proc


          !> get the coordinates of the new grdpt in the
          !> reference frame at t-dt such that we can
          !> estimate the procedure needed for the 
          !> computation of the new grid point with the
          !> grdpts_id configuration at t-dt
          do k=1,2
             bf_newgrdpt_coords0(k) = bf_newgrdpt_coords1(k) + bf_alignment1(k,1) - bf_alignment0(k,1)
          end do


          !> determine which procedure should be used
          !> to compute the new grid point if this is
          !> possible to estimate the procedure
          ierror_proc = get_newgrdpt_procedure_from_bf_layer(
     $         bf_localization,
     $         bf_alignment0,
     $         bf_grdpts_id0,
     $         bf_newgrdpt_coords0,
     $         bf_can_exchange_with_neighbor1,
     $         bf_can_exchange_with_neighbor2,
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
     $            bf_newgrdpt_coords0,
     $            bf_newgrdpt_coords1,
     $            nb_procedures,
     $            procedure_type,
     $            gradient_type,
     $            grdpts_available,
     $            data_needed_bounds0)

             !> if there were not enough grid points to
             !> compute the new gridpoint, the error
             !> message is initialized
             if(ierror_proc.neqv.BF_SUCCESS) then

                ierror = NEWGRDPT_DATA_NBGRDPTS_ERROR

             !> otherwise, the error is initialized
             !> to success
             else

                ierror = NEWGRDPT_NO_ERROR

             end if

          end if

        end function ask_bf_layer_to_compute_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> ask the buffer layer to determine which procedures
        !> are needed to compute the new grid-point or whether
        !> it is compulsory to create a temporary array to
        !> gather information from several buffer layers
        !
        !> @date
        !> 13_03_2015 - initial version - J.L. Desmarais
        !
        !>@param bf_localization
        !> cardinal coordinate identifying the position of the
        !> buffer layer
        !
        !>@param bf_alignment0
        !> relative position of the buffer layer compared to
        !> the interior domain at t-dt
        !
        !>@param bf_grdpts_id0
        !> integer array identifying the configuration of the
        !> grid-points at t-dt
        !
        !>@param bf_newgrdpt_coords0
        !> coordinates of the new grid-point in the reference
        !> frame of the buffer layer at t-dt
        !
        !>@param bf_can_exchange_with_neighbor1
        !> determine whether the buffer layer can exchange with
        !> the buffer layer of type neighbor1
        !
        !>@param bf_can_exchange_with_neighbor2
        !> determine whether the buffer layer can exchange with
        !> the buffer layer of type neighbor2
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
        !>@return ierror
        !> logical determining whether the buffer layer managed to
        !> determein which procedures are needed to compute the new grid
        !> point
        !--------------------------------------------------------------
        function get_newgrdpt_procedure_from_bf_layer(
     $       bf_localization,
     $       bf_alignment0,
     $       bf_grdpts_id0,
     $       bf_newgrdpt_coords0,
     $       bf_can_exchange_with_neighbor1,
     $       bf_can_exchange_with_neighbor2,
     $       nb_procedures,
     $       procedure_type,
     $       gradient_type)
     $       result(ierror)

          implicit none

          integer                       , intent(in)  :: bf_localization
          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment0
          integer       , dimension(:,:), intent(in)  :: bf_grdpts_id0
          integer(ikind), dimension(2)  , intent(in)  :: bf_newgrdpt_coords0
          logical                       , intent(in)  :: bf_can_exchange_with_neighbor1
          logical                       , intent(in)  :: bf_can_exchange_with_neighbor2
          integer                       , intent(out) :: nb_procedures
          integer       , dimension(4)  , intent(out) :: procedure_type
          integer       , dimension(4)  , intent(out) :: gradient_type
          logical                                     :: ierror

          integer                        :: size_x
          integer                        :: size_y
          logical                        :: ierror_proc
          logical                        :: grdpts_available

          integer       , dimension(3,3) :: grdpts_id0_tmp
          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2,2) :: gen_coords


          size_x = size(bf_grdpts_id0,1)
          size_y = size(bf_grdpts_id0,2)

          nb_procedures     = 0
          procedure_type(1) = no_bc_procedure_type
          gradient_type(1)  = no_gradient_type


          !> if the grid points [i-1,i+1]x[j-1,j+1] are 
          !> available around the new gridpoint [i,j], the
          !> procedure needed for the new grdpt can be
          !> immediately determined
          !------------------------------------------------------------
          if(
     $         (bf_newgrdpt_coords0(1).gt.1).and.
     $         (bf_newgrdpt_coords0(1).lt.size_x).and.
     $         (bf_newgrdpt_coords0(2).gt.1).and.
     $         (bf_newgrdpt_coords0(2).lt.size_y) ) then

             ierror_proc = get_newgrdpt_procedure(
     $            bf_newgrdpt_coords0(1),
     $            bf_newgrdpt_coords0(2),
     $            bf_grdpts_id0,
     $            nb_procedures,
     $            procedure_type,
     $            gradient_type)

             if(ierror_proc.neqv.BF_SUCCESS) then
                print '(''bf_compute_newgrdpt_class'')'
                print '(''get_newgrdpt_procedure_from_grdpts_id0'')'
                print '(''error when evaluating the newgrdpt_procedure'')'
                stop ''
             end if

             ierror = BF_SUCCESS

          !> if the grid points [i-1,i+1]x[j-1,j+1] are not
          !> available around the new gridpoint [i,j], we
          !> determine whether the grid points available are
          !> enough to estimate the procedure (by combining
          !> with the information from the interior domain)
          !------------------------------------------------------------
          else
             
             !> determine whether there are enough grdpts_id to
             !> determine the procedure for the new gridpoint
             grdpts_available = are_grdpts_available_to_get_newgrdpt_proc(
     $            bf_localization,
     $            size_x, size_y,
     $            bf_can_exchange_with_neighbor1,
     $            bf_can_exchange_with_neighbor2,
     $            bf_newgrdpt_coords0)


             !> if there are enough grid points to evaluate the procedure
             !> a temporary array is created using the grdpts_id from the
             !> interior domain and the procedure is evaluated
             if(grdpts_available) then

                !> determine the general coordinates for the
                !> corners of the grdpts_id extracted
                match_table = get_bf_layer_match_table(bf_alignment0)

                gen_coords(1,1) = bf_newgrdpt_coords0(1)-1 + match_table(1)
                gen_coords(1,2) = bf_newgrdpt_coords0(1)+1 + match_table(1)
                gen_coords(2,1) = bf_newgrdpt_coords0(2)-1 + match_table(2)
                gen_coords(2,2) = bf_newgrdpt_coords0(2)+1 + match_table(2)
                

                !> extract the grdpts_id to determine the new grdpt
                !> procedure

                !> extract the grdpts_id from the interior domain
                call get_grdpts_id_from_interior(
     $               grdpts_id0_tmp,
     $               gen_coords)

                !> extract the grdpts_id from the grdpts_id0
                call get_grdpts_id_from_bf_layer(
     $               grdpts_id0_tmp,
     $               gen_coords,
     $               bf_alignment0,
     $               bf_grdpts_id0)

                !> determine the procedure for the computation
                !> of the new grid point
                ierror_proc = get_newgrdpt_procedure(
     $               2,2,
     $               grdpts_id0_tmp,
     $               nb_procedures,
     $               procedure_type,
     $               gradient_type)

                if(ierror_proc.neqv.BF_SUCCESS) then
                   print '(''bf_compute_newgrdpt_class'')'
                   print '(''get_newgrdpt_procedure_from_grdpts_id0'')'
                   print '(''error extracting procedure for newgrdpt'')'
                   print '(''with temporary array'')'
                   stop ''
                end if

                ierror = BF_SUCCESS

             !> otherwise, the error message is initialized with
             !> not_enough_grdpts_proc_error
             else
                ierror = .not.BF_SUCCESS
             end if
                
          end if

        end function get_newgrdpt_procedure_from_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> verify if the grid points needed to compute the
        !> new grid point are available and try to compute
        !> the new grid point
        !
        !> @date
        !> 13_03_2015 - initial version - J.L. Desmarais
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
        !>@param bf_alignment0
        !> relative position of the buffer layer compared to
        !> the interior domain at t-dt
        !
        !>@param bf_grdpts_id0
        !> integer array identifying the configuration of the
        !> grid-points at t-dt
        !
        !>@param bf_x_map0
        !> x-coordinates for the buffer layer at t-dt
        !
        !>@param bf_y_map0
        !> y-coordinates for the buffer layer at t-dt
        !
        !>@param bf_nodes0
        !> nodes for the buffer layer at t-dt
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
        !>@param bf_newgrdpt_coords0
        !> coordinates of the new grid-point in the reference
        !> frame of the buffer layer at t-dt
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
        function compute_newgrdpt_from_bf_layer(
     $     p_model,
     $     t,
     $     dt,
     $     bf_alignment0,
     $     bf_grdpts_id0,
     $     bf_x_map0,
     $     bf_y_map0,
     $     bf_nodes0,
     $     bf_alignment1,
     $     bf_x_map1,
     $     bf_y_map1,
     $     bf_nodes1,
     $     bf_newgrdpt_coords0,
     $     bf_newgrdpt_coords1,
     $     nb_procedures,
     $     procedure_type,
     $     gradient_type,
     $     grdpts_available,
     $     data_needed_bounds0)
     $     result(ierror)

          implicit none

          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind)                     , intent(in)    :: t
          real(rkind)                     , intent(in)    :: dt
          integer(ikind), dimension(2,2)  , intent(in)    :: bf_alignment0
          integer       , dimension(:,:)  , intent(in)    :: bf_grdpts_id0
          real(rkind)   , dimension(:)    , intent(in)    :: bf_x_map0
          real(rkind)   , dimension(:)    , intent(in)    :: bf_y_map0
          real(rkind)   , dimension(:,:,:), intent(in)    :: bf_nodes0
          integer(ikind), dimension(2,2)  , intent(in)    :: bf_alignment1
          real(rkind)   , dimension(:)    , intent(in)    :: bf_x_map1
          real(rkind)   , dimension(:)    , intent(in)    :: bf_y_map1
          real(rkind)   , dimension(:,:,:), intent(inout) :: bf_nodes1
          integer(ikind), dimension(2)    , intent(in)    :: bf_newgrdpt_coords0
          integer(ikind), dimension(2)    , intent(in)    :: bf_newgrdpt_coords1
          integer                         , intent(in)    :: nb_procedures
          integer       , dimension(4)    , intent(in)    :: procedure_type
          integer       , dimension(4)    , intent(in)    :: gradient_type
          logical       , dimension(4)    , intent(out)   :: grdpts_available
          integer(ikind), dimension(2,2)  , intent(out)   :: data_needed_bounds0
          integer                                         :: ierror


          logical               :: all_grdpts_available
          logical               :: tmp_array_needed
          integer               :: k
          type(bf_newgrdpt)     :: bf_newgrdpt_used


          !> for each procedure_type verify that there are
          !> enough grid points available
          all_grdpts_available = .true.

          do k=1, nb_procedures

             call are_grdpts_available_to_get_newgrdpt_data(
     $            bf_grdpts_id0,
     $            bf_newgrdpt_coords0,
     $            procedure_type(k),
     $            gradient_type(k),
     $            data_needed_bounds0,
     $            tmp_array_needed,
     $            grdpts_available(k))

             if(.not.grdpts_available(k)) then
                all_grdpts_available = .false.
             end if

          end do

          
          !> if all grid points are available to compute the
          !> new grid point, the new grdpt is directly computed
          !> from the buffer layer
          if(all_grdpts_available) then

             bf_nodes1(
     $            bf_newgrdpt_coords1(1),
     $            bf_newgrdpt_coords1(2),
     $            :) =
     $       bf_newgrdpt_used%compute_newgrdpt(
     $            p_model, t, dt,
     $            bf_alignment0, bf_x_map0, bf_y_map0, bf_nodes0,
     $            bf_alignment1, bf_x_map1, bf_y_map1, bf_nodes1,
     $            bf_newgrdpt_coords1(1),bf_newgrdpt_coords1(2),
     $            nb_procedures, procedure_type, gradient_type)

             ierror = BF_SUCCESS


          !> otherwise we will need to create a temporary array to compute
          !> the new grdpt: if only one temporary array is needed,
          !> the overcost of having one big array for all data is negligeable
          !> compared to using a temporary array from the interior domain +
          !> the buffer layer, so we prefer to compute everything from the same
          !> big temporary array for all procedures for k in [1,nb_procedures]
          else

             ierror = .not.BF_SUCCESS

          end if

        end function compute_newgrdpt_from_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine whether there are enough grid points to
        !> compute the new grid point if the procedure is 
        !> identified by procedure_type and the gradient by
        !> gradient_type
        !
        !> @date
        !> 13_03_2015 - initial version - J.L. Desmarais
        !
        !>@param bf_grdpts_id0
        !> integer array identifying the configuration of the
        !> grid-points
        !
        !>@param bf_newgrdpt_coords0
        !> coordinates of the new grid-point in the reference
        !> frame of the buffer layer at t-dt
        !
        !>@param procedure_type
        !> integer identifying the procedure used to compute
        !> the new grid point
        !
        !>@param gradient_type
        !> integer identifying the gradient used to compute
        !> the new grid point
        !
        !>@param data_needed_bounds0
        !> identification of the grid-points needed around the
        !> new grid point to compute its values
        !
        !>@param tmp_array_needed
        !> logical determining whether a temporary array will be
        !> needed to gather the data required to compute the new
        !> grid point
        !
        !>@return grdpts_available
        !> logical determining whether there are enough grid points
        !> to compute the new grid point
        !--------------------------------------------------------------
        subroutine are_grdpts_available_to_get_newgrdpt_data(
     $     bf_grdpts_id0,
     $     bf_newgrdpt_coords0,
     $     procedure_type,
     $     gradient_type,
     $     data_needed_bounds0,
     $     tmp_array_needed,
     $     grdpts_available)

          implicit none
          
          integer       , dimension(:,:), intent(in)    :: bf_grdpts_id0
          integer(ikind), dimension(2)  , intent(in)    :: bf_newgrdpt_coords0
          integer                       , intent(in)    :: procedure_type
          integer                       , intent(in)    :: gradient_type
          integer       , dimension(2,2), intent(inout) :: data_needed_bounds0
          logical                       , intent(out)   :: tmp_array_needed
          logical                       , intent(out)   :: grdpts_available


          integer(ikind)                              :: size_x
          integer(ikind)                              :: size_y
          integer                                     :: nb_bounds
          integer       , dimension(2,2,2)            :: bounds
          integer                                     :: k
          logical                                     :: all_grdpts_exists

          
          !> estimate the bounds for the grid points to be checked
          call get_newgrdpt_verification_bounds(
     $         procedure_type,
     $         gradient_type,
     $         nb_bounds,
     $         bounds)


          !> update the bounds for the data needed
          do k=1, nb_bounds
             data_needed_bounds0(1,1) = min(data_needed_bounds0(1,1),bounds(1,1,k))
             data_needed_bounds0(1,2) = max(data_needed_bounds0(1,2),bounds(1,2,k))
             data_needed_bounds0(2,1) = min(data_needed_bounds0(2,1),bounds(2,1,k))
             data_needed_bounds0(2,2) = max(data_needed_bounds0(2,2),bounds(2,2,k))
          end do


          !> determine whether an intermediate array is needed to
          !> evaluate the data available
          tmp_array_needed = .false.

          size_x = size(bf_grdpts_id0,1)
          size_y = size(bf_grdpts_id0,2)


          do k=1, nb_bounds
             tmp_array_needed = tmp_array_needed.or.(
     $            ((bf_newgrdpt_coords0(1)+bounds(1,1,k)).lt.1).or. 
     $            ((bf_newgrdpt_coords0(1)+bounds(1,2,k)).gt.size_x).or.
     $            ((bf_newgrdpt_coords0(2)+bounds(2,1,k)).lt.1).or. 
     $            ((bf_newgrdpt_coords0(2)+bounds(2,2,k)).gt.size_y)
     $            )

          end do


          !> if no tmp_array is needed, the data needed for the new grdpt
          !> are directly checked on grdpts_id
          if(.not.tmp_array_needed) then

             all_grdpts_exists = .true.

             do k=1,nb_bounds

                if(.not.are_grdpts_available(
     $               bf_grdpts_id0,
     $               reshape((/
     $                  bf_newgrdpt_coords0(1)+bounds(1,1,k),
     $                  bf_newgrdpt_coords0(2)+bounds(2,1,k),
     $                  bf_newgrdpt_coords0(1)+bounds(1,2,k),
     $                  bf_newgrdpt_coords0(2)+bounds(2,2,k)/),
     $                  (/2,2/))
     $               )) then
                   
                   all_grdpts_exists = .false.

                end if

             end do

             if(.not.all_grdpts_exists) then

                print '(''bf_newgrdpt_dispatch_module'')'
                print '(''are_grdpts_available_to_get_newgrdpt_data'')'
                print '(''the data needed by the newgrdpt are in the buffer layer'')'
                print '(''but it does not seem they are all available'')'
                stop ''

             else
                grdpts_available = .true.
             end if

          else
             grdpts_available = .false.
          end if

        end subroutine are_grdpts_available_to_get_newgrdpt_data


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine whether there are enough grid points to
        !> identify the procedure computing the new grid point
        !
        !> @date
        !> 13_03_2015 - initial version - J.L. Desmarais
        !
        !>@param bf_localization
        !> cardinal coordinate identifying the mainlayer to which the
        !> buffer layer belongs
        !
        !>@param size_x
        !> size along the x-direction of the buffer layer
        !
        !>@param size_y
        !> size along the y-direction of the buffer layer
        !
        !>@param bf_can_exchange_with_neighbor1
        !> logical identifying whether the buffer layer can exchange
        !> grid points with its neighbor1
        !
        !>@param bf_can_exchange_with_neighbor2
        !> logical identifying whether the buffer layer can exchange
        !> grid points with its neighbor2
        !
        !>@param bf_newgrdpt_coords
        !> local coordinates of the new grid point to be computed
        !
        !>@return grdpts_available
        !> logical determining whether there are enough grid points
        !> to identify the procedure computing the new grid point
        !--------------------------------------------------------------
        function are_grdpts_available_to_get_newgrdpt_proc(
     $       bf_localization,
     $       size_x, size_y,
     $       bf_can_exchange_with_neighbor1,
     $       bf_can_exchange_with_neighbor2,
     $       bf_newgrdpt_coords)
     $       result(grdpts_available)
        
          implicit none

          integer                     , intent(in) :: bf_localization
          integer                     , intent(in) :: size_x
          integer                     , intent(in) :: size_y
          logical                     , intent(in) :: bf_can_exchange_with_neighbor1
          logical                     , intent(in) :: bf_can_exchange_with_neighbor2
          integer(ikind), dimension(2), intent(in) :: bf_newgrdpt_coords
          logical                                  :: grdpts_available


          select case(bf_localization)
            case(N)

               !([i-1] or [i+1]) and (j-1)
               grdpts_available = 
     $              ( ((bf_newgrdpt_coords(1)-1).ge.1).or.
     $                ((bf_newgrdpt_coords(1)+1).le.size_x)).and.
     $              ( ((bf_newgrdpt_coords(2)-1).ge.1) )

            case(S)

               !([i-1] or [i+1]) and (j+1)
               grdpts_available = 
     $              ( ((bf_newgrdpt_coords(1)-1).ge.1).or.
     $                ((bf_newgrdpt_coords(1)+1).le.size_x)).and.
     $              ( ((bf_newgrdpt_coords(2)+1).le.size_y) )
               
            case(E,W)

               if(bf_localization.eq.E) then
                  grdpts_available = (bf_newgrdpt_coords(1)-1).ge.1
               else
                  grdpts_available = (bf_newgrdpt_coords(1)+1).le.size_x
               end if

               if(bf_can_exchange_with_neighbor1) then

               !([i+1]) and (j-1,j+1)
                  if(bf_can_exchange_with_neighbor2) then
                     grdpts_available = 
     $                    grdpts_available.and.
     $                    ( (bf_newgrdpt_coords(2)-1).ge.1 ).and.
     $                    ( (bf_newgrdpt_coords(2)+1).le.size_y )

               !([i+1]) and (j-1)
                  else
                     grdpts_available = 
     $                    grdpts_available.and.
     $                    ( (bf_newgrdpt_coords(2)-1).ge.1 )
                  end if

               else

               !([i+1]) and (j+1)
                  if(bf_can_exchange_with_neighbor2) then
                     grdpts_available = 
     $                    grdpts_available.and.
     $                    ( (bf_newgrdpt_coords(2)+1).le.size_y )

               !([i+1]) and ((j-1) or (j+1))
                  else
                     grdpts_available = 
     $                    grdpts_available.and.(
     $                    ( (bf_newgrdpt_coords(2)-1).ge.1 ).or.
     $                    ( (bf_newgrdpt_coords(2)+1).le.size_y ))

                  end if                     

               end if

            case default
               call error_mainlayer_id(
     $              'bf_compute_newgrdpt_class',
     $              'are_grdpts_available_to_get_newgrdpt_proc',
     $              bf_localization)

          end select

        end function are_grdpts_available_to_get_newgrdpt_proc

      end module bf_newgrdpt_dispatch_module

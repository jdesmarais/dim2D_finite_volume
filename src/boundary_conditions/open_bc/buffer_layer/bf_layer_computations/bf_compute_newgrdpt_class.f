      !> @file
      !> extension of bf_compute to compute the new grid points
      !> belonging to the buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> extension of bf_compute to compute the new grid points
      !> belonging to the buffer layer
      !
      !> @date
      !> 03_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_compute_newgrdpt_class

        use bf_compute_class, only :
     $       bf_compute

        use bf_layer_extract_module, only :
     $       get_grdpts_id_from_interior,
     $       get_grdpts_id_from_bf_layer

        use bf_newgrdpt_procedure_module, only :
     $       get_newgrdpt_procedure

        use bf_newgrdpt_verification_module, only :
     $       verify_data_for_newgrdpt

        use bf_newgrdpt_class, only : 
     $       bf_newgrdpt

        use parameters_bf_layer, only :
     $       BF_SUCCESS,
     $       
     $       NEWGRDPT_PROC_NBGRDPTS_ERROR,
     $       NEWGRDPT_DATA_NBGRDPTS_ERROR,
     $       NEWGRDPT_SUCCESS

        use parameters_input, only :
     $       ne

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
        !> @param get_data_for_newgrdpt
        !> extract the data needed for the computation of the new
        !> gridpoint (from the previous time step and the grdpts_id)
        !
        !> @param compute_newgrdpt
        !> attempt to compute the new grid point
        !---------------------------------------------------------------
        type, extends(bf_compute) :: bf_compute_newgrdpt

          contains

          procedure, pass :: compute_newgrdpt

        end type bf_compute_newgrdpt

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
        function compute_newgrdpt(
     $       this,
     $       p_model, t, dt,
     $       newgrdpt_coords1,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes0,
     $       interior_nodes1,
     $       bf_localization,
     $       bf_alignment1,
     $       bf_nodes1,
     $       can_exchange_with_neighbor1,
     $       can_exchange_with_neighbor2,
     $       nb_procedures,
     $       procedure_type,
     $       gradient_type,
     $       grdpts_available,
     $       tmp_array_needed,
     $       data_needed_bounds0)
     $     result(ierror)

          implicit none

          class(bf_compute_newgrdpt)         , intent(in)  :: this
          type(pmodel_eq)                    , intent(in)  :: p_model
          real(rkind)                        , intent(in)  :: t
          real(rkind)                        , intent(in)  :: dt
          integer(ikind), dimension(2)       , intent(in)  :: newgrdpt_coords1
          real(rkind)   , dimension(nx)      , intent(in)  :: interior_x_map
          real(rkind)   , dimension(ny)      , intent(in)  :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne), intent(in)  :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne), intent(in)  :: interior_nodes1
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


        function get_newgrdpt_procedure_from_grdpts_id0(
     $     bf_localization,
     $     grdpts_id0,
     $     newgrdpt_coords0,
     $     can_exchange_with_neighbor1,
     $     can_exchange_with_neighbor2,
     $     nb_procedures,
     $     procedure_type,
     $     gradient_type)
     $     result(ierror)

          implicit none

          integer                       , intent(in)  :: bf_localization
          integer       , dimension(:,:), intent(in)  :: grdpts_id0
          integer(ikind), dimension(2)  , intent(in)  :: newgrdpt_coords0
          logical                       , intent(in)  :: can_exchange_with_neighbor1
          logical                       , intent(in)  :: can_exchange_with_neighbor2
          integer                       , intent(out) :: nb_procedures
          integer       , dimension(4)  , intent(out) :: procedure_type
          integer       , dimension(4)  , intent(out) :: gradient_type
          logical                                     :: ierror

          logical                 :: ierror_proc
          logical                 :: grdpts_available

          integer       , dimension(3,3) :: grdpts_id0_tmp
          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2,2) :: gen_coords
          integer                        :: i,j
          integer                        :: size_x
          integer                        :: size_y
          integer                        :: i_recv
          integer                        :: i_send
          integer                        :: j_recv
          integer                        :: j_send

          nb_procedures = 0
          procedure_type(1) = no_bc_procedure_type
          gradient_type(1)  = no_gradient_type


          !> if the grid points [i-1,i+1]x[j-1,j+1] are 
          !> available around the new gridpoint [i,j], the
          !> procedure needed for the new grdpt can be
          !> immediately determined
          !------------------------------------------------------------
          if(
     $         (newgrdpt_coords0(1).gt.1).and.
     $         (newgrdpt_coords0(1).lt.size(grdpts_id0,1)).and.
     $         (newgrdpt_coords0(2).gt.1).and.
     $         (newgrdpt_coords0(2).lt.size(grdpts_id0,2)) ) then

             ierror_proc = get_newgrdpt_procedure(
     $            newgrdpt_coords0(1),
     $            newgrdpt_coords0(2),
     $            grdpts_id0,
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
          !> enough to estimate the procedure
          !------------------------------------------------------------
          else
             
             !> determine whether there are enough grdpts_id to determine
             !> the procedure for the new gridpoint
             grdpts_available = are_grdpts_available_to_get_newgrdpt_proc(
     $            bf_localization,
     $            size(newgrdpt_coords0,1),
     $            size(newgrdpt_coords0,2),
     $            bf_grdpts_id0)


             !> if there are enough grid points to evaluate the procedure
             !> a temporary array is created using the grdpts_id from the
             !> interior domain and the procedure is evaluated
             if(grdpts_available) then

                !> determine the general coordinates for the
                !> corners of the grdpts_id extracted
                match_table = get_bf_layer_match_table(bf_alignment0)

                gen_coords(1,1) = newgrdpt_coords0(1)-1 + match_table(1)
                gen_coords(1,2) = newgrdpt_coords0(1)+1 + match_table(1)
                gen_coords(2,1) = newgrdpt_coords0(2)-1 + match_table(2)
                gen_coords(2,2) = newgrdpt_coords0(2)+1 + match_table(2)
                

                !> extract the grdpts_id to determine the new grdpt
                !> procedure

                !> extract the grdpts_id from the interior domain
                call get_grdpts_id_from_interior(
     $               grdpts_id0_tmp,
     $               gen_coords)

                !> extract the grdpts_id from the grdpts_id0
                call extract_grdpts_id_from_bf_layer(
     $               tmp_grdpts_id0,
     $               gen_coords,
     $               bf_alignment0,
     $               this%grdpts_id_tmp)

                !> determine the procedure for the computation
                !> of the new grid point
                ierror_proc = get_newgrdpt_procedure(
     $               2,2,
     $               tmp_grdpts_id0,
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

        end function get_newgrdpt_procedure_from_grdpts_id0



        function compute_newgrdpt_from_bf_layer(
     $     p_model,
     $     t,
     $     dt,
     $     bf_alignment0,
     $     grdpts_id0,
     $     bf_x_map0,
     $     bf_y_map0,
     $     bf_nodes0,
     $     bf_alignment1,
     $     bf_x_map1,
     $     bf_y_map1,
     $     bf_nodes1,
     $     newgrdpt_coords0,
     $     newgrdpt_coords1,
     $     nb_procedures,
     $     procedure_type,
     $     gradient_type,
     $     grdpts_available,
     $     tmp_array_needed,
     $     data_needed_bounds0)
     $     result(ierror)

          implicit none

          type(pmodel_eq)                    , intent(in) :: p_model
          real(rkind)                        , intent(in) :: t
          real(rkind)                        , intent(in) :: dt
          integer(ikind), dimension(2,2)     , intent(in) :: bf_alignment0
          integer       , dimension(:,:)     , intent(in) :: grdpts_id0
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map0
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map0
          real(rkind)   , dimension(:)       , intent(in) :: bf_nodes0
          integer(ikind), dimension(2,2)     , intent(in) :: bf_alignment1
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map1
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map1
          real(rkind)   , dimension(:)       , intent(in) :: bf_nodes1
          integer(ikind), dimension(2)       , intent(in) :: newgrdpt_coords0
          integer(ikind), dimension(2)       , intent(in) :: newgrdpt_coords1
          integer                            , intent(in) :: nb_procedures
          integer       , dimension(4)       , intent(in) :: procedure_type
          integer       , dimension(4)       , intent(in) :: gradient_type
          logical       , dimension(4)       , intent(out):: grdpts_available
          logical       , dimension(4)       , intent(out):: tmp_array_needed
          integer(ikind), dimension(2,2)     , intent(out):: data_needed_bounds0
          integer                                         :: ierror


          logical               :: all_grdpts_available
          integer               :: k
          type(bf_newgrdpt)     :: bf_newgrdpt_used


          !> for each procedure_type verify that there are
          !> enough grid points available          
          do k=1, nb_procedures

             call are_grdpts_available_to_get_newgrdpt_data(
     $            bf_localization,
     $            bf_grdpts_id0,
     $            bf_newgrdpt_coords0,
     $            procedure_type,
     $            gradient_type,
     $            data_needed_bounds0,
     $            tmp_array_needed(k),
     $            grdpts_available(k))

             all_grdpts_available = all_grdpts_available.and.grdpts_available(k)

          end do

          
          !> if all grid points are available to compute the
          !> new grid point, the new grdpt is directly computed
          if(all_grdpts_available) then

             bf_nodes1(
     $            newgrdpt_coords1(1),
     $            newgrdpt_coords1(2),
     $            :) =
     $       bf_newgrdpt_used%compute_newgrdpt(
     $            p_model, t, dt,
     $            bf_alignment0, bf_x_map0, bf_y_map0, bf_nodes0,
     $            bf_alignment1, bf_x_map1, bf_y_map1, bf_nodes1,
     $            newgrdpt_coords1(1),newgrdpt_coords1(2),
     $            nb_procedures, procedure_type, gradient_type)

             ierror = BF_SUCCESS

          !> otherwise we will need to create a temporary array to compute
          !> the new grdpt: if only one temporary array is needed,
          !> the overcost of having one big array for all data is negligeable
          !> compared to using a temporary array from the interior doamin +
          !> the buffer layer
          end if

             ierror = .not.BF_SUCCESS

          end if

        end function compute_newgrdpt_from_bf_layer


        !> determine whether the data are available to compute
        !> the new grid point with the procedure of type procedure_type
        !> using the gradient gradient_type
        subroutine are_grdpts_available_to_get_newgrdpt_data(
     $     bf_localization,
     $     bf_grdpts_id0,
     $     bf_newgrdpt_coords0,
     $     procedure_type,
     $     gradient_type,
     $     data_needed_bounds0,
     $     tmp_array_needed,
     $     grdpts_available)

          implicit none
          
          integer                       , intent(in)    :: bf_localization
          integer       , dimension(2,2), intent(in)    :: bf_grdpts_id0
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
          logical                                     :: tmp_array_needed
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
             tmp_array_needed = tmp_array_needed.and.(
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
     $                  bf_newgrdpt_coords0(2)+bounds(2,1,k)/),
     $                  (/2,2/))
     $               )) then
                   
                   all_grdpts_exists = .false.

                end if

             end do

             if(.not.all_grdpts_exists) then

                print '(''bf_compute_newgrdpt_class'')'
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

        end function are_grdpts_available_to_get_newgrdpt_data

      end module bf_compute_newgrdpt_class

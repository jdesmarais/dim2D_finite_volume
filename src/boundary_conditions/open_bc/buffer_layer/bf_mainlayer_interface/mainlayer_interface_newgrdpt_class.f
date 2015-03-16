      !> @file
      !> mainlayer_interface_dyn enhanced with procedures enabling
      !> to compute new grid points
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> mainlayer_interface_dyn enhanced with procedures enabling
      !> to compute new grid points
      !
      !> @date
      ! 14_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mainlayer_interface_newgrdpt_class

        use bf_layer_extract_module, only :
     $       get_map_from_interior,
     $       get_grdpts_id_from_interior,
     $       get_grdpts_id_from_bf_layer,
     $       get_nodes_from_interior

        use bf_newgrdpt_class, only :
     $       bf_newgrdpt

        use bf_newgrdpt_procedure_module, only :
     $       get_newgrdpt_procedure

        use bf_newgrdpt_verification_module, only :
     $       are_grdpts_available,
     $       get_newgrdpt_verification_bounds

        use bf_sublayer_class, only :
     $       bf_sublayer

        use mainlayer_interface_dyn_class, only :
     $       mainlayer_interface_dyn

        use parameters_bf_layer, only :
     $       NEWGRDPT_NO_ERROR,
     $       NEWGRDPT_NO_PREVIOUS_DATA_ERROR,
     $       NEWGRDPT_PROC_NBGRDPTS_ERROR,
     $       NEWGRDPT_DATA_NBGRDPTS_ERROR,
     $       BF_SUCCESS

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none


        private
        public :: mainlayer_interface_newgrdpt


        !> @class mainlayer_interface_newgrdpt
        !> mainlayer_interface_dyn enhanced with procedures enabling
        !> to compute new grid points
        !
        !> @param compute_newgrdpt
        !> compute the new grid point identified by its general
        !> coordinate
        !
        !> @param are_grdpts_available_to_compute_newgrdpt
        !> check whether grid-points are available to compute the
        !> new grid-point
        !
        !> @param collect_data_to_compute_newgrdpt
        !> collect the data needed to compute the new grid-point
        !------------------------------------------------------------
        type, extends(mainlayer_interface_dyn) :: mainlayer_interface_newgrdpt

          contains

          procedure,   pass :: compute_newgrdpt

          procedure, nopass :: are_grdpts_available_to_compute_newgrdpt
          procedure,   pass :: collect_data_to_compute_newgrdpt

        end type mainlayer_interface_newgrdpt


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid point identified by its coordinates in the
        !> general reference frame: if it is possible using the data available
        !> in the buffer layer or the data in the interior domain and
        !> its neighbors
        !
        !> @date
        !> 14_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> mainlayer_interface_newgrdpt encapsulating the links to the
        !> buffer layers at the interface b/w main layers
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
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes0
        !> nodes of the interior domain at t-dt
        !
        !>@param interior_nodes1
        !> nodes of the interior domain at t
        !
        !>@param bf_sublayer_ptr
        !> pointer to the buffer layer where the
        !> new grid point is determined
        !
        !>@param gen_newgrdpt_coords
        !> coordinates of the new grid-point computed in the
        !> general reference frame
        !
        !>@return ierror
        !> integer identifying whether the computation was successful
        !--------------------------------------------------------------
        function compute_newgrdpt(
     $       this,
     $       p_model,
     $       t,
     $       dt,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes0,
     $       interior_nodes1,
     $       bf_sublayer_ptr,
     $       gen_newgrdpt_coords)
     $       result(ierror)

          implicit none

          class(mainlayer_interface_newgrdpt)   , intent(inout) :: this
          type(pmodel_eq)                       , intent(in)    :: p_model
          real(rkind)                           , intent(in)    :: t
          real(rkind)                           , intent(in)    :: dt
          real(rkind)      , dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind)      , dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind)      , dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind)      , dimension(nx,ny,ne), intent(in)    :: interior_nodes1
          type(bf_sublayer), pointer            , intent(inout) :: bf_sublayer_ptr
          integer(ikind)   , dimension(2)       , intent(in)    :: gen_newgrdpt_coords
          integer                                               :: ierror

          integer(ikind), dimension(2)   :: bf_newgrdpt_coords1

          integer                        :: ierror_bf
          logical                        :: ierror_proc

          integer                        :: nb_procedures
          integer       , dimension(4)   :: procedure_type
          integer       , dimension(4)   :: gradient_type
          logical       , dimension(4)   :: grdpts_available
          integer(ikind), dimension(2,2) :: data_needed_bounds0

          integer       , parameter                       :: size_tmp=2*(bc_size+1)+1
          integer       , dimension(size_tmp,size_tmp)    :: grdpts_id_tmp
          real(rkind)   , dimension(size_tmp)             :: x_map_tmp
          real(rkind)   , dimension(size_tmp)             :: y_map_tmp
          real(rkind)   , dimension(size_tmp,size_tmp,ne) :: nodes0_tmp
          real(rkind)   , dimension(size_tmp,size_tmp,ne) :: nodes1_tmp
          integer(ikind), dimension(2)                    :: newgrdpt_coords_tmp
          integer(ikind), dimension(2,2)                  :: align_tmp

          integer                    :: k
          logical                    :: all_grdpts_available
          real(rkind), dimension(ne) :: new_grdpt

          type(bf_newgrdpt) :: bf_newgrdpt_used

          
          !------------------------------------------------------------
          !1) ask the buffer layer to which the new grid point is
          !   belonging whether it is possible to compute the new grid
          !   point using only the data available in this buffer layer
          !------------------------------------------------------------
          !1.1) compute the coordinates of the new grid point in the
          !     local reference frame of the buffer layer
          !------------------------------------------------------------
          bf_newgrdpt_coords1 = bf_sublayer_ptr%get_local_coord(gen_newgrdpt_coords)
          
          !------------------------------------------------------------
          !1.2) ask the buffer layer to compute the new grid point
          !------------------------------------------------------------
          ierror_bf = bf_sublayer_ptr%compute_newgrdpt(
     $         p_model,
     $         t,
     $         dt,
     $         bf_newgrdpt_coords1,
     $         nb_procedures,
     $         procedure_type,
     $         gradient_type,
     $         grdpts_available,
     $         data_needed_bounds0)


          !2) check whether the buffer layer was able to compute the
          !   new grid point by itself
          !------------------------------------------------------------
          ! - NEWGRDPT_NO_ERROR:
          !      the buffer layer could determine the procedure needed
          !      to compute the new grid-point and all grid points
          !      were available to compute the new grid point
          !
          ! - NEWGRDPT_NO_PREVIOUS_DATA_ERROR:
          !      the buffer layer does not contain the data at the
          !      previous time step and so it was impossible to
          !      compute the new grid-point
          !
          ! - NEWGRDPT_PROC_NBGRDPTS_ERROR:
          !      the buffer layer did not contain enough grid points
          !      to determine the procedure needed to compute the new
          !      grid-point and so it was impossible to compute
          !      the new grid-point
          !
          ! - NEWGRDPT_DATA_NBGRDPTS_ERROR:
          !      the buffer layer was able to determine the procedure
          !      needed to compute the new grid-point but the grid 
          !      points needed to compute the new grid point were not
          !      all available
          !------------------------------------------------------------
          if(
     $         (ierror_bf.eq.NEWGRDPT_NO_PREVIOUS_DATA_ERROR).or.
     $         (ierror_bf.eq.NEWGRDPT_PROC_NBGRDPTS_ERROR).or.
     $         (ierror_bf.eq.NEWGRDPT_DATA_NBGRDPTS_ERROR)) then


             ! if it was impossible to determine the procedure,
             ! it is not possible to estimate the restricted
             ! amount of grid points to compute the new grid
             ! point so the amount of data gathered corresponds
             ! to the maximum potentially [2*(bc_size+1)]
             ! around the new grid point
             if(
     $            (ierror_bf.eq.NEWGRDPT_NO_PREVIOUS_DATA_ERROR).or.
     $            (ierror_bf.eq.NEWGRDPT_PROC_NBGRDPTS_ERROR)) then

                data_needed_bounds0 = reshape((/
     $               -(bc_size+1), -(bc_size+1),
     $                (bc_size+1),  (bc_size+1)/),
     $               (/2,2/))

             ! for ierror.eq.NEWGRDPT_DATA_NBGRDPTS_ERROR, the
             ! data_needed_bounds has already been initialized
             ! by bf_sublayer_ptr%compute_newgrdpt
             end if

             
             ! we now gather the data needed to compute the new
             ! grid point: grdpts_id0, x_map0, y_map0, x_map1,
             ! y_map1, nodes0, nodes1
             call collect_data_to_compute_newgrdpt(
     $            this,
     $            bf_sublayer_ptr,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            gen_newgrdpt_coords,
     $            data_needed_bounds0,
     $            grdpts_id_tmp,
     $            x_map_tmp,
     $            y_map_tmp,
     $            nodes0_tmp,
     $            nodes1_tmp)


             ! in this frame, the coordinates of the new grdpt
             ! to be computed are identified by:
             ! [-data_needed_bounds0(1,1)+1,-data_needed_bounds0(2,1)+1]
             newgrdpt_coords_tmp = [
     $            -data_needed_bounds0(1,1)+1,
     $            -data_needed_bounds0(2,1)+1]


             ! if it was impossible to determine the procedure,
             ! we will now determine its characteristics
             if( 
     $            (ierror_bf.eq.NEWGRDPT_NO_PREVIOUS_DATA_ERROR).or.
     $            (ierror_bf.eq.NEWGRDPT_PROC_NBGRDPTS_ERROR)) then

                ierror_proc = get_newgrdpt_procedure(
     $               newgrdpt_coords_tmp(1),
     $               newgrdpt_coords_tmp(2),
     $               grdpts_id_tmp,
     $               nb_procedures,
     $               procedure_type,
     $               gradient_type)

                do k=1,nb_procedures
                   grdpts_available(k) = .false.
                end do

                if(ierror_proc.neqv.BF_SUCCESS) then
                   print '(''mainlayer_interface_newgrdpt'')'
                   print '(''compute_newgrdpt'')'
                   print '(''unable to find the procedure characteristics'')'
                   stop ''
                end if

             ! for ierror_bf.eq.NEWGRDPT_DATA_NBGRDPTS_ERROR
             ! the characteristics of the procedure have been
             ! initialized by bf_sublayer_ptr%compute_newgrdpt
             end if

             
             ! for each procedure needed to compute the new
             ! grid-point, verify that the grid points are
             ! effectively available
             all_grdpts_available = are_grdpts_available_to_compute_newgrdpt(
     $            grdpts_available,
     $            nb_procedures,
     $            procedure_type,
     $            gradient_type,
     $            newgrdpt_coords_tmp,
     $            grdpts_id_tmp)

             if(.not.all_grdpts_available) then
                print '(''mainlayer_interface_newgrdpt'')'
                print '(''compute_newgrdpt'')'
                print '(''error when checking the grdpts for newgrdpt'')'
                stop ''
             end if


             ! if all the grdpts are available, then the new
             ! grid point is computed
             align_tmp(1,1) = gen_newgrdpt_coords(1)+data_needed_bounds0(1,1)+bc_size
             align_tmp(1,2) = align_tmp(1,1)+2
             align_tmp(2,1) = gen_newgrdpt_coords(2)+data_needed_bounds0(2,1)+bc_size
             align_tmp(2,2) = align_tmp(2,1)+2

             new_grdpt = bf_newgrdpt_used%compute_newgrdpt(
     $            p_model, t, dt,
     $            align_tmp, x_map_tmp, y_map_tmp, nodes0_tmp,
     $            align_tmp, x_map_tmp, y_map_tmp, nodes1_tmp,
     $            newgrdpt_coords_tmp(1),
     $            newgrdpt_coords_tmp(2),
     $            nb_procedures, procedure_type, gradient_type)


             ! this new grdpt is then set in the buffer layer
             ! to which it belongs
             do k=1,ne
                call bf_sublayer_ptr%set_nodes_pt(
     $               bf_newgrdpt_coords1(1),
     $               bf_newgrdpt_coords1(2),
     $               k,
     $               new_grdpt(k))
             end do

             ierror = BF_SUCCESS


          else

             ierror = BF_SUCCESS

          end if

        end function compute_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine whether all grid-points needed to compute
        !> the new grid-point are available
        !
        !> @date
        !> 16_03_2015 - initial version - J.L. Desmarais
        !
        !> @param grdpts_available
        !> logical array determining the procedures already checked
        !
        !> @param nb_procedures
        !> number of porcedures needed to compute the new grid-point
        !
        !> @param procedure_type
        !> procedures needed to compute the new grid-point
        !
        !> @param gradient_type
        !> gradients needed to compute the new grid-point
        !
        !> @param interior_nodes0
        !> nodes of the interior domain at t-dt
        !
        !> @param tmp_newgrdpt_coords
        !> coordinates of the new grid point in the frame of the
        !> tmp_grdpts_id array
        !
        !> @param tmp_grdpts_id
        !> temporary array gathering the configuration of the grid
        !> points needed to compute the new grid-point
        !
        !> @return all_grdpts_available
        !> logical identifying whether all grid-points are available
        !--------------------------------------------------------------
        function are_grdpts_available_to_compute_newgrdpt(
     $     grdpts_available,
     $     nb_procedures,
     $     procedure_type,
     $     gradient_type,
     $     tmp_newgrdpt_coords,
     $     tmp_grdpts_id)
     $     result(all_grdpts_available)

          implicit none

          logical, dimension(4)  , intent(in) :: grdpts_available
          integer                , intent(in) :: nb_procedures
          integer, dimension(4)  , intent(in) :: procedure_type
          integer, dimension(4)  , intent(in) :: gradient_type
          integer, dimension(2)  , intent(in) :: tmp_newgrdpt_coords
          integer, dimension(:,:), intent(in) :: tmp_grdpts_id
          logical                             :: all_grdpts_available

          integer                   :: k,m
          integer                   :: nb_bounds
          integer, dimension(2,2,2) :: bounds


          all_grdpts_available = .true.


          ! loop over the procedure that are applied to compute
          ! the new grid-point
          do k=1, nb_procedures

          !1) for each procedure, determine the grdpts that should
          !   be checked
             call get_newgrdpt_verification_bounds(
     $            procedure_type(k),
     $            gradient_type(k),
     $            nb_bounds,
     $            bounds)

          !2) for each procedure, check if the grdpts asked to
          !   compute the new grid point are available
             if(.not.grdpts_available(k)) then

                do m=1, nb_bounds
                   if(.not.are_grdpts_available(
     $                  tmp_grdpts_id,
     $                  reshape((/
     $                  tmp_newgrdpt_coords(1)+bounds(1,1,m),
     $                  tmp_newgrdpt_coords(2)+bounds(2,1,m),
     $                  tmp_newgrdpt_coords(1)+bounds(1,2,m),
     $                  tmp_newgrdpt_coords(2)+bounds(2,2,m)/),
     $                  (/2,2/))
     $                  )) then
                      all_grdpts_available = .false.
                   end if
                end do

             end if

          end do

        end function are_grdpts_available_to_compute_newgrdpt

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> collect the data to compute the new grid-point
        !
        !> @date
        !> 16_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object with references to the buffer layers at the interface
        !> b/w the main layers
        !
        !> @param bf_sublayer_ptr
        !> pointer to the buffer layer whose new grid point is computed
        !
        !> @param interior_x_map0
        !> x-coordinates of the interior domain
        !
        !> @param interior_y_map0
        !> y-coordinates of the interior domain
        !
        !> @param interior_nodes0
        !> nodes of the interior domain at t-dt
        !
        !> @param interior_nodes1
        !> nodes of the interior domain at t
        !
        !> @param gen_newgrdpt_coords
        !> coordinates of the new grid point for which data are
        !> collected in the general reference frame
        !
        !> @param tmp_grdpts_id
        !> temporary array gathering the grdpts_id needed to compute
        !> the new grid-point
        !
        !> @param tmp_x_map
        !> temporary array gathering the x-coordinates needed to compute
        !> the new grid-point
        !
        !> @param tmp_y_map
        !> temporary array gathering the y-coordinates needed to compute
        !> the new grid-point
        !
        !> @param tmp_nodes0
        !> temporary array gathering the nodes needed to compute the
        !> new grid-point at t-dt
        !
        !> @param tmp_nodes1
        !> temporary array gathering the nodes needed to compute the
        !> new grid-point at t
        !--------------------------------------------------------------
        subroutine collect_data_to_compute_newgrdpt(
     $     this,
     $     bf_sublayer_ptr,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1,
     $     gen_newgrdpt_coords,
     $     data_needed_bounds0,
     $     tmp_grdpts_id,
     $     tmp_x_map,
     $     tmp_y_map,
     $     tmp_nodes0,
     $     tmp_nodes1)

          implicit none

          class(mainlayer_interface_newgrdpt), intent(in)  :: this
          type(bf_sublayer), pointer         , intent(in)  :: bf_sublayer_ptr
          real(rkind), dimension(nx)         , intent(in)  :: interior_x_map
          real(rkind), dimension(ny)         , intent(in)  :: interior_y_map
          real(rkind), dimension(nx,ny,ne)   , intent(in)  :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)   , intent(in)  :: interior_nodes1
          integer(ikind), dimension(2)       , intent(in)  :: gen_newgrdpt_coords
          integer(ikind), dimension(2,2)     , intent(in)  :: data_needed_bounds0
          integer       , dimension(:,:)     , intent(out) :: tmp_grdpts_id
          real(rkind)   , dimension(:)       , intent(out) :: tmp_x_map
          real(rkind)   , dimension(:)       , intent(out) :: tmp_y_map
          real(rkind)   , dimension(:,:,:)   , intent(out) :: tmp_nodes0
          real(rkind)   , dimension(:,:,:)   , intent(out) :: tmp_nodes1


          integer(ikind)   , dimension(2,2) :: gen_coords
          integer(ikind)   , dimension(6)   :: extract_param
          type(bf_sublayer), pointer        :: bf_neighbor_ptr


          ! determination of the coordinates of the SW and NE borders
          ! of the data to be extracted in the general frame
          gen_coords(1,1) = gen_newgrdpt_coords(1)+data_needed_bounds0(1,1)
          gen_coords(1,2) = gen_newgrdpt_coords(1)+data_needed_bounds0(1,2)
          gen_coords(2,1) = gen_newgrdpt_coords(2)+data_needed_bounds0(2,1)
          gen_coords(2,2) = gen_newgrdpt_coords(2)+data_needed_bounds0(2,2)


          ! extraction of the coordinate maps
          !============================================================
          ! extraction of the x_map_tmp
          call get_map_from_interior(
     $         tmp_x_map,
     $         gen_coords(1,:),
     $         interior_x_map)


          ! extraction of the y_map_tmp
          call get_map_from_interior(
     $         tmp_y_map,
     $         gen_coords(2,:),
     $         interior_y_map)


          ! extraction from the interior domain
          !============================================================
          ! grdpts_id
          !------------------------------------------------------------
          call get_grdpts_id_from_interior(
     $         tmp_grdpts_id,
     $         gen_coords)

          ! nodes0
          !------------------------------------------------------------
          call get_nodes_from_interior(
     $         tmp_nodes0,
     $         gen_coords,
     $         interior_nodes0,
     $         extract_param_out=extract_param)

          ! nodes1
          !------------------------------------------------------------
          call get_nodes_from_interior(
     $         tmp_nodes1,
     $         gen_coords,
     $         interior_nodes1,
     $         extract_param_in=extract_param)


          ! extraction from the neighbor1 if any
          !============================================================
          if(bf_sublayer_ptr%can_exchange_with_neighbor1()) then
             
             bf_neighbor_ptr => this%get_neighbor_sublayer_ptr(
     $            bf_sublayer_ptr%get_localization(),1)

          ! grdpts_id
          !------------------------------------------------------------
             call bf_neighbor_ptr%extract_grdpts_id(
     $            tmp_grdpts_id,
     $            gen_coords,
     $            extract_param_out=extract_param,
     $            previous_step=.true.)


          ! nodes0
          !------------------------------------------------------------
             call bf_neighbor_ptr%extract_nodes(
     $            tmp_nodes0,
     $            gen_coords,
     $            extract_param_in=extract_param,
     $            previous_step=.true.)

          ! nodes1
          !------------------------------------------------------------
             call bf_neighbor_ptr%extract_nodes(
     $            tmp_nodes1,
     $            gen_coords,
     $            previous_step=.false.)

          end if


          ! extraction from the neighbor2 if any
          !============================================================
          if(bf_sublayer_ptr%can_exchange_with_neighbor2()) then
             
             bf_neighbor_ptr => this%get_neighbor_sublayer_ptr(
     $            bf_sublayer_ptr%get_localization(),2)

          ! grdpts_id
          !------------------------------------------------------------
             call bf_neighbor_ptr%extract_grdpts_id(
     $            tmp_grdpts_id,
     $            gen_coords,
     $            extract_param_out=extract_param,
     $            previous_step=.true.)

          ! nodes0
          !------------------------------------------------------------
             call bf_neighbor_ptr%extract_nodes(
     $            tmp_nodes0,
     $            gen_coords,
     $            extract_param_in=extract_param,
     $            previous_step=.true.)

          ! nodes1
          !------------------------------------------------------------
             call bf_neighbor_ptr%extract_nodes(
     $            tmp_nodes1,
     $            gen_coords,
     $            previous_step=.false.)

          end if


          ! extraction from the current buffer layer
          !============================================================
          ! grdpts_id
          !------------------------------------------------------------
          call bf_sublayer_ptr%extract_grdpts_id(
     $         tmp_grdpts_id,
     $         gen_coords,
     $         extract_param_out=extract_param,
     $         previous_step=.true.)

          ! nodes0
          !------------------------------------------------------------
          call bf_sublayer_ptr%extract_nodes(
     $         tmp_nodes0,
     $         gen_coords,
     $         extract_param_in=extract_param,
     $         previous_step=.true.)

          ! nodes1
          !------------------------------------------------------------
          call bf_sublayer_ptr%extract_nodes(
     $         tmp_nodes1,
     $         gen_coords,
     $         previous_step=.false.)


        end subroutine collect_data_to_compute_newgrdpt
        
      end module mainlayer_interface_newgrdpt_class

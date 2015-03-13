      program test_bf_newgrdpt_dispatch

        use bf_newgrdpt_dispatch_module, only :
     $     get_newgrdpt_procedure_from_bf_layer,
     $     compute_newgrdpt_from_bf_layer,
     $     are_grdpts_available_to_get_newgrdpt_data,
     $     are_grdpts_available_to_get_newgrdpt_proc

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated,
     $       is_int_vector_validated,
     $       is_int_matrix_validated,
     $       is_boolean_vector_validated

        use dim2d_parameters, only :
     $       cv_r

        use parameters_bf_layer, only :
     $       align_N, align_E, align_W,
     $       
     $       NE_corner_type,
     $       SW_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       
     $       gradient_I_type,
     $       gradient_xI_yLR0_type,
     $       
     $       BF_SUCCESS,
     $       no_pt

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only : 
     $       pmodel_eq

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled      = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_are_grdpts_available_to_get_newgrdpt_proc(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_are_grdpts_available_to_get_newgrdpt_proc: '',L1)', test_loc
        print '()'


        test_loc = test_are_grdpts_available_to_get_newgrdpt_data(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_are_grdpts_available_to_get_newgrdpt_data: '',L1)', test_loc
        print '()'


        test_loc = test_compute_newgrdpt_from_buffer_layer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_newgrdpt_from_buffer_layer: '',L1)', test_loc
        print '()'


        test_loc = test_get_newgrdpt_procedure_from_buffer_layer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_newgrdpt_procedure_from_buffer_layer: '',L1)', test_loc
        print '()'


        contains


        function test_get_newgrdpt_procedure_from_buffer_layer(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer :: k
          logical :: test_loc

          integer                        :: bf_localization
          integer(ikind), dimension(2,2) :: bf_alignment0
          integer       , dimension(5,6) :: bf_grdpts_id0
          integer(ikind), dimension(2)   :: bf_newgrdpt_coords0
          logical                        :: bf_can_exchange_with_neighbor1
          logical                        :: bf_can_exchange_with_neighbor2

          integer                        :: test_nb_procedures
          integer       , dimension(4)   :: test_procedure_type
          integer       , dimension(4)   :: test_gradient_type
          logical                        :: test_ierror

          integer                        :: nb_procedures
          integer       , dimension(4)   :: procedure_type
          integer       , dimension(4)   :: gradient_type
          logical                        :: ierror


          test_validated = .true.


          do k=1,3

             !input
             call get_param_test_newgrdpt_proc_bf_layer(
     $            k,
     $            bf_localization,
     $            bf_alignment0,
     $            bf_grdpts_id0,
     $            bf_newgrdpt_coords0,
     $            bf_can_exchange_with_neighbor1,
     $            bf_can_exchange_with_neighbor2,
     $            test_nb_procedures,
     $            test_procedure_type,
     $            test_gradient_type,
     $            test_ierror)

             !output
             ierror = get_newgrdpt_procedure_from_bf_layer(
     $            bf_localization,
     $            bf_alignment0,
     $            bf_grdpts_id0,
     $            bf_newgrdpt_coords0,
     $            bf_can_exchange_with_neighbor1,
     $            bf_can_exchange_with_neighbor2,
     $            nb_procedures,
     $            procedure_type,
     $            gradient_type)

             !validation
             test_loc = ierror.eqv.test_ierror
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''ierror failed'')'
             end if

             if(ierror.eqv.BF_SUCCESS) then

                test_loc = nb_procedures.eq.test_nb_procedures
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''nb_procedures failed'')'
                end if

                test_loc = is_int_vector_validated(
     $               procedure_type(1:nb_procedures),
     $               test_procedure_type(1:nb_procedures),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''procedure_type failed'')'
                end if
                
                test_loc = is_int_vector_validated(
     $               gradient_type(1:nb_procedures),
     $               test_gradient_type(1:nb_procedures),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''gradient_type failed'')'
                end if

             end if

          end do

        end function test_get_newgrdpt_procedure_from_buffer_layer


        subroutine get_param_test_newgrdpt_proc_bf_layer(
     $     test_id,
     $     bf_localization,
     $     bf_alignment0,
     $     bf_grdpts_id0,
     $     bf_newgrdpt_coords0,
     $     bf_can_exchange_with_neighbor1,
     $     bf_can_exchange_with_neighbor2,
     $     test_nb_procedures,
     $     test_procedure_type,
     $     test_gradient_type,
     $     test_ierror)

          implicit none

          integer                       , intent(in)  :: test_id
          integer                       , intent(out) :: bf_localization
          integer(ikind), dimension(2,2), intent(out) :: bf_alignment0
          integer       , dimension(:,:), intent(out) :: bf_grdpts_id0
          integer(ikind), dimension(2)  , intent(out) :: bf_newgrdpt_coords0
          logical                       , intent(out) :: bf_can_exchange_with_neighbor1
          logical                       , intent(out) :: bf_can_exchange_with_neighbor2
          integer                       , intent(out) :: test_nb_procedures
          integer       , dimension(4)  , intent(out) :: test_procedure_type
          integer       , dimension(4)  , intent(out) :: test_gradient_type
          logical                       , intent(out) :: test_ierror


          bf_localization = E

          bf_alignment0 = reshape((/
     $         align_E,align_N-2,align_E,align_N-1/),
     $         (/2,2/))

          bf_can_exchange_with_neighbor1 = .false.
          bf_can_exchange_with_neighbor2 = .false.

          test_nb_procedures     = 1
          test_procedure_type(1) = E_edge_type
          test_gradient_type(1)  = gradient_I_type


          select case(test_id)

            case(1)
               bf_grdpts_id0 = reshape((/
     $              1,2,3,0,0,
     $              1,2,3,0,0,
     $              1,2,3,0,0,
     $              1,2,3,0,0,
     $              1,2,3,0,0,
     $              1,3,3,0,0/),
     $              (/5,6/))

               bf_newgrdpt_coords0 = [4,4]

               test_ierror = BF_SUCCESS

            case(2)
               bf_grdpts_id0 = reshape((/
     $              2,3,3,3,3,
     $              2,2,2,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              2,2,2,2,3,
     $              2,3,3,3,3/),
     $              (/5,6/))

               bf_newgrdpt_coords0 = [6,4]

               test_ierror = BF_SUCCESS

            case(3)
               bf_grdpts_id0 = reshape((/
     $              2,3,3,3,3,
     $              2,2,2,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              2,2,2,2,3,
     $              2,3,3,3,3/),
     $              (/5,6/))

               bf_can_exchange_with_neighbor2 = .true.

               bf_newgrdpt_coords0 = [7,4]

               test_ierror = .not.BF_SUCCESS

          end select

        end subroutine get_param_test_newgrdpt_proc_bf_layer


        function test_compute_newgrdpt_from_buffer_layer(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: k
          logical :: test_loc

          type(pmodel_eq)                  :: p_model
          real(rkind)                      :: t
          real(rkind)                      :: dt
          integer(ikind), dimension(2,2)   :: bf_alignment0
          integer       , dimension(6,5)   :: bf_grdpts_id0
          real(rkind)   , dimension(6)     :: bf_x_map0
          real(rkind)   , dimension(5)     :: bf_y_map0
          real(rkind)   , dimension(6,5,ne):: bf_nodes0
          integer(ikind), dimension(2,2)   :: bf_alignment1
          real(rkind)   , dimension(8)     :: bf_x_map1
          real(rkind)   , dimension(6)     :: bf_y_map1
          real(rkind)   , dimension(8,6,ne):: bf_nodes1
          integer(ikind), dimension(2)     :: bf_newgrdpt_coords0
          integer(ikind), dimension(2)     :: bf_newgrdpt_coords1
          integer                          :: nb_procedures
          integer       , dimension(4)     :: procedure_type
          integer       , dimension(4)     :: gradient_type

          logical       , dimension(4)     :: test_grdpts_available
          integer(ikind), dimension(2,2)   :: test_data_needed_bounds0
          logical                          :: test_ierror
          real(rkind)   , dimension(ne)    :: test_newgrdpt

          logical       , dimension(4)     :: grdpts_available
          integer(ikind), dimension(2,2)   :: data_needed_bounds0
          logical                          :: ierror


          test_validated = .true.

          call p_model%initial_conditions%ini_far_field()


          do k=1,4

             !input
             call get_param_test_compute_from_bf_layer(
     $            k,
     $            t,dt,
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
     $            test_grdpts_available,
     $            test_data_needed_bounds0,
     $            test_ierror,
     $            test_newgrdpt)

             !output
             ierror = compute_newgrdpt_from_bf_layer(
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

             !validation
             !ierror
             test_loc = ierror.eqv.test_ierror
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''ierror failed'')'
             end if

             !grdpts_available
             test_loc = is_boolean_vector_validated(
     $            grdpts_available(1:nb_procedures),
     $            test_grdpts_available(1:nb_procedures),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''grdpts_available failed'')'
             end if

             !data_needed_bounds0
             test_loc = is_int_matrix_validated(
     $            data_needed_bounds0,
     $            test_data_needed_bounds0,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''data_needed_bounds failed'')'
             end if

             !newgrdpt_data
             if(ierror.eqv.BF_SUCCESS) then
                test_loc = is_real_vector_validated(
     $               bf_nodes1(bf_newgrdpt_coords1(1),bf_newgrdpt_coords1(2),:),
     $               test_newgrdpt,
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''newgrdpt_data failed'')'
                end if
             end if             

          end do

        end function test_compute_newgrdpt_from_buffer_layer


        subroutine get_param_test_compute_from_bf_layer(
     $     test_id,
     $     t,dt,
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
     $     test_grdpts_available,
     $     test_data_needed_bounds0,
     $     test_ierror,
     $     test_newgrdpt)

          implicit none

          
          integer                         , intent(in)  :: test_id

          real(rkind)                     , intent(out) :: t
          real(rkind)                     , intent(out) :: dt
          integer(ikind), dimension(2,2)  , intent(out) :: bf_alignment0
          integer       , dimension(:,:)  , intent(out) :: bf_grdpts_id0
          real(rkind)   , dimension(:)    , intent(out) :: bf_x_map0
          real(rkind)   , dimension(:)    , intent(out) :: bf_y_map0
          real(rkind)   , dimension(:,:,:), intent(out) :: bf_nodes0
          integer(ikind), dimension(2,2)  , intent(out) :: bf_alignment1
          real(rkind)   , dimension(:)    , intent(out) :: bf_x_map1
          real(rkind)   , dimension(:)    , intent(out) :: bf_y_map1
          real(rkind)   , dimension(:,:,:), intent(out) :: bf_nodes1
          integer(ikind), dimension(2)    , intent(out) :: bf_newgrdpt_coords0
          integer(ikind), dimension(2)    , intent(out) :: bf_newgrdpt_coords1
          integer                         , intent(out) :: nb_procedures
          integer       , dimension(4)    , intent(out) :: procedure_type
          integer       , dimension(4)    , intent(out) :: gradient_type

          logical       , dimension(4)    , intent(out) :: test_grdpts_available
          integer(ikind), dimension(2,2)  , intent(out) :: test_data_needed_bounds0
          logical                         , intent(out) :: test_ierror
          real(rkind)   , dimension(ne)   , intent(out) :: test_newgrdpt

          real(rkind)   , dimension(ne) :: test_newgrdpt_W
          real(rkind)   , dimension(ne) :: test_newgrdpt_NE

          integer(ikind) :: i,j
          integer        :: k

          t  = 0.0d0
          dt = 0.25d0

          bf_alignment0 = reshape((/align_W+2,align_N,align_W+3,align_N/),(/2,2/))
          bf_x_map0 = [0.5d0, 1.50d0, 2.50d0, 3.50d0, 4.50d0, 5.50d0]
          bf_y_map0 = [0.0d0, 0.25d0, 0.50d0, 0.75d0, 1.00d0]

          bf_alignment1 = reshape((/align_W+1,align_N,align_W+4, align_N+1/), (/2,2/))
          bf_x_map1 = [-0.5d0, 0.50d0, 1.50d0, 2.50d0, 3.50d0, 4.50d0, 5.50d0, 6.50d0]
          bf_y_map1 = [ 0.0d0, 0.25d0, 0.50d0, 0.75d0, 1.00d0, 1.25d0]

          
          nb_procedures     = 2
          procedure_type(1) = NE_corner_type
          procedure_type(2) = W_edge_type
          gradient_type(2)  = gradient_I_type

          test_data_needed_bounds0 = reshape((/-3,-3,2,1/),(/2,2/))

          select case(test_id)
            case(1)
               bf_newgrdpt_coords0 = [4,4]
               bf_newgrdpt_coords1 = [5,4]

               bf_grdpts_id0 = reshape((/
     $              ((no_pt,i=1,6),j=1,5)/),(/6,5/))

               bf_nodes0 = reshape((/
     $              (((-99.0d0,i=1,6),j=1,5),k=1,ne)/),
     $              (/6,5,ne/))

               bf_nodes1 = reshape((/
     $              (((-99.0d0,i=1,8),j=1,6),k=1,ne)/),
     $              (/8,6,ne/))

               !data for the NE corner
               bf_nodes0(1:3,1:3,:) = reshape((/
     $              1.48d0, 1.30d0, 1.35d0,
     $              1.26d0, 1.45d0, 1.40d0,
     $              1.46d0, 1.27d0, 1.47d0,
     $              
     $              0.128d0, 0.127d0, 0.142d0,
     $              1.138d0, 0.148d0, 0.132d0,
     $              0.146d0, 0.143d0, 0.145d0,
     $              
     $              0.0050d0, 0.020d0, 0.060d0,
     $              0.0025d0, 0.001d0, 0.015d0, 
     $              0.0100d0, 0.002d0, 0.050d0,
     $              
     $              4.88d0, 4.870d0, 4.855d0,
     $              4.85d0, 4.865d0, 4.845d0,
     $              4.89d0, 4.870d0, 4.860d0/),
     $              (/3,3,ne/))

               bf_nodes1(2:4,1:3,:) = reshape((/
     $              1.50d0, 1.455d0, 1.48d0,
     $              1.20d0, 1.350d0, 1.25d0,
     $              1.49d0, 1.250d0, 1.40d0,
     $              
     $              0.128d0, 0.450d0, 0.135d0,
     $              0.148d0, 0.150d0, 0.122d0,
     $              0.142d0, 1.152d0, 0.236d0,
     $              
     $              0.006d0, 0.0600d0, 0.020d0,
     $              0.000d0, 0.0028d0, 0.035d0,
     $              0.020d0, 0.0030d0, 0.040d0,
     $              
     $              4.876d0, 4.825d0, 4.862d0,
     $              4.890d0, 4.871d0, 4.892d0,
     $              4.865d0, 4.757d0, 4.895d0/),
     $              (/3,3,ne/))

               bf_grdpts_id0(1:3,1:3) = reshape((/
     $              2,2,3,
     $              2,2,3,
     $              3,3,3/),
     $              (/3,3/))

               !data for the W_edge
               bf_nodes0(5:6,3:5,:) = reshape((/
     $              1.45d0,  1.26d0,
     $              1.27d0,  1.46d0,
     $              1.26d0,  1.48d0,
     $              
     $            -0.148d0, -1.138d0, 
     $            -0.143d0, -0.146d0, 
     $            -0.129d0, -0.123d0, 
     $              
     $             0.001d0,  0.0025d0,
     $             0.002d0,  0.0100d0,
     $             0.015d0,  0.0800d0,
     $              
     $             4.865d0,  4.85d0, 
     $             4.870d0,  4.89d0, 
     $             4.950d0,  4.83d0/),
     $              (/2,3,ne/))

               bf_nodes1(6:7,3:5,:) = reshape((/
     $              1.35d0, 1.20d0, 
     $              1.25d0, 1.49d0, 
     $              1.28d0, 1.52d0, 
     $              
     $             -0.150d0, -0.148d0, 
     $             -1.152d0, -0.142d0, 
     $             -0.185d0, -0.136d0,
     $              
     $              0.0028d0, 0.00d0,
     $              0.0030d0, 0.02d0,
     $              0.0850d0, 0.06d0,
     $              
     $              4.871d0, 4.890d0,
     $              4.757d0, 4.865d0,
     $              4.785d0, 4.625d0/),
     $              (/2,3,ne/))

               bf_grdpts_id0(5:6,3:5) = reshape((/
     $              3,2,
     $              3,2,
     $              3,2/),
     $              (/2,3/))

               test_newgrdpt_W = [
     $              1.10558898620299d0,
     $             -0.21493341075159d0,
     $              0.49142653305574d0,
     $              4.65867937589912d0]

               test_newgrdpt_NE = [
     $              1.55742091189301d0,
     $              0.56224898725725d0,
     $             -0.14504765609975d0,
     $              4.63647450016455d0]

               do k=1,ne
                  test_newgrdpt(k) = (test_newgrdpt_W(k)+test_newgrdpt_NE(k))/2.0d0
               end do

               test_grdpts_available = [.true.,.true.,.false.,.false.]
               test_ierror = BF_SUCCESS

            case(2)
               bf_newgrdpt_coords0 = [5,4]
               bf_newgrdpt_coords1 = [6,4]

               bf_grdpts_id0(2:4,1:3) = reshape((/
     $              2,2,3,
     $              2,2,3,
     $              3,3,3/),
     $              (/3,3/))

               bf_grdpts_id0(5:5,3:5) = reshape((/
     $              3,
     $              3,
     $              3/),
     $              (/1,3/))

               test_grdpts_available = [.true.,.false.,.false.,.false.]
               test_ierror = .not.BF_SUCCESS

            case(3)
               bf_newgrdpt_coords0 = [3,4]
               bf_newgrdpt_coords1 = [4,4]

               bf_grdpts_id0(1:2,2:4) = reshape((/
     $              2,3,
     $              2,3,
     $              3,3/),
     $              (/2,3/))

               bf_grdpts_id0(4:5,3:5) = reshape((/
     $              3,2,
     $              3,2,
     $              3,2/),
     $              (/2,3/))
               

               test_grdpts_available = [.false.,.true.,.false.,.false.]
               test_ierror = .not.BF_SUCCESS

            case(4)
               bf_newgrdpt_coords0 = [3,5]
               bf_newgrdpt_coords1 = [4,5]

               bf_grdpts_id0(1:2,2:4) = reshape((/
     $              2,3,
     $              2,3,
     $              3,3/),
     $              (/2,3/))

               bf_grdpts_id0(3:4,4:5) = reshape((/
     $              3,2,
     $              3,2/),
     $              (/2,2/))

               test_grdpts_available = [.false.,.false.,.false.,.false.]
               test_ierror = .not.BF_SUCCESS

          end select

        end subroutine get_param_test_compute_from_bf_layer



        function test_are_grdpts_available_to_get_newgrdpt_data(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer                        :: k
          integer       , dimension(7,4) :: bf_grdpts_id0
          integer(ikind), dimension(2)   :: bf_newgrdpt_coords0
          integer                        :: procedure_type
          integer                        :: gradient_type
          integer(ikind), dimension(2,2) :: test_data_needed_bounds0
          logical                        :: test_tmp_array_needed
          logical                        :: test_grdpts_available
          integer(ikind), dimension(2,2) :: data_needed_bounds0
          logical                        :: tmp_array_needed
          logical                        :: grdpts_available

          logical :: test_loc


          test_validated = .true.

          
          data_needed_bounds0 = reshape((/0,1,1,1/),(/2,2/))


          do k=1,5
             
             !input
             call get_param_test_are_grdpts_data(
     $            k,
     $            bf_grdpts_id0,
     $            bf_newgrdpt_coords0,
     $            procedure_type,
     $            gradient_type,
     $            test_data_needed_bounds0,
     $            test_tmp_array_needed,
     $            test_grdpts_available)

             !output
             call are_grdpts_available_to_get_newgrdpt_data(
     $            bf_grdpts_id0,
     $            bf_newgrdpt_coords0,
     $            procedure_type,
     $            gradient_type,
     $            data_needed_bounds0,
     $            tmp_array_needed,
     $            grdpts_available)

             !validation
             test_loc = tmp_array_needed.eqv.test_tmp_array_needed
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''tmp_array_needed('',I2,'') failed'')',k
             end if

             test_loc = grdpts_available.eqv.test_grdpts_available
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''grdpts_available('',I2,'') failed'')',k
             end if

             test_loc = is_int_matrix_validated(
     $            data_needed_bounds0,
     $            test_data_needed_bounds0,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''data_needed_bounds0('',I2,'') failed'')',k
             end if
             
          end do

        end function test_are_grdpts_available_to_get_newgrdpt_data


        subroutine get_param_test_are_grdpts_data(
     $     test_id,
     $     bf_grdpts_id0,
     $     bf_newgrdpt_coords0,
     $     procedure_type,
     $     gradient_type,
     $     test_data_needed_bounds0,
     $     test_tmp_array_needed,
     $     test_grdpts_available)

          implicit none

          integer                       , intent(in)  :: test_id
          integer       , dimension(7,4), intent(out) :: bf_grdpts_id0
          integer(ikind), dimension(2)  , intent(out) :: bf_newgrdpt_coords0
          integer                       , intent(out) :: procedure_type
          integer                       , intent(out) :: gradient_type
          integer(ikind), dimension(2,2), intent(out) :: test_data_needed_bounds0
          logical                       , intent(out) :: test_tmp_array_needed
          logical                       , intent(out) :: test_grdpts_available

          !  _____________
          ! |  _______    |
          ! | |2 2 2 2|   |
          ! | |2_3 3 2|   |
          ! |____0|3_2|___|
          !------------------------------------------------------------
          bf_grdpts_id0 = reshape((/
     $         0,0,0,3,2,0,0,
     $         0,2,3,3,2,0,0,
     $         0,2,2,2,2,0,0,
     $         0,0,0,0,0,0,0/),
     $         (/7,4/))

          procedure_type = SW_edge_type
          gradient_type  = gradient_xI_yLR0_type

          test_data_needed_bounds0 = reshape((/-1,0,2,2/),(/2,2/))

          select case(test_id)
          !  _________
          ! |   ___   |
          ! |  |   |  |
          ! |  |___|  |
          ! |_________|

            case(1)
               bf_newgrdpt_coords0   = [3,1]
               test_tmp_array_needed = .false.
               test_grdpts_available = .true.

          !     _________
          !   _|_        |
          !  | | |       |
          !  |_|_|       |
          !    |_________|

            case(2)
               bf_newgrdpt_coords0   = [1,1]
               test_tmp_array_needed = .true.
               test_grdpts_available = .false.

          !  _________
          ! |        _|_
          ! |       | | |
          ! |       |_|_|
          ! |_________|

            case(3)
               bf_newgrdpt_coords0   = [6,1]
               test_tmp_array_needed = .true.
               test_grdpts_available = .false.
          !     ___
          !  __|___|__ 
          ! |  |___|  |
          ! |         |
          ! |         |
          ! |_________|

            case(4)
               bf_newgrdpt_coords0   = [3,3]
               test_tmp_array_needed = .true.
               test_grdpts_available = .false.

           
          !  _________ 
          ! |         |
          ! |         |
          ! |   ___   |
          ! |__|___|__|
          !    |___|

            case(5)
               bf_newgrdpt_coords0   = [3,0]
               test_tmp_array_needed = .true.
               test_grdpts_available = .false.

          end select

        end subroutine get_param_test_are_grdpts_data


        function test_are_grdpts_available_to_get_newgrdpt_proc(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer               :: k
          integer               :: bf_localization
          integer               :: size_x
          integer               :: size_y
          logical               :: bf_can_exchange_with_neighbor1
          logical               :: bf_can_exchange_with_neighbor2
          integer, dimension(2) :: bf_newgrdpt_coords
          logical               :: test_grdpts_available
          logical               :: grdpts_available
          logical               :: test_loc

          
          test_validated = .true.


          do k=1,210

             !input
             call get_param_test_are_grdpts_proc(
     $            k,
     $            bf_localization,
     $            size_x,size_y,
     $            bf_can_exchange_with_neighbor1,
     $            bf_can_exchange_with_neighbor2,
     $            bf_newgrdpt_coords,
     $            test_grdpts_available)

             !output
             grdpts_available = are_grdpts_available_to_get_newgrdpt_proc(
     $            bf_localization,
     $            size_x, size_y,
     $            bf_can_exchange_with_neighbor1,
     $            bf_can_exchange_with_neighbor2,
     $            bf_newgrdpt_coords)

             !validation
             test_loc = grdpts_available.eqv.test_grdpts_available
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''bf_localization      : '',I1)' , bf_localization
                print '(''size_x               : '',I1)' , size_x
                print '(''size_y               : '',I1)' , size_y
                print '(''neighbor1            : '',L1)' , bf_can_exchange_with_neighbor1
                print '(''neighbor2            : '',L1)' , bf_can_exchange_with_neighbor2
                print '(''bf_newgrdpt_coords   : '',2I2)', bf_newgrdpt_coords
                print '(''test_grdpts_available: '',L1)' , test_grdpts_available
                print '()'
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_are_grdpts_available_to_get_newgrdpt_proc


        subroutine get_param_test_are_grdpts_proc(
     $     test_id,
     $     bf_localization,
     $     size_x,size_y,
     $     bf_can_exchange_with_neighbor1,
     $     bf_can_exchange_with_neighbor2,
     $     bf_newgrdpt_coords,
     $     test_grdpts_available)

          implicit none

          integer               , intent(in)  :: test_id
          integer               , intent(out) :: bf_localization
          integer               , intent(out) :: size_x
          integer               , intent(out) :: size_y
          logical               , intent(out) :: bf_can_exchange_with_neighbor1
          logical               , intent(out) :: bf_can_exchange_with_neighbor2
          integer, dimension(2) , intent(out) :: bf_newgrdpt_coords
          logical               , intent(out) :: test_grdpts_available

          integer :: test_id_loc
          integer :: test_id_card

          !    _       _       _
          ! 8 |_|_____|_|_____|_|
          ! 7 |_|     |_|     |_|
          ! 6 |_|     |_|     |_|
          !    _|      _      |_
          ! 4 |_|     |_|     |_|   21 pts tested = 21 bf_newgrdpt_coords tested
          !    _|      _      |_        |
          ! 2 |_|     |_|     |_|       |___\  test_id_loc [0,20]
          ! 1 |_|_____|_|_____|_|           /
          ! 0 |_|     |_|     |_|
          !    0       3       6
          !
          !            X
          !
          ! 10 localizations tested:  _____\  test_id_card [1,10]
          !   - N                          /
          !   - S
          !   - E ( neighbor1=T, neighbor2=T)
          !   - E ( neighbor1=T, neighbor2=F)
          !   - E ( neighbor1=F, neighbor2=T)
          !   - E ( neighbor1=F, neighbor2=F)
          !   - W ( neighbor1=T, neighbor2=T)
          !   - W ( neighbor1=T, neighbor2=F)
          !   - W ( neighbor1=F, neighbor2=T)
          !   - W ( neighbor1=F, neighbor2=F)
          !------------------------------------------------------------

          size_x = 5
          size_y = 7

          test_id_loc  = mod(test_id-1,21)
          test_id_card = ((test_id - test_id_loc)/21)+1


          !bf_newgrdpt_coords(1)
          select case(mod(test_id_loc,3))
            case(0)
               bf_newgrdpt_coords(1) = 0
            case(1)
               bf_newgrdpt_coords(1) = 3
            case(2)
               bf_newgrdpt_coords(1) = 6
            case default
               print '(''get_param_test_are_grdpts_proc'')'
               print '(''wrong bf_newgrdpt_coords(1)'')'
               stop ''
          end select

          !bf_newgrdpt_coords(2)
          select case((test_id_loc-mod(test_id_loc,3))/3)
            case(0)
               bf_newgrdpt_coords(2) = 0
            case(1)
               bf_newgrdpt_coords(2) = 1
            case(2)
               bf_newgrdpt_coords(2) = 2
            case(3)
               bf_newgrdpt_coords(2) = 4
            case(4)
               bf_newgrdpt_coords(2) = 6
            case(5)
               bf_newgrdpt_coords(2) = 7
            case(6)
               bf_newgrdpt_coords(2) = 8
            case default
               print '(''get_param_test_are_grdpts_proc'')'
               print '(''wrong bf_newgrdpt_coords(2)'')'
               stop ''
          end select
            
          !bf_neighbors
          bf_can_exchange_with_neighbor1 = .true.
          bf_can_exchange_with_neighbor2 = .true.

          select case(test_id_card)
            case(1,2,3,7)
               bf_can_exchange_with_neighbor1 = .true.
               bf_can_exchange_with_neighbor2 = .true.
            case(4)
               bf_can_exchange_with_neighbor2 = .false.
            case(5)
               bf_can_exchange_with_neighbor1 = .false.
            case(6)
               bf_can_exchange_with_neighbor1 = .false.
               bf_can_exchange_with_neighbor2 = .false.
            case(8)
               bf_can_exchange_with_neighbor2 = .false.
            case(9)
               bf_can_exchange_with_neighbor1 = .false.
            case(10)
               bf_can_exchange_with_neighbor1 = .false.
               bf_can_exchange_with_neighbor2 = .false.
            case default
               print '(''get_param_test_are_grdpts_proc'')'
               print '(''wrong can_exchange_with_neighbor1'')'
               print '(''wrong can_exchange_with_neighbor2'')'
               stop ''
          end select

          !bf_localization
          select case(test_id_card)
            case(1)
               bf_localization = N
            case(2)
               bf_localization = S
            case(3,4,5,6)
               bf_localization = E
            case(7,8,9,10)
               bf_localization = W
            case default
               print '(''get_param_test_are_grdpts_proc'')'
               print '(''wrong bf_localization'')'
               stop
          end select

          select case(test_id_card)
            !N
            case(1)
               select case(test_id_loc+1)
                 case(1,2,3,4,5,6)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !S
            case(2)
               select case(test_id_loc+1)
                 case(16,17,18,19,20,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !E(T,T)
            case(3)
               select case(test_id_loc+1)
                 case(1,2,3,4,5,6,7,10,13,16,17,18,19,20,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !E(T,F)
            case(4)
               select case(test_id_loc+1)
                 case(1,2,3,4,5,6,7,10,13,16,19)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !E(F,T)
            case(5)
               select case(test_id_loc+1)
                 case(1,4,7,10,13,16,17,18,19,20,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !E(F,F)
            case(6)
               select case(test_id_loc+1)
                 case(1,4,7,10,13,16,19)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !W(T,T)
            case(7)
               select case(test_id_loc+1)
                 case(1,2,3,4,5,6,9,12,15,16,17,18,19,20,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !W(T,F)
            case(8)
               select case(test_id_loc+1)
                 case(1,2,3,4,5,6,9,12,15,18,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !W(F,T)
            case(9)
               select case(test_id_loc+1)
                 case(3,6,9,12,15,16,17,18,19,20,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !W(F,F)
            case(10)
               select case(test_id_loc+1)
                 case(3,6,9,12,15,18,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            case default
               print '(''get_param_test_are_grdpts_proc'')'
               print '(''wrong test_id_card'')'

          end select

        end subroutine get_param_test_are_grdpts_proc


        subroutine check_inputs()

          implicit none

          real(rkind), dimension(4) :: far_field
          type(pmodel_eq)           :: p_model

          if(ne.ne.4) then
             stop 'the test requires ne=4: DIM2D model'
          end if

          if(.not.is_real_validated(cv_r,2.5d0,.false.)) then
             stop 'the test requires c_v_r=2.5'
          end if

          call p_model%initial_conditions%ini_far_field()

          far_field = p_model%get_far_field(0.0d0,1.0d0,1.0d0)

          if(.not.is_real_vector_validated(
     $         far_field,
     $         [1.46510213931996d0,0.146510214d0,0.0d0,2.84673289046992d0],
     $         .true.)) then
             print '(''the test requires p_model%get_far_field(t,x,y)='')'
             print '(''[1.465102139d0,0.14651021d0,0.0d0,2.84673289d0]'')'
             print '()'
             print '(''T0 should be 0.95'')'
             print '(''flow_direction should be x-direction'')'
             print '(''ic_choice should be newgrdpt_test'')'
             print '()'
             stop ''
             
          end if

        end subroutine check_inputs
        
      end program test_bf_newgrdpt_dispatch

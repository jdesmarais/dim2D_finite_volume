      program test_mainlayer_interface_newgrdpt

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated,
     $       is_real_matrix3D_validated,
     $       is_int_matrix_validated

        use dim2d_parameters, only :
     $       cv_r

        use mainlayer_interface_newgrdpt_class, only :
     $       mainlayer_interface_newgrdpt

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W,
     $       SE_interface_type,
     $       NE_interface_type,
     $       SW_edge_type,
     $       NE_edge_type,
     $       gradient_xI_yLR0_type,
     $       gradient_xLR0_yI_type

        use parameters_constant, only :
     $       N,S,E,W,
     $       obc_eigenqties_bc,
     $       obc_eigenqties_lin

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       obc_eigenqties_strategy

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


        test_loc = test_collect_data_to_compute_newgrdpt(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_collect_data_to_compute_newgrdpt: '',L1)', test_loc
        print '()'


        test_loc = test_are_grdpts_available_to_compute_newgrdpt(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_are_grdpts_available_to_compute_newgrdpt: '',L1)', test_loc
        print '()'


        test_loc = test_compute_newgrdpt(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_newgrdpt: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_compute_newgrdpt(detailled)
     $       result(test_validated)

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_newgrdpt)   :: mainlayer_interface_used
          type(pmodel_eq)                      :: p_model
          real(rkind)                          :: t
          real(rkind)                          :: dt
          real(rkind)   , dimension(nx)        :: interior_x_map
          real(rkind)   , dimension(ny)        :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne)  :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne)  :: interior_nodes1
          type(bf_sublayer), pointer           :: bf_sublayer_ptr
          integer(ikind), dimension(2)         :: gen_newgrdpt_coords
          logical                              :: ierror
         
          integer(ikind), dimension(2)         :: loc_newgrdpt_coords
          real(rkind)   , dimension(ne)        :: test_newgrdpt

          integer :: k

          test_validated = .true.


          call ini_mainlayer_interface_for_test(
     $         mainlayer_interface_used,
     $         bf_sublayer_ptr,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)

          call p_model%initial_conditions%ini_far_field()

          
          do k=1,4

             !input
             call get_param_test_compute_newgrdpt(
     $            k,
     $            mainlayer_interface_used,
     $            t,dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            bf_sublayer_ptr,
     $            gen_newgrdpt_coords,
     $            loc_newgrdpt_coords,
     $            test_newgrdpt)
             
             !output
             ierror = mainlayer_interface_used%compute_newgrdpt(
     $            p_model,
     $            t,dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            bf_sublayer_ptr,
     $            gen_newgrdpt_coords)

             !validation
             test_loc = is_real_vector_validated(
     $            bf_sublayer_ptr%nodes(
     $               loc_newgrdpt_coords(1),
     $               loc_newgrdpt_coords(2),
     $               :),
     $            test_newgrdpt,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_compute_newgrdpt


        subroutine get_param_test_compute_newgrdpt(
     $     test_id,
     $     mainlayer_interface_used,
     $     t,dt,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1,
     $     bf_sublayer_ptr,
     $     gen_newgrdpt_coords,
     $     loc_newgrdpt_coords,
     $     test_newgrdpt)

          implicit none

          integer                            , intent(in)    :: test_id
          type(mainlayer_interface_newgrdpt) , intent(in)    :: mainlayer_interface_used
          real(rkind)                        , intent(out)   :: t
          real(rkind)                        , intent(out)   :: dt
          real(rkind)   , dimension(nx)      , intent(inout) :: interior_x_map
          real(rkind)   , dimension(ny)      , intent(inout) :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne), intent(inout) :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne), intent(inout) :: interior_nodes1
          type(bf_sublayer), pointer         , intent(inout) :: bf_sublayer_ptr
          integer(ikind), dimension(2)       , intent(out)   :: gen_newgrdpt_coords
          integer(ikind), dimension(2)       , intent(out)   :: loc_newgrdpt_coords
          real(rkind)   , dimension(ne)      , intent(out)   :: test_newgrdpt


          type(bf_sublayer), pointer :: nbf_sublayer_ptr
          integer(ikind)             :: i,j


          t  = 0.0d0
          dt = 0.25d0

          select case(test_id)

            !NEWGRDPT_NO_ERROR
            case(1)

               bf_sublayer_ptr%x_map(1)   = 0.00d0
               bf_sublayer_ptr%x_map(2)   = 0.50d0
               bf_sublayer_ptr%x_map(5:8) = [0.5d0, 1.5d0, 2.5d0, 3.5d0]

               bf_sublayer_ptr%y_map(1)   = 0.00d0
               bf_sublayer_ptr%y_map(2)   = 0.25d0
               bf_sublayer_ptr%y_map(6:9) = [0.0d0,0.25d0,0.5d0,0.75d0]

               bf_sublayer_ptr%nodes(5:8,6:9,:) = reshape((/
     $              1.50d0, 1.455d0, 1.48d0, 0.0d0,
     $              1.20d0, 1.350d0, 1.25d0, 0.0d0,
     $              1.49d0, 1.250d0, 1.40d0, 0.0d0,
     $              0.00d0, 0.000d0, 0.00d0, 0.0d0,
     $              
     $              0.128d0, 0.450d0, 0.135d0, 0.0d0,
     $              0.148d0, 0.150d0, 0.122d0, 0.0d0,
     $              0.142d0, 1.152d0, 0.236d0, 0.0d0,
     $              0.000d0, 0.000d0, 0.000d0, 0.0d0,
     $              
     $              0.006d0, 0.0600d0, 0.020d0, 0.0d0,
     $              0.000d0, 0.0028d0, 0.035d0, 0.0d0,
     $              0.020d0, 0.0030d0, 0.040d0, 0.0d0,
     $              0.000d0, 0.0000d0, 0.000d0, 0.0d0,
     $              
     $              4.876d0, 4.825d0, 4.862d0, 0.0d0,
     $              4.890d0, 4.871d0, 4.892d0, 0.0d0,
     $              4.865d0, 4.757d0, 4.895d0, 0.0d0,
     $              0.000d0, 0.000d0, 0.000d0, 0.0d0
     $              /),
     $              (/4,4,ne/))

               bf_sublayer_ptr%bf_compute_used%x_map_tmp(5:7) = [0.5d0, 1.5d0 , 2.5d0]
               bf_sublayer_ptr%bf_compute_used%y_map_tmp(6:8) = [0.0d0, 0.25d0, 0.5d0]
               bf_sublayer_ptr%bf_compute_used%nodes_tmp(5:7,6:8,:) = reshape((/
     $              1.48d0, 1.30d0, 1.35d0,
     $              1.26d0, 1.45d0, 1.4d0,
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

               gen_newgrdpt_coords = [align_E+5,align_S+5]
               loc_newgrdpt_coords = [8,7]


               select case(obc_eigenqties_strategy)

                 case(obc_eigenqties_bc)
                    test_newgrdpt = [
     $                   1.22383078395524d0,
     $                   0.39531842478603d0,
     $                   -0.21050217290879d0,
     $                   4.19684181018688d0]

                 case(obc_eigenqties_lin)
                    test_newgrdpt = [
     $                   1.21167555521982d0,
     $                   0.35901827468671d0,
     $                   -0.20475732388332d0,
     $                   4.20002914561351d0]
                    
                 case default
                    print '(''test_bf_newgrdpt_prim'')'
                    print '(''test_sym_compute_newgrdpt_x'')'
                    print '(''obc_eigenqties_strategy not recognized'')'
                    stop ''
                    
               end select


            !NEWGRDPT_PROC_NBGRDPTS_ERROR
            case(2)

               interior_x_map = (/ ((i-1)*1.00d0,i=1,nx) /)
               interior_y_map = (/ ((j-1)*0.25d0,j=1,ny) /)

               !East buffer layer
               bf_sublayer_ptr%x_map(1)       = 0.00d0
               bf_sublayer_ptr%x_map(2)       = 0.50d0
               bf_sublayer_ptr%x_map(5:8)     = [0.5d0, 1.5d0, 2.5d0, 3.5d0]

               bf_sublayer_ptr%y_map(1)       = 0.00d0
               bf_sublayer_ptr%y_map(2)       = 0.25d0
               bf_sublayer_ptr%y_map(ny-1:ny) = [0.0d0,0.25d0] !,0.5d0,0.75d0]

               bf_sublayer_ptr%nodes(5:8,ny-1:ny,:) = reshape((/
     $              1.50d0, 1.455d0, 1.48d0, 0.0d0,
     $              1.20d0, 1.350d0, 1.25d0, 0.0d0,
     $              
     $              0.128d0, 0.450d0, 0.135d0, 0.0d0,
     $              0.148d0, 0.150d0, 0.122d0, 0.0d0,
     $              
     $              0.006d0, 0.0600d0, 0.020d0, 0.0d0,
     $              0.000d0, 0.0028d0, 0.035d0, 0.0d0,
     $              
     $              4.876d0, 4.825d0, 4.862d0, 0.0d0,
     $              4.890d0, 4.871d0, 4.892d0, 0.0d0
     $              /),
     $              (/4,2,ne/))

               bf_sublayer_ptr%bf_compute_used%x_map_tmp(5:7) = [0.5d0, 1.5d0 , 2.5d0]
               bf_sublayer_ptr%bf_compute_used%y_map_tmp(ny-1:ny) = [0.0d0, 0.25d0] !, 0.5d0]
               bf_sublayer_ptr%bf_compute_used%nodes_tmp(5:7,ny-1:ny,:) = reshape((/
     $              1.48d0, 1.30d0, 1.35d0,
     $              1.26d0, 1.45d0, 1.4d0,
     $              
     $              0.128d0, 0.127d0, 0.142d0,
     $              1.138d0, 0.148d0, 0.132d0,
     $              
     $              0.0050d0, 0.020d0, 0.060d0,
     $              0.0025d0, 0.001d0, 0.015d0,
     $              
     $              4.88d0, 4.870d0, 4.855d0,
     $              4.85d0, 4.865d0, 4.845d0/),
     $              (/3,2,ne/))

               gen_newgrdpt_coords = [align_E+5,align_N+1]
               loc_newgrdpt_coords = [8,ny]


               select case(obc_eigenqties_strategy)

                 case(obc_eigenqties_bc)
                    test_newgrdpt = [
     $                   1.22383078395524d0,
     $                   0.39531842478603d0,
     $                   -0.21050217290879d0,
     $                   4.19684181018688d0]

                 case(obc_eigenqties_lin)
                    test_newgrdpt = [
     $                   1.21167555521982d0,
     $                   0.35901827468671d0,
     $                   -0.20475732388332d0,
     $                   4.20002914561351d0]
                    
                 case default
                    print '(''test_bf_newgrdpt_prim'')'
                    print '(''test_sym_compute_newgrdpt_x'')'
                    print '(''obc_eigenqties_strategy not recognized'')'
                    stop ''
                    
               end select

               !North buffer layer
               nbf_sublayer_ptr => mainlayer_interface_used%NE_interface_N_ptr
               
               nbf_sublayer_ptr%y_map(2:5) = [0.0d0,0.25d0,0.5d0,0.75d0]

               nbf_sublayer_ptr%nodes(12:15,3:6,:) = reshape((/
     $              1.50d0, 1.455d0, 1.48d0, 0.0d0,
     $              1.20d0, 1.350d0, 1.25d0, 0.0d0,
     $              1.49d0, 1.250d0, 1.40d0, 0.0d0,
     $              0.00d0, 0.000d0, 0.00d0, 0.0d0,
     $              
     $              0.128d0, 0.450d0, 0.135d0, 0.0d0,
     $              0.148d0, 0.150d0, 0.122d0, 0.0d0,
     $              0.142d0, 1.152d0, 0.236d0, 0.0d0,
     $              0.000d0, 0.000d0, 0.000d0, 0.0d0,
     $              
     $              0.006d0, 0.0600d0, 0.020d0, 0.0d0,
     $              0.000d0, 0.0028d0, 0.035d0, 0.0d0,
     $              0.020d0, 0.0030d0, 0.040d0, 0.0d0,
     $              0.000d0, 0.0000d0, 0.000d0, 0.0d0,
     $              
     $              4.876d0, 4.825d0, 4.862d0, 0.0d0,
     $              4.890d0, 4.871d0, 4.892d0, 0.0d0,
     $              4.865d0, 4.757d0, 4.895d0, 0.0d0,
     $              0.000d0, 0.000d0, 0.000d0, 0.0d0
     $              /),
     $              (/4,4,ne/))

               nbf_sublayer_ptr%bf_compute_used%y_map_tmp(3:5) = [0.0d0, 0.25d0, 0.5d0]
               nbf_sublayer_ptr%bf_compute_used%nodes_tmp(12:14,3:5,:) = reshape((/
     $              1.48d0, 1.30d0, 1.35d0,
     $              1.26d0, 1.45d0, 1.4d0,
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

            !NEWGRDPT_DATA_NBGRDPTS_ERROR
            case(3)

               interior_x_map = (/ ((i-1)*1.00d0,i=1,nx) /)
               interior_y_map = (/ ((j-1)*0.25d0,j=1,ny) /)

               !East buffer layer
               bf_sublayer_ptr%bf_compute_used%grdpts_id_tmp(:,1:4) = reshape((/
     $              1,1,1,1,1,2,3,3,
     $              1,1,1,1,1,2,3,0,
     $              1,1,1,1,1,2,3,0,
     $              1,1,1,1,1,2,3,0/),
     $              (/8,4/))

               bf_sublayer_ptr%bf_compute_used%nodes_tmp(6:8,1:3,:) = reshape((/
     $              1.26d0, 1.45d0, 1.40d0,
     $              1.46d0, 1.27d0, 1.47d0,
     $              1.48d0, 1.26d0, 1.41d0,
     $              
     $              1.138d0, 0.148d0, 0.132d0,
     $              0.146d0, 0.143d0, 0.145d0,
     $              0.123d0, 0.129d0, 0.124d0,
     $              
     $              0.0025d0, 0.001d0, 0.015d0,
     $                0.01d0, 0.002d0,  0.05d0,
     $                0.08d0, 0.015d0,  0.09d0,
     $              
     $              4.85d0, 4.865d0, 4.845d0,
     $              4.89d0, 4.870d0, 4.860d0,
     $              4.83d0, 4.950d0, 4.62d0/),
     $              (/3,3,ne/))

               bf_sublayer_ptr%nodes(6:8,1:3,:) = reshape((/
     $               1.20d0, 1.35d0, 1.25d0,
     $               1.49d0, 1.25d0, 1.40d0,
     $               1.52d0, 1.28d0, 1.45d0,
     $              
     $              0.148d0, 0.150d0, 0.122d0,
     $              0.142d0, 1.152d0, 0.236d0,
     $              0.136d0, 0.185d0, 0.296d0,
     $              
     $               0.0d0, 0.0028d0, 0.035d0,
     $              0.02d0, 0.0030d0, 0.040d0,
     $              0.06d0,  0.085d0, 0.015d0,
     $              
     $              4.89d0,  4.871d0, 4.892d0,
     $             4.865d0,  4.757d0, 4.895d0,
     $             4.625d0,  4.785d0, 4.825d0/),
     $              (/3,3,ne/))

               gen_newgrdpt_coords = [align_E+5,align_S]
               loc_newgrdpt_coords = [8,2]

               !south buffer layer
               nbf_sublayer_ptr => mainlayer_interface_used%SE_interface_S_ptr

               nbf_sublayer_ptr%grdpts_id(8:11,4:9) = reshape((/
     $              1,1,2,3,
     $              1,2,2,3,
     $              1,2,3,3,
     $              1,2,3,0,
     $              1,2,3,0,
     $              1,2,3,0/),
     $              (/4,6/))

               nbf_sublayer_ptr%bf_compute_used%nodes_tmp(9:11,5:8,:) = reshape((/
     $              1.48d0, 1.30d0, 1.35d0,
     $              1.26d0, 1.45d0, 1.40d0,
     $              1.46d0, 1.27d0, 1.47d0,
     $              1.48d0, 1.26d0, 1.41d0,
     $              
     $              0.128d0, 0.127d0, 0.142d0,
     $              1.138d0, 0.148d0, 0.132d0,
     $              0.146d0, 0.143d0, 0.145d0,
     $              0.123d0, 0.129d0, 0.124d0,
     $              
     $              0.0050d0, 0.020d0, 0.060d0,
     $              0.0025d0, 0.001d0, 0.015d0,
     $              0.0100d0, 0.002d0, 0.050d0,
     $              0.0800d0, 0.015d0, 0.090d0,
     $              
     $              4.88d0, 4.870d0, 4.855d0,
     $              4.85d0, 4.865d0, 4.845d0,
     $              4.89d0, 4.870d0, 4.860d0,
     $              4.83d0, 4.950d0, 4.62d0/),
     $              (/3,4,ne/))

               nbf_sublayer_ptr%nodes(9:11,5:8,:) = reshape((/
     $              1.50d0, 1.455d0, 1.48d0,
     $              1.20d0, 1.350d0, 1.25d0,
     $              1.49d0, 1.250d0, 1.40d0,
     $              1.52d0, 1.28d0, 1.45d0,
     $              
     $              0.128d0, 0.450d0, 0.135d0,
     $              0.148d0, 0.150d0, 0.122d0,
     $              0.142d0, 1.152d0, 0.236d0,
     $              0.136d0, 0.185d0, 0.296d0,
     $              
     $              0.006d0, 0.0600d0, 0.02d0,
     $              0.000d0, 0.0028d0, 0.035d0,
     $              0.020d0, 0.0030d0, 0.04d0,
     $               0.06d0, 0.0850d0, 0.015d0,
     $              
     $              4.876d0, 4.825d0, 4.862d0,
     $              4.890d0, 4.871d0, 4.892d0,
     $              4.865d0, 4.757d0, 4.895d0,
     $              4.625d0, 4.785d0, 4.825d0/),
     $              (/3,4,ne/))

               !test_newgrdpt
               test_newgrdpt = [
     $              1.09198841545674d0,
     $              0.34791358035954d0,
     $              0.22873420436383d0,
     $              3.76108958832707d0]


            !NEWGRDPT_NO_PREVIOUS_DATA_ERROR
            case(4)

               interior_x_map = (/ ((i-1)*1.00d0,i=1,nx) /)
               interior_y_map = (/ ((j-1)*0.25d0,j=1,ny) /)
               
               call bf_sublayer_ptr%deallocate_after_timeInt()
               deallocate(bf_sublayer_ptr%grdpts_id)
               deallocate(bf_sublayer_ptr%nodes)

               bf_sublayer_ptr%alignment = reshape((/
     $              align_E,align_S+6,align_E,align_S+6/),
     $              (/2,2/)) 

               allocate(bf_sublayer_ptr%grdpts_id(5,5))
               bf_sublayer_ptr%grdpts_id = reshape((/
     $              1,1,2,3,0,
     $              1,1,2,3,0,
     $              1,1,2,3,0,
     $              1,1,2,3,0,
     $              1,1,2,3,0/),
     $              (/5,5/))

               allocate(bf_sublayer_ptr%nodes(5,5,ne))
               bf_sublayer_ptr%nodes(2:5,2:5,:) = reshape((/
     $              1.50d0, 1.455d0, 1.48d0, 0.0d0,
     $              1.20d0, 1.350d0, 1.25d0, 0.0d0,
     $              1.49d0, 1.250d0, 1.40d0, 0.0d0,
     $              0.00d0, 0.000d0, 0.00d0, 0.0d0,
     $              
     $              0.128d0, 0.450d0, 0.135d0, 0.0d0,
     $              0.148d0, 0.150d0, 0.122d0, 0.0d0,
     $              0.142d0, 1.152d0, 0.236d0, 0.0d0,
     $              0.000d0, 0.000d0, 0.000d0, 0.0d0,
     $              
     $              0.006d0, 0.0600d0, 0.020d0, 0.0d0,
     $              0.000d0, 0.0028d0, 0.035d0, 0.0d0,
     $              0.020d0, 0.0030d0, 0.040d0, 0.0d0,
     $              0.000d0, 0.0000d0, 0.000d0, 0.0d0,
     $              
     $              4.876d0, 4.825d0, 4.862d0, 0.0d0,
     $              4.890d0, 4.871d0, 4.892d0, 0.0d0,
     $              4.865d0, 4.757d0, 4.895d0, 0.0d0,
     $              0.000d0, 0.000d0, 0.000d0, 0.0d0
     $              /),
     $              (/4,4,ne/))

               gen_newgrdpt_coords = [align_E+2,align_S+6]
               loc_newgrdpt_coords = [5,3]

               call bf_sublayer_ptr%set_neighbor1_share(.false.)
               call bf_sublayer_ptr%set_neighbor2_share(.false.)


               interior_nodes0(nx-2:nx,align_S+5:align_S+7,:) = reshape((/
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

               interior_nodes1(nx-2:nx,align_S+5:align_S+7,:) = reshape((/
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


               select case(obc_eigenqties_strategy)

                 case(obc_eigenqties_bc)
                    test_newgrdpt = [
     $                   1.22383078395524d0,
     $                   0.39531842478603d0,
     $                   -0.21050217290879d0,
     $                   4.19684181018688d0]

                 case(obc_eigenqties_lin)
                    test_newgrdpt = [
     $                   1.21167555521982d0,
     $                   0.35901827468671d0,
     $                   -0.20475732388332d0,
     $                   4.20002914561351d0]
                    
                 case default
                    print '(''test_bf_newgrdpt_prim'')'
                    print '(''test_sym_compute_newgrdpt_x'')'
                    print '(''obc_eigenqties_strategy not recognized'')'
                    stop ''
                    
               end select

          end select

        end subroutine get_param_test_compute_newgrdpt


        function test_are_grdpts_available_to_compute_newgrdpt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(mainlayer_interface_newgrdpt)   :: mainlayer_interface_used
          logical, dimension(4)                :: grdpts_available
          integer                              :: nb_procedures
          integer, dimension(4)                :: procedure_type
          integer, dimension(4)                :: gradient_type
          integer, dimension(2)                :: tmp_newgrdpt_coords
          integer, dimension(6,6)              :: tmp_grdpts_id
          logical                              :: all_grdpts_available_test
          logical                              :: all_grdpts_available


          test_validated = .true.


          !test 1
          grdpts_available = [.false.,.false.,.false.,.false.]
          nb_procedures = 2
          procedure_type(1) = SW_edge_type
          procedure_type(2) = NE_edge_type
          gradient_type(1)  = gradient_xI_yLR0_type
          gradient_type(2)  = gradient_xLR0_yI_type
          
          tmp_newgrdpt_coords = [3,3]
          tmp_grdpts_id = reshape((/
     $         1,1,1,0,0,0,
     $         1,1,1,0,0,0,
     $         1,1,1,1,1,0,
     $         1,1,1,1,1,0,
     $         0,1,1,1,1,0,
     $         0,0,0,0,0,0/),
     $         (/6,6/))

          all_grdpts_available_test = .true.

          all_grdpts_available = mainlayer_interface_used%are_grdpts_available_to_compute_newgrdpt(
     $         grdpts_available,
     $         nb_procedures,
     $         procedure_type,
     $         gradient_type,
     $         tmp_newgrdpt_coords,
     $         tmp_grdpts_id)
          test_loc = all_grdpts_available.eqv.all_grdpts_available_test
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test 1 failed'')'
          end if


          !test 2
          tmp_grdpts_id(1,1) = 0
          all_grdpts_available_test = .false.

          all_grdpts_available = mainlayer_interface_used%are_grdpts_available_to_compute_newgrdpt(
     $         grdpts_available,
     $         nb_procedures,
     $         procedure_type,
     $         gradient_type,
     $         tmp_newgrdpt_coords,
     $         tmp_grdpts_id)
          test_loc = all_grdpts_available.eqv.all_grdpts_available_test
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test 2 failed'')'
          end if

        end function test_are_grdpts_available_to_compute_newgrdpt


        function test_collect_data_to_compute_newgrdpt(detailled)
     $       result(test_validated)

          implicit none
          
          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_newgrdpt)   :: mainlayer_interface_used
          type(bf_sublayer), pointer           :: bf_sublayer_ptr
          real(rkind)   , dimension(nx)        :: interior_x_map
          real(rkind)   , dimension(ny)        :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne)  :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne)  :: interior_nodes1
          integer(ikind), dimension(2)         :: gen_newgrdpt_coords
          integer(ikind), dimension(2,2)       :: data_needed_bounds0
          
          integer       , dimension(8,ny+6)    :: tmp_grdpts_id
          real(rkind)   , dimension(8)         :: tmp_x_map
          real(rkind)   , dimension(ny+6)      :: tmp_y_map
          real(rkind)   , dimension(8,ny+6,ne) :: tmp_nodes0
          real(rkind)   , dimension(8,ny+6,ne) :: tmp_nodes1

          integer       , dimension(8,ny+6)    :: tmp_grdpts_id_test
          real(rkind)   , dimension(8)         :: tmp_x_map_test
          real(rkind)   , dimension(ny+6)      :: tmp_y_map_test
          real(rkind)   , dimension(8,ny+6,ne) :: tmp_nodes0_test
          real(rkind)   , dimension(8,ny+6,ne) :: tmp_nodes1_test

          integer(ikind) :: i,j
          integer        :: k
          logical        :: test_loc


          test_validated = .true.


          !input
          gen_newgrdpt_coords(1) = align_E-2
          gen_newgrdpt_coords(2) = align_S-1
          
          data_needed_bounds0(1,1) = -1
          data_needed_bounds0(2,1) = -3
          data_needed_bounds0(1,2) =  6
          data_needed_bounds0(2,2) = -align_S+align_N+5

          call ini_mainlayer_interface_for_test(
     $         mainlayer_interface_used,
     $         bf_sublayer_ptr,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)

          !grdpts_id
          do j=1,4
             tmp_grdpts_id_test(:,j) = 
     $            [1,1,1,1,1,1,1,2]
          end do

          tmp_grdpts_id_test(:,5) =
     $         [1,1,1,1,1,1,2,2]

          do j=6,ny+6
             tmp_grdpts_id_test(:,j) = 
     $            [1,1,1,1,1,1,2,3]
          end do

          !x_map
          tmp_x_map_test = (/(align_E-3+(i-1),i=1,8)/)

          !y_map
          tmp_y_map_test = (/(align_S-4+(j-1),j=1,ny+6)/)

          !nodes0
          tmp_nodes0_test = reshape((/
     $         (((-200*(k-1)-20*(align_S-5+j-1)-(align_E-4+i-1),i=1,8),j=1,ny+6),k=1,ne)/),
     $         (/8,ny+6,ne/))

          !nodes1
          tmp_nodes1_test = reshape((/
     $         ((( 200*(k-1)+20*(align_S-5+j-1)+(align_E-4+i-1),i=1,8),j=1,ny+6),k=1,ne)/),
     $         (/8,ny+6,ne/))


          !output
          call mainlayer_interface_used%collect_data_to_compute_newgrdpt(
     $         bf_sublayer_ptr,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         gen_newgrdpt_coords,
     $         data_needed_bounds0,
     $         tmp_grdpts_id,
     $         tmp_x_map,
     $         tmp_y_map,
     $         tmp_nodes0,
     $         tmp_nodes1)


          !validation
          test_loc = is_real_vector_validated(
     $         tmp_x_map,
     $         tmp_x_map_test,
     $         detailled)
          if(detailled.and.(.not.test_loc)) then
             print '(''test x_map failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         tmp_y_map,
     $         tmp_y_map_test,
     $         detailled)
          if(detailled.and.(.not.test_loc)) then
             print '(''test y_map failed'')'
          end if

          test_loc = is_int_matrix_validated(
     $         tmp_grdpts_id,
     $         tmp_grdpts_id_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test grdpts_id failed'')'
          end if

          test_loc = is_real_matrix3D_validated(
     $         tmp_nodes0,
     $         tmp_nodes0_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nodes0 failed'')'
          end if
          
          test_loc = is_real_matrix3D_validated(
     $         tmp_nodes1,
     $         tmp_nodes1_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nodes1 failed'')'
          end if

        end function test_collect_data_to_compute_newgrdpt


        subroutine ini_mainlayer_interface_for_test(
     $     mainlayer_interface_used,
     $     bf_sublayer_ptr_tested,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1)

          implicit none

          type(mainlayer_interface_newgrdpt)  , intent(inout) :: mainlayer_interface_used
          type(bf_sublayer), pointer          , intent(out)   :: bf_sublayer_ptr_tested
          real(rkind)   , dimension(nx)       , intent(out)   :: interior_x_map
          real(rkind)   , dimension(ny)       , intent(out)   :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne) , intent(out)   :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne) , intent(out)   :: interior_nodes1


          type(bf_sublayer), pointer :: bf_sublayer_ptr
          integer(ikind)             :: i,j
          integer                    :: k


          interior_x_map  = (/ (i,i=1,nx) /)
          interior_y_map  = (/ (j,j=1,ny) /)
          interior_nodes0 = reshape((/
     $         (((-200*(k-1)-20*(j-1)-(i-1),i=1,nx),j=1,ny),k=1,ne) /),
     $         (/nx,ny,ne/))
          interior_nodes1 = reshape((/
     $         (((200*(k-1)+20*(j-1)+(i-1),i=1,nx),j=1,ny),k=1,ne) /),
     $         (/nx,ny,ne/))


          !South buffer layer
          !============================================================
          allocate(bf_sublayer_ptr)

          !alignment1
          !------------------------------------------------------------
          call bf_sublayer_ptr%ini(S)
          bf_sublayer_ptr%alignment = reshape((/
     $         align_E-3,align_S-4,align_E+3,align_S/),
     $         (/2,2/))

          !grdpts_id1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%grdpts_id(11,9))
          bf_sublayer_ptr%grdpts_id=reshape((/
     $         3,3,3,3,3,3,3,3,3,3,3,
     $         3,2,2,2,2,2,2,2,2,2,3,
     $         3,2,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,2,3,
     $         2,2,1,1,1,1,1,1,2,2,3,
     $         1,1,1,1,1,1,1,1,2,3,3,
     $         1,1,1,1,1,1,1,1,2,3,0/),
     $         (/11,9/))

          !x_map1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%x_map(11))
          bf_sublayer_ptr%x_map = (/ (align_E-5+(i-1),i=1,11) /)
          
          !y_map1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%y_map(9))
          bf_sublayer_ptr%y_map = (/ (align_S-6+(j-1),j=1,9) /)

          !nodes1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%nodes(11,9,ne))
          bf_sublayer_ptr%nodes = reshape((/
     $         (((200*(k-1)+20*(align_S-7+j-1)+(align_E-6+i-1),i=1,11),j=1,9),k=1,ne)/),
     $         (/11,9,ne/))
          

          !grdpts_id0
          !------------------------------------------------------------
          call bf_sublayer_ptr%allocate_before_timeInt()

          !nodes0
          !------------------------------------------------------------
          bf_sublayer_ptr%bf_compute_used%nodes_tmp = reshape((/
     $         (((-200*(k-1)-20*(align_S-7+j-1)-(align_E-6+i-1),i=1,11),j=1,9),k=1,ne)/),
     $         (/11,9,ne/))


          !neighbor1
          !------------------------------------------------------------
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         bf_sublayer_ptr)


          !East buffer layer
          !============================================================
          allocate(bf_sublayer_ptr)
          bf_sublayer_ptr_tested => bf_sublayer_ptr

          !alignment1
          !------------------------------------------------------------
          call bf_sublayer_ptr%ini(E)
          bf_sublayer_ptr%alignment = reshape((/
     $         align_E,align_S+1,align_E+3,align_N-1/),
     $         (/2,2/))

          !grdpts_id1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%grdpts_id(8,ny))
          bf_sublayer_ptr%grdpts_id(:,1:3) = reshape((/
     $         1,1,1,1,1,1,2,3,
     $         1,1,1,1,1,2,2,3,
     $         1,1,1,1,1,2,3,3/),
     $         (/8,3/))
          do j=4,ny
             bf_sublayer_ptr%grdpts_id(:,j) = [1,1,1,1,1,2,3,0]
          end do

          !x_map1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%x_map(8))
          bf_sublayer_ptr%x_map = (/ (align_E-2+(i-1),i=1,8) /)
          
          !y_map1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%y_map(ny))
          bf_sublayer_ptr%y_map = (/ (align_S-1+(j-1),j=1,ny) /)

          !nodes1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%nodes(8,ny,ne))
          bf_sublayer_ptr%nodes = reshape((/
     $         (((200*(k-1)+20*(align_S-2+j-1)+(align_E-3+i-1),i=1,8),j=1,ny),k=1,ne)/),
     $         (/8,ny,ne/))
          

          !grdpts_id0
          !------------------------------------------------------------
          call bf_sublayer_ptr%allocate_before_timeInt()

          !nodes0
          !------------------------------------------------------------
          bf_sublayer_ptr%bf_compute_used%nodes_tmp = reshape((/
     $         (((-200*(k-1)-20*(align_S-2+j-1)-(align_E-3+i-1),i=1,8),j=1,ny),k=1,ne)/),
     $         (/8,ny,ne/))


          !neighbor1
          !------------------------------------------------------------
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         bf_sublayer_ptr)

          !neighbor2
          !------------------------------------------------------------
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         NE_interface_type,
     $         bf_sublayer_ptr)


          !North buffer layer
          !============================================================
          allocate(bf_sublayer_ptr)

          !alignment1
          !------------------------------------------------------------
          call bf_sublayer_ptr%ini(N)
          bf_sublayer_ptr%alignment = reshape((/
     $         align_E-7,align_N,align_E+3,align_N+4/),
     $         (/2,2/))

          !grdpts_id1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%grdpts_id(15,9))
          bf_sublayer_ptr%grdpts_id = reshape((/
     $         1,1,1,1,1,1,1,1,1,1,1,1,2,3,0,
     $         1,1,1,1,1,1,1,1,1,1,1,1,2,3,0,
     $         3,3,1,1,1,1,1,1,1,1,1,1,2,3,0,
     $         2,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
     $         3,2,2,2,2,2,2,2,2,2,2,2,2,3,0,
     $         3,3,3,3,3,3,3,3,3,3,3,3,3,3,0/),
     $         (/15,9/))

          !x_map1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%x_map(15))
          bf_sublayer_ptr%x_map = (/ (align_E-9+(i-1),i=1,15) /)
          
          !y_map1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%y_map(9))
          bf_sublayer_ptr%y_map = (/ (align_N-2+(j-1),j=1,9) /)

          !nodes1
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%nodes(15,9,ne))
          bf_sublayer_ptr%nodes = reshape((/
     $         (((200*(k-1)+20*(align_N-3+j-1)+(align_E-10+i-1),i=1,15),j=1,9),k=1,ne)/),
     $         (/15,9,ne/))
          

          !grdpts_id0
          !------------------------------------------------------------
          call bf_sublayer_ptr%allocate_before_timeInt()

          !nodes0
          !------------------------------------------------------------
          bf_sublayer_ptr%bf_compute_used%nodes_tmp = reshape((/
     $         (((-200*(k-1)-20*(align_N-3+j-1)-(align_E-10+i-1),i=1,15),j=1,9),k=1,ne)/),
     $         (/15,9,ne/))


          !neighbor2
          !------------------------------------------------------------
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         NE_interface_type,
     $         bf_sublayer_ptr)

        end subroutine ini_mainlayer_interface_for_test


        subroutine check_inputs()

          implicit none

          type(pmodel_eq)            :: p_model
          real(rkind), dimension(ne) :: far_field

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


      end program test_mainlayer_interface_newgrdpt

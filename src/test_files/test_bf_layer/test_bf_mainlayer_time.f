      program test_bf_mainlayer_time

        use bc_operators_class, only :
     $       bc_operators

        use bf_mainlayer_time_class, only :
     $       bf_mainlayer_time

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_int_vector_validated,
     $       is_int_matrix_validated,
     $       is_real_matrix3D_validated

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W,
     $       
     $       NE_edge_type,
     $       NW_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       N_edge_type,
     $       no_overlap,
     $       NS_overlap

        use parameters_constant, only :
     $       N

        use parameters_input, only :
     $       nx,ny,ne,
     $       x_min,y_min

        use parameters_kind, only :
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use rk3tvd_steps_module, only :
     $       compute_1st_step,
     $       compute_1st_step_nopt
        
        use sd_operators_class, only :
     $       sd_operators
        
        use td_operators_class, only :
     $       td_operators

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_initialize_before_timeInt(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_initialize_before_timeInt: '',L1)', test_loc
        print '()'


        test_loc = test_finalize_after_timeInt(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_finalize_after_timeInt: '',L1)', test_loc
        print '()'
        
        
        test_loc = test_compute_time_dev(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_time_dev: '',L1)', test_loc
        print '()'


        test_loc = test_compute_integration_step(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_integration_step: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_initialize_before_timeInt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_mainlayer_time) :: bf_mainlayer_used
          integer, dimension(2)   :: x_borders_test
          integer, dimension(2)   :: y_borders_test
          integer, dimension(5,7) :: bc_sections_test

          integer                    :: nb_sublayers
          type(bf_sublayer), pointer :: current_sublayer

          integer :: k
          logical :: test_loc


          test_validated = .true.

          
          !input
          call initialize_bf_mainlayer_for_tests(bf_mainlayer_used)

          x_borders_test = [1,6]
          y_borders_test = [3,6]
          bc_sections_test = reshape((/
     $         NW_edge_type  , 1,1, no_overlap, NS_overlap,
     $         NE_edge_type  , 5,1, no_overlap, NS_overlap,
     $         W_edge_type   , 1,3, 4         , no_overlap,
     $         E_edge_type   , 5,3, 4         , no_overlap,
     $         NW_corner_type, 1,5, no_overlap, no_overlap,
     $         N_edge_type   , 3,5, 4         , no_overlap,
     $         NE_corner_type, 5,5, no_overlap, no_overlap/),
     $         (/5,7/))


          !output
          call bf_mainlayer_used%initialize_before_timeInt()


          !validation
          nb_sublayers = bf_mainlayer_used%get_nb_sublayers()
          current_sublayer => bf_mainlayer_used%get_head_sublayer()

          do k=1, nb_sublayers

             !integration borders
             test_loc = is_int_vector_validated(
     $            current_sublayer%x_borders,
     $            x_borders_test,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''x_borders failed'')'
             end if

             test_loc = is_int_vector_validated(
     $            current_sublayer%y_borders,
     $            y_borders_test,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''y_borders failed'')'
             end if

             !bc_sections
             test_loc = is_int_matrix_validated(
     $            current_sublayer%bc_sections,
     $            bc_sections_test,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''bc_sections failed'')'
             end if

             !arrays for time integration allocated
             test_loc = allocated(current_sublayer%bf_compute_used%alignment_tmp)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''alignment_tmp failed'')'
             end if

             test_loc = allocated(current_sublayer%bf_compute_used%grdpts_id_tmp)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''grdpts_id_tmp failed'')'
             end if

             test_loc = allocated(current_sublayer%bf_compute_used%nodes_tmp)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''nodes_tmp failed'')'
             end if

             test_loc = allocated(current_sublayer%bf_compute_used%time_dev)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''timedev_tmp failed'')'
             end if

             !next sublayer
             current_sublayer => current_sublayer%get_next()

          end do

        end function test_initialize_before_timeInt


        function test_finalize_after_timeInt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_mainlayer_time) :: bf_mainlayer_used
          integer, dimension(2)   :: x_borders_test
          integer, dimension(2)   :: y_borders_test
          integer, dimension(5,7) :: bc_sections_test

          integer                    :: nb_sublayers
          type(bf_sublayer), pointer :: current_sublayer

          integer :: k
          logical :: test_loc


          test_validated = .true.

          
          !input
          call initialize_bf_mainlayer_for_tests(bf_mainlayer_used)
          call bf_mainlayer_used%initialize_before_timeInt()

          x_borders_test = [1,6]
          y_borders_test = [3,6]
          bc_sections_test = reshape((/
     $         NW_edge_type  , 1,1, no_overlap, NS_overlap,
     $         NE_edge_type  , 5,1, no_overlap, NS_overlap,
     $         W_edge_type   , 1,3, 4         , no_overlap,
     $         E_edge_type   , 5,3, 4         , no_overlap,
     $         NW_corner_type, 1,5, no_overlap, no_overlap,
     $         N_edge_type   , 3,5, 4         , no_overlap,
     $         NE_corner_type, 5,5, no_overlap, no_overlap/),
     $         (/5,7/))

          !output
          call bf_mainlayer_used%finalize_after_timeInt()


          !validation
          nb_sublayers = bf_mainlayer_used%get_nb_sublayers()
          current_sublayer => bf_mainlayer_used%get_head_sublayer()

          do k=1, nb_sublayers

             !integration borders
             test_loc = is_int_vector_validated(
     $            current_sublayer%x_borders,
     $            x_borders_test,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''x_borders failed'')'
             end if

             test_loc = is_int_vector_validated(
     $            current_sublayer%y_borders,
     $            y_borders_test,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''y_borders failed'')'
             end if

             !bc_sections
             test_loc = is_int_matrix_validated(
     $            current_sublayer%bc_sections,
     $            bc_sections_test,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''bc_sections failed'')'
             end if

             !arrays for time integration deallocated
             test_loc = .not.allocated(current_sublayer%bf_compute_used%alignment_tmp)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''alignment_tmp failed'')'
             end if

             test_loc = .not.allocated(current_sublayer%bf_compute_used%grdpts_id_tmp)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''grdpts_id_tmp failed'')'
             end if

             test_loc = .not.allocated(current_sublayer%bf_compute_used%nodes_tmp)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''nodes_tmp failed'')'
             end if

             test_loc = .not.allocated(current_sublayer%bf_compute_used%time_dev)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''timedev_tmp failed'')'
             end if

             !next sublayer
             current_sublayer => current_sublayer%get_next()

          end do

        end function test_finalize_after_timeInt


        function test_compute_time_dev(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_mainlayer_time)             :: bf_mainlayer_used
          type(td_operators)                  :: td_operators_used
          type(sd_operators)                  :: s
          type(pmodel_eq)                     :: p_model
          type(bc_operators)                  :: bc_used
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
          
          
          real(rkind), dimension(6,6,ne) :: timedev_test
          real(rkind), dimension(6,6,ne) :: nodes_test

          integer                    :: nb_sublayers
          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: k
          logical                    :: test_loc


          test_validated = .true.

          
          !input
          call initialize_bf_mainlayer_for_tests(bf_mainlayer_used)
          call bf_mainlayer_used%initialize_before_timeInt()

          !output
          call bf_mainlayer_used%compute_time_dev(
     $         td_operators_used,
     $         0.0d0,s,p_model,bc_used,
     $         interior_nodes)

          !validation
          call get_arrays_for_time_integration_tests(
     $         timedev_test,
     $         nodes_test)

          nb_sublayers = bf_mainlayer_used%get_nb_sublayers()
          current_sublayer => bf_mainlayer_used%get_head_sublayer()

          do k=1,nb_sublayers

             test_loc = is_real_matrix3D_validated(
     $            current_sublayer%bf_compute_used%time_dev(1:6,3:6,:),
     $            timedev_test(1:6,3:6,:),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test timedev failed'')'
             end if

             !next sublayer
             current_sublayer => current_sublayer%get_next()

          end do          

        end function test_compute_time_dev


        function test_compute_integration_step(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_mainlayer_time)             :: bf_mainlayer_used
          type(td_operators)                  :: td_operators_used
          type(sd_operators)                  :: s
          type(pmodel_eq)                     :: p_model
          type(bc_operators)                  :: bc_used
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
          
          
          real(rkind), dimension(6,6,ne) :: timedev_test
          real(rkind), dimension(6,6,ne) :: nodes_test

          integer                    :: nb_sublayers
          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: k
          logical                    :: test_loc


          test_validated = .true.

          
          !input
          call initialize_bf_mainlayer_for_tests(bf_mainlayer_used)
          call bf_mainlayer_used%initialize_before_timeInt()

          !output
          call bf_mainlayer_used%compute_time_dev(
     $         td_operators_used,
     $         0.0d0,s,p_model,bc_used,
     $         interior_nodes)

          call bf_mainlayer_used%compute_integration_step(
     $         0.1d0,compute_1st_step_nopt)

          !validation
          call get_arrays_for_time_integration_tests(
     $         timedev_test,
     $         nodes_test)

          nb_sublayers = bf_mainlayer_used%get_nb_sublayers()
          current_sublayer => bf_mainlayer_used%get_head_sublayer()

          do k=1,nb_sublayers

             test_loc = is_real_matrix3D_validated(
     $            current_sublayer%nodes(1:6,3:6,:),
     $            nodes_test(1:6,3:6,:),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test timedev failed'')'
             end if

             !next sublayer
             current_sublayer => current_sublayer%get_next()

          end do          

        end function test_compute_integration_step


        subroutine initialize_bf_mainlayer_for_tests(bf_mainlayer_used)

          implicit none

          type(bf_mainlayer_time), intent(inout) :: bf_mainlayer_used


          integer, dimension(2,2,3)        :: bf_alignment
          integer    , dimension(6,6)      :: bf_grdpts_id
          real(rkind), dimension(6)        :: bf_x_map
          real(rkind), dimension(6)        :: bf_y_map
          real(rkind), dimension(6,6,ne)   :: bf_nodes

          type(bf_sublayer), pointer       :: added_sublayer
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          integer :: i,j,k


          bf_alignment(:,:,1) = reshape((/
     $         align_W-1, align_N, align_W, align_N+1/),
     $         (/2,2/))
          
          bf_alignment(:,:,2) = reshape((/
     $         align_W+9, align_N, align_W+10, align_N+1/),
     $         (/2,2/))

          bf_alignment(:,:,3) = reshape((/
     $         align_E, align_N, align_E+1, align_N+1/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         2,2,1,1,2,2,
     $         3,2,1,1,2,3,
     $         3,2,1,1,2,3,
     $         3,2,1,1,2,3,
     $         3,2,2,2,2,3,
     $         3,3,3,3,3,3/),
     $         (/6,6/))

          bf_x_map = (/ (x_min+(i-1)*0.1d0,i=1,6) /)
          bf_y_map = (/ (y_min+(j-1)*0.2d0,j=1,6) /)

          bf_nodes = reshape((/
     $         1.48d0, 1.30d0, 1.35d0, 1.31d0, 1.43d0, 1.31d0,
     $         1.26d0, 1.45d0,  1.4d0, 1.29d0, 1.37d0, 1.41d0,
     $         1.46d0, 1.27d0, 1.47d0, 1.28d0, 1.25d0, 1.43d0,
     $         1.48d0, 1.26d0, 1.41d0, 1.34d0, 1.31d0, 1.39d0,
     $         1.42d0, 1.46d0, 1.38d0, 1.26d0, 1.37d0, 1.33d0,
     $         1.41d0, 1.22d0, 1.42d0, 1.23d0, 1.21d0, 1.40d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0, 0.129d0, 0.136d0, 0.124d0,
     $         1.138d0, 0.148d0, 0.132d0, 0.125d0, 0.175d0, 0.123d0,
     $         0.146d0, 0.143d0, 0.145d0, 0.182d0, 0.135d0, 0.154d0,
     $         0.123d0, 0.129d0, 0.124d0, 0.162d0, 0.152d0, 0.142d0,
     $         0.168d0, 0.198d0, 0.186d0, 0.163d0, 0.126d0, 0.168d0,
     $         0.164d0, 0.134d0, 0.154d0, 0.128d0, 0.153d0, 0.145d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0, 0.056d0, 0.062d0, 0.062d0,
     $         0.0025d0, 0.001d0, 0.015d0,  0.07d0, 0.085d0, 0.011d0,
     $         0.0100d0, 0.002d0, 0.050d0,  0.08d0, 0.015d0, 0.057d0,
     $           0.08d0, 0.015d0,  0.09d0, 0.065d0, 0.042d0, 0.067d0,
     $          0.026d0,  0.03d0, 0.045d0, 0.052d0, 0.023d0, 0.051d0,
     $           0.02d0, 0.012d0, 0.098d0, 0.056d0, 0.024d0, 0.090d0,
     $         
     $         4.88d0, 4.870d0,	4.855d0, 4.834d0, 4.592d0, 4.834d0,
     $         4.85d0, 4.865d0, 4.845d0, 4.875d0, 4.815d0, 4.875d0,
     $         4.89d0, 4.870d0, 4.860d0, 4.826d0, 4.723d0, 4.826d0,
     $         4.83d0,  4.95d0,  4.62d0, 4.952d0, 4.852d0, 4.952d0,
     $         4.81d0, 4.758d0, 4.762d0,  4.95d0, 4.703d0, 4.95d0,
     $         4.98d0, 4.780d0, 4.608d0, 4.628d0, 4.237d0, 4.862d0
     $         /),
     $         (/6,6,ne/))


          call bf_mainlayer_used%ini(N)

          do k=1,3

             added_sublayer => bf_mainlayer_used%add_sublayer(
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            bf_alignment(:,:,k))
             added_sublayer%grdpts_id = bf_grdpts_id
             added_sublayer%x_map     = bf_x_map
             added_sublayer%y_map     = bf_y_map
             added_sublayer%nodes     = bf_nodes

          end do

        end subroutine initialize_bf_mainlayer_for_tests


        subroutine get_arrays_for_time_integration_tests(
     $     timedev,
     $     nodes)

          implicit none

          real(rkind), dimension(6,6,ne), intent(out) :: timedev
          real(rkind), dimension(6,6,ne), intent(out) :: nodes

          timedev = reshape((/
     $ -1.81324100585423d0,  0.60371361535018d0,  0.22725509073847d0, -0.01516507648066d0, -0.39267644142951d0,  0.78273426206585d0,
     $  5.03134472386253d0, -0.67003318244995d0,  0.27303140728232d0, -0.29994618105248d0, -0.67263996432621d0, -0.13101865506719d0,
     $ -1.15811574067973d0,  1.12311177707608d0, -0.38250000000000d0,  0.06250000000000d0,  0.60052276187758d0, -1.18870018265479d0,
     $ -1.02864307661362d0,  0.19700429558412d0, -0.15250000000000d0, -0.07000000000000d0,  0.31560947574760d0, -0.45412605745646d0,
     $  0.40805927971097d0, -0.80088630420635d0,  0.27048271888657d0,  0.50154931441799d0,  0.09766264811703d0, -0.15079486587804d0,
     $ -0.98264669761148d0,  1.24806263360275d0, -0.11611849274841d0,  0.29174417088895d0,  1.28936995304316d0, -1.70510436256785d0,
     $         
     $  2.60472632401217d0, -0.26924356708494d0,  0.05971427344415d0, -0.11973434160404d0, -0.61373581128368d0,  0.72390169362333d0,
     $ -3.91593075166454d0,  0.63670773311468d0,  2.08405901287141d0,  0.49456430227676d0, -1.28471481501711d0, -0.08091581906056d0,
     $  2.57644819808738d0, -2.52537309130274d0,  0.03811528724298d0,  3.22542982375567d0,  1.14048749664020d0, -2.59353316892192d0,
     $  2.69510728870803d0, -0.59404880070705d0, -0.77044485062662d0,  0.33599455400759d0,  0.81419441553750d0, -1.15339032734526d0,
     $ -0.09882904164233d0,  1.05192430004862d0,  1.82111769857176d0,  0.40566524349533d0,  0.43429241960583d0, -1.04478770117858d0,
     $  2.28349148178139d0, -1.17695976774214d0,  0.46213584285444d0,  3.28072840906496d0,  1.38626548143010d0, -3.80544677450882d0,
     $         
     $  2.35783121834333d0, -1.08165923763885d0, -0.59316468804325d0,  0.10792952200530d0,  0.17456794794897d0, -0.94780911241018d0,
     $ -2.15695569716296d0,  1.12365527686775d0, -0.43682383404828d0,  0.14339392863515d0,  0.23186352437052d0,  0.30223839343039d0,
     $ -2.30452326972487d0,  1.01189941403061d0,  0.28313451475454d0, -0.37048886481365d0,  0.34516507684081d0, -0.08294216734036d0,
     $  0.44664820279739d0, -0.98827841199167d0,  0.89046763230823d0, -0.12115348134101d0, -0.49817545660221d0,  0.40446845636165d0,
     $  0.93089499666899d0, -0.99470562254562d0,  0.21504753748582d0,  0.46310241294008d0,  0.09499466502237d0,  0.42812640862865d0,
     $ -0.22076183646944d0,  1.39570724539266d0, -0.37993357715650d0,  0.67965693329700d0,  1.52719765921909d0, -0.70591953134663d0,
     $                                                                                                
     $ -8.75180775755818d0,  3.11460294375328d0,  1.16377421960177d0,  0.17143007719475d0, -1.01919413772504d0,  3.33173243062430d0,
     $ 17.76022596044600d0, -3.25493190362516d0,  1.18639803028245d0, -1.52876924636906d0, -2.85520692644194d0, -0.53358624558689d0,
     $ -5.56960164311805d0,  5.69410602387996d0, -1.98776588829396d0,  0.20058149836543d0,  3.16783535205124d0, -5.02996274073871d0,
     $ -4.81560124148740d0,  1.03793539026216d0, -0.35265170025370d0, -0.89256825113173d0,  1.65862665333175d0, -2.14020075281347d0,
     $  1.98942665472833d0, -3.57627077636379d0,  0.86586517339179d0,  2.36820969873451d0,  1.33459499857288d0, -1.52420504277944d0,
     $ -4.80513012738320d0,  6.51768123603877d0, -0.25428019996007d0,  1.41053382625555d0,  6.55766680279622d0, -8.15016650575661d0
     $         /),
     $         (/6,6,ne/))

          nodes=reshape((/
     $  1.29867589941458d0, 1.36037136153502d0,  1.37272550907385d0, 1.30848349235193d0, 1.39073235585705d0,  1.38827342620659d0,
     $  1.76313447238625d0, 1.38299668175500d0,  1.42730314072823d0, 1.26000538189475d0, 1.30273600356738d0,  1.39689813449328d0,
     $  1.34418842593203d0, 1.38231117770761d0,  1.43175000000000d0, 1.28625000000000d0, 1.31005227618776d0,  1.31112998173452d0,
     $  1.37713569233864d0, 1.27970042955841d0,  1.39475000000000d0, 1.33300000000000d0, 1.34156094757476d0,  1.34458739425435d0,
     $  1.46080592797110d0, 1.37991136957937d0,  1.40704827188866d0, 1.31015493144180d0, 1.37976626481170d0,  1.31492051341220d0,
     $  1.31173533023885d0, 1.34480626336027d0,  1.40838815072516d0, 1.25917441708890d0, 1.33893699530432d0,  1.22948956374321d0,
     $                                                                                                                        
     $  0.38847263240122d0, 0.10007564329151d0,  0.14797142734441d0, 0.11702656583960d0, 0.07462641887163d0,  0.19639016936233d0,
     $  0.74640692483355d0, 0.21167077331147d0,  0.34040590128714d0, 0.17445643022768d0, 0.04652851849829d0,  0.11490841809394d0,
     $  0.40364481980874d0,-0.10953730913027d0,  0.14881152872430d0, 0.50454298237557d0, 0.24904874966402d0, -0.10535331689219d0,
     $  0.39251072887080d0, 0.06959511992930d0,  0.04695551493734d0, 0.19559945540076d0, 0.23341944155375d0,  0.02666096726547d0,
     $  0.15811709583577d0, 0.30319243000486d0,  0.36811176985718d0, 0.20356652434953d0, 0.16942924196058d0,  0.06352122988214d0,
     $  0.39234914817814d0, 0.01630402322579d0,  0.20021358428544d0, 0.45607284090650d0, 0.29162654814301d0, -0.23554467745088d0,
     $                                                                                                                        
     $  0.24078312183433d0,-0.08816592376389d0,  0.00068353119568d0, 0.06679295220053d0, 0.07945679479490d0, -0.03278091124102d0,
     $ -0.21319556971630d0, 0.11336552768678d0, -0.02868238340483d0, 0.08433939286351d0, 0.10818635243705d0,  0.04122383934304d0,
     $ -0.22045232697249d0, 0.10318994140306d0,  0.07831345147545d0, 0.04295111351863d0, 0.04951650768408d0,  0.04870578326596d0,
     $  0.12466482027974d0,-0.08382784119917d0,  0.17904676323082d0, 0.05288465186590d0,-0.00781754566022d0,  0.10744684563617d0,
     $  0.11908949966690d0,-0.06947056225456d0,  0.06650475374858d0, 0.09831024129401d0, 0.03249946650224d0,  0.09381264086287d0,
     $ -0.00207618364694d0, 0.15157072453927d0,  0.06000664228435d0, 0.12396569332970d0, 0.17671976592191d0,  0.01940804686534d0,
     $                                                                                                                        
     $  4.00481922424418d0, 5.18146029437533d0,  4.97137742196018d0, 4.85114300771947d0, 4.49008058622750d0,  5.16717324306243d0,
     $  6.62602259604460d0, 4.53950680963748d0,  4.96363980302824d0, 4.72212307536309d0, 4.52947930735581d0,  4.82164137544131d0,
     $  4.33303983568819d0, 5.43941060238800d0,  4.66122341117060d0, 4.84605814983654d0, 5.03978353520512d0,  4.32300372592613d0,
     $  4.34843987585126d0, 5.05379353902622d0,  4.58473482997463d0, 4.86274317488683d0, 5.01786266533318d0,  4.73797992471865d0,
     $  5.00894266547283d0, 4.40037292236362d0,  4.84858651733918d0, 5.18682096987345d0, 4.83645949985729d0,  4.79757949572206d0,
     $  4.49948698726168d0, 5.43176812360388d0,  4.58257198000399d0, 4.76905338262556d0, 4.89276668027962d0,  4.04698334942434d0
     $         /),
     $         (/6,6,ne/))

        end subroutine get_arrays_for_time_integration_tests

      
        subroutine check_inputs()

          implicit none

          if(.not.(
     $         (nx.eq.20).and.
     $         (ny.eq.25).and.
     $         (ne.eq.4))) then

             print '(''the test requires: '')'
             print '(''nx=20'')'
             print '(''ny=25'')'
             print '(''ne=4'')'
             print '(''use DIM2D model'')'
             print '(''use hedstrom b.c.'')'
             print '(''use newgrdpt_test i.c.'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_bf_mainlayer_time

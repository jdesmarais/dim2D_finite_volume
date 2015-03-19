      program test_bf_increase_coords

        use bf_increase_coords_module, only :
     $       get_mainlayer_id

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W

        use parameters_constant, only :
     $       N,S,E,W,interior

        use parameters_kind, only :
     $       ikind

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        test_loc = test_get_mainlayer_id(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_mainlayer_id: '',L1)', test_loc
        print '()'


        contains


        function test_get_mainlayer_id(detailled)
     $       result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          integer       , parameter             :: nb_tests = 9
          integer(ikind), dimension(2,nb_tests) :: test_gen_coords
          integer       , dimension(nb_tests)   :: test_mainlayer_id
          integer                               :: mainlayer_id

          integer :: k
          logical :: test_loc


          test_validated = .true.


          !input
          test_gen_coords(:,1) = [align_W            , align_S]
          test_gen_coords(:,2) = [(align_W+align_E)/2, align_S]
          test_gen_coords(:,3) = [align_E            , align_S]

          test_gen_coords(:,4) = [align_W            , (align_S+align_N)/2]
          test_gen_coords(:,5) = [(align_W+align_E)/2, (align_S+align_N)/2]
          test_gen_coords(:,6) = [align_E            , (align_S+align_N)/2]

          test_gen_coords(:,7) = [align_W            , align_N]
          test_gen_coords(:,8) = [(align_W+align_E)/2, align_N]
          test_gen_coords(:,9) = [align_E            , align_N]


          test_mainlayer_id(1:3) = [S,S,S]
          test_mainlayer_id(4:6) = [W,interior,E]
          test_mainlayer_id(7:9) = [N,N,N]


          do k=1, nb_tests

             !output
             mainlayer_id = get_mainlayer_id(test_gen_coords(:,k))

             !validation
             test_loc = mainlayer_id.eq.test_mainlayer_id(k)
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_get_mainlayer_id

      end program test_bf_increase_coords

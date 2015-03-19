      program test_bf_mainlayer

        use bf_mainlayer_class, only :
     $       bf_mainlayer

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_bf_layer, only :
     $       align_N,
     $       align_E,
     $       align_W

        use parameters_constant, only :
     $       N

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_get_sublayer_from_gen_coords(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_sublayer_from_gen_coords: '',L1)', test_loc
        print '()'


        contains


        function test_get_sublayer_from_gen_coords(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_mainlayer)         :: bf_mainlayer_used
          type(bf_sublayer), pointer :: bf_sublayer1
          type(bf_sublayer), pointer :: bf_sublayer2
          type(bf_sublayer), pointer :: bf_sublayer_match

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          integer       , parameter             :: nb_tests=8
          integer(ikind), dimension(2,nb_tests) :: gen_coords
          integer       , dimension(nb_tests)   :: tolerance
          integer       , dimension(nb_tests)   :: bf_sublayer_match_i

          integer :: k


          test_validated = .true.


          !input
          call bf_mainlayer_used%ini(N)

          bf_sublayer1 => bf_mainlayer_used%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $         align_W,align_N,align_W+10,align_N/),
     $         (/2,2/)))

          bf_sublayer2 => bf_mainlayer_used%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $         align_E-5,align_N,align_E,align_N/),
     $         (/2,2/)))

          
          gen_coords(:,1)        = [align_W-3,align_N]
          tolerance(1)           = 0
          bf_sublayer_match_i(1) = 0

          gen_coords(:,2)        = [align_W-3,align_N]
          tolerance(2)           = 1
          bf_sublayer_match_i(2) = 1

          gen_coords(:,3)        = [align_W+13,align_N]
          tolerance(3)           = 0
          bf_sublayer_match_i(3) = 0

          gen_coords(:,4)        = [align_W+13,align_N]
          tolerance(4)           = 1
          bf_sublayer_match_i(4) = 1


          gen_coords(:,5)        = [align_E-8,align_N]
          tolerance(5)           = 0
          bf_sublayer_match_i(5) = 0

          gen_coords(:,6)        = [align_E-8,align_N]
          tolerance(6)           = 1
          bf_sublayer_match_i(6) = 2

          gen_coords(:,7)        = [align_E+3,align_N]
          tolerance(7)           = 0
          bf_sublayer_match_i(7) = 0

          gen_coords(:,8)        = [align_E+3,align_N]
          tolerance(8)           = 1
          bf_sublayer_match_i(8) = 2


          do k=1,nb_tests


             !output
             bf_sublayer_match => bf_mainlayer_used%get_sublayer_from_gen_coords(
     $            gen_coords(:,k),
     $            tolerance(k))


             !validation
             select case(bf_sublayer_match_i(k))
               case(0)
                  test_loc = .not.associated(bf_sublayer_match)
               case(1)
                  test_loc = associated(bf_sublayer_match,bf_sublayer1)
               case(2)
                  test_loc = associated(bf_sublayer_match,bf_sublayer2)
             end select
             test_validated = test_validated.and.test_loc


             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if


          end do

        end function test_get_sublayer_from_gen_coords


        subroutine check_inputs()

          implicit none

          if(nx.le.24) then
             print '(''the testt requires:'')'
             print '(''nx>24'')'
             stop ''
          end if

        end subroutine check_inputs


      end program test_bf_mainlayer

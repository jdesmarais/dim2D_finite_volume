      program test_bf_interface_coords

        use bf_interface_coords_class, only :
     $       bf_interface_coords

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W

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


        call check_inputs()


        detailled = .true.
        test_validated = .true.


        test_loc = test_get_bf_layer_from_gen_coords(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_bf_layer_from_gen_coords: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_get_bf_layer_from_gen_coords(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_coords)  :: bf_interface_used

          type(bf_sublayer), pointer :: bf_sublayerN1
          type(bf_sublayer), pointer :: bf_sublayerN2
         
          type(bf_sublayer), pointer :: bf_sublayer_match

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          integer       , parameter             :: nb_tests=8
          integer(ikind), dimension(2,nb_tests) :: gen_coords
          integer       , dimension(nb_tests)   :: tolerance
          integer       , dimension(nb_tests)   :: bf_sublayer_match_i

          integer(ikind), dimension(2,2) :: alignment_tmp

          integer :: k


          test_validated = .true.


          !input
          call bf_interface_used%ini(interior_x_map,interior_y_map)

          alignment_tmp = reshape((/
     $         align_W,align_N,align_W+10,align_N/),
     $         (/2,2/))

          bf_sublayerN1 => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         alignment_tmp)

          alignment_tmp = reshape((/
     $         align_E-5,align_N,align_E,align_N/),
     $         (/2,2/))

          bf_sublayerN2 => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         alignment_tmp)

          
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
             bf_sublayer_match => bf_interface_used%get_bf_layer_from_gen_coords(
     $            gen_coords(:,k),
     $            tolerance=tolerance(k))


             !validation
             select case(bf_sublayer_match_i(k))
               case(0)
                  test_loc = .not.associated(bf_sublayer_match)
               case(1)
                  test_loc = associated(bf_sublayer_match,bf_sublayerN1)
               case(2)
                  test_loc = associated(bf_sublayer_match,bf_sublayerN2)
             end select
             test_validated = test_validated.and.test_loc


             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if


          end do

        end function test_get_bf_layer_from_gen_coords


        subroutine check_inputs()

          implicit none

          if(nx.le.24) then
             print '(''the test requires:'')'
             print '(''nx>24'')'
             stop ''
          end if

        end subroutine check_inputs

      end program test_bf_interface_coords

      program test_mpi_tag

        use mpi_tag_module, only :
     $     compute_mpi_tag

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        
        detailled = .true.
        test_validated = .true.


        test_loc = test_compute_mpi_tag(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_mpi_tag: '',L1)', test_loc
        print '()'

        print '(''test_validated: '',L1)', test_validated


        contains


        function test_compute_mpi_tag(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, parameter           :: nb_tests=4
          integer                      :: k
          logical                      :: test_loc
          integer                      :: mpi_tag

          integer, dimension(nb_tests) :: rank_send
          integer, dimension(nb_tests) :: rank_recv
          integer, dimension(nb_tests) :: nb_procs
          integer, dimension(nb_tests) :: test_mpi_tag


          test_validated = .true.


          rank_send = [0,1,2,3]
          rank_recv = [7,2,5,4]
          nb_procs  = [2,4,6,7]

          test_mpi_tag = [7,6,17,25]


          do k=1, nb_tests

             mpi_tag = compute_mpi_tag(
     $            rank_send(k),
     $            rank_recv(k),
     $            nb_procs(k))

             test_loc = mpi_tag.eq.test_mpi_tag(k)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do


        end function test_compute_mpi_tag

      end program test_mpi_tag

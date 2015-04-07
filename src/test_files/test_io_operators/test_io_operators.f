      program test_io_operators

        use check_data_module, only :
     $     is_char_validated

        use io_operators_module, only :
     $       get_filename,
     $       get_timestep

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        test_loc = test_get_filename(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_filename: '',L1)', test_loc
        print '()'


        test_loc = test_get_timestep(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_timestep: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated


        contains


        function test_get_filename(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated



          character(len=16)               :: filename
          character(len=16), dimension(4) :: test_filename
          integer          , dimension(4) :: test_timestep
          integer          , dimension(2) :: test_rank


          integer :: k
          integer :: test_loc

          
          test_validated = .true.


          test_filename = [
     $         'data0.nc',
     $         'data10.nc',
     $         'data0_0.nc',
     $         'data10_20.nc']

          test_timestep = [0,10,0,10]
          test_rank     = [0,20]          


          do k=1,2

             call get_filename(filename, test_timestep(k))

             test_loc = is_char_validated(
     $            filename,
     $            test_filename(k),
     $            detailled)

             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do


          do k=3,4

             call get_filename(filename, test_timestep(k),test_rank(k-2))

             test_loc = is_char_validated(
     $            filename,
     $            test_filename(k),
     $            detailled)

             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_get_filename


        function test_get_timestep(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated



          integer                         :: timestep
          character(len=16), dimension(3) :: test_filename
          integer          , dimension(3) :: test_timestep


          integer :: k
          integer :: test_loc

          
          test_validated = .true.


          test_filename = [
     $         'data0.nc',
     $         'data10.nc',
     $         'data1000.nc']

          test_timestep = [0,10,1000]


          do k=1,3

             timestep = get_timestep(trim(test_filename(k)))

             test_loc = timestep.eq.test_timestep(k)

             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_get_timestep

      end program test_io_operators

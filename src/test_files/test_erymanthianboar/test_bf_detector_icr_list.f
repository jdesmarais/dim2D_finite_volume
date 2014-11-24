      program test_bf_detector_icr_list

        use bf_detector_icr_list_class, only :
     $       bf_detector_icr_list

        use parameters_kind, only :
     $       ikind, rkind

        implicit none

        type(bf_detector_icr_list) :: bf_detectors_used
        logical                    :: test_validated
        logical                    :: detailled


        detailled = .true.

        test_validated = test_ini(bf_detectors_used)
        print '(''test_ini: '',L1)', test_validated
        print '()'

        test_validated = test_add_new_detector(bf_detectors_used,detailled)
        print '(''test_add_new_detector: '', L1)', test_validated
        print '()'

        test_validated = test_fill_new_detector_table(bf_detectors_used,detailled)
        print '(''test_fill_new_detector_table: '',L1)', test_validated
        print '()'

        test_validated = test_destroy(bf_detectors_used)
        print '(''test_destroy: '',L1)', test_validated
        print '()'


        contains

        function test_ini(bf_detectors_used)
     $       result(test_validated)
        
          implicit none

          type(bf_detector_icr_list), intent(inout) :: bf_detectors_used
          logical                                   :: test_validated

          integer :: mainlayer_id
          integer :: size_detectors_list

          mainlayer_id        = 1
          size_detectors_list = 4
          
          call bf_detectors_used%ini(mainlayer_id,size_detectors_list)

          test_validated = 
     $         (bf_detectors_used%mainlayer_id.eq.mainlayer_id).and.
     $         (size(bf_detectors_used%icoords,2).eq.size_detectors_list).and.
     $         (size(bf_detectors_used%rcoords,2).eq.size_detectors_list).and.
     $         (bf_detectors_used%get_nb_detectors().eq.0)

        end function test_ini


        function test_add_new_detector(bf_detectors_used,detailled)
     $       result(test_validated)
        
          implicit none

          type(bf_detector_icr_list), intent(inout) :: bf_detectors_used
          logical                   , intent(in)    :: detailled
          logical                                   :: test_validated

          integer                         :: k
          integer(ikind), dimension(2,6)  :: icoords_input
          real(rkind)   , dimension(2,6)  :: rcoords_input
          integer(ikind), dimension(2,12) :: icoords_test
          real(rkind)   , dimension(2,12) :: rcoords_test

          logical :: test_loc


          call get_test_data(
     $         icoords_input,
     $         rcoords_input,
     $         icoords_test,
     $         rcoords_test)          

          do k=1, size(icoords_input,2)

             call bf_detectors_used%add_new_detector(
     $            icoords_input(:,k),
     $            rcoords_input(:,k))

          end do

          test_validated = .true.

          do k=1, bf_detectors_used%get_nb_detectors()

             test_loc =
     $            icoords_test(1,k).eq.bf_detectors_used%icoords(1,k)

             test_loc = test_loc.and.(
     $            icoords_test(2,k).eq.bf_detectors_used%icoords(2,k))

             test_loc = test_loc.and.is_test_validated(
     $            rcoords_test(1,k),
     $            bf_detectors_used%rcoords(1,k),
     $            detailled)

             test_loc = test_loc.and.is_test_validated(
     $            rcoords_test(2,k),
     $            bf_detectors_used%rcoords(2,k),
     $            detailled)

          
             if(.not.test_loc) then
                print '(''**test failed at '',I2,'' **'')', k
             end if

             test_validated = test_validated.and.test_loc

          end do


          do k=1, bf_detectors_used%get_nb_detectors()

             print '(''icoords: '', 2I2, ''rcoords: '',2F8.2)', 
     $            bf_detectors_used%icoords(:,k),
     $            bf_detectors_used%rcoords(:,k)

          end do
          print '()'
          
        end function test_add_new_detector


        function test_fill_new_detector_table(bf_detectors_used,detailled)
     $       result(test_validated)
        
          implicit none

          type(bf_detector_icr_list), intent(inout) :: bf_detectors_used
          logical                   , intent(in)    :: detailled
          logical                                   :: test_validated

          integer                         :: k
          integer(ikind), dimension(2,6)  :: icoords_input
          real(rkind)   , dimension(2,6)  :: rcoords_input
          integer(ikind), dimension(2,12) :: icoords_test
          real(rkind)   , dimension(2,12) :: rcoords_test
          integer(ikind), dimension(2,13) :: new_icoords
          real(rkind)   , dimension(2,13) :: new_rcoords

          integer :: start_i
          logical :: test_loc


          call get_test_data(
     $         icoords_input,
     $         rcoords_input,
     $         icoords_test,
     $         rcoords_test)

          start_i = 2

          call bf_detectors_used%fill_new_detector_table(
     $         start_i,
     $         new_icoords,
     $         new_rcoords)


          test_validated = .true.

          do k=1, size(rcoords_test,2)

             test_loc = is_test_validated(
     $            rcoords_test(1,k),
     $            new_rcoords(1,start_i+k-1),
     $            detailled)

             test_loc = test_loc.and.is_test_validated(
     $            rcoords_test(2,k),
     $            new_rcoords(2,start_i+k-1),
     $            detailled)

          
             if(.not.test_loc) then
                print '(''**test failed at '',I2,'' **'')', k
             end if

             test_validated = test_validated.and.test_loc

          end do

        end function test_fill_new_detector_table


        function test_destroy(bf_detectors_used)
     $       result(test_validated)
        
          implicit none

          type(bf_detector_icr_list), intent(inout) :: bf_detectors_used
          logical                                   :: test_validated


          call bf_detectors_used%destroy()

          test_validated = 
     $         (.not.(allocated(bf_detectors_used%icoords))).and.
     $         (.not.(allocated(bf_detectors_used%rcoords)))

        end function test_destroy        


        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical    , intent(in) :: detailled
          logical                 :: test_validated

          if(detailled) then
             print *, int(var*1e5)
             print *, int(cst*1e5)
          end if
          
          test_validated=abs(
     $         int(var*1e5)-
     $         int(cst*1e5)).le.1
          
        end function is_test_validated


        subroutine get_test_data(
     $     icoords_input,
     $     rcoords_input,
     $     icoords_test,
     $     rcoords_test)

          implicit none

          integer(ikind), dimension(:,:), intent(inout) :: icoords_input
          real(rkind)   , dimension(:,:), intent(inout) :: rcoords_input
          integer(ikind), dimension(:,:), intent(inout) :: icoords_test
          real(rkind)   , dimension(:,:), intent(inout) :: rcoords_test


          icoords_input(:,1) = [1,1]
          rcoords_input(:,1) = [0.0d0,0.0d0]

          icoords_input(:,2) = [4,1]
          rcoords_input(:,2) = [0.75d0,0.0d0]

          icoords_input(:,3) = [4,3]
          rcoords_input(:,3) = [0.75d0,1.0d0]

          icoords_input(:,4) = [5,3]
          rcoords_input(:,4) = [1.0d0,1.0d0]

          icoords_input(:,5) = [5,3]
          rcoords_input(:,5) = [1.0d0,1.0d0]

          icoords_input(:,6) = [10,1]
          rcoords_input(:,6) = [2.0d0,0.0d0]


          icoords_test(:,1)  = [1,1]
          rcoords_test(:,1)  = [0.0d0,0.0d0]
                            
          icoords_test(:,2)  = [2,1]
          rcoords_test(:,2)  = [0.25d0,0.0d0]
          icoords_test(:,3)  = [3,1]
          rcoords_test(:,3)  = [0.5d0,0.0d0]
          icoords_test(:,4)  = [4,1]
          rcoords_test(:,4)  = [0.75d0,0.0d0]
                            
          icoords_test(:,5)  = [4,2]
          rcoords_test(:,5)  = [0.75d0,0.5d0]
          icoords_test(:,6)  = [4,3]
          rcoords_test(:,6)  = [0.75d0,1.0d0]
          
          icoords_test(:,7)  = [5,3]
          rcoords_test(:,7)  = [1.0d0, 1.0d0]
          icoords_test(:,8)  = [6,3]
          rcoords_test(:,8)  = [1.20d0,0.80d0]
          icoords_test(:,9)  = [7,2]
          rcoords_test(:,9)  = [1.40d0,0.60d0]
          icoords_test(:,10) = [8,2]
          rcoords_test(:,10) = [1.60d0,0.40d0]
          icoords_test(:,11) = [9,2]
          rcoords_test(:,11) = [1.80d0,0.20d0]
          icoords_test(:,12) = [10,1]
          rcoords_test(:,12) = [2.00d0,0.00d0]

        end subroutine get_test_data

      end program test_bf_detector_icr_list

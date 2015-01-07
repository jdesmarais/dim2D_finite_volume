      program test_bf_detector

        use bf_detector_module, only :
     $       determine_local_map_coordinates,
     $       get_inter_detector_param

        use check_data_module, only :
     $       is_vector_validated

        use parameters_input, only :
     $       nx,ny

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled      = .true.
        test_validated = .true.

        call test_inputs()

        test_loc = test_determine_local_map_coordinates(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_determine_local_map_coordinates: '',L1)', test_loc
        print '()'

        test_loc = test_get_inter_detector_param(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_inter_detector_param: '',L1)', test_loc
        print '()'
        
        contains


        function test_get_inter_detector_param(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer(ikind), dimension(2)              :: prev_icoord
          integer(ikind), dimension(2)              :: next_icoord
          real(rkind)   , dimension(nx)             :: interior_x_map
          real(rkind)   , dimension(ny)             :: interior_y_map

          real(rkind)   , dimension(2)              :: icoord_icr
          integer                                   :: inter_nb
          real(rkind)   , dimension(:), allocatable :: x_map_icr
          real(rkind)   , dimension(:), allocatable :: y_map_icr

          real(rkind)   , dimension(2)              :: test_icoord_icr
          integer                                   :: test_inter_nb
          real(rkind)   , dimension(:), allocatable :: test_x_map_icr
          real(rkind)   , dimension(:), allocatable :: test_y_map_icr

          integer :: test_id
          logical :: validated_icoord_icr
          logical :: validated_inter_nb
          logical :: validated_x_map_icr
          logical :: validated_y_map_icr

          interior_x_map = [0.1d0,0.2d0,0.3d0,0.4d0,0.45d0]
          interior_y_map = [0.0d0,0.2d0,0.4d0,0.6d0,0.8d0,0.9d0]
          
          do test_id=1,1

             call make_test_get_inter_detector_param(
     $            test_id,
     $            prev_icoord,
     $            next_icoord,
     $            test_icoord_icr,
     $            test_inter_nb,
     $            test_x_map_icr,
     $            test_y_map_icr)

             call get_inter_detector_param(
     $            prev_icoord,
     $            next_icoord,
     $            interior_x_map,
     $            interior_y_map,
     $            icoord_icr,
     $            inter_nb,
     $            x_map_icr,
     $            y_map_icr)

             validated_icoord_icr = is_vector_validated(
     $            icoord_icr,
     $            test_icoord_icr,
     $            detailled)

             validated_inter_nb =
     $            test_inter_nb.eq.inter_nb

             validated_x_map_icr = is_vector_validated(
     $            x_map_icr,
     $            test_x_map_icr,
     $            detailled)

             validated_y_map_icr = is_vector_validated(
     $            y_map_icr,
     $            test_y_map_icr,
     $            detailled)

             test_validated =
     $            validated_icoord_icr.and.
     $            validated_inter_nb.and.
     $            validated_x_map_icr.and.
     $            validated_y_map_icr

             if(detailled) then

                print '(''icoord_icr: '')', validated_icoord_icr
                if(.not.validated_icoord_icr) then
                   print '(''   - icoord_icr(1): '',I2,'' -> '',I2)',
     $                  icoord_icr(1), test_icoord_icr(1)
                   print '(''   - icoord_icr(2): '',I2,'' -> '',I2)',
     $                  icoord_icr(2), test_icoord_icr(2)
                   print '()'
                end if

                print '(''inter_nb: '')', validated_inter_nb
                if(.not.validated_inter_nb) then
                   print '(''   - '',I2,'' -> '',I2)',
     $                  inter_nb, test_inter_nb
                   print '()'
                end if

                print '(''x_map_icr: '')', validated_x_map_icr
                if(.not.validated_x_map_icr) then
                   print '('' -  x_map_icr: '')'
                   print *, x_map_icr
                   print '(''    -> '')'
                   print '('' -  test_x_map_icr: '')'
                   print *, test_x_map_icr
                   print '()'
                end if

                print '(''y_map_icr: '')', validated_y_map_icr
                if(.not.validated_y_map_icr) then
                   print '('' -  y_map_icr: '')'
                   print *, y_map_icr
                   print '(''    -> '')'
                   print '('' -  test_y_map_icr: '')'
                   print *, test_y_map_icr
                   print '()'
                end if

             end if

             deallocate(test_x_map_icr)
             deallocate(test_y_map_icr)
             deallocate(x_map_icr)
             deallocate(y_map_icr)

          end do

        end function test_get_inter_detector_param


        subroutine make_test_get_inter_detector_param(
     $     test_id,
     $     prev_icoord,
     $     next_icoord,
     $     test_icoord_icr,
     $     test_inter_nb,
     $     test_x_map_icr,
     $     test_y_map_icr)

          implicit none

          integer                                  , intent(in)  :: test_id
          integer(ikind), dimension(2)             , intent(out) :: prev_icoord
          integer(ikind), dimension(2)             , intent(out) :: next_icoord
          real(rkind)   , dimension(2)             , intent(out) :: test_icoord_icr
          integer                                  , intent(out) :: test_inter_nb
          real(rkind)   , dimension(:), allocatable, intent(out) :: test_x_map_icr
          real(rkind)   , dimension(:), allocatable, intent(out) :: test_y_map_icr


          select case(test_id)
            case(1)
               prev_icoord(1) = 1
               prev_icoord(2) = 1

               next_icoord(1) = 4
               next_icoord(2) = 1

               test_icoord_icr(1) = 1.0d0
               test_icoord_icr(2) = 0.0d0

               test_inter_nb = 2

               allocate(test_x_map_icr(3))
               allocate(test_y_map_icr(1))

               test_x_map_icr = [0.1d0,0.2d0,0.3d0]
               test_y_map_icr = [0.0d0]

            case default
               print '(''make_test_get_inter_detector_param'')'
               print '(''test not implemented: '',I2)', test_id
               stop ''
               
          end select


        end subroutine make_test_get_inter_detector_param

        
        function test_determine_local_map_coordinates(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer                      :: size_interior_map
          integer                      :: size_local_map
          real(rkind)   , dimension(5) :: interior_map
          real(rkind)   , dimension(3) :: local_map
          real(rkind)   , dimension(3) :: test_local_map
          integer(ikind)               :: first_icoord

          integer :: test_id

          test_validated = .true.

          size_interior_map = size(interior_map,1)
          interior_map = [0.1d0,0.2d0,0.3d0,0.4d0,0.45d0]

          size_local_map = 3

          do test_id=1,4

             !get the test presets
             call make_test_determine_local_map_coordinates(
     $            test_id, first_icoord, test_local_map)

             !compute the local map using the tested function
             call determine_local_map_coordinates(
     $            interior_map,
     $            size_interior_map,
     $            size_local_map,
     $            first_icoord,
     $            local_map)

             !compare the results
             test_loc = is_vector_validated(
     $            local_map,
     $            test_local_map,
     $            detailled)
             test_validated = test_validated.and.test_loc

             !print the potential differences if detailled
             if((.not.test_loc).and.detailled) then
                print '(''**test '',I2,'' failed**'')', test_id
                print '(''local_map     : '',3F8.3)', local_map
                print '(''test_local_map: '',3F8.3)', test_local_map
                print ''                
             end if

          end do          

        end function test_determine_local_map_coordinates        

        
        subroutine make_test_determine_local_map_coordinates(
     $     test_id, first_icoord, test_local_map)

          implicit none

          integer                  , intent(in)  :: test_id
          integer                  , intent(out) :: first_icoord
          real(rkind), dimension(3), intent(out) :: test_local_map

          select case(test_id)

            case(1)
               first_icoord = -2
               test_local_map = [-0.2d0,-0.1d0,0.0d0]

            case(2)
               first_icoord = 1
               test_local_map = [0.1d0,0.2d0,0.3d0]

            case(3)
               first_icoord = 3
               test_local_map = [0.3d0,0.4d0,0.45d0]

            case(4)
               first_icoord = 6
               test_local_map = [0.5d0,0.55d0,0.6d0]
               
            case default
               print '(''make_test_determine_locla_map_coordinates'')'
               print '(''test case not implemented: '',I2)', test_id
               stop ''

          end select

        end subroutine make_test_determine_local_map_coordinates


        subroutine test_inputs()

          implicit none

          logical :: inputs_validated

          inputs_validated = nx.eq.5
          inputs_validated = inputs_validated.and.(ny.eq.6)

          if(.not.inputs_validated) then
             
             print '(''the test requires:'')'
             print '(''nx=5'')', nx.eq.5
             print '(''ny=6'')', ny.eq.6
             print '()'
             stop ''

          end if

        end subroutine test_inputs

      end program test_bf_detector

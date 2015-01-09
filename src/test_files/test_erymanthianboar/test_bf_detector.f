      program test_bf_detector

        use bf_detector_module, only :
     $       determine_local_map_coordinates,
     $       get_inter_detector_param,
     $       get_inter_detector_coords,
     $     
     $       get_rot_coords,
     $       get_minimum_block_nb,
     $       get_inter_detector_rot_param,
     $       get_inter_detector_rot_coords

        use check_data_module, only :
     $       is_vector_validated,
     $       is_matrix_validated

        use parameters_input, only :
     $       nx,ny

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        test_validated = .true.

        detailled      = .false.

        call test_inputs()


        test_loc = test_determine_local_map_coordinates(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_determine_local_map_coordinates: '',L1)', test_loc
        print '()'

        test_loc = test_get_inter_detector_param(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_inter_detector_param: '',L1)', test_loc
        print '()'

c$$$        test_loc = test_get_inter_detector_param_sym(detailled)
c$$$        test_validated = test_validated.and.test_loc
c$$$        print '(''test_get_inter_detector_param_sym: '',L1)', test_loc
c$$$        print '()'

        test_loc = test_get_inter_detector_coords(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_inter_detector_coords: '',L1)', test_loc
        print '()'

        detailled      = .false.

        test_loc = test_get_rot_coords(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_rot_coords: '',L1)', test_loc
        print '()'
        
        test_loc = test_get_minimum_block_nb(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_minimum_block_nb: '',L1)', test_loc
        print '()'

        test_loc = test_get_inter_detector_rot_param(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_inter_detector_rot_param: '',L1)', test_loc
        print '()'

        test_loc = test_get_inter_detector_rot_coords(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_inter_detector_rot_coords: '',L1)', test_loc
        print '()'

        contains

        function test_get_inter_detector_rot_coords(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer(ikind), dimension(2)                :: prev_icoords
          integer(ikind), dimension(2)                :: next_icoords
          real(rkind)   , dimension(nx)               :: interior_x_map
          real(rkind)   , dimension(ny)               :: interior_y_map
                                                      
          real(rkind)   , dimension(2)                :: rot_icoords_r
          real(rkind)   , dimension(2)                :: icoords_icr
          integer(ikind)                              :: inter_nb
          real(rkind)   , dimension(:)  , allocatable :: x_map_icr
          real(rkind)   , dimension(:)  , allocatable :: y_map_icr

          integer(ikind), dimension(:,:), allocatable :: icoords_inter
          integer(ikind), dimension(:,:), allocatable :: test_icoords_inter
          real(rkind)   , dimension(:,:), allocatable :: rcoords_inter
          real(rkind)   , dimension(:,:), allocatable :: test_rcoords_inter

          integer :: test_id

          logical :: validated_icoords_inter
          logical :: validated_rcoords_inter

          integer :: k

          test_validated = .true.


          interior_x_map = [0.1d0,0.2d0,0.3d0,0.4d0,0.45d0]
          interior_y_map = [0.0d0,0.2d0,0.4d0,0.6d0,0.8d0,0.9d0]

          do test_id=1,24

             call make_test_get_inter_detector_rot_coords(
     $            test_id,
     $            interior_x_map,
     $            interior_y_map,
     $            prev_icoords,
     $            next_icoords,
     $            test_icoords_inter,
     $            test_rcoords_inter)

             call get_inter_detector_rot_param(
     $            prev_icoords,
     $            next_icoords,
     $            interior_x_map,
     $            interior_y_map,
     $            rot_icoords_r,
     $            icoords_icr,
     $            inter_nb,
     $            x_map_icr,
     $            y_map_icr)

             allocate(icoords_inter(2,inter_nb))
             allocate(rcoords_inter(2,inter_nb))             

             do k=1, inter_nb

                call get_inter_detector_rot_coords(
     $               prev_icoords,
     $               rot_icoords_r,
     $               icoords_icr,
     $               inter_nb,
     $               x_map_icr,
     $               y_map_icr,
     $               k,
     $               icoords_inter(:,k),
     $               rcoords_inter(:,k))

             end do

             validated_icoords_inter = .true.

             do k=1,inter_nb
                validated_icoords_inter = 
     $               validated_icoords_inter.and.(
     $               (icoords_inter(1,k).eq.test_icoords_inter(1,k)).and.
     $               (icoords_inter(2,k).eq.test_icoords_inter(2,k)))
             end do
             
             validated_rcoords_inter = is_matrix_validated(
     $            rcoords_inter,
     $            test_rcoords_inter,
     $            detailled)
             
             test_loc =
     $            validated_icoords_inter.and.
     $            validated_rcoords_inter

             test_validated = test_validated.and.test_loc

             if(detailled) then
                print '(''test '',I2)', test_id

                if(.not.validated_icoords_inter) then
                   print '(''icoords_inter:'')'
                   print *, icoords_inter
                   print '(''->'')'
                   print *,test_icoords_inter
                end if

                if(.not.validated_rcoords_inter) then
                   print '(''rcoords_inter:'')'
                   print *, rcoords_inter
                   print '(''->'')'
                   print *,test_rcoords_inter
                end if
             end if

             deallocate(icoords_inter)
             deallocate(rcoords_inter)

             deallocate(test_icoords_inter)
             deallocate(test_rcoords_inter)

             deallocate(x_map_icr)
             deallocate(y_map_icr)

          end do

        end function test_get_inter_detector_rot_coords

        
        subroutine make_test_get_inter_detector_rot_coords(
     $     test_id,
     $     interior_x_map,
     $     interior_y_map,
     $     prev_icoords,
     $     next_icoords,
     $     test_icoords_inter,
     $     test_rcoords_inter)

          implicit none

          integer                                    , intent(in)  :: test_id
          real(rkind)   , dimension(:)               , intent(in)  :: interior_x_map
          real(rkind)   , dimension(:)               , intent(in)  :: interior_y_map
          integer(ikind), dimension(2)               , intent(out) :: prev_icoords
          integer(ikind), dimension(2)               , intent(out) :: next_icoords
          integer(ikind), dimension(:,:), allocatable, intent(out) :: test_icoords_inter
          real(rkind)   , dimension(:,:), allocatable, intent(out) :: test_rcoords_inter

          integer :: k

          select case(test_id)
            case(1)
               prev_icoords = [1,1]
               next_icoords = [2,3]
               
               allocate(test_icoords_inter(2,2))

               test_icoords_inter = reshape((/
     $              1,2,
     $              2,2
     $              /),
     $              (/2,2/)
     $              )

            case(2)
               prev_icoords = [1,3]
               next_icoords = [2,1]
               
               allocate(test_icoords_inter(2,2))

               test_icoords_inter = reshape((/
     $              1,2,
     $              2,2
     $              /),
     $              (/2,2/)
     $              )

            case(3)
               prev_icoords = [2,1]
               next_icoords = [1,3]
               
               allocate(test_icoords_inter(2,2))

               test_icoords_inter = reshape((/
     $              2,2,
     $              1,2
     $              /),
     $              (/2,2/)
     $              )

            case(4)
               prev_icoords = [2,1]
               next_icoords = [1,3]
               
               allocate(test_icoords_inter(2,2))

               test_icoords_inter = reshape((/
     $              2,2,
     $              1,2
     $              /),
     $              (/2,2/)
     $              )

            case(5)
               prev_icoords = [1,1]
               next_icoords = [3,3]
               
               allocate(test_icoords_inter(2,1))

               test_icoords_inter = reshape((/
     $              2,2
     $              /),
     $              (/2,1/)
     $              )

            case(6)
               prev_icoords = [1,3]
               next_icoords = [3,1]
               
               allocate(test_icoords_inter(2,1))

               test_icoords_inter = reshape((/
     $              2,2
     $              /),
     $              (/2,1/)
     $              )

            case(7)
               prev_icoords = [3,1]
               next_icoords = [1,3]
               
               allocate(test_icoords_inter(2,1))

               test_icoords_inter = reshape((/
     $              2,2
     $              /),
     $              (/2,1/)
     $              )

            case(8)
               prev_icoords = [3,3]
               next_icoords = [1,1]
               
               allocate(test_icoords_inter(2,1))

               test_icoords_inter = reshape((/
     $              2,2
     $              /),
     $              (/2,1/)
     $              )

             case(9)
               prev_icoords = [1,1]
               next_icoords = [4,3]
               
               allocate(test_icoords_inter(2,2))

               test_icoords_inter = reshape((/
     $              2,2,
     $              3,2
     $              /),
     $              (/2,2/)
     $              )               

            case(10)
               prev_icoords = [1,3]
               next_icoords = [4,1]
               
               allocate(test_icoords_inter(2,2))

               test_icoords_inter = reshape((/
     $              2,2,
     $              3,2
     $              /),
     $              (/2,2/)
     $              )

            case(11)
               prev_icoords = [4,1]
               next_icoords = [1,3]
               
               allocate(test_icoords_inter(2,2))

               test_icoords_inter = reshape((/
     $              3,2,
     $              2,2
     $              /),
     $              (/2,2/)
     $              )

            case(12)
               prev_icoords = [4,3]
               next_icoords = [1,1]
               
               allocate(test_icoords_inter(2,2))

               test_icoords_inter = reshape((/
     $              3,2,
     $              2,2
     $              /),
     $              (/2,2/)
     $              )

            case(13)
               prev_icoords = [1,1]
               next_icoords = [5,7]
               
               allocate(test_icoords_inter(2,5))

               test_icoords_inter = reshape((/
     $              2,2,
     $              2,3,
     $              3,4,
     $              4,5,
     $              4,6
     $              /),
     $              (/2,5/)
     $              )

            case(14)
               prev_icoords = [5,1]
               next_icoords = [1,7]
               
               allocate(test_icoords_inter(2,5))

               test_icoords_inter = reshape((/
     $              4,2,
     $              4,3,
     $              3,4,
     $              2,5,
     $              2,6
     $              /),
     $              (/2,5/)
     $              )

            case(15)
               prev_icoords = [1,6]
               next_icoords = [5,0]
               
               allocate(test_icoords_inter(2,5))

               test_icoords_inter = reshape((/
     $              2,5,
     $              2,4,
     $              3,3,
     $              4,2,
     $              4,1
     $              /),
     $              (/2,5/)
     $              )

            case(16)
               prev_icoords = [5,6]
               next_icoords = [1,0]
               
               allocate(test_icoords_inter(2,5))

               test_icoords_inter = reshape((/
     $              4,5,
     $              4,4,
     $              3,3,
     $              2,2,
     $              2,1
     $              /),
     $              (/2,5/)
     $              )

            case(17)
               prev_icoords = [1,1]
               next_icoords = [3,1]
               
               allocate(test_icoords_inter(2,1))

               test_icoords_inter = reshape((/
     $              2,1
     $              /),
     $              (/2,1/)
     $              )

            case(18)
               prev_icoords = [3,1]
               next_icoords = [1,1]
               
               allocate(test_icoords_inter(2,1))

               test_icoords_inter = reshape((/
     $              2,1
     $              /),
     $              (/2,1/)
     $              )

            case(19)
               prev_icoords = [1,1]
               next_icoords = [1,3]
               
               allocate(test_icoords_inter(2,1))

               test_icoords_inter = reshape((/
     $              1,2
     $              /),
     $              (/2,1/)
     $              )

            case(20)
               prev_icoords = [1,3]
               next_icoords = [1,1]
               
               allocate(test_icoords_inter(2,1))

               test_icoords_inter = reshape((/
     $              1,2
     $              /),
     $              (/2,1/)
     $              )

            case(21)
               prev_icoords = [1,1]
               next_icoords = [6,1]
               
               allocate(test_icoords_inter(2,4))

               test_icoords_inter = reshape((/
     $              2,1,
     $              3,1,
     $              4,1,
     $              5,1
     $              /),
     $              (/2,4/)
     $              )

            case(22)
               prev_icoords = [6,1]
               next_icoords = [1,1]
               
               allocate(test_icoords_inter(2,4))

               test_icoords_inter = reshape((/
     $              5,1,
     $              4,1,
     $              3,1,
     $              2,1
     $              /),
     $              (/2,4/)
     $              )

            case(23)
               prev_icoords = [1,1]
               next_icoords = [1,6]
               
               allocate(test_icoords_inter(2,4))

               test_icoords_inter = reshape((/
     $              1,2,
     $              1,3,
     $              1,4,
     $              1,5
     $              /),
     $              (/2,4/)
     $              )

            case(24)
               prev_icoords = [1,6]
               next_icoords = [1,1]
               
               allocate(test_icoords_inter(2,4))

               test_icoords_inter = reshape((/
     $              1,5,
     $              1,4,
     $              1,3,
     $              1,2
     $              /),
     $              (/2,4/)
     $              )

            case default
               print '(''make_test_get_inter_detector_rot_coords'')'
               print '(''test not implemented: '',I2)', test_id
          end select


          allocate(test_rcoords_inter(
     $         size(test_icoords_inter,1),
     $         size(test_icoords_inter,2)))

          do k=1, size(test_icoords_inter,2)

             test_rcoords_inter(1,k) = interior_x_map(test_icoords_inter(1,k))
             test_rcoords_inter(2,k) = interior_y_map(test_icoords_inter(2,k))

          end do

        end subroutine make_test_get_inter_detector_rot_coords


        function test_get_inter_detector_rot_param(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer(ikind), dimension(2)              :: prev_icoords
          integer(ikind), dimension(2)              :: next_icoords
          real(rkind)   , dimension(nx)             :: interior_x_map
          real(rkind)   , dimension(ny)             :: interior_y_map

          real(rkind)   , dimension(2)              :: rot_icoords_r
          real(rkind)   , dimension(2)              :: icoords_icr
          integer(ikind)                            :: inter_nb
          real(rkind)   , dimension(:), allocatable :: x_map_icr
          real(rkind)   , dimension(:), allocatable :: y_map_icr

          real(rkind)   , dimension(2)              :: test_rot_icoords_r
          real(rkind)   , dimension(2)              :: test_icoords_icr
          integer(ikind)                            :: test_inter_nb
          real(rkind)   , dimension(:), allocatable :: test_x_map_icr
          real(rkind)   , dimension(:), allocatable :: test_y_map_icr

          integer :: test_id

          logical :: validated_rot_icoords_r
          logical :: validated_icoords_icr
          logical :: validated_inter_nb
          logical :: validated_x_map_icr
          logical :: validated_y_map_icr

          test_validated = .true.


          interior_x_map = [0.1d0,0.2d0,0.3d0,0.4d0,0.45d0]
          interior_y_map = [0.0d0,0.2d0,0.4d0,0.6d0,0.8d0,0.9d0]

          do test_id=1,11
             
             call make_test_get_inter_detector_rot_param(
     $            test_id,
     $            prev_icoords,
     $            next_icoords,
     $            test_rot_icoords_r,
     $            test_icoords_icr,
     $            test_inter_nb,
     $            test_x_map_icr,
     $            test_y_map_icr)

             call get_inter_detector_rot_param(
     $            prev_icoords,
     $            next_icoords,
     $            interior_x_map,
     $            interior_y_map,
     $            rot_icoords_r,
     $            icoords_icr,
     $            inter_nb,
     $            x_map_icr,
     $            y_map_icr)

             validated_rot_icoords_r = is_vector_validated(
     $            rot_icoords_r,
     $            test_rot_icoords_r,
     $            detailled)
             
             validated_icoords_icr = is_vector_validated(
     $            icoords_icr,
     $            test_icoords_icr,
     $            detailled)
             
             validated_inter_nb = inter_nb.eq.test_inter_nb

             validated_x_map_icr = is_vector_validated(
     $            x_map_icr,
     $            test_x_map_icr,
     $            detailled)

             validated_y_map_icr = is_vector_validated(
     $            y_map_icr,
     $            test_y_map_icr,
     $            detailled)

             test_loc =
     $            validated_rot_icoords_r.and.
     $            validated_icoords_icr.and.
     $            validated_inter_nb.and.
     $            validated_x_map_icr.and.
     $            validated_y_map_icr

             test_validated = test_validated.and.test_loc

             if(detailled) then

                print '(''test '',I2)', test_id

                if(.not.validated_rot_icoords_r) then
                   print '(''rot_icoords_r: '',2F8.5,''->'',2F8.5)',
     $                  rot_icoords_r,
     $                  test_rot_icoords_r
                end if
                
                if(.not.validated_icoords_icr) then
                   print '(''icoords_icr: '',2F8.5,''->'',2F8.5)',
     $               icoords_icr,
     $                  test_icoords_icr
                end if
                
                if(.not.validated_inter_nb) then
                   print '(''inter_nb: '',I2,''->'',I2)',
     $                  inter_nb,
     $                  test_inter_nb
                end if
                
                if(.not.validated_x_map_icr) then
                   print '(''x_map_icr'')'
                   print *, x_map_icr
                   print '(''->'')'
                   print *, test_x_map_icr
                end if
                
                if(.not.validated_y_map_icr) then
                   print '(''y_map_icr'')'
                   print *, y_map_icr
                   print '(''->'')'
                   print *, test_y_map_icr
                end if
             end if

             deallocate(x_map_icr)
             deallocate(y_map_icr)

          end do

        end function test_get_inter_detector_rot_param


        subroutine make_test_get_inter_detector_rot_param(
     $     test_id,
     $     prev_icoords,
     $     next_icoords,
     $     rot_icoords_r,
     $     icoords_icr,
     $     inter_nb,
     $     x_map_icr,
     $     y_map_icr)

          implicit none

          integer                                  , intent(in)  :: test_id
          integer(ikind), dimension(2)             , intent(out) :: prev_icoords
          integer(ikind), dimension(2)             , intent(out) :: next_icoords
          real(rkind)   , dimension(2)             , intent(out) :: rot_icoords_r
          real(rkind)   , dimension(2)             , intent(out) :: icoords_icr
          integer(ikind)                           , intent(out) :: inter_nb
          real(rkind)   , dimension(:), allocatable, intent(out) :: x_map_icr
          real(rkind)   , dimension(:), allocatable, intent(out) :: y_map_icr


          select case(test_id)
            case(1)
               prev_icoords = [1,1]
               next_icoords = [2,3]
               
               rot_icoords_r = [1.5d0,2.0d0]
               
               icoords_icr = [1.0d0/3.0d0,2.0d0/3.0d0]

               inter_nb = 2

               allocate(x_map_icr(2))
               allocate(y_map_icr(2))
               
               x_map_icr = [0.1d0,0.2d0]
               y_map_icr = [0.0d0,0.2d0]

            case(2)
               prev_icoords = [1,3]
               next_icoords = [2,1]
               
               rot_icoords_r = [1.5d0,2.0d0]
               
               icoords_icr = [1.0d0/3.0d0,-2.0d0/3.0d0]

               inter_nb = 2

               allocate(x_map_icr(2))
               allocate(y_map_icr(2))
               
               x_map_icr = [0.1d0,0.2d0]
               y_map_icr = [0.4d0,0.2d0]

            case(3)
               prev_icoords = [2,1]
               next_icoords = [1,3]
               
               rot_icoords_r = [1.5d0,2.0d0]
               
               icoords_icr = [-1.0d0/3.0d0,2.0d0/3.0d0]

               inter_nb = 2

               allocate(x_map_icr(2))
               allocate(y_map_icr(2))
               
               x_map_icr = [0.2d0,0.1d0]
               y_map_icr = [0.0d0,0.2d0]

            case(4)
               prev_icoords = [2,3]
               next_icoords = [1,1]
               
               rot_icoords_r = [1.5d0,2.0d0]
               
               icoords_icr = [-1.0d0/3.0d0,-2.0d0/3.0d0]

               inter_nb = 2

               allocate(x_map_icr(2))
               allocate(y_map_icr(2))
               
               x_map_icr = [0.2d0,0.1d0]
               y_map_icr = [0.4d0,0.2d0]

            case(5)
               prev_icoords = [1,1]
               next_icoords = [3,3]
               
               rot_icoords_r = [2.0d0,2.0d0]
               
               icoords_icr = [1.0d0,1.0d0]

               inter_nb = 1

               allocate(x_map_icr(2))
               allocate(y_map_icr(2))
               
               x_map_icr = [0.1d0,0.2d0]
               y_map_icr = [0.0d0,0.2d0]

            case(6)
               prev_icoords = [1,3]
               next_icoords = [3,1]
               
               rot_icoords_r = [2.0d0,2.0d0]
               
               icoords_icr = [1.0d0,-1.0d0]

               inter_nb = 1

               allocate(x_map_icr(2))
               allocate(y_map_icr(2))
               
               x_map_icr = [0.1d0,0.2d0]
               y_map_icr = [0.4d0,0.2d0]

            case(7)
               prev_icoords = [3,1]
               next_icoords = [1,3]
               
               rot_icoords_r = [2.0d0,2.0d0]
               
               icoords_icr = [-1.0d0,1.0d0]

               inter_nb = 1

               allocate(x_map_icr(2))
               allocate(y_map_icr(2))
               
               x_map_icr = [0.3d0,0.2d0]
               y_map_icr = [0.0d0,0.2d0]

            case(8)
               prev_icoords = [3,3]
               next_icoords = [1,1]
               
               rot_icoords_r = [2.0d0,2.0d0]
               
               icoords_icr = [-1.0d0,-1.0d0]

               inter_nb = 1

               allocate(x_map_icr(2))
               allocate(y_map_icr(2))
               
               x_map_icr = [0.3d0,0.2d0]
               y_map_icr = [0.4d0,0.2d0]

            case(9)
               prev_icoords = [1,1]
               next_icoords = [4,3]
               
               rot_icoords_r = [2.5d0,2.0d0]
               
               icoords_icr = [1.0d0,2.0d0/3.0d0]

               inter_nb = 2

               allocate(x_map_icr(3))
               allocate(y_map_icr(2))
               
               x_map_icr = [0.1d0,0.2d0,0.3d0]
               y_map_icr = [0.0d0,0.2d0]

            case(10)
               prev_icoords = [1,1]
               next_icoords = [6,5]
               
               rot_icoords_r = [3.5d0,3.0d0]
               
               icoords_icr = [1.0d0,4.0d0/5.0d0]

               inter_nb = 4

               allocate(x_map_icr(5))
               allocate(y_map_icr(4))
               
               x_map_icr = [0.1d0,0.2d0,0.3d0,0.4d0,0.45d0]
               y_map_icr = [0.0d0,0.2d0,0.4d0,0.6d0]

            case(11)
               prev_icoords = [1,1]
               next_icoords = [5,7]
               
               rot_icoords_r = [3.0d0,4.0d0]
               
               icoords_icr = [2.0d0/3.0d0,1.0d0]

               inter_nb = 5

               allocate(x_map_icr(4))
               allocate(y_map_icr(6))
               
               x_map_icr = [0.1d0,0.2d0,0.3d0,0.4d0]
               y_map_icr = [0.0d0,0.2d0,0.4d0,0.6d0,0.8d0,0.9d0]

            case default
               print '(''make_test_get_inter_detector_rot_param'')'
               print '(''test not implemented: '',I2)', test_id
          end select

        end subroutine make_test_get_inter_detector_rot_param


        function test_get_minimum_block_nb(detailled)
     $       result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          integer(ikind), dimension(2) :: prev_icoords
          integer(ikind), dimension(2) :: next_icoords
          integer(ikind)               :: block_nb
          integer(ikind)               :: test_block_nb

          integer :: test_id
          logical :: test_loc
          

          test_validated = .true.

          do test_id=1,20
              
             call make_test_get_minimum_block_nb(
     $             test_id,
     $             prev_icoords,
     $             next_icoords,
     $             test_block_nb)

             block_nb = get_minimum_block_nb(
     $            prev_icoords,
     $            next_icoords)

             test_loc = block_nb.eq.test_block_nb
             test_validated = test_validated.and.test_loc
             

             if((.not.test_loc).and.detailled) then
                print '(''test '',I2,'': '',I2,''->'',I2)',
     $               test_id, block_nb, test_block_nb
             end if

          end do

        end function test_get_minimum_block_nb


        subroutine make_test_get_minimum_block_nb(
     $     test_id,
     $     prev_icoords,
     $     next_icoords,
     $     test_block_nb)

          implicit none

          integer                     , intent(in)  :: test_id
          integer(ikind), dimension(2), intent(out) :: prev_icoords
          integer(ikind), dimension(2), intent(out) :: next_icoords
          integer(ikind)              , intent(out) :: test_block_nb


          select case(test_id)
            case(1)
               prev_icoords  = [0,0]
               next_icoords  = [-1,-1]
               test_block_nb = 0

            case(2)
               prev_icoords  = [0,0]
               next_icoords  = [0,-1]
               test_block_nb = 0

            case(3)
               prev_icoords  = [0,0]
               next_icoords  = [1,-1]
               test_block_nb = 0

            case(4)
               prev_icoords  = [0,0]
               next_icoords  = [-1,0]
               test_block_nb = 0

            case(5)
               prev_icoords  = [0,0]
               next_icoords  = [1,0]
               test_block_nb = 0

            case(6)
               prev_icoords  = [0,0]
               next_icoords  = [-1,1]
               test_block_nb = 0

            case(7)
               prev_icoords  = [0,0]
               next_icoords  = [0,1]
               test_block_nb = 0

            case(8)
               prev_icoords  = [0,0]
               next_icoords  = [1,1]
               test_block_nb = 0

            case(9)
               prev_icoords  = [0,0]
               next_icoords  = [-2,-2]
               test_block_nb = 1

            case(10)
               prev_icoords  = [0,0]
               next_icoords  = [2,-2]
               test_block_nb = 1

            case(11)
               prev_icoords  = [0,0]
               next_icoords  = [-2,2]
               test_block_nb = 1

            case(12)
               prev_icoords  = [0,0]
               next_icoords  = [2,2]
               test_block_nb = 1

            case(13)
               prev_icoords  = [0,0]
               next_icoords  = [0,-2]
               test_block_nb = 1

            case(14)
               prev_icoords  = [0,0]
               next_icoords  = [-2,0]
               test_block_nb = 1

            case(15)
               prev_icoords  = [0,0]
               next_icoords  = [2,0]
               test_block_nb = 1

            case(16)
               prev_icoords  = [0,0]
               next_icoords  = [0,2]
               test_block_nb = 1

            case(17)
               prev_icoords  = [0,0]
               next_icoords  = [-3,-1]
               test_block_nb = 2

            case(18)
               prev_icoords  = [0,0]
               next_icoords  = [3,-1]
               test_block_nb = 2

            case(19)
               prev_icoords  = [0,0]
               next_icoords  = [-3,1]
               test_block_nb = 2

            case(20)
               prev_icoords  = [0,0]
               next_icoords  = [3,1]
               test_block_nb = 2

            case default
               print '(''make_test_get_minimum_block_nb'')'
               print '(''test not implemented: '',I2)', test_id

          end select

        end subroutine make_test_get_minimum_block_nb

        
        function test_get_rot_coords(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer(ikind), dimension(2) :: prev_icoords
          integer(ikind), dimension(2) :: next_icoords
          integer(ikind), dimension(2) :: rot_icoords
          real(rkind)   , dimension(2) :: rot_rcoords
          integer(ikind), dimension(2) :: test_rot_icoords
          real(rkind)   , dimension(2) :: test_rot_rcoords

          integer :: test_id
          logical :: test_loc

          logical :: validated_rot_icoords
          logical :: validated_rot_rcoords

          test_validated = .true.

          do test_id=1,16

             call make_test_get_rot_coords(
     $            test_id,
     $            prev_icoords,
     $            next_icoords,
     $            test_rot_icoords,
     $            test_rot_rcoords)

             call get_rot_coords(
     $            prev_icoords,
     $            next_icoords,
     $            rot_icoords,
     $            rot_rcoords)

             validated_rot_icoords =
     $            (test_rot_icoords(1).eq.rot_icoords(1)).and.
     $            (test_rot_icoords(2).eq.rot_icoords(2))
             
             validated_rot_rcoords = is_vector_validated(
     $            rot_rcoords,
     $            test_rot_rcoords,
     $            detailled)

             test_loc = validated_rot_icoords.and.
     $            validated_rot_rcoords
          
             test_validated = test_validated.and.
     $            test_loc

             if(.not.validated_rot_icoords) then
                print '(''rot_icoords: '',2I2,''->'',2I2)',
     $               rot_icoords,
     $               test_rot_icoords
             end if

             if(.not.validated_rot_rcoords) then
                print '(''rot_icoords: '',2F8.5,''->'',2F8.5)',
     $               rot_icoords,
     $               test_rot_icoords                
             end if
                
          end do

        end function test_get_rot_coords


        subroutine make_test_get_rot_coords(
     $     test_id,
     $     prev_icoord,
     $     next_icoord,
     $     test_rot_icoords,
     $     test_rot_rcoords)

          implicit none

          integer                     , intent(in)  :: test_id
          integer(ikind), dimension(2), intent(out) :: prev_icoord
          integer(ikind), dimension(2), intent(out) :: next_icoord
          integer(ikind), dimension(2), intent(out) :: test_rot_icoords
          real(rkind)   , dimension(2), intent(out) :: test_rot_rcoords

          select case(test_id)
            case(1)
               prev_icoord(1) = 1
               prev_icoord(2) = 1
               
               next_icoord(1) = 5
               next_icoord(2) = 3

               test_rot_icoords(1) = 3
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 3.0d0
               test_rot_rcoords(2) = 2.0d0

            case(2)
               prev_icoord(1) = 5
               prev_icoord(2) = 1
               
               next_icoord(1) = 1
               next_icoord(2) = 3

               test_rot_icoords(1) = 3
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 3.0d0
               test_rot_rcoords(2) = 2.0d0

            case(3)
               prev_icoord(1) = 5
               prev_icoord(2) = 3
               
               next_icoord(1) = 1
               next_icoord(2) = 1

               test_rot_icoords(1) = 3
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 3.0d0
               test_rot_rcoords(2) = 2.0d0

            case(4)
               prev_icoord(1) = 1
               prev_icoord(2) = 3
               
               next_icoord(1) = 5
               next_icoord(2) = 1

               test_rot_icoords(1) = 3
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 3.0d0
               test_rot_rcoords(2) = 2.0d0

            case(5)
               prev_icoord(1) = 1
               prev_icoord(2) = 1
               
               next_icoord(1) = 4
               next_icoord(2) = 3

               test_rot_icoords(1) = 2
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 2.5d0
               test_rot_rcoords(2) = 2.0d0

            case(6)
               prev_icoord(1) = 4
               prev_icoord(2) = 1
               
               next_icoord(1) = 1
               next_icoord(2) = 3

               test_rot_icoords(1) = 3
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 2.5d0
               test_rot_rcoords(2) = 2.0d0

            case(7)
               prev_icoord(1) = 4
               prev_icoord(2) = 3
               
               next_icoord(1) = 1
               next_icoord(2) = 1

               test_rot_icoords(1) = 3
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 2.5d0
               test_rot_rcoords(2) = 2.0d0
              
            case(8)
               prev_icoord(1) = 1
               prev_icoord(2) = 3
               
               next_icoord(1) = 4
               next_icoord(2) = 1

               test_rot_icoords(1) = 2
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 2.5d0
               test_rot_rcoords(2) = 2.0d0

            case(9)
               prev_icoord(1) = 1
               prev_icoord(2) = 1
               
               next_icoord(1) = 3
               next_icoord(2) = 4

               test_rot_icoords(1) = 2
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 2.0d0
               test_rot_rcoords(2) = 2.5d0

            case(10)
               prev_icoord(1) = 3
               prev_icoord(2) = 1
               
               next_icoord(1) = 1
               next_icoord(2) = 4

               test_rot_icoords(1) = 2
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 2.0d0
               test_rot_rcoords(2) = 2.5d0

            case(11)
               prev_icoord(1) = 3
               prev_icoord(2) = 4
               
               next_icoord(1) = 1
               next_icoord(2) = 1

               test_rot_icoords(1) = 2
               test_rot_icoords(2) = 3

               test_rot_rcoords(1) = 2.0d0
               test_rot_rcoords(2) = 2.5d0

            case(12)
               prev_icoord(1) = 1
               prev_icoord(2) = 4
               
               next_icoord(1) = 3
               next_icoord(2) = 1

               test_rot_icoords(1) = 2
               test_rot_icoords(2) = 3

               test_rot_rcoords(1) = 2.0d0
               test_rot_rcoords(2) = 2.5d0

            case(13)
               prev_icoord(1) = 1
               prev_icoord(2) = 1
               
               next_icoord(1) = 4
               next_icoord(2) = 4

               test_rot_icoords(1) = 2
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 2.5d0
               test_rot_rcoords(2) = 2.5d0

            case(14)
               prev_icoord(1) = 4
               prev_icoord(2) = 1
               
               next_icoord(1) = 1
               next_icoord(2) = 4

               test_rot_icoords(1) = 3
               test_rot_icoords(2) = 2

               test_rot_rcoords(1) = 2.5d0
               test_rot_rcoords(2) = 2.5d0

            case(15)
               prev_icoord(1) = 4
               prev_icoord(2) = 4
               
               next_icoord(1) = 1
               next_icoord(2) = 1

               test_rot_icoords(1) = 3
               test_rot_icoords(2) = 3

               test_rot_rcoords(1) = 2.5d0
               test_rot_rcoords(2) = 2.5d0

            case(16)
               prev_icoord(1) = 1
               prev_icoord(2) = 4
               
               next_icoord(1) = 4
               next_icoord(2) = 1

               test_rot_icoords(1) = 2
               test_rot_icoords(2) = 3

               test_rot_rcoords(1) = 2.5d0
               test_rot_rcoords(2) = 2.5d0

            case default
               print '(''make_test_get_rot_coords'')'
               print '(''test_not_implemented: '', I2)', test_id
               print '()'

          end select

        end subroutine make_test_get_rot_coords


        function test_get_inter_detector_coords(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer(ikind), dimension(2) :: prev_icoord
          real(rkind)   , dimension(2) :: icoord_icr
          integer                      :: k
          real(rkind)   , dimension(3) :: x_map_icr
          real(rkind)   , dimension(2) :: y_map_icr
          
          integer(ikind), dimension(2) :: icoord_inter
          real(rkind)   , dimension(2) :: rcoord_inter
          integer(ikind), dimension(2) :: test_icoord_inter
          real(rkind)   , dimension(2) :: test_rcoord_inter
          
          integer :: test_id

          logical :: validated_icoord_inter
          logical :: validated_rcoord_inter

          do test_id=1,2

             !test presets
             call make_test_get_inter_detector_coords(
     $            test_id,
     $            prev_icoord,
     $            icoord_icr,
     $            k,
     $            x_map_icr,
     $            y_map_icr,
     $            test_icoord_inter,
     $            test_rcoord_inter)

             !compute the output with the tested function
             call get_inter_detector_coords(
     $            prev_icoord,
     $            icoord_icr,
     $            k,
     $            x_map_icr,
     $            y_map_icr,
     $            icoord_inter,
     $            rcoord_inter)
             
             !compare results
             validated_icoord_inter = 
     $            (icoord_inter(1).eq.test_icoord_inter(1)).and.
     $            (icoord_inter(2).eq.test_icoord_inter(2))

             validated_rcoord_inter = is_vector_validated(
     $            rcoord_inter,
     $            test_rcoord_inter,
     $            detailled)

             test_validated = validated_icoord_inter.and.
     $                        validated_rcoord_inter

             !display comparison results
             if(detailled) then

                print '(''test '',I2)', test_id

                print '(''  icoord_inter: '',L1)', validated_icoord_inter
                if(.not.validated_icoord_inter) then
                   print '(''     - ['',2I2,''] -> ['',2I2,'']'')',
     $                  icoord_inter, test_icoord_inter
                   print '()'
                end if

                print '(''  rcoord_inter: '',L1)', validated_rcoord_inter
                if(.not.validated_rcoord_inter) then
                   print '(''     - ['',2F8.5,''] -> ['',2F8.5,'']'')',
     $                  rcoord_inter, test_rcoord_inter
                   print '()'
                end if

             end if

          end do

        end function test_get_inter_detector_coords

      
        subroutine make_test_get_inter_detector_coords(
     $     test_id,
     $     prev_icoord,
     $     icoord_icr,
     $     k,
     $     x_map_icr,
     $     y_map_icr,
     $     test_icoord_inter,
     $     test_rcoord_inter)

          implicit none

          integer                     , intent(in)  :: test_id
          integer(ikind), dimension(2), intent(out) :: prev_icoord
          real(rkind)   , dimension(2), intent(out) :: icoord_icr
          integer(ikind)              , intent(out) :: k
          real(rkind)   , dimension(3), intent(out) :: x_map_icr
          real(rkind)   , dimension(2), intent(out) :: y_map_icr
          integer(ikind), dimension(2), intent(out) :: test_icoord_inter
          real(rkind)   , dimension(2), intent(out) :: test_rcoord_inter


          select case(test_id)
            case(1)
               prev_icoord(1) = 1
               prev_icoord(2) = 1
               
               icoord_icr(1) = 1.0d0
               icoord_icr(2) = 0.5d0
               
               k=2

               x_map_icr = [0.0d0,0.1d0,0.2d0]
               y_map_icr = [0.1d0,0.3d0]

               test_icoord_inter = [3,2]
               test_rcoord_inter = [0.2d0,0.3d0]

            case(2)
               prev_icoord(1) = 4
               prev_icoord(2) = 3
               
               icoord_icr(1) =-1.0d0
               icoord_icr(2) =-0.5d0
               
               k=2

               x_map_icr = [0.3d0,0.2d0,0.1d0]
               y_map_icr = [0.6d0,0.4d0]

               test_icoord_inter = [2,2]
               test_rcoord_inter = [0.1d0,0.4d0]

            case default
               print '(''make_test_get_inter_detector_coords'')'
               print '(''test not implemeted: '',I2)', test_id
               stop ''
          end select        

        end subroutine make_test_get_inter_detector_coords


        
c$$$        function test_get_inter_detector_param_sym(detailled)
c$$$     $       result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$          logical             :: test_validated
c$$$
c$$$          integer(ikind), dimension(2)              :: prev_icoord
c$$$          integer(ikind), dimension(2)              :: next_icoord
c$$$          real(rkind)   , dimension(nx)             :: interior_x_map
c$$$          real(rkind)   , dimension(ny)             :: interior_y_map
c$$$
c$$$          real(rkind)   , dimension(2)              :: icoord_icr
c$$$          integer                                   :: inter_nb
c$$$          real(rkind)   , dimension(:), allocatable :: x_map_icr
c$$$          real(rkind)   , dimension(:), allocatable :: y_map_icr
c$$$
c$$$          integer(ikind), dimension(2)              :: icoord_inter
c$$$          real(rkind)   , dimension(2)              :: rcoord_inter
c$$$
c$$$          integer :: k
c$$$
c$$$          interior_x_map = [0.1d0,0.2d0,0.3d0,0.4d0,0.45d0]
c$$$          interior_y_map = [0.0d0,0.2d0,0.4d0,0.6d0,0.8d0,0.9d0]
c$$$          
c$$$          !  ___ ___
c$$$          ! |XXX|   |
c$$$          ! |XXX|___|
c$$$          ! |   |   |
c$$$          ! |___|___|
c$$$          ! |   |XXX|
c$$$          ! |___|XXX|
c$$$          !------------------------------
c$$$          prev_icoord=[1,3]
c$$$          next_icoord=[2,1]
c$$$
c$$$          call get_inter_detector_param(
c$$$     $         prev_icoord,
c$$$     $         next_icoord,
c$$$     $         interior_x_map,
c$$$     $         interior_y_map,
c$$$     $         icoord_icr,
c$$$     $         inter_nb,
c$$$     $         x_map_icr,
c$$$     $         y_map_icr)
c$$$
c$$$          print '(''prev_icoord: '',2I2)', prev_icoord
c$$$          print '(''next_icoord: '',2I2)', next_icoord
c$$$          print '(''--------------------'')'
c$$$          print '(''icoord_icr     :'',2F9.5)', icoord_icr
c$$$          print '(''inter_icr      : '',I2)', inter_nb
c$$$
c$$$          do k=1, inter_nb
c$$$             
c$$$             call get_inter_detector_coords(
c$$$     $            prev_icoord,
c$$$     $            icoord_icr,
c$$$     $            k,
c$$$     $            x_map_icr,
c$$$     $            y_map_icr,
c$$$     $            icoord_inter,
c$$$     $            rcoord_inter)
c$$$
c$$$             print '(''icoord_inter('',I1,''): '',2I2)', k,icoord_inter
c$$$
c$$$          end do
c$$$
c$$$          deallocate(x_map_icr)
c$$$          deallocate(y_map_icr)
c$$$
c$$$          print '()'
c$$$
c$$$
c$$$          !  ___ ___
c$$$          ! |   |XXX|
c$$$          ! |___|XXX|
c$$$          ! |   |   |
c$$$          ! |___|___|
c$$$          ! |XXX|   |
c$$$          ! |XXX|___|
c$$$          !------------------------------
c$$$          call get_inter_detector_param(
c$$$     $         prev_icoord,
c$$$     $         next_icoord,
c$$$     $         interior_x_map,
c$$$     $         interior_y_map,
c$$$     $         icoord_icr,
c$$$     $         inter_nb,
c$$$     $         x_map_icr,
c$$$     $         y_map_icr)
c$$$
c$$$          print '(''prev_icoord: '',2I2)', prev_icoord
c$$$          print '(''next_icoord: '',2I2)', next_icoord
c$$$          print '(''--------------------'')'
c$$$          print '(''icoord_icr     :'',2F9.5)', icoord_icr
c$$$          print '(''inter_icr      : '',I2)', inter_nb
c$$$
c$$$          do k=1, inter_nb
c$$$             
c$$$             call get_inter_detector_coords(
c$$$     $            prev_icoord,
c$$$     $            icoord_icr,
c$$$     $            k,
c$$$     $            x_map_icr,
c$$$     $            y_map_icr,
c$$$     $            icoord_inter,
c$$$     $            rcoord_inter)
c$$$
c$$$             print '(''icoord_inter('',I1,''): '',2I2)', k,icoord_inter
c$$$
c$$$          end do
c$$$
c$$$          deallocate(x_map_icr)
c$$$          deallocate(y_map_icr)
c$$$          
c$$$          print '()'
c$$$
c$$$        end function test_get_inter_detector_param_sym


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
          
          do test_id=1,6

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

                print '(''test '',I2)', test_id

                print '(''  icoord_icr: '',L1)', validated_icoord_icr
                if(.not.validated_icoord_icr) then
                   print '(''     - icoord_icr(1): '',F8.5,''->'',F8.5)',
     $                  icoord_icr(1), test_icoord_icr(1)
                   print '(''     - icoord_icr(2): '',F8.5,''->'',F8.5)',
     $                  icoord_icr(2), test_icoord_icr(2)
                   print '()'
                end if

                print '(''  inter_nb: '',L1)', validated_inter_nb
                if(.not.validated_inter_nb) then
                   print '(''     - '',I2,'' -> '',I2)',
     $                  inter_nb, test_inter_nb
                   print '()'
                end if

                print '(''  x_map_icr: '',L1)', validated_x_map_icr
                if(.not.validated_x_map_icr) then
                   print '(''   -  x_map_icr: '')'
                   print *, x_map_icr
                   print '(''      -> '')'
                   print '(''   -  test_x_map_icr: '')'
                   print *, test_x_map_icr
                   print '()'
                end if

                print '(''  y_map_icr: '',L1)', validated_y_map_icr
                if(.not.validated_y_map_icr) then
                   print '(''   -  y_map_icr: '')'
                   print *, y_map_icr
                   print '(''      -> '')'
                   print '(''   -  test_y_map_icr: '')'
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
            !i_icr>0, j_icr=0
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

            !i_icr<0, j_icr=0
            case(2)
               prev_icoord(1) = 4
               prev_icoord(2) = 1

               next_icoord(1) = 1
               next_icoord(2) = 1

               test_icoord_icr(1) =-1.0d0
               test_icoord_icr(2) = 0.0d0

               test_inter_nb = 2

               allocate(test_x_map_icr(3))
               allocate(test_y_map_icr(1))

               test_x_map_icr = [0.4d0,0.3d0,0.2d0]
               test_y_map_icr = [0.0d0]

            !i_icr=0, j_icr>0
            case(3)
               prev_icoord(1) = 1
               prev_icoord(2) = 1

               next_icoord(1) = 1
               next_icoord(2) = 4

               test_icoord_icr(1) = 0.0d0
               test_icoord_icr(2) = 1.0d0

               test_inter_nb = 2

               allocate(test_x_map_icr(1))
               allocate(test_y_map_icr(3))

               test_x_map_icr = [0.1d0]
               test_y_map_icr = [0.0d0,0.2d0,0.4d0]

            !i_icr=0, j_icr<0
            case(4)
               prev_icoord(1) = 1
               prev_icoord(2) = 4

               next_icoord(1) = 1
               next_icoord(2) = 1

               test_icoord_icr(1) = 0.0d0
               test_icoord_icr(2) =-1.0d0

               test_inter_nb = 2

               allocate(test_x_map_icr(1))
               allocate(test_y_map_icr(3))

               test_x_map_icr = [0.1d0]
               test_y_map_icr = [0.6d0,0.4d0,0.2d0]

            !i_icr>0, j_icr>0
            case(5)
               prev_icoord(1) = 1
               prev_icoord(2) = 1

               next_icoord(1) = 4
               next_icoord(2) = 3

               test_icoord_icr(1) = 1.0d0
               test_icoord_icr(2) = 0.5d0

               test_inter_nb = 2

               allocate(test_x_map_icr(3))
               allocate(test_y_map_icr(2))

               test_x_map_icr = [0.1d0,0.2d0,0.3d0]
               test_y_map_icr = [0.0d0,0.2d0]

            !i_icr<0, j_icr>0
            case(6)
               prev_icoord(1) = 4
               prev_icoord(2) = 1

               next_icoord(1) = 1
               next_icoord(2) = 3

               test_icoord_icr(1) =-1.0d0
               test_icoord_icr(2) = 0.5d0

               test_inter_nb = 2

               allocate(test_x_map_icr(3))
               allocate(test_y_map_icr(2))

               test_x_map_icr = [0.4d0,0.3d0,0.2d0]
               test_y_map_icr = [0.0d0,0.2d0]

            !i_icr>0, j_icr<0
            case(7)
               prev_icoord(1) = 1
               prev_icoord(2) = 3

               next_icoord(1) = 4
               next_icoord(2) = 1

               test_icoord_icr(1) = 1.0d0
               test_icoord_icr(2) =-0.5d0

               test_inter_nb = 2

               allocate(test_x_map_icr(3))
               allocate(test_y_map_icr(2))

               test_x_map_icr = [0.1d0,0.2d0,0.3d0]
               test_y_map_icr = [0.4d0,0.2d0]

            !i_icr<0, j_icr<0
            case(8)
               prev_icoord(1) = 4
               prev_icoord(2) = 3

               next_icoord(1) = 1
               next_icoord(2) = 1

               test_icoord_icr(1) =-1.0d0
               test_icoord_icr(2) =-0.5d0

               test_inter_nb = 2

               allocate(test_x_map_icr(3))
               allocate(test_y_map_icr(2))

               test_x_map_icr = [0.4d0,0.3d0,0.2d0]
               test_y_map_icr = [0.4d0,0.2d0]

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
          logical                      :: sign

          integer :: test_id

          test_validated = .true.

          size_interior_map = size(interior_map,1)
          interior_map = [0.1d0,0.2d0,0.3d0,0.4d0,0.45d0]

          size_local_map = 3

          do test_id=1,12

             !get the test presets
             call make_test_determine_local_map_coordinates(
     $            test_id, first_icoord, sign, test_local_map)

             !compute the local map using the tested function
             call determine_local_map_coordinates(
     $            interior_map,
     $            size_interior_map,
     $            size_local_map,
     $            first_icoord,
     $            sign,
     $            local_map)

             !compare the results
             test_loc = is_vector_validated(
     $            local_map,
     $            test_local_map,
     $            detailled)
             test_validated = test_validated.and.test_loc

             !print the potential differences if detailled
             if(detailled) then

                if(.not.test_loc) then
                   print '(''**test '',I2,'' failed**'')', test_id
                   print '(''local_map     : '',3F8.3)', local_map
                   print '(''test_local_map: '',3F8.3)', test_local_map
                   print ''                
                else
                   print '(''test '',I2,'' validated'')', test_id
                end if

             end if

          end do          

        end function test_determine_local_map_coordinates        

        
        subroutine make_test_determine_local_map_coordinates(
     $     test_id, first_icoord, sign, test_local_map)

          implicit none

          integer                  , intent(in)  :: test_id
          integer                  , intent(out) :: first_icoord
          logical                  , intent(out) :: sign
          real(rkind), dimension(3), intent(out) :: test_local_map

          select case(test_id)

            case(1)
               first_icoord   = -2
               sign           = .true.
               test_local_map = [-0.2d0,-0.1d0,0.0d0]

            case(2)
               first_icoord   = 0
               sign           = .true.
               test_local_map = [0.0d0,0.1d0,0.2d0]

            case(3)
               first_icoord   = 1
               sign           = .true.
               test_local_map = [0.1d0,0.2d0,0.3d0]

            case(4)
               first_icoord   = 3
               sign           = .true.
               test_local_map = [0.3d0,0.4d0,0.45d0]

            case(5)
               first_icoord   = 5
               sign           = .true.
               test_local_map = [0.45d0,0.5d0,0.55d0]

            case(6)
               first_icoord   = 6
               sign           = .true.
               test_local_map = [0.5d0,0.55d0,0.6d0]

            case(7)
               first_icoord   = 8
               sign           = .false.
               test_local_map = [0.6d0,0.55d0,0.5d0]

            case(8)
               first_icoord   = 6
               sign           = .false.
               test_local_map = [0.5d0,0.45d0,0.4d0]
               
            case(9)
               first_icoord   = 5
               sign           = .false.
               test_local_map = [0.45d0,0.4d0,0.3d0]

            case(10)
               first_icoord   = 3
               sign           = .false.
               test_local_map = [0.3d0,0.2d0,0.1d0]

            case(11)
               first_icoord   = 1
               sign           = .false.
               test_local_map = [0.1d0,0.0d0,-0.1d0]

            case(12)
               first_icoord   = 0
               sign           = .false.
               test_local_map = [0.0d0,-0.1d0,-0.2d0]

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

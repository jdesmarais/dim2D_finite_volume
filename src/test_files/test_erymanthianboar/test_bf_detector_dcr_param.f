      program test_bf_detector_dcr_param

        use bf_detector_dcr_param_class, only :
     $     bf_detector_dcr_param,
     $     
     $     compute_new_list_param_g,
     $     finalize_new_list_g,
     $     pinpoint_detector_for_removal,
     $     get_border_detector_N,
     $     get_border_detector_S,
     $     get_border_detector_E,
     $     get_border_detector_W,
     $     get_detector_changes,
     $     get_segment_first_pt,
     $     replace_removed_detectors,
     $     should_be_removed_N,
     $     should_be_removed_S,
     $     should_be_removed_E,
     $     should_be_removed_W,
     $     get_rcoord
        
        use parameters_bf_layer, only :
     $       dct_icr_distance,
     $       dct_icr_N_default,
     $       dct_icr_S_default,
     $       dct_icr_E_default,
     $       dct_icr_W_default

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny

        use parameters_kind, only :
     $       ikind, rkind


        implicit none

        logical                     :: detailled
        logical                     :: test_validated


        if((nx.ne.20).or.(ny.ne.20).or.(dct_icr_distance.ne.4)) then
           print '(''the test requires nx=ny=20'')'
           print '(''the test requires dct_icr_distance=4'')'
           stop ''
        end if

        detailled = .false.
        test_validated = test_get_border_detector(detailled)
        print '(''test_get_border_detector: '',L1)', test_validated
        print '()'

        detailled = .false.
        test_validated = test_should_be_removed(detailled)
        print '(''test_should_be_removed: '',L1)', test_validated
        print '()'

        detailled = .false.
        test_validated = test_get_rcoords(detailled)
        print '(''test_get_rcoords: '',L1)', test_validated
        print '()'

        detailled = .false.
        test_validated = test_compute_new_list_param(detailled)
        print '(''test_compute_new_list_param: '',L1)', test_validated
        print '()'

        detailled = .true.
        test_validated = test_finalize_new_list(detailled)
        print '(''test_finalize_new_list: '',L1)', test_validated
        print '()'

        contains

        function test_get_border_detector(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          logical :: test_loc

          real(rkind) :: dx
          real(rkind) :: dy

          real(rkind)   , dimension(nx) :: interior_x_map
          real(rkind)   , dimension(ny) :: interior_y_map

          integer(ikind), dimension(2)  :: icoord
          real(rkind)   , dimension(2)  :: rcoord
          integer(ikind), dimension(2)  :: icoord_n
          real(rkind)   , dimension(2)  :: rcoord_n
          integer(ikind), dimension(2)  :: icoord_test
          real(rkind)   , dimension(2)  :: rcoord_test

          integer :: i


          test_validated = .true.

          dx = 0.1
          dy = 0.5


          !interior maps initialization
          do i=1, nx
             interior_x_map(i) = (i-1)*dx
          end do

          do i=1,ny
             interior_y_map(i) = (i-1)*dy
          end do


          !test_get_border_detector_N
          icoord = [nx,ny+1]
          rcoord = [(nx-1)*dx,ny*dy]

          icoord_test = [nx,dct_icr_N_default]
          rcoord_test = [(nx-1)*dx,interior_y_map(dct_icr_N_default)]

          call get_border_detector_N(
     $         icoord,
     $         rcoord,
     $         interior_y_map,
     $         icoord_n,
     $         rcoord_n)

          test_loc = icoord_n(1).eq.icoord_test(1)
          test_loc = test_loc.and.
     $         (icoord_n(2).eq.icoord_test(2))
          test_loc = test_loc.and.
     $         is_test_validated(rcoord_n(1),rcoord_test(1),detailled)
          test_loc = test_loc.and.
     $         is_test_validated(rcoord_n(2),rcoord_test(2),detailled)
          
          if(detailled) then
             print '(''test_get_border_detector_N: '',L1)', test_loc
          end if
          test_validated = test_validated.and.test_loc


          !test_get_border_detector_S
          icoord = [nx,0]
          rcoord = [(nx-1)*dx,-dy]

          icoord_test = [nx,dct_icr_S_default]
          rcoord_test = [(nx-1)*dx,interior_y_map(dct_icr_S_default)]

          call get_border_detector_S(
     $         icoord,
     $         rcoord,
     $         interior_y_map,
     $         icoord_n,
     $         rcoord_n)

          test_loc = icoord_n(1).eq.icoord_test(1)
          test_loc = test_loc.and.
     $         (icoord_n(2).eq.icoord_test(2))
          test_loc = test_loc.and.
     $         is_test_validated(rcoord_n(1),rcoord_test(1),detailled)
          test_loc = test_loc.and.
     $         is_test_validated(rcoord_n(2),rcoord_test(2),detailled)
          
          if(detailled) then
             print '(''test_get_border_detector_S: '',L1)', test_loc
          end if
          test_validated = test_validated.and.test_loc


          !test_get_border_detector_E
          icoord = [nx+1 ,5   ]
          rcoord = [nx*dx,4*dy]

          icoord_test = [dct_icr_E_default,5]
          rcoord_test = [interior_x_map(dct_icr_E_default),4*dy]

          call get_border_detector_E(
     $         icoord,
     $         rcoord,
     $         interior_x_map,
     $         icoord_n,
     $         rcoord_n)

          test_loc = icoord_n(1).eq.icoord_test(1)
          test_loc = test_loc.and.
     $         (icoord_n(2).eq.icoord_test(2))
          test_loc = test_loc.and.
     $         is_test_validated(rcoord_n(1),rcoord_test(1),detailled)
          test_loc = test_loc.and.
     $         is_test_validated(rcoord_n(2),rcoord_test(2),detailled)
          
          if(detailled) then
             print '(''test_get_border_detector_E: '',L1)', test_loc
          end if
          test_validated = test_validated.and.test_loc


          !test_get_border_detector_W
          icoord = [0  ,5   ]
          rcoord = [-dx,4*dy]

          icoord_test = [dct_icr_W_default,5]
          rcoord_test = [interior_x_map(dct_icr_W_default),4*dy]

          call get_border_detector_W(
     $         icoord,
     $         rcoord,
     $         interior_x_map,
     $         icoord_n,
     $         rcoord_n)

          test_loc = icoord_n(1).eq.icoord_test(1)
          test_loc = test_loc.and.
     $         (icoord_n(2).eq.icoord_test(2))
          test_loc = test_loc.and.
     $         is_test_validated(rcoord_n(1),rcoord_test(1),detailled)
          test_loc = test_loc.and.
     $         is_test_validated(rcoord_n(2),rcoord_test(2),detailled)
          
          if(detailled) then
             print '(''test_get_border_detector_W: '',L1)', test_loc
          end if
          test_validated = test_validated.and.test_loc

        end function test_get_border_detector


        function test_should_be_removed(detailled) result(test_validated)
          
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          logical                        :: test_loc
          integer(ikind), dimension(2,2) :: bf_align
          integer(ikind), dimension(2,5) :: test_gcoords
          logical       , dimension(5)   :: test_remove
          logical                        :: remove
          integer                        :: k
          logical                        :: test_N
          logical                        :: test_S
          logical                        :: test_E
          logical                        :: test_W
          

          !test should_be_removed_N
          !============================================================
          test_validated = .true.

          bf_align = reshape((/
     $         7,ny-1,
     $         10,ny+5/),
     $         (/2,2/))

          test_gcoords = reshape((/
     $         2 ,ny+1,
     $         5 ,ny+1,
     $         5 ,ny-2,
     $         12, ny,
     $         13, ny/),
     $         (/2,5/))

          test_remove = [.false.,.true.,.true.,.true.,.false.]

          do k=1,size(test_remove,1)
             
             remove = should_be_removed_N(
     $            bf_align,
     $            test_gcoords(:,k))

             test_loc = remove.eqv.test_remove(k)
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''** test '',I2,'' failed **'')', k
             end if

          end do

          test_N = test_validated

          if(detailled) then
             print '(''test_should_be_removed_N: '',L1)', test_N
          end if


          !test should_be_removed_N
          !============================================================
          test_validated = .true.
          bf_align = reshape((/
     $         7,-5,
     $         10,2/),
     $         (/2,2/))

          test_gcoords = reshape((/
     $         2 ,1,
     $         5 ,1,
     $         5 ,2,
     $         12, 2,
     $         13, 2/),
     $         (/2,5/))

          test_remove = [.false.,.true.,.true.,.true.,.false.]

          do k=1,size(test_remove,1)
             
             remove = should_be_removed_S(
     $            bf_align,
     $            test_gcoords(:,k))

             test_loc = remove.eqv.test_remove(k)
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''** test '',I2,'' failed **'')', k
             end if

          end do

          test_S = test_validated

          if(detailled) then
             print '(''test_should_be_removed_S: '',L1)', test_S
          end if


          !test should_be_removed_W
          !============================================================
          test_validated = .true.
          bf_align = reshape((/
     $         -5,7,
     $          2,20 /),
     $         (/2,2/))

          test_gcoords = reshape((/
     $         1, 2 ,
     $         1, 5 ,
     $         2, 5 ,
     $         2, 12,
     $         2, 20/),
     $         (/2,5/))

          test_remove = [.false.,.true.,.true.,.true.,.true.]

          do k=1,size(test_remove,1)
             
             remove = should_be_removed_W(
     $            bf_align,
     $            test_gcoords(:,k))

             test_loc = remove.eqv.test_remove(k)
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''** test '',I2,'' failed **'')', k
             end if

          end do

          test_W = test_validated

          if(detailled) then
             print '(''test_should_be_removed_W: '',L1)', test_W
          end if


          !test should_be_removed_E
          !============================================================
          test_validated = .true.
          bf_align = reshape((/
     $         19,7,
     $         25,20 /),
     $         (/2,2/))

          test_gcoords = reshape((/
     $         nx  , 2 ,
     $         nx  , 5 ,
     $         nx-1, 5 ,
     $         nx-1, 12,
     $         nx-1, 20/),
     $         (/2,5/))

          test_remove = [.false.,.true.,.true.,.true.,.true.]

          do k=1,size(test_remove,1)
             
             remove = should_be_removed_E(
     $            bf_align,
     $            test_gcoords(:,k))

             test_loc = remove.eqv.test_remove(k)
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''** test '',I2,'' failed **'')', k
             end if

          end do

          test_E = test_validated

          if(detailled) then
             print '(''test_should_be_removed_E: '',L1)', test_E
          end if
          

          test_validated = test_N.and.test_S.and.test_W.and.test_E

        end function test_should_be_removed


        function test_get_rcoords(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          logical                      :: test_loc
          integer                      :: k
          real(rkind)                  :: rcoord
          real(rkind)   , dimension(5) :: interior_map
          integer(ikind), dimension(3) :: test_i
          real(rkind)   , dimension(3) :: test_rcoord


          test_validated = .true.

          interior_map = [0.1d0, 0.5d0, 0.7d0, 0.8d0, 1.2d0]
          
          test_i       = [3    ,-2    ,6    ]
          test_rcoord  = [0.7d0,-1.1d0,1.6d0]
          
          do k=1, size(test_rcoord,1)

             rcoord = get_rcoord(test_i(k),interior_map)

             test_loc = is_test_validated(
     $            rcoord,
     $            test_rcoord(k),
     $            detailled)

             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''** test '',I2,'' failed **'')', k
             end if

          end do

        end function test_get_rcoords


        function test_compute_new_list_param(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          type(bf_detector_dcr_param)     :: bf_detector_dcr_param_used
          logical                         :: test_loc

          real(rkind)   , dimension(nx)   :: interior_x_map
          real(rkind)   , dimension(ny)   :: interior_y_map

          integer       , dimension(4)    :: bf_localization
          integer(ikind), dimension(2,11) :: icoords
          real(rkind)   , dimension(2,11) :: rcoords
          integer(ikind), dimension(2,2)  :: bf_align

          integer(ikind), dimension(2)    :: segment_borders_test
          integer(ikind), dimension(2)    :: first_icoord_test
          real(rkind)   , dimension(2)    :: first_rcoord_test
          integer(ikind), dimension(2)    :: last_icoord_test
          real(rkind)   , dimension(2)    :: last_rcoord_test 

          integer :: k,config


          test_validated = .true.


          !initialize the interior maps
          interior_x_map = [
     $         0.0d0,0.2d0,0.5d0,0.8d0,1.0d0,2.3d0,8.9d0,8.99d0,9.2d0,12.5d0,
     $         14.3d0,18.6d0,23.4d0,30.2d0,40.6d0,45.6d0,50.12d0,56.1d0,57.2d0,60.2d0]

          interior_y_map = [
     $         0.1d0,0.15d0,0.3d0,0.65d0,1.2d0,2.1d0,5.9d0,6.99d0,8.2d0,10.5d0,
     $         12.3d0,15.6d0,21.4d0,28.2d0,36.6d0,41.6d0,47.12d0,52.1d0,53.2d0,63.2d0]

          
          bf_localization = [N,S,E,W]


          do k=1,4

             if(detailled) then
                print '(''test bf_localization: '',I2)',
     $               bf_localization(k)
             end if

             do config=1,4

                !get the test data
                call get_data_test_compute_new_list_param(
     $               bf_localization(k), config,
     $               interior_x_map,
     $               interior_y_map,
     $               bf_align,
     $               icoords,
     $               rcoords,
     $               segment_borders_test,
     $               first_icoord_test,
     $               first_rcoord_test,
     $               last_icoord_test,
     $               last_rcoord_test)

                !apply compute_new_list_param
                call bf_detector_dcr_param_used%compute_new_list_param(
     $               bf_localization(k),
     $               bf_align,
     $               interior_x_map,
     $               interior_y_map,
     $               icoords,
     $               rcoords)
                
             !compare the results
             test_loc = compare_bf_detector_param(
     $            bf_detector_dcr_param_used,
     $            segment_borders_test,
     $            first_icoord_test,
     $            first_rcoord_test,
     $            last_icoord_test,
     $            last_rcoord_test,
     $            detailled)

             if(detailled.and.(.not.test_loc)) then
                print '(''** test '',I2,'' failed **'')', config
             end if

             !update test validation
             test_validated = test_validated.and.test_loc

          end do
       end do

      end function test_compute_new_list_param


      function test_finalize_new_list(detailled)
     $     result(test_validated)

        implicit none

        logical, intent(in) :: detailled
        logical             :: test_validated


        type(bf_detector_dcr_param)     :: bf_detector_dcr_param_used
        logical                         :: test_loc

        real(rkind)   , dimension(nx)   :: interior_x_map
        real(rkind)   , dimension(ny)   :: interior_y_map

        integer       , dimension(2)                :: bf_localization
        integer(ikind), dimension(2,2)              :: bf_align
        integer(ikind), dimension(:,:), allocatable :: icoords_input
        real(rkind)   , dimension(:,:), allocatable :: rcoords_input
        integer(ikind), dimension(:,:), allocatable :: icoords_output
        real(rkind)   , dimension(:,:), allocatable :: rcoords_output
        integer(ikind), dimension(:,:), allocatable :: icoords_output_test
        real(rkind)   , dimension(:,:), allocatable :: rcoords_output_test
        integer(ikind), dimension(:,:), allocatable :: tmp_icoords
        real(rkind)   , dimension(:,:), allocatable :: tmp_rcoords
        integer(ikind), dimension(2)                :: left_icoord_test
        real(rkind)   , dimension(2)                :: left_rcoord_test
        integer(ikind), dimension(2)                :: right_icoord_test
        real(rkind)   , dimension(2)                :: right_rcoord_test

        integer :: k,l
        integer :: config_detectors
        integer :: config_extra_pts

        logical, dimension(2,3) :: remove_test
        logical                 :: remove_first_dct
        logical                 :: remove_last_dct

        test_validated = .true.


        !initialize the interior maps
        interior_x_map = [
     $       0.0d0,0.2d0,0.5d0,0.8d0,1.0d0,2.3d0,8.9d0,8.99d0,9.2d0,12.5d0,
     $       14.3d0,18.6d0,23.4d0,30.2d0,40.6d0,45.6d0,50.12d0,56.1d0,57.2d0,60.2d0]

        interior_y_map = [
     $       0.1d0,0.15d0,0.3d0,0.65d0,1.2d0,2.1d0,5.9d0,6.99d0,8.2d0,10.5d0,
     $       12.3d0,15.6d0,21.4d0,28.2d0,36.6d0,41.6d0,47.12d0,52.1d0,53.2d0,63.2d0]

        bf_localization = [N,E]


        !============================================================
        !test config: no_extra_pts
        !============================================================
        if(detailled) then
           print '(''test config: no_extra_pts'')'
        end if        

        config_extra_pts = 1

        do k=1,2

           if(detailled) then
              print '(''test bf_localization: '',I2)',
     $             bf_localization(k)
           end if

           do config_detectors=1,4
              
              !allocate input detectors
              allocate(icoords_input(2,11))
              allocate(rcoords_input(2,11))

              !get the test data
              call get_data_test_finalize_new_list(
     $             bf_localization(k),
     $             config_detectors,
     $             config_extra_pts,
     $             interior_x_map,
     $             interior_y_map,
     $             bf_align,
     $             icoords_input,
     $             rcoords_input,
     $             icoords_output_test,
     $             rcoords_output_test,
     $             left_icoord_test,
     $             left_rcoord_test,
     $             right_icoord_test,
     $             right_rcoord_test)

              !use input
              allocate(icoords_output(2,11))
              allocate(rcoords_output(2,11))

              icoords_output = icoords_input
              rcoords_output = rcoords_input

              !compute_new_list_param
              call bf_detector_dcr_param_used%compute_new_list_param(
     $             bf_localization(k),
     $             bf_align,
     $             interior_x_map,
     $             interior_y_map,
     $             icoords_input,
     $             rcoords_input)

              !finalize the new list
              call bf_detector_dcr_param_used%finalize_new_list(
     $             bf_localization(k),
     $             interior_x_map,
     $             interior_y_map,
     $             icoords_output,
     $             rcoords_output,
     $             left_icoord_test,
     $             left_rcoord_test,
     $             right_icoord_test,
     $             right_rcoord_test,
     $             .false.,
     $             .false.)
              
              !compare the results
              test_loc = compare_bf_detector_new_lists(
     $             icoords_output,
     $             rcoords_output,
     $             icoords_output_test,
     $             rcoords_output_test,
     $             detailled)

              if(detailled.and.(.not.test_loc)) then
                 print '(''** test '',I2,'' failed **'')', config_detectors
              end if

              !update test validation
              test_validated = test_validated.and.test_loc

              !deallocate
              deallocate(icoords_input)
              deallocate(rcoords_input)
              deallocate(icoords_output)
              deallocate(rcoords_output)
              deallocate(icoords_output_test)
              deallocate(rcoords_output_test)

           end do
        end do


        !============================================================
        !test config: add_extra_pts
        !============================================================
        config_extra_pts = 2

        do k=1,2

           do config_detectors=1,1
              
              !allocate input detectors
              allocate(icoords_input(2,11))
              allocate(rcoords_input(2,11))

              !get the test data
              call get_data_test_finalize_new_list(
     $             bf_localization(k),
     $             config_detectors,
     $             config_extra_pts,
     $             interior_x_map,
     $             interior_y_map,
     $             bf_align,
     $             icoords_input,
     $             rcoords_input,
     $             icoords_output_test,
     $             rcoords_output_test,
     $             left_icoord_test,
     $             left_rcoord_test,
     $             right_icoord_test,
     $             right_rcoord_test)

              !use input
              allocate(icoords_output(2,11))
              allocate(rcoords_output(2,11))

              icoords_output = icoords_input
              rcoords_output = rcoords_input

              !compute_new_list_param
              call bf_detector_dcr_param_used%compute_new_list_param(
     $             bf_localization(k),
     $             bf_align,
     $             interior_x_map,
     $             interior_y_map,
     $             icoords_input,
     $             rcoords_input)

              !finalize the new list
              call bf_detector_dcr_param_used%finalize_new_list(
     $             bf_localization(k),
     $             interior_x_map,
     $             interior_y_map,
     $             icoords_output,
     $             rcoords_output,
     $             left_icoord_test,
     $             left_rcoord_test,
     $             right_icoord_test,
     $             right_rcoord_test,
     $             .false.,
     $             .false.)
              
              !compare the results
              test_loc = compare_bf_detector_new_lists(
     $             icoords_output,
     $             rcoords_output,
     $             icoords_output_test,
     $             rcoords_output_test,
     $             detailled)

              if(detailled.and.(.not.test_loc)) then
                 print '(''** test '',I2,'' failed **'')', config_detectors
              end if

              !update test validation
              test_validated = test_validated.and.test_loc

              !deallocate
              deallocate(icoords_input)
              deallocate(rcoords_input)
              deallocate(icoords_output)
              deallocate(rcoords_output)
              deallocate(icoords_output_test)
              deallocate(rcoords_output_test)

           end do
        end do


        !============================================================
        !test config: remove_dct
        !============================================================
        if(detailled) then
           print '(''test config: no_extra_pts'')'
        end if        

        config_extra_pts = 1

        remove_test = reshape((/
     $       .true.,.true.,
     $       .true.,.false.,
     $       .false.,.true./),
     $       (/2,3/))

        do l=1,3

           remove_first_dct = remove_test(1,l)
           remove_last_dct  = remove_test(2,l)

           do k=1,2
              
              if(detailled) then
                 print '(''test bf_localization: '',I2)',
     $                bf_localization(k)
              end if
              
              do config_detectors=1,4
                 
                 !allocate input detectors
                 allocate(icoords_input(2,11))
                 allocate(rcoords_input(2,11))

                 !get the test data
                 call get_data_test_finalize_new_list(
     $                bf_localization(k),
     $                config_detectors,
     $                config_extra_pts,
     $                interior_x_map,
     $                interior_y_map,
     $                bf_align,
     $                icoords_input,
     $                rcoords_input,
     $                icoords_output_test,
     $                rcoords_output_test,
     $                left_icoord_test,
     $                left_rcoord_test,
     $                right_icoord_test,
     $                right_rcoord_test)

                 !use input
                 allocate(icoords_output(2,11))
                 allocate(rcoords_output(2,11))

                 icoords_output = icoords_input
                 rcoords_output = rcoords_input

                 !reallocation for the test of remove_first_dct
                 if(remove_first_dct) then

                    if(remove_last_dct) then

                       allocate(tmp_icoords(2,size(icoords_output_test,2)-2))
                       allocate(tmp_rcoords(2,size(icoords_output_test,2)-2))
                       tmp_icoords = icoords_output_test(:,2:size(icoords_output_test,2)-1)
                       tmp_rcoords = rcoords_output_test(:,2:size(rcoords_output_test,2)-1)
                       call MOVE_ALLOC(tmp_icoords,icoords_output_test)
                       call MOVE_ALLOC(tmp_rcoords,rcoords_output_test)

                    else

                       allocate(tmp_icoords(2,size(icoords_output_test,2)-1))
                       allocate(tmp_rcoords(2,size(icoords_output_test,2)-1))
                       tmp_icoords = icoords_output_test(:,2:size(icoords_output_test,2))
                       tmp_rcoords = rcoords_output_test(:,2:size(rcoords_output_test,2))
                       call MOVE_ALLOC(tmp_icoords,icoords_output_test)
                       call MOVE_ALLOC(tmp_rcoords,rcoords_output_test)

                    end if

                 else

                    allocate(tmp_icoords(2,size(icoords_output_test,2)-1))
                    allocate(tmp_rcoords(2,size(icoords_output_test,2)-1))
                    tmp_icoords = icoords_output_test(:,1:size(icoords_output_test,2)-1)
                    tmp_rcoords = rcoords_output_test(:,1:size(rcoords_output_test,2)-1)
                    call MOVE_ALLOC(tmp_icoords,icoords_output_test)
                    call MOVE_ALLOC(tmp_rcoords,rcoords_output_test)
                    
                 end if
                 
                 !compute_new_list_param
                 call bf_detector_dcr_param_used%compute_new_list_param(
     $                bf_localization(k),
     $                bf_align,
     $                interior_x_map,
     $                interior_y_map,
     $                icoords_input,
     $                rcoords_input)

                 !finalize the new list
                 call bf_detector_dcr_param_used%finalize_new_list(
     $                bf_localization(k),
     $                interior_x_map,
     $                interior_y_map,
     $                icoords_output,
     $                rcoords_output,
     $                left_icoord_test,
     $                left_rcoord_test,
     $                right_icoord_test,
     $                right_rcoord_test,
     $                remove_first_dct,
     $                remove_last_dct)
                 
                 !compare the results
                 test_loc = compare_bf_detector_new_lists(
     $                icoords_output,
     $                rcoords_output,
     $                icoords_output_test,
     $                rcoords_output_test,
     $                detailled)

                 if(detailled.and.(.not.test_loc)) then
                    print '(''** test '',I2,'' failed **'')', config_detectors
                 end if

                 !update test validation
                 test_validated = test_validated.and.test_loc

                 !deallocate
                 deallocate(icoords_input)
                 deallocate(rcoords_input)
                 deallocate(icoords_output)
                 deallocate(rcoords_output)
                 deallocate(icoords_output_test)
                 deallocate(rcoords_output_test)

              end do
           end do
        end do

      end function test_finalize_new_list


      function compare_bf_detector_param(
     $     bf_detector_dcr_param_used,
     $     segment_borders_test,
     $     first_icoord_test,
     $     first_rcoord_test,
     $     last_icoord_test,
     $     last_rcoord_test,
     $     detailled)
     $     result(test_validated)

        implicit none

        type(bf_detector_dcr_param) , intent(in) :: bf_detector_dcr_param_used
        integer(ikind), dimension(2), intent(in) :: segment_borders_test
        integer(ikind), dimension(2), intent(in) :: first_icoord_test
        real(rkind)   , dimension(2), intent(in) :: first_rcoord_test
        integer(ikind), dimension(2), intent(in) :: last_icoord_test
        real(rkind)   , dimension(2), intent(in) :: last_rcoord_test
        logical                     , intent(in) :: detailled
        logical                                  :: test_validated

        integer(ikind), dimension(2) :: segment_borders
        integer(ikind), dimension(2) :: first_icoord
        real(rkind)   , dimension(2) :: first_rcoord
        integer(ikind), dimension(2) :: last_icoord
        real(rkind)   , dimension(2) :: last_rcoord

        logical :: test_loc

        test_validated = .true.

        call bf_detector_dcr_param_used%get_param(
     $         segment_borders,
     $         first_icoord,
     $         first_rcoord,
     $         last_icoord,
     $         last_rcoord)

        test_loc =
     $       (segment_borders(1).eq.segment_borders_test(1)).and.
     $       (segment_borders(2).eq.segment_borders_test(2))
        if(detailled.and.(.not.test_loc)) then
           print '('' - segment_borders: failed'')'
           print '(''   '',I2,'' -> '',I2)', segment_borders(1), segment_borders_test(1)
           print '(''   '',I2,'' -> '',I2)', segment_borders(2), segment_borders_test(2)
        end if
        test_validated = test_validated.and.test_loc
        
        
        test_loc =
     $       (first_icoord(1).eq.first_icoord_test(1)).and.
     $       (first_icoord(2).eq.first_icoord_test(2))
        if(detailled.and.(.not.test_loc)) then
           print '('' - first_icoord: failed'')'
           print '(''   '',I2,'' -> '',I2)', first_icoord(1), first_icoord_test(1)
           print '(''   '',I2,'' -> '',I2)', first_icoord(2), first_icoord_test(2)
        end if
        test_validated = test_validated.and.test_loc


        test_loc =
     $       (is_test_validated(first_rcoord(1),first_rcoord_test(1),.false.)).and.
     $       (is_test_validated(first_rcoord(2),first_rcoord_test(2),.false.))
        if(detailled.and.(.not.test_loc)) then
           print '('' - first_rcoord: failed'')'
           print '(''   '',F8.2,'' -> '',F8.2)', first_rcoord(1), first_rcoord_test(1)
           print '(''   '',F8.2,'' -> '',F8.2)', first_rcoord(2), first_rcoord_test(2)
        end if
        test_validated = test_validated.and.test_loc


        test_loc =
     $       (last_icoord(1).eq.last_icoord_test(1)).and.
     $       (last_icoord(2).eq.last_icoord_test(2))
        if(detailled.and.(.not.test_loc)) then
           print '('' - last_icoord: failed'')'
           print '(''   '',I2,'' -> '',I2)', last_icoord(1), last_icoord_test(1)
           print '(''   '',I2,'' -> '',I2)', last_icoord(2), last_icoord_test(2)
        end if
        test_validated = test_validated.and.test_loc


        test_loc =
     $       (is_test_validated(last_rcoord(1),last_rcoord_test(1),.false.)).and.
     $       (is_test_validated(last_rcoord(2),last_rcoord_test(2),.false.))
        if(detailled.and.(.not.test_loc)) then
           print '('' - last_rcoord: failed'')'
           print '(''   '',F8.2,'' -> '',F8.2)', last_rcoord(1), last_rcoord_test(1)
           print '(''   '',F8.2,'' -> '',F8.2)', last_rcoord(2), last_rcoord_test(2)
        end if
        test_validated = test_validated.and.test_loc

      end function compare_bf_detector_param


      subroutine get_data_test_compute_new_list_param(
     $     bf_localization,
     $     config,
     $     interior_x_map,
     $     interior_y_map,
     $     bf_align,
     $     icoords,
     $     rcoords,
     $     segment_borders_test,
     $     first_icoord_test,
     $     first_rcoord_test,
     $     last_icoord_test,
     $     last_rcoord_test)

         implicit none

         integer                        , intent(in)  :: bf_localization
         integer                        , intent(in)  :: config
         real(rkind)   , dimension(nx)  , intent(in)  :: interior_x_map
         real(rkind)   , dimension(ny)  , intent(in)  :: interior_y_map
         integer(ikind), dimension(2,2) , intent(out) :: bf_align
         integer(ikind), dimension(2,11), intent(out) :: icoords
         real(rkind)   , dimension(2,11), intent(out) :: rcoords
         integer(ikind), dimension(2)   , intent(out) :: segment_borders_test
         integer(ikind), dimension(2)   , intent(out) :: first_icoord_test
         real(rkind)   , dimension(2)   , intent(out) :: first_rcoord_test
         integer(ikind), dimension(2)   , intent(out) :: last_icoord_test
         real(rkind)   , dimension(2)   , intent(out) :: last_rcoord_test

         integer :: i

         integer, parameter :: no_detector_removed=1
         integer, parameter :: middle_detectors_removed=2
         integer, parameter :: middle_and_left_detectors_removed=3
         integer, parameter :: middle_and_right_detectors_removed=4
         integer, parameter :: all_detectors_removed=5

         
         select case(bf_localization)

           !North buffer layer
           !==========================================================
           case(N)

              icoords = reshape((/
     $             5,ny-2,
     $             5,ny-1,
     $             6,ny-1,
     $             6,ny,
     $             7,ny,
     $             8,ny,
     $             9,ny,
     $             10,ny,
     $             11,ny,
     $             11,ny-1,
     $             11,ny-2/),
     $             (/2,11/))

              do i=1, size(icoords,2)
                 rcoords(1,i) = interior_x_map(icoords(1,i))
                 rcoords(2,i) = interior_y_map(icoords(2,i))
              end do


              select case(config)

                case(no_detector_removed)

                   bf_align = reshape((/
     $                  14,ny-1,
     $                  16,ny+5/),
     $                  (/2,2/))

                   segment_borders_test = [0,0]
                   
                   first_icoord_test = icoords(:,1)
                   first_rcoord_test = rcoords(:,1)
                   last_icoord_test  = icoords(:,size(icoords,2))
                   last_rcoord_test  = rcoords(:,size(rcoords,2))


                case(middle_detectors_removed)
                   
                   bf_align = reshape((/
     $                  8,ny-1,
     $                  8,ny+5/),
     $                  (/2,2/))

                   segment_borders_test = [3,8]
                   
                   first_icoord_test = icoords(:,1)
                   first_rcoord_test = rcoords(:,1)
                   last_icoord_test  = icoords(:,size(icoords,2))
                   last_rcoord_test  = rcoords(:,size(rcoords,2))


                case(middle_and_left_detectors_removed)
                   
                   bf_align = reshape((/
     $                  7,ny-1,
     $                  8,ny+5/),
     $                  (/2,2/))

                   segment_borders_test = [1,8]
                   
                   first_icoord_test    = icoords(:,1)
                   first_rcoord_test    = rcoords(:,1)
                   first_icoord_test(2) = dct_icr_N_default
                   first_rcoord_test(2) = interior_y_map(dct_icr_N_default)
                   last_icoord_test     = icoords(:,size(icoords,2))
                   last_rcoord_test     = rcoords(:,size(rcoords,2))


               case(middle_and_right_detectors_removed)

                  bf_align = reshape((/
     $                 8,ny-1,
     $                 9,ny+5/),
     $                 (/2,2/))

                  segment_borders_test = [3,11]
                  
                  first_icoord_test   = icoords(:,1)
                  first_rcoord_test   = rcoords(:,1)
                  last_icoord_test    = icoords(:,size(icoords,2))
                  last_rcoord_test    = rcoords(:,size(rcoords,2))
                  last_icoord_test(2) = dct_icr_N_default
                  last_rcoord_test(2) = interior_y_map(dct_icr_N_default)


               case(all_detectors_removed)
                   
                  bf_align = reshape((/
     $                 7,ny-1,
     $                 9,ny+5/),
     $                 (/2,2/))

                  segment_borders_test = [1,11]
                  
                  first_icoord_test    = icoords(:,1)
                  first_rcoord_test    = rcoords(:,1)
                  first_icoord_test(2) = dct_icr_N_default
                  first_rcoord_test(2) = interior_y_map(dct_icr_N_default)

                  last_icoord_test     = icoords(:,size(icoords,2))
                  last_rcoord_test     = rcoords(:,size(rcoords,2))
                  last_icoord_test(2)  = dct_icr_N_default
                  last_rcoord_test(2)  = interior_y_map(dct_icr_N_default)
                   
               case default
                  print '(''test_bf_detector_dcr_param'')'
                  print '(''get_data_test_compute_new_list_param'')'
                  print '(''config not recognized: '',I2)', config
                  stop ''

             end select

           !South buffer layer
           !==========================================================
           case(S)

              icoords = reshape((/
     $             5,3,
     $             5,2,
     $             6,2,
     $             6,1,
     $             7,1,
     $             8,1,
     $             9,1,
     $             10,1,
     $             11,1,
     $             11,2,
     $             11,3/),
     $             (/2,11/))

              do i=1, size(icoords,2)
                 rcoords(1,i) = interior_x_map(icoords(1,i))
                 rcoords(2,i) = interior_y_map(icoords(2,i))
              end do


              select case(config)

                case(no_detector_removed)

                   bf_align = reshape((/
     $                  14,-5,
     $                  16, 2/),
     $                  (/2,2/))

                   segment_borders_test = [0,0]
                   
                   first_icoord_test = icoords(:,1)
                   first_rcoord_test = rcoords(:,1)
                   last_icoord_test  = icoords(:,size(icoords,2))
                   last_rcoord_test  = rcoords(:,size(rcoords,2))


                case(middle_detectors_removed)
                   
                   bf_align = reshape((/
     $                  8,-5,
     $                  8, 2/),
     $                  (/2,2/))

                   segment_borders_test = [3,8]
                   
                   first_icoord_test = icoords(:,1)
                   first_rcoord_test = rcoords(:,1)
                   last_icoord_test  = icoords(:,size(icoords,2))
                   last_rcoord_test  = rcoords(:,size(rcoords,2))


                case(middle_and_left_detectors_removed)
                   
                   bf_align = reshape((/
     $                  7,-5,
     $                  8, 2/),
     $                  (/2,2/))

                   segment_borders_test = [1,8]
                   
                   first_icoord_test    = icoords(:,1)
                   first_rcoord_test    = rcoords(:,1)
                   first_icoord_test(2) = dct_icr_S_default
                   first_rcoord_test(2) = interior_y_map(dct_icr_S_default)
                   last_icoord_test     = icoords(:,size(icoords,2))
                   last_rcoord_test     = rcoords(:,size(rcoords,2))


               case(middle_and_right_detectors_removed)

                  bf_align = reshape((/
     $                 8,-5,
     $                 9, 2/),
     $                 (/2,2/))

                  segment_borders_test = [3,11]
                  
                  first_icoord_test   = icoords(:,1)
                  first_rcoord_test   = rcoords(:,1)
                  last_icoord_test    = icoords(:,size(icoords,2))
                  last_rcoord_test    = rcoords(:,size(rcoords,2))
                  last_icoord_test(2) = dct_icr_S_default
                  last_rcoord_test(2) = interior_y_map(dct_icr_S_default)


               case(all_detectors_removed)
                   
                  bf_align = reshape((/
     $                 7,-5,
     $                 9, 2/),
     $                 (/2,2/))

                  segment_borders_test = [1,11]
                  
                  first_icoord_test    = icoords(:,1)
                  first_rcoord_test    = rcoords(:,1)
                  first_icoord_test(2) = dct_icr_S_default
                  first_rcoord_test(2) = interior_y_map(dct_icr_S_default)

                  last_icoord_test     = icoords(:,size(icoords,2))
                  last_rcoord_test     = rcoords(:,size(rcoords,2))
                  last_icoord_test(2)  = dct_icr_S_default
                  last_rcoord_test(2)  = interior_y_map(dct_icr_S_default)
                   
               case default
                  print '(''test_bf_detector_dcr_param'')'
                  print '(''get_data_test_compute_new_list_param'')'
                  print '(''config not recognized: '',I2)', config
                  stop ''

             end select


           !East buffer layer
           !==========================================================
           case(E)

              icoords = reshape((/
     $             nx-2, 5,
     $             nx-1, 5,
     $             nx-1, 6,
     $             nx  , 6,
     $             nx  , 7,
     $             nx  , 8,
     $             nx  , 9,
     $             nx  ,10,
     $             nx  ,11,
     $             nx-1,11,
     $             nx-2,11/),
     $             (/2,11/))

              do i=1, size(icoords,2)
                 rcoords(1,i) = interior_x_map(icoords(1,i))
                 rcoords(2,i) = interior_y_map(icoords(2,i))
              end do


              select case(config)

                case(no_detector_removed)

                   bf_align = reshape((/
     $                  nx-1 ,14,
     $                  nx+10,16/),
     $                  (/2,2/))

                   segment_borders_test = [0,0]
                   
                   first_icoord_test = icoords(:,1)
                   first_rcoord_test = rcoords(:,1)
                   last_icoord_test  = icoords(:,size(icoords,2))
                   last_rcoord_test  = rcoords(:,size(rcoords,2))


                case(middle_detectors_removed)
                   
                   bf_align = reshape((/
     $                  nx-1 ,8,
     $                  nx+10,8/),
     $                  (/2,2/))

                   segment_borders_test = [3,8]
                   
                   first_icoord_test = icoords(:,1)
                   first_rcoord_test = rcoords(:,1)
                   last_icoord_test  = icoords(:,size(icoords,2))
                   last_rcoord_test  = rcoords(:,size(rcoords,2))


                case(middle_and_left_detectors_removed)
                   
                   bf_align = reshape((/
     $                  nx-1 ,7,
     $                  nx+10,8/),
     $                  (/2,2/))

                   segment_borders_test = [1,8]
                   
                   first_icoord_test    = icoords(:,1)
                   first_rcoord_test    = rcoords(:,1)
                   first_icoord_test(1) = dct_icr_E_default
                   first_rcoord_test(1) = interior_x_map(dct_icr_E_default)
                   last_icoord_test     = icoords(:,size(icoords,2))
                   last_rcoord_test     = rcoords(:,size(rcoords,2))


               case(middle_and_right_detectors_removed)

                  bf_align = reshape((/
     $                 nx-1 ,8,
     $                 nx+10,9/),
     $                 (/2,2/))

                  segment_borders_test = [3,11]
                  
                  first_icoord_test   = icoords(:,1)
                  first_rcoord_test   = rcoords(:,1)
                  last_icoord_test    = icoords(:,size(icoords,2))
                  last_rcoord_test    = rcoords(:,size(rcoords,2))
                  last_icoord_test(1) = dct_icr_E_default
                  last_rcoord_test(1) = interior_x_map(dct_icr_E_default)


               case(all_detectors_removed)
                   
                  bf_align = reshape((/
     $                 nx-1 ,7,
     $                 nx+10,9/),
     $                 (/2,2/))

                  segment_borders_test = [1,11]
                  
                  first_icoord_test    = icoords(:,1)
                  first_rcoord_test    = rcoords(:,1)
                  first_icoord_test(1) = dct_icr_E_default
                  first_rcoord_test(1) = interior_x_map(dct_icr_E_default)

                  last_icoord_test     = icoords(:,size(icoords,2))
                  last_rcoord_test     = rcoords(:,size(rcoords,2))
                  last_icoord_test(1)  = dct_icr_E_default
                  last_rcoord_test(1)  = interior_x_map(dct_icr_E_default)
                   
               case default
                  print '(''test_bf_detector_dcr_param'')'
                  print '(''get_data_test_compute_new_list_param'')'
                  print '(''config not recognized: '',I2)', config
                  stop ''

               end select

           !West buffer layer
           !==========================================================
           case(W)

              icoords = reshape((/
     $             3,  5,
     $             2,  5,
     $             2,  6,
     $             1,  6,
     $             1,  7,
     $             1,  8,
     $             1,  9,
     $             1, 10,
     $             1, 11,
     $             2, 11,
     $             3, 11/),
     $             (/2,11/))

              do i=1, size(icoords,2)
                 rcoords(1,i) = interior_x_map(icoords(1,i))
                 rcoords(2,i) = interior_y_map(icoords(2,i))
              end do


              select case(config)

                case(no_detector_removed)

                   bf_align = reshape((/
     $                  -5,14,
     $                   2,16/),
     $                  (/2,2/))

                   segment_borders_test = [0,0]
                   
                   first_icoord_test = icoords(:,1)
                   first_rcoord_test = rcoords(:,1)
                   last_icoord_test  = icoords(:,size(icoords,2))
                   last_rcoord_test  = rcoords(:,size(rcoords,2))


                case(middle_detectors_removed)
                   
                   bf_align = reshape((/
     $                  -5,8,
     $                   2,8/),
     $                  (/2,2/))

                   segment_borders_test = [3,8]
                   
                   first_icoord_test = icoords(:,1)
                   first_rcoord_test = rcoords(:,1)
                   last_icoord_test  = icoords(:,size(icoords,2))
                   last_rcoord_test  = rcoords(:,size(rcoords,2))


                case(middle_and_left_detectors_removed)
                   
                   bf_align = reshape((/
     $                  -5,7,
     $                   2,8/),
     $                  (/2,2/))

                   segment_borders_test = [1,8]
                   
                   first_icoord_test    = icoords(:,1)
                   first_rcoord_test    = rcoords(:,1)
                   first_icoord_test(1) = dct_icr_W_default
                   first_rcoord_test(1) = interior_x_map(dct_icr_W_default)
                   last_icoord_test     = icoords(:,size(icoords,2))
                   last_rcoord_test     = rcoords(:,size(rcoords,2))


               case(middle_and_right_detectors_removed)

                  bf_align = reshape((/
     $                 -5,8,
     $                  2,9/),
     $                 (/2,2/))

                  segment_borders_test = [3,11]
                  
                  first_icoord_test   = icoords(:,1)
                  first_rcoord_test   = rcoords(:,1)
                  last_icoord_test    = icoords(:,size(icoords,2))
                  last_rcoord_test    = rcoords(:,size(rcoords,2))
                  last_icoord_test(1) = dct_icr_W_default
                  last_rcoord_test(1) = interior_x_map(dct_icr_W_default)


               case(all_detectors_removed)
                   
                  bf_align = reshape((/
     $                 -5,7,
     $                  2,9/),
     $                 (/2,2/))

                  segment_borders_test = [1,11]
                  
                  first_icoord_test    = icoords(:,1)
                  first_rcoord_test    = rcoords(:,1)
                  first_icoord_test(1) = dct_icr_W_default
                  first_rcoord_test(1) = interior_x_map(dct_icr_W_default)

                  last_icoord_test     = icoords(:,size(icoords,2))
                  last_rcoord_test     = rcoords(:,size(rcoords,2))
                  last_icoord_test(1)  = dct_icr_W_default
                  last_rcoord_test(1)  = interior_x_map(dct_icr_W_default)
                   
               case default
                  print '(''test_bf_detector_dcr_param'')'
                  print '(''get_data_test_compute_new_list_param'')'
                  print '(''config not recognized: '',I2)', config
                  stop ''

             end select


           case default
              print '(''test_bf_detctor_dcr_param'')'
              print '(''get_data_test_compute_new_list_param'')'
              print '(''bf_localization not recognized: '',I2)', bf_localization
              stop '' 
         end select 

        end subroutine get_data_test_compute_new_list_param


        subroutine get_data_test_finalize_new_list(
     $     bf_localization,
     $     config_detectors,
     $     config_extra_pts,
     $     interior_x_map,
     $     interior_y_map,
     $     bf_align,
     $     icoords_input,
     $     rcoords_input,
     $     icoords_output_test,
     $     rcoords_output_test,
     $     left_icoord_test,
     $     left_rcoord_test,
     $     right_icoord_test,
     $     right_rcoord_test)

         implicit none

         integer                                     , intent(in)  :: bf_localization
         integer                                     , intent(in)  :: config_detectors
         integer                                     , intent(in)  :: config_extra_pts
         real(rkind)   , dimension(nx)               , intent(in)  :: interior_x_map
         real(rkind)   , dimension(ny)               , intent(in)  :: interior_y_map
         integer(ikind), dimension(2,2)              , intent(out) :: bf_align
         integer(ikind), dimension(2,11)             , intent(out) :: icoords_input
         real(rkind)   , dimension(2,11)             , intent(out) :: rcoords_input
         integer(ikind), dimension(:,:) , allocatable, intent(out) :: icoords_output_test
         real(rkind)   , dimension(:,:) , allocatable, intent(out) :: rcoords_output_test
         integer(ikind), dimension(2)                , intent(out) :: left_icoord_test
         real(rkind)   , dimension(2)                , intent(out) :: left_rcoord_test
         integer(ikind), dimension(2)                , intent(out) :: right_icoord_test
         real(rkind)   , dimension(2)                , intent(out) :: right_rcoord_test

         integer :: i

         integer, parameter :: no_detector_removed=1
         integer, parameter :: middle_detectors_removed=2
         integer, parameter :: middle_and_left_detectors_removed=3
         integer, parameter :: middle_and_right_detectors_removed=4
         integer, parameter :: all_detectors_removed=5

         integer, parameter :: no_extra_pts=1
         integer, parameter :: add_extra_pts=2

         
         select case(bf_localization)

           !North buffer layer
           !==========================================================
           case(N)

              icoords_input = reshape((/
     $             5,ny-2,
     $             5,ny-1,
     $             6,ny-1,
     $             6,ny,
     $             7,ny,
     $             8,ny,
     $             9,ny,
     $             10,ny,
     $             11,ny,
     $             11,ny-1,
     $             11,ny-2/),
     $             (/2,11/))

              do i=1, size(icoords_input,2)
                 rcoords_input(1,i) = interior_x_map(icoords_input(1,i))
                 rcoords_input(2,i) = interior_y_map(icoords_input(2,i))
              end do


              select case(config_detectors)

                case(no_detector_removed)

                   bf_align = reshape((/
     $                  14,ny-1,
     $                  16,ny+5/),
     $                  (/2,2/))

                   select case(config_extra_pts)
                     case(no_extra_pts)
                        
                        allocate(icoords_output_test(2,11))
                        allocate(rcoords_output_test(2,11))
                        
                        icoords_output_test = icoords_input
                        rcoords_output_test = rcoords_input

                        left_icoord_test  = icoords_input(:,1)
                        left_rcoord_test  = rcoords_input(:,1)
                        right_icoord_test = icoords_input(:,11)
                        right_rcoord_test = rcoords_input(:,11)

                     case(add_extra_pts)

                        allocate(icoords_output_test(2,13))
                        allocate(rcoords_output_test(2,13))

                        icoords_output_test(:,1)    = [4,ny-2]
                        icoords_output_test(:,2:12) = icoords_input
                        icoords_output_test(:,13)   = [12,ny-2]

                        rcoords_output_test(1,1)    = 0.75d0
                        rcoords_output_test(2,1)    = interior_y_map(ny-2)
                        rcoords_output_test(:,2:12) = rcoords_input
                        rcoords_output_test(1,13)   = 18.85d0
                        rcoords_output_test(2,13)   = interior_y_map(ny-2)

                        left_icoord_test     = [3,ny-2]
                        left_rcoord_test(1)  = interior_x_map(3)
                        left_rcoord_test(2)  = interior_y_map(ny-2)
                        right_icoord_test    = [13,ny-2]
                        right_rcoord_test(1) = interior_x_map(13)
                        right_rcoord_test(2) = interior_y_map(ny-2)

                   end select


                case(middle_detectors_removed)
                   
                   bf_align = reshape((/
     $                  8,ny-1,
     $                  8,ny+5/),
     $                  (/2,2/))

                   allocate(icoords_output_test(2,10))
                   allocate(rcoords_output_test(2,10))
                        
                   icoords_output_test = reshape((/
     $                   5, ny-2,
     $                   5, ny-1,
     $                   6, ny-1,
     $                   7, ny-1,
     $                   8, ny-1,
     $                   9, ny-1,
     $                  10, ny-1,
     $                  11, ny,
     $                  11, ny-1,
     $                  11, ny-2/),
     $                  (/2,10/))

                   do i=1, size(icoords_output_test,2)
                      rcoords_output_test(1,i) = interior_x_map(icoords_output_test(1,i))
                      rcoords_output_test(2,i) = interior_y_map(icoords_output_test(2,i))
                   end do

                   left_icoord_test     = [4,ny-2]
                   left_rcoord_test(1)  = interior_x_map(4)
                   left_rcoord_test(2)  = interior_y_map(ny-2)
                   right_icoord_test    = [12,ny-2]
                   right_rcoord_test(1) = interior_x_map(12)
                   right_rcoord_test(2) = interior_y_map(ny-2)

                case(middle_and_left_detectors_removed)
                   
                   bf_align = reshape((/
     $                  7,ny-1,
     $                  8,ny+5/),
     $                  (/2,2/))

                   allocate(icoords_output_test(2,9))
                   allocate(rcoords_output_test(2,9))
                        
                   icoords_output_test = reshape((/
     $                   5, dct_icr_N_default,
     $                   6, dct_icr_N_default,
     $                   7, dct_icr_N_default,
     $                   8, dct_icr_N_default,
     $                   9, dct_icr_N_default,
     $                  10, dct_icr_N_default,
     $                  11, ny,
     $                  11, ny-1,
     $                  11, ny-2/),
     $                  (/2,9/))

                   do i=1, size(icoords_output_test,2)
                      rcoords_output_test(1,i) = interior_x_map(icoords_output_test(1,i))
                      rcoords_output_test(2,i) = interior_y_map(icoords_output_test(2,i))
                   end do

                   left_icoord_test     = [4,dct_icr_N_default]
                   left_rcoord_test(1)  = interior_x_map(4)
                   left_rcoord_test(2)  = interior_y_map(dct_icr_N_default)
                   right_icoord_test    = [12,ny-2]
                   right_rcoord_test(1) = interior_x_map(12)
                   right_rcoord_test(2) = interior_y_map(ny-2)


               case(middle_and_right_detectors_removed)

                  bf_align = reshape((/
     $                 8,ny-1,
     $                 9,ny+5/),
     $                 (/2,2/))

                  allocate(icoords_output_test(2,8))
                  allocate(rcoords_output_test(2,8))

                  icoords_output_test = reshape((/
     $                 5,ny-2,
     $                 5,ny-1,
     $                 6,ny-1,
     $                 7,ny-1,
     $                 8,ny-1,
     $                 9,ny-1,
     $                10,ny-1,
     $                11,ny-1/),
     $                 (/2,8/))
                  
                  do i=1, size(icoords_output_test,2)
                     rcoords_output_test(1,i) = interior_x_map(icoords_output_test(1,i))
                     rcoords_output_test(2,i) = interior_y_map(icoords_output_test(2,i))
                  end do

                  left_icoord_test     = [4,ny-2]
                  left_rcoord_test(1)  = interior_x_map(4)
                  left_rcoord_test(2)  = interior_y_map(ny-2)
                  right_icoord_test    = [12,dct_icr_N_default]
                  right_rcoord_test(1) = interior_x_map(12)
                  right_rcoord_test(2) = interior_y_map(dct_icr_N_default)


               case(all_detectors_removed)
                   
                  bf_align = reshape((/
     $                 7,ny-1,
     $                 9,ny+5/),
     $                 (/2,2/))

                  icoords_output_test = reshape((/
     $                 5,dct_icr_N_default,
     $                 6,dct_icr_N_default,
     $                 7,dct_icr_N_default,
     $                 8,dct_icr_N_default,
     $                 9,dct_icr_N_default,
     $                10,dct_icr_N_default,
     $                11,dct_icr_N_default/),
     $                 (/2,7/))
                  
                  do i=1, size(icoords_output_test,2)
                     rcoords_output_test(1,i) = interior_x_map(icoords_output_test(1,i))
                     rcoords_output_test(2,i) = interior_y_map(icoords_output_test(2,i))
                  end do

                  left_icoord_test     = [5,dct_icr_N_default]
                  left_rcoord_test(1)  = interior_x_map(5)
                  left_rcoord_test(2)  = interior_y_map(dct_icr_N_default)
                  right_icoord_test    = [11,dct_icr_N_default]
                  right_rcoord_test(1) = interior_x_map(11)
                  right_rcoord_test(2) = interior_y_map(dct_icr_N_default)

               case default
                  print '(''test_bf_detector_dcr_param'')'
                  print '(''get_data_test_compute_new_list_param'')'
                  print '(''config not recognized: '',I2)', config_detectors
                  stop ''

             end select


           !East buffer layer
           !==========================================================
           case(E)

              icoords_input = reshape((/
     $             nx-2, 5,
     $             nx-1, 5,
     $             nx-1, 6,
     $             nx  , 6,
     $             nx  , 7,
     $             nx  , 8,
     $             nx  , 9,
     $             nx  ,10,
     $             nx  ,11,
     $             nx-1,11,
     $             nx-2,11/),
     $             (/2,11/))

              do i=1, size(icoords_input,2)
                 rcoords_input(1,i) = interior_x_map(icoords_input(1,i))
                 rcoords_input(2,i) = interior_y_map(icoords_input(2,i))
              end do


              select case(config_detectors)

                case(no_detector_removed)

                   bf_align = reshape((/
     $                  nx-1,14,
     $                  nx+5,16/),
     $                  (/2,2/))

                   select case(config_extra_pts)
                     case(no_extra_pts)
                        
                        allocate(icoords_output_test(2,11))
                        allocate(rcoords_output_test(2,11))
                        
                        icoords_output_test = icoords_input
                        rcoords_output_test = rcoords_input

                        left_icoord_test  = icoords_input(:,1)
                        left_rcoord_test  = rcoords_input(:,1)
                        right_icoord_test = icoords_input(:,11)
                        right_rcoord_test = rcoords_input(:,11)

                     case(add_extra_pts)

                        allocate(icoords_output_test(2,13))
                        allocate(rcoords_output_test(2,13))

                        icoords_output_test(:,1)    = [nx-2,4]
                        icoords_output_test(:,2:12) = icoords_input
                        icoords_output_test(:,13)   = [nx-2,12]

                        rcoords_output_test(1,1)    = interior_x_map(nx-2)
                        rcoords_output_test(2,1)    = 0.75d0
                        rcoords_output_test(:,2:12) = rcoords_input
                        rcoords_output_test(1,13)   = interior_x_map(nx-2)
                        rcoords_output_test(2,13)   = 16.85d0

                        left_icoord_test     = [nx-2,3]
                        left_rcoord_test(1)  = interior_x_map(nx-2)
                        left_rcoord_test(2)  = interior_y_map(3)
                        right_icoord_test    = [nx-2,13]
                        right_rcoord_test(1) = interior_x_map(nx-2)
                        right_rcoord_test(2) = interior_y_map(13)

                   end select


                case(middle_detectors_removed)
                   
                   bf_align = reshape((/
     $                  nx-1,8,
     $                  nx+5,8/),
     $                  (/2,2/))

                   allocate(icoords_output_test(2,10))
                   allocate(rcoords_output_test(2,10))
                        
                   icoords_output_test = reshape((/
     $                  nx-2,  5, 
     $                  nx-1,  5, 
     $                  nx-1,  6, 
     $                  nx-1,  7, 
     $                  nx-1,  8, 
     $                  nx-1,  9, 
     $                  nx-1, 10, 
     $                  nx,   11, 
     $                  nx-1, 11, 
     $                  nx-2, 11 /),
     $                  (/2,10/))

                   do i=1, size(icoords_output_test,2)
                      rcoords_output_test(1,i) = interior_x_map(icoords_output_test(1,i))
                      rcoords_output_test(2,i) = interior_y_map(icoords_output_test(2,i))
                   end do

                   left_icoord_test     = [nx-2,4]
                   left_rcoord_test(1)  = interior_x_map(nx-2)
                   left_rcoord_test(2)  = interior_y_map(4)
                   right_icoord_test    = [nx-2,12]
                   right_rcoord_test(1) = interior_x_map(nx-2)
                   right_rcoord_test(2) = interior_y_map(12)

                case(middle_and_left_detectors_removed)
                   
                   bf_align = reshape((/
     $                  nx-1,7,
     $                  nx+5,8/),
     $                  (/2,2/))

                   allocate(icoords_output_test(2,9))
                   allocate(rcoords_output_test(2,9))
                        
                   icoords_output_test = reshape((/
     $                  dct_icr_E_default, 5,
     $                  dct_icr_E_default, 6,
     $                  dct_icr_E_default, 7,
     $                  dct_icr_E_default, 8,
     $                  dct_icr_E_default, 9,
     $                  dct_icr_E_default,10,
     $                  nx  , 11,
     $                  nx-1, 11,
     $                  nx-2, 11/),
     $                  (/2,9/))

                   do i=1, size(icoords_output_test,2)
                      rcoords_output_test(1,i) = interior_x_map(icoords_output_test(1,i))
                      rcoords_output_test(2,i) = interior_y_map(icoords_output_test(2,i))
                   end do

                   left_icoord_test     = [dct_icr_E_default,4]
                   left_rcoord_test(1)  = interior_x_map(dct_icr_E_default)
                   left_rcoord_test(2)  = interior_y_map(4)
                   right_icoord_test    = [nx-2,12]
                   right_rcoord_test(1) = interior_x_map(nx-2)
                   right_rcoord_test(2) = interior_y_map(12)


               case(middle_and_right_detectors_removed)

                  bf_align = reshape((/
     $                 nx-1,8,
     $                 nx+5,9/),
     $                 (/2,2/))

                  allocate(icoords_output_test(2,8))
                  allocate(rcoords_output_test(2,8))

                  icoords_output_test = reshape((/
     $                nx-2, 5,
     $                nx-1, 5,
     $                nx-1, 6,
     $                nx-1, 7,
     $                nx-1, 8,
     $                nx-1, 9,
     $                nx-1,10,
     $                nx-1,11/),
     $                 (/2,8/))
                  
                  do i=1, size(icoords_output_test,2)
                     rcoords_output_test(1,i) = interior_x_map(icoords_output_test(1,i))
                     rcoords_output_test(2,i) = interior_y_map(icoords_output_test(2,i))
                  end do

                  left_icoord_test     = [nx-2,4]
                  left_rcoord_test(1)  = interior_x_map(nx-2)
                  left_rcoord_test(2)  = interior_y_map(4)
                  right_icoord_test    = [dct_icr_E_default,12]
                  right_rcoord_test(1) = interior_x_map(dct_icr_E_default)
                  right_rcoord_test(2) = interior_y_map(12)


               case(all_detectors_removed)
                   
                  bf_align = reshape((/
     $                 nx-1,7,
     $                 nx+5,9/),
     $                 (/2,2/))

                  icoords_output_test = reshape((/
     $                dct_icr_E_default, 5,
     $                dct_icr_E_default, 6,
     $                dct_icr_E_default, 7,
     $                dct_icr_E_default, 8,
     $                dct_icr_E_default, 9,
     $                dct_icr_E_default,10,
     $                dct_icr_E_default,11/),
     $                (/2,7/))
                  
                  do i=1, size(icoords_output_test,2)
                     rcoords_output_test(1,i) = interior_x_map(icoords_output_test(1,i))
                     rcoords_output_test(2,i) = interior_y_map(icoords_output_test(2,i))
                  end do

                  left_icoord_test     = [dct_icr_E_default,5]
                  left_rcoord_test(1)  = interior_x_map(dct_icr_E_default)
                  left_rcoord_test(2)  = interior_y_map(5)
                  right_icoord_test    = [dct_icr_E_default,11]
                  right_rcoord_test(1) = interior_x_map(dct_icr_E_default)
                  right_rcoord_test(2) = interior_y_map(11)

               case default
                  print '(''test_bf_detector_dcr_param'')'
                  print '(''get_data_test_compute_new_list_param'')'
                  print '(''config not recognized: '',I2)', config_detectors
                  stop ''

             end select

             case default
                print '(''test_bf_detector_dcr_param'')'
                print '(''get_data_test_finalize_new_list'')'
                print '(''test not recognized: '',I2)', bf_localization
                stop ''                

          end select

        end subroutine get_data_test_finalize_new_list


        function compare_bf_detector_new_lists(
     $     icoords_output,
     $     rcoords_output,
     $     icoords_output_test,
     $     rcoords_output_test,
     $     detailled)
     $     result(test_validated)

          implicit none

          integer(ikind), dimension(:,:), intent(in) :: icoords_output
          real(rkind)   , dimension(:,:), intent(in) :: rcoords_output
          integer(ikind), dimension(:,:), intent(in) :: icoords_output_test
          real(rkind)   , dimension(:,:), intent(in) :: rcoords_output_test
          logical                       , intent(in) :: detailled
          logical                                    :: test_validated

          integer :: i
          logical :: test_loc


          test_validated = .true.

          test_loc = size(icoords_output,2).eq.size(icoords_output_test,2)
          if(detailled.and.(.not.test_loc)) then
             print '(''  - nb of detectors: '',I2,'' -> '',I2)',
     $            size(icoords_output,2), size(icoords_output_test,2)
          end if
          test_validated = test_validated.and.test_loc

          if(test_loc) then
             
             do i=1, size(icoords_output,2)

                !icoords
                test_loc =
     $               (icoords_output(1,i).eq.icoords_output_test(1,i)).and.
     $               (icoords_output(2,i).eq.icoords_output_test(2,i))
                if(detailled.and.(.not.test_loc)) then
                   print '(''   - icoords('',I2,''): '')',i
                   print '(''     '',I2,'' -> '',I2)', icoords_output(1,i), icoords_output_test(1,i)
                   print '(''     '',I2,'' -> '',I2)', icoords_output(2,i), icoords_output_test(2,i)
                end if
                test_validated = test_validated.and.test_loc

                !rcoords
                test_loc =
     $               (is_test_validated(rcoords_output(1,i),rcoords_output_test(1,i),.false.)).and.
     $               (is_test_validated(rcoords_output(2,i),rcoords_output_test(2,i),.false.))
                if(detailled.and.(.not.test_loc)) then
                   print '(''   - rcoords('',I2,''): '')',i
                   print '(''     '',F8.2,'' -> '',F8.2)', rcoords_output(1,i), rcoords_output_test(1,i)
                   print '(''     '',F8.2,'' -> '',F8.2)', rcoords_output(2,i), rcoords_output_test(2,i)
                end if
                test_validated = test_validated.and.test_loc                

             end do

          end if
          

        end function compare_bf_detector_new_lists


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
     $         int(var*10000.)-
     $         sign(int(abs(cst*10000.)),int(cst*10000.))).le.1
          
        end function is_test_validated        

      end program test_bf_detector_dcr_param

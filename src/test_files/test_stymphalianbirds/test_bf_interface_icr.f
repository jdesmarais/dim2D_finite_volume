      program test_bf_interface_icr

        use bf_interface_icr_class, only :
     $       bf_interface_icr

        use parameters_bf_layer, only :
     $       dct_icr_distance

        use parameters_input, only :
     $       nx,ny,dt

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        type(bf_interface_icr) :: bf_interface_icr_used
        logical                :: detailled
        logical                :: test_validated
        
        detailled = .true.

        test_validated = test_ini(detailled)
        print '(''test_ini: '',L1)', test_validated
        print '()'

        test_validated = test_get_central_grdpt(bf_interface_icr_used, detailled)
        print '(''test_bf_interface_icr: '',L1)', test_validated
        print '()'


        contains

        function test_ini(
     $       detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_interface_icr) :: bf_interface_icr_used
          
          real(rkind)   , dimension(nx) :: interior_x_map
          real(rkind)   , dimension(ny) :: interior_y_map

          integer(ikind), dimension(2,10) :: N_icoords_test
          integer(ikind), dimension(2,10) :: S_icoords_test
          integer(ikind), dimension(2,8)  :: E_icoords_test
          integer(ikind), dimension(2,8)  :: W_icoords_test
          real(rkind)   , dimension(2,10) :: N_rcoords_test
          real(rkind)   , dimension(2,10) :: S_rcoords_test
          real(rkind)   , dimension(2,8)  :: E_rcoords_test
          real(rkind)   , dimension(2,8)  :: W_rcoords_test

          integer(ikind), dimension(:,:), allocatable :: N_icoords
          integer(ikind), dimension(:,:), allocatable :: S_icoords
          integer(ikind), dimension(:,:), allocatable :: E_icoords
          integer(ikind), dimension(:,:), allocatable :: W_icoords
          real(rkind)   , dimension(:,:), allocatable :: N_rcoords
          real(rkind)   , dimension(:,:), allocatable :: S_rcoords
          real(rkind)   , dimension(:,:), allocatable :: E_rcoords
          real(rkind)   , dimension(:,:), allocatable :: W_rcoords

          integer :: i
          logical :: test_loc

          if((nx.ne.20).or.(ny.ne.20).or.(dct_icr_distance.ne.4)) then

             print '(''the test requires:'')'
             print '(''nx=20: '',L1)', (nx.ne.20)
             print '(''ny=20: '',L1)', (ny.ne.20)
             print '(''dct_icr_distance=4: '',L1)', (dct_icr_distance.ne.4)

          end if

          !initialize the interior maps
          interior_x_map = [
     $         0.0d0,0.2d0,0.5d0,0.8d0,1.0d0,
     $         .3d0,8.9d0,8.99d0,9.2d0,12.5d0,
     $         14.3d0,18.6d0,23.4d0,30.2d0,40.6d0,
     $         45.6d0,50.12d0,56.1d0,57.2d0,60.2d0]

          interior_y_map = [
     $         0.1d0,0.15d0,0.3d0,0.65d0,1.2d0,
     $         2.1d0,5.9d0,6.99d0,8.2d0,10.5d0,
     $         12.3d0,15.6d0,21.4d0,28.2d0,36.6d0,
     $         41.6d0,47.12d0,52.1d0,53.2d0,63.2d0]         

          N_icoords_test = reshape((/
     $         6,15,
     $         7,15,
     $         8,15,
     $         9,15,
     $         10,15,
     $         11,15,
     $         12,15,
     $         13,15,
     $         14,15,
     $         15,15/),
     $         (/2,10/))

          S_icoords_test = reshape((/
     $         6 ,6,
     $         7 ,6,
     $         8 ,6,
     $         9 ,6,
     $         10,6,
     $         11,6,
     $         12,6,
     $         13,6,
     $         14,6,
     $         15,6/),
     $         (/2,10/))

          E_icoords_test = reshape((/
     $         15,7 ,
     $         15,8 ,
     $         15,9 ,
     $         15,10,
     $         15,11,
     $         15,12,
     $         15,13,
     $         15,14/),
     $         (/2,8/))

          W_icoords_test = reshape((/
     $         6,7 ,
     $         6,8 ,
     $         6,9 ,
     $         6,10,
     $         6,11,
     $         6,12,
     $         6,13,
     $         6,14/),
     $         (/2,8/))

          do i=1,size(N_icoords_test,2)

             N_rcoords_test(1,i) = interior_x_map(N_icoords_test(1,i))
             S_rcoords_test(1,i) = interior_x_map(S_icoords_test(1,i))

             N_rcoords_test(2,i) = interior_y_map(N_icoords_test(2,i))
             S_rcoords_test(2,i) = interior_y_map(S_icoords_test(2,i))
             
          end do

          do i=1, size(E_icoords_test,2)

             E_rcoords_test(1,i) = interior_x_map(E_icoords_test(1,i))
             W_rcoords_test(1,i) = interior_x_map(W_icoords_test(1,i))

             E_rcoords_test(2,i) = interior_y_map(E_icoords_test(2,i))
             W_rcoords_test(2,i) = interior_y_map(W_icoords_test(2,i))

          end do


          !initialize the detector positions
          call bf_interface_icr_used%ini(interior_x_map, interior_y_map)

          
          !get the detectors initialized from the bf_interface
          call bf_interface_icr_used%get_dct_coords(
     $         N_icoords,
     $         S_icoords,
     $         E_icoords,
     $         W_icoords,
     $         N_rcoords,
     $         S_rcoords,
     $         E_rcoords,
     $         W_rcoords)

          !compare the detector initialized
          test_validated = .true.

          test_loc = compare_icoords(N_icoords,N_icoords_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = compare_icoords(S_icoords,S_icoords_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = compare_icoords(E_icoords,E_icoords_test,detailled)
          test_validated = test_validated.and.test_loc

c$$$          test_loc = compare_icoords(W_icoords,W_icoords_test,detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$
c$$$          test_loc = compare_rcoords(N_rcoords,N_rcoords_test,detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$
c$$$          test_loc = compare_rcoords(S_rcoords,S_rcoords_test,detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$
c$$$          test_loc = compare_rcoords(E_rcoords,E_rcoords_test,detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$
c$$$          test_loc = compare_rcoords(W_rcoords,W_rcoords_test,detailled)
c$$$          test_validated = test_validated.and.test_loc

        end function test_ini

      
        function compare_icoords(icoords, icoords_test, detailled)
     $     result(test_validated)

          implicit none

          integer(ikind), dimension(:,:), intent(in) :: icoords
          integer(ikind), dimension(:,:), intent(in) :: icoords_test
          logical                       , intent(in) :: detailled
          logical                                    :: test_validated

          integer :: i,j
          logical :: test_loc

          test_validated = .true.

          do j=1, size(icoords,2)
             do i=1, size(icoords,1)

                test_loc = icoords(i,j).eq.icoords_test(i,j)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test failed at: '',2I2)',i,j
                   print '(''  - '',I2,''->'',I2)', icoords(i,j), icoords_test(i,j)
                   
                end if

             end do
          end do

        end function compare_icoords


        function compare_rcoords(rcoords, rcoords_test, detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(:,:), intent(in) :: rcoords
          real(rkind), dimension(:,:), intent(in) :: rcoords_test
          logical                    , intent(in) :: detailled
          logical                                 :: test_validated

          integer :: i,j
          logical :: test_loc

          test_validated = .true.

          do j=1, size(rcoords,2)
             do i=1, size(rcoords,1)

                test_loc = is_test_validated(
     $               rcoords(i,j),
     $               rcoords_test(i,j),
     $               .false.)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test failed at: '',2I2)',i,j
                   print '(''  - '',F8.3,''->'',F8.3)', rcoords(i,j), rcoords_test(i,j)
                end if

             end do
          end do

        end function compare_rcoords


        function test_get_central_grdpt(
     $       bf_interface_icr_used,
     $       detailled)
     $       result(test_validated)

          implicit none

          class(bf_interface_icr), intent(in) :: bf_interface_icr_used
          logical                , intent(in) :: detailled
          logical                             :: test_validated


          integer(ikind), dimension(2)  :: d_icoord
          real(rkind)   , dimension(2)  :: d_rcoord
          real(rkind)   , dimension(2)  :: velocity
          real(rkind)   , dimension(nx) :: interior_x_map
          real(rkind)   , dimension(ny) :: interior_y_map
          integer(ikind), dimension(2)  :: d_icoord_n
          real(rkind)   , dimension(2)  :: d_rcoord_n
          integer(ikind), dimension(2)  :: cpt_coord
          integer(ikind), dimension(2)  :: d_icoord_n_test
          real(rkind)   , dimension(2)  :: d_rcoord_n_test
          integer(ikind), dimension(2)  :: cpt_coord_test

          logical :: test_loc

          if(.not.is_test_validated(dt,2.0d0,.false.)) then
             print '(''test_bf_interface_icr'')'
             print '(''test_get_central_grdpt'')'
             print '(''the test requires dt=2.0'')'
             stop ''
          end if
          
          d_icoord            = [2,1]
          d_rcoord            = [0.5d0,0.0d0]
          velocity            = [3.0d0,1.0d0]
          interior_x_map(1:2) = [0.0d0,6.0d0]
          interior_y_map(1:2) = [0.0d0,3.0d0]
          d_icoord_n_test     = [2,1]
          d_rcoord_n_test     = [6.5d0,2.0d0]
          cpt_coord_test      = [6,2]

          cpt_coord = bf_interface_icr_used%get_central_grdpt(
     $         d_icoord,
     $         d_rcoord,
     $         velocity,
     $         interior_x_map,
     $         interior_y_map,
     $         d_icoord_n,
     $         d_rcoord_n)

          test_validated = .true.

          test_loc = d_icoord_n(1).eq.d_icoord_n_test(1)
          if(detailled.and.(.not.test_loc)) then
             print '(''d_icoord_n(1): '',I2,'' -> '',I2)',
     $            d_icoord_n(1), d_icoord_n_test(1)
          end if
          test_validated = test_validated.and.test_loc

          test_loc = d_icoord_n(2).eq.d_icoord_n_test(2)
          if(detailled.and.(.not.test_loc)) then
             print '(''d_icoord_n(2): '',I2,'' -> '',I2)',
     $            d_icoord_n(2), d_icoord_n_test(2)
          end if
          test_validated = test_validated.and.test_loc

          test_loc = is_test_validated(d_rcoord_n(1),d_rcoord_n(1),detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_test_validated(d_rcoord_n(2),d_rcoord_n(2),detailled)
          test_validated = test_validated.and.test_loc          

        end function test_get_central_grdpt


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

      end program test_bf_interface_icr

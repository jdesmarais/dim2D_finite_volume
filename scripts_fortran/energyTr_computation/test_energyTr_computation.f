      program test_energyTr_computation

        use check_data_module, only :
     $       is_real_validated

        use energyTr_computation_module, only :
     $       compute_energyTr_across_x_edge,
     $       compute_energyTr_across_y_edge,
     $       compute_energyTr

        use parameters_kind, only :
     $       rkind


        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled = .true.
        test_validated = .true.

        
        test_loc = test_compute_energyTr_across_x_edge(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_energyTr_across_x_edge: '',L1)', test_loc
        print '()'


        test_loc = test_compute_energyTr_across_y_edge(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_energyTr_across_y_edge: '',L1)', test_loc
        print '()'


        test_loc = test_compute_energyTr(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_energyTr: '',L1)', test_loc
        print '()'

        
        contains


        function test_compute_energyTr_across_x_edge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(8)     :: x_map
          real(rkind), dimension(8,2,4) :: nodes
          integer                       :: i_min
          integer                       :: i_max
          integer                       :: i
          integer                       :: j
          integer                       :: k
          real(rkind)                   :: test_energyTr
          real(rkind)                   :: energyTr
          

          x_map = (/ (-99.0d0,i=1,8) /)
          nodes = reshape( (/ ( ( (-99.0d0,i=1,8) , j=1,2), k=1,4) /),
     $                     (/8,2,4/))

          x_map(2:6) = [ 0.31d0, 0.42d0, 0.62d0, 0.86d0, 0.94d0]
          nodes(3:5,1:2,1) = reshape((/ 0.1d0, 0.52d0, 1.23d0, 8.96d0, 5.12d0, 6.23d0/), (/3,2/))
          nodes(3:5,1:2,3) = reshape((/-1.2d0,  5.6d0, -7.8d0,  1.8d0,  8.2d0, -9.8d0/), (/3,2/))
          nodes(3:5,1:2,4) = reshape((/ 3.6d0,  4.5d0, 4.21d0,  9.6d0, 12.0d0, 7.56d0/), (/3,2/))

          i_min = 3
          i_max = 5
          j     = 1

          test_energyTr = 1.46510634298282d0


          energyTr = compute_EnergyTr_across_x_edge(
     $         x_map,
     $         nodes,
     $         i_min,
     $         i_max,
     $         j)

          test_validated = is_real_validated(
     $         energyTr,
     $         test_energyTr,
     $         detailled)

        end function test_compute_energyTr_across_x_edge


        function test_compute_energyTr_across_y_edge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(8)     :: y_map
          real(rkind), dimension(2,8,4) :: nodes
          integer                       :: j_min
          integer                       :: j_max
          integer                       :: i
          integer                       :: j
          integer                       :: k
          real(rkind)                   :: test_energyTr
          real(rkind)                   :: energyTr
          

          y_map = (/ (-99.0d0,i=1,8) /)
          nodes = reshape( (/ ( ( (-99.0d0,i=1,2) , j=1,8), k=1,4) /),
     $                     (/2,8,4/))

          y_map(2:6) = [0.31d0,0.42d0,0.62d0,0.86d0,0.94d0]
          nodes(1:2,3:5,1) = reshape((/ 0.1d0, 8.96d0, 0.52d0, 5.12d0, 1.23d0, 6.23d0/), (/2,3/))
          nodes(1:2,3:5,2) = reshape((/-1.2d0,  1.8d0,  5.6d0,  8.2d0, -7.8d0, -9.8d0/), (/2,3/))
          nodes(1:2,3:5,4) = reshape((/ 3.6d0,  9.6d0,  4.5d0, 12.0d0, 4.21d0, 7.56d0/), (/2,3/))

          j_min = 3
          j_max = 5
          i     = 1

          test_energyTr = 1.46510634298282d0


          energyTr = compute_EnergyTr_across_y_edge(
     $         y_map,
     $         nodes,
     $         i,
     $         j_min,
     $         j_max)

          test_validated = is_real_validated(
     $         energyTr,
     $         test_energyTr,
     $         detailled)

        end function test_compute_energyTr_across_y_edge



        function test_compute_energyTr(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(8)     :: x_map
          real(rkind), dimension(7)     :: y_map
          real(rkind), dimension(8,7,4) :: nodes
          integer                       :: i
          integer                       :: j
          integer                       :: k
          real(rkind), dimension(4)     :: test_energyTr_tab
          real(rkind)                   :: energyTr
          logical    , dimension(4)     :: borders_tab

          logical :: test_loc

          
          x_map = (/ (-99.0d0,i=1,8) /)
          y_map = (/ (-99.0d0,i=1,7) /)
          nodes = reshape( (/ ( ( (-99.0d0,i=1,8) , j=1,7), k=1,4) /),
     $                     (/8,7,4/))

          x_map(2:7) = [0.31d0, 0.56d0, 0.69d0, 0.75d0, 0.84d0, 0.92d0]
          y_map(2:6) = [1.23d0, 2.53d0, 4.12d0, 6.12d0, 7.52d0]


          nodes(3:6,2:3,1) = reshape((/0.1d0, 0.5d0, 7.8d0, 5.2d0, 5.2d0, 1.23d0, 9.56d0, 4.18d0/), (/4,2/))
          nodes(3:6,5:6,1) = reshape((/6.12d0, 5.12d0, 8.15d0, 9.6d0, 4.52d0, 1.23d0, 9.65d0, 1.23d0/), (/4,2/))
          nodes(2:3,3:5,1) = reshape((/8.25d0, 5.2d0, 6.12d0, 9.63d0, 7.45d0, 6.12d0/), (/2,3/))
          nodes(6:7,3:5,1) = reshape((/4.18d0, 1.45d0, 1.02d0, 7.56d0, 9.6d0, 4.89d0/), (/2,3/))

          nodes(2:3,3:5,2) = reshape((/9.26d0, 1.25d0,-6.12d0, 1.25d0,  6.6d0, 2.36d0/)   , (/2,3/))
          nodes(6:7,3:5,2) = reshape((/5.23d0,-5.69d0,-2.36d0,-6.12d0, 4.12d0, 8.56d0/), (/2,3/))

          nodes(3:6,2:3,3) = reshape((/-5.24d0, 7.62d0, 15.23d0, -5.69d0, -2.36d0, -6.12d0, 5.23d0,  8.56d0/), (/4,2/))
          nodes(3:6,5:6,3) = reshape((/ 2.36d0, 1.25d0, -2.13d0,  4.12d0,  8.26d0,  8.56d0, 1.25d0, 12.03d0/), (/4,2/))

          nodes(3:6,2:3,4) = reshape((/ 0.1d0,  0.5d0,  7.8d0,  5.2d0,   5.2d0, 1.23d0, 9.56d0, 4.18d0/), (/4,2/))
          nodes(3:6,5:6,4) = reshape((/2.12d0, 6.12d0, 8.15d0,  1.6d0, 14.52d0, 3.23d0, 7.65d0, 8.23d0/), (/4,2/))
          nodes(2:3,3:5,4) = reshape((/2.25d0,  3.2d0, 3.12d0, 7.63d0,  3.45d0, 2.12d0/), (/2,3/))
          nodes(6:7,3:5,4) = reshape((/4.18d0, 1.45d0, 1.02d0, 7.56d0,   1.6d0, 4.89d0/), (/2,3/))


          test_energyTr_tab = [
     $         5.44244206184023d0,
     $         6.85939144910197d0,
     $       -11.44960127521820d0,
     $        -1.49525113626922d0]

          
          test_validated = .true.

          do i=1,4
             
             borders_tab = (/ (.false.,i=1,4) /)

             borders_tab(i) = .true.

             energyTr = compute_energyTr(
     $            x_map,
     $            y_map,
     $            nodes,
     $            borders=borders_tab)

             test_loc = is_real_validated(
     $            energyTr,
     $            test_energyTr_tab(i),
     $            detailled)

             test_validated = test_validated.and.test_loc

             if(.not.test_loc .and. detailled) then
                print '(''test('',I2,'') failed'')', i
             end if

          end do

          energyTr = compute_energyTr(
     $         x_map,
     $         y_map,
     $         nodes)

          test_energyTr_tab(1) = -0.64301890054519d0

          test_loc = is_real_validated(
     $         energyTr,
     $         test_energyTr_tab(1),
     $         detailled)

          test_validated = test_validated.and.test_loc

          if(.not.test_loc .and. detailled) then
             print '(''test_all_borders failed'')'
          end if          

        end function test_compute_energyTr

      end program test_energyTr_computation

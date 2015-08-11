      program test_energyTr_computation

        use check_data_module, only :
     $     is_real_validated

        use energyTr_computation_module, only :
     $     compute_energyTr_across_x_edge,
     $     compute_energyTr_across_y_edge

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


c$$$        test_loc = test_compute_energyTr_across_y_edge(detailled)
c$$$        test_validated = test_validated.and.test_loc
c$$$        print '(''test_compute_energyTr_across_y_edge: '',L1)', test_loc
c$$$        print '()'

        
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

          x_map(2:6) = [0.31,0.42,0.62,0.86,0.94]
          nodes(3:5,1:2,1) = reshape((/0.1,0.52,1.23,8.96,5.12,6.23/), (/3,2/))
          nodes(3:5,1:2,3) = reshape((/-1.2,5.6,-7.8,1.8,8.2,-9.8/)  , (/3,2/))
          nodes(3:5,1:2,4) = reshape((/3.6,4.5,4.21,9.6,12.0,7.56/)  , (/3,2/))

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

          y_map(2:6) = [0.31,0.42,0.62,0.86,0.94]
          nodes(1:2,3:5,1) = reshape((/0.1,8.96,0.52,5.12,1.23,6.23/), (/2,3/))
          nodes(1:2,3:5,2) = reshape((/-1.2,1.8,5.6,8.2,-7.8,-9.8/)  , (/2,3/))
          nodes(1:2,3:5,4) = reshape((/3.6,9.6,12.0,4.5,4.21,7.56/)  , (/2,3/))

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

      end program test_energyTr_computation

      program test_bc_operators_wall_x_simplified

        use bc_operators_class, only :
     $     bc_operators
        
        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        test_loc = test_apply_bc_on_nodes(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_nodes: '',L1)', test_loc
        print '()'

        print '(''test_validated: '',L1)', test_validated

        contains


        function test_apply_bc_on_nodes(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(nx)       :: x_map
          real(rkind), dimension(ny)       :: y_map
          real(rkind), dimension(nx,ny,ne) :: nodes
          
          integer :: i
          
          x_map = reshape(
     $         (x_min+(i-1)*),
     $         (/nx/))



        end function test_apply_bc_on_nodes


      end program test_bc_operators_wall_x_simplified

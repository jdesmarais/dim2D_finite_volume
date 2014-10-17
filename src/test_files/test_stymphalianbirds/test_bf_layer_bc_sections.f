      program test_bf_layer_bc_sections

        use bf_layer_bc_procedure_module, only : 
     $     N_edge_type,
     $     S_edge_type,
     $     E_edge_type,
     $     W_edge_type,
     $     NE_corner_type,
     $     NW_corner_type,
     $     SE_corner_type,
     $     SW_corner_type,
     $     NE_edge_type,
     $     NW_edge_type,
     $     SE_edge_type,
     $     SW_edge_type

        use bf_layer_bc_sections_class, only :
     $     bf_layer_bc_sections

        implicit none

        type(bf_layer_bc_sections) :: bc_sections

        integer, dimension(5) :: new_bc_section
        integer               :: k
        
        !test ini()
        print '(''test ini()'')'
        call bc_sections%ini()
        print '(''nb_ele_temp : '',I2)', bc_sections%get_nb_ele_temp()
        print '(''nb_ele_final: '',I2)', bc_sections%get_nb_ele_final()
        print '()'

        !test add_to_temporary_bc_sections() and
        !add_to_final_bc_sections()
        print '(''test add_to_temporary_bc_sections()'')'
        
        do k=1,7
           new_bc_section = get_test_bc_section(k)
           call bc_sections%add_to_temporary_bc_sections(new_bc_section)
           call bc_sections%add_to_final_bc_sections(new_bc_section)
        end do
        call bc_sections%print_bc_sections()


        contains

        function get_test_bc_section(k)
     $       result(bc_section)

          implicit none

          integer, intent(in)   :: k
          integer, dimension(5) :: bc_section

          integer, dimension(12) :: proc_type

          proc_type = [
     $         N_edge_type,
     $         S_edge_type,
     $         E_edge_type,
     $         W_edge_type,
     $         NE_edge_type,
     $         NW_edge_type,
     $         SE_edge_type,
     $         SW_edge_type,
     $         NE_corner_type,
     $         NW_corner_type,
     $         SE_corner_type,
     $         SW_corner_type]

          bc_section = [
     $         proc_type(k),
     $         k,
     $         k,
     $         k+1,
     $         mod(k,3)]

        end function get_test_bc_section

      end program test_bf_layer_bc_sections

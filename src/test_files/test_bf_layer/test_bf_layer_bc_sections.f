      program test_bf_layer_bc_sections

        use bf_layer_bc_sections_class, only :
     $       bf_layer_bc_sections

        use check_data_module, only :
     $       is_int_vector_validated,
     $       is_int_matrix_validated

        use parameters_bf_layer, only : 
     $       align_N, align_S,
     $       align_E, align_W,
     $       
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $     
     $       no_overlap,
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap,
     $       NE_overlap,
     $       NW_overlap,
     $       SE_overlap,
     $       SW_overlap,
     $       NS_overlap,
     $       EW_overlap,
     $       
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt,
     $       no_pt,
     $       
     $       BF_SUCCESS,
     $       
     $       cpt1normal_and_cpt4normal,   
     $       cpt1normal_and_cpt4not,  
     $       cpt1normal_and_cpt4overlap,
     $       cpt1not_and_cpt4normal,
     $       cpt1not_and_cpt4not,
     $       cpt1not_and_cpt4overlap,
     $       cpt1overlap_and_cpt4normal,
     $       cpt1overlap_and_cpt4not,
     $       cpt1overlap_and_cpt4overlap,
     $       
     $       cpt2normal_and_cpt3normal,
     $       cpt2normal_and_cpt3not,  
     $       cpt2normal_and_cpt3overlap,
     $       cpt2not_and_cpt3normal, 
     $       cpt2not_and_cpt3not,     
     $       cpt2not_and_cpt3overlap,
     $       cpt2overlap_and_cpt3normal,
     $       cpt2overlap_and_cpt3not,
     $       cpt2overlap_and_cpt3overlap


        use parameters_input, only :
     $       nx,ny

        use parameters_kind, only :
     $       ikind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled      = .true.
        test_validated = .true.


        test_loc = test_ini(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_ini: '',L1)', test_loc
        print '()'        

        test_loc = test_deallocate_tables(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_deallocate_tables: '',L1)', test_loc
        print '()'

        test_loc = test_add_to_temporary_bc_section(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_add_to_temporary_bc_section: '',L1)', test_loc
        print '()'

        test_loc = test_add_to_final_bc_section(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_add_to_final_bc_section: '',L1)', test_loc
        print '()'

        test_loc = test_add_to_bc_sections(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_add_to_bc_sections: '',L1)', test_loc
        print '()'

        test_loc = test_remove_from_bc_sections_temp(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_remove_from_bc_sections_temp: '',L1)', test_loc
        print '()'

        test_loc = test_remove_from_bc_sections_buffer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_remove_from_bc_sections_buffer: '',L1)', test_loc
        print '()'

        test_loc = test_get_bc_section(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_bc_section: '',L1)', test_loc
        print '()'

        test_loc = test_analyse_grdpt_with_bc_section(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyse_grdpt_with_bc_section: '',L1)', test_loc
        print '()'

        test_loc = test_create_tmp_grdpts_id_for_analyse(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_create_tmp_grdpts_id_for_analyse: '',L1)', test_loc
        print '()'

        test_loc = test_finalize_bc_sections(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_finalize_bc_sections: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains

        
        function test_ini(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_bc_sections) :: bc_sections


          !output
          call bc_sections%ini()

          !validation
          test_validated = bc_sections%get_nb_ele_temp().eq.0
          test_validated = test_validated.and.(bc_sections%get_nb_ele_final().eq.0)

          !detailled
          if(detailled.and.(.not.test_validated)) then
             print '(''nb_ele_temp: '',I2,'' -> '',I2)', bc_sections%get_nb_ele_temp(), 0
             print '(''nb_ele_final: '',I2,'' -> '',I2)', bc_sections%get_nb_ele_final(), 0
          end if

        end function test_ini


        function test_deallocate_tables(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_bc_sections) :: bc_sections

          !output
          call bc_sections%deallocate_tables()

          !validation
          test_validated = .not.allocated(bc_sections%bc_sections_temp)
          test_validated = test_validated.and.(.not.allocated(bc_sections%bc_sections_buffer))
          test_validated = test_validated.and.(.not.allocated(bc_sections%bc_sections_final))

          !detailled
          if(detailled.and.(.not.test_validated)) then
             print '(''bc_sections_temp:   '', L1)', .not.allocated(bc_sections%bc_sections_temp)
             print '(''bc_sections_buffer: '', L1)', .not.allocated(bc_sections%bc_sections_buffer)
             print '(''bc_sections_final:  '', L1)', .not.allocated(bc_sections%bc_sections_final)
          end if

        end function test_deallocate_tables


        function test_add_to_temporary_bc_section(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_bc_sections)           :: bc_sections
          integer                              :: k
          integer, dimension(5,12)             :: test_bc_section
          integer, dimension(5,6)              :: test_bc_sections_temp
          integer, dimension(5,6)              :: test_bc_sections_buffer
          integer, dimension(:,:), allocatable :: bc_sections_temp
          integer, dimension(:,:), allocatable :: bc_sections_buffer
          logical                              :: test_loc


          test_validated = .true.


          !input
          call bc_sections%ini()

          test_bc_section = reshape((/
     $         N_edge_type   ,1,1,2,0, 
     $         S_edge_type   ,2,2,3,1, 
     $         E_edge_type   ,3,3,4,2, 
     $         W_edge_type   ,4,4,5,0, 
     $         NE_edge_type  ,5,5,0,2, 
     $         NW_edge_type  ,6,6,0,0, 
     $         SE_edge_type  ,7,7,0,1, 
     $         SW_edge_type  ,8,8,0,2, 
     $         NE_corner_type,9,9,0,0, 
     $         NW_corner_type,9,9,0,0, 
     $         SE_corner_type,9,9,0,0,  
     $         SW_corner_type,9,9,0,0/),
     $         (/5,12/))
          
          test_bc_sections_temp = reshape(
     $         (/
     $         N_edge_type ,1,1,2,0,
     $         S_edge_type ,2,2,3,1,
     $         E_edge_type ,3,3,4,2,
     $         W_edge_type ,4,4,5,0,
     $         NE_edge_type,5,5,0,2,
     $         NW_edge_type,6,6,0,0
     $         /),
     $         (/5,6/))

          test_bc_sections_buffer = reshape((/
     $         SE_edge_type  ,7,7,0,1,
     $         SW_edge_type  ,8,8,0,2,
     $         NE_corner_type,9,9,0,0,
     $         NW_corner_type,9,9,0,0,
     $         SE_corner_type,9,9,0,0, 
     $         SW_corner_type,9,9,0,0
     $         /),
     $         (/5,6/))

          !output
          do k=1,size(test_bc_section,2)
             call bc_sections%add_to_temporary_bc_sections(test_bc_section(:,k))
          end do

          !validation
          call bc_sections%get_bc_sections_temp(bc_sections_temp)
          call bc_sections%get_bc_sections_buffer(bc_sections_buffer)

          test_loc = is_int_matrix_validated(bc_sections_temp,test_bc_sections_temp,detailled)
          test_loc = test_loc.and.is_int_matrix_validated(bc_sections_buffer,test_bc_sections_buffer,detailled)
          test_validated = test_validated.and.test_loc

          !detailled
          if(detailled.and.(.not.test_loc)) then
             print '(''bc_sections_temp  : '',L1)', is_int_matrix_validated(bc_sections_temp  ,test_bc_sections_temp)
             print '(''bc_sections_buffer: '',L1)', is_int_matrix_validated(bc_sections_buffer,test_bc_sections_buffer)
          end if

          deallocate(bc_sections_temp)
          deallocate(bc_sections_buffer)

        end function test_add_to_temporary_bc_section


        function test_add_to_final_bc_section(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_bc_sections)           :: bc_sections
          integer                              :: k
          integer, dimension(5,12)             :: test_bc_section
          integer, dimension(5,12)             :: test_bc_sections_final
          integer, dimension(:,:), allocatable :: bc_sections_final
          logical                              :: test_loc

          test_validated = .true.


          !input
          call bc_sections%ini()

          test_bc_section = reshape((/
     $         N_edge_type   ,1,1,2,0, 
     $         S_edge_type   ,2,2,3,1, 
     $         E_edge_type   ,3,3,4,2, 
     $         W_edge_type   ,4,4,5,0, 
     $         NE_edge_type  ,5,5,0,2, 
     $         NW_edge_type  ,6,6,0,0, 
     $         SE_edge_type  ,7,7,0,1, 
     $         SW_edge_type  ,8,8,0,2, 
     $         NE_corner_type,9,9,0,0, 
     $         NW_corner_type,9,9,0,0, 
     $         SE_corner_type,9,9,0,0,  
     $         SW_corner_type,9,9,0,0/),
     $         (/5,12/))
          
          test_bc_sections_final = reshape(
     $         (/
     $         N_edge_type   ,1,1,2,0, 
     $         S_edge_type   ,2,2,3,1, 
     $         E_edge_type   ,3,3,4,2, 
     $         W_edge_type   ,4,4,5,0, 
     $         NE_edge_type  ,5,5,0,2, 
     $         NW_edge_type  ,6,6,0,0, 
     $         SE_edge_type  ,7,7,0,1, 
     $         SW_edge_type  ,8,8,0,2, 
     $         NE_corner_type,9,9,0,0, 
     $         NW_corner_type,9,9,0,0, 
     $         SE_corner_type,9,9,0,0,  
     $         SW_corner_type,9,9,0,0/),
     $         (/5,12/))


          !output
          do k=1,size(test_bc_section,2)
             call bc_sections%add_to_final_bc_sections(test_bc_section(:,k))
          end do

          !validation
          call bc_sections%get_bc_sections_final(bc_sections_final)
          test_loc = is_int_matrix_validated(bc_sections_final,test_bc_sections_final)
          test_validated = test_validated.and.test_loc

          !detailled
          if(detailled.and.(.not.test_loc)) then
             print '(''bc_sections_final: '',L1)', is_int_matrix_validated(bc_sections_final,test_bc_sections_final)
          end if

          deallocate(bc_sections_final)

        end function test_add_to_final_bc_section


        function test_add_to_bc_sections(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_bc_sections)           :: bc_sections
          integer                              :: k
          integer, dimension(5,12)             :: test_bc_section
          integer, dimension(5,6)              :: test_bc_sections_temp
          integer, dimension(5,2)              :: test_bc_sections_buffer
          integer, dimension(5,4)              :: test_bc_sections_final
          integer, dimension(:,:), allocatable :: bc_sections_temp
          integer, dimension(:,:), allocatable :: bc_sections_buffer
          integer, dimension(:,:), allocatable :: bc_sections_final
          logical                              :: test_loc

          test_validated = .true.


          !input
          call bc_sections%ini()

          test_bc_section = reshape((/
     $         N_edge_type   ,1,1,2,0, 
     $         S_edge_type   ,2,2,3,1, 
     $         E_edge_type   ,3,3,4,2, 
     $         W_edge_type   ,4,4,5,0, 
     $         NE_edge_type  ,5,5,0,2, 
     $         NW_edge_type  ,6,6,0,0, 
     $         SE_edge_type  ,7,7,0,1, 
     $         SW_edge_type  ,8,8,0,2, 
     $         NE_corner_type,9,9,0,0, 
     $         NW_corner_type,9,9,0,0, 
     $         SE_corner_type,9,9,0,0,  
     $         SW_corner_type,9,9,0,0/),
     $         (/5,12/))

          test_bc_sections_temp = reshape((/
     $         N_edge_type   ,1,1,2,0, 
     $         S_edge_type   ,2,2,3,1, 
     $         E_edge_type   ,3,3,4,2, 
     $         W_edge_type   ,4,4,5,0, 
     $         NE_edge_type  ,5,5,0,2, 
     $         NW_edge_type  ,6,6,0,0/),
     $         (/5,6/))

          test_bc_sections_buffer = reshape((/
     $         SE_edge_type  ,7,7,0,1, 
     $         SW_edge_type  ,8,8,0,2/),
     $         (/5,2/))

          test_bc_sections_final = reshape((/
     $         NE_corner_type,9,9,0,0, 
     $         NW_corner_type,9,9,0,0, 
     $         SE_corner_type,9,9,0,0,  
     $         SW_corner_type,9,9,0,0/),
     $         (/5,4/))

          !output
          do k=1,size(test_bc_section,2)
             call bc_sections%add_to_bc_sections(test_bc_section(:,k))
          end do

          !validation
          call bc_sections%get_bc_sections_temp(bc_sections_temp)
          test_loc = is_int_matrix_validated(bc_sections_temp,test_bc_sections_temp)
          test_validated = test_validated.and.test_loc

          call bc_sections%get_bc_sections_buffer(bc_sections_buffer)
          test_loc = is_int_matrix_validated(bc_sections_buffer,test_bc_sections_buffer)
          test_validated = test_validated.and.test_loc

          call bc_sections%get_bc_sections_final(bc_sections_final)
          test_loc = is_int_matrix_validated(bc_sections_final(:,1:4),test_bc_sections_final)
          test_validated = test_validated.and.test_loc

          !detailled
          if(detailled.and.(.not.test_loc)) then
             print '(''bc_sections_temp: '',L1)', is_int_matrix_validated(bc_sections_temp,test_bc_sections_temp)
             print '(''bc_sections_buffer: '',L1)', is_int_matrix_validated(bc_sections_buffer,test_bc_sections_buffer)
             print '(''bc_sections_final: '',L1)', is_int_matrix_validated(bc_sections_final(:,1:4),test_bc_sections_final)
          end if

          deallocate(bc_sections_final)

        end function test_add_to_bc_sections


        function test_remove_from_bc_sections_temp(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_bc_sections)           :: bc_sections
          integer                              :: k
          integer, dimension(5,8)              :: test_bc_section
          integer, dimension(5,6)              :: test_bc_sections_temp
          integer, dimension(5,1)              :: test_bc_sections_buffer
          integer, dimension(:,:), allocatable :: bc_sections_temp
          integer, dimension(:,:), allocatable :: bc_sections_buffer
          logical                              :: test_loc

          test_validated = .true.


          !input
          call bc_sections%ini()

          test_bc_section = reshape((/
     $         N_edge_type   ,1,1,2,0, 
     $         S_edge_type   ,2,2,3,1, 
     $         E_edge_type   ,3,3,4,2, 
     $         W_edge_type   ,4,4,5,0, 
     $         NE_edge_type  ,5,5,0,2, 
     $         NW_edge_type  ,6,6,0,0, 
     $         SE_edge_type  ,7,7,0,1, 
     $         SW_edge_type  ,8,8,0,2/),
     $         (/5,8/))

          test_bc_sections_temp = reshape((/
     $         N_edge_type   ,1,1,2,0,
     $         E_edge_type   ,3,3,4,2, 
     $         W_edge_type   ,4,4,5,0, 
     $         NE_edge_type  ,5,5,0,2, 
     $         NW_edge_type  ,6,6,0,0,
     $         SE_edge_type  ,7,7,0,1/),
     $         (/5,6/))

          test_bc_sections_buffer = reshape((/
     $         SW_edge_type  ,8,8,0,2/),
     $         (/5,1/))

          !output
          do k=1,size(test_bc_section,2)
             call bc_sections%add_to_bc_sections(test_bc_section(:,k))
          end do
          call bc_sections%remove_from_bc_sections_temp(2)

          !validation
          call bc_sections%get_bc_sections_temp(bc_sections_temp)
          test_loc = is_int_matrix_validated(bc_sections_temp,test_bc_sections_temp)
          test_validated = test_validated.and.test_loc

          call bc_sections%get_bc_sections_buffer(bc_sections_buffer)
          test_loc = is_int_matrix_validated(bc_sections_buffer(:,1:1),test_bc_sections_buffer)
          test_validated = test_validated.and.test_loc

          !detailled
          if(detailled.and.(.not.test_loc)) then
             print '(''bc_sections_temp: '',L1)', is_int_matrix_validated(bc_sections_temp,test_bc_sections_temp)
             print '(''bc_sections_buffer: '',L1)', is_int_matrix_validated(bc_sections_buffer(:,1:1),test_bc_sections_buffer)
          end if

          deallocate(bc_sections_temp)
          deallocate(bc_sections_buffer)

        end function test_remove_from_bc_sections_temp


        function test_remove_from_bc_sections_buffer(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_bc_sections)           :: bc_sections
          integer                              :: k
          integer, dimension(5,12)             :: test_bc_section
          integer, dimension(5,6)              :: test_bc_sections_temp
          integer, dimension(5,5)              :: test_bc_sections_buffer
          integer, dimension(:,:), allocatable :: bc_sections_temp
          integer, dimension(:,:), allocatable :: bc_sections_buffer
          logical                              :: test_loc

          test_validated = .true.


          !input
          call bc_sections%ini()

          test_bc_section = reshape((/
     $         N_edge_type   ,1,1,2,0, 
     $         S_edge_type   ,2,2,3,1, 
     $         E_edge_type   ,3,3,4,2, 
     $         W_edge_type   ,4,4,5,0, 
     $         NE_edge_type  ,5,5,0,2, 
     $         NW_edge_type  ,6,6,0,0, 
     $         SE_edge_type  ,7,7,0,1, 
     $         SW_edge_type  ,8,8,0,2,
     $         NE_corner_type,9,9,0,0, 
     $         NW_corner_type,9,9,0,0, 
     $         SE_corner_type,9,9,0,0,  
     $         SW_corner_type,9,9,0,0/),
     $         (/5,12/))

          test_bc_sections_temp = reshape((/
     $         N_edge_type   ,1,1,2,0,
     $         S_edge_type   ,2,2,3,1, 
     $         E_edge_type   ,3,3,4,2, 
     $         W_edge_type   ,4,4,5,0, 
     $         NE_edge_type  ,5,5,0,2, 
     $         NW_edge_type  ,6,6,0,0/),
     $         (/5,6/))

          test_bc_sections_buffer = reshape((/
     $         SE_edge_type  ,7,7,0,1,
     $         SW_edge_type  ,8,8,0,2,
     $         NE_corner_type,9,9,0,0, 
     $         NW_corner_type,9,9,0,0, 
     $         SW_corner_type,9,9,0,0/),
     $         (/5,5/))

          !output
          do k=1,size(test_bc_section,2)
             call bc_sections%add_to_temporary_bc_sections(test_bc_section(:,k))
          end do
          call bc_sections%remove_from_bc_sections_buffer(11)

          !validation
          call bc_sections%get_bc_sections_temp(bc_sections_temp)
          test_loc = is_int_matrix_validated(bc_sections_temp,test_bc_sections_temp)
          test_validated = test_validated.and.test_loc

          call bc_sections%get_bc_sections_buffer(bc_sections_buffer)
          test_loc = is_int_matrix_validated(bc_sections_buffer(:,1:5),test_bc_sections_buffer)
          test_validated = test_validated.and.test_loc

          !detailled
          if(detailled.and.(.not.test_loc)) then
             print '(''bc_sections_temp: '',L1)', is_int_matrix_validated(bc_sections_temp,test_bc_sections_temp)
             print '(''bc_sections_buffer: '',L1)', is_int_matrix_validated(bc_sections_buffer(:,1:5),test_bc_sections_buffer)
          end if

          deallocate(bc_sections_temp)
          deallocate(bc_sections_buffer)

        end function test_remove_from_bc_sections_buffer


c$$$        subroutine test_add_bc_sections()
c$$$
c$$$          implicit none
c$$$
c$$$          type(bf_layer_bc_sections) :: bc_sections
c$$$          integer                    :: k
c$$$          integer, dimension(5)      :: new_bc_section
c$$$
c$$$          integer, dimension(5,6)              :: test_bc_sections_temp
c$$$          integer, dimension(5,2)              :: test_bc_sections_buffer
c$$$          integer, dimension(5,1)              :: test_bc_sections_final
c$$$
c$$$          logical :: same
c$$$
c$$$          print '(''test add_to_bc_sections()'')'
c$$$          call bc_sections%ini()
c$$$          do k=1,9
c$$$             new_bc_section = get_test_bc_section(k)
c$$$             call bc_sections%add_to_bc_sections(new_bc_section)
c$$$          end do
c$$$
c$$$          !bc_sections_temp
c$$$          test_bc_sections_temp = reshape(
c$$$     $         (/
c$$$     $         N_edge_type,1,1,2,0,
c$$$     $         S_edge_type,2,2,3,1,
c$$$     $         E_edge_type,3,3,4,2,
c$$$     $         W_edge_type,4,4,5,0,
c$$$     $         NE_edge_type,5,5,0,2,
c$$$     $         NW_edge_type,6,6,0,0
c$$$     $         /),
c$$$     $         (/5,6/))
c$$$
c$$$          !bc_sections_buffer
c$$$          test_bc_sections_buffer = reshape(
c$$$     $         (/
c$$$     $         SE_edge_type,7,7,0,1,
c$$$     $         SW_edge_type,8,8,0,2
c$$$     $         /),
c$$$     $         (/5,2/))
c$$$          
c$$$          !bc_sections_final
c$$$          test_bc_sections_final = reshape(
c$$$     $         (/
c$$$     $         NE_corner_type,9,9,0,0
c$$$     $         /),
c$$$     $         (/5,1/))
c$$$
c$$$
c$$$          !same bc_sections_temp and bc_sections_buffer ?
c$$$          same = compare_bc_sections_temporary(
c$$$     $         bc_sections,
c$$$     $         test_bc_sections_temp,
c$$$     $         test_bc_sections_buffer)
c$$$          
c$$$          print '(''bc_sections_temporary: '',L1)', same
c$$$
c$$$
c$$$          !same bc_sections_final ?
c$$$          same = compare_bc_sections_final(
c$$$     $         bc_sections,
c$$$     $         test_bc_sections_final)          
c$$$
c$$$          print '(''bc_sections_final    : '',L1)', same
c$$$          print '()'
c$$$          
c$$$        end subroutine test_add_bc_sections
c$$$
c$$$
c$$$        function compare_bc_sections_final(
c$$$     $     bc_sections,
c$$$     $     test_bc_sections_final)
c$$$     $     result(same)
c$$$
c$$$          implicit none
c$$$
c$$$          type(bf_layer_bc_sections), intent(in) :: bc_sections
c$$$          integer, dimension(:,:)   , intent(in) :: test_bc_sections_final
c$$$          logical                                :: same
c$$$
c$$$          integer                                :: nb_ele_final
c$$$          integer, dimension(:,:), allocatable   :: bc_sections_final
c$$$          integer                                :: k
c$$$
c$$$
c$$$          same = .true.
c$$$
c$$$          nb_ele_final = bc_sections%get_nb_ele_final()
c$$$          call bc_sections%get_bc_sections_final(bc_sections_final)
c$$$
c$$$          if(allocated(bc_sections_final)) then
c$$$
c$$$             do k=1, min(size(bc_sections_final,2),nb_ele_final)
c$$$
c$$$                same = compare_bc_procedure(
c$$$     $               bc_sections_final(:,k),
c$$$     $               test_bc_sections_final(:,k))
c$$$
c$$$                if(.not.same) then
c$$$                   print '(''   k                    : '', I3)' , k
c$$$                   print '(''   bc_section_final     : '', 5I3)', bc_sections_final(:,k)
c$$$                   print '(''   test_bc_section_final: '', 5I3)', test_bc_sections_final(:,k)
c$$$                   exit
c$$$                end if
c$$$
c$$$             end do
c$$$
c$$$          end if
c$$$
c$$$        end function compare_bc_sections_final
c$$$
c$$$      
c$$$        function compare_bc_sections_sorted(
c$$$     $     bc_sections_sorted,
c$$$     $     test_bc_sections_sorted)
c$$$     $     result(same)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:)   , intent(in) :: bc_sections_sorted
c$$$          integer, dimension(:,:)   , intent(in) :: test_bc_sections_sorted
c$$$          logical                                :: same
c$$$
c$$$          integer                                :: k
c$$$
c$$$
c$$$          same = .true.
c$$$
c$$$          do k=1, size(bc_sections_sorted,2)
c$$$
c$$$             same = compare_bc_layer_sorted(
c$$$     $            bc_sections_sorted(:,k),
c$$$     $            test_bc_sections_sorted(:,k))
c$$$
c$$$             if(.not.same) then
c$$$                print '(''   k                     : '',  I3)', k
c$$$                print '(''   bc_section_sorted     : '', 4I3)', bc_sections_sorted(:,k)
c$$$                print '(''   test_bc_section_sorted: '', 4I3)', test_bc_sections_sorted(:,k)
c$$$                exit
c$$$             end if
c$$$             
c$$$          end do
c$$$
c$$$        end function compare_bc_sections_sorted
c$$$
c$$$
c$$$        function compare_bc_layer_sorted(bc_sorted,test_bc_sorted)
c$$$     $     result(same)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(4), intent(in) :: bc_sorted
c$$$          integer, dimension(4), intent(in) :: test_bc_sorted
c$$$          logical                           :: same
c$$$
c$$$          same = bc_sorted(1).eq.test_bc_sorted(1)          
c$$$
c$$$          select case(bc_sorted(1))
c$$$
c$$$            case(N_edge_type,S_edge_type,E_edge_type,W_edge_type)
c$$$               same = same.and.(bc_sorted(2).eq.test_bc_sorted(2))
c$$$               same = same.and.(bc_sorted(3).eq.test_bc_sorted(3))
c$$$               same = same.and.(bc_sorted(4).eq.test_bc_sorted(4))
c$$$
c$$$            case(
c$$$     $              NE_corner_type,NW_corner_type,
c$$$     $              SE_corner_type,SW_corner_type,
c$$$     $              NE_edge_type,NW_edge_type,
c$$$     $              SE_edge_type,SW_edge_type
c$$$     $              )
c$$$               same = same.and.(bc_sorted(2).eq.test_bc_sorted(2))
c$$$               same = same.and.(bc_sorted(3).eq.test_bc_sorted(3))
c$$$
c$$$            case default
c$$$               print '(''test_bf_layer_bc_sections'')'
c$$$               print '(''compare_bc_layer_sorted'')'
c$$$               print '(''case not recognized: '',I2)', bc_sorted(1)
c$$$               stop ''
c$$$
c$$$          end select
c$$$
c$$$        end function compare_bc_layer_sorted
c$$$
c$$$
c$$$        subroutine test_remove_bc_section()
c$$$
c$$$          implicit none
c$$$
c$$$          type(bf_layer_bc_sections)           :: bc_sections
c$$$          integer                              :: k
c$$$          integer, dimension(5)                :: new_bc_section
c$$$          integer, dimension(5,6)              :: test_bc_sections_temp
c$$$          integer, dimension(:,:), allocatable :: test_bc_sections_buffer
c$$$          logical                              :: same
c$$$
c$$$
c$$$          !initialization
c$$$          call bc_sections%ini()
c$$$          
c$$$          do k=1,9
c$$$             new_bc_section = get_test_bc_section(k)
c$$$             call bc_sections%add_to_temporary_bc_sections(new_bc_section)
c$$$          end do
c$$$
c$$$
c$$$          !test remove_from_bc_sections_temp()
c$$$          print '(''test remove_from_bc_sections_temp()'')'
c$$$          call bc_sections%remove_from_bc_sections_temp(1)
c$$$          test_bc_sections_temp = reshape(
c$$$     $         (/
c$$$     $         S_edge_type,2,2,3,0,
c$$$     $         E_edge_type,3,3,4,0,
c$$$     $         W_edge_type,4,4,5,0,
c$$$     $         NE_edge_type,5,5,0,2,
c$$$     $         NW_edge_type,6,6,0,0,
c$$$     $         SE_edge_type,7,7,0,1
c$$$     $         /),
c$$$     $         (/5,6/))
c$$$          allocate(test_bc_sections_buffer(5,2))
c$$$          test_bc_sections_buffer = reshape(
c$$$     $         (/
c$$$     $         SW_edge_type,8,8,0,2,
c$$$     $         NE_corner_type,9,9,0,0
c$$$     $         /),
c$$$     $         (/5,2/))
c$$$
c$$$          same = compare_bc_sections_temporary(
c$$$     $         bc_sections,
c$$$     $         test_bc_sections_temp,
c$$$     $         test_bc_sections_buffer)
c$$$
c$$$          deallocate(test_bc_sections_buffer)
c$$$
c$$$          print '(''remove_from_bc_sections_temp(1): '',L1)', same
c$$$
c$$$
c$$$          call bc_sections%remove_from_bc_sections_temp(4)
c$$$          test_bc_sections_temp = reshape(
c$$$     $         (/
c$$$     $         S_edge_type,2,2,3,0,
c$$$     $         E_edge_type,3,3,4,0,
c$$$     $         W_edge_type,4,4,5,0,
c$$$     $         NW_edge_type,6,6,0,0,
c$$$     $         SE_edge_type,7,7,0,1,
c$$$     $         SW_edge_type,8,8,0,2
c$$$     $         /),
c$$$     $         (/5,6/))
c$$$          allocate(test_bc_sections_buffer(5,1))
c$$$          test_bc_sections_buffer = reshape(
c$$$     $         (/
c$$$     $         NE_corner_type,9,9,0,0
c$$$     $         /),
c$$$     $         (/5,1/))
c$$$
c$$$          same = compare_bc_sections_temporary(
c$$$     $         bc_sections,
c$$$     $         test_bc_sections_temp,
c$$$     $         test_bc_sections_buffer)
c$$$
c$$$          deallocate(test_bc_sections_buffer)
c$$$
c$$$          print '(''remove_from_bc_sections_temp(4): '',L1)', same
c$$$
c$$$
c$$$          call bc_sections%remove_from_bc_sections_temp(6)
c$$$          test_bc_sections_temp = reshape(
c$$$     $         (/
c$$$     $         S_edge_type,2,2,3,0,
c$$$     $         E_edge_type,3,3,4,0,
c$$$     $         W_edge_type,4,4,5,0,
c$$$     $         NW_edge_type,6,6,0,0,
c$$$     $         SE_edge_type,7,7,0,1,
c$$$     $         NE_corner_type,9,9,0,0
c$$$     $         /),
c$$$     $         (/5,6/))
c$$$
c$$$          same = compare_bc_sections_temporary(
c$$$     $         bc_sections,
c$$$     $         test_bc_sections_temp,
c$$$     $         test_bc_sections_buffer)
c$$$
c$$$          print '(''remove_from_bc_sections_temp(6): '',L1)', same
c$$$
c$$$          print '()'
c$$$
c$$$
c$$$          do k=10,12
c$$$             new_bc_section = get_test_bc_section(k)
c$$$             call bc_sections%add_to_temporary_bc_sections(new_bc_section)
c$$$          end do
c$$$
c$$$          call bc_sections%remove_from_bc_sections_buffer(1)
c$$$          
c$$$          allocate(test_bc_sections_buffer(5,2))
c$$$
c$$$          test_bc_sections_buffer = reshape(
c$$$     $         (/
c$$$     $         SE_corner_type,11,11,0,2,
c$$$     $         SW_corner_type,12,12,0,0
c$$$     $         /),
c$$$     $         (/5,2/))
c$$$
c$$$          same = compare_bc_sections_temporary(
c$$$     $         bc_sections,
c$$$     $         test_bc_sections_temp,
c$$$     $         test_bc_sections_buffer)
c$$$
c$$$          deallocate(test_bc_sections_buffer)
c$$$
c$$$          print '(''remove_from_bc_sections_buffer(1): '',L1)', same
c$$$
c$$$
c$$$          call bc_sections%remove_from_bc_sections_buffer(2)
c$$$
c$$$          allocate(test_bc_sections_buffer(5,1))
c$$$
c$$$          test_bc_sections_buffer = reshape(
c$$$     $         (/
c$$$     $          SE_corner_type,11,11,0,2
c$$$     $         /),
c$$$     $         (/5,1/))
c$$$
c$$$          same = compare_bc_sections_temporary(
c$$$     $         bc_sections,
c$$$     $         test_bc_sections_temp,
c$$$     $         test_bc_sections_buffer)
c$$$
c$$$          deallocate(test_bc_sections_buffer)
c$$$
c$$$          print '(''remove_from_bc_sections_temp(2): '',L1)', same
c$$$
c$$$          print '()'
c$$$
c$$$      end subroutine test_remove_bc_section
c$$$
c$$$
c$$$      function compare_bc_sections_temporary(
c$$$     $     bc_sections,
c$$$     $     test_bc_sections_temp,
c$$$     $     test_bc_sections_buffer)
c$$$     $     result(same)
c$$$
c$$$        implicit none
c$$$
c$$$        class(bf_layer_bc_sections), intent(in) :: bc_sections
c$$$        integer, dimension(:,:)    , intent(in) :: test_bc_sections_temp
c$$$        integer, dimension(:,:)    , intent(in) :: test_bc_sections_buffer
c$$$        logical                                 :: same
c$$$
c$$$
c$$$        integer :: k
c$$$        integer :: nb_ele_temp
c$$$        integer, dimension(:,:), allocatable :: bc_sections_temp
c$$$        integer, dimension(:,:), allocatable :: bc_sections_buffer
c$$$
c$$$
c$$$        nb_ele_temp = bc_sections%get_nb_ele_temp()
c$$$        call bc_sections%get_bc_sections_temp(bc_sections_temp)
c$$$        call bc_sections%get_bc_sections_buffer(bc_sections_buffer)
c$$$
c$$$
c$$$        same = .true.
c$$$        if(allocated(bc_sections_temp)) then
c$$$
c$$$           do k=1, min(size(bc_sections_temp,2),nb_ele_temp)
c$$$
c$$$              same = compare_bc_procedure(
c$$$     $             bc_sections_temp(:,k),
c$$$     $             test_bc_sections_temp(:,k))
c$$$
c$$$              if(.not.same) then
c$$$                 print '(''   k                   : '', I3)' , k
c$$$                 print '(''   bc_section_temp     : '', 5I3)', bc_sections_temp(:,k)
c$$$                 print '(''   test_bc_section_temp: '', 5I3)', test_bc_sections_temp(:,k)
c$$$                 exit
c$$$              end if
c$$$
c$$$           end do
c$$$
c$$$           deallocate(bc_sections_temp)
c$$$
c$$$        end if
c$$$
c$$$        if(same) then
c$$$           if(allocated(bc_sections_buffer)) then
c$$$              
c$$$              do k=1, min(
c$$$     $             size(bc_sections_buffer,2),
c$$$     $             nb_ele_temp-size(bc_sections_temp,2))
c$$$                 
c$$$                 same = compare_bc_procedure(
c$$$     $                bc_sections_buffer(:,k),
c$$$     $                test_bc_sections_buffer(:,k))
c$$$
c$$$                 if(.not.same) then
c$$$                    print '(''   k                   : '', I3)' , k
c$$$                    print '(''   bc_section_buff     : '', 5I3)', bc_sections_buffer(:,k)
c$$$                    print '(''   test_bc_section_buff: '', 5I3)', test_bc_sections_buffer(:,k)
c$$$                    exit
c$$$                 end if
c$$$                 
c$$$              end do
c$$$
c$$$              deallocate(bc_sections_buffer)
c$$$           
c$$$           end if
c$$$        end if
c$$$
c$$$      end function compare_bc_sections_temporary
c$$$
c$$$
c$$$      function compare_bc_procedure(bc_sections,test_bc_sections)
c$$$     $     result(same)
c$$$
c$$$        implicit none
c$$$
c$$$        integer, dimension(:), intent(in) :: bc_sections
c$$$        integer, dimension(:), intent(in) :: test_bc_sections
c$$$        logical                           :: same
c$$$
c$$$        same = .true.
c$$$
c$$$        same = bc_sections(1).eq.test_bc_sections(1)
c$$$
c$$$        select case(bc_sections(1))
c$$$          case(N_edge_type,S_edge_type,E_edge_type,W_edge_type)
c$$$
c$$$             same = same.and.(bc_sections(2).eq.test_bc_sections(2))
c$$$             same = same.and.(bc_sections(3).eq.test_bc_sections(3))
c$$$             same = same.and.(bc_sections(4).eq.test_bc_sections(4))
c$$$
c$$$          case(NE_corner_type,NW_corner_type,SE_corner_type,SW_corner_type)
c$$$
c$$$             same = same.and.(bc_sections(2).eq.test_bc_sections(2))
c$$$             same = same.and.(bc_sections(3).eq.test_bc_sections(3))
c$$$
c$$$          case(NE_edge_type,NW_edge_type,SE_edge_type,SW_edge_type)
c$$$
c$$$             same = same.and.(bc_sections(2).eq.test_bc_sections(2))
c$$$             same = same.and.(bc_sections(3).eq.test_bc_sections(3))
c$$$             same = same.and.(bc_sections(5).eq.test_bc_sections(5))
c$$$
c$$$          case default
c$$$             print '(''test_bf_layer_bc_sections'')'
c$$$             print '(''compare_bc_procedure'')'
c$$$             print '(''case: '',I2)', bc_sections(1)
c$$$             stop 'case not recognized'
c$$$        end select
c$$$
c$$$      end function compare_bc_procedure


      function test_get_bc_section(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_bc_sections) :: bf_layer_bc_sections_used
          integer                    :: test_i
          integer                    :: test_j
          integer, dimension(3,3,12) :: test_grdpts_id
          integer, dimension(5,12)   :: test_bc_section
          integer, dimension(5)      :: bc_section
          integer                    :: k
          logical                    :: test_loc
          logical                    :: ierror

          integer, dimension(3,3) :: tmp_grdpts_id
          logical                 :: use_tmp_grdpts_id


          test_validated = .true.


          !input
          test_i = 2
          test_j = 2

          !N_edge
          test_grdpts_id(:,:,1) = reshape((/
     $         interior_pt,interior_pt,interior_pt,
     $         bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $         bc_pt,bc_pt,bc_pt/),
     $         (/3,3/))

          test_bc_section(:,1) = [N_edge_type,2,2,2,1]

          !S_edge
          test_grdpts_id(:,:,2) = reshape((/
     $         bc_pt,bc_pt,bc_pt,
     $         bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $         interior_pt,interior_pt,interior_pt/),
     $         (/3,3/))

          test_bc_section(:,2) = [S_edge_type,2,2,2,1]

          !E_edge
          test_grdpts_id(:,:,3) = reshape((/
     $         interior_pt,bc_interior_pt,bc_pt,
     $         interior_pt,bc_interior_pt,bc_pt,
     $         interior_pt,bc_interior_pt,bc_pt/),
     $         (/3,3/))

          test_bc_section(:,3) = [E_edge_type,2,2,2,1]

          !W_edge
          test_grdpts_id(:,:,4) = reshape((/
     $         bc_pt,bc_interior_pt,interior_pt,
     $         bc_pt,bc_interior_pt,interior_pt,
     $         bc_pt,bc_interior_pt,interior_pt/),
     $         (/3,3/))

          test_bc_section(:,4) = [W_edge_type,2,2,2,1]


          !NE_corner
          test_grdpts_id(:,:,5) = reshape((/
     $         interior_pt,bc_interior_pt,bc_pt,
     $         bc_interior_pt,bc_interior_pt,bc_pt,
     $         bc_pt,bc_pt,bc_pt/),
     $         (/3,3/))

          test_bc_section(:,5) = [NE_corner_type,2,2,2,1]

          !NW_corner
          test_grdpts_id(:,:,6) = reshape((/
     $         bc_pt,bc_interior_pt,interior_pt,
     $         bc_pt,bc_interior_pt,bc_interior_pt,
     $         bc_pt,bc_pt,bc_pt/),
     $         (/3,3/))

          test_bc_section(:,6) = [NW_corner_type,1,2,2,1]

          !SE_corner
          test_grdpts_id(:,:,7) = reshape((/
     $         bc_pt,bc_pt,bc_pt,
     $         bc_interior_pt,bc_interior_pt,bc_pt,
     $         interior_pt,bc_interior_pt,bc_pt/),
     $         (/3,3/))

          test_bc_section(:,7) = [SE_corner_type,2,1,2,1]

          !SW_corner
          test_grdpts_id(:,:,8) = reshape((/
     $         bc_pt,bc_pt,bc_pt,
     $         bc_pt,bc_interior_pt,bc_interior_pt,
     $         bc_pt,bc_interior_pt,interior_pt/),
     $         (/3,3/))

          test_bc_section(:,8) = [SW_corner_type,1,1,2,1]


          !NE_edge
          test_grdpts_id(:,:,9) = reshape((/
     $         interior_pt,interior_pt,interior_pt,
     $         interior_pt,bc_interior_pt,bc_interior_pt,
     $         interior_pt,bc_interior_pt,bc_pt/),
     $         (/3,3/))

          test_bc_section(:,9) = [NE_edge_type,2,2,2,1]

          !NW_edge
          test_grdpts_id(:,:,10) = reshape((/
     $         interior_pt,interior_pt,interior_pt,
     $         bc_interior_pt,bc_interior_pt,interior_pt,
     $         bc_pt,bc_interior_pt,interior_pt/),
     $         (/3,3/))

          test_bc_section(:,10) = [NW_edge_type,1,2,2,1]

          !SE_edge
          test_grdpts_id(:,:,11) = reshape((/
     $         interior_pt,bc_interior_pt,bc_pt,
     $         interior_pt,bc_interior_pt,bc_interior_pt,
     $         interior_pt,interior_pt,interior_pt/),
     $         (/3,3/))

          test_bc_section(:,11) = [SE_edge_type,2,1,2,1]

          !SW_edge
          test_grdpts_id(:,:,12) = reshape((/
     $         bc_pt,bc_interior_pt,interior_pt,
     $         bc_interior_pt,bc_interior_pt,interior_pt,
     $         interior_pt,interior_pt,interior_pt/),
     $         (/3,3/))

          test_bc_section(:,12) = [SW_edge_type,1,1,2,1]
          
          
          !test the edge detection
          !============================================================
          ! test w/o the tmp_grdpts_id
          !------------------------------------------------------------
          use_tmp_grdpts_id = .false.
          do k=1,4

             !output
             bc_section = bf_layer_bc_sections_used%get_bc_section(
     $            test_i,
     $            test_j,
     $            test_grdpts_id(:,:,k),
     $            tmp_grdpts_id, use_tmp_grdpts_id,
     $            ierror)

             !validation
             if(ierror.eqv.BF_SUCCESS) then
                test_loc = is_int_vector_validated(
     $               bc_section(1:4),
     $               test_bc_section(1:4,k),
     $               detailled)                
             else
                test_loc = .false.
             end if

             test_validated = test_validated.and.test_loc
             
             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test ('',I2,'',1) failed'')', k
             end if

          end do

          ! test w/ the tmp_grdpts_id
          !------------------------------------------------------------
          use_tmp_grdpts_id = .true.
          do k=1,4

             !input
             tmp_grdpts_id = test_grdpts_id(:,:,k)

             !output
             bc_section = bf_layer_bc_sections_used%get_bc_section(
     $            test_i,
     $            test_j,
     $            test_grdpts_id(:,:,1),
     $            tmp_grdpts_id, use_tmp_grdpts_id,
     $            ierror)

             !validation
             if(ierror.eqv.BF_SUCCESS) then
                test_loc = is_int_vector_validated(
     $               bc_section(1:4),
     $               test_bc_section(1:4,k),
     $               detailled)                
             else
                test_loc = .false.
             end if

             test_validated = test_validated.and.test_loc
             
             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test ('',I2,'',2) failed'')', k
             end if

          end do


          ! test the corner and anti-corner detection
          !============================================================
          !test w/o tmp_grdpts_id
          !------------------------------------------------------------
          use_tmp_grdpts_id = .false.
          do k=5,size(test_grdpts_id,3)

             !output
             bc_section = bf_layer_bc_sections_used%get_bc_section(
     $            test_i,
     $            test_j,
     $            test_grdpts_id(:,:,k),
     $            tmp_grdpts_id, use_tmp_grdpts_id,
     $            ierror)

             !validation
             if(ierror.eqv.BF_SUCCESS) then
                test_loc = is_int_vector_validated(
     $               bc_section(1:3),
     $               test_bc_section(1:3,k),
     $               detailled)
             else
                test_loc = .false.
             end if

             test_validated = test_validated.and.test_loc
             
             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test ('',I2,'',1) failed'')', k
             end if

          end do

          !test w/ tmp_grdpts_id
          !------------------------------------------------------------
          use_tmp_grdpts_id = .true.
          do k=5,size(test_grdpts_id,3)

             !input
             tmp_grdpts_id = test_grdpts_id(:,:,k)

             !output
             bc_section = bf_layer_bc_sections_used%get_bc_section(
     $            test_i,
     $            test_j,
     $            test_grdpts_id(:,:,1),
     $            tmp_grdpts_id, use_tmp_grdpts_id,
     $            ierror)

             !validation
             if(ierror.eqv.BF_SUCCESS) then
                test_loc = is_int_vector_validated(
     $               bc_section(1:3),
     $               test_bc_section(1:3,k),
     $               detailled)
             else
                test_loc = .false.
             end if

             test_validated = test_validated.and.test_loc
             
             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test ('',I2,'',2) failed'')', k
             end if

          end do

        end function test_get_bc_section


        function test_analyse_grdpt_with_bc_section(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_bc_sections) :: bf_layer_bc_sections_used
          integer                    :: test_i
          integer                    :: test_j
          integer, dimension(3,3,36) :: test_grdpts_id
          integer, dimension(5,36)   :: test_bc_section
          integer, dimension(5,36)   :: test_bc_section_after
          logical, dimension(36)     :: test_compatible
          logical, dimension(36)     :: test_remove_ele
                                     
          integer                    :: k
          logical                    :: test_loc
          logical                    :: compatible
          logical                    :: remove_ele

          integer, dimension(5)      :: bc_section_in
          integer, dimension(3,3)    :: tmp_grdpts_id
          logical                    :: use_tmp_grdpts_id

          test_i = 2
          test_j = 2

          test_validated = .true.
          

          !test N_edge

          !compatible
          !  -------    
          ! | 3 3 3 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
          test_grdpts_id(:,:,1) = reshape((/
     $         interior_pt,interior_pt,interior_pt,
     $         bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $         bc_pt,bc_pt,bc_pt/),
     $         (/3,3/))
          
          test_bc_section(:,1)       = [N_edge_type,1,1,2,0]
          test_bc_section_after(:,1) = [N_edge_type,1,2,2,0]
          
          test_compatible(1)   = .true.
          test_remove_ele(1)   = .false.

          !incompatible
          !  -------    
          ! | 3 3 2 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
          test_grdpts_id(:,:,2) = reshape((/
     $         interior_pt,interior_pt,interior_pt,
     $         bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $         bc_pt,bc_pt,bc_interior_pt
     $         /),
     $         (/3,3/))
          
          test_bc_section(:,2)       = test_bc_section(:,1)
          test_bc_section_after(:,2) = test_bc_section(:,1)
          
          test_compatible(2)   = .false.
          test_remove_ele(2)   = .false.

          !remove_ele
          !  -------    
          ! | 3 3 2 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
          test_grdpts_id(:,:,3) = reshape((/
     $         interior_pt,interior_pt,interior_pt,
     $         bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $         bc_pt,bc_pt,bc_interior_pt
     $         /),
     $         (/3,3/))
          
          test_bc_section(:,3)       = [N_edge_type,1,1,1,0]
          test_bc_section_after(:,3) = test_bc_section(:,3)
          
          test_compatible(3)   = .false.
          test_remove_ele(3)   = .true.


          !test S_edge

          !compatible
          !  -------    
          ! | 1 1 1 |
          ! | 2 2*2 |
          ! | 3 3 3 |
          !  -------
          test_grdpts_id(:,:,4) = reshape((/
     $         bc_pt,bc_pt,bc_pt,
     $         bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $         interior_pt,interior_pt,interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,4)       = [S_edge_type,1,1,2,0]
          test_bc_section_after(:,4) = [S_edge_type,1,2,2,0]
          
          test_compatible(4)   = .true.
          test_remove_ele(4)   = .false.

          !incompatible
          !  -------    
          ! | 1 1 1 |
          ! | 2 2*2 |
          ! | 3 3 2 |
          !  -------
          test_grdpts_id(:,:,5) = reshape((/
     $         bc_pt,bc_pt,bc_interior_pt,
     $         bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $         interior_pt,interior_pt,interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,5)       = test_bc_section(:,4)
          test_bc_section_after(:,5) = test_bc_section(:,4)
          
          test_compatible(5)   = .false.
          test_remove_ele(5)   = .false.

          !remove_ele
          !  -------    
          ! | 1 1 1 |
          ! | 2 2*2 |
          ! | 3 3 2 |
          !  -------
          test_grdpts_id(:,:,6) = reshape((/
     $         bc_pt,bc_pt,bc_interior_pt,
     $         bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $         interior_pt,interior_pt,interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,6)       = [S_edge_type,1,1,1,0]
          test_bc_section_after(:,6) = test_bc_section(:,6)
          
          test_compatible(6)   = .false.
          test_remove_ele(6)   = .true.


          !test E_edge

          !compatible
          !  -------    
          ! | 1 2 3 |
          ! | 1 2*3 |
          ! | 1 2 3 |
          !  -------
          test_grdpts_id(:,:,7) = reshape((/
     $         interior_pt,bc_interior_pt,bc_pt,
     $         interior_pt,bc_interior_pt,bc_pt,
     $         interior_pt,bc_interior_pt,bc_pt/),
     $         (/3,3/))
          
          test_bc_section(:,7)       = [E_edge_type,1,1,2,0]
          test_bc_section_after(:,7) = [E_edge_type,1,2,2,0]
          
          test_compatible(7)   = .true.
          test_remove_ele(7)   = .false.

          !incompatible
          !  -------    
          ! | 1 2 2 |
          ! | 1 2*3 |
          ! | 1 2 3 |
          !  -------
          test_grdpts_id(:,:,8) = reshape((/
     $         interior_pt,bc_interior_pt,bc_pt,
     $         interior_pt,bc_interior_pt,bc_pt,
     $         interior_pt,bc_interior_pt,bc_interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,8)       = test_bc_section(:,7)
          test_bc_section_after(:,8) = test_bc_section(:,7)
          
          test_compatible(8)   = .false.
          test_remove_ele(8)   = .false.

          !remove_ele
          !  -------
          ! | 1 1 1 |
          ! | 2 2*2 |
          ! | 3 3 2 |
          !  -------
          test_grdpts_id(:,:,9) = reshape((/
     $         interior_pt,bc_interior_pt,bc_pt,
     $         interior_pt,bc_interior_pt,bc_pt,
     $         interior_pt,bc_interior_pt,bc_interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,9)       = [E_edge_type,0,0,2,0]
          test_bc_section_after(:,9) = test_bc_section(:,9)
          
          test_compatible(9)   = .false.
          test_remove_ele(9)   = .true.

          !test W_edge

          !compatible
          !  -------    
          ! | 3 2 1 |
          ! | 3 2*1 |
          ! | 3 2 1 |
          !  -------
          test_grdpts_id(:,:,10) = reshape((/
     $         bc_pt,bc_interior_pt,interior_pt,
     $         bc_pt,bc_interior_pt,interior_pt,
     $         bc_pt,bc_interior_pt,interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,10)      = [W_edge_type,1,1,2,0]
          test_bc_section_after(:,10)= [W_edge_type,1,2,2,0]
          
          test_compatible(10)   = .true.
          test_remove_ele(10)   = .false.

          !incompatible
          !  -------    
          ! | 2 2 1 |
          ! | 3 2*1 |
          ! | 3 2 1 |
          !  -------
          test_grdpts_id(:,:,11) = reshape((/
     $         bc_pt,bc_interior_pt,interior_pt,
     $         bc_pt,bc_interior_pt,interior_pt,
     $         bc_interior_pt,bc_interior_pt,interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,11)      = test_bc_section(:,10)
          test_bc_section_after(:,11)= test_bc_section(:,10)
          
          test_compatible(11)  = .false.
          test_remove_ele(11)  = .false.

          !remove_ele
          !  -------
          ! | 2 2 1 |
          ! | 3 2*1 |
          ! | 3 2 1 |
          !  -------
          test_grdpts_id(:,:,12) = reshape((/
     $         bc_pt,bc_interior_pt,interior_pt,
     $         bc_pt,bc_interior_pt,interior_pt,
     $         bc_interior_pt,bc_interior_pt,interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,12)       = [W_edge_type,0,0,2,0]
          test_bc_section_after(:,12) = test_bc_section(:,12)
          
          test_compatible(12)  = .false.
          test_remove_ele(12)  = .true.

          !test NE_edge

          !compatible
          !  -------    
          ! | 2 3 3 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
          test_grdpts_id(:,:,13) = reshape((/
     $         interior_pt, interior_pt, interior_pt,
     $         bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $         bc_interior_pt, bc_pt, bc_pt/),
     $         (/3,3/))
          
          test_bc_section(:,13)      = [NW_edge_type,1,2,0,0]
          test_bc_section_after(:,13)= test_bc_section(:,13)
          
          test_compatible(13)   = .true.
          test_remove_ele(13)   = .false.

          !incompatible
          !  -------
          ! | 3 3 3 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
          test_grdpts_id(:,:,14) = reshape((/
     $         interior_pt, interior_pt, interior_pt,
     $         bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $         bc_pt, bc_pt, bc_pt/),
     $         (/3,3/))
          
          test_bc_section(:,14)      = [NW_edge_type,0,2,0,0]
          test_bc_section_after(:,14)= test_bc_section(:,14)
          
          test_compatible(14)   = .false.
          test_remove_ele(14)   = .false.

          !remove_ele
          !  -------
          ! | 3 3 3 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
          test_grdpts_id(:,:,15) = reshape((/
     $         interior_pt, interior_pt, interior_pt,
     $         bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $         bc_pt, bc_pt, bc_pt/),
     $         (/3,3/))
          
          test_bc_section(:,15)      = [NW_edge_type,0,0,0,0]
          test_bc_section_after(:,15)= test_bc_section(:,15)
          
          test_compatible(15)   = .false.
          test_remove_ele(15)   = .true.


          !test NW_edge

          !compatible
          !  -------    
          ! | 3 2 1 |
          ! | 2 2*1 |
          ! | 1 1 1 |
          !  -------
          test_grdpts_id(:,:,16) = reshape((/
     $         interior_pt, interior_pt, interior_pt,
     $         bc_interior_pt, bc_interior_pt, interior_pt,
     $         bc_pt, bc_interior_pt, interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,16)      = [NW_edge_type,1,2,0,0]
          test_bc_section_after(:,16)= test_bc_section(:,16)
          
          test_compatible(16)   = .true.
          test_remove_ele(16)   = .false.

          !incompatible
          !  -------
          ! | 3 3 3 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
          test_grdpts_id(:,:,17) = reshape((/
     $         interior_pt, interior_pt, interior_pt,
     $         bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $         bc_pt, bc_pt, bc_pt/),
     $         (/3,3/))
          
          test_bc_section(:,17)      = [NW_edge_type,0,2,0,0]
          test_bc_section_after(:,17)= test_bc_section(:,17)
          
          test_compatible(14)   = .false.
          test_remove_ele(14)   = .false.

          !remove_ele
          !  -------
          ! | 3 3 3 |
          ! | 2 2*2 |
          ! | 1 1 1 |
          !  -------
          test_grdpts_id(:,:,18) = reshape((/
     $         interior_pt, interior_pt, interior_pt,
     $         bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $         bc_pt, bc_pt, bc_pt/),
     $         (/3,3/))
          
          test_bc_section(:,18)      = [NW_edge_type,0,0,0,0]
          test_bc_section_after(:,18)= test_bc_section(:,18)
          
          test_compatible(18)   = .false.
          test_remove_ele(18)   = .true.


          !test SW_edge

          !compatible
          !  -------    
          ! | 1 1 1 |
          ! | 2 2*1 |
          ! | 3 2 1 |
          !  -------
          test_grdpts_id(:,:,19) = reshape((/
     $         bc_pt, bc_interior_pt, interior_pt,
     $         bc_interior_pt, bc_interior_pt, interior_pt,
     $         interior_pt, interior_pt, interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,19)      = [SW_edge_type,1,1,0,0]
          test_bc_section_after(:,19)= test_bc_section(:,19)
          
          test_compatible(19)   = .true.
          test_remove_ele(19)   = .false.

          !incompatible
          !  -------
          ! | 2 1 1 |
          ! | 2 2*2 |
          ! | 3 3 3 |
          !  -------
          test_grdpts_id(:,:,20) = reshape((/
     $         bc_pt, bc_pt, bc_pt,
     $         bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $         bc_interior_pt, interior_pt, interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,20)      = [SW_edge_type,0,2,0,0]
          test_bc_section_after(:,20)= test_bc_section(:,20)
          
          test_compatible(20) = .false.
          test_remove_ele(20) = .false.

          !remove_ele
          !  -------
          ! | 1 1 1 |
          ! | 2 2*2 |
          ! | 3 3 3 |
          !  -------
          test_grdpts_id(:,:,21) = reshape((/
     $         bc_pt, bc_pt, bc_pt,
     $         bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $         interior_pt, interior_pt, interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,21)      = [SW_edge_type,0,0,0,0]
          test_bc_section_after(:,21)= test_bc_section(:,21)
          
          test_compatible(21)   = .false.
          test_remove_ele(21)   = .true.

          !test SW_edge

          !compatible
          !  -------    
          ! | 1 1 1 |
          ! | 1 2*2 |
          ! | 2 2 3 |
          !  -------
          test_grdpts_id(:,:,22) = reshape((/
     $         bc_interior_pt, bc_interior_pt, bc_pt,
     $         interior_pt, bc_interior_pt, bc_interior_pt,
     $         interior_pt, interior_pt, interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,22)      = [SE_edge_type,2,1,0,0]
          test_bc_section_after(:,22)= test_bc_section(:,22)
          
          test_compatible(22)   = .true.
          test_remove_ele(22)   = .false.

          !incompatible
          !  -------
          ! | 1 1 2 |
          ! | 2 2*2 |
          ! | 3 3 3 |
          !  -------
          test_grdpts_id(:,:,23) = reshape((/
     $         bc_pt, bc_pt, bc_pt,
     $         bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $         interior_pt, interior_pt, bc_interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,23)      = [SE_edge_type,3,2,0,0]
          test_bc_section_after(:,23)= test_bc_section(:,23)
          
          test_compatible(23) = .false.
          test_remove_ele(23) = .false.

          !remove_ele
          !  -------
          ! | 1 1 2 |
          ! | 2 2*2 |
          ! | 3 3 3 |
          !  -------
          test_grdpts_id(:,:,24) = reshape((/
     $         bc_pt, bc_pt, bc_pt,
     $         bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $         interior_pt, interior_pt, bc_interior_pt/),
     $         (/3,3/))
          
          test_bc_section(:,24)      = [SE_edge_type,0,0,0,0]
          test_bc_section_after(:,24)= test_bc_section(:,24)
          
          test_compatible(24)   = .false.
          test_remove_ele(24)   = .true.


          ! test w/o tmp_grdpts_id
          !============================================================
          use_tmp_grdpts_id = .false.
          do k=1,24

             !input
             bc_section_in = test_bc_section(:,k)

             !output
             compatible = bf_layer_bc_sections_used%analyse_grdpt_with_bc_section(
     $            test_i,
     $            test_j,
     $            test_grdpts_id(:,:,k),
     $            tmp_grdpts_id, use_tmp_grdpts_id,
     $            bc_section_in,
     $            remove_ele)

             !validation
             test_loc = is_int_vector_validated(bc_section_in,test_bc_section_after(:,k))
             test_loc = test_loc.and.(compatible.eqv.test_compatible(k))
             test_loc = test_loc.and.(remove_ele.eqv.test_remove_ele(k))
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'',1) failed'')',k
                print '('' - bc_section: '',L1)', is_int_vector_validated(test_bc_section(:,k),test_bc_section_after(:,k))
                print '('' - compatible: '',L1)', compatible.eqv.test_compatible(k)
                print '('' - remove_ele: '',L1)', remove_ele.eqv.test_remove_ele(k)
                print '()'
             end if

          end do


          ! test w/ tmp_grdpts_id
          !============================================================
          use_tmp_grdpts_id = .true.
          do k=1,24

             !input
             tmp_grdpts_id = test_grdpts_id(:,:,k)
             bc_section_in = test_bc_section(:,k)

             !output
             compatible = bf_layer_bc_sections_used%analyse_grdpt_with_bc_section(
     $            test_i,
     $            test_j,
     $            test_grdpts_id(:,:,1),
     $            tmp_grdpts_id, use_tmp_grdpts_id,
     $            bc_section_in,
     $            remove_ele)

             !validation
             test_loc = is_int_vector_validated(bc_section_in,test_bc_section_after(:,k))
             test_loc = test_loc.and.(compatible.eqv.test_compatible(k))
             test_loc = test_loc.and.(remove_ele.eqv.test_remove_ele(k))
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'',2) failed'')',k
                print '('' - bc_section: '',L1)', is_int_vector_validated(test_bc_section(:,k),test_bc_section_after(:,k))
                print '('' - compatible: '',L1)', compatible.eqv.test_compatible(k)
                print '('' - remove_ele: '',L1)', remove_ele.eqv.test_remove_ele(k)
                print '()'
             end if

          end do

        end function test_analyse_grdpt_with_bc_section

c$$$        subroutine test_analyse_grdpt(detailled)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$
c$$$          type(bf_layer_bc_sections)           :: bc_sections
c$$$          integer                              :: i,j,k
c$$$          integer                              :: test_i_min
c$$$          integer                              :: test_i_max
c$$$          integer                              :: test_j_min
c$$$          integer                              :: test_j_max
c$$$          integer, dimension(:,:), allocatable :: grdpts_id
c$$$          integer, dimension(:,:), allocatable :: test_bc_sections_sorted
c$$$
c$$$          integer, dimension(:,:), allocatable :: bc_sections_sorted
c$$$
c$$$          logical :: same
c$$$          logical :: test_global
c$$$          
c$$$          logical :: ierror
c$$$
c$$$          print '(''test_analyse_grdpt()'')'          
c$$$
c$$$          test_global = .true.
c$$$
c$$$          do k=1,3
c$$$
c$$$             !initialize object for boundary layers
c$$$             call bc_sections%ini()
c$$$
c$$$             !create test
c$$$             call make_test_analyse_grdpt(
c$$$     $            k,
c$$$     $            test_i_min, test_i_max,
c$$$     $            test_j_min, test_j_max,
c$$$     $            grdpts_id,
c$$$     $            test_bc_sections_sorted)
c$$$
c$$$             !make test
c$$$             do j=test_j_min, test_j_max
c$$$                do i=test_i_min,test_i_max
c$$$
c$$$                   if(grdpts_id(i,j).eq.bc_interior_pt) then
c$$$                      call bc_sections%analyse_grdpt(i,j,grdpts_id,ierror)
c$$$                   end if
c$$$
c$$$                end do
c$$$             end do
c$$$
c$$$             call bc_sections%sort_bc_sections(bc_sections_sorted)
c$$$
c$$$             !compare results
c$$$             same = compare_bc_sections_sorted(
c$$$     $            bc_sections_sorted,
c$$$     $            test_bc_sections_sorted)
c$$$
c$$$             if(detailled) then
c$$$                print '(''test '',I2,'': '',L1)', k, same
c$$$                if(.not.same) then
c$$$                   call print_bc_sections(bc_sections_sorted)
c$$$                end if
c$$$             end if
c$$$
c$$$             test_global = test_global.and.same
c$$$
c$$$             !deallocate intermediate variables
c$$$             if(allocated(test_bc_sections_sorted)) then
c$$$                deallocate(test_bc_sections_sorted)
c$$$             end if
c$$$
c$$$             !deallocate the object with the boundary
c$$$             !layers
c$$$             call bc_sections%deallocate_tables()
c$$$
c$$$          end do
c$$$
c$$$          print '(''test_validated: '', L1)', test_global
c$$$          print '()'
c$$$
c$$$        end subroutine test_analyse_grdpt
c$$$
c$$$        
c$$$        subroutine print_bc_sections(bc_sections)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, dimension(:,:), intent(in) :: bc_sections
c$$$
c$$$          integer :: i
c$$$
c$$$          do i=1, size(bc_sections,2)
c$$$
c$$$             select case(bc_sections(1,i))
c$$$
c$$$               case(N_edge_type)
c$$$                  print '(I2,'' N_edge    [i,j]: '',2I3,'' i_max='',I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i),
c$$$     $                 bc_sections(4,i)
c$$$
c$$$               case(S_edge_type)
c$$$                  print '(I2,'' S_edge    [i,j]: '',2I3,'' i_max='',I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i),
c$$$     $                 bc_sections(4,i)
c$$$
c$$$               case(E_edge_type)
c$$$                  print '(I2,'' E_edge    [i,j]: '',2I3,'' j_max='',I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i),
c$$$     $                 bc_sections(4,i)
c$$$
c$$$               case(W_edge_type)
c$$$                  print '(I2,'' W_edge    [i,j]: '',2I3,'' j_max='',I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i),
c$$$     $                 bc_sections(4,i)
c$$$
c$$$               case(NE_corner_type)
c$$$                  print '(I2,'' NE_corner [i,j]: '',2I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i)
c$$$
c$$$               case(NW_corner_type)
c$$$                  print '(I2,'' NW_corner [i,j]: '',2I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i)
c$$$
c$$$               case(SW_corner_type)
c$$$                  print '(I2,'' SW_corner [i,j]: '',2I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i)
c$$$
c$$$               case(SE_corner_type)
c$$$                  print '(I2,'' SE_corner [i,j]: '',2I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i)
c$$$
c$$$               case(NE_edge_type)
c$$$                  print '(I2,'' NE_edge   [i,j]: '',2I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i)
c$$$
c$$$               case(NW_edge_type)
c$$$                  print '(I2,'' NW_edge   [i,j]: '',2I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i)
c$$$
c$$$               case(SW_edge_type)
c$$$                  print '(I2,'' SW_edge   [i,j]: '',2I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i)
c$$$
c$$$               case(SE_edge_type)
c$$$                  print '(I2,'' SE_edge   [i,j]: '',2I3)',
c$$$     $                 i,
c$$$     $                 bc_sections(2,i), bc_sections(3,i)
c$$$
c$$$             end select
c$$$
c$$$          end do
c$$$
c$$$        end subroutine print_bc_sections
c$$$
c$$$
c$$$        subroutine make_test_analyse_grdpt(
c$$$     $     test_id,
c$$$     $     test_i_min, test_i_max,
c$$$     $     test_j_min, test_j_max,
c$$$     $     grdpts_id,
c$$$     $     test_bc_sections_sorted)
c$$$
c$$$          implicit none
c$$$
c$$$          integer                             , intent(in)  :: test_id
c$$$          integer                             , intent(out) :: test_i_min
c$$$          integer                             , intent(out) :: test_i_max
c$$$          integer                             , intent(out) :: test_j_min
c$$$          integer                             , intent(out) :: test_j_max
c$$$          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
c$$$          integer, dimension(:,:), allocatable, intent(out) :: test_bc_sections_sorted
c$$$
c$$$
c$$$          integer :: i,j
c$$$
c$$$          
c$$$          select case(test_id)
c$$$
c$$$            !  ---- | ----- | ----
c$$$            ! |     |       |     |
c$$$            ! | 2 2 | 2 2 2 | 2 2 |
c$$$            ! | 1 1 | 1 1 1 | 1 1 |
c$$$            ! | 0 0 | 0 0 0 | 0 0 |
c$$$            ! ---------------------
c$$$            ! | 0 0 | 0 0 0 | 0 0 |
c$$$            ! | 0 0 | 0 0 0 | 0 0 |
c$$$            !  ---- | ----- | ----
c$$$            !
c$$$            case(1)
c$$$               test_i_min = 3
c$$$               test_i_max = 5
c$$$               test_j_min = 3
c$$$               test_j_max = 5
c$$$
c$$$               allocate(grdpts_id(7,6))
c$$$               
c$$$               do j=1,3
c$$$                  do i=1,7
c$$$                     grdpts_id(i,j) = interior_pt
c$$$                  end do
c$$$               end do
c$$$
c$$$               j=4
c$$$               do i=1,7
c$$$                  grdpts_id(i,j) = bc_interior_pt
c$$$               end do
c$$$
c$$$               j=5
c$$$               do i=1,7
c$$$                  grdpts_id(i,j) = bc_pt
c$$$               end do
c$$$
c$$$               j=6
c$$$               do i=1,7
c$$$                  grdpts_id(i,j) = no_pt
c$$$               end do
c$$$
c$$$               allocate(test_bc_sections_sorted(4,1))
c$$$               test_bc_sections_sorted = reshape(
c$$$     $              (/
c$$$     $              N_edge_type,3,4,5
c$$$     $              /),
c$$$     $              (/4,1/))
c$$$
c$$$
c$$$            !  ---- | ------- | ----
c$$$            ! |     |         |     |
c$$$            ! |     | 2 2 2 2 | 2 2 |
c$$$            ! |     | 2 1 1 1 | 1 1 |
c$$$            ! |     | 2 1 0 0 | 0 0 |
c$$$            ! -----------------------
c$$$            ! |     | 2 1 0 0 | 0 0 |
c$$$            ! |     | 2 1 0 0 | 0 0 |
c$$$            !  ---- | ------- | ----
c$$$            !
c$$$            case(2)
c$$$               test_i_min = 3
c$$$               test_i_max = 6
c$$$               test_j_min = 3
c$$$               test_j_max = 5
c$$$
c$$$               allocate(grdpts_id(8,6))
c$$$
c$$$               do j=1,3
c$$$                  do i=1,2
c$$$                     grdpts_id(i,j) = no_pt
c$$$                  end do
c$$$                  grdpts_id(3,j) = bc_pt
c$$$                  grdpts_id(4,j) = bc_interior_pt
c$$$                  do i=5,8
c$$$                     grdpts_id(i,j) = interior_pt
c$$$                  end do
c$$$               end do
c$$$
c$$$               j=4
c$$$               do i=1,2
c$$$                  grdpts_id(i,j) = no_pt
c$$$               end do
c$$$               grdpts_id(3,j) = bc_pt
c$$$               do i=4,8
c$$$                  grdpts_id(i,j) = bc_interior_pt
c$$$               end do
c$$$               
c$$$               j=5
c$$$               do i=1,2
c$$$                  grdpts_id(i,j) = no_pt
c$$$               end do
c$$$               do i=3,8
c$$$                  grdpts_id(i,j) = bc_pt
c$$$               end do
c$$$
c$$$               j=6
c$$$               do i=1,8
c$$$                  grdpts_id(i,j) = no_pt
c$$$               end do
c$$$
c$$$               allocate(test_bc_sections_sorted(4,3))
c$$$               test_bc_sections_sorted = reshape(
c$$$     $              (/
c$$$     $              W_edge_type,3,3,3,
c$$$     $              NW_corner_type,3,4,0,
c$$$     $              N_edge_type,5,4,6
c$$$     $              /),
c$$$     $              (/4,3/))
c$$$
c$$$            !    |     |         2 2 2 2 2 2 2 2 2                |    |
c$$$            !    |     |         2 1 1 1 1 1 1 1 2                |    |
c$$$            !    |     |         2 1 0 0 0 0 0 1 2                |    |
c$$$            !    |     |         2 1 0       0 1 2                |    |
c$$$            !18- |     |         2 1 0       0 1 2 2 2 2 2        |    |
c$$$            !17- |     | 2 2 2 2 2 1 0       0 1 1 1 1 1 2        |    |
c$$$            !    |     | 2 1 1 1 1 1 0       0 0 0 0 0 1 2        |    |
c$$$            !15- |     | 2 1 0 0 0 0 0               0 1 2        |    |
c$$$            !    |     | 2 1 0                       0 1 2        |    |
c$$$            !13- |     | 2 1 0                       0 1 2 2 2 2  |    |
c$$$            !    |     | 2 1 0                       0 1 1 1 1 2  |    |
c$$$            !11- |     | 2 1 0 0 0                   0 0 0 0 1 2  |    |
c$$$            !10- |     | 2 1 1 1 0                         0 1 2  |    |
c$$$            !    |     | 2 2 2 1 0             0 0 0 0 0 0 0 1 2  |    |
c$$$            ! 8- |     |     2 1 0 0           0 1 1 1 1 1 1 1 2  |    |
c$$$            !    |     |     2 1 1 0           0 1 2 2 2 2 2 2 2  |    |
c$$$            !    |     |     2 2 1 0           0 1 2              |    |
c$$$            ! 5- |     |       2 1 0           0 1 2              |    |
c$$$            !    |     | 2 2 2 2 1 0           0 1 2              |    |
c$$$            !    |     | 2 1 1 1 1 0           0 1 2              |    |
c$$$            !    -------------------------------------------------|----|
c$$$            !    |     | 2 1 0 0               0 1 2              |    |
c$$$            !    |     | 2 1 0 0               0 1 2              |    |
c$$$            !     ---- | ----------------------------------------- ---  
c$$$            !            |     | | |           |     |         |        
c$$$            !            3     6 7 8          14    17        22
c$$$            ! --------------------------------------------------
c$$$            case(3)
c$$$
c$$$               test_i_min = 3
c$$$               test_i_max = 23
c$$$               test_j_min = 3
c$$$               test_j_max = 22
c$$$
c$$$               allocate(grdpts_id(25,22))
c$$$
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c$$$     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
c$$$     $              0, 0, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
c$$$     $              0, 0, 3, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
c$$$     $              0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
c$$$     $              0, 0, 0, 0, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c$$$     $              0, 0, 0, 0, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0,
c$$$     $              0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 0, 0, 0,
c$$$     $              0, 0, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0,
c$$$     $              0, 0, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0,
c$$$     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0,
c$$$     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 0, 0, 0,
c$$$     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 0, 0, 0,
c$$$     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0,
c$$$     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0,
c$$$     $              0, 0, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0,
c$$$     $              0, 0, 3, 3, 3, 3, 3, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0, 0, 0,
c$$$     $              0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0,
c$$$     $              0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c$$$     $              0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c$$$     $              0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c$$$     $              0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
c$$$     $              /),
c$$$     $              (/25,22/))
c$$$
c$$$               allocate(test_bc_sections_sorted(4,30))
c$$$               test_bc_sections_sorted = reshape(
c$$$     $              (/
c$$$     $              NW_corner_type,3,3,0,
c$$$     $              N_edge_type,5,3,5,
c$$$     $              NW_edge_type,6,3,0,
c$$$     $              E_edge_type,15,3,6,
c$$$     $              W_edge_type,6,5,5,
c$$$     $              SW_corner_type,5,6,0,
c$$$     $              SW_edge_type,6,6,0,
c$$$     $              SE_edge_type,15,7,0,
c$$$     $              S_edge_type,17,7,20,
c$$$     $              SE_corner_type,21,7,0,
c$$$     $              W_edge_type,5,8,8,
c$$$     $              SW_corner_type,3,9,0,
c$$$     $              SW_edge_type,5,9,0,
c$$$     $              E_edge_type,21,9,11,
c$$$     $              W_edge_type,3,11,15,
c$$$     $              NE_edge_type,18,12,0,
c$$$     $              N_edge_type,20,12,20,
c$$$     $              NE_corner_type,21,12,0,
c$$$     $              E_edge_type,18,14,16,
c$$$     $              NW_corner_type,3,16,0,
c$$$     $              N_edge_type,5,16,6,
c$$$     $              NW_edge_type,7,16,0,
c$$$     $              NE_edge_type,14,17,0,
c$$$     $              N_edge_type,16,17,17,
c$$$     $              NE_corner_type,18,17,0,
c$$$     $              W_edge_type,7,18,20,
c$$$     $              E_edge_type,14,19,20,
c$$$     $              NW_corner_type,7,21,0,
c$$$     $              N_edge_type,9,21,13,
c$$$     $              NE_corner_type,14,21,0
c$$$     $              /),
c$$$     $              (/4,30/))
c$$$               
c$$$            case default
c$$$               print '(''test_bf_layer_bc_sections'')'
c$$$               print '(''make_test_analyse_grdpt'')'
c$$$               print '(''case : '', I2)', test_id
c$$$               stop
c$$$          end select
c$$$
c$$$c$$$          allocate(test_bc_sections_buffer(1,1))
c$$$c$$$          deallocate(test_bc_sections_buffer)
c$$$
c$$$        end subroutine make_test_analyse_grdpt
c$$$
c$$$
c$$$        subroutine test_add_overlap_between_corners_and_anti_corners(detailled)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$
c$$$
c$$$          type(bf_layer_bc_sections) :: bf_layer_bc_sections_used
c$$$          integer, dimension(4,29)   :: bc_sections_sorted
c$$$          integer, dimension(4,29)   :: test_bc_sections_sorted
c$$$
c$$$          logical :: test_validated
c$$$
c$$$
c$$$          !initialization of bc_sections_sorted
c$$$          bc_sections_sorted(:, 1) = [NE_edge_type  ,1, 1,0] !-> no_overlap
c$$$          bc_sections_sorted(:, 2) = [E_edge_type   ,1, 1,3]
c$$$          bc_sections_sorted(:, 3) = [NW_edge_type  ,1, 2,0] !-> N_overlap
c$$$          bc_sections_sorted(:, 4) = [N_edge_type   ,2, 2,3]
c$$$          bc_sections_sorted(:, 5) = [NW_corner_type,1, 3,0] !-> for N_overlap
c$$$          bc_sections_sorted(:, 6) = [S_edge_type   ,2, 3,3]
c$$$          bc_sections_sorted(:, 7) = [SW_edge_type  ,3, 3,0] !-> NE_overlap
c$$$          bc_sections_sorted(:, 8) = [SE_corner_type,4, 3,0] !-> for E_overlap
c$$$          bc_sections_sorted(:, 9) = [SW_corner_type,3, 4,0] !-> for N_overlap
c$$$          bc_sections_sorted(:,10) = [W_edge_type   ,2, 5,3]
c$$$          bc_sections_sorted(:,11) = [SW_corner_type,6, 5,0] !-> for W_overlap
c$$$          bc_sections_sorted(:,12) = [SE_edge_type  ,7, 5,0] !-> NW_overlap
c$$$          bc_sections_sorted(:,13) = [SE_corner_type,7, 6,0] !-> for N_overlap
c$$$
c$$$          bc_sections_sorted(:,14) = [SW_corner_type,6, 8,0] !-> for W_overlap
c$$$          bc_sections_sorted(:,15) = [SE_edge_type  ,7, 8,0] !-> W_overlap
c$$$          bc_sections_sorted(:,16) = [SE_edge_type  ,8,10,0] !-> E_overlap
c$$$          bc_sections_sorted(:,17) = [W_edge_type   ,9,10,9]
c$$$          bc_sections_sorted(:,18) = [SW_corner_type,9,10,0] !-> for E_overlap
c$$$
c$$$          bc_sections_sorted(:,19) = [NW_corner_type,1,12,0] !-> for S_overlap
c$$$          bc_sections_sorted(:,20) = [N_edge_type   ,2,12,3]
c$$$          bc_sections_sorted(:,21) = [NW_edge_type  ,1,13,0] !-> S_overlap
c$$$          bc_sections_sorted(:,22) = [S_edge_type   ,2,13,3]
c$$$          bc_sections_sorted(:,23) = [SW_corner_type,3,15,0] !-> for S_overlap
c$$$          bc_sections_sorted(:,24) = [SW_edge_type  ,3,16,0] !-> SE_overlap
c$$$          bc_sections_sorted(:,25) = [SE_corner_type,4,16,0] !-> for E_overlap
c$$$          bc_sections_sorted(:,26) = [W_edge_type   ,2,17,3]
c$$$          bc_sections_sorted(:,27) = [SE_corner_type,7,17,0] !-> for S_overlap
c$$$          bc_sections_sorted(:,28) = [SW_corner_type,6,18,0] !-> for W_overlap
c$$$          bc_sections_sorted(:,29) = [SE_edge_type  ,7,18,0] !-> SW_overlap
c$$$
c$$$          
c$$$
c$$$          !initialization of test_bc_sections_sorted
c$$$          test_bc_sections_sorted       = bc_sections_sorted
c$$$
c$$$          test_bc_sections_sorted(4,3)  = N_overlap
c$$$          test_bc_sections_sorted(4,7)  = NE_overlap
c$$$          test_bc_sections_sorted(4,12) = NW_overlap
c$$$
c$$$          test_bc_sections_sorted(4,15) = W_overlap
c$$$          test_bc_sections_sorted(4,16) = E_overlap
c$$$
c$$$          test_bc_sections_sorted(4,21) = S_overlap
c$$$          test_bc_sections_sorted(4,24) = SE_overlap
c$$$          test_bc_sections_sorted(4,29) = SW_overlap
c$$$          
c$$$
c$$$          !computation of the overlap
c$$$          call bf_layer_bc_sections_used%add_overlap_between_corners_and_anti_corners(
c$$$     $         bc_sections_sorted)
c$$$
c$$$
c$$$          !check the result of the computation of the overlap
c$$$          test_validated = is_int_matrix_validated(
c$$$     $         bc_sections_sorted,
c$$$     $         test_bc_sections_sorted,
c$$$     $         detailled)
c$$$
c$$$          print '(''add_overlap_between_corners_and_anti_corners'')'
c$$$          print '(''test_validated: '',L1)', test_validated
c$$$          print '()'
c$$$
c$$$        end subroutine test_add_overlap_between_corners_and_anti_corners


        function test_create_tmp_grdpts_id_for_analyse(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_layer_bc_sections)       :: bf_layer_bc_sections_used
          integer(ikind), dimension(2,2)   :: bf_alignment
          integer       , dimension(6,5)   :: grdpts_id
          integer(ikind), dimension(2,2)   :: central_coords
          integer       , dimension(3,3,2) :: test_grdpts_id
          integer       , dimension(3,3)   :: tmp_grdpts_id

          integer :: k
          logical :: test_loc

          
          test_validated = .true.


          bf_alignment = reshape((/
     $         align_E-2, align_N,
     $         align_E-1, align_N/),
     $         (/2,2/))

          grdpts_id = reshape((/
     $         1,1,1,1,2,3,
     $         1,1,1,1,2,3,
     $         2,2,2,2,2,3,
     $         3,2,2,2,3,3,
     $         3,3,3,3,3,0/),
     $         (/6,5/))

          
          central_coords(:,1) = [1,3]

          test_grdpts_id(:,:,1) = reshape((/
     $         1,1,1,
     $         2,2,2,
     $         3,3,2/),
     $         (/3,3/))

          central_coords(:,2) = [1,4]

          test_grdpts_id(:,:,2) = reshape((/
     $         2,2,2,
     $         3,3,2,
     $         0,3,3/),
     $         (/3,3/))

          do k=1,2

             !output
             tmp_grdpts_id = bf_layer_bc_sections_used%create_tmp_grdpts_id_for_analyse(
     $            bf_alignment,
     $            central_coords(:,k),
     $            grdpts_id)

             !validation
             test_loc = is_int_matrix_validated(
     $            tmp_grdpts_id,
     $            test_grdpts_id(:,:,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test '',I2,'' failed'')', k
             end if
                
          end do

        end function test_create_tmp_grdpts_id_for_analyse


        function test_finalize_bc_sections(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_layer_bc_sections)           :: bf_layer_bc_sections_used
          integer, dimension(26,22)            :: grdpts_id
          integer                              :: i_min
          integer                              :: i_max
          integer                              :: j_min
          integer                              :: j_max
          integer, dimension(2)                :: x_borders
          integer, dimension(2)                :: y_borders
          integer, dimension(5,60)             :: test_bc_sections
          integer                              :: i
          integer                              :: j
          logical                              :: ierror          
          integer, dimension(:,:), allocatable :: sorted_bc_sections

          integer(ikind), dimension(2,2) :: bf_alignment


          !22- |     |                    3 3 3 3 3 3                     |    |
          !    |     |            3 3 3 3 3 2 2 2 2 3 3 3 3 3             |    |
          !20- |     |            3 2 2 2 2 2     2 2 2 2 2 3             |    |
          !    |     |    3 3 3 3 3 2                     2 3 3 3 3 3     |    |
          !18- |     |  3 3 2 2 2 2 2                     2 2 2 2 2 3 3   |    |
          !17- |     |3 3 2 2                                     2 2 3 3 |    |
          !    |     |3 2 2                                         2 2 3 |    |
          !15- |     |3 2                                             2 3 |    |
          !    |     |3 2                                             2 3 |    |
          !13- |     |3 2 2                                         2 2 3 |    |
          !    |     |3 3 2 2                                     2 2 3 3 |    |
          !11- |     |  3 3 2 2 2 2 2                     2 2 2 2 2 3 3   |    |
          !10- |     |    3 3 3 3 3 2                     2 3 3 3 3 3     |    |
          !    |     |            3 2                     2 3             |    |
          ! 8- |     |            3 2 2 2 2 2     2 2 2 2 2 3             |    |
          !    |     |            3 3 3 3 3 2     2 3 3 3 3 3             |    |
          !    |     |                    3 2     2 3                     |    |
          !    |     |        3 3 3 3 3 3 3 2     2 3 3 3 3 3 3 3         |    |
          ! 4- |     |        3 2 2 2 2 2 2 2     2 2 2 2 2 2 2 3         |    |
          !    |     |        3 2                             2 3         |    |
          ! 2- |     |        3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3         |    |
          !    |     |        3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3         |    |
          !     ---- | ---------------------------------------------------------
          !           | | |   |   |       |       |         | |     | | |
          !           1 2 3   5   7       11      15        2021    242526
          ! --------------------------------------------------
          i_min = 2
          i_max = 25
          j_min = 2
          j_max = 21

          x_borders = [3,26]
          y_borders = [3,22]

          grdpts_id = reshape(
     $         (/
     $         0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0,
     $         0, 0, 0, 0, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0,
     $         0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0,
     $         0, 0, 0, 0, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0,
     $         0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 2, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 2, 1, 1, 2, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0,
     $         0, 0, 3, 3, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 3, 0, 0,
     $         0, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 0,
     $         3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3,
     $         3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3,
     $         3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3,
     $         3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3,
     $         3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3,
     $         3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3,
     $         0, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 0,
     $         0, 0, 3, 3, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 3, 0, 0,
     $         0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/),
     $         (/26,22/))
          
          test_bc_sections = reshape((/
     $         SW_corner_type, 5 ,  1, no_overlap                , NS_overlap,
     $         S_edge_type   , 7 ,  1, 20                        , NS_overlap,
     $         SE_corner_type, 21,  1, no_overlap                , NS_overlap,
     $         W_edge_type   , 5 ,  3, 3                         , no_overlap,
     $         E_edge_type   , 21,  3, 3                         , no_overlap,
     $         NW_corner_type, 5 ,  4, no_overlap                , no_overlap,
     $         N_edge_type   , 7 ,  4, 10                        , no_overlap,
     $         NW_edge_type  , 11,  4, no_overlap                , no_overlap,
     $         NE_edge_type  , 15,  4, no_overlap                , no_overlap,
     $         N_edge_type   , 17,  4, 20                        , no_overlap,
     $         NE_corner_type, 21,  4, no_overlap                , no_overlap,
     $         W_edge_type   , 11,  6, 6                         , no_overlap,
     $         E_edge_type   , 15,  6, 6                         , no_overlap,
     $         SW_corner_type, 7 ,  7, no_overlap                , no_overlap,
     $         S_edge_type   , 9 ,  7, 10                        , no_overlap,
     $         SW_edge_type  , 11,  7, no_overlap                , no_overlap,
     $         SE_edge_type  , 15,  7, no_overlap                , no_overlap,
     $         S_edge_type   , 17,  7, 18                        , no_overlap,
     $         SE_corner_type, 19,  7, no_overlap                , no_overlap,
     $         W_edge_type   , 7 ,  9, 9                         , no_overlap,
     $         E_edge_type   , 19,  9, 9                         , no_overlap,
     $         SW_corner_type, 3 , 10, cpt2normal_and_cpt3not    , no_overlap,
     $         S_edge_type   , 5 , 10, 6                         , no_overlap,
     $         SW_edge_type  , 7 , 10, no_overlap                , no_overlap,
     $         SE_edge_type  , 19, 10, no_overlap                , no_overlap,
     $         S_edge_type   , 21, 10, 22                        , no_overlap,
     $         SE_corner_type, 23, 10, cpt1normal_and_cpt4not    , no_overlap,
     $         SW_corner_type, 2 , 11, cpt2overlap_and_cpt3not   , W_overlap,
     $         SW_edge_type  , 3 , 11, cpt2normal_and_cpt3not    , SW_overlap,
     $         SE_edge_type  , 23, 11, cpt1normal_and_cpt4not    , SE_overlap,
     $         SE_corner_type, 24, 11, cpt1overlap_and_cpt4not   , no_overlap,
     $         SW_corner_type, 1 , 12, cpt2overlap_and_cpt3normal, EW_overlap,
     $         SW_edge_type  , 2 , 12, cpt2overlap_and_cpt3normal, SW_overlap,
     $         SE_edge_type  , 24, 12, cpt1overlap_and_cpt4normal, SE_overlap,
     $         SE_corner_type, 25, 12, cpt1overlap_and_cpt4normal, no_overlap,
     $         W_edge_type   , 1 , 14, 15                        , EW_overlap,
     $         E_edge_type   , 25, 14, 15                        , no_overlap,
     $         NW_corner_type, 1 , 16, cpt1normal_and_cpt4not    , EW_overlap,
     $         NW_edge_type  , 2,  16, cpt1normal_and_cpt4not    , NW_overlap,
     $         NE_edge_type  , 24, 16, cpt2normal_and_cpt3not    , NE_overlap,
     $         NE_corner_type, 25, 16, cpt2normal_and_cpt3not    , no_overlap,
     $         NW_corner_type, 2 , 17, cpt1overlap_and_cpt4not   , W_overlap,
     $         NW_edge_type  , 3 , 17, cpt1overlap_and_cpt4normal, NW_overlap,
     $         NE_edge_type  , 23, 17, cpt2overlap_and_cpt3normal, NE_overlap,
     $         NE_corner_type, 24, 17, cpt2overlap_and_cpt3not   , no_overlap,
     $         NW_corner_type,  3, 18, cpt1overlap_and_cpt4normal, no_overlap,
     $         N_edge_type   ,  5, 18, 6                         , no_overlap,
     $         NW_edge_type  ,  7, 18, no_overlap                , no_overlap,
     $         NE_edge_type  , 19, 18, no_overlap                , no_overlap,
     $         N_edge_type   , 21, 18, 22                        , no_overlap,
     $         NE_corner_type, 23, 18, cpt2overlap_and_cpt3normal, no_overlap,
     $         NW_corner_type, 7 , 20, no_overlap                , no_overlap,
     $         N_edge_type   , 9 , 20, 10                        , no_overlap,
     $         NW_edge_type  , 11, 20, no_overlap                , N_overlap,
     $         NE_edge_type  , 15, 20, no_overlap                , N_overlap,
     $         N_edge_type   , 17, 20, 18                        , no_overlap,
     $         NE_corner_type, 19, 20, no_overlap                , no_overlap,
     $         NW_corner_type, 11, 21, no_overlap                , no_overlap,
     $         N_edge_type   , 13, 21, 14                        , no_overlap,
     $         NE_corner_type, 15, 21, no_overlap                , no_overlap
     $         /),
     $         (/5,60/))

          !output
          call bf_layer_bc_sections_used%ini()
          do j=j_min,j_max
             do i=i_min,i_max
                if(grdpts_id(i,j).eq.bc_interior_pt) then
                   call bf_layer_bc_sections_used%analyse_grdpt(
     $                  bf_alignment,
     $                  i,j,grdpts_id,
     $                  ierror)
                end if
             end do
          end do

          call bf_layer_bc_sections_used%finalize_bc_sections(
     $         x_borders,
     $         y_borders,
     $         sorted_bc_sections)

          test_validated = is_int_matrix_validated(
     $         sorted_bc_sections,
     $         test_bc_sections,
     $         detailled)

        end function test_finalize_bc_sections


        subroutine check_inputs()

          implicit none

          if((nx.le.10).or.
     $       (ny.le.10)) then

             print '(''the test requires: '')'
             print '(''nx>10: '',L1)', (nx.gt.10)
             print '(''ny>10: '',L1)', (ny.gt.10)
             stop ''

          end if

        end subroutine check_inputs

      end program test_bf_layer_bc_sections

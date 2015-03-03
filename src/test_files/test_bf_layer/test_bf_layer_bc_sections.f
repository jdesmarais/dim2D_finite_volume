      program test_bf_layer_bc_sections

        use bf_layer_bc_sections_class, only :
     $       bf_layer_bc_sections

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only : 
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
     $       
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt,
     $       no_pt


        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        test_loc = test_ini(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_ini: '',L1)', test_loc
        print '()'        

        test_loc = test_deallocate_tables(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_deallocate_tables: '',L1)', test_loc
        print '()'

c$$$        test_loc = test_add_temporary_bc_section()
c$$$        call test_add_final_bc_section()
c$$$        call test_add_bc_sections()
c$$$
c$$$        !test remove_from_bc_sections_temp()
c$$$        call test_remove_bc_section()
c$$$
c$$$        !test deallocate_tables()
c$$$        call bc_sections%deallocate_tables()
c$$$
c$$$
c$$$        !test get_bc_section()
c$$$        call test_get_bc_section(.false.)
c$$$
c$$$        !test analyse_grdpt_with_bc_section()
c$$$        call test_analyse_grdpt_with_bc_section(.false.)
c$$$
c$$$        !test analyse_grdpt()
c$$$        call test_analyse_grdpt(.false.)
c$$$
c$$$        !test add_overlap_between_corners_and_anti_corners()
c$$$        call test_add_overlap_between_corners_and_anti_corners(.true.)
c$$$        
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


c$$$        function get_test_bc_section(k)
c$$$     $       result(bc_section)
c$$$
c$$$          implicit none
c$$$
c$$$          integer, intent(in)   :: k
c$$$          integer, dimension(5) :: bc_section
c$$$
c$$$          integer, dimension(12) :: proc_type
c$$$
c$$$          proc_type = [
c$$$     $         N_edge_type,
c$$$     $         S_edge_type,
c$$$     $         E_edge_type,
c$$$     $         W_edge_type,
c$$$     $         NE_edge_type,
c$$$     $         NW_edge_type,
c$$$     $         SE_edge_type,
c$$$     $         SW_edge_type,
c$$$     $         NE_corner_type,
c$$$     $         NW_corner_type,
c$$$     $         SE_corner_type,
c$$$     $         SW_corner_type]
c$$$
c$$$          bc_section = [
c$$$     $         proc_type(k),
c$$$     $         k,
c$$$     $         k,
c$$$     $         k+1,
c$$$     $         mod(k,3)]
c$$$
c$$$        end function get_test_bc_section
c$$$
c$$$      
c$$$        subroutine test_add_temporary_bc_section()
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
c$$$          print '(''test add_to_temporary_bc_sections()'')'
c$$$
c$$$          !initialization
c$$$          call bc_sections%ini()
c$$$          
c$$$          do k=1,9
c$$$             new_bc_section = get_test_bc_section(k)
c$$$             call bc_sections%add_to_temporary_bc_sections(new_bc_section)
c$$$          end do
c$$$
c$$$          !test
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
c$$$          allocate(test_bc_sections_buffer(5,3))
c$$$          test_bc_sections_buffer = reshape(
c$$$     $         (/
c$$$     $         SE_edge_type,7,7,0,1,
c$$$     $         SW_edge_type,8,8,0,2,
c$$$     $         NE_corner_type,9,9,0,0
c$$$     $         /),
c$$$     $         (/5,3/))
c$$$
c$$$          same = compare_bc_sections_temporary(
c$$$     $         bc_sections,
c$$$     $         test_bc_sections_temp,
c$$$     $         test_bc_sections_buffer)
c$$$
c$$$          deallocate(test_bc_sections_buffer)
c$$$
c$$$          print '(''add_temporary_bc_sections: '',L1)', same
c$$$          print '()'
c$$$
c$$$        end subroutine test_add_temporary_bc_section
c$$$
c$$$
c$$$        subroutine test_add_final_bc_section()
c$$$
c$$$          implicit none
c$$$
c$$$          type(bf_layer_bc_sections)           :: bc_sections
c$$$          integer                              :: k
c$$$          integer                              :: nb_ele_final
c$$$          integer, dimension(:,:), allocatable :: bc_sections_final
c$$$          integer, dimension(5,9)              :: test_bc_sections_final
c$$$          logical                              :: same
c$$$
c$$$          integer, dimension(5) :: new_bc_section
c$$$
c$$$          print '(''test add_to_final_bc_sections()'')'
c$$$
c$$$          !initialization
c$$$          call bc_sections%ini()
c$$$          
c$$$          do k=1,9
c$$$             new_bc_section = get_test_bc_section(k)
c$$$             call bc_sections%add_to_final_bc_sections(new_bc_section)
c$$$          end do
c$$$
c$$$          !test
c$$$          test_bc_sections_final = reshape(
c$$$     $         (/
c$$$     $         N_edge_type,1,1,2,0,
c$$$     $         S_edge_type,2,2,3,1,
c$$$     $         E_edge_type,3,3,4,2,
c$$$     $         W_edge_type,4,4,5,0,
c$$$     $         NE_edge_type,5,5,0,2,
c$$$     $         NW_edge_type,6,6,0,0,
c$$$     $         SE_edge_type,7,7,0,1,
c$$$     $         SW_edge_type,8,8,0,2,
c$$$     $         NE_corner_type,9,9,0,0
c$$$     $         /),
c$$$     $         (/5,9/))
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
c$$$             deallocate(bc_sections_final)
c$$$
c$$$          end if
c$$$
c$$$          print '(''add_final_bc_sections: '',L1)', same
c$$$          print '()'
c$$$
c$$$        end subroutine test_add_final_bc_section
c$$$
c$$$
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
c$$$
c$$$
c$$$      subroutine make_test_get_bc_section(
c$$$     $       test_id,
c$$$     $       grdpts_id,
c$$$     $       test_i,
c$$$     $       test_j,
c$$$     $       test_bc_section)
c$$$        
c$$$          implicit none
c$$$
c$$$          integer                             , intent(in)  :: test_id
c$$$          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
c$$$          integer                             , intent(out) :: test_i
c$$$          integer                             , intent(out) :: test_j
c$$$          integer, dimension(5)               , intent(out) :: test_bc_section
c$$$
c$$$          
c$$$          allocate(grdpts_id(3,3))
c$$$          test_i = 2
c$$$          test_j = 2
c$$$
c$$$          select case(test_id)
c$$$            case(1)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,interior_pt,
c$$$     $              interior_pt,interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = SW_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case(2)
c$$$               grdpts_id =  reshape(
c$$$     $              (/
c$$$     $              interior_pt,bc_interior_pt,bc_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_pt,
c$$$     $              bc_pt,bc_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = NE_corner_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$
c$$$            case(3)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = NW_corner_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 2
c$$$
c$$$            case(4)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              interior_pt,bc_interior_pt,bc_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              interior_pt,interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = SE_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case(5)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              bc_interior_pt,bc_interior_pt,bc_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = E_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$
c$$$            case(6)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_interior_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = NW_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case(7)
c$$$               grdpts_id = reshape( (/
c$$$     $              interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = NE_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case(8)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = W_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$               
c$$$            case(9)
c$$$               grdpts_id = reshape( (/
c$$$     $              interior_pt,bc_interior_pt,bc_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = E_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$
c$$$            case(10)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = SW_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case(11)
c$$$               grdpts_id = reshape( (/
c$$$     $              interior_pt,bc_interior_pt,bc_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = SE_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case(12)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = W_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$
c$$$            case(13)
c$$$               grdpts_id = reshape( (/
c$$$     $              interior_pt,bc_interior_pt,bc_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = E_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$
c$$$            case(14)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = W_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$
c$$$            case(15)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_pt,bc_pt,bc_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = SE_corner_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 1
c$$$
c$$$            case(16)
c$$$               grdpts_id = reshape( (/
c$$$     $              interior_pt,interior_pt,interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = NW_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case(17)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_interior_pt,bc_pt,bc_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              interior_pt,interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = SE_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case(18)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_interior_pt,interior_pt,interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = N_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$
c$$$            case(19)
c$$$               grdpts_id = reshape( (/
c$$$     $              interior_pt,interior_pt,bc_interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = N_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$
c$$$            case(20)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_pt,bc_pt,bc_interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              interior_pt,interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = SW_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case(21)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_pt,bc_pt,bc_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_interior_pt,interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = S_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$
c$$$            case(22)
c$$$               grdpts_id = reshape( (/
c$$$     $              interior_pt,interior_pt,interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_interior_pt,bc_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = NE_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case(23)
c$$$               grdpts_id = reshape( (/
c$$$     $              interior_pt,interior_pt,interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = NW_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case(24)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_pt,bc_pt,bc_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              interior_pt,interior_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = S_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$
c$$$            case(25)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_pt,bc_pt,bc_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              interior_pt,interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = S_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$
c$$$            case(26)
c$$$               grdpts_id = reshape( (/
c$$$     $              interior_pt,interior_pt,interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = N_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 2
c$$$
c$$$            case(27)
c$$$               grdpts_id = reshape( (/
c$$$     $              bc_pt,bc_pt,bc_pt,
c$$$     $              bc_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = SW_corner_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$
c$$$            case(28)
c$$$               grdpts_id = reshape( (/
c$$$     $              interior_pt,interior_pt,interior_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$               test_bc_section(1) = NE_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(5) = 1
c$$$
c$$$            case default
c$$$               print '(''test_bf_layer_bc_procedure'')'
c$$$               print '(''make_test_bf_layer_bc_procedure'')'
c$$$               print '(''test case not implemented: '', I2)', test_id
c$$$               stop ''
c$$$
c$$$          end select
c$$$
c$$$        end subroutine make_test_get_bc_section
c$$$
c$$$
c$$$        subroutine make_test_analyse_grdpt_with_bc_section(
c$$$     $       test_id,
c$$$     $       grdpts_id,
c$$$     $       test_i,
c$$$     $       test_j,
c$$$     $       test_bc_section,
c$$$     $       test_compatible,
c$$$     $       test_remove_ele)
c$$$        
c$$$          implicit none
c$$$
c$$$          integer                             , intent(in)  :: test_id
c$$$          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
c$$$          integer                             , intent(out) :: test_i
c$$$          integer                             , intent(out) :: test_j
c$$$          integer, dimension(5)               , intent(out) :: test_bc_section
c$$$          logical                             , intent(out) :: test_compatible
c$$$          logical                             , intent(out) :: test_remove_ele
c$$$
c$$$          
c$$$          allocate(grdpts_id(3,3))
c$$$
c$$$          select case(test_id)
c$$$
c$$$            !test W_edge: compatible
c$$$
c$$$            !  -------    
c$$$            ! | 2 1 0 |
c$$$            ! | 2 1*0 |
c$$$            ! | 2 1 0 |
c$$$            !  -------
c$$$            case(1)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt,
c$$$     $              bc_pt,bc_interior_pt,interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 2
c$$$               test_j             = 2
c$$$
c$$$               test_bc_section(1) = W_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(4) = 2
c$$$               test_bc_section(5) = 0
c$$$               test_compatible    = .true.
c$$$               test_remove_ele    = .false.
c$$$
c$$$            !test E_edge: compatible
c$$$
c$$$            !  -------    
c$$$            ! | 1 2   |
c$$$            ! | 1*2   |
c$$$            ! | 1 2   |
c$$$            !  -------
c$$$            case(2)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              bc_interior_pt,bc_pt,no_pt,
c$$$     $              bc_interior_pt,bc_pt,no_pt,
c$$$     $              bc_interior_pt,bc_pt,no_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 1
c$$$               test_j             = 2
c$$$
c$$$               test_bc_section(1) = E_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(4) = 1
c$$$               test_bc_section(5) = 0
c$$$               test_compatible    = .true.
c$$$               test_remove_ele    = .false.
c$$$
c$$$            !test N_edge: compatible
c$$$
c$$$            !  -------    
c$$$            ! |       |
c$$$            ! | 2 2 2 |
c$$$            ! | 1 1*1 |
c$$$            !  -------
c$$$            case(3)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_pt,bc_pt,
c$$$     $              no_pt,no_pt,no_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 2
c$$$               test_j             = 1
c$$$
c$$$               test_bc_section(1) = N_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(4) = 1
c$$$               test_bc_section(5) = 0
c$$$               test_compatible    = .true.
c$$$               test_remove_ele    = .false.
c$$$
c$$$            !test S_edge: compatible
c$$$
c$$$            !  -------    
c$$$            ! | 1 1*1 |
c$$$            ! | 2 2 2 |
c$$$            ! |       |
c$$$            !  -------
c$$$            case(4)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              no_pt,no_pt,no_pt,
c$$$     $              bc_pt,bc_pt,bc_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 2
c$$$               test_j             = 3
c$$$
c$$$               test_bc_section(1) = S_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(4) = 3
c$$$               test_bc_section(5) = 0
c$$$               test_compatible    = .true.
c$$$               test_remove_ele    = .false.
c$$$
c$$$            !test W_edge: not compatible
c$$$
c$$$            !  -------
c$$$            ! | 1 1 1 |
c$$$            ! | 2 2 1*|
c$$$            ! |   2 1 |
c$$$            !  -------
c$$$            case(5)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              no_pt,bc_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_pt,bc_interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 3
c$$$               test_j             = 2
c$$$
c$$$               test_bc_section(1) = W_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(4) = 3
c$$$               test_bc_section(5) = 0
c$$$               test_compatible    = .false.
c$$$               test_remove_ele    = .false.
c$$$
c$$$            !test E_edge: not compatible
c$$$
c$$$            !  -------    
c$$$            ! | 1 1 1 |
c$$$            ! | 1*2 2 |
c$$$            ! | 1 2   |
c$$$            !  -------
c$$$            case(6)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              bc_interior_pt,bc_pt,no_pt,
c$$$     $              bc_interior_pt,bc_pt,bc_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 1
c$$$               test_j             = 2
c$$$
c$$$               test_bc_section(1) = E_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(4) = 1
c$$$               test_bc_section(5) = 0
c$$$               test_compatible    = .false.
c$$$               test_remove_ele    = .false.
c$$$
c$$$            !test N_edge: not compatible
c$$$
c$$$            !  -------    
c$$$            ! |   2 1 |
c$$$            ! | 2 2 1 |
c$$$            ! | 1 1*1 |
c$$$            !  -------
c$$$            case(7)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_pt,bc_interior_pt,
c$$$     $              no_pt,bc_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 2
c$$$               test_j             = 1
c$$$
c$$$               test_bc_section(1) = N_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(4) = 1
c$$$               test_bc_section(5) = 0
c$$$               test_compatible    = .false.
c$$$               test_remove_ele    = .false.
c$$$
c$$$            !test S_edge: not compatible
c$$$
c$$$            !  -------    
c$$$            ! | 1 1*1 |
c$$$            ! | 2 2 1 |
c$$$            ! |   2 1 |
c$$$            !  -------
c$$$            case(8)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              no_pt,bc_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_pt,bc_interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 2
c$$$               test_j             = 3
c$$$
c$$$               test_bc_section(1) = S_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(4) = 3
c$$$               test_bc_section(5) = 0
c$$$               test_compatible    = .false.
c$$$               test_remove_ele    = .false.
c$$$
c$$$            !test SW_edge: compatible
c$$$
c$$$            !  -------    
c$$$            ! | 1 1 1 |
c$$$            ! | 2 2 1*|
c$$$            ! |   2 1 |
c$$$            !  -------
c$$$            case(9)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              no_pt,bc_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_pt,bc_interior_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 3
c$$$               test_j             = 2
c$$$
c$$$               test_bc_section(1) = SW_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 0
c$$$               test_bc_section(5) = 1
c$$$               test_compatible    = .true.
c$$$               test_remove_ele    = .false.
c$$$
c$$$            !test SE_edge: compatible
c$$$
c$$$            !  -------    
c$$$            ! | 1 1 1 |
c$$$            ! | 1*2 2 |
c$$$            ! | 1 2   |
c$$$            !  -------
c$$$            case(10)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              bc_interior_pt,bc_pt,no_pt,
c$$$     $              bc_interior_pt,bc_pt,bc_pt,
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 1
c$$$               test_j             = 2
c$$$
c$$$               test_bc_section(1) = SE_edge_type
c$$$               test_bc_section(2) = 1
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 0
c$$$               test_bc_section(5) = 2
c$$$               test_compatible    = .true.
c$$$               test_remove_ele    = .true.
c$$$
c$$$            !test NW_edge: compatible
c$$$
c$$$            !  -------    
c$$$            ! |   2 1 |
c$$$            ! | 2 2 1 |
c$$$            ! | 1 1*1 |
c$$$            !  -------
c$$$            case(11)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              bc_pt,bc_pt,bc_interior_pt,
c$$$     $              no_pt,bc_pt,bc_interior_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 2
c$$$               test_j             = 1
c$$$
c$$$               test_bc_section(1) = NW_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 1
c$$$               test_bc_section(4) = 0
c$$$               test_bc_section(5) = 2
c$$$               test_compatible    = .true.
c$$$               test_remove_ele    = .true.
c$$$
c$$$            !test NE_edge: compatible
c$$$
c$$$            !  -------    
c$$$            ! | 0 1 2 |
c$$$            ! | 0 1*1 |
c$$$            ! | 0 0 0 |
c$$$            !  -------
c$$$            case(12)
c$$$               grdpts_id = reshape(
c$$$     $              (/
c$$$     $              interior_pt,interior_pt,interior_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_interior_pt,
c$$$     $              interior_pt,bc_interior_pt,bc_pt
c$$$     $              /),
c$$$     $              (/3,3/))
c$$$
c$$$               test_i             = 2
c$$$               test_j             = 2
c$$$
c$$$               test_bc_section(1) = NE_edge_type
c$$$               test_bc_section(2) = 2
c$$$               test_bc_section(3) = 2
c$$$               test_bc_section(4) = 0
c$$$               test_bc_section(5) = 2
c$$$               test_compatible    = .true.
c$$$               test_remove_ele    = .true.
c$$$
c$$$            end select
c$$$
c$$$        end subroutine make_test_analyse_grdpt_with_bc_section
c$$$
c$$$
c$$$        subroutine test_get_bc_section(detailled)
c$$$
c$$$          implicit none
c$$$
c$$$          logical                 , intent(in) :: detailled
c$$$
c$$$          type(bf_layer_bc_sections)           :: bc_sections
c$$$
c$$$          integer, dimension(5)                :: new_bc_section
c$$$          integer                              :: k
c$$$
c$$$          integer, dimension(:,:), allocatable :: grdpts_id
c$$$          integer                              :: test_i
c$$$          integer                              :: test_j
c$$$          integer, dimension(5)                :: test_bc_section
c$$$          logical                              :: test_validated
c$$$
c$$$          logical                              :: test_global
c$$$
c$$$          logical                              :: ierror
c$$$
c$$$
c$$$          test_global = .true.
c$$$          
c$$$          print '(''test get_bc_section()'')'
c$$$
c$$$          do k=1,26
c$$$
c$$$             call make_test_get_bc_section(
c$$$     $            k,
c$$$     $            grdpts_id,
c$$$     $            test_i,
c$$$     $            test_j,
c$$$     $            test_bc_section)
c$$$
c$$$             new_bc_section = bc_sections%get_bc_section(
c$$$     $            test_i,
c$$$     $            test_j,
c$$$     $            grdpts_id,
c$$$     $            ierror)
c$$$
c$$$
c$$$             test_validated = test_bc_section(1).eq.new_bc_section(1)
c$$$
c$$$             select case(test_bc_section(1))
c$$$
c$$$               !if this is an edge procedure, only the 2:4
c$$$               !elements should be tested
c$$$               case(N_edge_type,S_edge_type,E_edge_type,W_edge_type)
c$$$                  
c$$$                  test_validated = test_validated.and.(
c$$$     $                 test_bc_section(2).eq.new_bc_section(2))
c$$$                  test_validated = test_validated.and.(
c$$$     $                 test_bc_section(3).eq.new_bc_section(3))
c$$$                  test_validated = test_validated.and.(
c$$$     $                 test_bc_section(4).eq.new_bc_section(4))
c$$$
c$$$                  if(detailled) then
c$$$                     print '(''test '',I2,'':'',L1)', k, test_validated
c$$$                  end if
c$$$
c$$$                  if(.not.test_validated) then
c$$$                     print '(''procedure_type: '',I2, 2X,I2)', test_bc_section(1), new_bc_section(1)
c$$$                     print '(''edge_min      : '',I2, 2X,I2)', test_bc_section(2), new_bc_section(2)
c$$$                     print '(''edge_max      : '',I2, 2X,I2)', test_bc_section(3), new_bc_section(3)
c$$$                     print '(''coord         : '',I2, 2X,I2)', test_bc_section(4), new_bc_section(4)
c$$$                  end if
c$$$
c$$$               !if this is a corner procedure, only the 2:3
c$$$               !elements should be tested
c$$$               case(NE_corner_type, NW_corner_type, SE_corner_type, SW_corner_type)
c$$$
c$$$                  test_validated = test_validated.and.(
c$$$     $                 test_bc_section(2).eq.new_bc_section(2))
c$$$                  test_validated = test_validated.and.(
c$$$     $                 test_bc_section(3).eq.new_bc_section(3))
c$$$
c$$$                  if(detailled) then
c$$$                     print '(''test '',I2,'':'',L1)', k, test_validated
c$$$                  end if
c$$$
c$$$                  if(.not.test_validated) then
c$$$                     print '(''procedure_type: '',I2, 2X,I2)', test_bc_section(1), new_bc_section(1)
c$$$                     print '(''i_min         : '',I2, 2X,I2)', test_bc_section(2), new_bc_section(2)
c$$$                     print '(''j_min         : '',I2, 2X,I2)', test_bc_section(3), new_bc_section(3)
c$$$                  end if
c$$$
c$$$               !if this is a special edge procedure, only the 2:3 and 5
c$$$               !elements should be tested
c$$$               case(NE_edge_type,NW_edge_type,SE_edge_type,SW_edge_type)
c$$$
c$$$                  test_validated = test_validated.and.(
c$$$     $                 test_bc_section(2).eq.new_bc_section(2))
c$$$                  test_validated = test_validated.and.(
c$$$     $                 test_bc_section(3).eq.new_bc_section(3))
c$$$                  test_validated = test_validated.and.(
c$$$     $                 test_bc_section(5).eq.new_bc_section(5))
c$$$
c$$$                  if(detailled) then
c$$$                     print '(''test '',I2,'':'',L1)', k, test_validated
c$$$                  end if
c$$$
c$$$                  if(.not.test_validated) then
c$$$                     print '(''procedure_type: '',I2, 2X,I2)', test_bc_section(1), new_bc_section(1)
c$$$                     print '(''i_min         : '',I2, 2X,I2)', test_bc_section(2), new_bc_section(2)
c$$$                     print '(''j_min         : '',I2, 2X,I2)', test_bc_section(3), new_bc_section(3)
c$$$                     print '(''match_nb      : '',I2, 2X,I2)', test_bc_section(5), new_bc_section(5)
c$$$                  end if
c$$$
c$$$               case default
c$$$                  
c$$$                  print '(''test_bf_layer_bc_sections'')'
c$$$                  print '(''test_get_bc-sections()'')'
c$$$                  print '(''case '', I2, ''not recognized'')', k
c$$$                  stop ''
c$$$
c$$$             end select
c$$$
c$$$             deallocate(grdpts_id)
c$$$
c$$$             test_global = test_global.and.test_validated
c$$$             
c$$$          end do
c$$$
c$$$          print '(''test_validated: '',L1)', test_global
c$$$          print '()'
c$$$
c$$$        end subroutine test_get_bc_section
c$$$
c$$$
c$$$        subroutine test_analyse_grdpt_with_bc_section(detailled)
c$$$
c$$$          implicit none
c$$$
c$$$          logical                 , intent(in) :: detailled
c$$$
c$$$          type(bf_layer_bc_sections)           :: bc_sections
c$$$
c$$$          integer                              :: k
c$$$          integer, dimension(:,:), allocatable :: grdpts_id
c$$$          logical                              :: compatible
c$$$          logical                              :: remove_ele
c$$$
c$$$          integer                              :: test_i
c$$$          integer                              :: test_j
c$$$          integer, dimension(5)                :: test_bc_section
c$$$          logical                              :: test_validated
c$$$          logical                              :: test_compatible
c$$$          logical                              :: test_remove_ele
c$$$
c$$$          logical :: test_global
c$$$
c$$$          test_global = .true.          
c$$$
c$$$          
c$$$          print '(''test analyse_grdpt_with_bc_section()'')'
c$$$          do k=1,12
c$$$
c$$$             call make_test_analyse_grdpt_with_bc_section(
c$$$     $            k,
c$$$     $            grdpts_id,
c$$$     $            test_i,
c$$$     $            test_j,
c$$$     $            test_bc_section,
c$$$     $            test_compatible,
c$$$     $            test_remove_ele)
c$$$
c$$$             compatible = bc_sections%analyse_grdpt_with_bc_section(
c$$$     $            test_i,test_j,grdpts_id,test_bc_section,remove_ele)
c$$$
c$$$             test_validated = compatible.eqv.test_compatible
c$$$             test_validated = test_validated.and.(remove_ele.eqv.test_remove_ele)
c$$$
c$$$             if(detailled) then
c$$$                print '(''test '',I2,'':'',L1)', k, test_validated
c$$$             end if
c$$$
c$$$             if(.not.test_validated) then
c$$$                print '(''  compatible    : '',L1, 2X,L1)', test_compatible, compatible
c$$$                print '(''  remove_ele    : '',L1, 2X,L1)', test_remove_ele, remove_ele
c$$$             end if
c$$$
c$$$             deallocate(grdpts_id)
c$$$
c$$$             test_global = test_global.and.test_validated
c$$$
c$$$          end do
c$$$          print '(''test_validated: '',L1)', test_global
c$$$          print '()'
c$$$
c$$$        end subroutine test_analyse_grdpt_with_bc_section
c$$$
c$$$
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

      end program test_bf_layer_bc_sections

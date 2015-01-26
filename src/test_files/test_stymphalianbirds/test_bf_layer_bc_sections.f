      program test_bf_layer_bc_sections

        use bf_layer_bc_procedure_module, only : 
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
     $       SW_edge_type

        use bf_layer_bc_sections_class, only :
     $       bf_layer_bc_sections,
     $       no_overlap,
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap,
     $       NE_overlap,
     $       NW_overlap,
     $       SE_overlap,
     $       SW_overlap

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt,
     $       no_pt


        implicit none

        type(bf_layer_bc_sections)           :: bc_sections

        !test ini()
        call test_ini()

        !test add_to_temporary_bc_sections() and
        !add_to_final_bc_sections()
        call test_add_temporary_bc_section()
        call test_add_final_bc_section()
        call test_add_bc_sections()

        !test remove_from_bc_sections_temp()
        call test_remove_bc_section()

        !test deallocate_tables()
        call bc_sections%deallocate_tables()


        !test get_bc_section()
        call test_get_bc_section(.false.)

        !test analyse_grdpt_with_bc_section()
        call test_analyse_grdpt_with_bc_section(.false.)

        !test analyse_grdpt()
        call test_analyse_grdpt(.false.)

        !test add_overlap_between_corners_and_anti_corners()
        call test_add_overlap_between_corners_and_anti_corners(.true.)
        
        contains

        
        subroutine test_ini()

          implicit none

          type(bf_layer_bc_sections) :: bc_sections
          logical                    :: test_validated

          print '(''test ini()'')'
          call bc_sections%ini()

          test_validated = bc_sections%get_nb_ele_temp().eq.0
          test_validated = test_validated.and.(bc_sections%get_nb_ele_final().eq.0)

          print '(''test_ini: '',L1)', test_validated
          print '()'

        end subroutine test_ini


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

      
        subroutine test_add_temporary_bc_section()

          implicit none

          type(bf_layer_bc_sections)           :: bc_sections
          integer                              :: k
          integer, dimension(5)                :: new_bc_section
          integer, dimension(5,6)              :: test_bc_sections_temp
          integer, dimension(:,:), allocatable :: test_bc_sections_buffer
          logical                              :: same


          print '(''test add_to_temporary_bc_sections()'')'

          !initialization
          call bc_sections%ini()
          
          do k=1,9
             new_bc_section = get_test_bc_section(k)
             call bc_sections%add_to_temporary_bc_sections(new_bc_section)
          end do

          !test
          test_bc_sections_temp = reshape(
     $         (/
     $         N_edge_type,1,1,2,0,
     $         S_edge_type,2,2,3,1,
     $         E_edge_type,3,3,4,2,
     $         W_edge_type,4,4,5,0,
     $         NE_edge_type,5,5,0,2,
     $         NW_edge_type,6,6,0,0
     $         /),
     $         (/5,6/))
          allocate(test_bc_sections_buffer(5,3))
          test_bc_sections_buffer = reshape(
     $         (/
     $         SE_edge_type,7,7,0,1,
     $         SW_edge_type,8,8,0,2,
     $         NE_corner_type,9,9,0,0
     $         /),
     $         (/5,3/))

          same = compare_bc_sections_temporary(
     $         bc_sections,
     $         test_bc_sections_temp,
     $         test_bc_sections_buffer)

          deallocate(test_bc_sections_buffer)

          print '(''add_temporary_bc_sections: '',L1)', same
          print '()'

        end subroutine test_add_temporary_bc_section


        subroutine test_add_final_bc_section()

          implicit none

          type(bf_layer_bc_sections)           :: bc_sections
          integer                              :: k
          integer                              :: nb_ele_final
          integer, dimension(:,:), allocatable :: bc_sections_final
          integer, dimension(5,9)              :: test_bc_sections_final
          logical                              :: same

          integer, dimension(5) :: new_bc_section

          print '(''test add_to_final_bc_sections()'')'

          !initialization
          call bc_sections%ini()
          
          do k=1,9
             new_bc_section = get_test_bc_section(k)
             call bc_sections%add_to_final_bc_sections(new_bc_section)
          end do

          !test
          test_bc_sections_final = reshape(
     $         (/
     $         N_edge_type,1,1,2,0,
     $         S_edge_type,2,2,3,1,
     $         E_edge_type,3,3,4,2,
     $         W_edge_type,4,4,5,0,
     $         NE_edge_type,5,5,0,2,
     $         NW_edge_type,6,6,0,0,
     $         SE_edge_type,7,7,0,1,
     $         SW_edge_type,8,8,0,2,
     $         NE_corner_type,9,9,0,0
     $         /),
     $         (/5,9/))

          same = .true.

          nb_ele_final = bc_sections%get_nb_ele_final()
          call bc_sections%get_bc_sections_final(bc_sections_final)

          if(allocated(bc_sections_final)) then

             do k=1, min(size(bc_sections_final,2),nb_ele_final)

                same = compare_bc_procedure(
     $               bc_sections_final(:,k),
     $               test_bc_sections_final(:,k))

                if(.not.same) then
                   print '(''   k                    : '', I3)' , k
                   print '(''   bc_section_final     : '', 5I3)', bc_sections_final(:,k)
                   print '(''   test_bc_section_final: '', 5I3)', test_bc_sections_final(:,k)
                   exit
                end if

             end do

             deallocate(bc_sections_final)

          end if

          print '(''add_final_bc_sections: '',L1)', same
          print '()'

        end subroutine test_add_final_bc_section


        subroutine test_add_bc_sections()

          implicit none

          type(bf_layer_bc_sections) :: bc_sections
          integer                    :: k
          integer, dimension(5)      :: new_bc_section

          integer, dimension(5,6)              :: test_bc_sections_temp
          integer, dimension(5,2)              :: test_bc_sections_buffer
          integer, dimension(5,1)              :: test_bc_sections_final

          logical :: same

          print '(''test add_to_bc_sections()'')'
          call bc_sections%ini()
          do k=1,9
             new_bc_section = get_test_bc_section(k)
             call bc_sections%add_to_bc_sections(new_bc_section)
          end do

          !bc_sections_temp
          test_bc_sections_temp = reshape(
     $         (/
     $         N_edge_type,1,1,2,0,
     $         S_edge_type,2,2,3,1,
     $         E_edge_type,3,3,4,2,
     $         W_edge_type,4,4,5,0,
     $         NE_edge_type,5,5,0,2,
     $         NW_edge_type,6,6,0,0
     $         /),
     $         (/5,6/))

          !bc_sections_buffer
          test_bc_sections_buffer = reshape(
     $         (/
     $         SE_edge_type,7,7,0,1,
     $         SW_edge_type,8,8,0,2
     $         /),
     $         (/5,2/))
          
          !bc_sections_final
          test_bc_sections_final = reshape(
     $         (/
     $         NE_corner_type,9,9,0,0
     $         /),
     $         (/5,1/))


          !same bc_sections_temp and bc_sections_buffer ?
          same = compare_bc_sections_temporary(
     $         bc_sections,
     $         test_bc_sections_temp,
     $         test_bc_sections_buffer)
          
          print '(''bc_sections_temporary: '',L1)', same


          !same bc_sections_final ?
          same = compare_bc_sections_final(
     $         bc_sections,
     $         test_bc_sections_final)          

          print '(''bc_sections_final    : '',L1)', same
          print '()'
          
        end subroutine test_add_bc_sections


        function compare_bc_sections_final(
     $     bc_sections,
     $     test_bc_sections_final)
     $     result(same)

          implicit none

          type(bf_layer_bc_sections), intent(in) :: bc_sections
          integer, dimension(:,:)   , intent(in) :: test_bc_sections_final
          logical                                :: same

          integer                                :: nb_ele_final
          integer, dimension(:,:), allocatable   :: bc_sections_final
          integer                                :: k


          same = .true.

          nb_ele_final = bc_sections%get_nb_ele_final()
          call bc_sections%get_bc_sections_final(bc_sections_final)

          if(allocated(bc_sections_final)) then

             do k=1, min(size(bc_sections_final,2),nb_ele_final)

                same = compare_bc_procedure(
     $               bc_sections_final(:,k),
     $               test_bc_sections_final(:,k))

                if(.not.same) then
                   print '(''   k                    : '', I3)' , k
                   print '(''   bc_section_final     : '', 5I3)', bc_sections_final(:,k)
                   print '(''   test_bc_section_final: '', 5I3)', test_bc_sections_final(:,k)
                   exit
                end if

             end do

          end if

        end function compare_bc_sections_final

      
        function compare_bc_sections_sorted(
     $     bc_sections_sorted,
     $     test_bc_sections_sorted)
     $     result(same)

          implicit none

          integer, dimension(:,:)   , intent(in) :: bc_sections_sorted
          integer, dimension(:,:)   , intent(in) :: test_bc_sections_sorted
          logical                                :: same

          integer                                :: k


          same = .true.

          do k=1, size(bc_sections_sorted,2)

             same = compare_bc_layer_sorted(
     $            bc_sections_sorted(:,k),
     $            test_bc_sections_sorted(:,k))

             if(.not.same) then
                print '(''   k                     : '',  I3)', k
                print '(''   bc_section_sorted     : '', 4I3)', bc_sections_sorted(:,k)
                print '(''   test_bc_section_sorted: '', 4I3)', test_bc_sections_sorted(:,k)
                exit
             end if
             
          end do

        end function compare_bc_sections_sorted


        function compare_bc_layer_sorted(bc_sorted,test_bc_sorted)
     $     result(same)

          implicit none

          integer, dimension(4), intent(in) :: bc_sorted
          integer, dimension(4), intent(in) :: test_bc_sorted
          logical                           :: same

          same = bc_sorted(1).eq.test_bc_sorted(1)          

          select case(bc_sorted(1))

            case(N_edge_type,S_edge_type,E_edge_type,W_edge_type)
               same = same.and.(bc_sorted(2).eq.test_bc_sorted(2))
               same = same.and.(bc_sorted(3).eq.test_bc_sorted(3))
               same = same.and.(bc_sorted(4).eq.test_bc_sorted(4))

            case(
     $              NE_corner_type,NW_corner_type,
     $              SE_corner_type,SW_corner_type,
     $              NE_edge_type,NW_edge_type,
     $              SE_edge_type,SW_edge_type
     $              )
               same = same.and.(bc_sorted(2).eq.test_bc_sorted(2))
               same = same.and.(bc_sorted(3).eq.test_bc_sorted(3))

            case default
               print '(''test_bf_layer_bc_sections'')'
               print '(''compare_bc_layer_sorted'')'
               print '(''case not recognized: '',I2)', bc_sorted(1)
               stop ''

          end select

        end function compare_bc_layer_sorted


        subroutine test_remove_bc_section()

          implicit none

          type(bf_layer_bc_sections)           :: bc_sections
          integer                              :: k
          integer, dimension(5)                :: new_bc_section
          integer, dimension(5,6)              :: test_bc_sections_temp
          integer, dimension(:,:), allocatable :: test_bc_sections_buffer
          logical                              :: same


          !initialization
          call bc_sections%ini()
          
          do k=1,9
             new_bc_section = get_test_bc_section(k)
             call bc_sections%add_to_temporary_bc_sections(new_bc_section)
          end do


          !test remove_from_bc_sections_temp()
          print '(''test remove_from_bc_sections_temp()'')'
          call bc_sections%remove_from_bc_sections_temp(1)
          test_bc_sections_temp = reshape(
     $         (/
     $         S_edge_type,2,2,3,0,
     $         E_edge_type,3,3,4,0,
     $         W_edge_type,4,4,5,0,
     $         NE_edge_type,5,5,0,2,
     $         NW_edge_type,6,6,0,0,
     $         SE_edge_type,7,7,0,1
     $         /),
     $         (/5,6/))
          allocate(test_bc_sections_buffer(5,2))
          test_bc_sections_buffer = reshape(
     $         (/
     $         SW_edge_type,8,8,0,2,
     $         NE_corner_type,9,9,0,0
     $         /),
     $         (/5,2/))

          same = compare_bc_sections_temporary(
     $         bc_sections,
     $         test_bc_sections_temp,
     $         test_bc_sections_buffer)

          deallocate(test_bc_sections_buffer)

          print '(''remove_from_bc_sections_temp(1): '',L1)', same


          call bc_sections%remove_from_bc_sections_temp(4)
          test_bc_sections_temp = reshape(
     $         (/
     $         S_edge_type,2,2,3,0,
     $         E_edge_type,3,3,4,0,
     $         W_edge_type,4,4,5,0,
     $         NW_edge_type,6,6,0,0,
     $         SE_edge_type,7,7,0,1,
     $         SW_edge_type,8,8,0,2
     $         /),
     $         (/5,6/))
          allocate(test_bc_sections_buffer(5,1))
          test_bc_sections_buffer = reshape(
     $         (/
     $         NE_corner_type,9,9,0,0
     $         /),
     $         (/5,1/))

          same = compare_bc_sections_temporary(
     $         bc_sections,
     $         test_bc_sections_temp,
     $         test_bc_sections_buffer)

          deallocate(test_bc_sections_buffer)

          print '(''remove_from_bc_sections_temp(4): '',L1)', same


          call bc_sections%remove_from_bc_sections_temp(6)
          test_bc_sections_temp = reshape(
     $         (/
     $         S_edge_type,2,2,3,0,
     $         E_edge_type,3,3,4,0,
     $         W_edge_type,4,4,5,0,
     $         NW_edge_type,6,6,0,0,
     $         SE_edge_type,7,7,0,1,
     $         NE_corner_type,9,9,0,0
     $         /),
     $         (/5,6/))

          same = compare_bc_sections_temporary(
     $         bc_sections,
     $         test_bc_sections_temp,
     $         test_bc_sections_buffer)

          print '(''remove_from_bc_sections_temp(6): '',L1)', same

          print '()'


          do k=10,12
             new_bc_section = get_test_bc_section(k)
             call bc_sections%add_to_temporary_bc_sections(new_bc_section)
          end do

          call bc_sections%remove_from_bc_sections_buffer(1)
          
          allocate(test_bc_sections_buffer(5,2))

          test_bc_sections_buffer = reshape(
     $         (/
     $         SE_corner_type,11,11,0,2,
     $         SW_corner_type,12,12,0,0
     $         /),
     $         (/5,2/))

          same = compare_bc_sections_temporary(
     $         bc_sections,
     $         test_bc_sections_temp,
     $         test_bc_sections_buffer)

          deallocate(test_bc_sections_buffer)

          print '(''remove_from_bc_sections_buffer(1): '',L1)', same


          call bc_sections%remove_from_bc_sections_buffer(2)

          allocate(test_bc_sections_buffer(5,1))

          test_bc_sections_buffer = reshape(
     $         (/
     $          SE_corner_type,11,11,0,2
     $         /),
     $         (/5,1/))

          same = compare_bc_sections_temporary(
     $         bc_sections,
     $         test_bc_sections_temp,
     $         test_bc_sections_buffer)

          deallocate(test_bc_sections_buffer)

          print '(''remove_from_bc_sections_temp(2): '',L1)', same

          print '()'

      end subroutine test_remove_bc_section


      function compare_bc_sections_temporary(
     $     bc_sections,
     $     test_bc_sections_temp,
     $     test_bc_sections_buffer)
     $     result(same)

        implicit none

        class(bf_layer_bc_sections), intent(in) :: bc_sections
        integer, dimension(:,:)    , intent(in) :: test_bc_sections_temp
        integer, dimension(:,:)    , intent(in) :: test_bc_sections_buffer
        logical                                 :: same


        integer :: k
        integer :: nb_ele_temp
        integer, dimension(:,:), allocatable :: bc_sections_temp
        integer, dimension(:,:), allocatable :: bc_sections_buffer


        nb_ele_temp = bc_sections%get_nb_ele_temp()
        call bc_sections%get_bc_sections_temp(bc_sections_temp)
        call bc_sections%get_bc_sections_buffer(bc_sections_buffer)


        same = .true.
        if(allocated(bc_sections_temp)) then

           do k=1, min(size(bc_sections_temp,2),nb_ele_temp)

              same = compare_bc_procedure(
     $             bc_sections_temp(:,k),
     $             test_bc_sections_temp(:,k))

              if(.not.same) then
                 print '(''   k                   : '', I3)' , k
                 print '(''   bc_section_temp     : '', 5I3)', bc_sections_temp(:,k)
                 print '(''   test_bc_section_temp: '', 5I3)', test_bc_sections_temp(:,k)
                 exit
              end if

           end do

           deallocate(bc_sections_temp)

        end if

        if(same) then
           if(allocated(bc_sections_buffer)) then
              
              do k=1, min(
     $             size(bc_sections_buffer,2),
     $             nb_ele_temp-size(bc_sections_temp,2))
                 
                 same = compare_bc_procedure(
     $                bc_sections_buffer(:,k),
     $                test_bc_sections_buffer(:,k))

                 if(.not.same) then
                    print '(''   k                   : '', I3)' , k
                    print '(''   bc_section_buff     : '', 5I3)', bc_sections_buffer(:,k)
                    print '(''   test_bc_section_buff: '', 5I3)', test_bc_sections_buffer(:,k)
                    exit
                 end if
                 
              end do

              deallocate(bc_sections_buffer)
           
           end if
        end if

      end function compare_bc_sections_temporary


      function compare_bc_procedure(bc_sections,test_bc_sections)
     $     result(same)

        implicit none

        integer, dimension(:), intent(in) :: bc_sections
        integer, dimension(:), intent(in) :: test_bc_sections
        logical                           :: same

        same = .true.

        same = bc_sections(1).eq.test_bc_sections(1)

        select case(bc_sections(1))
          case(N_edge_type,S_edge_type,E_edge_type,W_edge_type)

             same = same.and.(bc_sections(2).eq.test_bc_sections(2))
             same = same.and.(bc_sections(3).eq.test_bc_sections(3))
             same = same.and.(bc_sections(4).eq.test_bc_sections(4))

          case(NE_corner_type,NW_corner_type,SE_corner_type,SW_corner_type)

             same = same.and.(bc_sections(2).eq.test_bc_sections(2))
             same = same.and.(bc_sections(3).eq.test_bc_sections(3))

          case(NE_edge_type,NW_edge_type,SE_edge_type,SW_edge_type)

             same = same.and.(bc_sections(2).eq.test_bc_sections(2))
             same = same.and.(bc_sections(3).eq.test_bc_sections(3))
             same = same.and.(bc_sections(5).eq.test_bc_sections(5))

          case default
             print '(''test_bf_layer_bc_sections'')'
             print '(''compare_bc_procedure'')'
             print '(''case: '',I2)', bc_sections(1)
             stop 'case not recognized'
        end select

      end function compare_bc_procedure


      subroutine make_test_get_bc_section(
     $       test_id,
     $       grdpts_id,
     $       test_i,
     $       test_j,
     $       test_bc_section)
        
          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          integer, dimension(5)               , intent(out) :: test_bc_section

          
          allocate(grdpts_id(3,3))
          test_i = 2
          test_j = 2

          select case(test_id)
            case(1)
               grdpts_id = reshape(
     $              (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SW_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(2)
               grdpts_id =  reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NE_corner_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2

            case(3)
               grdpts_id = reshape(
     $              (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NW_corner_type
               test_bc_section(2) = 1
               test_bc_section(3) = 2

            case(4)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SE_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(5)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = E_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(6)
               grdpts_id = reshape( (/
     $              bc_interior_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NW_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(7)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NE_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(8)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = W_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2
               
            case(9)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = E_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(10)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SW_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case(11)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SE_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case(12)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = W_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(13)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = E_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(14)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = W_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(15)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SE_corner_type
               test_bc_section(2) = 2
               test_bc_section(3) = 1

            case(16)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NW_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case(17)
               grdpts_id = reshape( (/
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SE_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(18)
               grdpts_id = reshape( (/
     $              bc_interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = N_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(19)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = N_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(20)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SW_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(21)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = S_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(22)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NE_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case(23)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NW_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case(24)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = S_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(25)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = S_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(26)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = N_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(27)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SW_corner_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1

            case(28)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NE_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case default
               print '(''test_bf_layer_bc_procedure'')'
               print '(''make_test_bf_layer_bc_procedure'')'
               print '(''test case not implemented: '', I2)', test_id
               stop ''

          end select

        end subroutine make_test_get_bc_section


        subroutine make_test_analyse_grdpt_with_bc_section(
     $       test_id,
     $       grdpts_id,
     $       test_i,
     $       test_j,
     $       test_bc_section,
     $       test_compatible,
     $       test_remove_ele)
        
          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          integer, dimension(5)               , intent(out) :: test_bc_section
          logical                             , intent(out) :: test_compatible
          logical                             , intent(out) :: test_remove_ele

          
          allocate(grdpts_id(3,3))

          select case(test_id)

            !test W_edge: compatible

            !  -------    
            ! | 2 1 0 |
            ! | 2 1*0 |
            ! | 2 1 0 |
            !  -------
            case(1)
               grdpts_id = reshape(
     $              (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 2

               test_bc_section(1) = W_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 2
               test_bc_section(5) = 0
               test_compatible    = .true.
               test_remove_ele    = .false.

            !test E_edge: compatible

            !  -------    
            ! | 1 2   |
            ! | 1*2   |
            ! | 1 2   |
            !  -------
            case(2)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_pt,no_pt
     $              /),
     $              (/3,3/))

               test_i             = 1
               test_j             = 2

               test_bc_section(1) = E_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 1
               test_bc_section(5) = 0
               test_compatible    = .true.
               test_remove_ele    = .false.

            !test N_edge: compatible

            !  -------    
            ! |       |
            ! | 2 2 2 |
            ! | 1 1*1 |
            !  -------
            case(3)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 1

               test_bc_section(1) = N_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 1
               test_bc_section(5) = 0
               test_compatible    = .true.
               test_remove_ele    = .false.

            !test S_edge: compatible

            !  -------    
            ! | 1 1*1 |
            ! | 2 2 2 |
            ! |       |
            !  -------
            case(4)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 3

               test_bc_section(1) = S_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 3
               test_bc_section(5) = 0
               test_compatible    = .true.
               test_remove_ele    = .false.

            !test W_edge: not compatible

            !  -------
            ! | 1 1 1 |
            ! | 2 2 1*|
            ! |   2 1 |
            !  -------
            case(5)
               grdpts_id = reshape(
     $              (/
     $              no_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 3
               test_j             = 2

               test_bc_section(1) = W_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 3
               test_bc_section(5) = 0
               test_compatible    = .false.
               test_remove_ele    = .false.

            !test E_edge: not compatible

            !  -------    
            ! | 1 1 1 |
            ! | 1*2 2 |
            ! | 1 2   |
            !  -------
            case(6)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 1
               test_j             = 2

               test_bc_section(1) = E_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 1
               test_bc_section(5) = 0
               test_compatible    = .false.
               test_remove_ele    = .false.

            !test N_edge: not compatible

            !  -------    
            ! |   2 1 |
            ! | 2 2 1 |
            ! | 1 1*1 |
            !  -------
            case(7)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              no_pt,bc_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 1

               test_bc_section(1) = N_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 1
               test_bc_section(5) = 0
               test_compatible    = .false.
               test_remove_ele    = .false.

            !test S_edge: not compatible

            !  -------    
            ! | 1 1*1 |
            ! | 2 2 1 |
            ! |   2 1 |
            !  -------
            case(8)
               grdpts_id = reshape(
     $              (/
     $              no_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 3

               test_bc_section(1) = S_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 3
               test_bc_section(5) = 0
               test_compatible    = .false.
               test_remove_ele    = .false.

            !test SW_edge: compatible

            !  -------    
            ! | 1 1 1 |
            ! | 2 2 1*|
            ! |   2 1 |
            !  -------
            case(9)
               grdpts_id = reshape(
     $              (/
     $              no_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 3
               test_j             = 2

               test_bc_section(1) = SW_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 0
               test_bc_section(5) = 1
               test_compatible    = .true.
               test_remove_ele    = .false.

            !test SE_edge: compatible

            !  -------    
            ! | 1 1 1 |
            ! | 1*2 2 |
            ! | 1 2   |
            !  -------
            case(10)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 1
               test_j             = 2

               test_bc_section(1) = SE_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 2
               test_bc_section(4) = 0
               test_bc_section(5) = 2
               test_compatible    = .true.
               test_remove_ele    = .true.

            !test NW_edge: compatible

            !  -------    
            ! |   2 1 |
            ! | 2 2 1 |
            ! | 1 1*1 |
            !  -------
            case(11)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              no_pt,bc_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 1

               test_bc_section(1) = NW_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 1
               test_bc_section(4) = 0
               test_bc_section(5) = 2
               test_compatible    = .true.
               test_remove_ele    = .true.

            !test NE_edge: compatible

            !  -------    
            ! | 0 1 2 |
            ! | 0 1*1 |
            ! | 0 0 0 |
            !  -------
            case(12)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,interior_pt,interior_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 2

               test_bc_section(1) = NE_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 0
               test_bc_section(5) = 2
               test_compatible    = .true.
               test_remove_ele    = .true.

            end select

        end subroutine make_test_analyse_grdpt_with_bc_section


        subroutine test_get_bc_section(detailled)

          implicit none

          logical                 , intent(in) :: detailled

          type(bf_layer_bc_sections)           :: bc_sections

          integer, dimension(5)                :: new_bc_section
          integer                              :: k

          integer, dimension(:,:), allocatable :: grdpts_id
          integer                              :: test_i
          integer                              :: test_j
          integer, dimension(5)                :: test_bc_section
          logical                              :: test_validated

          logical                              :: test_global

          logical                              :: ierror


          test_global = .true.
          
          print '(''test get_bc_section()'')'

          do k=1,26

             call make_test_get_bc_section(
     $            k,
     $            grdpts_id,
     $            test_i,
     $            test_j,
     $            test_bc_section)

             new_bc_section = bc_sections%get_bc_section(
     $            test_i,
     $            test_j,
     $            grdpts_id,
     $            ierror)


             test_validated = test_bc_section(1).eq.new_bc_section(1)

             select case(test_bc_section(1))

               !if this is an edge procedure, only the 2:4
               !elements should be tested
               case(N_edge_type,S_edge_type,E_edge_type,W_edge_type)
                  
                  test_validated = test_validated.and.(
     $                 test_bc_section(2).eq.new_bc_section(2))
                  test_validated = test_validated.and.(
     $                 test_bc_section(3).eq.new_bc_section(3))
                  test_validated = test_validated.and.(
     $                 test_bc_section(4).eq.new_bc_section(4))

                  if(detailled) then
                     print '(''test '',I2,'':'',L1)', k, test_validated
                  end if

                  if(.not.test_validated) then
                     print '(''procedure_type: '',I2, 2X,I2)', test_bc_section(1), new_bc_section(1)
                     print '(''edge_min      : '',I2, 2X,I2)', test_bc_section(2), new_bc_section(2)
                     print '(''edge_max      : '',I2, 2X,I2)', test_bc_section(3), new_bc_section(3)
                     print '(''coord         : '',I2, 2X,I2)', test_bc_section(4), new_bc_section(4)
                  end if

               !if this is a corner procedure, only the 2:3
               !elements should be tested
               case(NE_corner_type, NW_corner_type, SE_corner_type, SW_corner_type)

                  test_validated = test_validated.and.(
     $                 test_bc_section(2).eq.new_bc_section(2))
                  test_validated = test_validated.and.(
     $                 test_bc_section(3).eq.new_bc_section(3))

                  if(detailled) then
                     print '(''test '',I2,'':'',L1)', k, test_validated
                  end if

                  if(.not.test_validated) then
                     print '(''procedure_type: '',I2, 2X,I2)', test_bc_section(1), new_bc_section(1)
                     print '(''i_min         : '',I2, 2X,I2)', test_bc_section(2), new_bc_section(2)
                     print '(''j_min         : '',I2, 2X,I2)', test_bc_section(3), new_bc_section(3)
                  end if

               !if this is a special edge procedure, only the 2:3 and 5
               !elements should be tested
               case(NE_edge_type,NW_edge_type,SE_edge_type,SW_edge_type)

                  test_validated = test_validated.and.(
     $                 test_bc_section(2).eq.new_bc_section(2))
                  test_validated = test_validated.and.(
     $                 test_bc_section(3).eq.new_bc_section(3))
                  test_validated = test_validated.and.(
     $                 test_bc_section(5).eq.new_bc_section(5))

                  if(detailled) then
                     print '(''test '',I2,'':'',L1)', k, test_validated
                  end if

                  if(.not.test_validated) then
                     print '(''procedure_type: '',I2, 2X,I2)', test_bc_section(1), new_bc_section(1)
                     print '(''i_min         : '',I2, 2X,I2)', test_bc_section(2), new_bc_section(2)
                     print '(''j_min         : '',I2, 2X,I2)', test_bc_section(3), new_bc_section(3)
                     print '(''match_nb      : '',I2, 2X,I2)', test_bc_section(5), new_bc_section(5)
                  end if

               case default
                  
                  print '(''test_bf_layer_bc_sections'')'
                  print '(''test_get_bc-sections()'')'
                  print '(''case '', I2, ''not recognized'')', k
                  stop ''

             end select

             deallocate(grdpts_id)

             test_global = test_global.and.test_validated
             
          end do

          print '(''test_validated: '',L1)', test_global
          print '()'

        end subroutine test_get_bc_section


        subroutine test_analyse_grdpt_with_bc_section(detailled)

          implicit none

          logical                 , intent(in) :: detailled

          type(bf_layer_bc_sections)           :: bc_sections

          integer                              :: k
          integer, dimension(:,:), allocatable :: grdpts_id
          logical                              :: compatible
          logical                              :: remove_ele

          integer                              :: test_i
          integer                              :: test_j
          integer, dimension(5)                :: test_bc_section
          logical                              :: test_validated
          logical                              :: test_compatible
          logical                              :: test_remove_ele

          logical :: test_global

          test_global = .true.          

          
          print '(''test analyse_grdpt_with_bc_section()'')'
          do k=1,12

             call make_test_analyse_grdpt_with_bc_section(
     $            k,
     $            grdpts_id,
     $            test_i,
     $            test_j,
     $            test_bc_section,
     $            test_compatible,
     $            test_remove_ele)

             compatible = bc_sections%analyse_grdpt_with_bc_section(
     $            test_i,test_j,grdpts_id,test_bc_section,remove_ele)

             test_validated = compatible.eqv.test_compatible
             test_validated = test_validated.and.(remove_ele.eqv.test_remove_ele)

             if(detailled) then
                print '(''test '',I2,'':'',L1)', k, test_validated
             end if

             if(.not.test_validated) then
                print '(''  compatible    : '',L1, 2X,L1)', test_compatible, compatible
                print '(''  remove_ele    : '',L1, 2X,L1)', test_remove_ele, remove_ele
             end if

             deallocate(grdpts_id)

             test_global = test_global.and.test_validated

          end do
          print '(''test_validated: '',L1)', test_global
          print '()'

        end subroutine test_analyse_grdpt_with_bc_section


        subroutine test_analyse_grdpt(detailled)

          implicit none

          logical, intent(in) :: detailled

          type(bf_layer_bc_sections)           :: bc_sections
          integer                              :: i,j,k
          integer                              :: test_i_min
          integer                              :: test_i_max
          integer                              :: test_j_min
          integer                              :: test_j_max
          integer, dimension(:,:), allocatable :: grdpts_id
          integer, dimension(:,:), allocatable :: test_bc_sections_sorted

          integer, dimension(:,:), allocatable :: bc_sections_sorted

          logical :: same
          logical :: test_global
          
          logical :: ierror

          print '(''test_analyse_grdpt()'')'          

          test_global = .true.

          do k=1,3

             !initialize object for boundary layers
             call bc_sections%ini()

             !create test
             call make_test_analyse_grdpt(
     $            k,
     $            test_i_min, test_i_max,
     $            test_j_min, test_j_max,
     $            grdpts_id,
     $            test_bc_sections_sorted)

             !make test
             do j=test_j_min, test_j_max
                do i=test_i_min,test_i_max

                   if(grdpts_id(i,j).eq.bc_interior_pt) then
                      call bc_sections%analyse_grdpt(i,j,grdpts_id,ierror)
                   end if

                end do
             end do

             call bc_sections%sort_bc_sections(bc_sections_sorted)

             !compare results
             same = compare_bc_sections_sorted(
     $            bc_sections_sorted,
     $            test_bc_sections_sorted)

             if(detailled) then
                print '(''test '',I2,'': '',L1)', k, same
                if(.not.same) then
                   call print_bc_sections(bc_sections_sorted)
                end if
             end if

             test_global = test_global.and.same

             !deallocate intermediate variables
             if(allocated(test_bc_sections_sorted)) then
                deallocate(test_bc_sections_sorted)
             end if

             !deallocate the object with the boundary
             !layers
             call bc_sections%deallocate_tables()

          end do

          print '(''test_validated: '', L1)', test_global
          print '()'

        end subroutine test_analyse_grdpt

        
        subroutine print_bc_sections(bc_sections)

          implicit none

          integer, dimension(:,:), intent(in) :: bc_sections

          integer :: i

          do i=1, size(bc_sections,2)

             select case(bc_sections(1,i))

               case(N_edge_type)
                  print '(I2,'' N_edge    [i,j]: '',2I3,'' i_max='',I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i),
     $                 bc_sections(4,i)

               case(S_edge_type)
                  print '(I2,'' S_edge    [i,j]: '',2I3,'' i_max='',I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i),
     $                 bc_sections(4,i)

               case(E_edge_type)
                  print '(I2,'' E_edge    [i,j]: '',2I3,'' j_max='',I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i),
     $                 bc_sections(4,i)

               case(W_edge_type)
                  print '(I2,'' W_edge    [i,j]: '',2I3,'' j_max='',I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i),
     $                 bc_sections(4,i)

               case(NE_corner_type)
                  print '(I2,'' NE_corner [i,j]: '',2I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i)

               case(NW_corner_type)
                  print '(I2,'' NW_corner [i,j]: '',2I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i)

               case(SW_corner_type)
                  print '(I2,'' SW_corner [i,j]: '',2I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i)

               case(SE_corner_type)
                  print '(I2,'' SE_corner [i,j]: '',2I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i)

               case(NE_edge_type)
                  print '(I2,'' NE_edge   [i,j]: '',2I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i)

               case(NW_edge_type)
                  print '(I2,'' NW_edge   [i,j]: '',2I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i)

               case(SW_edge_type)
                  print '(I2,'' SW_edge   [i,j]: '',2I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i)

               case(SE_edge_type)
                  print '(I2,'' SE_edge   [i,j]: '',2I3)',
     $                 i,
     $                 bc_sections(2,i), bc_sections(3,i)

             end select

          end do

        end subroutine print_bc_sections


        subroutine make_test_analyse_grdpt(
     $     test_id,
     $     test_i_min, test_i_max,
     $     test_j_min, test_j_max,
     $     grdpts_id,
     $     test_bc_sections_sorted)

          implicit none

          integer                             , intent(in)  :: test_id
          integer                             , intent(out) :: test_i_min
          integer                             , intent(out) :: test_i_max
          integer                             , intent(out) :: test_j_min
          integer                             , intent(out) :: test_j_max
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer, dimension(:,:), allocatable, intent(out) :: test_bc_sections_sorted


          integer :: i,j

          
          select case(test_id)

            !  ---- | ----- | ----
            ! |     |       |     |
            ! | 2 2 | 2 2 2 | 2 2 |
            ! | 1 1 | 1 1 1 | 1 1 |
            ! | 0 0 | 0 0 0 | 0 0 |
            ! ---------------------
            ! | 0 0 | 0 0 0 | 0 0 |
            ! | 0 0 | 0 0 0 | 0 0 |
            !  ---- | ----- | ----
            !
            case(1)
               test_i_min = 3
               test_i_max = 5
               test_j_min = 3
               test_j_max = 5

               allocate(grdpts_id(7,6))
               
               do j=1,3
                  do i=1,7
                     grdpts_id(i,j) = interior_pt
                  end do
               end do

               j=4
               do i=1,7
                  grdpts_id(i,j) = bc_interior_pt
               end do

               j=5
               do i=1,7
                  grdpts_id(i,j) = bc_pt
               end do

               j=6
               do i=1,7
                  grdpts_id(i,j) = no_pt
               end do

               allocate(test_bc_sections_sorted(4,1))
               test_bc_sections_sorted = reshape(
     $              (/
     $              N_edge_type,3,4,5
     $              /),
     $              (/4,1/))


            !  ---- | ------- | ----
            ! |     |         |     |
            ! |     | 2 2 2 2 | 2 2 |
            ! |     | 2 1 1 1 | 1 1 |
            ! |     | 2 1 0 0 | 0 0 |
            ! -----------------------
            ! |     | 2 1 0 0 | 0 0 |
            ! |     | 2 1 0 0 | 0 0 |
            !  ---- | ------- | ----
            !
            case(2)
               test_i_min = 3
               test_i_max = 6
               test_j_min = 3
               test_j_max = 5

               allocate(grdpts_id(8,6))

               do j=1,3
                  do i=1,2
                     grdpts_id(i,j) = no_pt
                  end do
                  grdpts_id(3,j) = bc_pt
                  grdpts_id(4,j) = bc_interior_pt
                  do i=5,8
                     grdpts_id(i,j) = interior_pt
                  end do
               end do

               j=4
               do i=1,2
                  grdpts_id(i,j) = no_pt
               end do
               grdpts_id(3,j) = bc_pt
               do i=4,8
                  grdpts_id(i,j) = bc_interior_pt
               end do
               
               j=5
               do i=1,2
                  grdpts_id(i,j) = no_pt
               end do
               do i=3,8
                  grdpts_id(i,j) = bc_pt
               end do

               j=6
               do i=1,8
                  grdpts_id(i,j) = no_pt
               end do

               allocate(test_bc_sections_sorted(4,3))
               test_bc_sections_sorted = reshape(
     $              (/
     $              W_edge_type,3,3,3,
     $              NW_corner_type,3,4,0,
     $              N_edge_type,5,4,6
     $              /),
     $              (/4,3/))

            !    |     |         2 2 2 2 2 2 2 2 2                |    |
            !    |     |         2 1 1 1 1 1 1 1 2                |    |
            !    |     |         2 1 0 0 0 0 0 1 2                |    |
            !    |     |         2 1 0       0 1 2                |    |
            !18- |     |         2 1 0       0 1 2 2 2 2 2        |    |
            !17- |     | 2 2 2 2 2 1 0       0 1 1 1 1 1 2        |    |
            !    |     | 2 1 1 1 1 1 0       0 0 0 0 0 1 2        |    |
            !15- |     | 2 1 0 0 0 0 0               0 1 2        |    |
            !    |     | 2 1 0                       0 1 2        |    |
            !13- |     | 2 1 0                       0 1 2 2 2 2  |    |
            !    |     | 2 1 0                       0 1 1 1 1 2  |    |
            !11- |     | 2 1 0 0 0                   0 0 0 0 1 2  |    |
            !10- |     | 2 1 1 1 0                         0 1 2  |    |
            !    |     | 2 2 2 1 0             0 0 0 0 0 0 0 1 2  |    |
            ! 8- |     |     2 1 0 0           0 1 1 1 1 1 1 1 2  |    |
            !    |     |     2 1 1 0           0 1 2 2 2 2 2 2 2  |    |
            !    |     |     2 2 1 0           0 1 2              |    |
            ! 5- |     |       2 1 0           0 1 2              |    |
            !    |     | 2 2 2 2 1 0           0 1 2              |    |
            !    |     | 2 1 1 1 1 0           0 1 2              |    |
            !    -------------------------------------------------|----|
            !    |     | 2 1 0 0               0 1 2              |    |
            !    |     | 2 1 0 0               0 1 2              |    |
            !     ---- | ----------------------------------------- ---  
            !            |     | | |           |     |         |        
            !            3     6 7 8          14    17        22
            ! --------------------------------------------------
            case(3)

               test_i_min = 3
               test_i_max = 23
               test_j_min = 3
               test_j_max = 22

               allocate(grdpts_id(25,22))

               grdpts_id = reshape(
     $              (/
     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     $              0, 0, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     $              0, 0, 3, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     $              0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     $              0, 0, 0, 0, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $              0, 0, 0, 0, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0,
     $              0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 0, 0, 0,
     $              0, 0, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0,
     $              0, 0, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0,
     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0,
     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 0, 0, 0,
     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 0, 0, 0,
     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0,
     $              0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0,
     $              0, 0, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0,
     $              0, 0, 3, 3, 3, 3, 3, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0, 0, 0,
     $              0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0,
     $              0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $              0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $              0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $              0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     $              /),
     $              (/25,22/))

               allocate(test_bc_sections_sorted(4,30))
               test_bc_sections_sorted = reshape(
     $              (/
     $              NW_corner_type,3,3,0,
     $              N_edge_type,5,3,5,
     $              NW_edge_type,6,3,0,
     $              E_edge_type,15,3,6,
     $              W_edge_type,6,5,5,
     $              SW_corner_type,5,6,0,
     $              SW_edge_type,6,6,0,
     $              SE_edge_type,15,7,0,
     $              S_edge_type,17,7,20,
     $              SE_corner_type,21,7,0,
     $              W_edge_type,5,8,8,
     $              SW_corner_type,3,9,0,
     $              SW_edge_type,5,9,0,
     $              E_edge_type,21,9,11,
     $              W_edge_type,3,11,15,
     $              NE_edge_type,18,12,0,
     $              N_edge_type,20,12,20,
     $              NE_corner_type,21,12,0,
     $              E_edge_type,18,14,16,
     $              NW_corner_type,3,16,0,
     $              N_edge_type,5,16,6,
     $              NW_edge_type,7,16,0,
     $              NE_edge_type,14,17,0,
     $              N_edge_type,16,17,17,
     $              NE_corner_type,18,17,0,
     $              W_edge_type,7,18,20,
     $              E_edge_type,14,19,20,
     $              NW_corner_type,7,21,0,
     $              N_edge_type,9,21,13,
     $              NE_corner_type,14,21,0
     $              /),
     $              (/4,30/))
               
            case default
               print '(''test_bf_layer_bc_sections'')'
               print '(''make_test_analyse_grdpt'')'
               print '(''case : '', I2)', test_id
               stop
          end select

c$$$          allocate(test_bc_sections_buffer(1,1))
c$$$          deallocate(test_bc_sections_buffer)

        end subroutine make_test_analyse_grdpt


        subroutine test_add_overlap_between_corners_and_anti_corners(detailled)

          implicit none

          logical, intent(in) :: detailled


          type(bf_layer_bc_sections) :: bf_layer_bc_sections_used
          integer, dimension(4,18)   :: bc_sections_sorted
          integer, dimension(4,18)   :: test_bc_sections_sorted

          logical :: test_validated


          !initialization of bc_sections_sorted
          bc_sections_sorted(:, 1) = [NE_edge_type  ,1, 1,0] !-> no_overlap
          bc_sections_sorted(:, 2) = [E_edge_type   ,1, 1,3]
          bc_sections_sorted(:, 3) = [NW_edge_type  ,1, 2,0] !-> N_overlap
          bc_sections_sorted(:, 4) = [N_edge_type   ,2, 2,3]
          bc_sections_sorted(:, 5) = [NW_corner_type,1, 3,0] !-> for N_overlap
          bc_sections_sorted(:, 6) = [S_edge_type   ,2, 3,3]
          bc_sections_sorted(:, 7) = [SW_edge_type  ,3, 3,0] !-> NE_overlap
          bc_sections_sorted(:, 8) = [SE_corner_type,4, 3,0] !-> for E_overlap
          bc_sections_sorted(:, 9) = [SW_corner_type,3, 4,0] !-> for N_overlap
          bc_sections_sorted(:,10) = [W_edge_type   ,2, 5,3]
          bc_sections_sorted(:,11) = [SW_corner_type,6, 5,0] !-> for W_overlap
          bc_sections_sorted(:,12) = [SE_edge_type  ,7, 5,0] !-> NW_overlap
          bc_sections_sorted(:,13) = [SE_corner_type,7, 6,0] !-> for N_overlap
          bc_sections_sorted(:,14) = [SW_corner_type,6, 8,0] !-> for W_overlap
          bc_sections_sorted(:,15) = [SE_edge_type  ,7, 8,0] !-> W_overlap
          bc_sections_sorted(:,16) = [SE_edge_type  ,8,10,0] !-> E_overlap
          bc_sections_sorted(:,17) = [W_edge_type   ,9,10,9]
          bc_sections_sorted(:,18) = [SW_corner_type,9,10,0] !-> for E_overlap
          

          !initialization of test_bc_sections_sorted
          test_bc_sections_sorted       = bc_sections_sorted

          test_bc_sections_sorted(4,3)  = N_overlap
          test_bc_sections_sorted(4,7)  = NE_overlap
          test_bc_sections_sorted(4,12) = NW_overlap
          test_bc_sections_sorted(4,15) = W_overlap
          test_bc_sections_sorted(4,16) = E_overlap          
          

          !computation of the overlap
          call bf_layer_bc_sections_used%add_overlap_between_corners_and_anti_corners(
     $         bc_sections_sorted)


          !check the result of the computation of the overlap
          test_validated = is_int_matrix_validated(
     $         bc_sections_sorted,
     $         test_bc_sections_sorted,
     $         detailled)

          print '(''add_overlap_between_corners_and_anti_corners'')'
          print '(''test_validated: '',L1)', test_validated
          print '()'

        end subroutine test_add_overlap_between_corners_and_anti_corners

      end program test_bf_layer_bc_sections

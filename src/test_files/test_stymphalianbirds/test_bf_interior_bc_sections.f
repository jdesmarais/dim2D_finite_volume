      !> @file
      !> test the subroutines of bf_interior_bc_sections_module.f
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the subroutines of bf_interior_bc_sections_module.f
      !
      !> @date
      !> 29_10_2014 - initial version         - J.L. Desmarais
      !-----------------------------------------------------------------
       program test_bf_interior_bc_sections

        use bf_layer_bc_procedure_module, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       SW_corner_type,
     $       SE_corner_type,
     $       NW_corner_type,
     $       NE_corner_type

        use bf_interior_bc_sections_module, only :
     $       ini_interior_bc_sections,
     $       determine_interior_bc_sections,
     $       close_last_bc_section,
     $       set_full_interior_bc_section,
     $       minimize_interior_bc_section,
     $       process_bc_sections_into_bc_procedure

        use parameters_input, only :
     $       nx,
     $       ny,
     $       bc_size

        use parameters_kind, only : 
     $       ikind

        implicit none


        integer(ikind), dimension(:,:), allocatable :: bf_alignments
        integer(ikind), dimension(:,:), allocatable :: test_bc_sections
        integer(ikind), dimension(:,:), allocatable :: test_bc_procedures
        integer                                     :: test_interior_inf
        integer                                     :: test_interior_sup
        integer                                     :: test_case_id
        logical                                     :: detailled
        logical                                     :: test_validated

        
        if(  (nx.ne.10).or.
     $       (ny.ne.10)) then

           stop 'the test requires (nx,ny)=(10,10)'

        end if


        detailled = .false.

        do test_case_id=1,12

           call ini_bf_alignments(
     $          test_case_id,
     $          bf_alignments,
     $          test_bc_sections,
     $          test_bc_procedures,
     $          test_interior_inf,
     $          test_interior_sup)
           
           test_validated = make_test_bc_sections(
     $          bf_alignments,
     $          test_bc_sections,
     $          test_bc_procedures,
     $          test_interior_inf,
     $          test_interior_sup,
     $          detailled)

           print '(''test '',I3,'': '',L1)',
     $          test_case_id,
     $          test_validated

        end do

        print '()'

        contains

        subroutine ini_bf_alignments(
     $       test_case_id,
     $       bf_alignments,
     $       test_bc_sections,
     $       test_bc_procedures,
     $       test_interior_inf,
     $       test_interior_sup)

          implicit none

          integer                                    , intent(in)  :: test_case_id
          integer(ikind), dimension(:,:), allocatable, intent(out) :: bf_alignments
          integer(ikind), dimension(:,:), allocatable, intent(out) :: test_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(out) :: test_bc_procedures
          integer(ikind)                             , intent(out) :: test_interior_inf
          integer(ikind)                             , intent(out) :: test_interior_sup


          if(allocated(bf_alignments)) then
             deallocate(bf_alignments)
          end if

          if(allocated(test_bc_sections)) then
             deallocate(test_bc_sections)
          end if

          if(allocated(test_bc_procedures)) then
             deallocate(test_bc_procedures)
          end if
          

          select case(test_case_id)

             !no buffer layer 
             case(1)
                allocate(test_bc_sections(2,1))
                test_bc_sections(:,1) = [1,nx]
                
                allocate(test_bc_procedures(4,3))
                test_bc_procedures(:,1) = [NW_corner_type,1,ny-bc_size+1,bc_size]
                test_bc_procedures(:,2) = [N_edge_type,bc_size+1,ny-bc_size+1,nx-bc_size]
                test_bc_procedures(:,3) = [NE_corner_type,nx-bc_size+1,ny-bc_size+1,nx]

                test_interior_inf = 1
                test_interior_sup = nx

             !buffer layer in the interior region
             case(2)
                allocate(bf_alignments(2,1))
                bf_alignments(:,1) = [5,nx-4]
                
                allocate(test_bc_sections(2,2))
                test_bc_sections(:,1) = [1,2]
                test_bc_sections(:,2) = [nx-1,nx]
                
                allocate(test_bc_procedures(4,2))
                test_bc_procedures(:,1) = [NW_corner_type,1,ny-bc_size+1,2]
                test_bc_procedures(:,2) = [NE_corner_type,nx-1,ny-bc_size+1,nx]

                test_interior_inf = 1
                test_interior_sup = nx
                
             !buffer layer on the left side
             case(3)
                allocate(bf_alignments(2,1))
                bf_alignments(:,1) = [-3,nx-4]
                
                allocate(test_bc_sections(2,1))
                test_bc_sections(:,1) = [nx-1,nx]

                allocate(test_bc_procedures(4,1))
                test_bc_procedures(:,1) = [NE_corner_type,nx-1,ny-bc_size+1,nx]

                test_interior_inf = 1
                test_interior_sup = nx
                
             !buffer layer on the right side
             case(4)
                allocate(bf_alignments(2,1))
                bf_alignments(:,1) = [nx-5,nx+1]
                
                allocate(test_bc_sections(2,1))
                test_bc_sections(:,1) = [1,nx-8]

                allocate(test_bc_procedures(4,1))
                test_bc_procedures(:,1) = [NW_corner_type,1,ny-bc_size+1,bc_size]

                test_interior_inf = 1
                test_interior_sup = nx

             !buffer layer just on the left side
             case(5)
                allocate(bf_alignments(2,1))
                bf_alignments(:,1) = [3,nx-5]
                
                allocate(test_bc_sections(2,1))
                test_bc_sections(:,1) = [nx-2,nx]

                allocate(test_bc_procedures(4,2))
                test_bc_procedures(:,1) = [N_edge_type,nx-2,ny-bc_size+1,nx-2]
                test_bc_procedures(:,2) = [NE_corner_type,nx-1,ny-bc_size+1,nx]

                test_interior_inf = 1
                test_interior_sup = nx

             !buffer layer just on the right side
             case(6)
                allocate(bf_alignments(2,1))
                bf_alignments(:,1) = [nx-4,nx-2]
                
                allocate(test_bc_sections(2,1))
                test_bc_sections(:,1) = [1,nx-7]

                allocate(test_bc_procedures(4,2))
                test_bc_procedures(:,1) = [NW_corner_type,1,ny-bc_size+1,bc_size]
                test_bc_procedures(:,2) = [N_edge_type,bc_size+1,ny-bc_size+1,nx-7]

                test_interior_inf = 1
                test_interior_sup = nx

             !buffer layer covering the interior
             case(7)
                allocate(bf_alignments(2,1))
                bf_alignments(:,1) = [1,nx]
                
                test_interior_inf = 1
                test_interior_sup = nx

             !buffer layer just covering the interior
             case(8)
                allocate(bf_alignments(2,1))
                bf_alignments(:,1) = [3,nx-2]
                
                test_interior_inf = 1
                test_interior_sup = nx

             !2 buffer layers: one on the outside left and one inside
             case(9)
                allocate(bf_alignments(2,2))
                bf_alignments(:,1) = [-10,-2]
                bf_alignments(:,2) = [5,nx-4]
                
                allocate(test_bc_sections(2,2))
                test_bc_sections(:,1) = [1,2]
                test_bc_sections(:,2) = [nx-1,nx]

                allocate(test_bc_procedures(4,2))
                test_bc_procedures(:,1) = [NW_corner_type,1,ny-bc_size+1,bc_size]
                test_bc_procedures(:,2) = [NE_corner_type,nx-1,ny-bc_size+1,nx]

                test_interior_inf = 1
                test_interior_sup = nx

             !2 buffer layers: one on the outside left, one on the
             !outside right
             case(10)
                allocate(bf_alignments(2,2))
                bf_alignments(:,1) = [-10,-2]
                bf_alignments(:,2) = [nx+2,nx+6]
                
                allocate(test_bc_sections(2,1))
                test_bc_sections(:,1) = [1,nx]

                allocate(test_bc_procedures(4,3))
                test_bc_procedures(:,1) = [NW_corner_type,1,ny-bc_size+1,bc_size]
                test_bc_procedures(:,2) = [N_edge_type,bc_size+1,ny-bc_size+1,nx-bc_size]
                test_bc_procedures(:,3) = [NE_corner_type,nx-1,ny-bc_size+1,nx]

                test_interior_inf = 1
                test_interior_sup = nx

             !2 buffer layers in the interior: one on the left
             !one inside
             case(11)
                allocate(bf_alignments(2,2))
                bf_alignments(:,1) = [-10,0]
                bf_alignments(:,2) = [6,6]
                
                allocate(test_bc_sections(2,2))
                test_bc_sections(:,1) = [3,3]
                test_bc_sections(:,2) = [9,10]

                allocate(test_bc_procedures(4,2))
                test_bc_procedures(:,1) = [N_edge_type,3,ny-bc_size+1,3]
                test_bc_procedures(:,2) = [NE_corner_type,nx-1,ny-bc_size+1,nx]

                test_interior_inf = 1
                test_interior_sup = nx

             !2 buffer layers in the interior: one inside,
             !one in the outside right
             case(12)
                allocate(bf_alignments(2,2))
                bf_alignments(:,1) = [5,5]
                bf_alignments(:,2) = [11,16]
                
                allocate(test_bc_sections(2,2))
                test_bc_sections(:,1) = [1,2]
                test_bc_sections(:,2) = [8,8]

                allocate(test_bc_procedures(4,2))
                test_bc_procedures(:,1) = [NW_corner_type,1,ny-bc_size+1,bc_size]
                test_bc_procedures(:,2) = [N_edge_type,8,ny-bc_size+1,8]

                test_interior_inf = 1
                test_interior_sup = nx

             case default

                stop 'ini_bf_alignment: case not implemented'
           
          end select

        end subroutine ini_bf_alignments


        function make_test_bc_sections(
     $     bf_alignments,
     $     test_bc_sections,
     $     test_bc_procedures,
     $     interior_inf,
     $     interior_sup,
     $     detailled)
     $     result(test_validated)

          implicit none
          
          integer(ikind), dimension(:,:), allocatable, intent(in) :: bf_alignments
          integer(ikind), dimension(:,:), allocatable, intent(in) :: test_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(in) :: test_bc_procedures
          integer(ikind)                             , intent(in) :: interior_inf
          integer(ikind)                             , intent(in) :: interior_sup
          logical                                    , intent(in) :: detailled
          logical                                                 :: test_validated

          
          integer                                     :: k
          integer(ikind), dimension(2)                :: bf_alignment

          integer                                     :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable :: bc_sections
          integer(ikind), dimension(:,:), allocatable :: bc_sections_tmp
          integer(ikind), dimension(:,:), allocatable :: bc_procedures
          logical                                     :: min_initialized
          logical                                     :: max_initialized
          logical                                     :: no_bf_common_with_interior

          logical :: test_bc_sec
          logical :: test_bc_proc


          !initialize the interior_bc_sections
          call ini_interior_bc_sections(
     $         nb_bc_sections,
     $         min_initialized,
     $         max_initialized,
     $         no_bf_common_with_interior)

          if(allocated(bf_alignments)) then

             !determine the interior_bc_sections
             do k=1, size(bf_alignments,2)
             
                bf_alignment = bf_alignments(:,k)
                
                call determine_interior_bc_sections(
     $               bf_alignment,
     $               interior_inf,
     $               interior_sup,
     $               nb_bc_sections,
     $               bc_sections,
     $               min_initialized,
     $               max_initialized,
     $               no_bf_common_with_interior)
                
             end do


            !finalize the interior_bc_sections
             call close_last_bc_section(
     $            nb_bc_sections,
     $            bc_sections,
     $            interior_sup,
     $            min_initialized,
     $            max_initialized)

          end if

          if(no_bf_common_with_interior) then
             call set_full_interior_bc_section(
     $            nb_bc_sections,
     $            bc_sections,
     $            min_initialized,
     $            max_initialized,
     $            interior_inf,
     $            interior_sup)
          end if

          call minimize_interior_bc_section(
     $         nb_bc_sections,
     $         bc_sections)           

          !compare the bc_sections computed with the 
          !test_bc_sections of the test
          test_bc_sec = compare_bc_sections(
     $         nb_bc_sections,
     $         bc_sections,
     $         test_bc_sections,
     $         detailled)
          
          if(detailled) then
             print '(''test_bc_sections: '',L1)', test_bc_sec
          end if


          !compute the bc_procedures
          call process_bc_sections_into_bc_procedure(
     $         bc_sections,
     $         bc_sections_tmp,
     $         bc_sections_tmp,
     $         bc_sections_tmp,
     $         bc_procedures)

          test_bc_proc = compare_bc_procedures(
     $         test_bc_procedures,
     $         bc_procedures)

          if(detailled) then
             print '(''test_bc_procedures: '',L1)', test_bc_proc
          end if


          !print bc_sections
          if(detailled) then
             call print_bc_sections(
     $            nb_bc_sections,
     $            bc_sections)

             call print_bc_procedures(
     $            bc_procedures)

          end if

          test_validated = test_bc_sec.and.test_bc_proc

        end function make_test_bc_sections


        function compare_bc_sections(
     $     nb_bc_sections,
     $     bc_sections,
     $     test_bc_sections,
     $     detailled)
     $     result(test_validated)

          implicit none

          integer                                    , intent(in) :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(in) :: bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(in) :: test_bc_sections
          logical                                    , intent(in) :: detailled
          logical                                                 :: test_validated
          
          
          integer :: k
          logical :: loc


          if(allocated(test_bc_sections)) then
             
             test_validated = .true.

             do k=1, size(test_bc_sections,2)

                loc =
     $               (test_bc_sections(1,k).eq.bc_sections(1,k)).and.
     $               (test_bc_sections(2,k).eq.bc_sections(2,k))
                
                if(detailled.and.(.not.loc)) then
                   
                   print '(''  k='',I2,'','',1X,2I3,''|'',2I3)',
     $                  k,
     $                  test_bc_sections(:,k),
     $                  bc_sections(:,k)

                end if

                test_validated = test_validated.and.loc

             end do

          else

             test_validated = (nb_bc_sections.eq.0)

          end if
          
        end function compare_bc_sections


        function compare_bc_procedures(
     $     bc_procedures,
     $     test_bc_procedures)
     $     result(test_validated)

          implicit none

          integer(ikind), dimension(:,:), allocatable, intent(in) :: bc_procedures
          integer(ikind), dimension(:,:), allocatable, intent(in) :: test_bc_procedures
          logical                                                 :: test_validated


          integer :: k

          
          if(allocated(test_bc_procedures)) then
             
             test_validated = .true.

             do k=1, size(test_bc_procedures,2)

                test_validated = test_validated.and.
     $               (test_bc_procedures(1,k).eq.bc_procedures(1,k))

                select case(test_bc_procedures(1,k))
                  case(N_edge_type,S_edge_type,
     $                 E_edge_type,W_edge_type)
                  
                  test_validated = test_validated.and.
     $                 (test_bc_procedures(2,k).eq.bc_procedures(2,k))
                  test_validated = test_validated.and.
     $                 (test_bc_procedures(3,k).eq.bc_procedures(3,k))
                  test_validated = test_validated.and.
     $                 (test_bc_procedures(4,k).eq.bc_procedures(4,k))
                  
                  case(NE_corner_type,NW_corner_type,
     $                 SE_corner_type, SW_corner_type)

                  test_validated = test_validated.and.
     $                 (test_bc_procedures(2,k).eq.bc_procedures(2,k))
                  test_validated = test_validated.and.
     $                 (test_bc_procedures(3,k).eq.bc_procedures(3,k))

                  case default
                     print '(''test_bf_interior_bc_sections'')'
                     print '(''compare_bc_procedures'')'
                     stop 'case not recognized'

                end select

             end do

          else
             test_validated = .not.allocated(bc_procedures)
          end if

        end function compare_bc_procedures

      
        subroutine print_bc_sections(nb_bc_sections,bc_sections)

          implicit none

          integer                       , intent(in) :: nb_bc_sections
          integer(ikind), dimension(:,:), intent(in) :: bc_sections

          integer :: k

          print '()'

          if(nb_bc_sections.gt.0) then
             do k=1, nb_bc_sections
                print '(''bc_section '',I2, '' : '',2I3)',
     $               k,
     $               bc_sections(:,k)
             end do

          else
             print '(''no bc_sections'')'
          end if

          
        end subroutine print_bc_sections

        
        subroutine print_bc_procedures(bc_procedures)

          implicit none

          integer(ikind), dimension(:,:), allocatable, intent(in) :: bc_procedures

          integer :: k

          print '()'

          if(allocated(bc_procedures)) then
             
             do k=1, size(bc_procedures,2)

                select case(bc_procedures(1,k))

                  case(N_edge_type)
                     print '(''bc_proc '',I2,'': N_edge   '',2I2)',
     $                    k,
     $                    bc_procedures(2,k),
     $                    bc_procedures(4,k)
                     
                  case(S_edge_type)
                     print '(''bc_proc '',I2,'': S_edge   '',2I2)',
     $                    k,
     $                    bc_procedures(2,k),
     $                    bc_procedures(4,k)

                  case(E_edge_type)
                     print '(''bc_proc '',I2,'': E_edge   '',2I2)',
     $                    k,
     $                    bc_procedures(3,k),
     $                    bc_procedures(4,k)

                  case(W_edge_type)
                     print '(''bc_proc '',I2,'': W_edge   '',2I2)',
     $                    k,
     $                    bc_procedures(3,k),
     $                    bc_procedures(4,k)

                  case(NW_corner_type)
                     print '(''bc_proc '',I2,'': NW_corner'')', k
                     
                  case(NE_corner_type)
                     print '(''bc_proc '',I2,'': NE_corner'')', k

                  case(SW_corner_type)
                     print '(''bc_proc '',I2,'': SW_corner'')', k

                  case(SE_corner_type)
                     print '(''bc_proc '',I2,'': SE_corner'')', k

                  case default
                     print '(''test_bf_interior_bc_sections'')'
                     print '(''print_bc_procedures'')'
                     stop 'case not recognized'

                end select

             end do

          else
             print '(''no bc_procedures'')'
          end if
          
        end subroutine print_bc_procedures

      end program test_bf_interior_bc_sections

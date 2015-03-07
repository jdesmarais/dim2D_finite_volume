      program test_bf_mainlayer_print

        use bf_mainlayer_print_class, only :
     $     bf_mainlayer_print

        use bf_sublayer_class, only :
     $       bf_sublayer
        
        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W

        use parameters_constant, only :
     $       N

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       rkind

        implicit none


        call test_print_netcdf()
        print '(''check_files :'')'
        print '(''   - N_1_10.nc'')'
        print '(''   - N_2_10.nc'')'
        print '(''   - N_3_10.nc'')'
        print '()'


        contains


        subroutine test_print_netcdf()

          implicit none

          type(bf_mainlayer_print)   :: bf_mainlayer_used
          integer, dimension(2,2,3)  :: bf_alignment
          type(bf_sublayer), pointer :: added_sublayer
          real(rkind), dimension(nx) :: interior_x_map
          real(rkind), dimension(ny) :: interior_y_map

          real(rkind), dimension(9    ,6,ne) :: bf_nodes1
          real(rkind), dimension(nx-16,6,ne) :: bf_nodes2
          real(rkind), dimension(9    ,6,ne) :: bf_nodes3
          real(rkind), dimension(nx,ny,ne)   :: interior_nodes

          integer :: i,j,k

          
          !input
          bf_nodes1 = reshape((/
     $         (((500*(k-1) + 50*(align_N-3+j-1) + 5*(align_W-5+i-1),
     $         i=1,size(bf_nodes1,1)),j=1,size(bf_nodes1,2)),k=1,ne)/),
     $         (/size(bf_nodes1,1),size(bf_nodes1,2),ne/))

          bf_nodes2 = reshape((/
     $         (((500*(k-1) + 50*(align_N-3+j-1) + 5*(align_W+6+i-1),
     $         i=1,size(bf_nodes2,1)),j=1,size(bf_nodes2,2)),k=1,ne)/),
     $         (/size(bf_nodes2,1),size(bf_nodes2,2),ne/))

          bf_nodes3 = reshape((/
     $         (((500*(k-1) + 50*(align_N-3+j-1) + 5*(align_E-5+i-1),
     $         i=1,size(bf_nodes3,1)),j=1,size(bf_nodes3,2)),k=1,ne)/),
     $         (/size(bf_nodes3,1),size(bf_nodes3,2),ne/))
          
          interior_nodes = reshape((/
     $         (((500*(k-1) + 50*(j-1) + 5*(i-1),
     $         i=1,size(interior_nodes,1)),j=1,size(interior_nodes,2)),k=1,ne)/),
     $         (/size(interior_nodes,1),size(interior_nodes,2),ne/))


          call bf_mainlayer_used%ini(N)

          bf_alignment(:,:,1) = reshape((/
     $         align_W-2, align_N, align_W+2, align_N+1/),
     $         (/2,2/))
          
          bf_alignment(:,:,2) = reshape((/
     $         align_W+9, align_N, align_E-9, align_N+1/),
     $         (/2,2/))

          bf_alignment(:,:,3) = reshape((/
     $         align_E-2, align_N, align_E+2, align_N+1/),
     $         (/2,2/))

          added_sublayer => bf_mainlayer_used%add_sublayer(
     $             interior_x_map,
     $             interior_y_map,
     $             interior_nodes,
     $             bf_alignment(:,:,1))
          added_sublayer%nodes = bf_nodes1

          added_sublayer => bf_mainlayer_used%add_sublayer(
     $             interior_x_map,
     $             interior_y_map,
     $             interior_nodes,
     $             bf_alignment(:,:,2))
          added_sublayer%nodes = bf_nodes2

          added_sublayer => bf_mainlayer_used%add_sublayer(
     $             interior_x_map,
     $             interior_y_map,
     $             interior_nodes,
     $             bf_alignment(:,:,3))
          added_sublayer%nodes = bf_nodes3


          !output
          call bf_mainlayer_used%print_netcdf(
     $         10,
     $         ['md','qx','qy','Ed'],
     $         ['mass_density','momentum_x','momentum_y','total_energy'],
     $         ['(kg.m-3)/(kg.m-3)','(kg.m-2.s-1)/(kg.m-2.s-1)','(kg.m-2.s-1)/(kg.m-2.s-1)','(J.m-3)/(J.m-3)'],
     $         0.10d0)

        end subroutine test_print_netcdf


        subroutine check_inputs()

          implicit none

          if(.not.(
     $         (nx.eq.20).and.
     $         (ny.eq.25).and.
     $         (ne.eq.4))) then

             print '(''the test requires: '')'
             print '(''nx=20'')'
             print '(''ny=25'')'
             print '(''ne=4'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_bf_mainlayer_print

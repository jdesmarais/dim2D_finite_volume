      program test_bf_layer_path_prog

        use bf_layer_path_class      , only : bf_layer_path
        use bf_sublayer_class        , only : bf_sublayer
        use interface_abstract_class , only : interface_abstract
        use parameters_constant      , only : N,W
        use parameters_input         , only : nx,ny,ne,bc_size
        use parameters_kind          , only : ikind, rkind
        use test_bf_layer_module     , only : print_nodes,
     $                                        print_sizes

        
        implicit none


        type(bf_layer_path)                         :: path_tested
        integer(ikind), dimension(:,:), allocatable :: bc_interior_pt_table
        type(interface_abstract)                    :: interface_used
        real(rkind), dimension(nx,ny,ne)            :: nodes
        real(rkind)                                 :: path_id

        integer :: i,k

        !test of the procedure ini()
        call path_tested%ini()

        print '()'
        print '(''test ini() : initialization of bf_layer_path'')'
        print '(''--------------------------------------------'')'
        print '(''matching_layer associated? '', L1)', associated(path_tested%matching_sublayer)
        print '(''ends? '', L1)', path_tested%ends
        print '(''ends_with_corner? '', L1)', path_tested%ends_with_corner
        print '(''nb pts in path? '', I2)', path_tested%nb_pts
        print '(''neighbors for corresponding layer? '', 4L2)', path_tested%neighbors
        print '()'


        !test of the procedure analyze_pt()
        print '()'
        print '(''test analyze_pt()'')'
        print '(''--------------------------------------------'')'

        !initialization of the nodes
        call ini_nodes(nodes)
        call print_sizes(nodes,'interior_sizes.dat')

        !initialization of the interface
        call ini_interface(interface_used,nodes)

        !initialization of the bc_interior_pt
        !for the paths to be created
        call ini_bc_interior_pt(bc_interior_pt_table)

        !initialize the path
        call path_tested%ini()


        !create the paths        
        path_id = 0.1
        do k=1, size(bc_interior_pt_table,2)

           !analyze the points saved in bc_interior_pt_table
           call path_tested%analyze_pt(
     $          bc_interior_pt_table(:,k),
     $          interface_used)

           !if by any chance the path ends, it should be saved
           !in the nodes table written, reinitialize and the 
           !last grid point analyzed should be reinitialized
           if(path_tested%ends) then
              
              !write the content of the path in the nodes table
              do i=1, path_tested%nb_pts
                 nodes(path_tested%pts(1,i),path_tested%pts(2,i),1) = path_id
              end do

              !print the type of the path
              print '(''path created'')'
              print '(''nb_pts: '', I2)', path_tested%nb_pts
              print '(''alignment: '', 4I3)', path_tested%alignment
              if(associated(path_tested%matching_sublayer)) then
                 print '(''sublayer_matching: associated'')'
                 print '(''  + localization: '', I2)',
     $                path_tested%matching_sublayer%element%localization
                 print '(''  + alignment: '', 4I3)',
     $                path_tested%matching_sublayer%element%alignment
              else
                 print '(''sublayer_matching: non associated'')'
              end if
              print '('''')'


              !reinitialize the path
              call path_tested%reinitialize()

              !analyze the last grid point that led to the path end
              call path_tested%analyze_pt(
     $             bc_interior_pt_table(:,k),
     $             interface_used)
              
              path_id = path_id + 0.1

              !if the last point was a corner then print the corner
              if(path_tested%ends_with_corner) then
                 nodes(bc_interior_pt_table(1,k), bc_interior_pt_table(2,k),1)=path_id
                 path_id = path_id + 0.1
              end if


           end if

        end do

        call print_nodes(nodes,'interior_nodes.dat')

        contains
        
        subroutine ini_nodes(nodes)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

          integer(ikind) :: i,j,k

          do k=1, ne
             do j=1, ny
                do i=1, nx
                   nodes(i,j,k) = 1.0
                end do
             end do
          end do


        end subroutine ini_nodes


        subroutine ini_bc_interior_pt(bc_interior_pt_table)

          implicit none

          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_interior_pt_table

          
          integer :: nb_bc_interior_pt
          
          nb_bc_interior_pt = 16
          allocate(bc_interior_pt_table(2,nb_bc_interior_pt))

          bc_interior_pt_table(:,1)  = [bc_size+1,ny-1]
          bc_interior_pt_table(:,2)  = [bc_size+2,ny-1]
          bc_interior_pt_table(:,3)  = [bc_size+4,ny-1]
          bc_interior_pt_table(:,4)  = [bc_size+7,ny-1]
                                     
          bc_interior_pt_table(:,5)  = [nx-4,ny-1]
          bc_interior_pt_table(:,6)  = [nx-3,ny-1]
          bc_interior_pt_table(:,7)  = [nx-2,ny-1]

          bc_interior_pt_table(:,8)  = [nx-1,10]
          bc_interior_pt_table(:,9)  = [nx-1,7]
          bc_interior_pt_table(:,10) = [nx-1,6]

          bc_interior_pt_table(:,11) = [nx-1,bc_size]

          bc_interior_pt_table(:,12) = [bc_size, 4]
          bc_interior_pt_table(:,13) = [bc_size, 5]
          bc_interior_pt_table(:,14) = [bc_size, 8]
          bc_interior_pt_table(:,15) = [bc_size, 13]
          
          bc_interior_pt_table(:,16) = [bc_size, ny-1]

        end subroutine ini_bc_interior_pt


        subroutine ini_interface(interface_tested, nodes)

          implicit none

          type(interface_abstract)     , intent(inout) :: interface_tested
          real(rkind), dimension(:,:,:), intent(in)    :: nodes

          type(bf_sublayer), pointer :: added_sublayer
          integer, dimension(2,2)    :: alignment
          logical, dimension(4)      :: neighbors

          character(len=22)          :: sizes_filename
          character(len=22)          :: nodes_filename
          character(len=22)          :: grdid_filename

          
          !initialize the interface as the number of sublayers
          call print_nb_sublayers('sublayers_nb.dat',2)
          call interface_tested%ini()


          !add a sublayer to the North main layer
          alignment(1,1) = bc_size + 1
          alignment(1,2) = bc_size + 1
          alignment(2,1) = ny-1
          alignment(2,2) = ny-1

          neighbors = [.false.,.false.,.false.,.false.]

          added_sublayer => interface_tested%add_sublayer(
     $             N, alignment, nodes, neighbors)

          write(sizes_filename,'(A2,I1,''_sizes.dat'')') 'N_',1
          write(nodes_filename,'(A2,I1,''_nodes.dat'')') 'N_',1
          write(grdid_filename,'(A2,I1,''_grdpt_id.dat'')') 'N_',1
          call added_sublayer%element%print_sizes(sizes_filename)
          call added_sublayer%element%print_nodes(nodes_filename)
          call added_sublayer%element%print_grdpts_id(grdid_filename)


          !add two sublayers to the West main layer
          alignment(1,1) = bc_size + 1
          alignment(1,2) = bc_size + 1
          alignment(2,1) = ny-bc_size-4
          alignment(2,2) = ny-bc_size

          neighbors = [.false.,.false.,.false.,.false.]

          added_sublayer => interface_tested%add_sublayer(
     $             W, alignment, nodes, neighbors)      

          write(sizes_filename,'(A2,I1,''_sizes.dat'')') 'W_',1
          write(nodes_filename,'(A2,I1,''_nodes.dat'')') 'W_',1
          write(grdid_filename,'(A2,I1,''_grdpt_id.dat'')') 'W_',1
          call added_sublayer%element%print_sizes(sizes_filename)
          call added_sublayer%element%print_nodes(nodes_filename)
          call added_sublayer%element%print_grdpts_id(grdid_filename)


          alignment(1,1) = bc_size + 1
          alignment(1,2) = bc_size + 1
          alignment(2,1) = ny-bc_size-13
          alignment(2,2) = ny-bc_size-13

          neighbors = [.false.,.false.,.false.,.false.]

          added_sublayer => interface_tested%add_sublayer(
     $             W, alignment, nodes, neighbors)      

          write(sizes_filename,'(A2,I1,''_sizes.dat'')') 'W_',2
          write(nodes_filename,'(A2,I1,''_nodes.dat'')') 'W_',2
          write(grdid_filename,'(A2,I1,''_grdpt_id.dat'')') 'W_',2
          call added_sublayer%element%print_sizes(sizes_filename)
          call added_sublayer%element%print_nodes(nodes_filename)
          call added_sublayer%element%print_grdpts_id(grdid_filename)

        end subroutine ini_interface


        subroutine print_nb_sublayers(filename,nb_sublayers)

          implicit none

          character(*), intent(in) :: filename
          integer     , intent(in) :: nb_sublayers

          integer :: ios
          
          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) nb_sublayers
              close(unit=1)
           else
              stop 'file opening pb'
           end if

        end subroutine print_nb_sublayers

      end program test_bf_layer_path_prog

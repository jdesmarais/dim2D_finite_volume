      program test_bf_path_icr_prog

        use bf_path_icr_class , only : bf_path_icr
        use bf_sublayer_class   , only : bf_sublayer
        use bf_interface_class  , only : bf_interface
        use parameters_bf_layer , only : align_N, align_W
        use parameters_constant , only : N,W
        use parameters_input    , only : nx,ny,ne,bc_size
        use parameters_kind     , only : ikind, rkind
        use test_bf_layer_module, only : print_interior_data,
     $                                   ini_grdpts_id,
     $                                   ini_nodes

        
        implicit none


        type(bf_path_icr)                         :: path_tested
        integer(ikind), dimension(:,:), allocatable :: bc_interior_pt_table
        type(bf_interface)                          :: interface_used
        real(rkind), dimension(nx,ny,ne)            :: nodes
        real(rkind), dimension(nx,ny,ne)            :: nodes_for_path
        integer    , dimension(nx,ny)               :: grdpts_id
        real(rkind)                                 :: path_id
        integer(ikind), dimension(2)                :: pt

        integer :: i,k, file_index


        !call test_create_filenames()

        !test of the procedure ini()
        call path_tested%ini()

        !test of the procedure analyze_pt()
        print '()'
        print '(''test analyze_pt()'')'
        print '(''--------------------------------------------'')'

        !initialization of the nodes
        call ini_nodes(nodes)
        call ini_nodes_for_path(nodes_for_path)
        call ini_grdpts_id(grdpts_id)

        !print interior data
        call print_interior_data(nodes,
     $                           grdpts_id, 
     $                           'interior_nodes1.dat',
     $                           'interior_grdpts_id1.dat',
     $                           'interior_sizes1.dat')        

        !initialization of the interface
        call ini_interface(interface_used)!,nodes)

        !print the initial state
        call interface_used%print_binary(
     $       'nodes1.dat',
     $       'grdpt_id1.dat',
     $       'sizes1.dat',
     $       '1.dat')

        !initialization of the bc_interior_pt
        !for the paths to be created
        call ini_bc_interior_pt(bc_interior_pt_table)

        !initialize the path
        call path_tested%ini()


        !create the paths        
        path_id = 0.1
        file_index=2
        do k=1, size(bc_interior_pt_table,2)

           !analyze the points saved in bc_interior_pt_table
           call path_tested%analyze_pt(
     $          bc_interior_pt_table(:,k),
     $          interface_used)

           !if by any chance the path ends, it should be saved
           !in the nodes table written, reinitialize and the 
           !last grid point analyzed should be reinitialized
           if(path_tested%is_ended()) then
              
              !write the content of the path in the nodes table
              do i=1, path_tested%get_nb_pts()
                 pt = path_tested%get_pt(i)
                 nodes_for_path(pt(1),pt(2),1) = path_id
              end do

              !print the type of the path
              call path_tested%print_on_screen()

c$$$              !print path map + print interior + bf_interface
c$$$              call print_path_map_and_interface(
c$$$     $             nodes_for_path, nodes, grdpts_id, interface_used,
c$$$     $             file_index)

              !process the path
              call path_tested%process_path(interface_used, nodes)

c$$$              !print path map + print interior + bf_interface
c$$$              call print_path_map_and_interface(
c$$$     $             nodes_for_path, nodes, grdpts_id, interface_used,
c$$$     $             file_index)

              !reinitialize the path
              call path_tested%reinitialize()

              !analyze the last grid point that led to the path end
              call path_tested%analyze_pt(
     $             bc_interior_pt_table(:,k),
     $             interface_used)
              
              path_id = path_id + 0.1

           end if

        end do

        !print path map + print interior + bf_interface
        call print_path_map_and_interface(
     $       nodes_for_path, nodes, grdpts_id, interface_used,
     $       file_index)

        print *, file_index

        contains

        subroutine print_path_map_and_interface(
     $       nodes_for_path, nodes, grdpts_id, interface_used, index)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_for_path
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer    , dimension(nx,ny)   , intent(in)    :: grdpts_id
          class(bf_interface)             , intent(in)    :: interface_used
          integer                         , intent(inout) :: index
          

          character(len=40) :: interior_nodes_filename
          character(len=40) :: interior_grdpts_filename
          character(len=40) :: interior_sizes_filename
          character(len=15)  :: bf_nodes_filename
          character(len=15)  :: bf_grdpts_filename
          character(len=15)  :: bf_sizes_filename
          character(len=15)  :: bf_nb_sublayers_filename


          !create the names for the files
          call create_filenames(
     $         index,
     $         interior_nodes_filename,
     $         interior_grdpts_filename,
     $         interior_sizes_filename,
     $         bf_nodes_filename,
     $         bf_grdpts_filename,
     $         bf_sizes_filename,
     $         bf_nb_sublayers_filename)

          !print the nodes for path
          call print_interior_data(nodes_for_path,
     $                             grdpts_id,
     $                             interior_nodes_filename,
     $                             interior_grdpts_filename,
     $                             interior_sizes_filename)

          !print the interface at ths state
          call interface_used%print_binary(
     $         bf_nodes_filename,
     $         bf_grdpts_filename,
     $         bf_sizes_filename,
     $         bf_nb_sublayers_filename)

          index = index+1


          !create the names for the files
          call create_filenames(
     $         index,
     $         interior_nodes_filename,
     $         interior_grdpts_filename,
     $         interior_sizes_filename,
     $         bf_nodes_filename,
     $         bf_grdpts_filename,
     $         bf_sizes_filename,
     $         bf_nb_sublayers_filename)

          !print the interior nodes
          call print_interior_data(nodes,
     $                             grdpts_id,
     $                             interior_nodes_filename,
     $                             interior_grdpts_filename,
     $                             interior_sizes_filename)

          !print the interface at ths state
          call interface_used%print_binary(
     $         bf_nodes_filename,
     $         bf_grdpts_filename,
     $         bf_sizes_filename,
     $         bf_nb_sublayers_filename)

          index = index+1          

        end subroutine print_path_map_and_interface


        subroutine create_filenames(
     $         index,
     $         interior_nodes_filename,
     $         interior_grdpts_filename,
     $         interior_sizes_filename,
     $         bf_nodes_filename,
     $         bf_grdpts_filename,
     $         bf_sizes_filename,
     $         bf_nb_sublayers_filename)

          implicit none

          integer     , intent(in)  :: index
          character(*), intent(out) :: interior_nodes_filename
          character(*), intent(out) :: interior_grdpts_filename
          character(*), intent(out) :: interior_sizes_filename
          character(*), intent(out) :: bf_nodes_filename
          character(*), intent(out) :: bf_grdpts_filename
          character(*), intent(out) :: bf_sizes_filename
          character(*), intent(out) :: bf_nb_sublayers_filename
          
          character(len=11) :: interior_nodes_format
          character(len=11) :: interior_grdpts_format
          character(len=10) :: bf_nodes_format
          character(len=10) :: bf_grdpts_format
          character(len=7) :: bf_nb_sublayers_format


          !write formats
          write(interior_nodes_format,
     $         '(''(A14,I'',I1,'',A4)'')')
     $         floor(real(index)/10.0)+1
          write(interior_grdpts_format,
     $         '(''(A18,I'',I1,'',A4)'')')
     $         floor(real(index)/10.0)+1
          write(bf_nodes_format,
     $         '(''(A5,I'',I1,'',A4)'')')
     $         floor(real(index)/10.0)+1
          write(bf_grdpts_format,
     $         '(''(A8,I'',I1,'',A4)'')')
     $         floor(real(index)/10.0)+1
          write(bf_nb_sublayers_format,
     $         '(''(I'',I1,'',A4)'')')
     $         floor(real(index)/10.0)+1


          !write filenames
          write(interior_nodes_filename, interior_nodes_format)
     $         'interior_nodes',
     $         index,
     $         '.dat'

          write(interior_grdpts_filename, interior_grdpts_format)
     $         'interior_grdpts_id',
     $         index,
     $         '.dat'

          write(interior_sizes_filename, interior_nodes_format)
     $         'interior_sizes',
     $         index,
     $         '.dat'

          write(bf_nodes_filename, bf_nodes_format)
     $         'nodes',
     $         index,
     $         '.dat'

          write(bf_grdpts_filename, bf_grdpts_format)
     $         'grdpt_id',
     $         index,
     $         '.dat'

          write(bf_sizes_filename, bf_nodes_format)
     $         'sizes',
     $         index,
     $         '.dat'

          write(bf_nb_sublayers_filename, bf_nb_sublayers_format)
     $         index,
     $         '.dat'

        end subroutine create_filenames


        subroutine test_create_filenames()

          implicit none

          integer      :: i
          character(len=30) :: interior_nodes_filename
          character(len=30) :: interior_grdpts_filename
          character(len=30) :: interior_sizes_filename
          character(len=30) :: bf_nodes_filename
          character(len=30) :: bf_grdpts_filename
          character(len=30) :: bf_sizes_filename
          character(len=30) :: bf_nb_sublayers_filename

          do i=1, 10

             call create_filenames(
     $            i,
     $            interior_nodes_filename,
     $            interior_grdpts_filename,
     $            interior_sizes_filename,
     $            bf_nodes_filename,
     $            bf_grdpts_filename,
     $            bf_sizes_filename,
     $            bf_nb_sublayers_filename)
             
             print *, interior_nodes_filename
             print *, interior_grdpts_filename
             print *, interior_sizes_filename
             print *, bf_nodes_filename
             print *, bf_grdpts_filename
             print *, bf_sizes_filename
             print *, bf_nb_sublayers_filename

          end do

        end subroutine test_create_filenames


        subroutine ini_bc_interior_pt(bc_interior_pt_table)

          implicit none

          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_interior_pt_table

          
          integer :: nb_bc_interior_pt
          
          nb_bc_interior_pt = 27
          allocate(bc_interior_pt_table(2,nb_bc_interior_pt))

          bc_interior_pt_table(:,1)  = [bc_size,ny-2]
          bc_interior_pt_table(:,2)  = [bc_size,ny-1]
          bc_interior_pt_table(:,3)  = [bc_size+1,ny-1]
          bc_interior_pt_table(:,4)  = [bc_size+2,ny-1]
          bc_interior_pt_table(:,5)  = [bc_size+4,ny-1]
          bc_interior_pt_table(:,6)  = [bc_size+7,ny-1]
                                     
          bc_interior_pt_table(:,7)  = [nx-4,ny-1]
          bc_interior_pt_table(:,8)  = [nx-3,ny-1]
          bc_interior_pt_table(:,9)  = [nx-2,ny-1]
          bc_interior_pt_table(:,10)  = [nx-1,ny-1]
          bc_interior_pt_table(:,11)  = [nx-1,ny-2]

          bc_interior_pt_table(:,12)  = [nx-1,ny-2]
          bc_interior_pt_table(:,13)  = [nx-1,ny-3]
          bc_interior_pt_table(:,14) = [nx-1,ny-4]

          bc_interior_pt_table(:,15) = [nx-1,10]
          bc_interior_pt_table(:,16) = [nx-1,7]
          bc_interior_pt_table(:,17) = [nx-1,6]

          bc_interior_pt_table(:,18) = [nx-1,bc_size+1]
          bc_interior_pt_table(:,19) = [nx-1,bc_size]
          bc_interior_pt_table(:,20) = [nx-2,bc_size]
          bc_interior_pt_table(:,21) = [nx-3,bc_size]

          bc_interior_pt_table(:,22) = [bc_size+1, bc_size]
          bc_interior_pt_table(:,23) = [bc_size, bc_size]
          bc_interior_pt_table(:,24) = [bc_size, bc_size+1]
          bc_interior_pt_table(:,25) = [bc_size, bc_size+2]
          bc_interior_pt_table(:,26) = [bc_size, 8]
          bc_interior_pt_table(:,27) = [bc_size, 13]
         
        end subroutine ini_bc_interior_pt


        subroutine ini_interface(interface_tested)!, nodes)

          implicit none

          type(bf_interface)           , intent(inout) :: interface_tested
c$$$          real(rkind), dimension(:,:,:), intent(in)    :: nodes

c$$$          type(bf_sublayer), pointer :: added_sublayer
c$$$          integer, dimension(2,2)    :: alignment
          

          !initialize the interface as the number of sublayers
          call interface_tested%ini()


c$$$          !add a sublayer to the North main layer
c$$$          alignment(1,1) = bc_size + 1
c$$$          alignment(1,2) = bc_size + 1
c$$$          alignment(2,1) = ny-1
c$$$          alignment(2,2) = ny-1
c$$$
c$$$          added_sublayer => interface_tested%allocate_sublayer(
c$$$     $             N, nodes, alignment)
c$$$
c$$$
c$$$          !add two sublayers to the West main layer
c$$$          alignment(1,1) = align_W
c$$$          alignment(1,2) = align_W
c$$$          alignment(2,1) = ny-bc_size-4
c$$$          alignment(2,2) = align_N-1
c$$$
c$$$          added_sublayer => interface_tested%allocate_sublayer(
c$$$     $             W, nodes, alignment)
c$$$
c$$$
c$$$          alignment(1,1) = align_W
c$$$          alignment(1,2) = align_W
c$$$          alignment(2,1) = ny-bc_size-13
c$$$          alignment(2,2) = ny-bc_size-13
c$$$
c$$$          added_sublayer => interface_tested%allocate_sublayer(
c$$$     $             W, nodes, alignment)

        end subroutine ini_interface

      
        subroutine ini_nodes_for_path(nodes)

          implicit none

          real(rkind), dimension(:,:,:), intent(out) :: nodes

          integer(ikind) :: i,j,k


          do k=1, size(nodes,3)
             do j=1, size(nodes,2)
                do i=1, size(nodes,1)
                   nodes(i,j,k) = 1.0
                end do
             end do
          end do

        end subroutine ini_nodes_for_path

      end program test_bf_path_icr_prog

      program test_bf_layer_path_prog

        use bf_layer_path_class , only : bf_layer_path
        use bf_sublayer_class   , only : bf_sublayer
        use bf_interface_class  , only : bf_interface
        use parameters_constant , only : N,W
        use parameters_input    , only : nx,ny,ne,bc_size
        use parameters_kind     , only : ikind, rkind
        use test_bf_layer_module, only : print_interior_data,
     $                                   ini_grdpts_id

        
        implicit none


        type(bf_layer_path)                         :: path_tested
        integer(ikind), dimension(:,:), allocatable :: bc_interior_pt_table
        type(bf_interface)                          :: interface_used
        real(rkind), dimension(nx,ny,ne)            :: nodes
        integer    , dimension(nx,ny)               :: grdpts_id
        real(rkind)                                 :: path_id
        integer(ikind), dimension(2)                :: pt

        integer :: i,k

        !test of the procedure ini()
        call path_tested%ini()

        !test of the procedure analyze_pt()
        print '()'
        print '(''test analyze_pt()'')'
        print '(''--------------------------------------------'')'

        !initialization of the nodes
        call ini_nodes(nodes)
        call ini_grdpts_id(grdpts_id)
        

        !initialization of the interface
        call ini_interface(interface_used,nodes)

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
                 nodes(pt(1),pt(2),1) = path_id
              end do

              !print the type of the path
              call path_tested%print_on_screen()

              !reinitialize the path
              call path_tested%reinitialize()

              !analyze the last grid point that led to the path end
              call path_tested%analyze_pt(
     $             bc_interior_pt_table(:,k),
     $             interface_used)
              
              path_id = path_id + 0.1

           end if

        end do

        call print_interior_data(nodes,
     $                           grdpts_id, 
     $                           'interior_nodes.dat',
     $                           'interior_grdpts_id.dat',
     $                           'interior_sizes.dat')


        contains


        subroutine ini_bc_interior_pt(bc_interior_pt_table)

          implicit none

          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_interior_pt_table

          
          integer :: nb_bc_interior_pt
          
          nb_bc_interior_pt = 19
          allocate(bc_interior_pt_table(2,nb_bc_interior_pt))

          bc_interior_pt_table(:,1)  = [bc_size+1,ny-1]
          bc_interior_pt_table(:,2)  = [bc_size+2,ny-1]
          bc_interior_pt_table(:,3)  = [bc_size+4,ny-1]
          bc_interior_pt_table(:,4)  = [bc_size+7,ny-1]
                                     
          bc_interior_pt_table(:,5)  = [nx-4,ny-1]
          bc_interior_pt_table(:,6)  = [nx-3,ny-1]
          bc_interior_pt_table(:,7)  = [nx-2,ny-1]

          bc_interior_pt_table(:,8)  = [nx-1,ny-2]
          bc_interior_pt_table(:,9)  = [nx-1,ny-3]
          bc_interior_pt_table(:,10) = [nx-1,ny-4]

          bc_interior_pt_table(:,11) = [nx-1,10]
          bc_interior_pt_table(:,12) = [nx-1,7]
          bc_interior_pt_table(:,13) = [nx-1,6]

          bc_interior_pt_table(:,14) = [nx-1,bc_size]

          bc_interior_pt_table(:,15) = [bc_size, 4]
          bc_interior_pt_table(:,16) = [bc_size, 5]
          bc_interior_pt_table(:,17) = [bc_size, 8]
          bc_interior_pt_table(:,18) = [bc_size, 13]
         
          bc_interior_pt_table(:,19) = [bc_size, ny-1]

        end subroutine ini_bc_interior_pt


        subroutine ini_interface(interface_tested, nodes)

          implicit none

          type(bf_interface)           , intent(inout) :: interface_tested
          real(rkind), dimension(:,:,:), intent(in)    :: nodes

          type(bf_sublayer), pointer :: added_sublayer
          integer, dimension(2,2)    :: alignment
          

          !initialize the interface as the number of sublayers
          call interface_tested%ini()


          !add a sublayer to the North main layer
          alignment(1,1) = bc_size + 1
          alignment(1,2) = bc_size + 1
          alignment(2,1) = ny-1
          alignment(2,2) = ny-1

          added_sublayer => interface_tested%allocate_sublayer(
     $             N, nodes, alignment)


          !add two sublayers to the West main layer
          alignment(1,1) = bc_size + 1
          alignment(1,2) = bc_size + 1
          alignment(2,1) = ny-bc_size-4
          alignment(2,2) = ny-bc_size

          added_sublayer => interface_tested%allocate_sublayer(
     $             W, nodes, alignment)


          alignment(1,1) = bc_size + 1
          alignment(1,2) = bc_size + 1
          alignment(2,1) = ny-bc_size-13
          alignment(2,2) = ny-bc_size-13

          added_sublayer => interface_tested%allocate_sublayer(
     $             W, nodes, alignment)

        end subroutine ini_interface

      
        subroutine ini_nodes(nodes)

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

        end subroutine ini_nodes

      end program test_bf_layer_path_prog

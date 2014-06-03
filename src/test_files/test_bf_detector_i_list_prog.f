      program test_bf_detctor_i_list_prog

        use bf_detector_i_list_class, only : bf_detector_i_list
        use bf_interface_class      , only : bf_interface
        use parameters_bf_layer     , only : no_pt, interior_pt
        use parameters_constant     , only : N,S,E,W
        use parameters_input        , only : nx,ny,ne,bc_size
        use parameters_kind         , only : ikind, rkind
        use test_bf_layer_module    , only : print_interior_data

        implicit none

        type(bf_detector_i_list)                    :: detector_list
        type(bf_interface)                          :: interface_used
        real(rkind), dimension(nx,ny,ne)            :: matrix
        integer    , dimension(nx,ny)               :: grdpts_id
        integer(ikind), dimension(:,:), allocatable :: detectors_added
        integer(ikind), dimension(:,:), allocatable :: new_detectors

        integer :: mainlayer_id
        integer :: size_preallocated
        integer :: k

        character(len=30) :: nodes_filename
        character(len=30) :: grdpts_id_filename
        character(len=30) :: sizes_filename

        logical :: print_on_matrix =.false.

        do mainlayer_id=1,4

           !initialize the matrix where the position of the
           !detectors is pinned
           call ini(matrix, grdpts_id)        
           
           !initialize the bf_detector_i_list for the
           !main layer and 6 detectors preallocated
           size_preallocated = 6
           call detector_list%ini(mainlayer_id, size_preallocated)
           
           
           !initialize the detectors saved in the bf_detector_i_list
           call get_detector_test(mainlayer_id, detectors_added)
           do k=1, size(detectors_added,2)
              
              !add a new detector to the object saving the detector
              !coordinates
              call detector_list%add_new_detector(
     $             detectors_added(:,k))
           
              !pin-point the detectors asked
              grdpts_id(detectors_added(1,k),detectors_added(2,k)) =
     $             interior_pt
           end do
           
           if(print_on_matrix) then

              !print the detectors saved by the object bf_detector_i_list
              call detector_list%print_on_matrix(matrix(:,:,1))

           else

              !put the new detectors in another table
              allocate(new_detectors(2,detector_list%get_nb_detectors()))
              call detector_list%fill_new_detector_table(1,new_detectors)
              
              do k=1, size(new_detectors,2)
                 matrix(new_detectors(1,k), new_detectors(2,k),1) = 0.3d0
              end do

           end if
           

           !write the file names
           write(nodes_filename    , '(''interior_nodes'',I1,''.dat'')')
     $          mainlayer_id
           write(grdpts_id_filename, '(''interior_grdpts_id'',I1,''.dat'')')
     $          mainlayer_id
           write(sizes_filename    , '(''interior_sizes'',I1,''.dat'')')
     $          mainlayer_id


           !print the data
           call print_interior_data(matrix,
     $                              grdpts_id, 
     $                              nodes_filename,
     $                              grdpts_id_filename,
     $                              sizes_filename)

           !deallocate the grdpts
           deallocate(detectors_added)
           call detector_list%destroy()

           !deallocate the new detectors
           if(allocated(new_detectors)) then
              deallocate(new_detectors)
           end if

        end do


        contains


        subroutine get_detector_test(mainlayer_id, detectors_added)

          implicit none

          integer                             , intent(in) :: mainlayer_id
          integer, dimension(:,:), allocatable, intent(out) :: detectors_added


          select case(mainlayer_id)
            case(N)
               allocate(detectors_added(2,10))
               detectors_added(:,1)  = [bc_size, ny-8]
               detectors_added(:,2)  = [bc_size, ny-5]
               detectors_added(:,3)  = [bc_size, ny-3]
                                     
               detectors_added(:,4)  = [bc_size  , ny-1]
               detectors_added(:,5)  = [bc_size+5, ny-1]
               detectors_added(:,6)  = [bc_size+8, ny-1]
               detectors_added(:,7)  = [bc_size+12, ny-1]

               detectors_added(:,8)  = [nx-1,ny-1]
               detectors_added(:,9)  = [nx-1,ny-3]
               detectors_added(:,10) = [nx-1,ny-7]
               
            case(S)
               allocate(detectors_added(2,10))
               detectors_added(:,1)  = [bc_size, bc_size+8]
               detectors_added(:,2)  = [bc_size, bc_size+3]
               detectors_added(:,3)  = [bc_size, bc_size+1]
                                     
               detectors_added(:,4)  = [bc_size+3,  bc_size]
               detectors_added(:,5)  = [bc_size+5,  bc_size]
               detectors_added(:,6)  = [bc_size+8,  bc_size]
               detectors_added(:,7)  = [bc_size+12, bc_size]

               detectors_added(:,8)  = [nx-1, bc_size+2]
               detectors_added(:,9)  = [nx-1, bc_size+6]
               detectors_added(:,10) = [nx-1, bc_size+10]

            case(E)
               allocate(detectors_added(2,10))
               detectors_added(:,1)  = [nx-8, bc_size]
               detectors_added(:,2)  = [nx-5, bc_size]
               detectors_added(:,3)  = [nx-2, bc_size]

               detectors_added(:,4)  = [nx-1, bc_size+3]
               detectors_added(:,5)  = [nx-1, bc_size+5]
               detectors_added(:,6)  = [nx-1, bc_size+6]
               detectors_added(:,7)  = [nx-1, bc_size+10]

               detectors_added(:,8)  = [nx-1, ny-1]
               detectors_added(:,9)  = [nx-4, ny-1]
               detectors_added(:,10) = [nx-6, ny-1]

            case(W)
               allocate(detectors_added(2,10))
               detectors_added(:,1)  = [bc_size+8, bc_size]
               detectors_added(:,2)  = [bc_size+5, bc_size]
               detectors_added(:,3)  = [bc_size+2, bc_size]

               detectors_added(:,4)  = [bc_size, bc_size+3]
               detectors_added(:,5)  = [bc_size, bc_size+5]
               detectors_added(:,6)  = [bc_size, bc_size+6]
               detectors_added(:,7)  = [bc_size, bc_size+10]

               detectors_added(:,8)  = [bc_size+1, ny-1]
               detectors_added(:,9)  = [bc_size+4, ny-1]
               detectors_added(:,10) = [bc_size+6, ny-1]
          end select            

        end subroutine get_detector_test

        subroutine ini(matrix, grdpts_id)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: matrix
          integer    , dimension(nx,ny)   , intent(out) :: grdpts_id

          integer(ikind) :: i,j

          do j=1, ny
             do i=1, nx
                matrix(i,j,1)  = 1.0
                grdpts_id(i,j) = no_pt
             end do
          end do

        end subroutine ini

      end program test_bf_detctor_i_list_prog

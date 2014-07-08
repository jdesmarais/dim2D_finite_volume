      program test_bf_interface_icr_prog

        use bf_detector_icr_list_class, only : bf_detector_icr_list
        use bf_interface_icr_class  , only : bf_interface_icr
        use bf_sublayer_class       , only : bf_sublayer
        use parameters_bf_layer     , only : align_N, interior_pt, no_pt
        use parameters_constant     , only : N,S,E,W
        use parameters_input        , only : nx,ny,ne,bc_size
        use parameters_kind         , only : ikind, rkind
        use test_bf_layer_module    , only : print_interior_data,
     $                                       ini_grdpts_id

        implicit none


        type(bf_interface_icr)              :: interface_used
        real(rkind)   , dimension(nx,ny,ne) :: nodes
        integer       , dimension(nx,ny)    :: grdpts_id
        integer(ikind), dimension(7,7,2)    :: cpt_coords_tab
        integer(ikind), dimension(2)        :: cpt_coords
        integer(ikind), dimension(2)        :: cpt_coords_p
        integer       , dimension(3,3)      :: nbc_template
        integer :: i,j,k
        integer :: nb_mgrdpts
        integer, dimension(2,10) :: mgrdpts
        integer :: mgrdpts_i, mgrdpts_j
        integer, parameter :: grdpt_checked = no_pt
        type(bf_sublayer), pointer :: added_sublayer

        !for the recombination test
        type(bf_detector_icr_list), dimension(4) :: detector_list
        integer(ikind), dimension(:,:), allocatable :: detectors_added
        integer :: size_preallocated
        integer :: mainlayer_id


        !initialization of nodes and grdpts_id
        call ini_cst_nodes(nodes)
        call ini_grdpts_id(grdpts_id)

        !test of ini()
        call interface_used%ini()
        call initialize_sublayers_in_interface(
     $       interface_used, nodes, added_sublayer)

        call interface_used%print_idetectors_on(nodes(:,:,1))

        !print the detector positions
        call print_interior_data(nodes,
     $                           grdpts_id,
     $                           'interior_nodes1.dat',
     $                           'interior_grdpts_id1.dat',
     $                           'interior_sizes1.dat')

        !print the interface
        call interface_used%print_binary(
     $       'nodes1.dat',
     $       'grdpt_id1.dat',
     $       'sizes1.dat',
     $       '1.dat')

        !test of create_nbc_interior_pt_template()
        !and check_nbc_interior_pt_template
        call ini_cpt_coords_table(cpt_coords_tab)
        do j=1,7
           do i=1,7

              if((i.ne.4).and.(j.ne.4)) then

                 !test of create_nbc_interior_pt_template()
                 cpt_coords = cpt_coords_tab(i,j,:)
                 nbc_template = interface_used%create_nbc_interior_pt_template(
     $                cpt_coords)
                 !call print_nbc_template(i,j,nbc_template)
   
                 !test of check_nbc_interior_pt_template()
                 
                 nb_mgrdpts = 0
                 call interface_used%check_nbc_interior_pt_template(
     $                nbc_template,
     $                nx/2, ny/2,
     $                cpt_coords(1), cpt_coords(2),
     $                nb_mgrdpts,
     $                mgrdpts)

                 do k=1, nb_mgrdpts
                    mgrdpts_i = mgrdpts(1,k)
                    mgrdpts_j = mgrdpts(2,k)
                    grdpts_id(mgrdpts_i,mgrdpts_j) = grdpt_checked
                 end do

              end if

           end do
        end do

        !print the grdpt_checked
        call print_interior_data(nodes,
     $                           grdpts_id,
     $                           'interior_nodes2.dat',
     $                           'interior_grdpts_id2.dat',
     $                           'interior_sizes2.dat')

        !test the check_nbc_interior
        call ini_grdpts_id(grdpts_id)
        cpt_coords = cpt_coords_tab(2,7,:)
        nb_mgrdpts = 0
        call added_sublayer%check_neighboring_bc_interior_pts(
     $       nx/2, ny/2,
     $       cpt_coords(1), cpt_coords(2),
     $       nb_mgrdpts,
     $       mgrdpts)
        do k=1, nb_mgrdpts
           mgrdpts_i = mgrdpts(1,k)
           mgrdpts_j = mgrdpts(2,k)
           grdpts_id(mgrdpts_i,mgrdpts_j) = grdpt_checked
        end do

        !print the grid point tested
        call print_interior_data(nodes,
     $                           grdpts_id,
     $                           'interior_nodes3.dat',
     $                           'interior_grdpts_id3.dat',
     $                           'interior_sizes3.dat')


        !test the check_nbc_interior
        call ini_grdpts_id(grdpts_id)
        cpt_coords_p = cpt_coords_tab(1,7,:)
        cpt_coords   = cpt_coords_tab(2,7,:)
        nb_mgrdpts   = 0
        call added_sublayer%check_neighboring_bc_interior_pts(
     $       cpt_coords_p(1), cpt_coords_p(2),
     $       cpt_coords(1), cpt_coords(2),
     $       nb_mgrdpts,
     $       mgrdpts)
        do k=1, nb_mgrdpts
           mgrdpts_i = mgrdpts(1,k)
           mgrdpts_j = mgrdpts(2,k)
           grdpts_id(mgrdpts_i,mgrdpts_j) = grdpt_checked
        end do

        !print the grid point tested
        call print_interior_data(nodes,
     $                           grdpts_id,
     $                           'interior_nodes4.dat',
     $                           'interior_grdpts_id4.dat',
     $                           'interior_sizes4.dat')


        !test the recombination of detector lists

        !initialize the matrix where the position of the
        !detectors is pinned
        do j=1, ny
           do i=1, nx
              nodes(i,j,1)   = 1.0
              grdpts_id(i,j) = no_pt
           end do
        end do

        !1) initialize the bf_detector_lists that will be
        !   recombined
        do mainlayer_id=1,4
           
           !initialize the bf_detector_icr_list for the
           !main layer and 6 detectors preallocated
           size_preallocated = 6
           call detector_list(mainlayer_id)%ini(
     $          mainlayer_id, size_preallocated)
           
           
           !initialize the detectors saved in the bf_detector_icr_list
           call get_detector_test(mainlayer_id, detectors_added)
           do k=1, size(detectors_added,2)
              
              !add a new detector to the object saving the detector
              !coordinates
              call detector_list(mainlayer_id)%add_new_detector(
     $             detectors_added(:,k))

              grdpts_id(detectors_added(1,k), detectors_added(2,k)) =
     $             interior_pt

           end do

        end do

        !2) print the lists
        do mainlayer_id=1,4
           call detector_list(mainlayer_id)%print_on_matrix(
     $          nodes(:,:,1))
        end do
        call print_interior_data(nodes,
     $                           grdpts_id,
     $                           'interior_nodes5.dat',
     $                           'interior_grdpts_id5.dat',
     $                           'interior_sizes5.dat')

        !3) recombine the lists
        call interface_used%combine_bf_idetector_lists(
     $       detector_list(1), detector_list(2),
     $       detector_list(3), detector_list(4))           

        !4) write the new detector position after recombination\
        !   on the nodes table
        do j=1, ny
           do i=1, nx
              nodes(i,j,1)   = 1.0
           end do
        end do
        call interface_used%print_idetectors_on(nodes(:,:,1))

        !5) print the position of the new detectors
        call print_interior_data(nodes,
     $                           grdpts_id,
     $                           'interior_nodes6.dat',
     $                           'interior_grdpts_id6.dat',
     $                           'interior_sizes6.dat')

        
        contains

        !< initialize the position of the detectors that will be
        !> added to the detector_lists
        subroutine get_detector_test(mainlayer_id, detectors_added)

          implicit none

          integer                             , intent(in)  :: mainlayer_id
          integer, dimension(:,:), allocatable, intent(out) :: detectors_added

          allocate(detectors_added(2,5))

          select case(mainlayer_id)
            case(N)

               detectors_added(:,1)  = [bc_size  ,  ny-1]
               detectors_added(:,2)  = [bc_size+5,  ny-1]
               detectors_added(:,3)  = [bc_size+8,  ny-1]
               detectors_added(:,4)  = [bc_size+12, ny-1]
               detectors_added(:,5)  = [bc_size+14, ny]
               
            case(S)

               detectors_added(:,1)  = [bc_size+3,  bc_size]
               detectors_added(:,2)  = [bc_size+5,  bc_size]
               detectors_added(:,3)  = [bc_size+8,  bc_size]
               detectors_added(:,4)  = [bc_size+12, bc_size]
               detectors_added(:,5)  = [bc_size+14, 1]

            case(E)

               detectors_added(:,1)  = [nx-1, bc_size+3]
               detectors_added(:,2)  = [nx-1, bc_size+5]
               detectors_added(:,3)  = [nx-1, bc_size+6]
               detectors_added(:,4)  = [nx-1, bc_size+10]
               detectors_added(:,5)  = [nx  , bc_size+14]

            case(W)

               detectors_added(:,1)  = [bc_size, bc_size+3]
               detectors_added(:,2)  = [bc_size, bc_size+5]
               detectors_added(:,3)  = [bc_size, bc_size+6]
               detectors_added(:,4)  = [bc_size, bc_size+10]
               detectors_added(:,5)  = [1      , bc_size+14]

          end select            

        end subroutine get_detector_test


        !< initialize the sublayers for the interface
        subroutine initialize_sublayers_in_interface(
     $       interface_used, nodes, added_sublayer)

          implicit none

          type(bf_interface_icr)          , intent(inout) :: interface_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          type(bf_sublayer), pointer      , intent(out)   :: added_sublayer

          
          integer(ikind), dimension(2,2) :: alignment


          alignment(1,1) = bc_size+1
          alignment(1,2) = bc_size+4
          alignment(2,1) = align_N
          alignment(2,2) = align_N


          added_sublayer => interface_used%allocate_sublayer(
     $         N, nodes, alignment)

          !modify the sublayer to see if the copy works
          call added_sublayer%set_grdpts_id_pt(bc_size+1,2*bc_size+1,interior_pt)

        end subroutine initialize_sublayers_in_interface



        !< print the content of the matrix nbc_template
        !> on an output file
        subroutine print_nbc_template(i,j,nbc_template)

          implicit none

          integer                , intent(in) :: i
          integer                , intent(in) :: j
          integer, dimension(3,3), intent(in) :: nbc_template


          character(len=20) :: nbc_filename
          integer :: ios


          write(nbc_filename,
     $         '(''nbc_template'',I1,I1,''.dat'')') i,j

          open(unit=1,
     $          file=nbc_filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) nbc_template
              close(unit=1)
           else
              stop 'file opening pb'
           end if


        end subroutine print_nbc_template


        !< initialize the coordinates tested as cpt_coords
        subroutine ini_cpt_coords_table(cpt_coords)
        
          implicit none

          integer(ikind), dimension(7,7,2), intent(out) :: cpt_coords

          integer(ikind), dimension(7) :: x_coords
          integer(ikind), dimension(7) :: y_coords
          integer :: i,j

          x_coords = [ 1, 2, 3, 4, nx-2, nx-1, nx]
          y_coords = [ 1, 2, 3, 4, ny-2, ny-1, ny]

          do j=1,7
             do i=1,7
                cpt_coords(i,j,1) = x_coords(i)
                cpt_coords(i,j,2) = y_coords(j)
             end do
          end do

        end subroutine ini_cpt_coords_table


        !< initialize the nodes using a constant variable
        subroutine ini_cst_nodes(nodes)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes

          integer(ikind) :: i,j

          do j=1, ny
             do i=1, nx
                nodes(i,j,1)  = 1.0
             end do
          end do

        end subroutine ini_cst_nodes


      end program test_bf_interface_icr_prog

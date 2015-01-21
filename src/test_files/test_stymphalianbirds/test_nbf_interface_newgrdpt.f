      program test_nbf_interface

        use bf_sublayer_class, only :
     $       bf_sublayer

        use nbf_interface_newgrdpt_class, only :
     $       nbf_interface_newgrdpt,
     $       finalize_grdpts_around_new_interior_pt,
     $       finalize_grdpts_for_bc_pt_crenel

        use parameters_constant, only :
     $       N,E,S

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use wave2d_parameters, only :
     $       c


        implicit none


        logical :: detailled
        logical :: test_validated

        detailled = .true.

        
        !test compute_newgrdpt
        !------------------------------------------------------------
        
        !all the data needed for the computation of the new grid point
        !are located in the buffer layer
        test_validated = test_compute_newgrdpt(1, detailled)
        print '(''test_compute_newgrdpt - config 1: '',L1)', test_validated
        print '()'

        !the data needed for the computation of the new grid point need to
        !be extracted from the interior
        test_validated = test_compute_newgrdpt(2, detailled)
        print '(''test_compute_newgrdpt - config 2: '',L1)', test_validated
        print '()'

        !the data needed for the computation of the new grid point need to
        !be extracted from the interior and the neighboring buffer layers
        test_validated = test_compute_newgrdpt(3, detailled)
        print '(''test_compute_newgrdpt - config 3: '',L1)', test_validated
        print '()'
        
        !the buffer layer has just been created and no data are stored
        !at t=t-dt
        test_validated = test_compute_newgrdpt(4, detailled)
        print '(''test_compute_newgrdpt - config 4: '',L1)', test_validated
        print '()'

        !test update_bf_grdpts_after_increase
        !------------------------------------------------------------
        test_validated = test_update_bf_grdpts_after_increase(detailled)
        print '(''test_update_bf_grdpts_after_increase: '',L1)', test_validated
        print '()'

        !test finalize_grdpts_around_new_interior_pt
        !------------------------------------------------------------
        test_validated = test_finalize_grdpts_around_new_interior_pt(1,detailled)
        print '(''test_finalize_grdpts_around_new_interior_pt-1: '',L1)', test_validated
        print '()'

        test_validated = test_finalize_grdpts_around_new_interior_pt(2,detailled)
        print '(''test_finalize_grdpts_around_new_interior_pt-2: '',L1)', test_validated
        print '()'

        !test finalize_grdpts_for_bc_pt_crenel
        !------------------------------------------------------------
        test_validated = test_finalize_grdpts_for_bc_pt_crenel(1,detailled)
        print '(''test_finalize_grdpts_for_bc_pt_crenel-1: '',L1)', test_validated
        print '()'

        test_validated = test_finalize_grdpts_for_bc_pt_crenel(2,detailled)
        print '(''test_finalize_grdpts_for_bc_pt_crenel-2: '',L1)', test_validated
        print '()'


        contains


        function test_compute_newgrdpt(config, detailled)
     $       result(test_validated)

          implicit none
          
          integer, intent(in) :: config
          logical, intent(in) :: detailled
          logical             :: test_validated

          type(nbf_interface_newgrdpt)     :: nbf_interface_used
          type(bf_sublayer), pointer       :: bf_sublayer_used

          type(pmodel_eq)                  :: p_model
          real(rkind)                      :: t
          real(rkind)                      :: dt
          real(rkind)                      :: dx
          real(rkind)                      :: dy
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne) :: interior_nodes1
          integer(ikind)                   :: i1
          integer(ikind)                   :: j1
          real(rkind), dimension(ne)       :: newgrdpt_data
          real(rkind), dimension(ne)       :: newgrdpt
          logical                          :: test_loc

          integer :: k

          !test requirements
          if((nx.ne.10).or.(ny.ne.10).or.(ne.ne.3).or.
     $         (.not.(is_test_validated(c,0.5d0,.false.)))) then

             print '(''test_nbf_interface'')'
             print '(''test_compute_newgrdpt'')'
             print '(''the test requires:'')'
             print '(''nx=10 : '',L1)', nx.eq.10
             print '(''ny=10 : '',L1)', ny.eq.10
             print '(''ne= 3 : '',L1)', ne.eq.3
             print '('' -> the wave2d model should be used'')'
             print '(''c=0.5 : '',L1)', is_test_validated(c,0.5d0,.false.)
             print '()'

             stop ''

          end if

          test_validated = .true.

          allocate(bf_sublayer_used)


          !initialize the interior x_map and y_map
          select case(config)
             case(1)
                dx = 0.25d0
                dy = 0.25d0
             case(2)
                dx = 1.0d0
                dy = 0.25d0
             case(3)
                dx = 1.0d0
                dy = 0.25d0
             case(4)
                dx = 0.25d0
                dy = 0.25d0                
             case default
                print '(''test_nbf_interface'')'
                print '(''test_compute_newgrdpt'')'
                print '(''config not implemented: '',I2)',config
                stop ''
          end select

          do k=1, nx
             interior_x_map(k) = (k-1)*dx
          end do

          do k=1, ny
             interior_y_map(k) = (k-1)*dy
          end do


          !initialize the buffer layer(s), the interior
          !nodes and the links in the nbf_interface for
          !the test
          call initialize_data_for_test_compute_newgrpdt(
     $         config,
     $         dx,dy,
     $         interior_nodes0,
     $         interior_nodes1,
     $         bf_sublayer_used,
     $         nbf_interface_used,
     $         i1,j1,
     $         newgrdpt_data)
          
          
          !test the computation of the new grdpt
          ! - using only a buffer layer
          ! - using the temporary tables
          t  = 0.0d0
          dt = 0.25d0
          call nbf_interface_used%compute_newgrdpt(
     $         bf_sublayer_used,
     $         p_model, t,dt,
     $         i1,j1,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)

          !compare the results with what is expected
          newgrdpt = bf_sublayer_used%get_nodes([i1,j1])

          do k=1, ne
             test_loc = is_test_validated(
     $            newgrdpt(k), newgrdpt_data(k), detailled)
             test_validated = test_validated.and.test_loc
          end do

        end function test_compute_newgrdpt


        subroutine initialize_data_for_test_compute_newgrpdt(
     $     config,
     $     dx,dy,
     $     interior_nodes0,
     $     interior_nodes1,
     $     bf_sublayer_used,
     $     nbf_interface_used,
     $     i1,j1,
     $     newgrdpt_data)

          implicit none

          integer                           , intent(in)    :: config
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          real(rkind)  , dimension(nx,ny,ne), intent(inout) :: interior_nodes0
          real(rkind)  , dimension(nx,ny,ne), intent(inout) :: interior_nodes1
          type(bf_sublayer)  , pointer      , intent(inout) :: bf_sublayer_used
          type(nbf_interface_newgrdpt)      , intent(out)   :: nbf_interface_used
          integer(ikind)                    , intent(out)   :: i1
          integer(ikind)                    , intent(out)   :: j1
          real(rkind)        , dimension(ne), intent(out)   :: newgrdpt_data


          select case(config)
            case(1)
               call initialize_data_for_test_compute_newgrpdt_config1(
     $              dx,dy,
     $              bf_sublayer_used,
     $              nbf_interface_used,
     $              i1,j1,
     $              newgrdpt_data)

            case(2)
               call initialize_data_for_test_compute_newgrpdt_config2(
     $              dx,dy,
     $              interior_nodes0,
     $              interior_nodes1,
     $              bf_sublayer_used,
     $              nbf_interface_used,
     $              i1,j1,
     $              newgrdpt_data)

            case(3)
               call initialize_data_for_test_compute_newgrpdt_config2(
     $              dx,dy,
     $              interior_nodes0,
     $              interior_nodes1,
     $              bf_sublayer_used,
     $              nbf_interface_used,
     $              i1,j1,
     $              newgrdpt_data)

            case(4)
               call initialize_data_for_test_compute_newgrpdt_config4(
     $              dx,dy,
     $              interior_nodes0,
     $              interior_nodes1,
     $              bf_sublayer_used,
     $              nbf_interface_used,
     $              i1,j1,
     $              newgrdpt_data)

            case(5)
               call initialize_data_for_test_compute_newgrpdt_config5(
     $              dx,dy,
     $              interior_nodes0,
     $              interior_nodes1,
     $              bf_sublayer_used,
     $              nbf_interface_used,
     $              i1,j1,
     $              newgrdpt_data)

            case default
               print '(''test_nbf_interface'')'
               print '(''initialize_data_for-test_compute_newgrdpt'')'
               print '(''test not implemented: '')', config
               stop ''

          end select

        end subroutine initialize_data_for_test_compute_newgrpdt


        subroutine initialize_data_for_test_compute_newgrpdt_config1(
     $     dx,dy,
     $     bf_sublayer_used,
     $     nbf_interface_used,
     $     i1,j1,
     $     newgrdpt_data)

          implicit none

          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(bf_sublayer)  , pointer      , intent(inout) :: bf_sublayer_used
          type(nbf_interface_newgrdpt)      , intent(out)   :: nbf_interface_used
          integer(ikind)                    , intent(out)   :: i1
          integer(ikind)                    , intent(out)   :: j1
          real(rkind)        , dimension(ne), intent(out)   :: newgrdpt_data

          integer(ikind), dimension(:,:)  , allocatable :: align0
          real(rkind)   , dimension(:)    , allocatable :: x_map0
          real(rkind)   , dimension(:)    , allocatable :: y_map0
          real(rkind)   , dimension(:,:,:), allocatable :: nodes0
          integer       , dimension(:,:)  , allocatable :: grdpts_id0

          integer(ikind), dimension(2,2)                :: align1
          real(rkind)   , dimension(:)    , allocatable :: x_map1
          real(rkind)   , dimension(:)    , allocatable :: y_map1
          real(rkind)   , dimension(:,:,:), allocatable :: nodes1

          integer                                     :: i

          !attribute initialization of buffer layer at t-dt
          allocate(align0(2,2))
          align0(1,1) = 5
          align0(1,2) = 6
          align0(2,1) = ny-1
          align0(2,2) = ny

          allocate(x_map0(6))
          allocate(y_map0(6))

          do i=1,6
             x_map0(i) = (i+1)*dx
          end do

          do i=1,6
             y_map0(i) = (i+ny-3)*dy
          end do

          allocate(nodes0(6,6,ne))

          nodes0(5:6,5:6,:) = reshape((/
     $         2.0d0, -0.5d0,
     $         3.0d0, 1.25d0,
     $         
     $       -8.25d0, 7.85d0,
     $        3.26d0, 9.23d0,
     $         
     $       -0.75d0,-0.45d0,
     $        3.26d0, 6.15d0
     $         /),
     $         (/2,2,3/))

          allocate(grdpts_id0(6,6))
          
          grdpts_id0 = reshape((/
     $         1,1,1,1,1,1,
     $         1,1,1,1,1,1,
     $         2,2,1,1,2,2,
     $         3,2,1,1,2,3,
     $         3,2,2,2,2,3,
     $         3,3,3,3,3,3
     $         /),
     $         (/6,6/))


          !attribute initialization of buffer layer at t
          align1(1,1) = 5
          align1(1,2) = 6
          align1(2,1) = ny-1
          align1(2,2) = ny+1

          allocate(x_map1(6))
          allocate(y_map1(7))

          x_map1 = x_map0

          do i=1,7
             y_map1(i) = (i+ny-3)*dy
          end do

          allocate(nodes1(6,7,ne))

          nodes1(5:6,5:6,:) = reshape((/
     $         1.0d0,  0.5d0,
     $        2.45d0,-0.26d0,
     $         
     $        2.05d0, 9.26d0,
     $       -2.15d0, 7.85d0,
     $         
     $        0.25d0, 0.10d0,
     $       -0.75d0,-8.52d0
     $         /),
     $         (/2,2,3/))


          !set the nodes in bf_compute object
          call bf_sublayer_used%bf_compute_used%set_alignment(align0)
          call bf_sublayer_used%bf_compute_used%set_grdpts_id(grdpts_id0)
          call bf_sublayer_used%bf_compute_used%set_x_map(x_map0)
          call bf_sublayer_used%bf_compute_used%set_y_map(y_map0)
          call bf_sublayer_used%bf_compute_used%set_nodes(nodes0)

          !set the nodes in the bf_layer_object
          call bf_sublayer_used%ini(N)
          call bf_sublayer_used%set_alignment_tab(align1)
          call bf_sublayer_used%set_x_map(x_map1)
          call bf_sublayer_used%set_y_map(y_map1)
          call bf_sublayer_used%set_nodes(nodes1)


          !set the indices of the grid point computed
          i1 = 6
          j1 = 7

          
          !newgrdpt_data
          newgrdpt_data = [-2.77171875d0,9.87d0,4.370469d0]


          !initialize the nbf_interface with no links
          call nbf_interface_used%ini()

        end subroutine initialize_data_for_test_compute_newgrpdt_config1


        subroutine initialize_data_for_test_compute_newgrpdt_config2(
     $     dx,dy,
     $     interior_nodes0,
     $     interior_nodes1,
     $     bf_sublayer_used,
     $     nbf_interface_used,
     $     i1,j1,
     $     newgrdpt_data)

          implicit none

          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: interior_nodes1
          type(bf_sublayer), pointer        , intent(inout) :: bf_sublayer_used
          type(nbf_interface_newgrdpt)      , intent(out)   :: nbf_interface_used
          integer(ikind)                    , intent(out)   :: i1
          integer(ikind)                    , intent(out)   :: j1
          real(rkind)        , dimension(ne), intent(out)   :: newgrdpt_data

          integer(ikind), dimension(:,:)  , allocatable :: align0
          real(rkind)   , dimension(:)    , allocatable :: x_map0
          real(rkind)   , dimension(:)    , allocatable :: y_map0
          real(rkind)   , dimension(:,:,:), allocatable :: nodes0
          integer       , dimension(:,:)  , allocatable :: grdpts_id0

          integer(ikind), dimension(2,2)                :: align1
          real(rkind)   , dimension(:)    , allocatable :: x_map1
          real(rkind)   , dimension(:)    , allocatable :: y_map1
          real(rkind)   , dimension(:,:,:), allocatable :: nodes1

          integer                                     :: i

          !attribute initialization of buffer layer at t-dt
          allocate(align0(2,2))
          align0(1,1) = 5
          align0(1,2) = 6
          align0(2,1) = ny-1
          align0(2,2) = ny

          allocate(x_map0(6))
          allocate(y_map0(6))

          do i=1,6
             x_map0(i) = (i+1)*dx
          end do

          do i=1,6
             y_map0(i) = (i+ny-3)*dy
          end do

          allocate(nodes0(6,6,ne))

          nodes0(5:6,3:6,:) = reshape((/
     $         0.6d0,  0.2d0,
     $         1.0d0,  2.0d0,
     $         0.5d0, -0.5d0,
     $         2.3d0,-6.23d0,
     $         
     $       -3.25d0, 6.12d0,
     $        0.25d0,-0.75d0,
     $         0.1d0,-0.45d0,
     $       -2.56d0, 1.23d0,
     $         
     $        9.26d0, 3.25d0,
     $        2.05d0,-8.25d0,
     $        9.26d0, 7.85d0,
     $       -6.23d0, 2.36d0
     $         /),
     $         (/2,4,3/))

          allocate(grdpts_id0(6,6))
          
          grdpts_id0 = reshape((/
     $         1,1,1,1,1,1,
     $         1,1,1,1,1,1,
     $         2,2,1,1,2,2,
     $         3,2,1,1,2,3,
     $         3,2,2,2,2,3,
     $         3,3,3,3,3,3
     $         /),
     $         (/6,6/))


          !attribute initialization of buffer layer at t
          align1(1,1) = 5
          align1(1,2) = 7
          align1(2,1) = ny-1
          align1(2,2) = ny

          allocate(x_map1(7))
          allocate(y_map1(6))

          x_map1(1:6) = x_map0(:)
          x_map1(7)   = x_map1(6)+dx

          y_map1      = y_map0

          allocate(nodes1(7,6,ne))

          nodes1(5:7,3:6,:) = reshape((/
     $         2.3d0,  7.8d0, 1.2d0,
     $         8.9d0,  1.0d0,2.45d0,
     $         0.2d0,  0.5d0, 0.0d0,
     $        6.23d0,-5.15d0, 0.0d0,
     $         
     $        -5.2d0, 1.23d0, 7.15d0,
     $        9.26d0, 0.25d0,-0.75d0,
     $         2.3d0,  0.1d0,  0.0d0,
     $       7.154d0,  1.2d0,  0.0d0,
     $         
     $        9.63d0,  1.2d0, 7.32d0,
     $        1.25d0, 2.05d0,-2.15d0,
     $        7.26d0, 9.26d0,  0.0d0,
     $        5.26d0, 1.45d0,  0.0d0
     $         /),
     $         (/3,4,3/))


          !set the nodes in bf_compute object
          call bf_sublayer_used%bf_compute_used%set_alignment(align0)
          call bf_sublayer_used%bf_compute_used%set_grdpts_id(grdpts_id0)
          call bf_sublayer_used%bf_compute_used%set_x_map(x_map0)
          call bf_sublayer_used%bf_compute_used%set_y_map(y_map0)
          call bf_sublayer_used%bf_compute_used%set_nodes(nodes0)

          !set the nodes in the bf_layer_object
          call bf_sublayer_used%ini(N)
          call bf_sublayer_used%set_alignment_tab(align1)
          call bf_sublayer_used%set_x_map(x_map1)
          call bf_sublayer_used%set_y_map(y_map1)
          call bf_sublayer_used%set_nodes(nodes1)


          !set the indices of the grid point computed
          i1 = 7
          j1 = 5

          
          !newgrdpt_data
          newgrdpt_data = [-4.992695567d0,-3.40475565d0,12.39524435d0]


          !initialize the nbf_interface with no links
          call nbf_interface_used%ini()


          !initialize the interior nodes
          interior_nodes0(7:10,9:10,:) = reshape((/
     $         0.6d0,  0.2d0,  0.8d0, -2.3d0,
     $         1.0d0,  2.0d0,  3.0d0, 1.23d0,
     $         
     $       -3.25d0, 6.12d0, 7.15d0, 1.23d0,
     $        0.25d0,-0.75d0, 3.26d0, 5.86d0,
     $         
     $        9.26d0, 3.25d0, 4.15d0, 9.56d0,
     $        2.05d0,-8.25d0, 3.26d0, 1.23d0
     $         /),
     $         (/4,2,3/))

          interior_nodes1(7:10,9:10,:) = reshape((/
     $         2.3d0,  7.8d0,  1.2d0,  1.5d0,
     $         8.9d0,  1.0d0, 2.45d0,  3.0d0,
     $         
     $        -5.2d0, 1.23d0, 7.15d0,  6.2d0,
     $        9.26d0, 0.25d0,-0.75d0, 3.26d0,
     $         
     $        9.63d0,  1.2d0, 7.32d0, 1.52d0,
     $        1.25d0, 2.05d0,-2.15d0, 3.26d0
     $         /),
     $         (/4,2,3/))

        end subroutine initialize_data_for_test_compute_newgrpdt_config2


        subroutine initialize_data_for_test_compute_newgrpdt_config3(
     $     dx,dy,
     $     interior_nodes0,
     $     interior_nodes1,
     $     bf_sublayer_used,
     $     nbf_interface_used,
     $     i1,j1,
     $     newgrdpt_data)

          implicit none

          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: interior_nodes1
          type(bf_sublayer), pointer        , intent(inout) :: bf_sublayer_used
          type(nbf_interface_newgrdpt)      , intent(out)   :: nbf_interface_used
          integer(ikind)                    , intent(out)   :: i1
          integer(ikind)                    , intent(out)   :: j1
          real(rkind)        , dimension(ne), intent(out)   :: newgrdpt_data

          integer(ikind), dimension(:,:)  , allocatable :: align0
          real(rkind)   , dimension(:)    , allocatable :: x_map0
          real(rkind)   , dimension(:)    , allocatable :: y_map0
          real(rkind)   , dimension(:,:,:), allocatable :: nodes0
          integer       , dimension(:,:)  , allocatable :: grdpts_id0

          integer(ikind), dimension(2,2)                :: align1
          real(rkind)   , dimension(:)    , allocatable :: x_map1
          real(rkind)   , dimension(:)    , allocatable :: y_map1
          real(rkind)   , dimension(:,:,:), allocatable :: nodes1

          type(bf_sublayer), pointer :: bf_sublayer_n_used

          integer                                     :: i

          !attribute initialization of buffer layer at t-dt
          allocate(align0(2,2))
          align0(1,1) = 7
          align0(1,2) = 8
          align0(2,1) = ny-1
          align0(2,2) = ny

          allocate(x_map0(6))
          allocate(y_map0(6))

          do i=1,6
             x_map0(i) = (i+3)*dx
          end do

          do i=1,6
             y_map0(i) = (i+ny-3)*dy
          end do

          allocate(nodes0(6,6,ne))

          nodes0(5:6,3:6,:) = reshape((/
     $         0.6d0,  0.2d0,
     $         1.0d0,  2.0d0,
     $         0.5d0, -0.5d0,
     $         2.3d0,-6.23d0,
     $         
     $       -3.25d0, 6.12d0,
     $        0.25d0,-0.75d0,
     $         0.1d0,-0.45d0,
     $       -2.56d0, 1.23d0,
     $         
     $        9.26d0, 3.25d0,
     $        2.05d0,-8.25d0,
     $        9.26d0, 7.85d0,
     $       -6.23d0, 2.36d0
     $         /),
     $         (/2,4,3/))

          allocate(grdpts_id0(6,6))
          
          grdpts_id0 = reshape((/
     $         1,1,1,1,1,1,
     $         1,1,1,1,1,1,
     $         2,2,1,1,2,2,
     $         3,2,1,1,2,3,
     $         3,2,2,2,2,3,
     $         3,3,3,3,3,3
     $         /),
     $         (/6,6/))


          !attribute initialization of buffer layer at t
          align1(1,1) = 7
          align1(1,2) = 9
          align1(2,1) = ny-1
          align1(2,2) = ny

          allocate(x_map1(7))
          allocate(y_map1(6))

          x_map1(1:6) = x_map0(:)
          x_map1(7)   = x_map1(6)+dx

          y_map1      = y_map0

          allocate(nodes1(7,6,ne))

          nodes1(5:7,3:6,:) = reshape((/
     $         2.3d0,  7.8d0, 1.2d0,
     $         8.9d0,  1.0d0,2.45d0,
     $         0.2d0,  0.5d0, 0.0d0,
     $        6.23d0,-5.15d0, 0.0d0,
     $         
     $        -5.2d0, 1.23d0, 7.15d0,
     $        9.26d0, 0.25d0,-0.75d0,
     $         2.3d0,  0.1d0,  0.0d0,
     $       7.154d0,  1.2d0,  0.0d0,
     $         
     $        9.63d0,  1.2d0, 7.32d0,
     $        1.25d0, 2.05d0,-2.15d0,
     $        7.26d0, 9.26d0,  0.0d0,
     $        5.26d0, 1.45d0,  0.0d0
     $         /),
     $         (/3,4,3/))


          !set the nodes in bf_compute object
          call bf_sublayer_used%bf_compute_used%set_alignment(align0)
          call bf_sublayer_used%bf_compute_used%set_grdpts_id(grdpts_id0)
          call bf_sublayer_used%bf_compute_used%set_x_map(x_map0)
          call bf_sublayer_used%bf_compute_used%set_y_map(y_map0)
          call bf_sublayer_used%bf_compute_used%set_nodes(nodes0)

          !set the nodes in the bf_layer_object
          call bf_sublayer_used%ini(N)
          call bf_sublayer_used%set_alignment_tab(align1)
          call bf_sublayer_used%set_x_map(x_map1)
          call bf_sublayer_used%set_y_map(y_map1)
          call bf_sublayer_used%set_nodes(nodes1)


          !set the indices of the grid point computed
          i1 = 7
          j1 = 5

          
          !newgrdpt_data
          newgrdpt_data = [-4.992695567d0,-3.40475565d0,12.39524435d0]


          !initialize the nbf_interface with no links
          call nbf_interface_used%ini()

          
          !initialize the neighboring buffer layer
          allocate(align0(2,2))
          align0(1,1) = 9
          align0(1,2) = 11
          align0(2,1) = 7
          align0(2,2) = 8
          
          allocate(x_map0(7))
          allocate(y_map0(6))

          do i=1,7
             x_map0(i) = (i+5)*dx
          end do

          do i=1,6
             y_map0(i) = (i+3)*dy
          end do
          
          allocate(nodes0(7,6,ne))

          nodes0(2:5,5:6,:) = reshape((/
     $         0.6d0,  0.2d0,  0.8d0, -2.3d0,
     $         1.0d0,  2.0d0,  3.0d0, 1.23d0,
     $         
     $       -3.25d0, 6.12d0, 7.15d0, 1.23d0,
     $        0.25d0,-0.75d0, 3.26d0, 5.86d0,
     $         
     $        9.26d0, 3.25d0, 4.15d0, 9.56d0,
     $        2.05d0,-8.25d0, 3.26d0, 1.23d0
     $         /),
     $         (/4,2,3/))

          allocate(grdpts_id0(7,6))

          grdpts_id0 = reshape((/
     $         1,1,2,3,3,3,3,
     $         1,1,2,2,2,2,3,
     $         1,1,1,1,1,2,3,
     $         1,1,1,1,1,2,3,
     $         1,1,2,2,2,2,3,
     $         1,1,2,3,3,3,3
     $         /),
     $         (/7,6/))

          !attribute initialization of buffer layer at t
          align1(1,1) = 9
          align1(1,2) = 11
          align1(2,1) = 7
          align1(2,2) = 8

          allocate(x_map1(7))
          allocate(y_map1(6))

          x_map1 = x_map0
          y_map1 = y_map0

          allocate(nodes1(7,6,ne))

           nodes1(2:5,5:6,:) = reshape((/
     $         2.3d0,  7.8d0,  1.2d0,  1.5d0,
     $         8.9d0,  1.0d0, 2.45d0,  3.0d0,
     $         
     $        -5.2d0, 1.23d0, 7.15d0,  6.2d0,
     $        9.26d0, 0.25d0,-0.75d0, 3.26d0,
     $         
     $        9.63d0,  1.2d0, 7.32d0, 1.52d0,
     $        1.25d0, 2.05d0,-2.15d0, 3.26d0
     $         /),
     $         (/4,2,3/))

          !allocate the neighboring bufferlayer
          allocate(bf_sublayer_n_used)
          
          !set the nodes in bf_compute object
          call bf_sublayer_n_used%bf_compute_used%set_alignment(align0)
          call bf_sublayer_n_used%bf_compute_used%set_grdpts_id(grdpts_id0)
          call bf_sublayer_n_used%bf_compute_used%set_x_map(x_map0)
          call bf_sublayer_n_used%bf_compute_used%set_y_map(y_map0)
          call bf_sublayer_n_used%bf_compute_used%set_nodes(nodes0)

          !set the nodes in the bf_layer_object
          call bf_sublayer_n_used%ini(E)
          call bf_sublayer_n_used%set_alignment_tab(align1)
          call bf_sublayer_n_used%set_x_map(x_map1)
          call bf_sublayer_n_used%set_y_map(y_map1)
          call bf_sublayer_n_used%set_nodes(nodes1)

          !set as a potential neighbor in nbf_interface
          call nbf_interface_used%link_neighbor1_to_bf_sublayer(bf_sublayer_used)
          call nbf_interface_used%link_neighbor2_to_bf_sublayer(bf_sublayer_n_used)

          !initialize the interior nodes
          interior_nodes0(8:10,9:10,:) = reshape((/
     $         0.6d0,  0.2d0,  0.8d0,
     $         1.0d0,  2.0d0,  3.0d0,
     $         
     $       -3.25d0, 6.12d0, 7.15d0,
     $        0.25d0,-0.75d0, 3.26d0,
     $         
     $        9.26d0, 3.25d0, 4.15d0,
     $        2.05d0,-8.25d0, 3.26d0
     $         /),
     $         (/3,2,3/))

          interior_nodes1(8:10,9:10,:) = reshape((/
     $         2.3d0,  7.8d0,  1.2d0,
     $         8.9d0,  1.0d0, 2.45d0,
     $         
     $        -5.2d0, 1.23d0, 7.15d0,
     $        9.26d0, 0.25d0,-0.75d0,
     $         
     $        9.63d0,  1.2d0, 7.32d0,
     $        1.25d0, 2.05d0,-2.15d0
     $         /),
     $         (/3,2,3/))

        end subroutine initialize_data_for_test_compute_newgrpdt_config3


        subroutine initialize_data_for_test_compute_newgrpdt_config4(
     $     dx,dy,
     $     interior_nodes0,
     $     interior_nodes1,
     $     bf_sublayer_used,
     $     nbf_interface_used,
     $     i1,j1,
     $     newgrdpt_data)

          implicit none

          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: interior_nodes1
          type(bf_sublayer)  , pointer      , intent(inout) :: bf_sublayer_used
          type(nbf_interface_newgrdpt)      , intent(out)   :: nbf_interface_used
          integer(ikind)                    , intent(out)   :: i1
          integer(ikind)                    , intent(out)   :: j1
          real(rkind)        , dimension(ne), intent(out)   :: newgrdpt_data

          integer(ikind), dimension(2,2)                :: align1
          real(rkind)   , dimension(:)    , allocatable :: x_map1
          real(rkind)   , dimension(:)    , allocatable :: y_map1
          real(rkind)   , dimension(:,:,:), allocatable :: nodes1

          integer                                     :: i

          !attribute initialization of buffer layer at t
          align1(1,1) = 7
          align1(1,2) = 8
          align1(2,1) = ny-1
          align1(2,2) = ny-1

          allocate(x_map1(6))
          allocate(y_map1(5))

          do i=1,6
             x_map1(i) = (3+i)*dx
          end do

          do i=1,5
             y_map1(i) = (5+i)*dy
          end do

          allocate(nodes1(6,5,ne))

          nodes1(5:6,3:4,:) = reshape((/
     $         1.0d0,  0.5d0,
     $        2.45d0,-0.26d0,
     $         
     $        2.05d0, 9.26d0,
     $       -2.15d0, 7.85d0,
     $         
     $        0.25d0, 0.10d0,
     $       -0.75d0,-8.52d0
     $         /),
     $         (/2,2,3/))


          !set the nodes in the bf_layer_object
          call bf_sublayer_used%ini(N)
          call bf_sublayer_used%set_alignment_tab(align1)
          call bf_sublayer_used%set_x_map(x_map1)
          call bf_sublayer_used%set_y_map(y_map1)
          call bf_sublayer_used%set_nodes(nodes1)


          !set the indices of the grid point computed
          i1 = 6
          j1 = 5

          
          !newgrdpt_data
          newgrdpt_data = [-2.77171875d0,9.87d0,4.370469d0]


          !initialize the nbf_interface with no links
          call nbf_interface_used%ini()

          !initialize the interior_nodes
          interior_nodes0(9:10,9:10,:) = reshape((/
     $         2.0d0, -0.5d0,
     $         3.0d0, 1.25d0,
     $         
     $       -8.25d0, 7.85d0,
     $        3.26d0, 9.23d0,
     $         
     $       -0.75d0,-0.45d0,
     $        3.26d0, 6.15d0
     $         /),
     $         (/2,2,3/))

          interior_nodes1(9:10,9:10,:) = reshape((/
     $         1.0d0,  0.5d0,
     $        2.45d0,-0.26d0,
     $         
     $        2.05d0, 9.26d0,
     $       -2.15d0, 7.85d0,
     $         
     $        0.25d0, 0.10d0,
     $       -0.75d0,-8.52d0
     $         /),
     $         (/2,2,3/))

        end subroutine initialize_data_for_test_compute_newgrpdt_config4


        subroutine initialize_data_for_test_compute_newgrpdt_config5(
     $     dx,dy,
     $     interior_nodes0,
     $     interior_nodes1,
     $     bf_sublayer_used,
     $     nbf_interface_used,
     $     i1,j1,
     $     newgrdpt_data)

          implicit none

          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: interior_nodes1
          type(bf_sublayer)  , pointer      , intent(inout) :: bf_sublayer_used
          type(nbf_interface_newgrdpt)      , intent(out)   :: nbf_interface_used
          integer(ikind)                    , intent(out)   :: i1
          integer(ikind)                    , intent(out)   :: j1
          real(rkind)        , dimension(ne), intent(out)   :: newgrdpt_data

          integer(ikind), dimension(2,2)                :: align1
          real(rkind)   , dimension(:)    , allocatable :: x_map1
          real(rkind)   , dimension(:)    , allocatable :: y_map1
          real(rkind)   , dimension(:,:,:), allocatable :: nodes1
          integer(ikind), dimension(:,:)  , allocatable :: grdpts_id1

          integer :: i


          !attribute initialization of buffer layer at t
          align1(1,1) = 8
          align1(1,2) = 8
          align1(2,1) = 9
          align1(2,2) = 9

          allocate(x_map1(5))
          allocate(y_map1(5))

          do i=1,5
             x_map1(i) = (4+i)*dx
          end do

          do i=1,5
             y_map1(i) = (5+i)*dy
          end do

          allocate(grdpts_id1(5,5))

          grdpts_id1 = reshape((/
     $         1,1,1,2,3,
     $         1,1,1,2,3,
     $         2,2,2,2,3,
     $         3,3,3,3,3,
     $         0,0,0,0,0
     $         /),
     $         (/5,5/))         

          allocate(nodes1(5,5,ne))

          nodes1(4:5,3:4,:) = reshape((/
     $         1.0d0,  0.5d0,
     $        2.45d0,-0.26d0,
     $         
     $        2.05d0, 9.26d0,
     $       -2.15d0, 7.85d0,
     $         
     $        0.25d0, 0.10d0,
     $       -0.75d0,-8.52d0
     $         /),
     $         (/2,2,3/))


          !set the nodes in the bf_layer_object
          call bf_sublayer_used%ini(N)
          call bf_sublayer_used%set_alignment_tab(align1)
          call bf_sublayer_used%set_x_map(x_map1)
          call bf_sublayer_used%set_y_map(y_map1)
          call bf_sublayer_used%set_nodes(nodes1)
          call bf_sublayer_used%set_grdpts_id(grdpts_id1)


          !set the indices of the grid point computed
          i1 = 5
          j1 = 5

          
          !newgrdpt_data
          newgrdpt_data = [-2.77171875d0,9.87d0,4.370469d0]


          !initialize the nbf_interface with no links
          call nbf_interface_used%ini()

          !initialize the interior_nodes
          interior_nodes0(9:10,9:10,:) = reshape((/
     $         2.0d0, -0.5d0,
     $         3.0d0, 1.25d0,
     $         
     $       -8.25d0, 7.85d0,
     $        3.26d0, 9.23d0,
     $         
     $       -0.75d0,-0.45d0,
     $        3.26d0, 6.15d0
     $         /),
     $         (/2,2,3/))

          interior_nodes1(9:10,9:10,:) = reshape((/
     $         1.0d0,  0.5d0,
     $        2.45d0,-0.26d0,
     $         
     $        2.05d0, 9.26d0,
     $       -2.15d0, 7.85d0,
     $         
     $        0.25d0, 0.10d0,
     $       -0.75d0,-8.52d0
     $         /),
     $         (/2,2,3/))

        end subroutine initialize_data_for_test_compute_newgrpdt_config5


        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical    , intent(in) :: detailled
          logical                 :: test_validated

          if(detailled) then
             print *, int(var*1e5)
             print *, int(cst*1e5)
          end if
          
          test_validated=abs(
     $         int(var*1e5)-
     $         int(cst*1e5)).le.1
          
        end function is_test_validated


        function test_update_bf_grdpts_after_increase(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(nbf_interface_newgrdpt)                :: nbf_interface_used
          type(bf_sublayer), pointer                  :: bf_sublayer_used
                                                      
          type(pmodel_eq)                             :: p_model
          real(rkind)                                 :: t
          real(rkind)                                 :: dt
          real(rkind)                                 :: dx
          real(rkind)                                 :: dy
          real(rkind)   , dimension(nx)               :: interior_x_map
          real(rkind)   , dimension(ny)               :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne)         :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne)         :: interior_nodes1
          integer(ikind)                              :: i1
          integer(ikind)                              :: j1
          real(rkind)   , dimension(ne)               :: newgrdpt_data
          real(rkind)   , dimension(ne)               :: newgrdpt
          integer(ikind), dimension(2,1)              :: selected_grdpts
          integer       , dimension(5,5)              :: grdpts_id1_data
          integer       , dimension(:,:), allocatable :: bf_grdpts_id1

          logical :: test_loc
          logical :: test_node
          logical :: test_grdpts_id

          integer(ikind) :: i,j
          integer        :: k
          integer        :: config

          !test requirements
          if((nx.ne.10).or.(ny.ne.10).or.(ne.ne.3).or.
     $         (.not.(is_test_validated(c,0.5d0,.false.)))) then

             print '(''test_nbf_interface'')'
             print '(''test_compute_newgrdpt'')'
             print '(''the test requires:'')'
             print '(''nx=10 : '',L1)', nx.eq.10
             print '(''ny=10 : '',L1)', ny.eq.10
             print '(''ne= 3 : '',L1)', ne.eq.3
             print '('' -> the wave2d model should be used'')'
             print '(''c=0.5 : '',L1)', is_test_validated(c,0.5d0,.false.)
             print '()'

             stop ''

          end if

          config = 5
          test_validated = .true.

          allocate(bf_sublayer_used)


          !initialize the interior x_map and y_map
          dx = 0.25d0
          dy = 0.25d0

          do k=1, nx
             interior_x_map(k) = (k-1)*dx
          end do

          do k=1, ny
             interior_y_map(k) = (k-1)*dy
          end do


          !initialize the buffer layer(s), the interior
          !nodes and the links in the nbf_interface for
          !the test
          call initialize_data_for_test_compute_newgrpdt(
     $         config,
     $         dx,dy,
     $         interior_nodes0,
     $         interior_nodes1,
     $         bf_sublayer_used,
     $         nbf_interface_used,
     $         i1,j1,
     $         newgrdpt_data)
          
          
          !test the computation of the new grdpt and the
          !update the neighboring grid points
          t  = 0.0d0
          dt = 0.25d0

          selected_grdpts(:,1) = [8,9]

          call nbf_interface_used%update_bf_grdpts_after_increase(
     $         bf_sublayer_used,
     $         p_model, t,dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         selected_grdpts)


          !compare the results with what is expected for the nodes
          newgrdpt = bf_sublayer_used%get_nodes([i1,j1])

          do k=1, ne
             test_loc = is_test_validated(
     $            newgrdpt(k), newgrdpt_data(k), detailled)
             test_validated = test_validated.and.test_loc
          end do
          test_node = test_validated
          print '(''test_node: '',L1)', test_validated


          !compare the results with what is expected for the grdpts_id
          grdpts_id1_data = reshape((/
     $         1,1,1,2,3,
     $         1,1,1,2,3,
     $         2,2,1,2,3,
     $         3,2,2,2,3,
     $         3,3,3,3,3
     $         /),
     $         (/5,5/))

          call bf_sublayer_used%get_grdpts_id(bf_grdpts_id1)
          
          test_validated=.true.
          do j=1,5
             do i=1,5
                test_loc = grdpts_id1_data(i,j).eq.bf_grdpts_id1(i,j)
                test_validated = test_validated.and.test_loc
                if(.not.test_loc) then
                   print '(2I2,''->'',2I2)',
     $                  i,j,
     $                  grdpts_id1_data(i,j),
     $                  bf_grdpts_id1(i,j)
                end if
             end do
          end do
          test_grdpts_id = test_validated
          print '(''test_grdpts_id: '',L1)', test_grdpts_id


          test_validated = test_node.and.test_grdpts_id

        end function test_update_bf_grdpts_after_increase


        function test_finalize_grdpts_around_new_interior_pt(
     $     config,
     $     detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: config
          logical, intent(in) :: detailled
          logical             :: test_validated

          type(nbf_interface_newgrdpt)                :: nbf_interface_used
          type(bf_sublayer)                           :: bf_sublayer_used
          integer(ikind), dimension(2)                :: match_table
          integer(ikind)                              :: i_center
          integer(ikind)                              :: j_center
          integer       , dimension(:,:), allocatable :: grdpts_id
          integer       , dimension(:,:), allocatable :: test_grdpts_id


          !============================================================
          !initialize the inputs
          !============================================================
          call get_data_for_test_finalize_grdpts_around(
     $         config,
     $         nbf_interface_used,
     $         bf_sublayer_used,
     $         match_table,
     $         i_center,
     $         j_center,
     $         test_grdpts_id)          


          !============================================================
          !test finalize_grdpts_around_new_interior_pt
          !============================================================
          call finalize_grdpts_around_new_interior_pt(
     $         nbf_interface_used,
     $         bf_sublayer_used,
     $         [i_center,j_center],
     $         match_table)


          !============================================================
          !compare results
          !============================================================
          call bf_sublayer_used%get_grdpts_id(grdpts_id)
          test_validated = compare_grdpts_id(
     $         grdpts_id,
     $         test_grdpts_id,
     $         detailled)

        end function test_finalize_grdpts_around_new_interior_pt


        subroutine get_data_for_test_finalize_grdpts_around(
     $     config,
     $     nbf_interface_used,
     $     bf_sublayer_used,
     $     match_table,
     $     i_center,
     $     j_center,
     $     test_grdpts_id)
 
          implicit none

          integer                             , intent(in)    :: config
          type(nbf_interface_newgrdpt)        , intent(inout) :: nbf_interface_used
          type(bf_sublayer)                   , intent(inout) :: bf_sublayer_used
          integer(ikind), dimension(2)        , intent(out)   :: match_table
          integer(ikind)                      , intent(out)   :: i_center
          integer(ikind)                      , intent(out)   :: j_center
          integer, dimension(:,:), allocatable, intent(out)   :: test_grdpts_id
          

          integer(ikind), dimension(2,2)              :: bf_align
          integer       , dimension(:,:), allocatable :: grdpts_id


          select case(config)

            case(1)

               !initialize the nbf_interface
               call nbf_interface_used%ini()
          
               !initialize the buffer layer
               call bf_sublayer_used%ini(S)

               allocate(grdpts_id(7,5))
               
               grdpts_id = reshape((/
     $              3,3,3,3,3,3,3,
     $              2,2,2,2,2,2,2,
     $              1,1,2,1,1,1,1,
     $              1,1,1,1,1,1,1,
     $              1,1,1,1,1,1,1/),
     $              (/7,5/))
               
               bf_align = reshape((/
     $              6,2,
     $              8,2/),
     $              (/2,2/))
               
               call bf_sublayer_used%set_alignment_tab(bf_align)
               call bf_sublayer_used%set_grdpts_id(grdpts_id)
               
               match_table = bf_sublayer_used%get_general_to_local_coord_tab()
               i_center    = 4
               j_center    = 3
               
               allocate(test_grdpts_id(7,5))
               test_grdpts_id = reshape((/
     $              3,3,3,3,3,3,3,
     $              2,2,2,2,2,2,2,
     $              1,1,1,1,1,1,1,
     $              1,1,1,1,1,1,1,
     $              1,1,1,1,1,1,1/),
     $              (/7,5/))

            case(2)

               !initialize the nbf_interface
               call nbf_interface_used%ini()
          
               !initialize the buffer layer
               call bf_sublayer_used%ini(S)

               allocate(grdpts_id(9,5))
               
               grdpts_id = reshape((/
     $              3,3,3,3,3,3,3,3,3,
     $              2,2,2,2,2,2,2,2,2,
     $              1,1,2,2,1,2,2,1,1,
     $              1,1,1,1,1,1,1,1,1,
     $              1,1,1,1,1,1,1,1,1/),
     $              (/9,5/))
               
               bf_align = reshape((/
     $              3,2,
     $              7,2/),
     $              (/2,2/))
               
               call bf_sublayer_used%set_alignment_tab(bf_align)
               call bf_sublayer_used%set_grdpts_id(grdpts_id)
               
               match_table = bf_sublayer_used%get_general_to_local_coord_tab()
               i_center    = 5
               j_center    = 3
               
               allocate(test_grdpts_id(9,5))
               test_grdpts_id = reshape((/
     $              3,3,3,3,3,3,3,3,3,
     $              2,2,2,2,2,2,2,2,2,
     $              1,1,1,1,1,1,1,1,1,
     $              1,1,1,1,1,1,1,1,1,
     $              1,1,1,1,1,1,1,1,1/),
     $              (/9,5/))

            case default
               print '(''test_nbf_interface_newgrdpt'')'
               print '(''get_data_for_test_finalize_grdpts'')'
               print '(''config for test not recognized: '',I2)', config
               stop ''               

          end select               

        end subroutine get_data_for_test_finalize_grdpts_around


        function compare_grdpts_id(grdpts_id,test_grdpts_id,detailled)
     $     result(test_validated)

          implicit none

          integer, dimension(:,:), intent(in) :: grdpts_id
          integer, dimension(:,:), intent(in) :: test_grdpts_id
          logical                , intent(in) :: detailled
          logical                             :: test_validated


          logical :: test_loc
          integer :: i,j

          
          test_validated = .true.


          do j=1, size(grdpts_id,2)
             do i=1, size(grdpts_id,1)

                test_loc = grdpts_id(i,j).eq.test_grdpts_id(i,j)
                test_validated = test_loc.and.test_validated
                if(detailled.and.(.not.test_loc)) then

                   print '(''** test failed at '',2I2,''**'')', i,j
                   print '(I2,'' -> '',I2)', grdpts_id(i,j), test_grdpts_id(i,j)

                end if

             end do
          end do

        end function compare_grdpts_id


        function test_finalize_grdpts_for_bc_pt_crenel(
     $     config,
     $     detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: config
          logical, intent(in) :: detailled
          logical             :: test_validated

          type(nbf_interface_newgrdpt)                :: nbf_interface_used
          type(bf_sublayer), pointer                  :: bf_sublayer_used
          integer(ikind), dimension(2)                :: match_table
          integer(ikind)                              :: i_center
          integer(ikind)                              :: j_center
          integer       , dimension(:,:), allocatable :: grdpts_id
          integer       , dimension(:,:), allocatable :: test_grdpts_id


          allocate(bf_sublayer_used)

          !============================================================
          !initialize the inputs
          !============================================================
          call get_data_for_test_finalize_grdpts_crenel(
     $         config,
     $         nbf_interface_used,
     $         bf_sublayer_used,
     $         match_table,
     $         i_center,
     $         j_center,
     $         test_grdpts_id)          


          !============================================================
          !test finalize_grdpts_around_new_interior_pt
          !============================================================
          call finalize_grdpts_for_bc_pt_crenel(
     $         nbf_interface_used,
     $         bf_sublayer_used,
     $         [i_center,j_center],
     $         match_table)


          !============================================================
          !compare results
          !============================================================
          call bf_sublayer_used%get_grdpts_id(grdpts_id)
          test_validated = compare_grdpts_id(
     $         grdpts_id,
     $         test_grdpts_id,
     $         detailled)

        end function test_finalize_grdpts_for_bc_pt_crenel


        subroutine get_data_for_test_finalize_grdpts_crenel(
     $     config,
     $     nbf_interface_used,
     $     bf_sublayer_used,
     $     match_table,
     $     i_center,
     $     j_center,
     $     test_grdpts_id)
 
          implicit none

          integer                             , intent(in)    :: config
          type(nbf_interface_newgrdpt)        , intent(inout) :: nbf_interface_used
          type(bf_sublayer), pointer          , intent(inout) :: bf_sublayer_used
          integer(ikind), dimension(2)        , intent(out)   :: match_table
          integer(ikind)                      , intent(out)   :: i_center
          integer(ikind)                      , intent(out)   :: j_center
          integer, dimension(:,:), allocatable, intent(out)   :: test_grdpts_id
          

          integer(ikind), dimension(2,2)              :: bf_align
          integer       , dimension(:,:), allocatable :: grdpts_id

          type(bf_sublayer), pointer :: bf_sublayer_used2

          select case(config)

            case(1)

               !initialize the nbf_interface
               call nbf_interface_used%ini()
          
               !initialize the buffer layer
               call bf_sublayer_used%ini(E)

               allocate(grdpts_id(5,10))
               
               grdpts_id = reshape((/
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,2,2,3,
     $              1,1,2,3,3,
     $              1,1,2,3,3,
     $              1,1,2,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3
     $              /),
     $              (/5,10/))
               
               bf_align = reshape((/
     $              9,3,
     $              9,8/),
     $              (/2,2/))
               
               call bf_sublayer_used%set_alignment_tab(bf_align)
               call bf_sublayer_used%set_grdpts_id(grdpts_id)
               
               match_table = bf_sublayer_used%get_general_to_local_coord_tab()
               i_center    = 3
               j_center    = 4
               
               allocate(test_grdpts_id(5,10))
               test_grdpts_id = reshape((/
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3/),
     $              (/5,10/))

            case(2)

               !initialize the nbf_interface
               call nbf_interface_used%ini()
          
               !initialize the buffer layer
               call bf_sublayer_used%ini(E)

               allocate(grdpts_id(5,10))
               
               grdpts_id = reshape((/
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,2,2,3,
     $              1,1,2,3,3
     $              /),
     $              (/5,10/))
               
               bf_align = reshape((/
     $              9,3,
     $              9,8/),
     $              (/2,2/))
               
               call bf_sublayer_used%set_alignment_tab(bf_align)
               call bf_sublayer_used%set_grdpts_id(grdpts_id)
               call bf_sublayer_used%set_neighbor2_share(neighbor2_share=.true.)

               match_table = bf_sublayer_used%get_general_to_local_coord_tab()
               i_center    = 3
               j_center    = 8

               allocate(test_grdpts_id(5,10))
               test_grdpts_id = reshape((/
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3/),
     $              (/5,10/))

               !add a second buffer layer
               allocate(bf_sublayer_used2)
               call bf_sublayer_used2%ini(N)

               allocate(grdpts_id(10,9))
               
               grdpts_id = reshape((/
     $              1,1,1,1,1,1,1,1,2,3,
     $              1,1,1,1,1,1,1,1,2,3,
     $              2,2,1,1,1,1,1,2,2,3,
     $              3,2,1,1,1,1,1,2,3,3,
     $              3,2,1,1,1,1,1,2,3,3,
     $              3,2,1,1,1,1,1,2,2,3,
     $              3,2,1,1,1,1,1,1,2,3,
     $              3,2,2,2,2,2,2,2,2,3,
     $              3,3,3,3,3,3,3,3,3,3/),
     $              (/10,9/))
               
               bf_align = reshape((/
     $              4,9,
     $              9,13/),
     $              (/2,2/))
               
               call bf_sublayer_used2%set_alignment_tab(bf_align)
               call bf_sublayer_used2%set_grdpts_id(grdpts_id)
               call bf_sublayer_used2%set_neighbor2_share(neighbor2_share=.true.)

               call nbf_interface_used%link_neighbor2_to_bf_sublayer(
     $              bf_sublayer_used)

               call nbf_interface_used%link_neighbor2_to_bf_sublayer(
     $              bf_sublayer_used2)

            case default
               print '(''test_nbf_interface_newgrdpt'')'
               print '(''get_data_for_test_finalize_grdpts'')'
               print '(''config for test not recognized: '',I2)', config
               stop ''               

          end select               

        end subroutine get_data_for_test_finalize_grdpts_crenel

      end program test_nbf_interface

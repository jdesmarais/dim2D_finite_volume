      program test_nbf_interface

        use bf_sublayer_class, only :
     $     bf_sublayer

        use nbf_interface_class, only :
     $     nbf_interface

        use parameters_constant, only :
     $       N,E

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
        
        
        !all the data needed for the computation of the new grid point
        !are located in the buffer layer
        test_validated = test_compute_newgrdpt(1, detailled)
        print '(''test_compute_newgrdpt - config 1: '',L1)', test_validated
        
        !the data needed for the computation of the new grid point need to
        !be extracted from the interior
        test_validated = test_compute_newgrdpt(2, detailled)
        print '(''test_compute_newgrdpt - config 2: '',L1)', test_validated

        !the data needed for the computation of the new grid point need to
        !be extracted from the interior and the neighboring buffer layers
        test_validated = test_compute_newgrdpt(3, detailled)
        print '(''test_compute_newgrdpt - config 3: '',L1)', test_validated
        
        !the buffer layer has just been created and no data are stored
        !at t=t-dt
        test_validated = test_compute_newgrdpt(4, detailled)
        print '(''test_compute_newgrdpt - config 4: '',L1)', test_validated
        

        contains


        function test_compute_newgrdpt(config, detailled)
     $       result(test_validated)

          implicit none
          
          integer, intent(in) :: config
          logical, intent(in) :: detailled
          logical             :: test_validated

          type(nbf_interface)              :: nbf_interface_used
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
          type(nbf_interface)               , intent(out)   :: nbf_interface_used
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
          type(nbf_interface)               , intent(out)   :: nbf_interface_used
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
          type(nbf_interface)               , intent(out)   :: nbf_interface_used
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
          type(nbf_interface)               , intent(out)   :: nbf_interface_used
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
          type(nbf_interface)               , intent(out)   :: nbf_interface_used
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

      end program test_nbf_interface

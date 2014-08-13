      program test_lodi_poinsot_ns2d
      
        use sd_operators_fd_module, only :
     $       gradient_x_interior,
     $       gradient_y_interior

        use lodi_abstract_class, only :
     $       lodi_abstract

        use lodi_ns2d_class, only :
     $       lodi_ns2d

        use lodi_inflow_class , only :
     $       lodi_inflow

        use lodi_outflow_class, only :
     $       lodi_outflow

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none


        type(lodi_inflow)             :: lodi_inflow_tested
        type(lodi_outflow)            :: lodi_outflow_tested
        type(pmodel_eq)               :: p_model
        real(rkind)                   :: t        
        real(rkind), dimension(5,5,ne):: nodes
        real(rkind), dimension(5)     :: x_map
        real(rkind), dimension(5)     :: y_map
        integer(ikind)                :: i
        integer(ikind)                :: j
        real(rkind), dimension(ne,2)  :: test_data
        logical                       :: loc
        logical                       :: detailled
        logical                       :: test_validated


        !initialization of the nodes
        call initialize_nodes(nodes,x_map,y_map)


        !test parameters
        i = 3
        j = 3
        test_validated = .true.


        !test lodi inflow
        print '(''test lodi_inflow'')'
        print '(''---------------------------------'')'
        detailled = .false.

        test_data(1,1) = 0.0
        test_data(2,1) = 4.767578736
        test_data(3,1) = 7.151368103
        test_data(4,1) = 7.151368103

        test_data(1,2) = 0.0
        test_data(2,2) = 0.204258766
        test_data(3,2) = 0.306388148
        test_data(4,2) = 0.306388148
        
        call lodi_inflow_tested%ini()
        loc = test_lodi(
     $       lodi_inflow_tested, test_data, 
     $       p_model, t, nodes, x_map, y_map, i,j,
     $       detailled)

        test_validated = test_validated.and.loc
        if(.not.detailled) then
           print '(''---------------------------------'')'
           print '(''test_validated: '',L3)', loc
        end if
        print '()'
        print '()'


        !test lodi outflow
        print '(''test lodi_outflow'')'
        print '(''---------------------------------'')'

        detailled = .false.

        test_data(1,1) = 0.139321207
        test_data(2,1) = 0.160183735
        test_data(3,1) = 0.168993797
        test_data(4,1) = 5.625884172

        test_data(1,2) = 0.160849567
        test_data(2,2) =-1.205979506
        test_data(3,2) = 0.168993797
        test_data(4,2) =-1.272356148
        
        call lodi_outflow_tested%ini()
        loc = test_lodi(
     $       lodi_outflow_tested, test_data, 
     $       p_model, t, nodes, x_map, y_map, i,j,
     $       detailled)

        test_validated = test_validated.and.loc
        if(.not.detailled) then
           print '(''---------------------------------'')'
           print '(''test_validated: '',L3)', loc
        end if
        print '()'
        print '()'


        !test the computation of the time derivatives from
        !the LODI vector for the inflow b.c.
        print '(''test lodi_inflow: time_dev'')'
        print '(''---------------------------------'')'

        detailled = .false.

        test_data(1,1) = -29.7979273
        test_data(2,1) = -20.73328705
        test_data(3,1) = -22.04438499
        test_data(4,1) = -23.9488677

        test_data(1,2) =-1.276641286
        test_data(2,2) =-0.403835509
        test_data(3,2) =-1.428900765
        test_data(4,2) =-1.231197507
        
        call lodi_inflow_tested%ini()
        loc = test_lodi_timedev(
     $       lodi_inflow_tested, test_data, 
     $       p_model, t, nodes, x_map, y_map, i,j,
     $       detailled)

        test_validated = test_validated.and.loc
        if(.not.detailled) then
           print '(''---------------------------------'')'
           print '(''test_validated: '',L3)', loc
        end if
        print '()'
        print '()'


        contains

        subroutine initialize_nodes(nodes,x_map,y_map)

          implicit none

          real(rkind), dimension(:,:,:), intent(out) :: nodes
          real(rkind), dimension(:)    , intent(out) :: x_map
          real(rkind), dimension(:)    , intent(out) :: y_map

          integer(ikind) :: j
          real(rkind)    :: dx
          real(rkind)    :: dy


          !space steps
          dx = 0.5
          dy = 0.6

          do j=1, size(x_map,1)
             x_map(j) = (j-1)*dx
          end do

          do j=1, size(y_map,1)
             y_map(j) = (j-1)*dy
          end do


          !mass
          nodes(1,1,1) = 0.5
          nodes(2,1,1) = 0.2
          nodes(3,1,1) = 1.2
          nodes(4,1,1) = 5.0
          nodes(5,1,1) = 0.6

          nodes(1,2,1) = 3.0
          nodes(2,2,1) = 4.2
          nodes(3,2,1) = 11.0
          nodes(4,2,1) = 10.6
          nodes(5,2,1) = 5.2

          nodes(1,3,1) = -14.2
          nodes(2,3,1) = 23
          nodes(3,3,1) = 9.8
          nodes(4,3,1) = 3.4
          nodes(5,3,1) = 9.1

          nodes(1,4,1) = 2.45
          nodes(2,4,1) = 0.2
          nodes(3,4,1) = 9.0
          nodes(4,4,1) = 5.4
          nodes(5,4,1) =-2.3

          nodes(1,5,1) = 3.6
          nodes(2,5,1) = 0.1
          nodes(3,5,1) = 6.3
          nodes(4,5,1) = 8.9
          nodes(5,5,1) = -4.23


          !momentum-x
          nodes(1,1,2) = 7.012
          nodes(2,1,2) =-6.323
          nodes(3,1,2) = 3.012
          nodes(4,1,2) = 4.5
          nodes(5,1,2) = 9.6
                    
          nodes(1,2,2) = 4.26
          nodes(2,2,2) = 4.23
          nodes(3,2,2) = 4.5
          nodes(4,2,2) = 7.56
          nodes(5,2,2) = 7.21
                    
          nodes(1,3,2) = 0.23
          nodes(2,3,2) = 7.23
          nodes(3,3,2) = 3.1
          nodes(4,3,2) = 8.9
          nodes(5,3,2) = 9.3
                    
          nodes(1,4,2) = 8.23
          nodes(2,4,2) = -3.1
          nodes(3,4,2) = 6.03
          nodes(4,4,2) = 6.25
          nodes(5,4,2) = 5.12
                    
          nodes(1,5,2) = 3.2
          nodes(2,5,2) = 8.12
          nodes(3,5,2) = 8.9
          nodes(4,5,2) = 4.2
          nodes(5,5,2) = 7.8


          !momentum-y
          nodes(1,1,3) = 7.1
          nodes(2,1,3) = 1.052
          nodes(3,1,3) = 1.23
          nodes(4,1,3) = 7.89
          nodes(5,1,3) = 8.0
                    
          nodes(1,2,3) = 8.362
          nodes(2,2,3) = 4.56
          nodes(3,2,3) = 9.6
          nodes(4,2,3) = 8.96
          nodes(5,2,3) = -3.23
                    
          nodes(1,3,3) = 2.53
          nodes(2,3,3) = -3.23
          nodes(3,3,3) = 7.25
          nodes(4,3,3) = 1.02
          nodes(5,3,3) = 9.26
                    
          nodes(1,4,3) = 8.965
          nodes(2,4,3) = 4.789
          nodes(3,4,3) = 4.56
          nodes(4,4,3) = 3.012
          nodes(5,4,3) = -1.45
                    
          nodes(1,5,3) = 6.26
          nodes(2,5,3) = 5.201
          nodes(3,5,3) = 2.03
          nodes(4,5,3) = 7.89
          nodes(5,5,3) = 9.889


          !total energy
          nodes(1,1,4) = 6.23
          nodes(2,1,4) = 4.12
          nodes(3,1,4) = -3.6
          nodes(4,1,4) = -6.52
          nodes(5,1,4) = 9.57
                    
          nodes(1,2,4) = -0.12
          nodes(2,2,4) = 8.2
          nodes(3,2,4) = 1.2
          nodes(4,2,4) = 7.89
          nodes(5,2,4) = 5.62
                    
          nodes(1,3,4) = -6.23
          nodes(2,3,4) = 6.201
          nodes(3,3,4) = 6.7
          nodes(4,3,4) = 4.12
          nodes(5,3,4) = 1.29
                    
          nodes(1,4,4) = 1.2
          nodes(2,4,4) = 7.958
          nodes(3,4,4) = 1
          nodes(4,4,4) = -5.62
          nodes(5,4,4) = 0.36
                    
          nodes(1,5,4) = 9.6
          nodes(2,5,4) = 6.12
          nodes(3,5,4) = 8.9
          nodes(4,5,4) = 8.95
          nodes(5,5,4) = 6.3

          print '()'
          print '(''mass_density'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,1)
          end do
          print '()'

          print '()'
          print '(''momentum-x'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,2)
          end do
          print '()'

          print '()'
          print '(''momentum-y'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,3)
          end do
          print '()'

          print '()'
          print '(''total energy'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,4)
          end do
          print '()'
          print '()'

        end subroutine initialize_nodes


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
     $         int(var*10000.)-
     $         sign(int(abs(cst*10000.)),int(cst*10000.))).le.1
          
        end function is_test_validated


        function test_lodi_vector(lodi,test_data,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(ne), intent(in) :: lodi
          real(rkind), dimension(ne), intent(in) :: test_data
          logical                   , intent(in) :: detailled
          logical                                :: test_validated

          integer :: i
          logical :: loc

          test_validated = .true.

          do i=1,ne

             loc = is_test_validated(lodi(i), test_data(i), detailled)
             test_validated = test_validated.and.loc
             if(detailled) print '(''lodi('',I2,''): '',L3)', i, loc

          end do

        end function test_lodi_vector


        function test_lodi(
     $     lodi_tested, test_data, 
     $     p_model, t, nodes, x_map, y_map, i,j,
     $     detailled)
     $     result(test_validated)

          implicit none

          class(lodi_abstract)         , intent(in) :: lodi_tested
          real(rkind), dimension(4,2)  , intent(in) :: test_data
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated

          real(rkind), dimension(ne) :: lodi
          logical                    :: loc
          logical                    :: detailled_loc


          test_validated = .true.

          
          print '(''test compute_x_lodi'')'
          detailled_loc = detailled
          lodi = lodi_tested%compute_x_lodi(
     $         p_model, t, nodes, x_map, y_map, i,j, gradient_x_interior)
          loc  = test_lodi_vector(lodi,test_data(:,1),detailled_loc)
          test_validated = test_validated.and.loc
          if(.not.detailled_loc) print '(''test_validated: '',L3)', test_validated
          print '()'


          print '(''test compute_y_lodi'')'
          detailled_loc = detailled
          lodi = lodi_tested%compute_y_lodi(
     $         p_model, t, nodes, x_map, y_map, i,j, gradient_y_interior)
          loc  = test_lodi_vector(lodi,test_data(:,2),detailled_loc)
          test_validated = test_validated.and.loc
          if(.not.detailled_loc) print '(''test_validated: '',L3)', test_validated
          print '()'

        end function test_lodi


        function test_lodi_timedev(
     $     lodi_tested, test_data, 
     $     p_model, t, nodes, x_map, y_map, i,j,
     $     detailled)
     $     result(test_validated)

          implicit none

          class(lodi_ns2d)             , intent(in) :: lodi_tested
          real(rkind), dimension(4,2)  , intent(in) :: test_data
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated

          real(rkind), dimension(ne) :: timedev
          logical                    :: loc
          logical                    :: detailled_loc


          test_validated = .true.

          
          print '(''test compute_x_timedev'')'
          detailled_loc = detailled
          timedev = lodi_tested%compute_x_timedev(
     $         p_model, t, nodes, x_map, y_map, i,j, gradient_x_interior)
          loc  = test_lodi_vector(timedev,test_data(:,1),detailled_loc)
          test_validated = test_validated.and.loc
          if(.not.detailled_loc) print '(''test_validated: '',L3)', test_validated
          print '()'


          print '(''test compute_y_timedev'')'
          detailled_loc = detailled
          timedev = lodi_tested%compute_y_timedev(
     $         p_model, t, nodes, x_map, y_map, i,j, gradient_y_interior)
          loc  = test_lodi_vector(timedev,test_data(:,2),detailled_loc)
          test_validated = test_validated.and.loc
          if(.not.detailled_loc) print '(''test_validated: '',L3)', test_validated
          print '()'

        end function test_lodi_timedev

      end program test_lodi_poinsot_ns2d

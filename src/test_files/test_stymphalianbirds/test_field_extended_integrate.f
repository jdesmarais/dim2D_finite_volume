      !to test whether the integration on the two fields
      !leads to teh same results
      program test_field_extended_integration

        use field_class, only :
     $       field

        use field_extended_class, only : 
     $       field_extended

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        
        type(field)          :: field_used
        type(field_extended) :: field_extended_used

        real(rkind), dimension(nx,ny,ne) :: nodes_test
        real(rkind), dimension(nx,ny,ne) :: nodes_field
        real(rkind), dimension(nx,ny,ne) :: nodes_field_extended

        logical :: detailled
        logical :: test_validated

        real(rkind) :: dt

        
        detailled = .true.


        !verify the inputs
        if((nx.ne.10).or.(ny.ne.10).or.(ne.ne.4)) then
           print '(''test_field_extended_integrate'')'
           print '(''nx.eq.10: '',L1)', nx.eq.10
           print '(''ny.eq.10: '',L1)', ny.eq.10
           print '(''ne.eq.4 : '',L1)', ne.eq.4
           stop ''
        end if


        !initialize the fields
        call field_used%ini()
        call field_extended_used%ini()


        !set the nodes
        call initialize_nodes(nodes_test)
        call field_used%set_nodes(nodes_test)
        call field_extended_used%set_nodes(nodes_test)


        !integrate in time
        dt = 0.0001
        call field_used%integrate(dt)
        call field_extended_used%integrate(dt)


        !get the nodes
        nodes_field = field_used%get_nodes()
        nodes_field_extended = field_extended_used%get_nodes()


        !compare the nodes
        test_validated = compare_nodes(nodes_field,nodes_field_extended,detailled)
        print '(''test_field_extended_integrate: '',L1)', test_validated
        print '()'

        contains


        subroutine initialize_nodes(nodes)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes

          integer(ikind) :: i,j
          integer        :: k


          do k=1,ne
             do j=1,ny
                do i=1,nx
                   nodes(i,j,k) = 0.0d0
                end do
             end do
          end do
          
          !mass----------------
          nodes(1,1,1) =  2.3d0
          nodes(2,1,1) =  1.02d0
          nodes(3,1,1) =  3.2d0
                          
          nodes(1,2,1) =  1.2d0
          nodes(2,2,1) =  8.6d0
          nodes(3,2,1) =  6.13d0
                          
          nodes(1,3,1) =  0.23d0
          nodes(2,3,1) =  4.5d0
          nodes(3,3,1) =  7.13d0
                          
          nodes(1,4,1) =  8.5d0
          nodes(2,4,1) =  7.8d0
          nodes(3,4,1) =  1.5d0
                          
          nodes(1,5,1) =  0.2d0
          nodes(2,5,1) =  3.6d0
          nodes(3,5,1) =  9.23d0

          !momentum-x----------
          nodes(1,1,2) = -6.045d0
          nodes(2,1,2) =  8.125d0
          nodes(3,1,2) =  3.26d0
                    
          nodes(1,2,2) = -6.3d0
          nodes(2,2,2) =  7.98d0
          nodes(3,2,2) = -5.89d0
                    
          nodes(1,3,2) = -0.15d0
          nodes(2,3,2) =  6.213d0
          nodes(3,3,2) =  2.105d0
                    
          nodes(1,4,2) =  8.23d0
          nodes(2,4,2) = -3.012d0
          nodes(3,4,2) =  6.213d0
                    
          nodes(1,5,2) = -1.23d0
          nodes(2,5,2) =  7.8d0
          nodes(3,5,2) = -0.15d0

          !momentum-y----------
          nodes(1,1,3) =  2.01d0
          nodes(2,1,3) =  3.25d0
          nodes(3,1,3) =  6.2d0
                    
          nodes(1,2,3) =  7.135d0
          nodes(2,2,3) = -2.01d0
          nodes(3,2,3) = -3.06d0
                    
          nodes(1,3,3) =  9.46d0
          nodes(2,3,3) = -9.16d0
          nodes(3,3,3) =  4.12d0
                    
          nodes(1,4,3) =  2.13d0
          nodes(2,4,3) = -2.15d0
          nodes(3,4,3) = -3.25d0
                    
          nodes(1,5,3) =  6.1023d0
          nodes(2,5,3) =  5.23d0
          nodes(3,5,3) =  1.12d0

          !total_energy--------
          nodes(1,1,4) =  20.1d0
          nodes(2,1,4) =  895.26d0
          nodes(3,1,4) =  961.23d0

          nodes(1,2,4) =  78.256d0
          nodes(2,2,4) =  8.45d0
          nodes(3,2,4) =  7.4d0

          nodes(1,3,4) =  256.12d0
          nodes(2,3,4) =  163.48d0
          nodes(3,3,4) =  9.56d0
                    
          nodes(1,4,4) =  56.12d0
          nodes(2,4,4) =  7.89d0
          nodes(3,4,4) =  629.12d0
                    
          nodes(1,5,4) =  102.3d0
          nodes(2,5,4) =  231.02d0
          nodes(3,5,4) =  7.123d0

          
          !fill the nodes (i in [4,5])x(j in [1,5]) by using
          !the symmetry along the x-axis

          do k=1, ne
             do j=1,5
                do i=1,2
                   nodes(6-i,j,k) =   nodes(i,j,k)
                end do
             end do
          end do

        end subroutine initialize_nodes


        function compare_nodes(
     $     nodes_field,
     $     nodes_field_extended,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes_field
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes_field_extended
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated

          integer(ikind) :: i,j
          integer        :: k
          logical        :: test_loc

          test_validated = .true.

          do k=1,ne
             do j=1,ny
                do i=1,nx

                   test_loc = is_test_validated(
     $                  nodes_field(i,j,k),
     $                  nodes_field_extended(i,j,k),
     $                  detailled)

                   test_validated = test_validated.and.test_loc

                   if(detailled.and.(.not.test_loc)) then

                      print '(''['',3I2,'']'',F10.4,'' -> '',F10.4)', 
     $                     i,j,k,
     $                     nodes_field_extended(i,j,k),
     $                     nodes_field(i,j,k)

                   end if

                end do
             end do
          end do

        end function compare_nodes


        function is_test_validated(var,cst,detailled)
     $     result(test_validated)

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

      end program test_field_extended_integration

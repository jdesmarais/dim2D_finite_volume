      program test_field_extended_integration

        use field_class, only :
     $       field

        use field_extended_class, only : 
     $       field_extended

        use parameters_input, only :
     $       nx,ny,ne


        implicit none

        
        type(field)          :: field_used
        type(field_extended) :: field_extended_used

        real(rkind), dimension(nx,ny,ne) :: nodes_test


        !verify the inputs
        if((nx.ne.6).or.(ny.ne.6).or.(ne.ne.4)) then
           
        end if


        call field_used%ini()
        call field_extended_used%ini()


        !set the nodes
        call initialize_nodes(nodes_test)
        call field_used%set_nodes(nodes_test)
        call field_extended_used%set_nodes(nodes_test)


        !integrate in time
        


        contains


        subroutine initialize_nodes(nodes)

          implicit none

          
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


      end program test_field_extended_integration

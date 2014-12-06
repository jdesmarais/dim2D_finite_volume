      program test_hedstrom_xy_corners

        use hedstrom_xy_corners_module, only :
     $       compute_n_timedev_with_openbc

        use parameters_input, only :
     $       nx,ny,ne
        
        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_fd_ncoords_module, only :
     $     gradient_n1_oneside_L0,
     $     gradient_n1_oneside_L1,
     $     gradient_n1_oneside_R1,
     $     gradient_n1_oneside_R0,
     $     gradient_n2_oneside_L0,
     $     gradient_n2_oneside_L1,
     $     gradient_n2_oneside_R1,
     $     gradient_n2_oneside_R0


        implicit none

        real(rkind), dimension(nx,ny,ne)   :: nodes
        real(rkind)                        :: dx
        real(rkind)                        :: dy

        logical                            :: test_loc
        logical                            :: test_validated
        logical                            :: detailled

        test_validated = .true.
        
        if((nx.ne.7).or.(ny.ne.5).or.(ne.ne.3)) then
           stop 'the test requires (nx,ny,ne)=(7,5,3)'
        end if

        print '()'
        print '(''*************************'')'
        print '(''WARNING: use wave2d model'')'
        print '(''*************************'')'
        print '()'

        detailled = .false.
        call initialize_nodes(nodes,dx,dy,detailled)


        !test the n_gradients
        detailled = .true.
        test_loc = test_compute_gradients_at_edges(nodes,dx,dy,detailled)
        test_validated = test_validated.and.test_loc
        print '(''test gradients at the corners: '',L1)', test_validated
        print '()'


        !test the time derivatives
        detailled = .true.
        test_loc = test_apply_bc_on_timedev(nodes,dx,dy,detailled)
        test_validated = test_validated.and.test_loc
        print '(''test apply_bc_on_timedev on corners: '',L1)', test_validated
        print '()'           


        print '(''test_hedstrom_xy_corners: '',L1)', test_validated
        print '(''--------------------------------'')'
        print '()'

        contains


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
     $         nint(var*1e5)-
     $         nint(cst*1e5)).le.1
          
        end function is_test_validated        


        subroutine initialize_nodes(nodes,dx,dy,detailled)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          real(rkind)                     , intent(out) :: dx
          real(rkind)                     , intent(out) :: dy
          logical                         , intent(in)  :: detailled
          
          integer :: k
          integer(ikind) :: j

          dx=0.5
          dy=0.6

          !position
          nodes(1,1,1)=  0.5d0
          nodes(2,1,1)=  0.2d0
          nodes(3,1,1)=  1.2d0
          nodes(4,1,1)=  5.0d0
          nodes(5,1,1)=  0.6d0
          nodes(6,1,1)= -3.6d0
          nodes(7,1,1)=-6.52d0

          nodes(1,2,1)= 3.0d0
          nodes(2,2,1)= 4.2d0
          nodes(3,2,1)=11.0d0
          nodes(4,2,1)=10.6d0
          nodes(5,2,1)= 5.2d0
          nodes(6,2,1)= 1.2d0
          nodes(7,2,1)=7.89d0

          nodes(1,3,1)=-14.2d0
          nodes(2,3,1)=   23d0
          nodes(3,3,1)=  9.8d0
          nodes(4,3,1)=  3.4d0
          nodes(5,3,1)=  9.1d0
          nodes(6,3,1)=  6.7d0
          nodes(7,3,1)= 4.12d0

          nodes(1,4,1)= 2.45d0
          nodes(2,4,1)=  0.2d0
          nodes(3,4,1)=  9.0d0
          nodes(4,4,1)=  5.4d0
          nodes(5,4,1)= -2.3d0
          nodes(6,4,1)=  1.0d0
          nodes(7,4,1)=-5.62d0

          nodes(1,5,1)=  3.6d0
          nodes(2,5,1)=  0.1d0
          nodes(3,5,1)=  6.3d0
          nodes(4,5,1)=  8.9d0
          nodes(5,5,1)=-4.23d0
          nodes(6,5,1)=  8.9d0
          nodes(7,5,1)= 8.95d0


          !velocity_x
          nodes(1,1,2)= 7.012d0
          nodes(2,1,2)=-6.323d0
          nodes(3,1,2)= 3.012d0
          nodes(4,1,2)=   4.5d0
          nodes(5,1,2)=   9.6d0
          nodes(6,1,2)=  9.57d0
          nodes(7,1,2)= 6.012d0
                    
          nodes(1,2,2)=4.26d0
          nodes(2,2,2)=4.23d0
          nodes(3,2,2)= 4.5d0
          nodes(4,2,2)=7.56d0
          nodes(5,2,2)=7.21d0
          nodes(6,2,2)=5.62d0
          nodes(7,2,2)=3.62d0
                    
          nodes(1,3,2)=0.23d0
          nodes(2,3,2)=7.23d0
          nodes(3,3,2)= 3.1d0
          nodes(4,3,2)= 8.9d0
          nodes(5,3,2)= 9.3d0
          nodes(6,3,2)=1.29d0
          nodes(7,3,2)=7.26d0
                    
          nodes(1,4,2)=8.23d0
          nodes(2,4,2)=-3.1d0
          nodes(3,4,2)=6.03d0
          nodes(4,4,2)=6.25d0
          nodes(5,4,2)=5.12d0
          nodes(6,4,2)=0.36d0
          nodes(7,4,2)=6.89d0
                    
          nodes(1,5,2)= 3.2d0
          nodes(2,5,2)=8.12d0
          nodes(3,5,2)= 8.9d0
          nodes(4,5,2)= 4.2d0
          nodes(5,5,2)= 7.8d0
          nodes(6,5,2)= 6.3d0
          nodes(7,5,2)= 1.2d0


          !velocity_y
          nodes(1,1,3)=  7.1d0
          nodes(2,1,3)=1.052d0
          nodes(3,1,3)= 1.23d0
          nodes(4,1,3)= 7.89d0
          nodes(5,1,3)=  8.0d0
          nodes(6,1,3)= 6.23d0
          nodes(7,1,3)= 4.12d0
                    
          nodes(1,2,3)=8.362d0
          nodes(2,2,3)= 4.56d0
          nodes(3,2,3)=  9.6d0
          nodes(4,2,3)= 8.96d0
          nodes(5,2,3)=-3.23d0
          nodes(6,2,3)=-0.12d0
          nodes(7,2,3)=  8.2d0
                    
          nodes(1,3,3)= 2.53d0
          nodes(2,3,3)=-3.23d0
          nodes(3,3,3)= 7.25d0
          nodes(4,3,3)= 1.02d0
          nodes(5,3,3)= 9.26d0
          nodes(6,3,3)=-6.23d0
          nodes(7,3,3)=6.201d0
                    
          nodes(1,4,3)=8.965d0
          nodes(2,4,3)=4.789d0
          nodes(3,4,3)= 4.56d0
          nodes(4,4,3)=3.012d0
          nodes(5,4,3)=-1.45d0
          nodes(6,4,3)=  1.2d0
          nodes(7,4,3)=7.958d0
                    
          nodes(1,5,3)= 6.26d0
          nodes(2,5,3)=5.201d0
          nodes(3,5,3)= 2.03d0
          nodes(4,5,3)= 7.89d0
          nodes(5,5,3)=9.889d0
          nodes(6,5,3)=  9.6d0
          nodes(7,5,3)= 6.12d0

          if(detailled) then
             print '()'
             do k=1,ne
                print '(''nodes(:,:,'',I1,'')'')', k
                print '(''-------------------------------'')'
                do j=1,5
                   print '(7F8.3)', nodes(:,6-j,k)
                end do
                print '()'
             end do
             print '()'
          end if

        end subroutine initialize_nodes


        function test_compute_gradients_at_edges(
     $     nodes,dx,dy,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)  :: nodes
          real(rkind)                     , intent(in)  :: dx
          real(rkind)                     , intent(in)  :: dy
          logical                         , intent(in)  :: detailled
          logical                                       :: test_validated

          real(rkind), dimension(nx,ny,ne) :: gradients_n
          type(pmodel_eq)                  :: p_model
          integer(ikind)                   :: i,j

          !============================================================
          !computation of the gradients
          !============================================================
          !south layers
          j=1
          i=1
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n2_oneside_L0,dx,dy)

          i=2
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n2_oneside_L0,dx,dy)
          
          i=nx-1
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n1_oneside_R0,dx,dy)
          
          i=nx
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n1_oneside_R0,dx,dy)
          
          
          j=2
          i=1
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n2_oneside_L0,dx,dy)
          
          i=2
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n2_oneside_L1,dx,dy)
          
          i=nx-1
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n1_oneside_R1,dx,dy)
          
          i=nx
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n1_oneside_R0,dx,dy)


          !north layers
          j=ny-1
          i=1
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n1_oneside_L0,dx,dy)
          
          i=2
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n1_oneside_L1,dx,dy)
          
          i=nx-1
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n2_oneside_R1,dx,dy)
          
          i=nx
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n2_oneside_R0,dx,dy)


          j=ny
          i=1
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n1_oneside_L0,dx,dy)
          
          i=2
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n1_oneside_L0,dx,dy)
          
          i=nx-1
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n2_oneside_R0,dx,dy)
          
          i=nx
          gradients_n(i,j,:) = p_model%compute_n_gradient(
     $         nodes,i,j,gradient_n2_oneside_R0,dx,dy)

          !============================================================
          !test of the gradients
          !============================================================
          test_validated = test_gradients(gradients_n,detailled)

        end function test_compute_gradients_at_edges


        function test_gradients(gradients_n,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: gradients_n
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated


          real(rkind), dimension(nx,ny,ne)             :: test_data

          !position
          test_data(1,1,1) = 4.73736456d0
          test_data(2,1,1) =   13.82798d0
          test_data(6,1,1) = -11.267245d0
          test_data(7,1,1) =   -9.88445d0

          test_data(1,2,1) = 25.607376d0
          test_data(2,2,1) =  5.953715d0
          test_data(6,2,1) =  -9.99968d0
          test_data(7,2,1) =  1.523639d0

          test_data(1,4,1) = 26.3115788d0
          test_data(2,4,1) =   3.969143d0
          test_data(6,4,1) =   -0.09603d0
          test_data(7,4,1) =-15.7741436d0

          test_data(1,5,1) =-4.35325392d0
          test_data(2,5,1) =   11.39528d0
          test_data(6,5,1) =   14.34013d0
          test_data(7,5,1) =   10.17893d0


          !velocity_x
          test_data(1,1,2) =-3.561986d0
          test_data(2,1,2) = 13.85743d0
          test_data(6,1,2) =  3.02167d0
          test_data(7,1,2) = 0.501905d0

          test_data(1,2,2) = 3.80269533d0
          test_data(2,2,2) =    -2.5044d0
          test_data(6,2,2) =   -2.10493d0
          test_data(7,2,2) =   2.983259d0

          test_data(1,4,2) = -1.2803688d0
          test_data(2,4,2) =   -0.06402d0
          test_data(6,4,2) =   -5.18549d0
          test_data(7,4,2) =   7.170065d0

          test_data(1,5,2) =-8.0663234d0
          test_data(2,5,2) =-  2.67597d0
          test_data(6,5,2) =  1.510835d0
          test_data(7,5,2) =   1.07551d0


          !velocity_y
          test_data(1,1,3) = -3.25214d0
          test_data(2,1,3) = 10.94459d0
          test_data(6,1,3) = 12.11229d0
          test_data(7,1,3) = 5.428764d0

          test_data(1,2,3) =-14.842034d0
          test_data(2,2,3) = 0.096028d0
          test_data(6,2,3) = -3.29055d0
          test_data(7,2,3) = 18.47572d0

          test_data(1,4,3) = -15.6141d0
          test_data(2,4,3) = 0.633783d0
          test_data(6,4,3) = -2.01018d0
          test_data(7,4,3) = 18.16587d0

          test_data(1,5,3) = -1.88342d0
          test_data(2,5,3) = -0.82072d0
          test_data(6,5,3) = 14.14808d0
          test_data(7,5,3) = 6.299414d0


          !test the gradients          
          test_validated = check_corners(gradients_n,test_data,detailled)

        end function test_gradients


        function test_apply_bc_on_timedev(nodes,dx,dy,detailled) result(test_validated)

          implicit none

          
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated

          type(bc_operators)                 :: bc_used
          type(pmodel_eq)                    :: p_model
          real(rkind), dimension(nx)         :: x_map
          real(rkind), dimension(ny)         :: y_map
          real(rkind)                        :: t
          real(rkind), dimension(nx+1,ny,ne) :: flux_x
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          real(rkind), dimension(nx,ny,ne)   :: timedev
          real(rkind), dimension(nx,ny,ne)   :: test_data
          integer(ikind)                     :: i,j
        

          !initialization of the maps
          do i=1, nx
             x_map(i) = (i-1)*dx
          end do
          
          do j=1,ny
             y_map(j) = (j-1)*dy
          end do
          
          !initialization of the b.c.
          call bc_used%ini(p_model)
          
          !computation of the b.c.
          call bc_used%apply_bc_on_timedev(
     $         p_model,
     $         t,nodes,x_map,y_map,
     $         flux_x,flux_y,
     $         timedev)

          !position
          test_data(1,1,1) =-0.010118481d0
          test_data(2,1,1) =   3.9207078d0
          test_data(6,1,1) = 0.604900933d0
          test_data(7,1,1) = 0.800078953d0
          
          test_data(1,2,1) = 2.225172995d0
          test_data(2,2,1) = 0.531342192d0
          test_data(6,2,1) = 1.354755166d0
          test_data(7,2,1) =-1.559808018d0

          test_data(1,4,1) =  4.555881949d0
          test_data(2,4,1) =  0.434465433d0
          test_data(6,4,1) = -0.624010158d0
          test_data(7,4,1) =  4.211169629d0

          test_data(1,5,1) = -1.090653137d0
          test_data(2,5,1) =  1.260427419d0
          test_data(6,5,1) = -0.408451101d0
          test_data(7,5,1) = -0.620509123d0

          !velocity_x
          test_data(1,1,2) = -0.007154846d0
          test_data(2,1,2) =  2.772359072d0
          test_data(6,1,2) = -0.427729552d0
          test_data(7,1,2) = -0.565741253d0

          test_data(1,2,2) =  1.573434914d0
          test_data(2,2,2) =  0.375715667d0
          test_data(6,2,2) = -0.957956565d0
          test_data(7,2,2) =  1.102950827d0

          test_data(1,4,2) =  3.221495021d0
          test_data(2,4,2) =  0.307213454d0
          test_data(6,4,2) =  0.441241815d0
          test_data(7,4,2) = -2.977746602d0

          test_data(1,5,2) = -0.771208229d0
          test_data(2,5,2) =  0.891256775d0
          test_data(6,5,2) =  0.288818544d0
          test_data(7,5,2) =  0.438766209d0

          !velocity_y
          test_data(1,1,3) = -0.007154846d0
          test_data(2,1,3) =  2.772359072d0
          test_data(6,1,3) =  0.427729552d0
          test_data(7,1,3) =  0.565741253d0

          test_data(1,2,3) =  1.573434914d0
          test_data(2,2,3) =  0.375715667d0
          test_data(6,2,3) =  0.957956565d0
          test_data(7,2,3) = -1.102950827d0

          test_data(1,4,3) = -3.221495021d0
          test_data(2,4,3) = -0.307213454d0
          test_data(6,4,3) =  0.441241815d0
          test_data(7,4,3) = -2.977746602d0

          test_data(1,5,3) =  0.771208229d0
          test_data(2,5,3) = -0.891256775d0
          test_data(6,5,3) =  0.288818544d0
          test_data(7,5,3) =  0.438766209d0

          !test the gradients          
          test_validated = check_corners(timedev,test_data,detailled)

        end function test_apply_bc_on_timedev

      
        function check_corners(gradients_n, test_data, detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: gradients_n
          real(rkind), dimension(nx,ny,ne), intent(in) :: test_data
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated


          integer(ikind) :: i,j,k


          test_validated = .true.
        

          do k=1,3
             do j=1,2
                do i=1,2
                   call check_data(
     $                  gradients_n,test_data,
     $                  i,j,k,
     $                  test_validated,
     $                  detailled)
                end do

                do i=6,7
                   call check_data(
     $                  gradients_n,test_data,
     $                  i,j,k,
     $                  test_validated,
     $                  detailled)                   
                end do
             end do
             
             do j=4,5
                do i=1,2
                   call check_data(
     $                  gradients_n,test_data,
     $                  i,j,k,
     $                  test_validated,
     $                  detailled)
                end do

                do i=6,7
                   call check_data(
     $                  gradients_n,test_data,
     $                  i,j,k,
     $                  test_validated,
     $                  detailled)                   
                end do
             end do
          end do

          if(.not.detailled) print '(''test_validated: '',L3)',
     $         test_validated

        end function check_corners


        subroutine check_data(
     $     gradients,
     $     test_data,
     $     i,j,k,
     $     test_validated,
     $     detailled)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)    :: gradients
          real(rkind), dimension(nx,ny,ne), intent(in)    :: test_data
          integer(ikind)                  , intent(in)    :: i,j,k
          logical                         , intent(inout) :: test_validated
          logical                         , intent(in)    :: detailled
          
          logical :: test_loc

          test_loc = is_test_validated(
     $         gradients(i,j,k),
     $         test_data(i,j,k),
     $         .false.)

          if(detailled.and.(.not.test_loc)) then
             print '(''['',3I2'']: '',F8.3,'' -> '',F8.3)',
     $            i,j,k,
     $            gradients(i,j,k),
     $            test_data(i,j,k)
          end if

          test_validated = test_validated.and.test_loc
          
        end subroutine check_data

      end program test_hedstrom_xy_corners

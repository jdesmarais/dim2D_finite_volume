      program test_hedstrom_xy_corners

        use bc_operators_class, only :
     $     bc_operators

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

        type(bc_operators)                 :: bc_used
        type(pmodel_eq)                    :: p_model
        real(rkind), dimension(nx,ny,ne)   :: nodes
        real(rkind), dimension(nx)         :: x_map
        real(rkind), dimension(ny)         :: y_map
        real(rkind)                        :: dx
        real(rkind)                        :: dy
        real(rkind)                        :: t
        real(rkind), dimension(nx,ny,ne)   :: gradients_n
        real(rkind), dimension(nx+1,ny,ne) :: flux_x
        real(rkind), dimension(nx,ny+1,ne) :: flux_y
        real(rkind), dimension(nx,ny,ne)   :: timedev
        logical                            :: detailled
        integer(ikind) :: i,j
        
        if((nx.ne.7).or.(ny.ne.5).or.(ne.ne.3)) then
           stop 'the test requires (nx,ny,ne)=(7,5,3)'
        end if

        print '()'
        print '(''*************************'')'
        print '(''WARNING: use wave2d model'')'
        print '(''*************************'')'
        print '()'


        call initialize_nodes(nodes,dx,dy)


        !test the n_gradients
        call compute_gradients_at_edges(
     $       gradients_n,p_model,nodes,dx,dy)

        print '(''test gradients at the corners'')'
        print '(''-----------------------------'')'
        detailled = .false.
        call test_gradients(gradients_n,detailled)
        print '()'
        print '()'


        !test the time derivatives
        do i=1, nx
           x_map(i) = (i-1)*dx
        end do

        do j=1,ny
           y_map(j) = (j-1)*dy
        end do


        call bc_used%ini(p_model)
        call bc_used%apply_bc_on_timedev(
     $       p_model,
     $       t,nodes,x_map,y_map,
     $       flux_x,flux_y,
     $       timedev)
        print '(''test time_dev at the corners '')'
        print '(''-----------------------------'')'
        detailled = .false.
        call test_timedev(timedev,detailled)
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
     $         int(var*10000.)-
     $         sign(int(abs(cst*10000.)),int(cst*10000.))).le.1
          
        end function is_test_validated        


        subroutine initialize_nodes(nodes,dx,dy)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          real(rkind)                     , intent(out) :: dx
          real(rkind)                     , intent(out) :: dy
          
          integer :: k
          integer(ikind) :: j

          dx=0.5
          dy=0.6

          !position
          nodes(1,1,1)=0.5
          nodes(2,1,1)=0.2
          nodes(3,1,1)=1.2
          nodes(4,1,1)=5.0
          nodes(5,1,1)=0.6
          nodes(6,1,1)=-3.6
          nodes(7,1,1)=-6.52

          nodes(1,2,1)=3.0
          nodes(2,2,1)=4.2
          nodes(3,2,1)=11.0
          nodes(4,2,1)=10.6
          nodes(5,2,1)=5.2
          nodes(6,2,1)=1.2
          nodes(7,2,1)=7.89

          nodes(1,3,1)=-14.2
          nodes(2,3,1)=23
          nodes(3,3,1)=9.8
          nodes(4,3,1)=3.4
          nodes(5,3,1)=9.1
          nodes(6,3,1)=6.7
          nodes(7,3,1)=4.12

          nodes(1,4,1)=2.45
          nodes(2,4,1)=0.2
          nodes(3,4,1)=9.0
          nodes(4,4,1)=5.4
          nodes(5,4,1)=-2.3
          nodes(6,4,1)=1.0
          nodes(7,4,1)=-5.62

          nodes(1,5,1)=3.6
          nodes(2,5,1)=0.1
          nodes(3,5,1)=6.3
          nodes(4,5,1)=8.9
          nodes(5,5,1)=-4.23
          nodes(6,5,1)=8.9
          nodes(7,5,1)=8.95


          !velocity_x
          nodes(1,1,2)=7.012
          nodes(2,1,2)=-6.323
          nodes(3,1,2)=3.012
          nodes(4,1,2)=4.5
          nodes(5,1,2)=9.6
          nodes(6,1,2)=9.57
          nodes(7,1,2)=6.012
                    
          nodes(1,2,2)=4.26
          nodes(2,2,2)=4.23
          nodes(3,2,2)=4.5
          nodes(4,2,2)=7.56
          nodes(5,2,2)=7.21
          nodes(6,2,2)=5.62
          nodes(7,2,2)=3.62
                    
          nodes(1,3,2)=0.23
          nodes(2,3,2)=7.23
          nodes(3,3,2)=3.1
          nodes(4,3,2)=8.9
          nodes(5,3,2)=9.3
          nodes(6,3,2)=1.29
          nodes(7,3,2)=7.26
                    
          nodes(1,4,2)=8.23
          nodes(2,4,2)=-3.1
          nodes(3,4,2)=6.03
          nodes(4,4,2)=6.25
          nodes(5,4,2)=5.12
          nodes(6,4,2)=0.36
          nodes(7,4,2)=6.89
                    
          nodes(1,5,2)=3.2
          nodes(2,5,2)=8.12
          nodes(3,5,2)=8.9
          nodes(4,5,2)=4.2
          nodes(5,5,2)=7.8
          nodes(6,5,2)=6.3
          nodes(7,5,2)=1.2


          !velocity_y
          nodes(1,1,3)=7.1
          nodes(2,1,3)=1.052
          nodes(3,1,3)=1.23
          nodes(4,1,3)=7.89
          nodes(5,1,3)=8.0
          nodes(6,1,3)=6.23
          nodes(7,1,3)=4.12
                    
          nodes(1,2,3)=8.362
          nodes(2,2,3)=4.56
          nodes(3,2,3)=9.6
          nodes(4,2,3)=8.96
          nodes(5,2,3)=-3.23
          nodes(6,2,3)=-0.12
          nodes(7,2,3)=8.2
                    
          nodes(1,3,3)=2.53
          nodes(2,3,3)=-3.23
          nodes(3,3,3)=7.25
          nodes(4,3,3)=1.02
          nodes(5,3,3)=9.26
          nodes(6,3,3)=-6.23
          nodes(7,3,3)=6.201
                    
          nodes(1,4,3)=8.965
          nodes(2,4,3)=4.789
          nodes(3,4,3)=4.56
          nodes(4,4,3)=3.012
          nodes(5,4,3)=-1.45
          nodes(6,4,3)=1.2
          nodes(7,4,3)=7.958
                    
          nodes(1,5,3)=6.26
          nodes(2,5,3)=5.201
          nodes(3,5,3)=2.03
          nodes(4,5,3)=7.89
          nodes(5,5,3)=9.889
          nodes(6,5,3)=9.6
          nodes(7,5,3)=6.12

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

        end subroutine initialize_nodes


        subroutine compute_gradients_at_edges(
     $     gradients_n, p_model, nodes, dx, dy)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: gradients_n
          type(pmodel_eq)                 , intent(in)  :: p_model
          real(rkind), dimension(nx,ny,ne), intent(in)  :: nodes
          real(rkind)                     , intent(in)  :: dx
          real(rkind)                     , intent(in)  :: dy

          integer(ikind) :: i,j

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

        end subroutine compute_gradients_at_edges


        subroutine test_gradients(gradients_n,detailled)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: gradients_n
          logical                         , intent(in)  :: detailled
          real(rkind), dimension(nx,ny,ne)              :: test_data


          !position
          test_data(1,1,1) = 4.73736456
          test_data(2,1,1) = 13.82798
          test_data(6,1,1) =-11.2672
          test_data(7,1,1) =-9.88445

          test_data(1,2,1) = 25.607376
          test_data(2,2,1) = 5.953715
          test_data(6,2,1) = -9.99968
          test_data(7,2,1) = 1.523639

          test_data(1,4,1) = 26.3115788
          test_data(2,4,1) = 3.969143
          test_data(6,4,1) = -0.09603
          test_data(7,4,1) =-15.7741

          test_data(1,5,1) =-4.35325392
          test_data(2,5,1) = 11.39528
          test_data(6,5,1) = 14.34013
          test_data(7,5,1) = 10.17893


          !velocity_x
          test_data(1,1,2) =-3.561986
          test_data(2,1,2) = 13.85743
          test_data(6,1,2) = 3.02167
          test_data(7,1,2) = 0.501905

          test_data(1,2,2) = 3.80269533
          test_data(2,2,2) = -2.5044
          test_data(6,2,2) = -2.10493
          test_data(7,2,2) = 2.983259

          test_data(1,4,2) = -1.2803688
          test_data(2,4,2) = -0.06402
          test_data(6,4,2) = -5.18549
          test_data(7,4,2) = 7.170065

          test_data(1,5,2) =-8.0663234
          test_data(2,5,2) =-2.67597
          test_data(6,5,2) = 1.510835
          test_data(7,5,2) = 1.07551


          !velocity_y
          test_data(1,1,3) =-3.25214
          test_data(2,1,3) = 10.94459
          test_data(6,1,3) = 12.11229
          test_data(7,1,3) = 5.428764

          test_data(1,2,3) = -14.842
          test_data(2,2,3) =  0.096028
          test_data(6,2,3) = -3.29055
          test_data(7,2,3) = 18.47572

          test_data(1,4,3) = -15.6141
          test_data(2,4,3) = 0.633783
          test_data(6,4,3) = -2.01018
          test_data(7,4,3) = 18.16587

          test_data(1,5,3) = -1.88342
          test_data(2,5,3) =-0.82072
          test_data(6,5,3) = 14.14808
          test_data(7,5,3) = 6.299414


          !test the gradients          
          call check_corners(gradients_n,test_data,detailled)          

        end subroutine test_gradients


        subroutine test_timedev(timedev,detailled)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: timedev
          logical                         , intent(in)  :: detailled
          real(rkind), dimension(nx,ny,ne)              :: test_data


          !position
          test_data(1,1,1) =-0.01430969
          test_data(2,1,1) = 5.544718
          test_data(6,1,1) = 0.855459
          test_data(7,1,1) = 1.131483

          test_data(1,2,1) = 3.14686983
          test_data(2,2,1) = 0.751431
          test_data(6,2,1) = 1.915913
          test_data(7,2,1) =-2.2059

          test_data(1,4,1) = 6.44299004
          test_data(2,4,1) = 0.614427
          test_data(6,4,1) =-0.88248
          test_data(7,4,1) = 5.955493

          test_data(1,5,1) =-1.54241646
          test_data(2,5,1) = 1.782514
          test_data(6,5,1) =-0.57764
          test_data(7,5,1) =-0.87753


          !velocity_x
          test_data(1,1,2) =-0.0101185
          test_data(2,1,2) = 3.920708
          test_data(6,1,2) =-0.6049
          test_data(7,1,2) =-0.80008 

          test_data(1,2,2) = 2.225173 
          test_data(2,2,2) = 0.531342
          test_data(6,2,2) = -1.35476
          test_data(7,2,2) = 1.559808

          test_data(1,4,2) = 4.55588195 
          test_data(2,4,2) = 0.434465
          test_data(6,4,2) = 0.62401
          test_data(7,4,2) =-4.21117  

          test_data(1,5,2) =-1.0906531 
          test_data(2,5,2) = 1.260427
          test_data(6,5,2) = 0.408451
          test_data(7,5,2) = 0.620509


          !velocity_y
          test_data(1,1,3) =-0.01012 
          test_data(2,1,3) = 3.920708
          test_data(6,1,3) = 0.604901
          test_data(7,1,3) = 0.800079

          test_data(1,2,3) = 2.225173
          test_data(2,2,3) = 0.531342
          test_data(6,2,3) = 1.354755
          test_data(7,2,3) =-1.55981 

          test_data(1,4,3) = -4.55588
          test_data(2,4,3) = -0.43447
          test_data(6,4,3) =  0.62401
          test_data(7,4,3) = -4.21117

          test_data(1,5,3) = 1.090653
          test_data(2,5,3) =-1.26043
          test_data(6,5,3) = 0.408451
          test_data(7,5,3) = 0.620509


          !test the gradients          
          call check_corners(timedev,test_data,detailled)

        end subroutine test_timedev

      
        subroutine check_corners(gradients_n, test_data, detailled)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: gradients_n
          real(rkind), dimension(nx,ny,ne), intent(in) :: test_data
          logical                         , intent(in) :: detailled


          integer(ikind) :: i,j,k
          logical        :: test_validated


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

        end subroutine check_corners


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
          
          logical :: loc

          loc = is_test_validated(
     $         gradients(i,j,k),
     $         test_data(i,j,k),detailled)
          if(detailled) print '(''test_grad('',3I2,''):'',L2)', 
     $         i,j,k, loc
          test_validated = test_validated.and.loc
          
        end subroutine check_data

      end program test_hedstrom_xy_corners

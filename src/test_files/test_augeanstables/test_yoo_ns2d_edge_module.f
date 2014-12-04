      module test_yoo_ns2d_edge_module

        use parameters_constant, only :
     $       vector_x

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq


        implicit none


        private
        public ::
     $       initialize_nodes,
     $       initialize_lodi_intermediate,
     $       get_test_data_for_lodi_outflow_x,
     $       get_test_data_for_lodi_inflow_x,
     $       get_test_data_for_lodi_outflow_timedevx,
     $       get_test_data_for_lodi_inflow_timedevx

        contains


        subroutine initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)

          implicit none

          type(pmodel_eq)                 , intent(in)  :: p_model
          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          real(rkind), dimension(nx)      , intent(out) :: x_map
          real(rkind), dimension(ny)      , intent(out) :: y_map
          real(rkind)                     , intent(out) :: dx
          real(rkind)                     , intent(out) :: dy

          integer, dimension(ne) :: var_type
          integer(ikind) :: i,j
          integer        :: k


          !fill the nodes (i in [1,3])x(j in [1,5])
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
          nodes(3,1,2) =  0.0d0
                    
          nodes(1,2,2) = -6.3d0
          nodes(2,2,2) =  7.98d0
          nodes(3,2,2) =  0.0d0
                    
          nodes(1,3,2) = -0.15d0
          nodes(2,3,2) = -6.213d0
          nodes(3,3,2) =  0.0d0
                    
          nodes(1,4,2) =  8.23d0
          nodes(2,4,2) =  3.012d0
          nodes(3,4,2) =  0.0d0
                    
          nodes(1,5,2) = -1.23d0
          nodes(2,5,2) =  7.8d0
          nodes(3,5,2) =  0.0d0

          !momentum-y----------
          nodes(1,1,3) =  2.01d0
          nodes(2,1,3) =  3.25d0
          nodes(3,1,3) =  6.2d0
                    
          nodes(1,2,3) =  7.135d0
          nodes(2,2,3) = -2.01d0
          nodes(3,2,3) =  3.06d0
                    
          nodes(1,3,3) =  9.46d0
          nodes(2,3,3) =  9.16d0
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
          var_type = p_model%get_var_type()

          do k=1, ne
             do j=1,5
                do i=1,2
                   if(var_type(k).eq.vector_x) then
                      nodes(6-i,j,k) = - nodes(i,j,k)
                   else
                      nodes(6-i,j,k) =   nodes(i,j,k)
                   end if
                end do
             end do
          end do


          !set to zero the nodes (i in [3])x(j in [1,5]) for
          !variables of type vector_x
          i=3
          do k=1,ne
             if(var_type(k).eq.vector_x) then
                do j=1,5
                   nodes(i,j,k)=0.0
                end do
             end if
          end do


          !initialize dx
          dx = 0.4d0

          !intiialize dy
          dy = 0.6d0

          !initialize the x_map
          do i=1,5
             x_map(i) = (i-1)*dx
          end do

          !initialize the y_map
          do i=1,5
             y_map(i) = (i-1)*dy
          end do
          
       end subroutine initialize_nodes


       subroutine initialize_lodi_intermediate(transverse_lodi, viscous_lodi)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: transverse_lodi
          real(rkind), dimension(nx,ny,ne), intent(out) :: viscous_lodi


          integer(ikind) :: i,j
          integer        :: k

          
          do k=1, ne
             do j=1, ny
                do i=1, nx
                   transverse_lodi(i,j,k) = 0.1d0
                   viscous_lodi(i,j,k) = 0.2d0
                end do
             end do
          end do

       end subroutine initialize_lodi_intermediate


               subroutine get_test_data_for_lodi_outflow_x(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          !L1
          test_data(1,5,1) =  446.77785d0
          test_data(1,4,1) = -1.2737843d0
          test_data(1,3,1) = 63.7416509d0
          test_data(1,2,1) = 81.1066497d0
          test_data(1,1,1) = -15.193722d0
                                       
          test_data(2,5,1) = -82.306674d0
          test_data(2,4,1) = -1.1667903d0
          test_data(2,3,1) = 69.9870196d0
          test_data(2,2,1) = -6.3174795d0
          test_data(2,1,1) = 10.5902500d0
                                       
          test_data(4,5,1) = -82.306674d0
          test_data(4,4,1) = -1.1667903d0
          test_data(4,3,1) = 69.9870196d0
          test_data(4,2,1) = -6.3174795d0
          test_data(4,1,1) = 10.5902500d0
                                       
          test_data(5,5,1) = 446.777854d0
          test_data(5,4,1) = -1.2737843d0
          test_data(5,3,1) = 63.7416509d0
          test_data(5,2,1) = 81.1066497d0
          test_data(5,1,1) = -15.193722d0
                                       
          !L2                          
          test_data(1,5,2) = 612.011517d0
          test_data(1,4,2) = 60.8978772d0
          test_data(1,3,2) = -1973.1928d0
          test_data(1,2,2) = -3957.7409d0
          test_data(1,1,2) = 3753.61415d0
                                       
          test_data(2,5,2) = 1648.38519d0
          test_data(2,4,2) = -187.99854d0
          test_data(2,3,2) = -501.81569d0
          test_data(2,2,2) = 29.5245331d0
          test_data(2,1,2) = 2106.98704d0
                                       
          test_data(4,5,2) = 1648.38519d0
          test_data(4,4,2) = -187.99854d0
          test_data(4,3,2) = -501.81569d0
          test_data(4,2,2) = 29.5245331d0
          test_data(4,1,2) = 2106.98704d0
                                       
          test_data(5,5,2) = 612.011517d0
          test_data(5,4,2) = 60.8978772d0
          test_data(5,3,2) = -1973.1928d0
          test_data(5,2,2) = -3957.7409d0
          test_data(5,1,2) = 3753.61415d0
                                       
          !L3                      
          test_data(1,5,3) = -3872.8478d0
          test_data(1,4,3) = 69.6020753d0
          test_data(1,3,3) = -2763.6543d0
          test_data(1,2,3) = 1973.37048d0
          test_data(1,1,3) = -6295.0801d0
                                       
          test_data(2,5,3) = 1367.715676d0
          test_data(2,4,3) = -299.003693d0
          test_data(2,3,3) = 497.1727008d0
          test_data(2,2,3) = -11.7202103d0
          test_data(2,1,3) = -15463.5656d0
                       
          test_data(4,5,3) = 98.40374684d0
          test_data(4,4,3) = -7.49140142d0
          test_data(4,3,3) = 63.94362383d0
          test_data(4,2,3) = -8.81893645d0
          test_data(4,1,3) = 417.8827844d0
                             
          test_data(5,5,3) = -8.8057865d0
          test_data(5,4,3) = 14.9460083d0
          test_data(5,3,3) = 19.5734597d0
          test_data(5,2,3) = 9.17385277d0
          test_data(5,1,3) = -5.4297627d0
                                       
          !L4              
          test_data(1,5,4) = -8.8057865d0
          test_data(1,4,4) = 14.9460083d0
          test_data(1,3,4) = 19.5734597d0
          test_data(1,2,4) = 9.17385277d0
          test_data(1,1,4) = -5.4297627d0
                                       
          test_data(2,5,4) = 98.40374684d0
          test_data(2,4,4) = -7.49140142d0
          test_data(2,3,4) = 63.94362383d0
          test_data(2,2,4) = -8.81893645d0
          test_data(2,1,4) = 417.8827844d0
                                       
          test_data(4,5,4) = 1367.715676d0
          test_data(4,4,4) = -299.003693d0
          test_data(4,3,4) = 497.1727008d0
          test_data(4,2,4) = -11.7202103d0
          test_data(4,1,4) = -15463.5656d0
                       
          test_data(5,5,4) = -3872.8478d0
          test_data(5,4,4) = 69.6020753d0
          test_data(5,3,4) = -2763.6543d0
          test_data(5,2,4) = 1973.37048d0
          test_data(5,1,4) = -6295.0801d0

        end subroutine get_test_data_for_lodi_outflow_x


        subroutine get_test_data_for_lodi_outflow_timedevx(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data


          !mass
          test_data(1,5,1) = 44.10934585d0
          test_data(1,4,1) = -15.2165587d0
          test_data(1,3,1) = 11.25514976d0
          test_data(1,2,1) = 79.09267117d0
          test_data(1,1,1) = -110.744439d0

          test_data(2,5,1) = -35.2692175d0
          test_data(2,4,1) = 341.6286828d0
          test_data(2,3,1) = 5.979210546d0
          test_data(2,2,1) = -33.0248749d0
          test_data(2,1,1) = 5.796465574d0
                                      
          test_data(4,5,1) = -35.2692175d0
          test_data(4,4,1) = 341.6286828d0
          test_data(4,3,1) = 5.979210546d0
          test_data(4,2,1) = -33.0248749d0
          test_data(4,1,1) = 5.796465574d0
                                      
          test_data(5,5,1) = 44.10934585d0
          test_data(5,4,1) = -15.2165587d0
          test_data(5,3,1) = 11.25514976d0
          test_data(5,2,1) = 79.09267117d0
          test_data(5,1,1) = -110.744439d0

                                 
          !momentum-x
          test_data(1,5,2) = -623.274248d0
          test_data(1,4,2) = -4.23812689d0 
          test_data(1,3,2) = -88.0603233d0
          test_data(1,2,2) = -254.873883d0
          test_data(1,1,2) = -1056.25085d0

          test_data(2,5,2) = 0.818692465d0
          test_data(2,4,2) = -13.9166020d0
          test_data(2,3,2) = 27.35376184d0
          test_data(2,2,2) = -32.5438124d0
          test_data(2,1,2) = -213.608694d0
                             
          test_data(4,5,2) = -0.81869246d0
          test_data(4,4,2) = 13.91660206d0  
          test_data(4,3,2) = -27.3537618d0
          test_data(4,2,2) = 32.54381247d0  
          test_data(4,1,2) = 213.6086941d0
                             
          test_data(5,5,2) = 623.274248d0
          test_data(5,4,2) = 4.23812689d0
          test_data(5,3,2) = 88.0603233d0
          test_data(5,2,2) = 254.873883d0
          test_data(5,1,2) = 1056.25085d0

                                 
          !momentum-y
          test_data(1,5,3) = 1256.48673d0
          test_data(1,4,3) = 7.01407643d0 
          test_data(1,3,3) = 448.268623d0
          test_data(1,2,3) = 372.943861d0
          test_data(1,1,3) = -61.835447d0
                             
          test_data(2,5,3) = 245.0656909d0
          test_data(2,4,3) = -85.0659158d0
          test_data(2,3,3) = -302.770573d0
          test_data(2,2,3) = 62.04892821d0
          test_data(2,1,3) = 7.667075478d0
                             
          test_data(4,5,3) = 245.0656909d0
          test_data(4,4,3) = -85.0659158d0
          test_data(4,3,3) = -302.770573d0
          test_data(4,2,3) = 62.04892821d0
          test_data(4,1,3) = 7.667075478d0
                             
          test_data(5,5,3) = 1256.48673d0 
          test_data(5,4,3) = 7.01407643d0
          test_data(5,3,3) = 448.268623d0 
          test_data(5,2,3) = 372.943861d0 
          test_data(5,1,3) = -61.835447d0

          !total energy
          test_data(1,5,4) = 23715.67526d0
          test_data(1,4,4) = -58.1465572d0
          test_data(1,3,4) = 11030.34145d0
          test_data(1,2,4) = -419.433190d0
          test_data(1,1,4) = 7872.233274d0
                             
          test_data(2,5,4) = -621.786095d0
          test_data(2,4,4) = 209.4959552d0
          test_data(2,3,4) = -1092.99628d0
          test_data(2,2,4) = -14.1760078d0
          test_data(2,1,4) = 9393.828732d0
                             
          test_data(4,5,4) = -621.786095d0
          test_data(4,4,4) = 209.4959552d0
          test_data(4,3,4) = -1092.99628d0
          test_data(4,2,4) = -14.1760078d0
          test_data(4,1,4) = 9393.828732d0
                             
          test_data(5,5,4) = 23715.67526d0
          test_data(5,4,4) = -58.1465572d0
          test_data(5,3,4) = 11030.34145d0
          test_data(5,2,4) = -419.433190d0
          test_data(5,1,4) = 7872.233274d0

        end subroutine get_test_data_for_lodi_outflow_timedevx


        subroutine get_test_data_for_lodi_inflow_x(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data


          !L1
          test_data(1,5,1) = 540.4347568d0
          test_data(1,4,1) = 0.562369759d0
          test_data(1,3,1) = 306.7471758d0
          test_data(1,2,1) = 24.23568908d0
          test_data(1,1,1) = 3.421977833d0
                             
          test_data(2,5,1) = 1.70950990d0
          test_data(2,4,1) = -0.1563713d0
          test_data(2,3,1) = 2.83158308d0
          test_data(2,2,1) = -0.7406068d0
          test_data(2,1,1) = 3.06662209d0
                             
          test_data(4,5,1) = 1.70950990d0
          test_data(4,4,1) = -0.1563713d0
          test_data(4,3,1) = 2.83158308d0
          test_data(4,2,1) = -0.7406068d0
          test_data(4,1,1) = 3.06662209d0
                             
          test_data(5,5,1) = 540.4347568d0
          test_data(5,4,1) = 0.562369759d0
          test_data(5,3,1) = 306.7471758d0
          test_data(5,2,1) = 24.23568908d0
          test_data(5,1,1) = 3.421977833d0

          !L2
          test_data(1,5,2) = 54.23088137d0
          test_data(1,4,2) = -12.8598311d0
          test_data(1,3,2) = 1217.924473d0
          test_data(1,2,2) = 30.54118812d0
          test_data(1,1,2) = -43.3197689d0
                                      
          test_data(2,5,2) =  25.57906023d0
          test_data(2,4,2) = -21.11020646d0
          test_data(2,3,2) =   9.36048316d0
          test_data(2,2,2) = -57.19780745d0
          test_data(2,1,2) =  478.8211866d0 
                                         
          test_data(4,5,2) =  25.57906023d0
          test_data(4,4,2) = -21.11020646d0
          test_data(4,3,2) =   9.36048316d0
          test_data(4,2,2) = -57.19780745d0
          test_data(4,1,2) =  478.8211866d0
                                      
          test_data(5,5,2) = 54.23088137d0
          test_data(5,4,2) = -12.8598311d0
          test_data(5,3,2) = 1217.924473d0
          test_data(5,2,2) = 30.54118812d0
          test_data(5,1,2) = -43.3197689d0

          !L3
          test_data(1,5,3) = -3872.8478d0
          test_data(1,4,3) = 69.6020753d0 
          test_data(1,3,3) = -2763.6543d0
          test_data(1,2,3) = 1973.37048d0
          test_data(1,1,3) = -6295.0801d0
                                         
          test_data(2,5,3) = 1367.715676d0 
          test_data(2,4,3) = -299.003693d0 
          test_data(2,3,3) = 497.1727008d0 
          test_data(2,2,3) = -11.7202103d0
          test_data(2,1,3) = -15463.5656d0   

          test_data(4,5,3) =-11.60674622d0
          test_data(4,4,3) =-4.945547454d0
          test_data(4,3,3) = 1.687066552d0
          test_data(4,2,3) =-7.054968428d0
          test_data(4,1,3) =-33.34939081d0
                                          
          test_data(5,5,3) =  19.04542009d0
          test_data(5,4,3) = -7.119291704d0
          test_data(5,3,3) = -1.242953254d0
          test_data(5,2,3) =  16.10798298d0
          test_data(5,1,3) =  6.287316017d0


          !L4
          test_data(1,5,4) =-27.0795799d0
          test_data(1,4,4) =0.142473002d0
          test_data(1,3,4) =-6.13425760d0
          test_data(1,2,4) =-23.2670170d0
          test_data(1,1,4) =-13.4246405d0
                                         
          test_data(2,5,4) =4.643253785d0
          test_data(2,4,4) =-2.04939360d0
          test_data(2,3,4) =-8.66793344d0
          test_data(2,2,4) =-0.09566610d0
          test_data(2,1,4) =26.39325625d0
                                           
          test_data(4,5,4) =1367.715676d0
          test_data(4,4,4) =-299.003693d0
          test_data(4,3,4) =497.1727008d0
          test_data(4,2,4) =-11.7202103d0
          test_data(4,1,4) =-15463.5656d0 
                                         
          test_data(5,5,4) =-3872.8478d0 
          test_data(5,4,4) =69.6020753d0 
          test_data(5,3,4) =-2763.6543d0 
          test_data(5,2,4) =1973.37048d0 
          test_data(5,1,4) =-6295.0801d0 

        end subroutine get_test_data_for_lodi_inflow_x


        subroutine get_test_data_for_lodi_inflow_timedevx(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          !mass
          test_data(1,5,1) = 62.9278839d0
          test_data(1,4,1) = -3.2465581d0
          test_data(1,3,1) = 0.56177567d0
          test_data(1,2,1) = -26.811349d0
          test_data(1,1,1) = 586.903034d0
                                 
          test_data(2,5,1) = -10.5411492d0
          test_data(2,4,1) = 171.8291839d0
          test_data(2,3,1) = -6.85357331d0
          test_data(2,2,1) = 108.2349372d0
          test_data(2,1,1) = 7.748555568d0
                   
          test_data(4,5,1) = -10.42081796d0
          test_data(4,4,1) =  173.2788844d0
          test_data(4,3,1) = -6.993488833d0
          test_data(4,2,1) =  114.2030122d0
          test_data(4,1,1) =  7.780526165d0                                     
                                     
          test_data(5,5,1) =  62.16233605d0
          test_data(5,4,1) =  -2.71104871d0
          test_data(5,3,1) =  0.553547195d0
          test_data(5,2,1) = -27.33626217d0
          test_data(5,1,1) =  585.0940041d0
               
                    
          !momentum-x        
          test_data(1,5,2) = -737.343573d0
          test_data(1,4,2) = 10.19423164d0
          test_data(1,3,2) = -80.3408007d0
          test_data(1,2,2) = 303.7707912d0
          test_data(1,1,2) = -2888.13781d0
                             
          test_data(2,5,2) = 60.10132975d0
          test_data(2,4,2) = -82.2078604d0
          test_data(2,3,2) = 51.03983127d0
          test_data(2,2,2) = 92.82001573d0
          test_data(2,1,2) = -191.655149d0

          test_data(4,5,2) = -61.35083065d0
          test_data(4,4,2) =  80.19916456d0
          test_data(4,3,2) = -50.38188365d0
          test_data(4,2,2) = -102.9148968d0
          test_data(4,1,2) =  190.4232389d0
                            
          test_data(5,5,2) =  736.8372929d0     
          test_data(5,4,2) = -12.10713806d0     
          test_data(5,3,2) =  80.47729343d0     
          test_data(5,2,2) = -303.3118948d0
          test_data(5,1,2) =  2887.605745d0


          !momentum-y
          test_data(1,5,3) = 1811.93718d0
          test_data(1,4,3) = -5.5936922d0
          test_data(1,3,3) = -47.445772d0
          test_data(1,2,3) = -188.49864d0
          test_data(1,1,3) = 505.031667d0
                             
          test_data(2,5,3) = -21.4681830d0
          test_data(2,4,3) = -46.1434760d0
          test_data(2,3,3) = -26.6929531d0
          test_data(2,2,3) = -18.9275509d0
          test_data(2,1,3) = 21.56107056d0
                             
          test_data(4,5,3) = -21.29336841d0
          test_data(4,4,3) = -46.54307295d0
          test_data(4,3,3) = -26.97775891d0
          test_data(4,2,3) = -20.32241503d0
          test_data(4,1,3) =  21.66293765d0
                             
          test_data(5,5,3) =  1788.579165d0
          test_data(5,4,3) = -5.459499861d0
          test_data(5,3,3) = -47.78421365d0
          test_data(5,2,3) = -191.6196857d0
          test_data(5,1,3) =  503.4507328d0


          !total energy
          test_data(1,5,4) = 32263.07675d0
          test_data(1,4,4) = -42.2159859d0
          test_data(1,3,4) = -297.028576d0
          test_data(1,2,4) = -3334.73191d0
          test_data(1,1,4) = 10512.30477d0
                                    
          test_data(2,5,4) = -894.371848d0
          test_data(2,4,4) = 187.4252157d0
          test_data(2,3,4) = -470.451434d0
          test_data(2,2,4) = 49.86206488d0
          test_data(2,1,4) = 9834.749701d0
                             
          test_data(4,5,4) =  -879.632557d0
          test_data(4,4,4) =  190.3199834d0
          test_data(4,3,4) = -477.4657912d0
          test_data(4,2,4) =  62.04236643d0
          test_data(4,1,4) =  9888.517688d0

          test_data(5,5,4) =  31883.50279d0
          test_data(5,4,4) = -35.15171993d0
          test_data(5,3,4) = -307.5663991d0
          test_data(5,2,4) =  -3363.89861d0
          test_data(5,1,4) =  10501.67974d0

        end subroutine get_test_data_for_lodi_inflow_timedevx

      end module test_yoo_ns2d_edge_module

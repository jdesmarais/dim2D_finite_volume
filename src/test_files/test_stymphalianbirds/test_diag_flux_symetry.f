      program test_diag_flux_symetry

        use parameters_constant, only :
     $     vector_x,
     $     vector_y

        use parameters_input, only :
     $     ne

        use parameters_kind, only : 
     $     ikind,
     $     rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_n1_oneside_R1_class, only :
     $       sd_operators_n1_oneside_R1

        use sd_operators_n2_oneside_R1_class, only :
     $       sd_operators_n2_oneside_R1

        implicit none

        real(rkind)                      :: dx
        real(rkind)                      :: dy
        real(rkind)                      :: dn
        real(rkind), dimension(5,5,ne)   :: nodes
        real(rkind), dimension(5,5,ne)   :: nodes_sym
        real(rkind), dimension(5,5,ne)   :: nodes_n
        real(rkind), dimension(5,5,ne)   :: nodes_n_sym

        type(sd_operators_n1_oneside_R1) :: sd_n1_R1
        type(sd_operators_n2_oneside_R1) :: sd_n2_R1

        type(pmodel_eq)                  :: p_model

        real(rkind), dimension(ne)       :: flux_diag1
        real(rkind), dimension(ne)       :: flux_diag2
        real(rkind), dimension(ne)       :: flux_diag1_sym
        real(rkind), dimension(ne)       :: flux_diag2_sym

        integer(ikind) :: i,j,k
        integer(ikind) :: i_c,j_c


        !initialize the nodes
        call initialize_nodes(dx,dy,nodes)
        dn = Sqrt(dx**2+dy**2)

        i_c = 3
        j_c = 3
        
        !make the symmetry comparison for NE_edge and SE_edge
        !compute the fluxes for the SE_edge
        do j=1,size(nodes,2)
           do i=1, size(nodes,1)
              nodes_n(i,j,:) = p_model%compute_xy_to_n_var(nodes(i,j,:))
           end do
        end do

        flux_diag1_sym = p_model%compute_flux_y_oneside(
     $       nodes_n,dn,dn, i_c  , j_c  , sd_n1_R1)
        flux_diag2_sym = p_model%compute_flux_y_oneside(
     $       nodes_n,dn,dn, i_c+1, j_c+1, sd_n1_R1)

        !create the symmetry for the nodes
        !and compute the fluxes for the NE_edge
        nodes_sym = make_symetry_y(p_model,nodes)

        do j=1,size(nodes,2)
           do i=1, size(nodes,1)
              nodes_n_sym(i,j,:) = p_model%compute_xy_to_n_var(nodes_sym(i,j,:))
           end do
        end do

        flux_diag1 = p_model%compute_flux_x_oneside(
     $       nodes_n_sym,dn,dn, i_c  , j_c  , sd_n2_R1)
        flux_diag2 = p_model%compute_flux_x_oneside(
     $       nodes_n_sym,dn,dn, i_c+1, j_c-1, sd_n2_R1)

        !compare the fluxes
        !print *, 'flux_diag1: '    , flux_diag1
        !print *, 'flux_diag1_sym: ', flux_diag1_sym
        !print '()'
        !print *, 'flux_diag2: '    , flux_diag2
        !print *, 'flux_diag2_sym: ', flux_diag2_sym
        

        print '(A20,'' | '',A20)', 'flux_diag1', 'flux_diag1_sym'
        print '()'
        do k=1, ne
           print '(F20.15,'' | '',F20.15)', flux_diag1(k), flux_diag1_sym(k)
           print '(A20,'' | '',A20)', ' ', ' '
        end do
        print '()'

        print '(A20,'' | '',A20)', 'flux_diag2', 'flux_diag2_sym'
        print '()'
        do k=1, ne
           print '(F20.15,'' | '',F20.15)', flux_diag2(k), flux_diag2_sym(k)
           print '(A20,'' | '',A20)', ' ', ' '
        end do
        print '()'
        

        contains

        subroutine initialize_nodes(dx,dy,nodes)

          implicit none

          real(rkind), dimension(5,5,ne), intent(out) :: nodes
          real(rkind)                   , intent(out) :: dx
          real(rkind)                   , intent(out) :: dy

          dx = 0.1
          dy = 0.2

          !   _ _ _ _ _
          !  |_|_|1|_|1|
          !  |_|1|_|1|_|
          !  |1|_|2|_|1|
          !  |_|1|_|1|_|
          !  |1|_|1|_|_|
          !------------------
          nodes(1,1,1) = 1.3d0
          nodes(3,1,1) = 1.4d0
          nodes(2,2,1) = 1.5d0
          nodes(4,2,1) = 1.2d0
          nodes(1,3,1) = 1.25d0
          nodes(3,3,1) = 1.45d0
          nodes(5,3,1) = 1.1d0
          nodes(2,4,1) = 1.2d0
          nodes(4,4,1) = 1.35d0
          nodes(3,5,1) = 1.45d0
          nodes(5,5,1) = 1.55d0
          
          do j=1,5
             do i=1,4
                nodes(i,j,2) = nodes(i,j,1)*0.1
                nodes(i,j,3) = nodes(i,j,1)*-0.02
             end do
          end do

          j=3
          do i=1,5
             nodes(i,j,3) = 0.0
          end do
          
          nodes(1,1,4) = 5.0d0
          nodes(3,1,4) = 4.98d0
          nodes(2,2,4) = 4.8d0
          nodes(4,2,4) = 4.93d0
          nodes(1,3,4) = 4.95d0
          nodes(3,3,4) = 4.88d0
          nodes(5,3,4) = 4.91d0
          nodes(2,4,4) = 4.93d0
          nodes(4,4,4) = 5.0d0
          nodes(3,5,4) = 4.93d0
          nodes(5,5,4) = 4.87d0

        end subroutine initialize_nodes


        function make_symetry_x(p_model,nodes) result(tmp_nodes)
        
          implicit none

          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind)    , dimension(5,5,ne), intent(inout) :: nodes
          real(rkind)    , dimension(5,5,ne)                :: tmp_nodes

          integer       , dimension(ne)     :: var_type
          integer(ikind)                    :: i,j
          integer                           :: k
          real(rkind)                       :: coeff

          var_type = p_model%get_var_type()


          do k=1,ne

             coeff = 1.0d0

             if(var_type(k).eq.vector_x) then
                coeff = -1.0d0
             end if

             do j=1, size(tmp_nodes,2)
                do i=1, size(tmp_nodes,1)
                
                   tmp_nodes(i,j,k) = coeff*nodes(5-i+1,j,k)
             
                end do
             end do
          end do

          nodes = tmp_nodes

        end function make_symetry_x


        function make_symetry_y(p_model,nodes) result(tmp_nodes)
        
          implicit none

          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind)    , dimension(5,5,ne), intent(inout) :: nodes
          real(rkind)    , dimension(5,5,ne)                :: tmp_nodes


          integer       , dimension(ne)     :: var_type
          integer(ikind)                    :: i,j
          integer                           :: k
          real(rkind)                       :: coeff


          var_type = p_model%get_var_type()


          do k=1,ne

             coeff = 1.0d0

             if(var_type(k).eq.vector_y) then
                coeff = -1.0d0
             end if

             do j=1, size(tmp_nodes,2)
                do i=1, size(tmp_nodes,1)
                
                   tmp_nodes(i,j,k) = coeff*nodes(i,5-j+1,k)
             
                end do
             end do
          end do          

        end function make_symetry_y

      end program test_diag_flux_symetry

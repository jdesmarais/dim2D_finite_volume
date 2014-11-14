      module bf_newgrdpt_class

        use bf_layer_newgrdpt_procedure_module, only :
     $       no_gradient_type,
     $       gradient_L0_type,
     $       gradient_R0_type

        use parameters_constant, only :
     $       left, right, x_direction, y_direction

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use interface_primary, only :
     $       gradient_y_proc

        implicit none

        private
        public :: bf_newgrdpt

        !object encapsulating the subroutines for the computation of
        !the newgrid points
        type :: bf_newgrdpt

          contains

          procedure, nopass :: compute_newgrdpt_x
          procedure, nopass :: get_interpolation_coeff_1D
          procedure, nopass :: interpolate_1D
          procedure, nopass :: compute_NewtonCotes_integration

        end type bf_newgrdpt


        contains


        !p_model       : physical model
        !dt            : time step
        !              
        !bf_align0     : alignment of the buffer layer at t=t-dt
        !bf_x_map0     : x-coordinates of the buffer layer at t=t-dt
        !bf_y_map0     : y-coordinates of the buffer layer at t=t-dt
        !bf_nodes0     : nodes of the buffer layer at t=t-dt
        !
        !bf_align1     : alignment of the buffer layer at t=t
        !bf_x_map1     : x-coordinates of the buffer layer at t=t
        !bf_y_map1     : y-coordinates of the buffer layer at t=t
        !bf_nodes1     : nodes of the buffer layer at t=t
        !              
        !i1            : x-index identifying the new grdpt at t=t
        !j1            : y-index identifying the new grdpt at t=t
        !              
        !side_x        : logical identifying the type of boundary
        !                (E or W)
        !gradient_type : integer identifying the type of gradient
        !                to apply for the transverse terms
        !-----------------------------------------------------------------
        subroutine compute_newgrdpt_x(
     $       p_model, t, dt,
     $       bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $       bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $       i1,j1, side_x, gradient_y)

          implicit none

          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)                        , intent(in)    :: t
          real(rkind)                        , intent(in)    :: dt
          integer(ikind), dimension(2,2)     , intent(in)    :: bf_align0
          real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map0
          real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map0
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes0
          integer(ikind), dimension(2,2)     , intent(in)    :: bf_align1
          real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map1
          real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map1
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: bf_nodes1
          integer(ikind)                     , intent(in)    :: i1
          integer(ikind)                     , intent(in)    :: j1
          logical                            , intent(in)    :: side_x
          procedure(gradient_y_proc)                         :: gradient_y

          integer                       :: k
          
          integer                       :: dir, dir2
          integer(ikind)                :: i_eigen
          real(rkind), dimension(ne)    :: eigenvalues_x
          real(rkind), dimension(ne,ne) :: left_eigenM
          real(rkind), dimension(ne,ne) :: right_eigenM
          real(rkind), dimension(ne)    :: n_amp0
          real(rkind), dimension(ne)    :: t_amp0
          real(rkind), dimension(ne)    :: t_map1
          real(rkind), dimension(ne)    :: amp
          real(rkind), dimension(ne)    :: char_amp

          real(rkind)                   :: x0,x1
          integer(ikind)                :: i0_inter1, i0_inter2, j0_inter
          integer(ikind)                :: i1_inter1, i1_inter2, j1_inter
          

          !0) determine the direction
          dir  = x_direction
          dir2 = y_direction


          !1) determine the x-coordinate of the new grid point computed
          x1   = bf_x_map1(i1)


          !2) determine where the eigenvalues are evaluated
          !   and which indices are needed for the interpolation
          !   of the grid points
          if(side_x.eqv.right) then

             !x-index for the evaluation of the eigenvalues
             i_eigen = i1-1
             
             !indices for the interpolation of the data at t
             i0_inter1 = bf_align1(dir,1) -bf_align0(dir,1) + i_eigen-1
             i0_inter2 = bf_align1(dir,1) -bf_align0(dir,1) + i_eigen
             j0_inter  = bf_align1(dir2,1)-bf_align0(dir2,1)+ j1
             
             !indices for the interpolation of the data at t+dt
             i1_inter1 = i_eigen-1
             i1_inter2 = i_eigen
             j1_inter  = j1

          else
             
             !x-index for the evaluation of the eigenvalues
             i_eigen = i1+1

             !indices for the interpolation of the data at t
             i0_inter1 = bf_align1(dir,1) -bf_align0(dir,1) +i_eigen
             i0_inter2 = bf_align1(dir,1) -bf_align0(dir,1) +i_eigen+1
             j0_inter  = bf_align1(dir2,1)-bf_align0(dir2,1)+j1
             
             !indices for the interpolation of the data at t+dt
             i1_inter1 = i_eigen
             i1_inter2 = i_eigen+1
             j1_inter  = j1
             
          end if

          dy = bf_y_map0(2)-bf_y_map0(1)
          dy = bf_y_map1(2)-bf_y_map1(1)


          !3) create the interpolation coefficients for the data at t

          !3.1) create the interpolation coefficients for the nodes
          x_map_inter(1) = bf_x_map0(i0_inter1)
          x_map_inter(2) = bf_x_map0(i0_inter2)
          
          nodes_inter(1,:) = bf_nodes0(i0_inter1,j0_inter,:)
          nodes_inter(2,:) = bf_nodes0(i0_inter1,j0_inter,:)

          inter_nodes0 = get_interpolation_coeff_1D(x_map_inter,nodes_inter)


          !3.2) create the interpolation coefficients for the transverse terms
          nodes_inter(1,:) = MATMUL(
     $         p_model%compute_y_gradient(bf_nodes0,i0_inter1,j0_inter,gradient_y,dy),
     $         p_model%compute_x_transM(bf_nodes0(i0_inter1,j0_inter,:)))

          nodes_inter(2,:) = MATMUL(
     $         p_model%compute_y_gradient(bf_nodes0,i0_inter2,j0_inter,gradient_y,dy),
     $         p_model%compute_x_transM(bf_nodes0(i0_inter2,j0_inter,:)))

          inter_trans0 = get_interpolation_coeff_1D(x_map_inter,nodes_inter)


          !4) create the interpolation coefficients for the data at t+dt
          x_map_inter(1) = bf_x_map1(i1_inter1)
          x_map_inter(2) = bf_x_map1(i1_inter2)

          nodes_inter(1,:) = MATMUL(
     $         p_model%compute_y_gradient(bf_nodes1,i1_inter1,j1_inter,gradient_y,dy),
     $         p_model%compute_x_transM(bf_nodes1(i1_inter1,j1_inter,:)))

          nodes_inter(2,:) = MATMUL(
     $         p_model%compute_y_gradient(bf_nodes1,i1_inter2,j1_inter,gradient_y,dy),
     $         p_model%compute_x_transM(bf_nodes1(i1_inter2,j1_inter,:)))

          inter_trans1 = get_interpolation_coeff_1D(x_map_inter,nodes_inter)

          t_amp1 = interpolate_1D(x1,inter_trans1)


          !5) evaluate the eigenvalues at t+dt
          eigenvalues_x = p_model%compute_x_eigenvalues(bf_nodes1(i_eigen,j1,:))


          !6) determine the left eigenvector corresponding to the eigenvalue
          left_eigenM = p_model%compute_x_lefteigenvector(bf_nodes1(i_eigen,j1,:))
             

          !4) determine the characteristic amplitude
          do k=1,ne


             !4.1) determine the position where the characteristic
             !     amplitude should be estimated
             x0 = x1 - eigenvalues_x(k)*dt


             !4.2) determine the normal and transverse contributions of
             !     the hyperbolic terms to the characteristic amplitude
             if(side_x.eq.right) then

                if(eigenvalues_x(k).ge.0) then
                   
                   n_amp0 = interpolate_1D(x0,inter_nodes0)
                   t_amp0 = interpolate_1D(x0,inter_trans0)
                   
                   amp =
     $                  n_map0 -
     $                  compute_NewtonCotes_integration(t_amp0, t_amp1, dt)

                else

                   n_amp0 = p_model%get_far_field(t,bf_x_map1(i1),bf_y_map1(j1))
                   t_amp0 = [0.0d0,0.0d0,0.0d0]

                end if

             else
                
                if(eigenvalues_x(k).gt.0) then

                   n_amp0 = p_model%get_far_field(t,x0,bf_y_map0(j0_inter))
                   t_amp0 = [0.0d0,0.0d0,0.0d0]

                else

                   n_amp0 = interpolate_1D(x0,inter_nodes0)
                   t_amp0 = interpolate_1D(x0,inter_trans0)                   

                end if
             end if


             !4.3) combine the information on the nodes at t-dt and the approximation
             !     of the integration of the transverse terms from t-dt to t
             amp =
     $            n_amp0 -
     $            compute_NewtonCotes_integration(t_amp0, t_amp1, dt)

             
             !4.4) compute the scalar product of the left eigenvector corresponding
             !     to the eigenvalue with the characteristic amplitude
             char_amp(k) = DOT_PRODUCT(char_amp,left_eigenM(:,k))
             
          end do


          !5) determine the right eigenmatrix
          right_eigenM = p_model%compute_x_righteigenvector(bf_nodes1(i_eigen,j1,:))


          !6) determine the new grid point
          bf_nodes1(i1,j1,:) = MATMUL(char_amp,right_eigenM)


        end subroutine compute_newgrdpt_x


        !get the interpolation coefficients for a 1st order polynomial
        !fit
        function get_interpolation_coeff_1D(x_map,nodes)
     $     result(inter_coeff)
        
          implicit none

          real(rkind), dimension(2)   , intent(in) :: x_map
          real(rkind), dimension(2,ne), intent(in) :: nodes
          real(rkind), dimension(2,ne)             :: inter_coeff

          do k=1, ne

             inter_coeff(1,k) = (nodes(2,k) - nodes(1,k))/(x_map(2)-x_map(1))
             inter_coeff(2,k) = nodes(1,k) - inter_coeff(1,k)*x_map(1)

          end do

        end function get_interpolation_coeff_1D


        !interpolate in 1D using a 1st order approximation
        function interpolate_1D(
     $     x,
     $     inter_coeff)
     $     result(nodes_inter)

          implicit none

          real(rkind)                 , intent(in) :: x
          real(rkind), dimension(2,ne), intent(in) :: inter_coeff
          real(rkind), dimension(ne)               :: nodes_inter

          integer :: k

          do k=1,ne
             nodes_inter(k) = x*inter_coeff(1,k) + inter_coeff(2,k)
          end do

        end function interpolate_1D


        !integrate a function between two points using Newton-Cotes
        !approximation
        function compute_NewtonCotes_integration(data0,data1,dt)
     $     result(data_integrated)

          implicit none

          real(rkind), dimension(ne), intent(in) :: data0
          real(rkind), dimension(ne), intent(in) :: data1
          real(rkind)               , intent(in) :: dt
          real(rkind), dimension(ne)             :: data_integrated

          integer :: k

          do k=1,ne

             data_integrated(k) = (data0(k) + 0.5d0*(data1(k)-data0(k)))*dt

          end do

        end function compute_NewtonCotes_integration

      end module bf_newgrdpt_class

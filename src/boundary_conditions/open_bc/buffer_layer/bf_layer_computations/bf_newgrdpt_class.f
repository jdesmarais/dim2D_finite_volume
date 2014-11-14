      !> @file
      !> main subroutines to compute the new grid points after extension
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the main subroutines computing the new grid
      !> points after the computational domain extension
      !
      !> @date
      !> 14_11_2014 - initial version         - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_newgrdpt_class

        use bf_layer_newgrdpt_procedure_module, only :
     $       no_gradient_type,
     $       gradient_L0_type,
     $       gradient_R0_type

        use n_coords_module, only :
     $       get_x_coord,
     $       get_y_coord,
     $       get_n1_coord,
     $       get_n2_coord

        use parameters_constant, only :
     $       left, right, x_direction, y_direction

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        implicit none

        private
        public :: bf_newgrdpt


        !>@class bf_compute
        !> class encapsulating the main subroutines to compute the
        !> new gridpoints after the computational domain extension
        !
        !>@param compute_newgrdpt_x
        !> compute a new grid point obtained by extension in the
        !> x-direction
        !
        !>@param get_interpolation_coeff_1D
        !> get the interpolation coefficients for a 1st order polynomial
        !> function
        !
        !>@param interpolate_1D
        !> make use of the interpolation coefficients to interpolate the
        !> grid points data
        !
        !>@param compute_NewtonCotes_integration
        !> perform Newton-Cotes integration to integrate the data between
        !> two grid points
        !---------------------------------------------------------------
        type :: bf_newgrdpt

          contains

          procedure, nopass :: compute_newgrdpt_x
          procedure, nopass :: compute_newgrdpt_y
          procedure, nopass :: compute_newgrdpt_xy

          procedure, nopass :: get_interpolation_coeff_1D
          procedure, nopass :: interpolate_1D

          procedure, nopass :: get_interpolation_coeff_2D
          procedure, nopass :: interpolate_2D

          procedure, nopass :: compute_NewtonCotes_integration

        end type bf_newgrdpt


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid point obtained by extension of the
        !> computational domain in the x-direction
        !
        !> @date
        !> 14_11_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> physical model
        !
        !>@param dt
        !> time step
        !              
        !>@param bf_align0
        !> alignment of the buffer layer at t=t-dt
        !
        !>@param bf_x_map0
        !> x-coordinates of the buffer layer at t=t-dt
        !
        !>@param bf_y_map0
        !> y-coordinates of the buffer layer at t=t-dt
        !
        !>@param bf_nodes0
        !> nodes of the buffer layer at t=t-dt
        !
        !>@param bf_align1
        !> alignment of the buffer layer at t=t
        !
        !>@param bf_x_map1
        !> x-coordinates of the buffer layer at t=t
        !
        !>@param bf_y_map1
        !> y-coordinates of the buffer layer at t=t
        !
        !>@param bf_nodes1
        !> nodes of the buffer layer at t=t
        !              
        !>@param i1
        !> x-index identifying the new grdpt at t=t
        !
        !>@param j1
        !> y-index identifying the new grdpt at t=t
        !              
        !>@param side_x
        !> logical identifying the type of boundary (E or W)
        !
        !>@param gradient_y
        !> gradient procedure applied to compute the
        !> the transverse terms
        !--------------------------------------------------------------
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
          real(rkind), dimension(ne)    :: t_amp1
          real(rkind), dimension(ne)    :: amp
          real(rkind), dimension(ne)    :: char_amp

          real(rkind)                   :: dy
          real(rkind)                   :: x0,x1,y0,y1
          integer(ikind)                :: i0_inter1, i0_inter2, j0_inter
          integer(ikind)                :: i1_inter1, i1_inter2, j1_inter
          real(rkind), dimension(2)     :: x_map_inter
          real(rkind), dimension(2,ne)  :: nodes_inter
          real(rkind), dimension(2,ne)  :: inter_nodes0
          real(rkind), dimension(2,ne)  :: inter_trans0
          real(rkind), dimension(2,ne)  :: inter_trans1

          !0) determine the direction
          dir  = x_direction
          dir2 = y_direction


          !1) determine the x-coordinate of the new grid point computed
          x1   = bf_x_map1(i1)
          y1   = bf_y_map1(j1)


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
          nodes_inter(2,:) = bf_nodes0(i0_inter2,j0_inter,:)

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
             

          !7) determine the characteristic amplitude
          y0 = y1

          do k=1,ne


             !7.1) determine the position where the characteristic
             !     amplitude should be estimated
             x0 = x1 - eigenvalues_x(k)*dt


             !7.2) determine the normal and transverse contributions of
             !     the hyperbolic terms to the characteristic amplitude
             if(side_x.eq.right) then

                if(eigenvalues_x(k).ge.0) then
                   
                   n_amp0 = interpolate_1D(x0,inter_nodes0)
                   t_amp0 = interpolate_1D(x0,inter_trans0)
                   
                else

                   n_amp0 = p_model%get_far_field(t,x0,y0)
                   t_amp0 = [0.0d0,0.0d0,0.0d0]

                end if

             else
                
                if(eigenvalues_x(k).gt.0) then

                   n_amp0 = p_model%get_far_field(t,x0,y0)
                   t_amp0 = [0.0d0,0.0d0,0.0d0]

                else

                   n_amp0 = interpolate_1D(x0,inter_nodes0)
                   t_amp0 = interpolate_1D(x0,inter_trans0)                   

                end if
             end if


             !7.3) combine the information on the nodes at t-dt and the approximation
             !     of the integration of the transverse terms from t-dt to t
             amp =
     $            n_amp0 -
     $            compute_NewtonCotes_integration(t_amp0, t_amp1, dt)

             
             !7.4) compute the scalar product of the left eigenvector corresponding
             !     to the eigenvalue with the characteristic amplitude
             char_amp(k) = DOT_PRODUCT(amp,left_eigenM(:,k))
             
          end do


          !8) determine the right eigenmatrix
          right_eigenM = p_model%compute_x_righteigenvector(bf_nodes1(i_eigen,j1,:))


          !9) determine the new grid point
          bf_nodes1(i1,j1,:) = MATMUL(char_amp,right_eigenM)


        end subroutine compute_newgrdpt_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid point obtained by extension of the
        !> computational domain in the x-direction
        !
        !> @date
        !> 14_11_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> physical model
        !
        !>@param dt
        !> time step
        !              
        !>@param bf_align0
        !> alignment of the buffer layer at t=t-dt
        !
        !>@param bf_x_map0
        !> x-coordinates of the buffer layer at t=t-dt
        !
        !>@param bf_y_map0
        !> y-coordinates of the buffer layer at t=t-dt
        !
        !>@param bf_nodes0
        !> nodes of the buffer layer at t=t-dt
        !
        !>@param bf_align1
        !> alignment of the buffer layer at t=t
        !
        !>@param bf_x_map1
        !> x-coordinates of the buffer layer at t=t
        !
        !>@param bf_y_map1
        !> y-coordinates of the buffer layer at t=t
        !
        !>@param bf_nodes1
        !> nodes of the buffer layer at t=t
        !              
        !>@param i1
        !> x-index identifying the new grdpt at t=t
        !
        !>@param j1
        !> y-index identifying the new grdpt at t=t
        !              
        !>@param side_x
        !> logical identifying the type of boundary (E or W)
        !
        !>@param gradient_y
        !> gradient procedure applied to compute the
        !> the transverse terms
        !--------------------------------------------------------------
        subroutine compute_newgrdpt_y(
     $       p_model, t, dt,
     $       bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $       bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $       i1,j1, side_y, gradient_x)

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
          logical                            , intent(in)    :: side_y
          procedure(gradient_x_proc)                         :: gradient_x

          integer                       :: k
          
          integer                       :: dir, dir2
          integer(ikind)                :: j_eigen
          real(rkind), dimension(ne)    :: eigenvalues_y
          real(rkind), dimension(ne,ne) :: left_eigenM
          real(rkind), dimension(ne,ne) :: right_eigenM
          real(rkind), dimension(ne)    :: n_amp0
          real(rkind), dimension(ne)    :: t_amp0
          real(rkind), dimension(ne)    :: t_amp1
          real(rkind), dimension(ne)    :: amp
          real(rkind), dimension(ne)    :: char_amp

          real(rkind)                   :: dx
          real(rkind)                   :: x0,x1,y0,y1
          integer(ikind)                :: j0_inter1, j0_inter2, i0_inter
          integer(ikind)                :: j1_inter1, j1_inter2, i1_inter
          real(rkind), dimension(2)     :: y_map_inter
          real(rkind), dimension(2,ne)  :: nodes_inter
          real(rkind), dimension(2,ne)  :: inter_nodes0
          real(rkind), dimension(2,ne)  :: inter_trans0
          real(rkind), dimension(2,ne)  :: inter_trans1

          !0) determine the direction
          dir  = y_direction
          dir2 = x_direction


          !1) determine the x-coordinate of the new grid point computed
          x1   = bf_x_map1(i1)
          y1   = bf_y_map1(j1)


          !2) determine where the eigenvalues are evaluated
          !   and which indices are needed for the interpolation
          !   of the grid points
          if(side_y.eqv.right) then

             !x-index for the evaluation of the eigenvalues
             j_eigen = j1-1
             
             !indices for the interpolation of the data at t
             j0_inter1 = bf_align1(dir,1) -bf_align0(dir,1) + j_eigen-1
             j0_inter2 = bf_align1(dir,1) -bf_align0(dir,1) + j_eigen
             i0_inter  = bf_align1(dir2,1)-bf_align0(dir2,1)+ i1
             
             !indices for the interpolation of the data at t+dt
             j1_inter1 = j_eigen-1
             j1_inter2 = j_eigen
             i1_inter  = i1

          else
             
             !x-index for the evaluation of the eigenvalues
             j_eigen = j1+1

             !indices for the interpolation of the data at t
             j0_inter1 = bf_align1(dir,1) -bf_align0(dir,1) +j_eigen
             j0_inter2 = bf_align1(dir,1) -bf_align0(dir,1) +j_eigen+1
             i0_inter  = bf_align1(dir2,1)-bf_align0(dir2,1)+i1
             
             !indices for the interpolation of the data at t+dt
             j1_inter1 = j_eigen
             j1_inter2 = j_eigen+1
             i1_inter  = i1
             
          end if

          dx = bf_x_map0(2)-bf_x_map0(1)
          dx = bf_x_map1(2)-bf_x_map1(1)


          !3) create the interpolation coefficients for the data at t

          !3.1) create the interpolation coefficients for the nodes
          y_map_inter(1) = bf_y_map0(j0_inter1)
          y_map_inter(2) = bf_y_map0(j0_inter2)
          
          nodes_inter(1,:) = bf_nodes0(i0_inter,j0_inter1,:)
          nodes_inter(2,:) = bf_nodes0(i0_inter,j0_inter2,:)

          inter_nodes0 = get_interpolation_coeff_1D(y_map_inter,nodes_inter)


          !3.2) create the interpolation coefficients for the transverse terms
          nodes_inter(1,:) = MATMUL(
     $         p_model%compute_x_gradient(bf_nodes0,i0_inter,j0_inter1,gradient_x,dx),
     $         p_model%compute_y_transM(bf_nodes0(i0_inter,j0_inter1,:)))

          nodes_inter(2,:) = MATMUL(
     $         p_model%compute_x_gradient(bf_nodes0,i0_inter,j0_inter2,gradient_x,dx),
     $         p_model%compute_y_transM(bf_nodes0(i0_inter,j0_inter2,:)))

          inter_trans0 = get_interpolation_coeff_1D(y_map_inter,nodes_inter)


          !4) create the interpolation coefficients for the data at t+dt
          y_map_inter(1) = bf_y_map1(j1_inter1)
          y_map_inter(2) = bf_y_map1(j1_inter2)

          nodes_inter(1,:) = MATMUL(
     $         p_model%compute_x_gradient(bf_nodes1,i1_inter,j1_inter1,gradient_x,dx),
     $         p_model%compute_y_transM(bf_nodes1(i1_inter,j1_inter1,:)))

          nodes_inter(2,:) = MATMUL(
     $         p_model%compute_x_gradient(bf_nodes1,i1_inter,j1_inter2,gradient_x,dx),
     $         p_model%compute_y_transM(bf_nodes1(i1_inter,j1_inter2,:)))

          inter_trans1 = get_interpolation_coeff_1D(y_map_inter,nodes_inter)

          t_amp1 = interpolate_1D(y1,inter_trans1)


          !5) evaluate the eigenvalues at t+dt
          eigenvalues_y = p_model%compute_y_eigenvalues(bf_nodes1(i1,j_eigen,:))


          !6) determine the left eigenvector corresponding to the eigenvalue
          left_eigenM = p_model%compute_y_lefteigenvector(bf_nodes1(i1,j_eigen,:))
             

          !7) determine the characteristic amplitude
          x0 = x1
          do k=1,ne


             !7.1) determine the position where the characteristic
             !     amplitude should be estimated
             y0 = y1 - eigenvalues_y(k)*dt


             !7.2) determine the normal and transverse contributions of
             !     the hyperbolic terms to the characteristic amplitude
             if(side_y.eq.right) then

                if(eigenvalues_y(k).ge.0) then
                   
                   n_amp0 = interpolate_1D(y0,inter_nodes0)
                   t_amp0 = interpolate_1D(y0,inter_trans0)
                   
                else

                   n_amp0 = p_model%get_far_field(t,x0,y0)
                   t_amp0 = [0.0d0,0.0d0,0.0d0]

                end if

             else
                
                if(eigenvalues_y(k).gt.0) then

                   n_amp0 = p_model%get_far_field(t,x0,y0)
                   t_amp0 = [0.0d0,0.0d0,0.0d0]

                else

                   n_amp0 = interpolate_1D(y0,inter_nodes0)
                   t_amp0 = interpolate_1D(y0,inter_trans0)

                end if
             end if


             !7.3) combine the information on the nodes at t-dt and the approximation
             !     of the integration of the transverse terms from t-dt to t
             amp =
     $            n_amp0 -
     $            compute_NewtonCotes_integration(t_amp0, t_amp1, dt)

             
             !7.4) compute the scalar product of the left eigenvector corresponding
             !     to the eigenvalue with the characteristic amplitude
             char_amp(k) = DOT_PRODUCT(amp,left_eigenM(:,k))
             
          end do


          !8) determine the right eigenmatrix
          right_eigenM = p_model%compute_y_righteigenvector(bf_nodes1(i1,j_eigen,:))


          !9) determine the new grid point
          bf_nodes1(i1,j1,:) = MATMUL(char_amp,right_eigenM)

        end subroutine compute_newgrdpt_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid point obtained by extension of the
        !> computational domain in the n1-direction
        !
        !> @date
        !> 14_11_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> physical model
        !
        !>@param dt
        !> time step
        !              
        !>@param bf_align0
        !> alignment of the buffer layer at t=t-dt
        !
        !>@param bf_x_map0
        !> x-coordinates of the buffer layer at t=t-dt
        !
        !>@param bf_y_map0
        !> y-coordinates of the buffer layer at t=t-dt
        !
        !>@param bf_nodes0
        !> nodes of the buffer layer at t=t-dt
        !
        !>@param bf_align1
        !> alignment of the buffer layer at t=t
        !
        !>@param bf_x_map1
        !> x-coordinates of the buffer layer at t=t
        !
        !>@param bf_y_map1
        !> y-coordinates of the buffer layer at t=t
        !
        !>@param bf_nodes1
        !> nodes of the buffer layer at t=t
        !              
        !>@param i1
        !> x-index identifying the new grdpt at t=t
        !
        !>@param j1
        !> y-index identifying the new grdpt at t=t
        !              
        !>@param side_x
        !> logical identifying the type of boundary (E or W)
        !
        !>@param gradient_y
        !> gradient procedure applied to compute the
        !> the transverse terms
        !--------------------------------------------------------------
        subroutine compute_newgrdpt_xy(
     $     p_model, t, dt,
     $     bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $     bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $     i1,j1,
     $     n_direction,
     $     incoming_proc,
     $     interpolation_indices)

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
          integer                            , intent(in)    :: n_direction
          procedure(incoming_proc)                           :: is_incoming
          procedure(gradient_n_proc)         , intent(in)    :: gradient_n_index1 !procedure for computing the gradient at the grid point of index1
          procedure(gradient_n_proc)         , intent(in)    :: gradient_n_index2 !procedure for computing the gradient at the grid point of index2
          procedure(gradient_n_proc)         , intent(in)    :: gradient_n_index3 !procedure for computing the gradient at the grid point of index3
          integer(ikind), dimension(2)       , intent(in)    :: eigen_indices
          integer(ikind), dimension(2,3)     , intent(in)    :: inter_indices1

          integer                       :: k          
          integer                       :: dir, dir2
          real(rkind), dimension(ne)    :: eigenvalues_n
          real(rkind), dimension(ne,ne) :: left_eigenM
          real(rkind), dimension(ne,ne) :: right_eigenM
          real(rkind), dimension(ne)    :: n_amp0
          real(rkind), dimension(ne)    :: t_amp0
          real(rkind), dimension(ne)    :: t_amp1
          real(rkind), dimension(ne)    :: amp
          real(rkind), dimension(ne)    :: char_amp

          real(rkind)                     :: n0,n1
          integer(ikind), dimension(2,3)  :: inter_indices0
          real(rkind)   , dimension(3)    :: n1_inter     !intermediate array where the n1-coordinates of the interpolation points are saved
          real(rkind)   , dimension(3)    :: n2_inter     !intermediate array where the n2-coordinates of the interpolation points are saved
          real(rkind)   , dimension(3,ne) :: nodes_inter  !intermediate array where the data of the interpolation points are saved
          real(rkind)   , dimension(3,ne) :: inter_nodes0 !interpolation coefficient for the plane by the nodes at t
          real(rkind)   , dimension(3,ne) :: inter_trans0 !interpolation coefficient for the plane by the transverse terms at t
          real(rkind)   , dimension(3,ne) :: inter_trans1 !interpolation coefficient for the plane by the transverse terms at t+dt

          !0) determine the direction
          dir  = x_direction
          dir2 = y_direction


          !1) determine the (x,y)-coordinates of the new grid point computed
          x1 = bf_x_map1(i1)
          y1 = bf_y_map(j1)


          !2) convert them into (n1,n2) coordinates
          n1_1 = get_n1_coord(x1,y1)
          n2_1 = get_n2_coord(x1,y1)


          !3) determine where the eigenvalues are evaluated
          !   and which indices are needed for the interpolation
          !   of the grid points
          !   i.e. convert the interpolation indices at t into
          !   interpolation indices at t-dt
          do k=1, 3
             
             inter_indices0(1,k) = bf_align1(dir,1) - bf_align0(dir,1) + inter_indices1(1,k)
             inter_indices0(2,k) = bf_align1(dir,2) - bf_align0(dir,2) + inter_indices1(2,k)

          end do


          !3) create the interpolation coefficients for the data at t-dt

          !3.1) create the interpolation coefficients for the nodes
          do k=1,3
             
             i_x_inter        = inter_indices0(1,k)
             i_y_inter        = inter_indices0(2,k)

             x_inter          = bf_x_map0(i_x_inter)
             y_inter          = bf_y_map0(i_y_inter)

             n1_inter(k)      = get_n1_coords(x_inter,y_inter)
             n2_inter(k)      = get_n2_coords(x_inter,y_inter)
             nodes_inter(k,:) = bf_nodes0(i_x_inter,i_y_inter,:)

          end do
          inter_nodes0 = get_interpolation_coeff_2D(n1_inter,n2_inter,nodes_inter)


          !3.2) create the interpolation coefficients for the transverse terms

          !3.2.1) create the data at the interpolation grid points
          call get_n_transverse_data_for_interpolation(
     $         p_model,
     $         bf_nodes0,
     $         inter_indices0,
     $         n_direction,
     $         gradient_n_index1,
     $         gradient_n_index2,
     $         gradient_n_index3,
     $         nodes_inter)
          
          !3.2.2) create the interpolation plane for the
          !       contribution of the transverse term at t-dt
          inter_trans0 = get_interpolation_coeff_2D(
     $         n1_map_inter,
     $         n2_map_inter,
     $         nodes_inter)


          !4) create the interpolation coefficients for the data at t+dt

          !4.1) create the coordinate maps (n1,n2) identifying the
          !     position of the interpolation points
          do k=1,3
             
             i_x_inter        = inter_indices1(1,k)
             i_y_inter        = inter_indices1(2,k)

             x_inter          = bf_x_map1(i_x_inter)
             y_inter          = bf_y_map1(i_y_inter)

             n1_inter(k)      = get_n1_coords(x_inter,y_inter)
             n2_inter(k)      = get_n2_coords(x_inter,y_inter)

          end do

          !4.2) compute the transverse terms at the interpolation
          !     points
          call get_n_transverse_data_for_interpolation(
     $         p_model,
     $         bf_nodes1,
     $         inter_indices1,
     $         n_direction,
     $         gradient_n_index1,
     $         gradient_n_index2,
     $         gradient_n_index3,
     $         dx,dy,
     $         nodes_inter1)

          !4.3) create the interpolation plane for the
          !     contribution of the transverse term at t
          inter_trans1 = get_interpolation_coeff_2D(
     $         n1_map_inter,
     $         n2_map_inter,
     $         nodes_inter)

          t_amp1 = interpolate_2D(n1_1,n2_1,inter_trans1)


          !5) evaluate the eigen data at t
          select case(n_direction)

            case(n1_direction)
               eigenvalues_n = p_model%compute_n1_eigenvalues(
     $              bf_nodes1(indices_eigen(1),indices_eigen(2),:))

               left_eigenM = p_model%compute_n1_lefteigenvector(
     $              bf_nodes1(indices_eigen(1),indices_eigen(2),:))

               right_eigenM = p_model%compute_n1_righteigenvector(
     $              bf_nodes1(indices_eigen(1),indices_eigen(2),:))

            case(n2_direction)
               eigenvalues_n = p_model%compute_n2_eigenvalues(
     $              bf_nodes1(indices_eigen(1),indices_eigen(2),:))

               left_eigenM = p_model%compute_n2_lefteigenvector(
     $              bf_nodes1(indices_eigen(1),indices_eigen(2),:))

               right_eigenM = p_model%compute_n2_righteigenvector(
     $              bf_nodes1(indices_eigen(1),indices_eigen(2),:))

            case default
               print '(''bf_newgrdpt_class'')'
               print '(''compute_newgrdpt_xy'')'
               print '(''direction not recognized: '', I2)', n_direction
               stop ''

          end select


          !7) determine the characteristic amplitude
          select case(n_direction)
            case(n1_direction)

               n2_0 = n2_1

               do k=1,ne
                  
                 !7.1) determine the position where the characteristic
                 !     amplitude should be estimated
                 n1_0 = n1_1 - eigenvalues_n(k)*dt


                 !7.2) determine the normal and transverse contributions of
                 !     the hyperbolic terms to the characteristic amplitude
                 if(side_x.eq.right) then

                    if(eigenvalues_n(k).ge.0) then
                   
                       n_amp0 = interpolate_2D(n1_0, n2_0, inter_nodes0)
                       t_amp0 = interpolate_2D(n1_0, n2_0, inter_trans0)
                       
                    else

                       x0 = get_x_coord(n1_0,n2_0)
                       y0 = get_y_coord(n1_0,n2_0)

                       n_amp0 = p_model%get_far_field(t,x0,y0)
                       t_amp0 = [0.0d0,0.0d0,0.0d0]
                       
                    end if
                    
                 else
                    
                    if(eigenvalues_x(k).gt.0) then

                       x0 = get_x_coord(n1_0,n2_0)
                       y0 = get_y_coord(n1_0,n2_0)

                       n_amp0 = p_model%get_far_field(t,x0,y0)
                       t_amp0 = [0.0d0,0.0d0,0.0d0]

                    else

                       n_amp0 = interpolate_2D(n1_0, n2_0, inter_nodes0)
                       t_amp0 = interpolate_2D(n1_0, n2_0, inter_trans0)

                    end if
                 end if


                 !7.3) combine the information on the nodes at t-dt and the approximation
                 !     of the integration of the transverse terms from t-dt to t
                 amp =
     $                n_amp0 -
     $                compute_NewtonCotes_integration(t_amp0, t_amp1, dt)

                 
                 !7.4) compute the scalar product of the left eigenvector corresponding
                 !     to the eigenvalue with the characteristic amplitude
                 char_amp(k) = DOT_PRODUCT(amp,left_eigenM(:,k))
                 
              end do

           case(n2_direction)

              n1_0 = n1_1

              do k=1,ne
                 
                 !7.1) determine the position where the characteristic
                 !     amplitude should be estimated
                 n2_0 = n2_1 - eigenvalues_n(k)*dt


                 !7.2) determine the normal and transverse contributions of
                 !     the hyperbolic terms to the characteristic amplitude
                 if(side_x.eq.right) then

                    if(eigenvalues_n(k).ge.0) then
                       
                       n_amp0 = interpolate_2D(n1_0, n2_0, inter_nodes0)
                       t_amp0 = interpolate_2D(n1_0, n2_0, inter_trans0)
                       
                    else

                       x0 = get_x_coord(n1_0,n2_0)
                       y0 = get_y_coord(n1_0,n2_0)

                       n_amp0 = p_model%get_far_field(t,x0,y0)
                       t_amp0 = [0.0d0,0.0d0,0.0d0]
                       
                    end if
                    
                 else
                    
                    if(eigenvalues_x(k).gt.0) then

                       x0 = get_x_coord(n1_0,n2_0)
                       y0 = get_y_coord(n1_0,n2_0)

                       n_amp0 = p_model%get_far_field(t,x0,y0)
                       t_amp0 = [0.0d0,0.0d0,0.0d0]

                    else

                       n_amp0 = interpolate_2D(n1_0, n2_0, inter_nodes0)
                       t_amp0 = interpolate_2D(n1_0, n2_0, inter_trans0)

                    end if
                 end if


                 !7.3) combine the information on the nodes at t-dt and the approximation
                 !     of the integration of the transverse terms from t-dt to t
                 amp =
     $                n_amp0 -
     $                compute_NewtonCotes_integration(t_amp0, t_amp1, dt)

                 
                 !7.4) compute the scalar product of the left eigenvector corresponding
                 !     to the eigenvalue with the characteristic amplitude
                 char_amp(k) = DOT_PRODUCT(amp,left_eigenM(:,k))
                 
              end do

             case default
                print '(''bf_newgrdpt_class'')'
                print '(''compute_newgrdpt_xy'')'
                print '(''direction not recognized: ''I2)', n_direction
                stop ''

           end select


          !9) determine the new grid point
          bf_nodes1(i1,j1,:) = MATMUL(char_amp,right_eigenM)


        end subroutine compute_newgrdpt_xy


        !compute the contribution of the transverse terms
        !at the location of the interpolation points
        subroutine get_n_transverse_data_for_interpolation(
     $     p_model,
     $     bf_nodes,
     $     inter_indices,
     $     n_direction,
     $     gradient_n_index1,
     $     gradient_n_index2,
     $     gradient_n_index3,
     $     dx,dy,
     $     nodes_inter)

          implicit none

          type(pmodel_eq)                  , intent(in)  :: p_model
          real(rkind)    , dimension(:,:,:), intent(in)  :: bf_nodes
          integer(ikind) , dimension(2,3)  , intent(in)  :: inter_indices
          integer                          , intent(in)  :: n_direction
          procedure(gradient_n_proc)       , intent(in)  :: gradient_n_index1
          procedure(gradient_n_proc)       , intent(in)  :: gradient_n_index2
          procedure(gradient_n_proc)       , intent(in)  :: gradient_n_index3
          real(rkind)                      , intent(in)  :: dx
          real(rkind)                      , intent(in)  :: dy
          real(rkind)    , dimension(3,ne) , intent(out) :: nodes_inter

          
          real(rkind) :: dn
          integer     :: k
          

          select case(n_direction)
            case(n1_direction)
               
               dn = get_dn2(dx,dy)

               !compute the transverse terms along the n2 direction at the grid point 1
               k=1
               nodes_inter(1,:) = MATMUL(
     $                 p_model%compute_n_gradient(
     $                        bf_nodes,
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),
     $                        gradient_n_index1,
     $                        dn),
     $                 p_model%compute_n1_transM(
     $                        bf_nodes(
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),:)))

               !compute the transverse terms along the n2 direction at the grid point 2
               k=2
               nodes_inter(1,:) = MATMUL(
     $                 p_model%compute_n_gradient(
     $                        bf_nodes,
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),
     $                        gradient_n_index2,
     $                        dn),
     $                 p_model%compute_n1_transM(
     $                        bf_nodes(
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),:)))

               !compute the transverse terms along the n3 direction at the grid point 3
               k=3
               nodes_inter(1,:) = MATMUL(
     $                 p_model%compute_n_gradient(
     $                        bf_nodes,
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),
     $                        gradient_n_index3,
     $                        dn),
     $                 p_model%compute_n1_transM(
     $                        bf_nodes(
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),:)))

             case(n2_direction)
               
               dn = get_dn1(dx,dy)

               !compute the transverse terms along the n1 direction at the grid point 1
               k=1
               nodes_inter(1,:) = MATMUL(
     $                 p_model%compute_n_gradient(
     $                        bf_nodes,
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),
     $                        gradient_n_index1,
     $                        dn),
     $                 p_model%compute_n2_transM(
     $                        bf_nodes(
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),:)))

               !compute the transverse terms along the n1 direction at the grid point 2
               k=2
               nodes_inter(1,:) = MATMUL(
     $                 p_model%compute_n_gradient(
     $                        bf_nodes,
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),
     $                        gradient_n_index2,
     $                        dn),
     $                 p_model%compute_n2_transM(
     $                        bf_nodes(
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),:)))

               !compute the transverse terms along the n1 direction at the grid point 3
               k=3
               nodes_inter(1,:) = MATMUL(
     $                 p_model%compute_n_gradient(
     $                        bf_nodes,
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),
     $                        gradient_n_index2,
     $                        dn),
     $                 p_model%compute_n2_transM(
     $                        bf_nodes(
     $                        inter_indices(1,k),
     $                        inter_indices(2,k),:)))

            case default
               print '(''compute_newgrdpt_xy'')'
               print '(''direction not recognized: '',I2)', n_direction
               stop ''

          end select

        end subroutine get_n_transverse_data_for_interpolation


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the interpolation coefficients for a 1st order
        !> polynomial fit: get (a,b) such that:
        !> a*x_map(1)+b = nodes(1,k)
        !
        !> @date
        !> 14_11_2014 - initial version - J.L. Desmarais
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param nodes
        !> interpolation points
        !              
        !>@return inter_coeff
        !> (a,b) for each governing variable
        !--------------------------------------------------------------
        function get_interpolation_coeff_1D(x_map,nodes)
     $     result(inter_coeff)
        
          implicit none

          real(rkind), dimension(2)   , intent(in) :: x_map
          real(rkind), dimension(2,ne), intent(in) :: nodes
          real(rkind), dimension(2,ne)             :: inter_coeff

          integer :: k

          do k=1, ne

             inter_coeff(1,k) = (nodes(2,k) - nodes(1,k))/(x_map(2)-x_map(1))
             inter_coeff(2,k) = nodes(1,k) - inter_coeff(1,k)*x_map(1)

          end do

        end function get_interpolation_coeff_1D


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> from the coefficients (a,b) for each governing variable
        !> compute ax+b
        !
        !> @date
        !> 14_11_2014 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !
        !>@param inter_coeff
        !> coefficients (a,b) for each governing variable
        !              
        !>@return nodes_inter
        !> nodes interpolated at x
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the interpolation coefficients for a 1st order
        !> polynomial fit: get (a,b,c) such that:
        !> a*x_map(1)+b*y_map(1) + c = nodes(1,k)
        !
        !> @date
        !> 14_11_2014 - initial version - J.L. Desmarais
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates        
        !
        !>@param nodes
        !> interpolation points
        !              
        !>@return inter_coeff
        !> (a,b,c) for each governing variable
        !--------------------------------------------------------------
        function get_interpolation_coeff_2D(x_map,y_map,nodes)
     $     result(inter_coeff)
        
          implicit none

          real(rkind), dimension(3)   , intent(in) :: x_map
          real(rkind), dimension(3)   , intent(in) :: y_map
          real(rkind), dimension(3,ne), intent(in) :: nodes
          real(rkind), dimension(3,ne)             :: inter_coeff

          integer                   :: k
          real(rkind), dimension(3) :: A
          real(rkind), dimension(3) :: B
          real(rkind), dimension(3) :: C
          real(rkind), dimension(3) :: n

          !create the vector identifying the
          !points on the plane: coordinates (x,y)
          A(1) = x_map(1)
          A(2) = y_map(1)

          B(1) = x_map(2)
          B(2) = y_map(2)

          C(1) = x_map(3)
          C(2) = y_map(3)

          do k=1, ne

             !create the vector identifying the
             !points on the plane: data (z)
             A(3) = nodes(1,k)
             B(3) = nodes(2,k)
             C(3) = nodes(3,k)

             !get the vector normal to the plane
             n = get_plane_normal_vector(A,B,C)

             !compute the coefficient a,b,c such that
             !the nodes are given by the equation:
             !nodes(1) = a*x + b*y + c
             inter_coeff(1,k) = -n(1)/n(3)
             inter_coeff(2,k) = -n(2)/n(3)
             inter_coeff(3,k) =  -A(1)*inter_coeff(1,k) - A(2)*inter_coeff(2,k) + A(3)

          end do

        end function get_interpolation_coeff_2D


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> from the coefficients (a,b,c) for each governing variable
        !> compute ax+by+c
        !
        !> @date
        !> 14_11_2014 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !
        !>@param y
        !> y-coordinate
        !
        !>@param inter_coeff
        !> coefficients (a,b,c) for each governing variable
        !              
        !>@return nodes_inter
        !> nodes interpolated at (x,y)
        !--------------------------------------------------------------
        function interpolate_2D(
     $     x,
     $     y,
     $     inter_coeff)
     $     result(nodes_inter)

          implicit none

          real(rkind)                 , intent(in) :: x
          real(rkind)                 , intent(in) :: y
          real(rkind), dimension(3,ne), intent(in) :: inter_coeff
          real(rkind), dimension(ne)               :: nodes_inter

          integer :: k

          do k=1,ne
             nodes_inter(k) = x*inter_coeff(1,k) + y*inter_coeff(2,k) + inter_coeff(3,k)
          end do

        end function interpolate_2D


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get a vector normal to the plane created by three points
        !
        !> @date
        !> 14_11_2014 - initial version - J.L. Desmarais
        !
        !>@param A
        !> first point (x_A,y_A,z_A)
        !
        !>@param B
        !> second point (x_B,y_B,z_B)
        !
        !>@param C
        !> third point (x_C,y_C,z_C)
        !
        !>@return n
        !> normal vector
        !--------------------------------------------------------------
        function get_plane_normal_vector(A,B,C) result(n)

          implicit none

          real(rkind), dimension(3), intent(in) :: A
          real(rkind), dimension(3), intent(in) :: B
          real(rkind), dimension(3), intent(in) :: C
          real(rkind), dimension(3)             :: n

          real(rkind), dimension(3) :: AB
          real(rkind), dimension(3) :: AC

          integer :: k

          do k=1,3
             AB(k) = B(k) - A(k)
             AC(k) = C(k) - A(k)
          end do

          n = cross_product(AB,AC)          

        end function get_plane_normal_vector


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the cross product between 2 vectors
        !
        !> @date
        !> 14_11_2014 - initial version - J.L. Desmarais
        !
        !>@param v1
        !> first vector
        !
        !>@param v2
        !> second vector
        !
        !>@return var
        !> cross product
        !--------------------------------------------------------------
        function cross_product(v1,v2) result(var)

          implicit none

          real(rkind), dimension(3), intent(in) :: v1
          real(rkind), dimension(3), intent(in) :: v2
          real(rkind), dimension(3)             :: var

          var(1) =  v1(2)*v2(3) - v1(3)*v2(2)
          var(2) = -v1(1)*v2(3) + v1(3)*v2(1)
          var(3) =  v1(1)*v2(2) - v1(2)*v2(1)

        end function cross_product


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> integrate a function between two points using Newton-Cotes
        !
        !> @date
        !> 14_11_2014 - initial version - J.L. Desmarais
        !
        !>@param data0
        !> data at t=t0
        !
        !>@param data1
        !> data at t=t1
        !              
        !>@param dt
        !> time step dt=t1-t0
        !              
        !>@return data_integrated
        !> integration from data0 to data1
        !--------------------------------------------------------------
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

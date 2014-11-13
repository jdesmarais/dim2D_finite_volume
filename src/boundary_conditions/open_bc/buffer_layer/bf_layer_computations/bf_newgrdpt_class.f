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

        implicit none

        private
        public :: bf_newgrdpt

        !object encapsulating the subroutines for the computation of
        !the newgrid points
        type :: bf_newgrdpt

          contains

          procedure, nopass :: compute_newgrdpt_x
          procedure, nopass :: interpolate_char_amp_x
          procedure, nopass :: interpolate1D
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
     $       i1,j1, side_x, gradient_type)

          implicit none

          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)                        , intent(in)    :: t
          real(rkind)                        , intent(in)    :: dt
          integer(ikind), dimension(4)       , intent(in)    :: bf_align0
          real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map0
          real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map0
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes0
          integer(ikind), dimension(4)       , intent(in)    :: bf_align1
          real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map1
          real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map1
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: bf_nodes1
          integer(ikind)                     , intent(in)    :: i1
          integer(ikind)                     , intent(in)    :: j1
          logical                            , intent(in)    :: side_x
          integer                            , intent(in)    :: gradient_type

          integer                       :: k
          
          integer                       :: dir, dir2
          integer(ikind)                :: i_eigen
          real(rkind), dimension(ne)    :: eigenvalues_x
          real(rkind), dimension(ne,ne) :: right_eigenM
          real(rkind), dimension(ne)    :: char_amp

          real(rkind)                  :: x0,x1
          integer(ikind)               :: i0_inter1, i0_inter2, j0_inter
          integer(ikind)               :: i1_inter1, i1_inter2, j1_inter
          real(rkind), dimension(2)    :: x_map_inter
          real(rkind), dimension(2,ne) :: nodes_inter
          

          !0) determine the direction
          dir  = x_direction
          dir2 = y_direction


          !1) determine the x-coordinate of the new grid point computed
          x1   = bf_x_map1(i1)


          !2) determine where the eigenvalues are evaluated
          if(side_x.eqv.right) then
             i_eigen = i1-1
          else
             i_eigen = i1+1
          end if


          !3) evaluate the eigenvalues
          eigenvalues_x = p_model%compute_x_eigenvalues(bf_nodes1(i_eigen,j,:))


          !4) determine the characteristic amplitude
          do k=1,ne


             !4.1) determine the position where the characteristic
             !     amplitude should be estimated
             x0 = bf_x_map1(i1) - eigenvalues_x(k)*dt


             !4.2) determine the normal and transverse contributions of
             !     the hyperbolic terms to the characteristic amplitude
             if(side_x.eq.right) then

                if(eigenvalues_x(k).ge.0) then
                   
                   i0_inter1 = bf_align1(dir,1) -bf-align0(dir,1) +i_eigen-2
                   i0_inter2 = bf_align1(dir,1) -bf-align0(dir,1) +i_eigen-1
                   j0_inter  = bf_align1(dir2,1)-bf_align0(dir2,1)+j1

                   i1_inter1 = i_eigen-2
                   i1_inter2 = i_eigen-1
                   j1_inter  = j1

                   char_map = interpolate_char_amp_x(
     $                  p_model,t,dt,
     $                  x0, i0_inter1, i0_inter2, j0_inter,
     $                  bf_x_map0, bf_y_map0, bf_nodes0,
     $                  x1, i1_inter1, i1_inter2, j1_inter,
     $                  bf_x_map1, bf_y_map1, bf_nodes1)

                else

                   char_amp = p_model%get_far_field(t,bf_x_map1(i1),bf_y_map1(j1))

                end if

             else
                
                if(eigenvalues_x(k).gt.0) then

                   char_map = p_model%get_far_field(t,bf_x_map1(i1),bf_y_map1(j1))

                else

                   i0_inter1 = bf_align1(dir,1) -bf-align0(dir,1) +i_eigen+1
                   i0_inter2 = bf_align1(dir,1) -bf-align0(dir,1) +i_eigen+2
                   j0_inter  = bf_align1(dir2,1)-bf_align0(dir2,1)+j1

                   i1_inter1 = i_eigen+1
                   i1_inter2 = i_eigen+2
                   j1_inter  = j1

                   char_map = interpolate_char_amp_x(
     $                  p_model,t,dt,
     $                  x0, i0_inter1, i0_inter2, j0_inter,
     $                  bf_x_map0, bf_y_map0, bf_nodes0,
     $                  x1, i1_inter1, i1_inter2, j1_inter,
     $                  bf_x_map1, bf_y_map1, bf_nodes1)

                end if
             end if
          end do


          !4) determine the right eigenmatrix
          right_eigenM = p_model%compute_x_righteigenvector(int_nodes1(i_eigen,j1,:))


          !5) determine the new grid point
          int_nodes1(i1,j1,:) = MATMUL(char_amp,right_eigenM)


        end subroutine compute_newgrdpt_x


        !interpolate the characteristic amplitude in the x-direction
        function interpolate_char_amp_x(
     $     p_model, dt,
     $     x0, i0_inter1, i0_inter2, j0_inter,
     $     bf_x_map0, bf_y_map0, bf_nodes0,
     $     x1, i1_inter1, i1_inter2, j1_inter,
     $     bf_x_map1, bf_y_map1, bf_nodes1,
     $     gradient_y)
     $     result(char_amp)
        
          implicit none

          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: dt

          real(rkind)                  , intent(in) :: x0
          integer(ikind)               , intent(in) :: i0_inter1
          integer(ikind)               , intent(in) :: i0_inter2
          integer(ikind)               , intent(in) :: j0_inter
          real(rkind)    , dimension(:), intent(in) :: bf_x_map0
          real(rkind)    , dimension(:), intent(in) :: bf_y_map0
          real(rkind)    , dimension(:), intent(in) :: bf_nodes0

          real(rkind)                  , intent(in) :: x1
          integer(ikind)               , intent(in) :: i1_inter1
          integer(ikind)               , intent(in) :: i1_inter2
          integer(ikind)               , intent(in) :: j1_inter
          real(rkind)    , dimension(:), intent(in) :: bf_x_map1
          real(rkind)    , dimension(:), intent(in) :: bf_y_map1
          real(rkind)    , dimension(:), intent(in) :: bf_nodes1

          procedure(gradient_y_proc)                :: gradient_y

          real(rkind)    , dimension(ne)            :: char_amp


          real(rkind), dimension(2)    :: x_map_inter
          real(rkind), dimension(2,ne) :: nodes_inter
          real(rkind), dimension(ne)   :: n_char_amp
          real(rkind), dimension(ne)   :: t_char_amp0
          real(rkind), dimension(ne)   :: t_char_amp1
          real(rkind), dimension(ne)   :: t_char_amp


          !x_map for the interpolation at t-dt
          x_map_inter(1) = bf_x_map0(i0_inter1)
          x_map_inter(2) = bf_x_map0(i0_inter2)


          !computation of the normal contribution
          nodes_inter(1,:) = bf_nodes0(i0_inter1,j_inter,:)
          nodes_inter(2,:) = bf_nodes0(i0_inter2,j_inter,:)

          n_char_amp = interpolate_1D(x0,x_map_inter,nodes_inter)


          !computation of the transverse contribution

          !   computation of the transverse component at t-dt
          nodes_inter(1,:) = MATMUL(
     $         p_model%compute_y_gradient(bf_nodes0,i0_inter1,j0_inter,gradient_y,dy),
     $         p_model%compute_x_transM(i0_inter1,j0_inter,bf_nodes0))

          nodes_inter(2,:) = MATMUL(
     $         p_model%compute_y_gradient(bf_nodes0,i0_inter2,j0_inter,gradient_y,dy),
     $         p_model%compute_x_transM(i0_inter2,j0_inter,bf_nodes0))
          
          t_char_amp0 = interpolate_1D(x0,x_map_inter,nodes_inter)

          
          !   computation of the transverse component at t
          x_map_inter(1) = bf_x_map1(i1_inter1)
          x_map_inter(2) = bf_x_map1(i1_inter2)

          nodes_inter(1,:) = MATMUL(
     $         p_model%compute_y_gradient(bf_nodes1,i1_inter1,j1_inter,gradient_y,dy),
     $         p_model%compute_x_transM(i1_inter1,j1_inter,bf_nodes1))

          nodes_inter(2,:) = MATMUL(
     $         p_model%compute_y_gradient(bf_nodes1,i1_inter2,j1_inter,gradient_y,dy),
     $         p_model%compute_x_transM(i1_inter2,j1_inter,bf_nodes1))

          t_char_amp1 = interpolate_1D(x1,x_map_inter,nodes_inter)

          !   time integration using Newton-Cotes
          t_char_amp = compute_NewtonCotes_integration(
     $         t_char_amp0, t_char_map1, dt)
          

          !adding the two contributions: normal + transverse
          do k=1,ne
             char_amp(k) = n_char_amp(k) + t_char_amp(k)
          end do

        end function interpolate_char_amp_x


        !interpolate in 1D using a 1st order approximation
        function interpolate1D(
     $     x,
     $     x_map,
     $     nodes)
     $     result(nodes_inter)

          implicit none

          real(rkind)                 , intent(in) :: x
          real(rkind), dimension(2)   , intent(in) :: x_map
          real(rkind), dimension(2,ne), intent(in) :: nodes
          real(rkind), dimension(ne)               :: nodes_inter

          real(rkind) :: a
          integer     :: k

          a = (x-x_map(1))/(x_map(2)-x_map(1))

          do k=1,ne
             nodes_inter(k) = nodes(1,k) + a*(nodes(2,k)-nodes(1,k))
          end do

        end function interpolate1D


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

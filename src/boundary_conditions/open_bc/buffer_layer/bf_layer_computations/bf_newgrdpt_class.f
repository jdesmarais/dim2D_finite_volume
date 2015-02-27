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
      !> 14_11_2014 - initial version           - J.L. Desmarais
      !> 06_02_2015 - using primitive variables - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_newgrdpt_class

        use bf_layer_newgrdpt_procedure_module, only :
     $       error_gradient_type

        use interface_primary, only :
     $       gradient_proc

        use n_coords_module, only :
     $       get_dn,
     $       get_x_coord,
     $       get_y_coord,
     $       get_n1_coord,
     $       get_n2_coord

        use openbc_operators_module, only :
     $       incoming_proc        

        use parameters_bf_layer, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       gradient_I_type,
     $       gradient_L0_type,
     $       gradient_R0_type,
     $       gradient_xLR0_yI_type,
     $       gradient_xI_yLR0_type,
     $       gradient_xLR0_yLR0_type

        use parameters_constant, only :
     $       left, right,
     $       x_direction, y_direction,
     $       n1_direction, n2_direction

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq        

        use sd_operators_fd_module, only :
     $       gradient_x_interior,
     $       gradient_y_interior,
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_R0

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

          procedure, nopass :: compute_newgrdpt
          
          procedure, nopass :: compute_newgrdpt_local
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
        !>@param t
        !> time
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
        !>@param procedure_type
        !> integer identifying the procedure that should be
        !> applied to compute the new grid point
        !
        !>@param gradient_type
        !> integer identifying the gradient procedure needed
        !> to compute the transverse terms
        !--------------------------------------------------------------
        function compute_newgrdpt(
     $       p_model, t, dt,
     $       bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $       bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $       i1,j1,
     $       nb_procedures, procedure_type, gradient_type)
     $       result(new_grdpt)

          implicit none

          type(pmodel_eq)                    , intent(in) :: p_model
          real(rkind)                        , intent(in) :: t
          real(rkind)                        , intent(in) :: dt
          integer(ikind), dimension(2,2)     , intent(in) :: bf_align0
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map0
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map0
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes0
          integer(ikind), dimension(2,2)     , intent(in) :: bf_align1
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map1
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map1
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes1
          integer(ikind)                     , intent(in) :: i1
          integer(ikind)                     , intent(in) :: j1
          integer                            , intent(in) :: nb_procedures
          integer       , dimension(4)       , intent(in) :: procedure_type
          integer       , dimension(4)       , intent(in) :: gradient_type
          real(rkind)   , dimension(ne)                   :: new_grdpt


          integer                      :: k
          integer                      :: l
          real(rkind), dimension(ne,4) :: new_grdpt_data


          if(nb_procedures.ge.1) then

             !compute the new grid point according to each procedure
             do k=1, nb_procedures

                new_grdpt_data(:,k) = compute_newgrdpt_local(
     $               p_model,t,dt,
     $               bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $               bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $               i1,j1,
     $               procedure_type(k), gradient_type(k))

             end do

             !compute the new grid point as an average of the
             !grid points given by the procedures
             do l=1,ne
                new_grdpt(l) = 0.0d0
                do k=1,nb_procedures
                   new_grdpt(l) = new_grdpt(l) + new_grdpt_data(l,k)
                end do
                new_grdpt(l) = new_grdpt(l)/nb_procedures
             end do

          end if

        end function compute_newgrdpt


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
        !>@param t
        !> time
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
        !>@param procedure_type
        !> integer identifying the procedure that should be
        !> applied to compute the new grid point
        !
        !>@param gradient_type
        !> integer identifying the gradient procedure needed
        !> to compute the transverse terms
        !--------------------------------------------------------------
        function compute_newgrdpt_local(
     $       p_model, t, dt,
     $       bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $       bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $       i1,j1,
     $       procedure_type, gradient_type)
     $       result(new_grdpt)

          implicit none

          type(pmodel_eq)                    , intent(in) :: p_model
          real(rkind)                        , intent(in) :: t
          real(rkind)                        , intent(in) :: dt
          integer(ikind), dimension(2,2)     , intent(in) :: bf_align0
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map0
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map0
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes0
          integer(ikind), dimension(2,2)     , intent(in) :: bf_align1
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map1
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map1
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes1
          integer(ikind)                     , intent(in) :: i1
          integer(ikind)                     , intent(in) :: j1
          integer                            , intent(in) :: procedure_type
          integer                            , intent(in) :: gradient_type
          real(rkind)   , dimension(ne)                   :: new_grdpt

          integer                        :: n_direction
          logical                        :: side
          integer(ikind), dimension(2)   :: eigen_indices
          integer(ikind), dimension(2,3) :: inter_indices1


          select case(procedure_type)
          
            case(SW_corner_type)
               
               n_direction         = n2_direction
               side                = left
               eigen_indices       = [i1+1,j1+1]
               inter_indices1(:,1) = [i1+1,j1+1]
               inter_indices1(:,2) = [i1+2,j1+1]
               inter_indices1(:,3) = [i1+1,j1+2]

               new_grdpt = compute_newgrdpt_xy(
     $              p_model, t, dt,
     $              bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $              bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $              i1,j1,
     $              n_direction,
     $              side,
     $              gradient_x_x_oneside_L0,
     $              gradient_y_y_oneside_L0,
     $              gradient_x_interior,
     $              gradient_y_y_oneside_L0,
     $              gradient_x_x_oneside_L0,
     $              gradient_y_interior,
     $              eigen_indices,
     $              inter_indices1)


            case(SE_corner_type)

               n_direction         = n1_direction
               side                = right
               eigen_indices       = [i1-1,j1+1]
               inter_indices1(:,1) = [i1-2,j1+1]
               inter_indices1(:,2) = [i1-1,j1+1]
               inter_indices1(:,3) = [i1-1,j1+2]

               new_grdpt = compute_newgrdpt_xy(
     $              p_model, t, dt,
     $              bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $              bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $              i1,j1,
     $              n_direction,
     $              side,
     $              gradient_x_interior,
     $              gradient_y_y_oneside_L0,
     $              gradient_x_x_oneside_R0,
     $              gradient_y_y_oneside_L0,
     $              gradient_x_x_oneside_R0,
     $              gradient_y_interior,
     $              eigen_indices,
     $              inter_indices1)


            case(NW_corner_type)

               n_direction         = n1_direction
               side                = left
               eigen_indices       = [i1+1,j1-1]
               inter_indices1(:,1) = [i1+1,j1-2]
               inter_indices1(:,2) = [i1+1,j1-1]
               inter_indices1(:,3) = [i1+2,j1-1]

               new_grdpt = compute_newgrdpt_xy(
     $              p_model, t, dt,
     $              bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $              bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $              i1,j1,
     $              n_direction,
     $              side,
     $              gradient_x_x_oneside_L0,
     $              gradient_y_interior,
     $              gradient_x_x_oneside_L0,
     $              gradient_y_y_oneside_R0,
     $              gradient_x_interior,
     $              gradient_y_y_oneside_R0,
     $              eigen_indices,
     $              inter_indices1)


            case(NE_corner_type)

               n_direction         = n2_direction
               side                = right
               eigen_indices       = [i1-1,j1-1]
               inter_indices1(:,1) = [i1-1,j1-2]
               inter_indices1(:,2) = [i1-2,j1-1]
               inter_indices1(:,3) = [i1-1,j1-1]

               new_grdpt = compute_newgrdpt_xy(
     $              p_model, t, dt,
     $              bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $              bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $              i1,j1,
     $              n_direction,
     $              side,
     $              gradient_x_x_oneside_R0,
     $              gradient_y_interior,
     $              gradient_x_interior,
     $              gradient_y_y_oneside_R0,
     $              gradient_x_x_oneside_R0,
     $              gradient_y_y_oneside_R0,
     $              eigen_indices,
     $              inter_indices1)


            case(S_edge_type)

               side = left

               select case(gradient_type)

                 case(gradient_I_type)

                    new_grdpt = compute_newgrdpt_y(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_x_interior)

                 case(gradient_L0_type)

                    new_grdpt = compute_newgrdpt_y(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_x_x_oneside_L0)

                 case(gradient_R0_type)

                    new_grdpt = compute_newgrdpt_y(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_x_x_oneside_R0)

                 case default
                    call error_gradient_type(
     $                   'bf_newgrdpt_class',
     $                   'compute_newgrdpt',
     $                   gradient_type)

               end select

                    
            case(E_edge_type)

               side = right

               select case(gradient_type)

                 case(gradient_I_type)

                    new_grdpt = compute_newgrdpt_x(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_y_interior)

                 case(gradient_L0_type)

                    new_grdpt = compute_newgrdpt_x(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_y_y_oneside_L0)

                 case(gradient_R0_type)

                    new_grdpt = compute_newgrdpt_x(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_y_y_oneside_R0)

                 case default
                    call error_gradient_type(
     $                   'bf_newgrdpt_class',
     $                   'compute_newgrdpt',
     $                   gradient_type)

               end select


            case(W_edge_type)

               side = left

               select case(gradient_type)

                 case(gradient_I_type)

                    new_grdpt = compute_newgrdpt_x(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_y_interior)

                 case(gradient_L0_type)

                    new_grdpt = compute_newgrdpt_x(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_y_y_oneside_L0)

                 case(gradient_R0_type)

                    new_grdpt = compute_newgrdpt_x(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_y_y_oneside_R0)

                 case default
                    call error_gradient_type(
     $                   'bf_newgrdpt_class',
     $                   'compute_newgrdpt',
     $                   gradient_type)

               end select


            case(N_edge_type)

               side = right

               select case(gradient_type)

                 case(gradient_I_type)

                    new_grdpt = compute_newgrdpt_y(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_x_interior)

                 case(gradient_L0_type)

                    new_grdpt = compute_newgrdpt_y(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_x_x_oneside_L0)

                 case(gradient_R0_type)

                    new_grdpt = compute_newgrdpt_y(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1, side, gradient_x_x_oneside_R0)

                 case default
                    call error_gradient_type(
     $                   'bf_newgrdpt_class',
     $                   'compute_newgrdpt',
     $                   gradient_type)

               end select

            case(SE_edge_type)
               
               n_direction         = n1_direction
               side                = right
               eigen_indices       = [i1-1,j1+1]
               inter_indices1(:,1) = [i1-1,j1]
               inter_indices1(:,2) = [i1-1,j1+1]
               inter_indices1(:,3) = [i1  ,j1+1]

               select case(gradient_type)
                 case(gradient_I_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_interior,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_interior,
     $                   gradient_y_y_oneside_L0,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xI_yLR0_type)
                    
                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_y_oneside_L0,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_interior,
     $                   gradient_y_y_oneside_L0,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xLR0_yI_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_interior,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_y_oneside_L0,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xLR0_yLR0_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_y_oneside_L0,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_y_oneside_L0,
     $                   eigen_indices,
     $                   inter_indices1)                    

                 case default
                    call error_gradient_type(
     $                   'bf_newgrdpt_class',
     $                   'compute_newgrdpt',
     $                   gradient_type)
              end select


            case(SW_edge_type)
               
               n_direction         = n2_direction
               side                = left
               eigen_indices       = [i1+1,j1+1]
               inter_indices1(:,1) = [i1+1,j1]
               inter_indices1(:,2) = [i1  ,j1+1]
               inter_indices1(:,3) = [i1+1,j1+1]


               select case(gradient_type)

                 case(gradient_I_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_interior,
     $                   gradient_x_interior,
     $                   gradient_y_y_oneside_L0,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xLR0_yI_type)
                    
                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_interior,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_y_oneside_L0,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xI_yLR0_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_y_oneside_L0,
     $                   gradient_x_interior,
     $                   gradient_y_y_oneside_L0,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xLR0_yLR0_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_y_oneside_L0,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_y_oneside_L0,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   eigen_indices,
     $                   inter_indices1)

                 case default

                    call error_gradient_type(
     $                   'bf_newgrdpt_class',
     $                   'compute_newgrdpt',
     $                   gradient_type)

               end select

            case(NE_edge_type)
               
               n_direction         = n2_direction
               side                = right
               eigen_indices       = [i1-1,j1-1]
               inter_indices1(:,1) = [i1-1,j1-1]
               inter_indices1(:,2) = [i1  ,j1-1]
               inter_indices1(:,3) = [i1-1,j1  ]

               select case(gradient_type)

                 case(gradient_I_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_interior,
     $                   gradient_y_y_oneside_R0,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_interior,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xLR0_yI_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_y_oneside_R0,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_interior,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xI_yLR0_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_interior,
     $                   gradient_y_y_oneside_R0,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_y_oneside_R0,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xLR0_yLR0_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_y_oneside_R0,
     $                   gradient_x_x_oneside_R0,
     $                   gradient_y_y_oneside_R0,
     $                   eigen_indices,
     $                   inter_indices1)

                 case default

                    call error_gradient_type(
     $                   'bf_newgrdpt_class',
     $                   'compute_newgrdpt',
     $                   gradient_type)

               end select

            case(NW_edge_type)
               
               n_direction         = n1_direction
               side                = left
               eigen_indices       = [i1+1,j1-1]
               inter_indices1(:,1) = [i1  ,j1-1]
               inter_indices1(:,2) = [i1+1,j1-1]
               inter_indices1(:,3) = [i1+1,j1  ]

               select case(gradient_type)

                 case(gradient_I_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_interior,
     $                   gradient_y_y_oneside_R0,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_interior,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xLR0_yI_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_y_oneside_R0,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_interior,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xI_yLR0_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_interior,
     $                   gradient_y_y_oneside_R0,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_y_oneside_R0,
     $                   eigen_indices,
     $                   inter_indices1)

                 case(gradient_xLR0_yLR0_type)

                    new_grdpt = compute_newgrdpt_xy(
     $                   p_model, t, dt,
     $                   bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $                   bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $                   i1,j1,
     $                   n_direction,
     $                   side,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_y_oneside_R0,
     $                   gradient_x_interior,
     $                   gradient_y_interior,
     $                   gradient_x_x_oneside_L0,
     $                   gradient_y_y_oneside_R0,
     $                   eigen_indices,
     $                   inter_indices1)

                 case default
                    
                    call error_gradient_type(
     $                   'bf_newgrdpt_class',
     $                   'compute_newgrdpt',
     $                   gradient_type)

               end select
                    

            case default
               print '(''bf_newgrdpt_class'')'
               print '(''compute_newgrdpt'')'
               print '(''procedure_type not recognized: '',I2)', procedure_type
               stop ''

          end select

        end function compute_newgrdpt_local


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
        function compute_newgrdpt_x(
     $     p_model, t, dt,
     $     bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $     bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $     i1,j1, side_x, gradient_y)
     $     result(new_grdpt)

          implicit none

          type(pmodel_eq)                    , intent(in) :: p_model
          real(rkind)                        , intent(in) :: t
          real(rkind)                        , intent(in) :: dt
          integer(ikind), dimension(2,2)     , intent(in) :: bf_align0
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map0
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map0
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes0
          integer(ikind), dimension(2,2)     , intent(in) :: bf_align1
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map1
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map1
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes1
          integer(ikind)                     , intent(in) :: i1
          integer(ikind)                     , intent(in) :: j1
          logical                            , intent(in) :: side_x
          procedure(gradient_proc)                        :: gradient_y
          real(rkind)   , dimension(ne)                   :: new_grdpt

          integer                       :: k,l
          
          integer                       :: dir, dir2
          integer(ikind)                :: i_eigen
          real(rkind), dimension(ne+1)  :: obc_prim_nodes
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

          real(rkind), dimension(ne)    :: new_grdpt_prim


          !--------------------------------------------------
          !0) determine the direction
          !--------------------------------------------------
          dir  = x_direction
          dir2 = y_direction


          !--------------------------------------------------
          !1) determine the x-coordinate of the new
          !   grid point computed
          !--------------------------------------------------
          x1   = bf_x_map1(i1)
          y1   = bf_y_map1(j1)


          !--------------------------------------------------
          !2) determine where the eigenvalues are
          !   evaluated and which indices are needed
          !   for the interpolation of the grid points
          !--------------------------------------------------
          if(side_x.eqv.right) then

             !x-index for the evaluation of the eigenvalues
             i_eigen = i1-1
             
             !indices for the interpolation of the data at t-dt
             i0_inter1 = bf_align1(dir,1) -bf_align0(dir,1) + i_eigen-1
             i0_inter2 = bf_align1(dir,1) -bf_align0(dir,1) + i_eigen
             j0_inter  = bf_align1(dir2,1)-bf_align0(dir2,1)+ j1
             
             !indices for the interpolation of the data at t
             i1_inter1 = i_eigen-1
             i1_inter2 = i_eigen
             j1_inter  = j1

          else
             
             !x-index for the evaluation of the eigenvalues
             i_eigen = i1+1

             !indices for the interpolation of the data at t-dt
             i0_inter1 = bf_align1(dir,1) -bf_align0(dir,1) +i_eigen
             i0_inter2 = bf_align1(dir,1) -bf_align0(dir,1) +i_eigen+1
             j0_inter  = bf_align1(dir2,1)-bf_align0(dir2,1)+j1
             
             !indices for the interpolation of the data at t
             i1_inter1 = i_eigen
             i1_inter2 = i_eigen+1
             j1_inter  = j1
             
          end if

          dy = bf_y_map0(2)-bf_y_map0(1)
          dy = bf_y_map1(2)-bf_y_map1(1)


          !--------------------------------------------------
          !3) create the interpolation coefficients for
          !   the data at t
          !--------------------------------------------------

          !3.1) create the interpolation coefficients for the nodes
          x_map_inter(1) = bf_x_map0(i0_inter1)
          x_map_inter(2) = bf_x_map0(i0_inter2)
          
          nodes_inter(1,:) = p_model%compute_prim_var(bf_nodes0(i0_inter1,j0_inter,:))
          nodes_inter(2,:) = p_model%compute_prim_var(bf_nodes0(i0_inter2,j0_inter,:))

          inter_nodes0 = get_interpolation_coeff_1D(x_map_inter,nodes_inter)


          !3.2) create the interpolation coefficients for the transverse terms

          !3.2.1) transverse terms at (i0_inter, j0_inter1) at t-dt
          obc_prim_nodes = p_model%get_prim_obc_eigenqties(
     $                        t-dt,x_map_inter(1),y0,
     $                        bf_nodes0(i0_inter1,j0_inter,:))

          nodes_inter(1,:) = MATMUL(
     $         p_model%compute_gradient_prim(bf_nodes0,i0_inter1,j0_inter,gradient_y,dy),
     $         p_model%compute_x_transM_prim(obc_prim_nodes))

          !3.2.2) transverse terms at (i0_inter, j0_inter2) at t-dt
          obc_prim_nodes = p_model%get_prim_obc_eigenqties(
     $                        t-dt,x_map_inter(2),y0,
     $                        bf_nodes0(i0_inter2,j0_inter,:))

          nodes_inter(2,:) = MATMUL(
     $         p_model%compute_gradient_prim(bf_nodes0,i0_inter2,j0_inter,gradient_y,dy),
     $         p_model%compute_x_transM_prim(obc_prim_nodes))

          !3.2.3) interpolation coefficients for the transverse terms at t
          inter_trans0 = get_interpolation_coeff_1D(x_map_inter,nodes_inter)


          !--------------------------------------------------
          !4) create the interpolation coefficients for
          !   the data at t+dt
          !--------------------------------------------------
          x_map_inter(1) = bf_x_map1(i1_inter1)
          x_map_inter(2) = bf_x_map1(i1_inter2)

          !4.1) transverse terms at (i1_inter, j1_inter1) at t
          obc_prim_nodes = p_model%get_prim_obc_eigenqties(
     $                        t,x_map_inter(1),y1,
     $                        bf_nodes1(i1_inter1,j1_inter,:))

          nodes_inter(1,:) = MATMUL(
     $         p_model%compute_gradient_prim(bf_nodes1,i1_inter1,j1_inter,gradient_y,dy),
     $         p_model%compute_x_transM_prim(obc_prim_nodes))

          !4.2) transverse terms at (i1_inter, j1_inter2) at t+dt
          obc_prim_nodes = p_model%get_prim_obc_eigenqties(
     $                         t,x_map_inter(2),y1,
     $                         bf_nodes1(i1_inter2,j1_inter,:))

          nodes_inter(2,:) = MATMUL(
     $         p_model%compute_gradient_prim(bf_nodes1,i1_inter2,j1_inter,gradient_y,dy),
     $         p_model%compute_x_transM_prim(obc_prim_nodes))

          !4.3) interpolation coefficients for the transverse terms at t+dt
          inter_trans1 = get_interpolation_coeff_1D(x_map_inter,nodes_inter)


          !--------------------------------------------------
          !5) interpolate the transverse terms
          !   at (x1,y1)
          !--------------------------------------------------
          t_amp1 = interpolate_1D(x1,inter_trans1)


          !--------------------------------------------------
          !6) determine the nodes for the computation
          !   of the eigenquantities at t+dt
          !--------------------------------------------------
          obc_prim_nodes = p_model%get_prim_obc_eigenqties(
     $                        t,bf_x_map1(i_eigen),y1,
     $                        bf_nodes1(i_eigen,j1,:))


          !--------------------------------------------------
          !7) evaluate the eigenvalues at t+dt
          !--------------------------------------------------
          eigenvalues_x  = p_model%compute_x_eigenvalues_prim(
     $                        obc_prim_nodes)


          !--------------------------------------------------
          !8) determine the left eigenvector
          !   corresponding to the eigenvalue
          !--------------------------------------------------
          left_eigenM    = p_model%compute_x_lefteigenvector_prim(
     $                        obc_prim_nodes)


          !--------------------------------------------------
          !9) determine the characteristic amplitude
          !--------------------------------------------------
          y0 = y1

          do k=1,ne


             !9.1) determine the position where the
             !     characteristic amplitude should
             !     be estimated
             x0 = x1 - eigenvalues_x(k)*dt


             !9.2) determine the normal and transverse
             !     contributions of the hyperbolic
             !     terms to the characteristic amplitude
             if(side_x.eqv.right) then

                if(eigenvalues_x(k).ge.0) then
                   
                   n_amp0 = interpolate_1D(x0,inter_nodes0)
                   t_amp0 = interpolate_1D(x0,inter_trans0)
                   
                else

                   n_amp0 = p_model%compute_prim_var(
     $                         p_model%get_far_field(t-dt,x0,y0))
                   do l=1,ne
                      t_amp0(l) = 0.0d0
                   end do

                end if

             else
                
                if(eigenvalues_x(k).gt.0) then

                   n_amp0 = p_model%compute_prim_var(
     $                         p_model%get_far_field(t-dt,x0,y0))
                   do l=1,ne
                      t_amp0(l) = 0.0d0
                   end do

                else

                   n_amp0 = interpolate_1D(x0,inter_nodes0)
                   t_amp0 = interpolate_1D(x0,inter_trans0)

                end if
             end if


             !9.3) combine the information on the nodes
             !     at t-dt and the approximation of the
             !     integration of the transverse terms
             !     from t-dt to t
             amp =
     $            n_amp0 -
     $            compute_NewtonCotes_integration(t_amp0, t_amp1, dt)

             
             !9.4) compute the scalar product of the
             !     left eigenvector corresponding to
             !     the eigenvalue with the
             !     characteristic amplitude
             char_amp(k) = DOT_PRODUCT(amp,left_eigenM(:,k))
             
          end do


          !--------------------------------------------------
          !10) determine the right eigenmatrix
          !--------------------------------------------------
          right_eigenM = p_model%compute_x_righteigenvector_prim(
     $                      obc_prim_nodes)


          !--------------------------------------------------
          !11) determine the new grid point
          !--------------------------------------------------
          new_grdpt_prim = MATMUL(char_amp,right_eigenM)
          new_grdpt      = p_model%compute_cons_var(new_grdpt_prim)
          
        end function compute_newgrdpt_x


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
        function compute_newgrdpt_y(
     $     p_model, t, dt,
     $     bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $     bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $     i1,j1, side_y, gradient_x)
     $     result(new_grdpt)

          implicit none

          type(pmodel_eq)                    , intent(in) :: p_model
          real(rkind)                        , intent(in) :: t
          real(rkind)                        , intent(in) :: dt
          integer(ikind), dimension(2,2)     , intent(in) :: bf_align0
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map0
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map0
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes0
          integer(ikind), dimension(2,2)     , intent(in) :: bf_align1
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map1
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map1
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes1
          integer(ikind)                     , intent(in) :: i1
          integer(ikind)                     , intent(in) :: j1
          logical                            , intent(in) :: side_y
          procedure(gradient_proc)                        :: gradient_x
          real(rkind)   , dimension(ne)                   :: new_grdpt

          integer                       :: k,l
          
          integer                       :: dir, dir2
          integer(ikind)                :: j_eigen
          real(rkind), dimension(ne+1)  :: obc_prim_nodes
          real(rkind), dimension(ne)    :: eigenvalues_y
          real(rkind), dimension(ne,ne) :: left_eigenM
          real(rkind), dimension(ne,ne) :: right_eigenM
          real(rkind), dimension(ne)    :: n_amp0
          real(rkind), dimension(ne)    :: t_amp0
          real(rkind), dimension(ne)    :: t_amp1
          real(rkind), dimension(ne)    :: amp
          real(rkind), dimension(ne)    :: char_amp
          real(rkind), dimension(ne)    :: new_grdpt_prim

          real(rkind)                   :: dx
          real(rkind)                   :: x0,x1,y0,y1
          integer(ikind)                :: j0_inter1, j0_inter2, i0_inter
          integer(ikind)                :: j1_inter1, j1_inter2, i1_inter
          real(rkind), dimension(2)     :: y_map_inter
          real(rkind), dimension(2,ne)  :: nodes_inter
          real(rkind), dimension(2,ne)  :: inter_nodes0
          real(rkind), dimension(2,ne)  :: inter_trans0
          real(rkind), dimension(2,ne)  :: inter_trans1


          !--------------------------------------------------
          !0) determine the direction
          !--------------------------------------------------
          dir  = y_direction
          dir2 = x_direction


          !--------------------------------------------------
          !1) determine the x-coordinate of the new
          !   grid point computed
          !--------------------------------------------------
          x1   = bf_x_map1(i1)
          y1   = bf_y_map1(j1)


          !--------------------------------------------------
          !2) determine where the eigenvalues are
          !   evaluated and which indices are needed
          !   for the interpolation of the grid points
          !--------------------------------------------------
          if(side_y.eqv.right) then

             !x-index for the evaluation of the eigenvalues
             j_eigen = j1-1
             
             !indices for the interpolation of the data at t-dt
             j0_inter1 = bf_align1(dir,1) -bf_align0(dir,1) + j_eigen-1
             j0_inter2 = bf_align1(dir,1) -bf_align0(dir,1) + j_eigen
             i0_inter  = bf_align1(dir2,1)-bf_align0(dir2,1)+ i1
             
             !indices for the interpolation of the data at t
             j1_inter1 = j_eigen-1
             j1_inter2 = j_eigen
             i1_inter  = i1

          else
             
             !x-index for the evaluation of the eigenvalues
             j_eigen = j1+1

             !indices for the interpolation of the data at t-dt
             j0_inter1 = bf_align1(dir,1) -bf_align0(dir,1) +j_eigen
             j0_inter2 = bf_align1(dir,1) -bf_align0(dir,1) +j_eigen+1
             i0_inter  = bf_align1(dir2,1)-bf_align0(dir2,1)+i1
             
             !indices for the interpolation of the data at t
             j1_inter1 = j_eigen
             j1_inter2 = j_eigen+1
             i1_inter  = i1
             
          end if

          dx = bf_x_map0(2)-bf_x_map0(1)
          dx = bf_x_map1(2)-bf_x_map1(1)


          !--------------------------------------------------
          !3) create the interpolation coefficients for
          !   the data at t
          !--------------------------------------------------

          !3.1) create the interpolation coefficients for the nodes
          y_map_inter(1) = bf_y_map0(j0_inter1)
          y_map_inter(2) = bf_y_map0(j0_inter2)
          
          nodes_inter(1,:) = p_model%compute_prim_var(bf_nodes0(i0_inter,j0_inter1,:))
          nodes_inter(2,:) = p_model%compute_prim_var(bf_nodes0(i0_inter,j0_inter2,:))

          inter_nodes0 = get_interpolation_coeff_1D(y_map_inter,nodes_inter)


          !3.2) create the interpolation coefficients for the transverse terms

          !3.2.1) transverse terms at (i0_inter, j0_inter1) at t-dt
          obc_prim_nodes = p_model%get_prim_obc_eigenqties(
     $                        t-dt,x0,y_map_inter(1),
     $                        bf_nodes0(i0_inter,j0_inter1,:))

          nodes_inter(1,:) = MATMUL(
     $         p_model%compute_gradient_prim(bf_nodes0,i0_inter,j0_inter1,gradient_x,dx),
     $         p_model%compute_y_transM_prim(obc_prim_nodes))


          !3.2.2) transverse terms at (i0_inter, j0_inter2) at t-dt
          obc_prim_nodes = p_model%get_prim_obc_eigenqties(
     $                        t-dt,x0,y_map_inter(2),
     $                        bf_nodes0(i0_inter,j0_inter2,:))

          nodes_inter(2,:) = MATMUL(
     $         p_model%compute_gradient_prim(bf_nodes0,i0_inter,j0_inter2,gradient_x,dx),
     $         p_model%compute_y_transM_prim(obc_prim_nodes))

          !3.2.3) interpolation coefficients for the transverse terms at t-dt
          inter_trans0 = get_interpolation_coeff_1D(y_map_inter,nodes_inter)


          !--------------------------------------------------
          !4) create the interpolation coefficients for
          !   the data at t
          !--------------------------------------------------
          y_map_inter(1) = bf_y_map1(j1_inter1)
          y_map_inter(2) = bf_y_map1(j1_inter2)

          !4.1) transverse terms at (i1_inter, j1_inter1) at t
          obc_prim_nodes = p_model%get_prim_obc_eigenqties(
     $                        t,x1,y_map_inter(1),
     $                        bf_nodes1(i1_inter,j1_inter1,:))

          nodes_inter(1,:) = MATMUL(
     $         p_model%compute_gradient_prim(bf_nodes1,i1_inter,j1_inter1,gradient_x,dx),
     $         p_model%compute_y_transM_prim(obc_prim_nodes))

          !4.2) transverse terms at (i1_inter, j1_inter2) at t
          obc_prim_nodes = p_model%get_prim_obc_eigenqties(
     $                        t,x1,y_map_inter(2),
     $                        bf_nodes1(i1_inter,j1_inter2,:))

          nodes_inter(2,:) = MATMUL(
     $         p_model%compute_gradient_prim(bf_nodes1,i1_inter,j1_inter2,gradient_x,dx),
     $         p_model%compute_y_transM_prim(obc_prim_nodes))

          !4.3) interpolation coefficients for the transverse terms at t
          inter_trans1 = get_interpolation_coeff_1D(y_map_inter,nodes_inter)


          !--------------------------------------------------
          !5) interpolate the transverse terms
          !   at (x1,y1)
          !--------------------------------------------------
          t_amp1 = interpolate_1D(y1,inter_trans1)


          !--------------------------------------------------
          !6) determine the nodes for the computation
          !   of the eigenquantities at t
          !--------------------------------------------------
          obc_prim_nodes = p_model%get_prim_obc_eigenqties(
     $                        t,x1,bf_y_map1(j_eigen),
     $                        bf_nodes1(i1,j_eigen,:))


          !--------------------------------------------------
          !7) evaluate the eigenvalues at t
          !--------------------------------------------------
          eigenvalues_y  = p_model%compute_y_eigenvalues_prim(
     $                        obc_prim_nodes)


          !--------------------------------------------------
          !8) determine the left eigenvector
          !   corresponding to the eigenvalue
          !--------------------------------------------------
          left_eigenM    = p_model%compute_y_lefteigenvector_prim(
     $                        obc_prim_nodes)
             

          !--------------------------------------------------
          !9) determine the characteristic amplitude
          !--------------------------------------------------
          x0 = x1

          do k=1,ne


             !9.1) determine the position where the
             !     characteristic amplitude should
             !     be estimated
             y0 = y1 - eigenvalues_y(k)*dt


             !9.2) determine the normal and transverse
             !     contributions of the hyperbolic
             !     terms to the characteristic amplitude
             if(side_y.eqv.right) then

                if(eigenvalues_y(k).ge.0) then
                   
                   n_amp0 = interpolate_1D(y0,inter_nodes0)
                   t_amp0 = interpolate_1D(y0,inter_trans0)
                   
                else

                   n_amp0 = p_model%compute_prim_var(
     $                         p_model%get_far_field(t-dt,x0,y0))
                   do l=1,ne
                      t_amp0(l) = 0.0d0
                   end do

                end if

             else
                
                if(eigenvalues_y(k).gt.0) then

                   n_amp0 = p_model%compute_prim_var(
     $                         p_model%get_far_field(t-dt,x0,y0))
                   do l=1,ne
                      t_amp0(l) = 0.0d0
                   end do

                else

                   n_amp0 = interpolate_1D(y0,inter_nodes0)
                   t_amp0 = interpolate_1D(y0,inter_trans0)

                end if
             end if


             !9.3) combine the information on the nodes
             !     at t-dt and the approximation of the
             !     integration of the transverse terms
             !     from t-dt to t
             amp =
     $            n_amp0 -
     $            compute_NewtonCotes_integration(t_amp0, t_amp1, dt)

             
             !9.4) compute the scalar product of the
             !     left eigenvector corresponding to
             !     the eigenvalue with the
             !     characteristic amplitude
             char_amp(k) = DOT_PRODUCT(amp,left_eigenM(:,k))
             
          end do


          !--------------------------------------------------
          !10) determine the right eigenmatrix
          !--------------------------------------------------
          right_eigenM = p_model%compute_y_righteigenvector_prim(
     $                      obc_prim_nodes)


          !--------------------------------------------------
          !11) determine the new grid point
          !--------------------------------------------------
          new_grdpt_prim = MATMUL(char_amp,right_eigenM)
          new_grdpt      = p_model%compute_cons_var(new_grdpt_prim)

        end function compute_newgrdpt_y


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
        !>@param j1
        !> y-index identifying the new grdpt at t=t
        !             
        !>@param n_direction
        !> integer identifying whether the vector normal to the
        !> boundary is located along the n1- or n2- direction
        !
        !>@param side_n
        !> logical determining in which direction the information
        !> is incoming
        !             
        !>@param gradient_n_index1
        !> procedure for computing the gradient at the grid point
        !> of index1
        !
        !>@param gradient_n_index2
        !> procedure for computing the gradient at the grid point
        !> of index2
        !
        !>@param gradient_n_index3
        !> procedure for computing the gradient at the grid point
        !> of index3
        !
        !>@param eigen_indices
        !> index where the eigen quantities  (eigenvalues, left
        !> eigenvector, right eigenvector) are evaluated
        !
        !>@param inter_indices
        !> indices identifying the interpolation points
        !--------------------------------------------------------------
        function compute_newgrdpt_xy(
     $     p_model, t, dt,
     $     bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $     bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $     i1,j1,
     $     n_direction,
     $     side_n,
     $     gradient_x_index1,
     $     gradient_y_index1,
     $     gradient_x_index2,
     $     gradient_y_index2,
     $     gradient_x_index3,
     $     gradient_y_index3,
     $     eigen_indices,
     $     inter_indices1)
     $     result(new_grdpt)

          implicit none

          type(pmodel_eq)                    , intent(in) :: p_model
          real(rkind)                        , intent(in) :: t
          real(rkind)                        , intent(in) :: dt
          integer(ikind), dimension(2,2)     , intent(in) :: bf_align0
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map0
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map0
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes0
          integer(ikind), dimension(2,2)     , intent(in) :: bf_align1
          real(rkind)   , dimension(:)       , intent(in) :: bf_x_map1
          real(rkind)   , dimension(:)       , intent(in) :: bf_y_map1
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes1
          integer(ikind)                     , intent(in) :: i1
          integer(ikind)                     , intent(in) :: j1
          integer                            , intent(in) :: n_direction
          logical                            , intent(in) :: side_n
          procedure(gradient_proc)                        :: gradient_x_index1
          procedure(gradient_proc)                        :: gradient_y_index1
          procedure(gradient_proc)                        :: gradient_x_index2
          procedure(gradient_proc)                        :: gradient_y_index2
          procedure(gradient_proc)                        :: gradient_x_index3
          procedure(gradient_proc)                        :: gradient_y_index3
          integer(ikind), dimension(2)       , intent(in) :: eigen_indices
          integer(ikind), dimension(2,3)     , intent(in) :: inter_indices1
          real(rkind)   , dimension(ne)                   :: new_grdpt

          
          !x1,y1
          !cartesian coordinates of the new grid point computed
          !
          !x0,y0
          !cartesian coordinates of the grid point where the
          !characteristic amplitude can be evaluated at t=t-dt
          !
          !n1_1, n2_1
          !(n1,n2)-coordinates of the new grid point computed
          !at t=t
          !
          !n1_0, n2_0
          !(n1,n2)-coordinates of the grid point where the
          !characteristic amplitude can be evaluated at t=t-dt
          !----------------------------------------------------
          real(rkind)                   :: x1,y1
          real(rkind)                   :: x0,y0
          real(rkind)                   :: n1_1,n2_1
          real(rkind)                   :: n1_0,n2_0
          integer                       :: k,l

       
          !nodes_eigenqties
          !grdpts used to evaluate the eigenquantities
          !
          !eigenvalues_n
          !eigenvalues in the direction dir
          !
          !left_eigenM
          !left eigenmatrix in the direction dir
          !
          !right_eigenM
          !right eigenmatrix in the direction dir
          !
          !n_amp0
          !normal contribution to the characteristic amplitude
          !at t=t-dt
          !
          !t_amp0
          !transverse contribution to the characteristic amplitude
          !at t=t-dt
          !
          !t_amp1
          !transverse contribution to the characteristic amplitude
          !at t=t
          !
          !amp
          !amplitude computed as the sum of the normal and the time
          !integration of the transverse contributions
          !
          !char_amp
          !characteristic wave which is assumed constant in time
          !----------------------------------------------------
          real(rkind), dimension(ne+1)  :: obc_prim_nodes
          real(rkind), dimension(ne)    :: eigenvalues_n
          real(rkind), dimension(ne,ne) :: left_eigenM
          real(rkind), dimension(ne,ne) :: right_eigenM
          real(rkind), dimension(ne)    :: n_amp0
          real(rkind), dimension(ne)    :: t_amp0
          real(rkind), dimension(ne)    :: t_amp1
          real(rkind), dimension(ne)    :: amp
          real(rkind), dimension(ne)    :: char_amp
          real(rkind), dimension(ne)    :: new_grdpt_prim_n
          real(rkind), dimension(ne)    :: new_grdpt_prim


          !inter_indices_0
          !integer identifying the interpolation points for the
          !computation of the new grid points on the grdpts_id
          !array at t=t-dt
          !
          !i_x_inter:
          !index identifying the x-coordinate of the interpolation
          !point
          !
          !i_y_inter:
          !index identifying the y-coordinate of the interpolation
          !point
          !
          !x_inter:
          !x-coordinates of the interpolation point
          !
          !y_inter:
          !y-coordinates of the interpolation point
          !
          !n1_inter:
          !intermediate array where the n1-coordinates of
          !the interpolation points are saved
          !
          !n2_inter:
          !intermediate array where the n2-coordinates of
          !the interpolation points are saved
          !
          !nodes_inter:
          !intermediate array where the data of the
          !interpolation points are saved
          !
          !inter_nodes0:
          !interpolation coefficient for the plane by the
          !nodes at t-dt
          !
          !inter_trans0:
          !interpolation coefficient for the plane by the
          !transverse terms at t-dt
          !
          !inter_trans1:
          !interpolation coefficient for the plane by the
          !transverse terms at t
          !------------------------------------------------
          integer(ikind), dimension(2,3)  :: inter_indices0
          integer(ikind)                  :: i_x_inter
          integer(ikind)                  :: i_y_inter
          real(rkind)                     :: x_inter
          real(rkind)                     :: y_inter
          real(rkind)   , dimension(3)    :: n1_inter     
          real(rkind)   , dimension(3)    :: n2_inter     
          real(rkind)   , dimension(3,ne) :: nodes_inter  
          real(rkind)   , dimension(3,ne) :: inter_nodes0
          real(rkind)   , dimension(3,ne) :: inter_trans0
          real(rkind)   , dimension(3,ne) :: inter_trans1
          real(rkind)                     :: dx,dy


          !--------------------------------------------------
          !1) determine the (x,y)-coordinates
          !   of the new grid point computed
          !--------------------------------------------------
          x1 = bf_x_map1(i1)
          y1 = bf_y_map1(j1)

          dx = bf_x_map1(2)-bf_x_map1(1)
          dy = bf_y_map1(2)-bf_y_map1(1)

          !--------------------------------------------------
          !2) convert them into (n1,n2)
          !   coordinates
          !--------------------------------------------------
          n1_1 = get_n1_coord(x1,y1)
          n2_1 = get_n2_coord(x1,y1)


          !--------------------------------------------------
          !3) determine where the eigenvalues
          !   are evaluated and which indices
          !   are needed for the interpolation
          !   of the grid points
          !   i.e. convert the interpolation
          !   indices at t into interpolation
          !   indices at t-dt
          !--------------------------------------------------
          do k=1, 3
             
             inter_indices0(1,k) = bf_align1(1,1) - bf_align0(1,1) + inter_indices1(1,k)
             inter_indices0(2,k) = bf_align1(2,1) - bf_align0(2,1) + inter_indices1(2,k)

          end do


          !--------------------------------------------------
          !4) create the interpolation
          !   coefficients for the data at
          !   t-dt
          !
          !   the interpolation data are the
          !   primitive variables converted
          !   into primitive variables in the
          !   (n1,n2) reference frame
          !--------------------------------------------------

          !4.1) extract the interpolation
          !     data for the nodes at t-dt
          do k=1,3
             
             i_x_inter        = inter_indices0(1,k)
             i_y_inter        = inter_indices0(2,k)

             x_inter          = bf_x_map0(i_x_inter)
             y_inter          = bf_y_map0(i_y_inter)

             n1_inter(k)      = get_n1_coord(x_inter,y_inter)
             n2_inter(k)      = get_n2_coord(x_inter,y_inter)

             nodes_inter(k,:) = p_model%compute_xy_to_n_var(
     $                             p_model%compute_prim_var(
     $                                bf_nodes0(i_x_inter,i_y_inter,:)))

          end do

          !4.2) create the interpolation
          !     coefficients for the nodes
          !     at t-dt
          inter_nodes0 = get_interpolation_coeff_2D(
     $         n1_inter,
     $         n2_inter,
     $         nodes_inter)


          !4.3) create the interpolation
          !     coefficients for the
          !     transverse terms at t-dt

          !4.3.1) create the data at the
          !       interpolation grid points
          call get_n_transverse_data_for_interpolation(
     $         t-dt,
     $         p_model,
     $         bf_x_map0,
     $         bf_y_map0,
     $         bf_nodes0,
     $         inter_indices0,
     $         n_direction,
     $         gradient_x_index1,
     $         gradient_y_index1,
     $         gradient_x_index2,
     $         gradient_y_index2,
     $         gradient_x_index3,
     $         gradient_y_index3,
     $         dx,dy,
     $         nodes_inter)
          
          !4.3.2) create the interpolation
          !       plane for the contribution
          !       of the transverse term at t-dt
          inter_trans0 = get_interpolation_coeff_2D(
     $         n1_inter,
     $         n2_inter,
     $         nodes_inter)


          !--------------------------------------------------
          !5) create the interpolation coefficients
          !   for the data at t+dt
          !
          !   the interpolation data are the
          !   primitive variables converted
          !   into primitive variables in the
          !   (n1,n2) reference frame
          !--------------------------------------------------
          
          !5.1) create the coordinate maps (n1,n2)
          !     identifying the position of the
          !     interpolation points
          do k=1,3
             
             i_x_inter        = inter_indices1(1,k)
             i_y_inter        = inter_indices1(2,k)

             x_inter          = bf_x_map1(i_x_inter)
             y_inter          = bf_y_map1(i_y_inter)

             n1_inter(k)      = get_n1_coord(x_inter,y_inter)
             n2_inter(k)      = get_n2_coord(x_inter,y_inter)

          end do

          !5.2) compute the transverse terms at
          !     the interpolation points
          call get_n_transverse_data_for_interpolation(
     $         t,
     $         p_model,
     $         bf_x_map1,
     $         bf_y_map1,
     $         bf_nodes1,
     $         inter_indices1,
     $         n_direction,
     $         gradient_x_index1,
     $         gradient_y_index1,
     $         gradient_x_index2,
     $         gradient_y_index2,
     $         gradient_x_index3,
     $         gradient_y_index3,
     $         dx,dy,
     $         nodes_inter)

          !5.3) create the interpolation plane
          !     for the contribution of the
          !     transverse term at t
          inter_trans1 = get_interpolation_coeff_2D(
     $         n1_inter,
     $         n2_inter,
     $         nodes_inter)


          !--------------------------------------------------
          !6) interpolate the transverse terms
          !   at t
          !--------------------------------------------------
          t_amp1 = interpolate_2D(n1_1,n2_1,inter_trans1)


          !--------------------------------------------------
          !7) determine the nodes for the
          !   computation of the eigenquantities
          !   at t
          !--------------------------------------------------
          obc_prim_nodes =  p_model%get_prim_obc_eigenqties(
     $                        t,x1,y1,
     $                        bf_nodes1(eigen_indices(1),
     $                                   eigen_indices(2),
     $                                   :))

          obc_prim_nodes(1:ne) = p_model%compute_xy_to_n_var(
     $                              obc_prim_nodes(1:ne))


          !--------------------------------------------------
          !8) evaluate the eigenquantities at t
          !--------------------------------------------------
          select case(n_direction)

            case(n1_direction)

               eigenvalues_n = p_model%compute_x_eigenvalues_prim(obc_prim_nodes)
               left_eigenM   = p_model%compute_x_lefteigenvector_prim(obc_prim_nodes)
               right_eigenM  = p_model%compute_x_righteigenvector_prim(obc_prim_nodes)

            case(n2_direction)

               eigenvalues_n = p_model%compute_y_eigenvalues_prim(obc_prim_nodes)
               left_eigenM   = p_model%compute_y_lefteigenvector_prim(obc_prim_nodes)
               right_eigenM  = p_model%compute_y_righteigenvector_prim(obc_prim_nodes)

            case default
               print '(''bf_newgrdpt_class'')'
               print '(''compute_newgrdpt_xy'')'
               print '(''direction not recognized: '', I2)', n_direction
               stop ''

          end select

          
          !--------------------------------------------------
          !9) determine the characteristic amplitude
          !--------------------------------------------------
           do k=1,ne
              
             !9.1) determine the position where the
             !     characteristic amplitude should
             !     be estimated
             select case(n_direction)
             
               case(n1_direction)
                  n1_0 = n1_1 - eigenvalues_n(k)*dt
                  n2_0 = n2_1
                  
               case(n2_direction)
                  n1_0 = n1_1
                  n2_0 = n2_1 - eigenvalues_n(k)*dt

               case default
                  print '(''bf_newgrdpt_class'')'
                  print '(''compute_newgrdpt_xy'')'
                  print '(''direction not recognized: '',I2)', n_direction
                  stop ''

             end select


             !9.2) determine the normal and
             !     transverse contributions of the
             !     hyperbolic terms to the
             !     characteristic amplitude
             if(side_n.eqv.right) then

                if(eigenvalues_n(k).ge.0) then
               
                   n_amp0 = interpolate_2D(n1_0, n2_0, inter_nodes0)
                   t_amp0 = interpolate_2D(n1_0, n2_0, inter_trans0)
                   
                else

                   x0 = get_x_coord(n1_0,n2_0)
                   y0 = get_y_coord(n1_0,n2_0)

                   n_amp0 = p_model%compute_xy_to_n_var(
     $                         p_model%compute_prim_var(
     $                            p_model%get_far_field(t,x0,y0)))

                   do l=1,ne
                      t_amp0(l) = 0.0d0
                   end do
                   
                end if
                
             else
                
                if(eigenvalues_n(k).gt.0) then

                   x0 = get_x_coord(n1_0,n2_0)
                   y0 = get_y_coord(n1_0,n2_0)

                   n_amp0 = p_model%compute_xy_to_n_var(
     $                         p_model%compute_prim_var(
     $                            p_model%get_far_field(t,x0,y0)))

                   do l=1,ne
                      t_amp0(l) = 0.0d0
                   end do

                else

                   n_amp0 = interpolate_2D(n1_0, n2_0, inter_nodes0)
                   t_amp0 = interpolate_2D(n1_0, n2_0, inter_trans0)

                end if
             end if


             !9.3) combine the information on the
             !     nodes at t-dt and the approximation
             !     of the integration of the transverse
             !     terms from t-dt to t
             amp =
     $            n_amp0 -
     $            compute_NewtonCotes_integration(t_amp0, t_amp1, dt)

             
             !9.4) compute the scalar product of the
             !     left eigenvector corresponding to
             !     the eigenvalue with the characteristic
             !     amplitude
             char_amp(k) = DOT_PRODUCT(amp,left_eigenM(:,k))
             
          end do


          !--------------------------------------------------
          !10) determine the new grid point
          !--------------------------------------------------
          new_grdpt_prim_n = MATMUL(char_amp,right_eigenM)
          new_grdpt_prim   = p_model%compute_n_to_xy_var(new_grdpt_prim_n)
          new_grdpt        = p_model%compute_cons_var(new_grdpt_prim)

        end function compute_newgrdpt_xy


        !compute the contribution of the transverse terms
        !at the location of the interpolation points
        subroutine get_n_transverse_data_for_interpolation(
     $     t,
     $     p_model,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     inter_indices,
     $     n_direction,
     $     gradient_x_index1,
     $     gradient_y_index1,
     $     gradient_x_index2,
     $     gradient_y_index2,
     $     gradient_x_index3,
     $     gradient_y_index3,
     $     dx,dy,
     $     nodes_inter)

          implicit none

          real(rkind)                      , intent(in)  :: t
          type(pmodel_eq)                  , intent(in)  :: p_model
          real(rkind)    , dimension(:)    , intent(in)  :: bf_x_map
          real(rkind)    , dimension(:)    , intent(in)  :: bf_y_map
          real(rkind)    , dimension(:,:,:), intent(in)  :: bf_nodes
          integer(ikind) , dimension(2,3)  , intent(in)  :: inter_indices
          integer                          , intent(in)  :: n_direction
          procedure(gradient_proc)                       :: gradient_x_index1
          procedure(gradient_proc)                       :: gradient_y_index1
          procedure(gradient_proc)                       :: gradient_x_index2
          procedure(gradient_proc)                       :: gradient_y_index2
          procedure(gradient_proc)                       :: gradient_x_index3
          procedure(gradient_proc)                       :: gradient_y_index3
          real(rkind)                      , intent(in)  :: dx
          real(rkind)                      , intent(in)  :: dy
          real(rkind)    , dimension(3,ne) , intent(out) :: nodes_inter


          integer                       :: k,l
          real(rkind), dimension(ne)    :: x_gradient      ! x-gradient evaluated with the primitive var
          real(rkind), dimension(ne)    :: y_gradient      ! y-gradient evaluated with the primitive var
          real(rkind), dimension(ne)    :: n_gradient_prim ! n-gradient evaluated with the primitive var
          real(rkind), dimension(ne)    :: n_gradient      ! n-gradient evaluated with the (n1,n2)-primitive var
          real(rkind), dimension(ne,ne) :: n_transM
          real(rkind), dimension(ne+1)  :: obc_prim_nodes


          !----------------------------------------
          ! compute the transverse terms at the
          ! grid point 1
          !----------------------------------------
          do k=1,3

             select case(k)

               case(1)
                  x_gradient = p_model%compute_gradient_prim(
     $                 bf_nodes,
     $                 inter_indices(1,k),
     $                 inter_indices(2,k),
     $                 gradient_x_index1,
     $                 dx)

                  y_gradient = p_model%compute_gradient_prim(
     $                 bf_nodes,
     $                 inter_indices(1,k),
     $                 inter_indices(2,k),
     $                 gradient_y_index1,
     $                 dy)

               case(2)
                  x_gradient = p_model%compute_gradient_prim(
     $                 bf_nodes,
     $                 inter_indices(1,k),
     $                 inter_indices(2,k),
     $                 gradient_x_index2,
     $                 dx)

                  y_gradient = p_model%compute_gradient_prim(
     $                 bf_nodes,
     $                 inter_indices(1,k),
     $                 inter_indices(2,k),
     $                 gradient_y_index2,
     $                 dy)
                  
               case(3)
                  x_gradient = p_model%compute_gradient_prim(
     $                 bf_nodes,
     $                 inter_indices(1,k),
     $                 inter_indices(2,k),
     $                 gradient_x_index3,
     $                 dx)

                  y_gradient = p_model%compute_gradient_prim(
     $                 bf_nodes,
     $                 inter_indices(1,k),
     $                 inter_indices(2,k),
     $                 gradient_y_index3,
     $                 dy)

             end select

             select case(n_direction)

               !the eigenvalues are evaluated in the n1-direction
               !so, the transverse terms (i.e. the gradient) is
               !evaluated in the n2-direction
               case(n1_direction)

                  do l=1,ne
                     n_gradient_prim(l) = get_n2_coord(x_gradient(l),y_gradient(l))
                  end do

               
               !the eigenvalues are evaluated in the n2-direction
               !so, the transverse terms (i.e. the gradient) is
               !evaluated in the n1-direction
               case(n2_direction)
                  
                  do l=1,ne
                     n_gradient_prim(l) = get_n1_coord(x_gradient(l),y_gradient(l))
                  end do

             end select

             n_gradient = p_model%compute_xy_to_n_var(n_gradient_prim)

             obc_prim_nodes = p_model%get_prim_obc_eigenqties(
     $                           t,
     $                           bf_x_map(inter_indices(1,k)),
     $                           bf_y_map(inter_indices(2,k)),
     $                           bf_nodes(inter_indices(1,k),
     $                                    inter_indices(2,k),
     $                                    :)
     $                           )

             obc_prim_nodes(1:ne) = p_model%compute_xy_to_n_var(
     $                                 obc_prim_nodes(1:ne))


             select case(n_direction)
               case(n1_direction)
                  n_transM = p_model%compute_x_transM_prim(obc_prim_nodes)
               case(n2_direction)
                  n_transM = p_model%compute_y_transM_prim(obc_prim_nodes)
               case default
                  print '(''compute_newgrdpt_xy'')'
                  print '(''direction not recognized: '',I2)', n_direction
                  stop ''
             end select
          
             nodes_inter(k,:) = MATMUL(n_gradient,n_transM)

          end do

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

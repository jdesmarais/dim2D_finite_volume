      !> @file
      !> module implementing subroutines to compute the gradients,
      !> the eigenvalues and the eigenvectors in the (x-y) and (x+y)
      !> directions for the dim2d governing equations
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing subroutines to compute the gradients,
      !> the eigenvalues and the eigenvectors in the (x-y) and (x+y)
      !> directions for the dim2d governing equations
      !
      !> @date
      !> 10_12_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_ncoords_module

        use dim2d_prim_module, only :
     $       compute_jacobian_prim_to_cons,
     $       compute_jacobian_cons_to_prim,
     $       speed_of_sound

        use parameters_input, only :
     $       ne
        
        use parameters_kind, only :
     $       ikind,
     $       rkind

        private
        public ::
     $       compute_n1_eigenvalues_dim2d,
     $       compute_n2_eigenvalues_dim2d,
     $       compute_n1_lefteigenvector_dim2d,
     $       compute_n1_righteigenvector_dim2d,
     $       compute_n2_lefteigenvector_dim2d,
     $       compute_n2_righteigenvector_dim2d,
     $       compute_n1_transM_dim2d,
     $       compute_n2_transM_dim2d


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the (x-y)-direction
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_n1_eigenvalues_dim2d(nodes) result(eigenvalues)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne)             :: eigenvalues

          real(rkind) :: ux
          real(rkind) :: uy
          real(rkind) :: c
          real(rkind) :: ratio
          real(rkind) :: u_av

          ux = nodes(2)/nodes(1)
          uy = nodes(3)/nodes(1)
          c  = speed_of_sound(nodes)

          if(rkind.eq.8) then
             ratio = 0.5d0*Sqrt(2.0d0)
          else
             ratio = 0.5*Sqrt(2.0)
          end if

          u_av  = ratio*(ux-uy)

          eigenvalues(1) =  u_av
          eigenvalues(2) =  u_av
          eigenvalues(3) =  u_av - c
          eigenvalues(4) =  u_av + c

        end function compute_n1_eigenvalues_dim2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the (x+y)-direction
        !
        !> @date
        !> 10_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_n2_eigenvalues_dim2d(nodes) result(eigenvalues)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne)             :: eigenvalues

          real(rkind) :: ux
          real(rkind) :: uy
          real(rkind) :: c
          real(rkind) :: ratio
          real(rkind) :: u_av

          ux = nodes(2)/nodes(1)
          uy = nodes(3)/nodes(1)
          c  = speed_of_sound(nodes)

          if(rkind.eq.8) then
             ratio = 0.5d0*Sqrt(2.0d0)
          else
             ratio = 0.5*Sqrt(2.0)
          end if

          u_av  = ratio*(ux+uy)

          eigenvalues(1) =  u_av
          eigenvalues(2) =  u_av
          eigenvalues(3) =  u_av - c
          eigenvalues(4) =  u_av + c

        end function compute_n2_eigenvalues_dim2d

      end module dim2d_ncoords_module

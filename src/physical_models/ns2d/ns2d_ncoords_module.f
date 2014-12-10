      !> @file
      !> module implementing subroutines to compute the gradients,
      !> the eigenvalues and the eigenvectors in the (x-y) and (x+y)
      !> directions for the ns2d governing equations
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing subroutines to compute the gradients,
      !> the eigenvalues and the eigenvectors in the (x-y) and (x+y)
      !> directions for the ns2d governing equations
      !
      !> @date
      !> 02_12_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
       module ns2d_ncoords_module

        use ns2d_prim_module, only :
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
     $       compute_n1_eigenvalues_ns2d,
     $       compute_n2_eigenvalues_ns2d,
     $       compute_n1_lefteigenvector_ns2d,
     $       compute_n1_righteigenvector_ns2d,
     $       compute_n2_lefteigenvector_ns2d,
     $       compute_n2_righteigenvector_ns2d,
     $       compute_n1_transM_ns2d,
     $       compute_n2_transM_ns2d


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the (x-y)-direction
        !
        !> @date
        !> 02_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_n1_eigenvalues_ns2d(nodes) result(eigenvalues)

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
          eigenvalues(3) =  -c + u_av
          eigenvalues(4) =   c + u_av

        end function compute_n1_eigenvalues_ns2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the (x+y)-direction
        !
        !> @date
        !> 02_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_n2_eigenvalues_ns2d(nodes) result(eigenvalues)

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
          eigenvalues(3) =  -c + u_av
          eigenvalues(4) =   c + u_av

        end function compute_n2_eigenvalues_ns2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenvector for the hyperbolic terms
        !> in the (x-y)-direction. By denoting L the left eigenmatrix, the
        !> result of the function is L[k,:]
        !
        !> @date
        !> 02_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_n1_lefteigenvector_ns2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind), dimension(ne,ne) :: jacPrimCons
          real(rkind), dimension(ne,ne) :: leftEigenMPrim
          real(rkind) :: c
          real(rkind) :: c2_inv
          real(rkind) :: ratio


          !computation of J, the jacobian matrix from primitive
          !to conservative variables
          jacPrimCons = compute_jacobian_prim_to_cons(nodes)


          !left eigenmatrix for the primitive
          !variables, L_p
          c = speed_of_sound(nodes)

          if(rkind.eq.8) then

             c2_inv = 1.0d0/c**2
             ratio = 0.25d0*nodes(1)*c*Sqrt(2.0d0)

             leftEigenMPrim(1,1) =  0.0d0
             leftEigenMPrim(2,1) =  0.5d0
             leftEigenMPrim(3,1) =  0.5d0
             leftEigenMPrim(4,1) =  0.0d0

             leftEigenMPrim(1,2) =  1.0d0
             leftEigenMPrim(2,2) =  0.0d0
             leftEigenMPrim(3,2) =  0.0d0
             leftEigenMPrim(4,2) = -c2_inv

             leftEigenMPrim(1,3) =  0.0d0
             leftEigenMPrim(2,3) = -ratio
             leftEigenMPrim(3,3) =  ratio
             leftEigenMPrim(4,3) =  0.5d0

             leftEigenMPrim(1,4) =  0.0d0
             leftEigenMPrim(2,4) =  ratio
             leftEigenMPrim(3,4) = -ratio
             leftEigenMPrim(4,4) =  0.5d0

          else

             c2_inv = 1.0d0/c**2
             ratio = 0.25d0*nodes(1)*c*Sqrt(2.0d0)

             leftEigenMPrim(1,1) =  0.0
             leftEigenMPrim(2,1) =  0.5
             leftEigenMPrim(3,1) =  0.5
             leftEigenMPrim(4,1) =  0.0

             leftEigenMPrim(1,2) =  1.0
             leftEigenMPrim(2,2) =  0.0
             leftEigenMPrim(3,2) =  0.0
             leftEigenMPrim(4,2) = -c2_inv

             leftEigenMPrim(1,3) =  0.0
             leftEigenMPrim(2,3) = -ratio
             leftEigenMPrim(3,3) =  ratio
             leftEigenMPrim(4,3) =  0.5

             leftEigenMPrim(1,4) =  0.0
             leftEigenMPrim(2,4) =  ratio
             leftEigenMPrim(3,4) = -ratio
             leftEigenMPrim(4,4) =  0.5

          end if

          eigenvect = MATMUL(jacPrimCons,leftEigenMPrim)          

        end function compute_n1_lefteigenvector_ns2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the right eigenvector for the hyperbolic terms
        !> in the (x-y)-direction. By denoting R the right eigenmatrix, the
        !> result of the function is R[k,:]
        !
        !> @date
        !> 02_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvector at the location of the grid point
        !--------------------------------------------------------------
        function compute_n1_righteigenvector_ns2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect

          real(rkind), dimension(ne,ne) :: jacConsPrim
          real(rkind), dimension(ne,ne) :: rightEigenMPrim
          real(rkind) :: c
          real(rkind) :: c2_inv
          real(rkind) :: ratio


          !computation of J, the jacobian matrix from primitive
          !to conservative variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)


          !right eigenmatrix for the primitive
          !variables, L_p
          c = speed_of_sound(nodes)

          if(rkind.eq.8) then

             c2_inv = 1.0d0/c**2
             ratio = 0.5d0*Sqrt(2.0d0)/(nodes(1)*c)

             rightEigenMPrim(1,1) =  0.0d0
             rightEigenMPrim(2,1) =  1.0d0
             rightEigenMPrim(3,1) =  c2_inv
             rightEigenMPrim(4,1) =  c2_inv

             rightEigenMPrim(1,2) =  1.0d0
             rightEigenMPrim(2,2) =  0.0d0
             rightEigenMPrim(3,2) = -ratio
             rightEigenMPrim(4,2) =  ratio

             rightEigenMPrim(1,3) =  1.0d0
             rightEigenMPrim(2,3) =  0.0d0
             rightEigenMPrim(3,3) =  ratio
             rightEigenMPrim(4,3) = -ratio

             rightEigenMPrim(1,4) =  0.0d0
             rightEigenMPrim(2,4) =  0.0d0
             rightEigenMPrim(3,4) =  1.0d0
             rightEigenMPrim(4,4) =  1.0d0

          else

             c2_inv = 1.0/c**2
             ratio = 0.5*Sqrt(2.0)/(nodes(1)*c)

             rightEigenMPrim(1,1) =  0.0
             rightEigenMPrim(2,1) =  1.0
             rightEigenMPrim(3,1) =  c2_inv
             rightEigenMPrim(4,1) =  c2_inv

             rightEigenMPrim(1,2) =  1.0
             rightEigenMPrim(2,2) =  0.0
             rightEigenMPrim(3,2) = -ratio
             rightEigenMPrim(4,2) =  ratio

             rightEigenMPrim(1,3) =  1.0
             rightEigenMPrim(2,3) =  0.0
             rightEigenMPrim(3,3) =  ratio
             rightEigenMPrim(4,3) = -ratio

             rightEigenMPrim(1,4) =  0.0
             rightEigenMPrim(2,4) =  0.0
             rightEigenMPrim(3,4) =  1.0
             rightEigenMPrim(4,4) =  1.0

          end if

          eigenvect = MATMUL(rightEigenMPrim,jacConsPrim)

        end function compute_n1_righteigenvector_ns2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenvector for the hyperbolic terms
        !> in the (x+y)-direction. By denoting L the left eigenmatrix, the
        !> result of the function is L[k,:]
        !
        !> @date
        !> 02_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_n2_lefteigenvector_ns2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind), dimension(ne,ne) :: jacPrimCons
          real(rkind), dimension(ne,ne) :: leftEigenMPrim
          real(rkind) :: c
          real(rkind) :: c2_inv
          real(rkind) :: ratio

          
          !computation of J, the jacobian matrix from primitive
          !to conservative variables
          jacPrimCons = compute_jacobian_prim_to_cons(nodes)

          !left eigenmatrix for the primitive
          !variables, L_p
          c = speed_of_sound(nodes)

          if(rkind.eq.8) then

             c2_inv = 1.0d0/c**2
             ratio = 0.25d0*nodes(1)*c*Sqrt(2.0d0)

             leftEigenMPrim(1,1) =  0.0d0
             leftEigenMPrim(2,1) = -0.5d0
             leftEigenMPrim(3,1) =  0.5d0
             leftEigenMPrim(4,1) =  0.0d0

             leftEigenMPrim(1,2) =  1.0d0
             leftEigenMPrim(2,2) =  0.0d0
             leftEigenMPrim(3,2) =  0.0d0
             leftEigenMPrim(4,2) = -c2_inv

             leftEigenMPrim(1,3) =  0.0d0
             leftEigenMPrim(2,3) = -ratio
             leftEigenMPrim(3,3) = -ratio
             leftEigenMPrim(4,3) =  0.5d0

             leftEigenMPrim(1,4) =  0.0d0
             leftEigenMPrim(2,4) =  ratio
             leftEigenMPrim(3,4) =  ratio
             leftEigenMPrim(4,4) =  0.5d0

          else

             c2_inv = 1.0/c**2
             ratio = 0.25*nodes(1)*c*Sqrt(2.0)

             leftEigenMPrim(1,1) =  0.0
             leftEigenMPrim(2,1) = -0.5
             leftEigenMPrim(3,1) =  0.5
             leftEigenMPrim(4,1) =  0.0

             leftEigenMPrim(1,2) =  1.0
             leftEigenMPrim(2,2) =  0.0
             leftEigenMPrim(3,2) =  0.0
             leftEigenMPrim(4,2) = -c2_inv

             leftEigenMPrim(1,3) =  0.0
             leftEigenMPrim(2,3) = -ratio
             leftEigenMPrim(3,3) = -ratio
             leftEigenMPrim(4,3) =  0.5

             leftEigenMPrim(1,4) =  0.0
             leftEigenMPrim(2,4) =  ratio
             leftEigenMPrim(3,4) =  ratio
             leftEigenMPrim(4,4) =  0.5

          end if

          !compute the left eigenmatrix by L = L_p.J
          eigenvect = MATMUL(jacPrimCons,leftEigenMPrim)

        end function compute_n2_lefteigenvector_ns2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the right eigenvector for the hyperbolic terms
        !> in the (x+y)-direction. By denoting R the right eigenmatrix, the
        !> result of the function is R[k,:]
        !
        !> @date
        !> 02_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvector at the location of the grid point
        !--------------------------------------------------------------
        function compute_n2_righteigenvector_ns2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect

          real(rkind), dimension(ne,ne) :: jacConsPrim
          real(rkind), dimension(ne,ne) :: rightEigenMPrim
          real(rkind) :: c
          real(rkind) :: c2_inv
          real(rkind) :: ratio

          
          !computation of J, the jacobian matrix from primitive
          !to conservative variables
          jacConsPrim = compute_jacobian_cons_to_prim(nodes)

          !right eigenmatrix for the primitive
          !variables, L_p
          c = speed_of_sound(nodes)

          if(rkind.eq.8) then

             c2_inv = 1.0d0/c**2
             ratio = 0.5d0*Sqrt(2.0d0)/(nodes(1)*c)

             rightEigenMPrim(1,1) =  0.0d0
             rightEigenMPrim(2,1) =  1.0d0
             rightEigenMPrim(3,1) =  c2_inv
             rightEigenMPrim(4,1) =  c2_inv

             rightEigenMPrim(1,2) = -1.0d0
             rightEigenMPrim(2,2) =  0.0d0
             rightEigenMPrim(3,2) = -ratio
             rightEigenMPrim(4,2) =  ratio

             rightEigenMPrim(1,3) =  1.0d0
             rightEigenMPrim(2,3) =  0.0d0
             rightEigenMPrim(3,3) = -ratio
             rightEigenMPrim(4,3) =  ratio

             rightEigenMPrim(1,4) =  0.0d0
             rightEigenMPrim(2,4) =  0.0d0
             rightEigenMPrim(3,4) =  1.0d0
             rightEigenMPrim(4,4) =  1.0d0

          else

             c2_inv = 1.0/c**2
             ratio = 0.5*Sqrt(2.0)/(nodes(1)*c)

             rightEigenMPrim(1,1) =  0.0
             rightEigenMPrim(2,1) =  1.0
             rightEigenMPrim(3,1) =  c2_inv
             rightEigenMPrim(4,1) =  c2_inv

             rightEigenMPrim(1,2) = -1.0
             rightEigenMPrim(2,2) =  0.0
             rightEigenMPrim(3,2) = -ratio
             rightEigenMPrim(4,2) =  ratio

             rightEigenMPrim(1,3) =  1.0
             rightEigenMPrim(2,3) =  0.0
             rightEigenMPrim(3,3) = -ratio
             rightEigenMPrim(4,3) =  ratio

             rightEigenMPrim(1,4) =  0.0
             rightEigenMPrim(2,4) =  0.0
             rightEigenMPrim(3,4) =  1.0
             rightEigenMPrim(4,4) =  1.0

          end if

          !compute the left eigenmatrix by L = L_p.J
          eigenvect = MATMUL(rightEigenMPrim,jacConsPrim)

        end function compute_n2_righteigenvector_ns2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the transverse matrix for the (x-y)-direction
        !
        !> @date
        !> 02_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvector at the location of the grid point
        !--------------------------------------------------------------
        function compute_n1_transM_ns2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect

          
          real(rkind), dimension(ne,ne) :: jacConsPrim
          real(rkind), dimension(ne,ne) :: jacPrimCons
          real(rkind), dimension(ne,ne) :: n1transMPrim
          real(rkind)                   :: ratio
          real(rkind)                   :: u_x
          real(rkind)                   :: u_y
          real(rkind)                   :: c
          real(rkind)                   :: u_av
          real(rkind)                   :: E_c
          real(rkind)                   :: rho_av
          real(rkind)                   :: V_av

          jacConsPrim = compute_jacobian_cons_to_prim(nodes)
          jacPrimCons = compute_jacobian_prim_to_cons(nodes)

          
          
          if(rkind.eq.8) then

             ratio = 0.5d0*Sqrt(2.0d0)
             u_x   = nodes(2)/nodes(1)
             u_y   = nodes(3)/nodes(1) 
             c     = speed_of_sound(nodes)
             
             u_av   = ratio*(u_x+u_y)
             E_c    = ratio*nodes(1)*c**2
             rho_av = ratio*nodes(1)
             V_av   = ratio/nodes(1)

             n1transMPrim(1,1) =  u_av
             n1transMPrim(2,1) =  rho_av
             n1transMPrim(3,1) =  rho_av
             n1transMPrim(4,1) =  0.0d0

             n1transMPrim(1,2) =  0.0d0
             n1transMPrim(2,2) =  u_av
             n1transMPrim(3,2) =  0.0d0
             n1transMPrim(4,2) =  V_av

             n1transMPrim(1,3) =  0.0d0
             n1transMPrim(2,3) =  0.0d0
             n1transMPrim(3,3) =  u_av
             n1transMPrim(4,3) =  V_av

             n1transMPrim(1,4) =  0.0d0
             n1transMPrim(2,4) =  E_c
             n1transMPrim(3,4) =  E_c
             n1transMPrim(4,4) =  u_av

          else

             ratio = 0.5*Sqrt(2.0)
             u_x   = nodes(2)/nodes(1)
             u_y   = nodes(3)/nodes(1) 
             c     = speed_of_sound(nodes)
             
             u_av   = ratio*(u_x+u_y)
             E_c    = ratio*nodes(1)*c**2
             rho_av = ratio*nodes(1)
             V_av   = ratio/nodes(1)

             n1transMPrim(1,1) =  u_av
             n1transMPrim(2,1) =  rho_av
             n1transMPrim(3,1) =  rho_av
             n1transMPrim(4,1) =  0.0

             n1transMPrim(1,2) =  0.0
             n1transMPrim(2,2) =  u_av
             n1transMPrim(3,2) =  0.0
             n1transMPrim(4,2) =  V_av

             n1transMPrim(1,3) =  0.0
             n1transMPrim(2,3) =  0.0
             n1transMPrim(3,3) =  u_av
             n1transMPrim(4,3) =  V_av

             n1transMPrim(1,4) =  0.0
             n1transMPrim(2,4) =  E_c
             n1transMPrim(3,4) =  E_c
             n1transMPrim(4,4) =  u_av

          end if

          eigenvect = MATMUL(MATMUL(jacPrimCons,n1transMPrim),jacConsPrim)

        end function compute_n1_transM_ns2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the transverse matrix for the (x-y)-direction
        !
        !> @date
        !> 02_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvector at the location of the grid point
        !--------------------------------------------------------------
        function compute_n2_transM_ns2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect

          
          real(rkind), dimension(ne,ne) :: jacConsPrim
          real(rkind), dimension(ne,ne) :: jacPrimCons
          real(rkind), dimension(ne,ne) :: n2transMPrim
          real(rkind)                   :: ratio
          real(rkind)                   :: u_x
          real(rkind)                   :: u_y
          real(rkind)                   :: c
          real(rkind)                   :: u_av
          real(rkind)                   :: E_c
          real(rkind)                   :: rho_av
          real(rkind)                   :: V_av


          jacConsPrim = compute_jacobian_cons_to_prim(nodes)
          jacPrimCons = compute_jacobian_prim_to_cons(nodes)
          
          
          if(rkind.eq.8) then

             ratio = 0.5d0*Sqrt(2.0d0)
             u_x   = nodes(2)/nodes(1)
             u_y   = nodes(3)/nodes(1) 
             c     = speed_of_sound(nodes)
             
             u_av   = ratio*(u_x-u_y)
             E_c    = ratio*nodes(1)*c**2
             rho_av = ratio*nodes(1)
             V_av   = ratio/nodes(1)

             n2transMPrim(1,1) =  u_av
             n2transMPrim(2,1) =  rho_av
             n2transMPrim(3,1) = -rho_av
             n2transMPrim(4,1) =  0.0d0

             n2transMPrim(1,2) =  0.0d0
             n2transMPrim(2,2) =  u_av
             n2transMPrim(3,2) =  0.0d0
             n2transMPrim(4,2) =  V_av

             n2transMPrim(1,3) =  0.0d0
             n2transMPrim(2,3) =  0.0d0
             n2transMPrim(3,3) =  u_av
             n2transMPrim(4,3) = -V_av

             n2transMPrim(1,4) =  0.0d0
             n2transMPrim(2,4) =  E_c
             n2transMPrim(3,4) = -E_c
             n2transMPrim(4,4) =  u_av

          else

             ratio = 0.5*Sqrt(2.0)
             u_x   = nodes(2)/nodes(1)
             u_y   = nodes(3)/nodes(1) 
             c     = speed_of_sound(nodes)
             
             u_av   = ratio*(u_x-u_y)
             E_c    = ratio*nodes(1)*c**2
             rho_av = ratio*nodes(1)
             V_av   = ratio/nodes(1)

             n2transMPrim(1,1) =  u_av
             n2transMPrim(2,1) =  rho_av
             n2transMPrim(3,1) = -rho_av
             n2transMPrim(4,1) =  0.0

             n2transMPrim(1,2) =  0.0
             n2transMPrim(2,2) =  u_av
             n2transMPrim(3,2) =  0.0
             n2transMPrim(4,2) =  V_av

             n2transMPrim(1,3) =  0.0
             n2transMPrim(2,3) =  0.0
             n2transMPrim(3,3) =  u_av
             n2transMPrim(4,3) = -V_av

             n2transMPrim(1,4) =  0.0
             n2transMPrim(2,4) =  E_c
             n2transMPrim(3,4) = -E_c
             n2transMPrim(4,4) =  u_av

          end if

          eigenvect = MATMUL(MATMUL(jacPrimCons,n2transMPrim),jacConsPrim)

        end function compute_n2_transM_ns2d

      end module ns2d_ncoords_module
      

      module ns_vdw2d_prim_module

        use dim2d_prim_module, only :
     $     classical_pressure_local

        use dim2d_parameters, only :
     $     cv_r

        use parameters_input, only :
     $     ne

        use parameters_kind, only :
     $     rkind

        implicit none

        private
        public ::
     $       compute_prim_var_ns_vdw2d,
     $       compute_cons_var_ns_vdw2d,
     $       compute_jacobian_prim_to_cons_ns_vdw2d,
     $       compute_jacobian_cons_to_prim_ns_vdw2d,
     $       compute_x_transM_ns_vdw2d,
     $       compute_y_transM_ns_vdw2d,
     $       compute_x_eigenvalues_ns_vdw2d,
     $       compute_y_eigenvalues_ns_vdw2d,
     $       compute_x_lefteigenvector_ns_vdw2d,
     $       compute_x_righteigenvector_ns_vdw2d,
     $       compute_y_lefteigenvector_ns_vdw2d,
     $       compute_y_righteigenvector_ns_vdw2d

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the primitive variables from the conservative
        !> variables
        !
        !> @date
        !> 04_02_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes_cons
        !> conservative variables
        !
        !>@return nodes_prim
        !> primitive variables
        !--------------------------------------------------------------
        function compute_prim_var_ns_vdw2d(nodes_cons) result(nodes_prim)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes_cons
          real(rkind), dimension(ne)             :: nodes_prim


          nodes_prim(1) = nodes_cons(1)                        !mass density
          nodes_prim(2) = nodes_cons(2)/nodes_cons(1)          !velocity-x
          nodes_prim(3) = nodes_cons(3)/nodes_cons(1)          !velocity-y
          nodes_prim(4) = classical_pressure_local(nodes_cons) !pressure

        end function compute_prim_var_ns_vdw2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the conservative variables from the
        !> primitive variables
        !
        !> @date
        !> 04_02_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes_prim
        !> primitive variables
        !
        !>@return nodes_cons
        !> conservative variables
        !--------------------------------------------------------------
        function compute_cons_var_ns_vdw2d(nodes_prim) result(nodes_cons)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes_prim
          real(rkind), dimension(ne)             :: nodes_cons

          real(rkind) :: T

          nodes_cons(1) = nodes_prim(1)
          nodes_cons(2) = nodes_prim(1)*nodes_prim(2)
          nodes_cons(3) = nodes_prim(1)*nodes_prim(3)

          T = (nodes_prim(4) + 3.0d0*nodes_prim(1)**2)*
     $        (3.0d0 - nodes_prim(1))/
     $        (8.0d0*nodes_prim(1))

          nodes_cons(4) = 0.5d0*nodes_prim(1)*(
     $                       nodes_prim(2)**2+
     $                       nodes_prim(3)**2) + 
     $                    nodes_prim(1)*(
     $                       8.0d0/3.0d0*cv_r*T -
     $                       3.0d0*nodes_prim(1))


        end function compute_cons_var_ns_vdw2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the Jacobian matrix for primitive to
        !> to conservative variables
        !
        !> @date
        !> 03_02_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density
        !
        !>@param ux
        !> velocity along the x-direction
        !
        !>@param uy
        !> velocity along the y-direction
        !
        !>@param P
        !> pressure
        !
        !>@return jac_matrix
        !> jacobian matrix for primitive to conservative
        !> variables \f$ J^p_v = \frac{\partial p}{\partial v} \f$
        !--------------------------------------------------------------
        function compute_jacobian_prim_to_cons_ns_vdw2d(md,ux,uy,P)
     $     result(jac_matrix)

          implicit none

          real(rkind)                  , intent(in) :: md
          real(rkind)                  , intent(in) :: ux
          real(rkind)                  , intent(in) :: uy
          real(rkind)                  , intent(in) :: P  
          real(rkind), dimension(ne,ne)             :: jac_matrix

          real(rkind) :: ratio
          real(rkind) :: dPdrho


          if(rkind.eq.8) then             

             ratio = 3.0d0/(cv_r*(-3.0d0+md))

             dPdrho =
     $            -(
     $               3.0d0*(ux**2+uy**2+12.0d0*md) +
     $               2.0d0*cv_r*(P+9.0d0*(md-2.0d0)*md)
     $            )/
     $            (
     $               2.d0*cv_r*(-3.0d0+md)
     $            )

             jac_matrix(1,1) = 1.0d0
             jac_matrix(2,1) = 0.0d0
             jac_matrix(3,1) = 0.0d0
             jac_matrix(4,1) = 0.0d0

             jac_matrix(1,2) = -  ux/md
             jac_matrix(2,2) = 1.0d0/md
             jac_matrix(3,2) = 0.0d0
             jac_matrix(4,2) = 0.0d0

             jac_matrix(1,3) = -  uy/md
             jac_matrix(2,3) = 0.0d0
             jac_matrix(3,3) = 1.0d0/md
             jac_matrix(4,3) = 0.0d0

             jac_matrix(1,4) =  dPdrho
             jac_matrix(2,4) =  ratio*ux
             jac_matrix(3,4) =  ratio*uy
             jac_matrix(4,4) = -ratio

          else

             ratio = 3.0/(cv_r*(-3.0+md))

             dPdrho =
     $            -(
     $               3.0*(ux**2+uy**2+12.0*md) +
     $               2.0*cv_r*(P+9.0*(md-2.0)*md)
     $            )/
     $            (
     $               2.*cv_r*(-3.0+md)
     $            )

             jac_matrix(1,1) = 1.0
             jac_matrix(2,1) = 0.0
             jac_matrix(3,1) = 0.0
             jac_matrix(4,1) = 0.0

             jac_matrix(1,2) = -ux/md
             jac_matrix(2,2) = 1.0/md
             jac_matrix(3,2) = 0.0
             jac_matrix(4,2) = 0.0

             jac_matrix(1,3) = -uy/md
             jac_matrix(2,3) = 0.0
             jac_matrix(3,3) = 1.0/md
             jac_matrix(4,3) = 0.0

             jac_matrix(1,4) =  dPdrho
             jac_matrix(2,4) =  ratio*ux
             jac_matrix(3,4) =  ratio*uy
             jac_matrix(4,4) = -ratio

          end if

        end function compute_jacobian_prim_to_cons_ns_vdw2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the Jacobian matrix for conservative to
        !> to primitive variables
        !
        !> @date
        !> 03_02_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density
        !
        !>@param ux
        !> velocity along the x-direction
        !
        !>@param uy
        !> velocity along the y-direction
        !
        !>@param P
        !> pressure
        !
        !>@return jac_matrix
        !> jacobian matrix for conservative to primitive
        !> variables \f$ J^v_p = \frac{\partial v}{\partial p} \f$
        !--------------------------------------------------------------
        function compute_jacobian_cons_to_prim_ns_vdw2d(md,ux,uy,P)
     $     result(jac_matrix)

          implicit none

          real(rkind)                  , intent(in) :: md
          real(rkind)                  , intent(in) :: ux
          real(rkind)                  , intent(in) :: uy
          real(rkind)                  , intent(in) :: P
          real(rkind), dimension(ne,ne)             :: jac_matrix


          real(rkind) :: dEdrho


          if(rkind.eq.8) then

             dEdrho =
     $            0.5d0*(ux**2+uy**2) -
     $            6.0d0*md -
     $            1.0d0/3.0d0*cv_r*(P + 9.0d0*(md-2.0d0)*md)
             
             jac_matrix(1,1) = 1.0d0
             jac_matrix(2,1) = 0.0d0
             jac_matrix(3,1) = 0.0d0
             jac_matrix(4,1) = 0.0d0
                        
             jac_matrix(1,2) = ux
             jac_matrix(2,2) = md
             jac_matrix(3,2) = 0.0d0
             jac_matrix(4,2) = 0.0d0
                        
             jac_matrix(1,3) = uy
             jac_matrix(2,3) = 0.0d0
             jac_matrix(3,3) = md
             jac_matrix(4,3) = 0.0d0
                        
             jac_matrix(1,4) = dEdrho
             jac_matrix(2,4) = md*ux
             jac_matrix(3,4) = md*uy
             jac_matrix(4,4) = cv_r*(1.0d0-md/3.0d0)

          else

             dEdrho =
     $            0.5*(ux**2+uy**2) -
     $            6.0*md -
     $            1.0/3.0*cv_r*(P + 9.0*(md-2.0)*md)
             
             jac_matrix(1,1) = 1.0
             jac_matrix(2,1) = 0.0
             jac_matrix(3,1) = 0.0
             jac_matrix(4,1) = 0.0
                        
             jac_matrix(1,2) = ux
             jac_matrix(2,2) = md
             jac_matrix(3,2) = 0.0
             jac_matrix(4,2) = 0.0
                        
             jac_matrix(1,3) = uy
             jac_matrix(2,3) = 0.0
             jac_matrix(3,3) = md
             jac_matrix(4,3) = 0.0
                        
             jac_matrix(1,4) = dEdrho
             jac_matrix(2,4) = md*ux
             jac_matrix(3,4) = md*uy
             jac_matrix(4,4) = cv_r*(1.0-md/3.0)

          end if

        end function compute_jacobian_cons_to_prim_ns_vdw2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of $A^2$, the hyperbolic matrix in the
        !> y-direction
        !
        !> @date
        !> 04_02_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density
        !
        !>@param uy
        !> velocity in the y-direction
        !
        !>@param c
        !> speed of sound
        !
        !>@return transM
        !> hyperbolic matrix in the y-direction
        !--------------------------------------------------------------
        function compute_x_transM_ns_vdw2d(md,uy,c) result(transM)

          implicit none

          real(rkind)                  , intent(in) :: md
          real(rkind)                  , intent(in) :: uy
          real(rkind)                  , intent(in) :: c
          real(rkind), dimension(ne,ne)             :: transM


          if(rkind.eq.8) then

             transM(1,1) = uy
             transM(2,1) = 0.0d0
             transM(3,1) = md
             transM(4,1) = 0.0d0

             transM(1,2) = 0.0d0
             transM(2,2) = uy
             transM(3,2) = 0.0d0
             transM(4,2) = 0.0d0

             transM(1,3) = 0.0d0
             transM(2,3) = 0.0d0
             transM(3,3) = uy
             transM(4,3) = 1.0d0/md

             transM(1,4) = 0.0d0
             transM(2,4) = 0.0d0
             transM(3,4) = c**2*md
             transM(4,4) = uy

          else
             
             transM(1,1) = uy
             transM(2,1) = 0.0
             transM(3,1) = md
             transM(4,1) = 0.0

             transM(1,2) = 0.0
             transM(2,2) = uy
             transM(3,2) = 0.0
             transM(4,2) = 0.0

             transM(1,3) = 0.0
             transM(2,3) = 0.0
             transM(3,3) = uy
             transM(4,3) = 1.0/md

             transM(1,4) = 0.0
             transM(2,4) = 0.0
             transM(3,4) = c**2*md
             transM(4,4) = uy

          end if

        end function compute_x_transM_ns_vdw2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the $A^1$, the hyperbolic matrix in the
        !> x-direction
        !
        !> @date
        !> 04_02_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density
        !
        !>@param ux
        !> velocity in the x-direction
        !
        !>@param c
        !> speed of sound
        !
        !>@return transM
        !> hyperbolic matrix in the x-direction
        !--------------------------------------------------------------
        function compute_y_transM_ns_vdw2d(md,ux,c) result(transM)

          implicit none

          real(rkind)                  , intent(in) :: md
          real(rkind)                  , intent(in) :: ux
          real(rkind)                  , intent(in) :: c
          real(rkind), dimension(ne,ne)             :: transM


          if(rkind.eq.8) then
             
             transM(1,1) = ux
             transM(2,1) = md
             transM(3,1) = 0.0d0
             transM(4,1) = 0.0d0

             transM(1,2) = 0.0d0
             transM(2,2) = ux
             transM(3,2) = 0.0d0
             transM(4,2) = 1.0d0/md

             transM(1,3) = 0.0d0
             transM(2,3) = 0.0d0
             transM(3,3) = ux
             transM(4,3) = 0.0d0

             transM(1,4) = 0.0d0
             transM(2,4) = c**2*md
             transM(3,4) = 0.0d0
             transM(4,4) = ux

          else
             
             transM(1,1) = ux
             transM(2,1) = md
             transM(3,1) = 0.0
             transM(4,1) = 0.0

             transM(1,2) = 0.0
             transM(2,2) = ux
             transM(3,2) = 0.0
             transM(4,2) = 1.0/md

             transM(1,3) = 0.0
             transM(2,3) = 0.0
             transM(3,3) = ux
             transM(4,3) = 0.0

             transM(1,4) = 0.0
             transM(2,4) = c**2*md
             transM(3,4) = 0.0
             transM(4,4) = ux

          end if

        end function compute_y_transM_ns_vdw2d
        

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the x-direction in primitive form
        !
        !> @date
        !> 03_02_2015 - initial version - J.L. Desmarais
        !
        !>@param uy
        !> velocity along the y-direction
        !
        !>@param c
        !> speed of sound
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_x_eigenvalues_ns_vdw2d(ux,c) result(eigenvalues)

          implicit none

          real(rkind)               , intent(in) :: ux
          real(rkind)               , intent(in) :: c
          real(rkind), dimension(ne)             :: eigenvalues

          eigenvalues(1) = ux
          eigenvalues(2) = ux
          eigenvalues(3) = ux-c
          eigenvalues(4) = ux+c

        end function compute_x_eigenvalues_ns_vdw2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the y-direction in primitive form
        !
        !> @date
        !> 03_02_2015 - initial version - J.L. Desmarais
        !
        !>@param uy
        !> velocity along the y-direction
        !
        !>@param c
        !> speed of sound
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_y_eigenvalues_ns_vdw2d(uy,c) result(eigenvalues)

          implicit none

          real(rkind)               , intent(in) :: uy
          real(rkind)               , intent(in) :: c
          real(rkind), dimension(ne)             :: eigenvalues

          eigenvalues(1) = uy
          eigenvalues(2) = uy
          eigenvalues(3) = uy-c
          eigenvalues(4) = uy+c

        end function compute_y_eigenvalues_ns_vdw2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenmatrix for the hyperbolic
        !> terms in the x-direction in primitive form
        !
        !> @date
        !> 03_02_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density
        !
        !>@param c
        !> speed of sound
        !
        !>@return eigenvect
        !> left eigenmatrix
        !--------------------------------------------------------------
        function compute_x_lefteigenvector_ns_vdw2d(md,c) result(eigenvect)

          implicit none

          real(rkind)                  , intent(in) :: md
          real(rkind)                  , intent(in) :: c
          real(rkind), dimension(ne,ne)             :: eigenvect
          
          real(rkind) :: q_c

          
          if(rkind.eq.8) then

             q_c = 0.5d0*md*c

             eigenvect(1,1) = 0.0d0
             eigenvect(2,1) = 0.0d0
             eigenvect(3,1) = 1.0d0
             eigenvect(4,1) = 0.0d0

             eigenvect(1,2) = 1.0d0
             eigenvect(2,2) = 0.0d0
             eigenvect(3,2) = 0.0d0
             eigenvect(4,2) =-1.0d0/c**2

             eigenvect(1,3) = 0.0d0
             eigenvect(2,3) =-q_c
             eigenvect(3,3) = 0.0d0
             eigenvect(4,3) = 0.5d0

             eigenvect(1,4) = 0.0d0
             eigenvect(2,4) = q_c
             eigenvect(3,4) = 0.0d0
             eigenvect(4,4) = 0.5d0

          else

             q_c = 0.5*md*c

             eigenvect(1,1) = 0.0
             eigenvect(2,1) = 0.0
             eigenvect(3,1) = 1.0
             eigenvect(4,1) = 0.0

             eigenvect(1,2) = 1.0
             eigenvect(2,2) = 0.0
             eigenvect(3,2) = 0.0
             eigenvect(4,2) =-1.0/c**2

             eigenvect(1,3) = 0.0
             eigenvect(2,3) =-q_c
             eigenvect(3,3) = 0.0
             eigenvect(4,3) = 0.5

             eigenvect(1,4) = 0.0
             eigenvect(2,4) = q_c
             eigenvect(3,4) = 0.0
             eigenvect(4,4) = 0.5

          end if

        end function compute_x_lefteigenvector_ns_vdw2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the right eigenmatrix for the hyperbolic
        !> terms in the x-direction in primitive form
        !
        !> @date
        !> 03_02_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density
        !
        !>@param c
        !> speed of sound
        !
        !>@return eigenvect
        !> left eigenmatrix
        !--------------------------------------------------------------
        function compute_x_righteigenvector_ns_vdw2d(md,c) result(eigenvect)

          implicit none

          real(rkind)                  , intent(in) :: md
          real(rkind)                  , intent(in) :: c
          real(rkind), dimension(ne,ne)             :: eigenvect


          real(rkind) :: inv_c
          real(rkind) :: inv_q_c
          

          if(rkind.eq.8) then

             inv_c   = 1.0d0/c**2
             inv_q_c = 1.0d0/(md*c)

             eigenvect(1,1) = 0.0d0
             eigenvect(2,1) = 1.0d0
             eigenvect(3,1) = inv_c
             eigenvect(4,1) = inv_c
             
             eigenvect(1,2) = 0.0d0
             eigenvect(2,2) = 0.0d0
             eigenvect(3,2) =-inv_q_c
             eigenvect(4,2) = inv_q_c

             eigenvect(1,3) = 1.0d0
             eigenvect(2,3) = 0.0d0
             eigenvect(3,3) = 0.0d0
             eigenvect(4,3) = 0.0d0

             eigenvect(1,4) = 0.0d0
             eigenvect(2,4) = 0.0d0
             eigenvect(3,4) = 1.0d0
             eigenvect(4,4) = 1.0d0

          else

             inv_c   = 1.0/c**2
             inv_q_c = 1.0/(md*c)

             eigenvect(1,1) = 0.0
             eigenvect(2,1) = 1.0
             eigenvect(3,1) = inv_c
             eigenvect(4,1) = inv_c
             
             eigenvect(1,2) = 0.0
             eigenvect(2,2) = 0.0
             eigenvect(3,2) =-inv_q_c
             eigenvect(4,2) = inv_q_c

             eigenvect(1,3) = 1.0
             eigenvect(2,3) = 0.0
             eigenvect(3,3) = 0.0
             eigenvect(4,3) = 0.0

             eigenvect(1,4) = 0.0
             eigenvect(2,4) = 0.0
             eigenvect(3,4) = 1.0
             eigenvect(4,4) = 1.0

          end if

        end function compute_x_righteigenvector_ns_vdw2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenmatrix for the hyperbolic
        !> terms in the y-direction in primitive form
        !
        !> @date
        !> 03_02_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density
        !
        !>@param c
        !> speed of sound
        !
        !>@return eigenvect
        !> left eigenmatrix
        !--------------------------------------------------------------
        function compute_y_lefteigenvector_ns_vdw2d(md,c) result(eigenvect)

          implicit none

          real(rkind)                  , intent(in) :: md
          real(rkind)                  , intent(in) :: c
          real(rkind), dimension(ne,ne)             :: eigenvect

          real(rkind) :: q_c


          if(rkind.eq.8) then

             q_c = 0.5d0*md*c

             eigenvect(1,1) = 0.0d0
             eigenvect(2,1) = 1.0d0
             eigenvect(3,1) = 0.0d0
             eigenvect(4,1) = 0.0d0

             eigenvect(1,2) = 1.0d0
             eigenvect(2,2) = 0.0d0
             eigenvect(3,2) = 0.0d0
             eigenvect(4,2) =-1.0d0/c**2

             eigenvect(1,3) = 0.0d0
             eigenvect(2,3) = 0.0d0
             eigenvect(3,3) =-q_c
             eigenvect(4,3) = 0.5d0

             eigenvect(1,4) = 0.0d0
             eigenvect(2,4) = 0.0d0
             eigenvect(3,4) = q_c
             eigenvect(4,4) = 0.5d0

          else

             q_c = 0.5*md*c

             eigenvect(1,1) = 0.0
             eigenvect(2,1) = 1.0
             eigenvect(3,1) = 0.0
             eigenvect(4,1) = 0.0

             eigenvect(1,2) = 1.0
             eigenvect(2,2) = 0.0
             eigenvect(3,2) = 0.0
             eigenvect(4,2) =-1.0/c**2

             eigenvect(1,3) = 0.0
             eigenvect(2,3) = 0.0
             eigenvect(3,3) =-q_c
             eigenvect(4,3) = 0.5

             eigenvect(1,4) = 0.0
             eigenvect(2,4) = 0.0
             eigenvect(3,4) = q_c
             eigenvect(4,4) = 0.5

          end if

        end function compute_y_lefteigenvector_ns_vdw2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the right eigenmatrix for the hyperbolic
        !> terms in the y-direction in primitive form
        !
        !> @date
        !> 03_02_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density
        !
        !>@param c
        !> speed of sound
        !
        !>@return eigenvect
        !> left eigenmatrix
        !--------------------------------------------------------------
        function compute_y_righteigenvector_ns_vdw2d(md,c) result(eigenvect)

          implicit none

          real(rkind)                  , intent(in) :: md
          real(rkind)                  , intent(in) :: c
          real(rkind), dimension(ne,ne)             :: eigenvect

          real(rkind) :: inv_c
          real(rkind) :: inv_q_c
          

          if(rkind.eq.8) then

             inv_c   = 1.0d0/c**2
             inv_q_c = 1.0d0/(md*c)

             eigenvect(1,1) = 0.0d0
             eigenvect(2,1) = 1.0d0
             eigenvect(3,1) = inv_c
             eigenvect(4,1) = inv_c

             eigenvect(1,2) = 1.0d0
             eigenvect(2,2) = 0.0d0
             eigenvect(3,2) = 0.0d0
             eigenvect(4,2) = 0.0d0

             eigenvect(1,3) = 0.0d0
             eigenvect(2,3) = 0.0d0
             eigenvect(3,3) =-inv_q_c
             eigenvect(4,3) = inv_q_c

             eigenvect(1,4) = 0.0d0
             eigenvect(2,4) = 0.0d0
             eigenvect(3,4) = 1.0d0
             eigenvect(4,4) = 1.0d0

          else

             inv_c   = 1.0/c**2
             inv_q_c = 1.0/(md*c)

             eigenvect(1,1) = 0.0
             eigenvect(2,1) = 1.0
             eigenvect(3,1) = inv_c
             eigenvect(4,1) = inv_c

             eigenvect(1,2) = 1.0
             eigenvect(2,2) = 0.0
             eigenvect(3,2) = 0.0
             eigenvect(4,2) = 0.0

             eigenvect(1,3) = 0.0
             eigenvect(2,3) = 0.0
             eigenvect(3,3) =-inv_q_c
             eigenvect(4,3) = inv_q_c

             eigenvect(1,4) = 0.0
             eigenvect(2,4) = 0.0
             eigenvect(3,4) = 1.0
             eigenvect(4,4) = 1.0

          end if

        end function compute_y_righteigenvector_ns_vdw2d

      end module ns_vdw2d_prim_module

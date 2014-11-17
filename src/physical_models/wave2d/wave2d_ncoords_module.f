      !> @file
      !> module implemeting subroutines to compute the gradients,
      !> the eigenvalues and the eigenvectors in the (x-y) and (x+y)
      !> directions for the wave2d governing equations
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implemeting subroutines to compute the gradients,
      !> the eigenvalues and the eigenvectors in the (x-y) and (x+y)
      !> directions for the wave2d governing equations
      !
      !> @date
      !> 05_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
       module wave2d_ncoords_module

        use interface_primary       , only : gradient_n_proc
        use wave2d_parameters       , only : c
        use wave2d_prim_module      , only : position,
     $                                       velocity_x,
     $                                       velocity_y
        use parameters_input        , only : nx,ny,ne
        use parameters_kind         , only : ikind, rkind


        private
        public :: 
     $       compute_n_eigenvalues_wave2d,
     $       compute_n1_lefteigenvector_wave2d,
     $       compute_n1_righteigenvector_wave2d,
     $       compute_n2_lefteigenvector_wave2d,
     $       compute_n2_righteigenvector_wave2d,
     $       compute_n1_transM_wave2d,
     $       compute_n2_transM_wave2d


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the eigenvalues for the hyperbolic terms
        !> in the (x-y)-direction
        !
        !> @date
        !> 05_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_n_eigenvalues_wave2d(nodes) result(eigenvalues)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne)             :: eigenvalues


          real(rkind) :: node_s

          node_s = nodes(1)

          if(rkind.eq.8) then
             eigenvalues(1) =  0.0d0
             eigenvalues(2) = -c**2
             eigenvalues(3) =  c**2
          else
             eigenvalues(1) =  0.0
             eigenvalues(2) = -c**2
             eigenvalues(3) =  c**2
          end if

        end function compute_n_eigenvalues_wave2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenvector for the hyperbolic terms
        !> in the (x-y)-direction. By denoting L the left eigenmatrix, the
        !> result of the function is L[k,:]
        !
        !> @date
        !> 05_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_n1_lefteigenvector_wave2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind) :: node_s

          node_s = nodes(1)


          if(rkind.eq.8) then
             eigenvect(1,1) =  0.0d0
             eigenvect(2,1) =  0.5d0
             eigenvect(3,1) =  0.5d0
                              
             eigenvect(1,2) = -0.25d0*Sqrt(2.0d0)
             eigenvect(2,2) = -0.25d0
             eigenvect(3,2) =  0.25d0
                              
             eigenvect(1,3) =  0.25d0*Sqrt(2.0d0)
             eigenvect(2,3) = -0.25d0
             eigenvect(3,3) =  0.25d0
                              
          else                
             eigenvect(1,1) =  0.0
             eigenvect(2,1) =  0.5
             eigenvect(3,1) =  0.5
                              
             eigenvect(1,2) = -0.25*Sqrt(2.0)
             eigenvect(2,2) = -0.25
             eigenvect(3,2) =  0.25
                              
             eigenvect(1,3) =  0.25*Sqrt(2.0)
             eigenvect(2,3) = -0.25
             eigenvect(3,3) =  0.25

          end if

        end function compute_n1_lefteigenvector_wave2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the right eigenvector for the hyperbolic terms
        !> in the (x-y)-direction. By denoting R the right eigenmatrix, the
        !> result of the function is R[k,:]
        !
        !> @date
        !> 05_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvector at the location of the grid point
        !--------------------------------------------------------------
        function compute_n1_righteigenvector_wave2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind) :: node_s

          node_s = nodes(1)

          if(rkind.eq.8) then
             eigenvect(1,1) =  0.0d0
             eigenvect(2,1) = -Sqrt(2.0d0)
             eigenvect(3,1) =  Sqrt(2.0d0)
             
             eigenvect(1,2) =  1.0d0
             eigenvect(2,2) = -1.0d0
             eigenvect(3,2) = -1.0d0
                              
             eigenvect(1,3) =  1.0d0
             eigenvect(2,3) =  1.0d0
             eigenvect(3,3) =  1.0d0

          else
             eigenvect(1,1) =  0.0
             eigenvect(2,1) = -Sqrt(2.0)
             eigenvect(3,1) =  Sqrt(2.0)
             
             eigenvect(1,2) =  1.0
             eigenvect(2,2) = -1.0
             eigenvect(3,2) = -1.0
                              
             eigenvect(1,3) =  1.0
             eigenvect(2,3) =  1.0
             eigenvect(3,3) =  1.0

          end if

        end function compute_n1_righteigenvector_wave2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the left eigenvector for the hyperbolic terms
        !> in the (x+y)-direction. By denoting L the left eigenmatrix, the
        !> result of the function is L[k,:]
        !
        !> @date
        !> 06_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvalues
        !> eigenvalues at the location of the grid point
        !--------------------------------------------------------------
        function compute_n2_lefteigenvector_wave2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind) :: node_s

          node_s = nodes(1)


          if(rkind.eq.8) then
             eigenvect(1,1) =  0.0d0
             eigenvect(2,1) = -0.5d0
             eigenvect(3,1) =  0.5d0
                              
             eigenvect(1,2) =  0.25d0*Sqrt(2.0d0)
             eigenvect(2,2) =  0.25d0
             eigenvect(3,2) =  0.25d0
                              
             eigenvect(1,3) = -0.25d0*Sqrt(2.0d0)
             eigenvect(2,3) =  0.25d0
             eigenvect(3,3) =  0.25d0

          else
             eigenvect(1,1) =  0.0
             eigenvect(2,1) = -0.5
             eigenvect(3,1) =  0.5
                              
             eigenvect(1,2) =  0.25*Sqrt(2.0)
             eigenvect(2,2) =  0.25
             eigenvect(3,2) =  0.25

             eigenvect(1,3) = -0.25*Sqrt(2.0)
             eigenvect(2,3) =  0.25
             eigenvect(3,3) =  0.25
          end if

        end function compute_n2_lefteigenvector_wave2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the right eigenvector for the hyperbolic terms
        !> in the (x+y)-direction. By denoting R the right eigenmatrix, the
        !> result of the function is R[k,:]
        !
        !> @date
        !> 06_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvector at the location of the grid point
        !--------------------------------------------------------------
        function compute_n2_righteigenvector_wave2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind) :: node_s

          node_s = nodes(1)

          if(rkind.eq.8) then
             eigenvect(1,1) =  0.0d0
             eigenvect(2,1) =  Sqrt(2.0d0)
             eigenvect(3,1) = -Sqrt(2.0d0)

             eigenvect(1,2) = -1.0d0
             eigenvect(2,2) =  1.0d0
             eigenvect(3,2) =  1.0d0
                            
             eigenvect(1,3) =  1.0d0
             eigenvect(2,3) =  1.0d0
             eigenvect(3,3) =  1.0d0

          else              
             eigenvect(1,1) =  0.0
             eigenvect(2,1) =  Sqrt(2.0)
             eigenvect(3,1) = -Sqrt(2.0)
                            
             eigenvect(1,2) = -1.0
             eigenvect(2,2) =  1.0
             eigenvect(3,2) =  1.0
                            
             eigenvect(1,3) =  1.0
             eigenvect(2,3) =  1.0
             eigenvect(3,3) =  1.0

          end if       

        end function compute_n2_righteigenvector_wave2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the transverse matrix for the direction n1
        !
        !> @date
        !> 17_11_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvector at the location of the grid point
        !--------------------------------------------------------------
        function compute_n1_transM_wave2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind) :: node_s
          real(rkind) :: a

          node_s = nodes(1)

          if(rkind.eq.8) then
             a = 0.5d0*Sqrt(2.0d0)

             eigenvect(1,1) =  0.0d0
             eigenvect(2,1) = -a*c**2
             eigenvect(3,1) =  a*c**2

             eigenvect(1,2) = -a*c**2
             eigenvect(2,2) =  0.0d0
             eigenvect(3,2) =  0.0d0
                            
             eigenvect(1,3) =  a*c**2
             eigenvect(2,3) =  0.0d0
             eigenvect(3,3) =  0.0d0

          else              
             a = 0.5*Sqrt(2.0)

             eigenvect(1,1) =  0.0
             eigenvect(2,1) = -a*c**2
             eigenvect(3,1) =  a*c**2

             eigenvect(1,2) = -a*c**2
             eigenvect(2,2) =  0.0
             eigenvect(3,2) =  0.0
                            
             eigenvect(1,3) =  a*c**2
             eigenvect(2,3) =  0.0
             eigenvect(3,3) =  0.0

          end if       

        end function compute_n1_transM_wave2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the transverse matrix for the direction n2
        !
        !> @date
        !> 17_11_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@return eigenvect
        !> eigenvector at the location of the grid point
        !--------------------------------------------------------------
        function compute_n2_transM_wave2d(nodes) result(eigenvect)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne,ne)          :: eigenvect


          real(rkind) :: node_s
          real(rkind) :: a

          node_s = nodes(1)

          if(rkind.eq.8) then
             a = 0.5d0*Sqrt(2.0d0)

             eigenvect(1,1) =  0.0d0
             eigenvect(2,1) = -a*c**2
             eigenvect(3,1) = -a*c**2

             eigenvect(1,2) = -a*c**2
             eigenvect(2,2) =  0.0d0
             eigenvect(3,2) =  0.0d0
                            
             eigenvect(1,3) = -a*c**2
             eigenvect(2,3) =  0.0d0
             eigenvect(3,3) =  0.0d0

          else              
             a = 0.5*Sqrt(2.0)

             eigenvect(1,1) =  0.0
             eigenvect(2,1) = -a*c**2
             eigenvect(3,1) = -a*c**2

             eigenvect(1,2) = -a*c**2
             eigenvect(2,2) =  0.0
             eigenvect(3,2) =  0.0
                            
             eigenvect(1,3) = -a*c**2
             eigenvect(2,3) =  0.0
             eigenvect(3,3) =  0.0

          end if       

        end function compute_n2_transM_wave2d

      end module wave2d_ncoords_module
      

      !> @file
      !> module encapsulating subroutines for the computation
      !> of the boundary layers using the reflection xy boundary
      !> conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines for the computation
      !> of the boundary layers using the reflection xy boundary
      !> conditions
      !
      !> @date
      ! 23_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module reflection_xy_module

        use cg_operators_class , only : cg_operators
        use dim2d_eq_class     , only : dim2d_eq
        use parameters_constant, only : vector_x, vector_y
        use parameters_input   , only : nx,ny,ne
        use parameters_kind    , only : ikind, rkind
        
        implicit none
        
        private
        public :: reflection_x_prefactor,
     $            reflection_y_prefactor

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the prefactor for the reflection
        !> along the x-axis: the prefactor is equal to -1 or +1
        !> depending if the variable type is a vector_x or not
        !
        !> @date
        !> 23_08_2013 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> physical model
        !
        !>@param prefactor
        !> table containing the prefactor for the reflection
        !--------------------------------------------------------------
        function reflection_x_prefactor(p_model) result(prefactor)

          implicit none

          type(dim2d_eq), intent(in) :: p_model
          integer, dimension(ne)     :: prefactor
            
          integer, dimension(ne)     :: var_type
          integer :: k

          var_type = p_model%get_var_type()

          do k=1,ne
             if(var_type(k).eq.vector_x) then
                prefactor(k)=-1
             else
                prefactor(k)= 1
             end if
          end do
          
        end function reflection_x_prefactor


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the prefactor for the reflection
        !> along the y-axis: the prefactor is equal to -1 or +1
        !> depending if the variable type is a vector_y or not
        !
        !> @date
        !> 23_08_2013 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> physical model
        !
        !>@param prefactor
        !> table containing the prefactor for the reflection
        !--------------------------------------------------------------
        function reflection_y_prefactor(p_model) result(prefactor)

          implicit none

          type(dim2d_eq), intent(in) :: p_model
          integer, dimension(ne)     :: prefactor
          
          integer, dimension(ne)     :: var_type
          integer :: k

          var_type = p_model%get_var_type()

          do k=1,ne
             if(var_type(k).eq.vector_y) then
                prefactor(k)=-1
             else
                prefactor(k)= 1
             end if
          end do
          
        end function reflection_y_prefactor

      end module reflection_xy_module

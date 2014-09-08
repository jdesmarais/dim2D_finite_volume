      !> @file
      !> module implementing the functions useful to distinguish the 
      !> acoustic components of the LODI vector
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the functions useful to distinguish the 
      !> acoustic components of the LODI vector
      !
      !> @date
      ! 05_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_component_module

        use parameters_constant, only :
     $     left

        implicit none

        private
        public ::
     $       get_incoming_acoustic_component,
     $       get_other_acoustic_component,
     $       get_sign_acoustic_component

        contains


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the incoming acoustic component corresponding to the
        !> side location of the boundary conditions
        !
        !> @date
        !> 04_09_2014 - initial version - J.L. Desmarais
        !
        !>@param side
        !> side location of the boundary condition
        !
        !>@return i
        !> LODI component corresponding to the incoming acoustic amplitude (3 or 4)
        !---------------------------------------------------------------
        function get_incoming_acoustic_component(side) result(i)

          implicit none

          logical, intent(in) :: side
          integer             :: i

          
          if(side.eqv.left) then
             i = 4
          else
             i = 3
          end if

        end function get_incoming_acoustic_component



        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the other acoustic amplitude (3->4, 4->3)
        !
        !> @date
        !> 04_09_2014 - initial version - J.L. Desmarais
        !
        !>@param i
        !> LODI component corresponding to an acoustic amplitude (3 or 4)
        !
        !>@return i_star
        !> LODI component corresponding to the other acoustic amplitude (3 or 4)
        !---------------------------------------------------------------
        function get_other_acoustic_component(i) result(i_star)

          implicit none

          integer, intent(in) :: i
          integer             :: i_star

          i_star = - i + 7

        end function get_other_acoustic_component


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the sign for the acoustic amplitude (3->-1, 4->1)
        !
        !> @date
        !> 04_09_2014 - initial version - J.L. Desmarais
        !
        !>@param i
        !> LODI component corresponding to an acoustic amplitude (3 or 4)
        !
        !>@return sign
        !> sign corresponding to the accoustic LODI component
        !---------------------------------------------------------------
        function get_sign_acoustic_component(i) result(sign)

          implicit none

          integer, intent(in) :: i
          integer             :: sign

          sign = 2*i-7

        end function get_sign_acoustic_component

      end module lodi_component_module

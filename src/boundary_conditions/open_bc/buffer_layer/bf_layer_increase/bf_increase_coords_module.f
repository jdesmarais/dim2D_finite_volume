      !> @file
      !> module encapsulating the subroutines for the determination of
      !> the main layer the activated gridpoint belongs to
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the subroutines for the determination of
      !> the main layer the activated gridpoint belongs to
      !
      !> @date
      !> 19_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_increase_coords_module

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W

        use parameters_constant, only :
     $       N,S,E,W,interior
        
        use parameters_input, only :
     $       nx,ny

        use parameters_kind, only :
     $       ikind

        
        implicit none


        private
        public ::
     $       get_mainlayer_coord


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> find the main layer to which to the grid-point with
        !> coordinates expressed in the general frame as well as
        !> its neighbors [-bc_size+gen_coords(1),gen_coords(1)+bc_size]x
        !> [-bc_size+gen_coords(2),gen_coords(2)+bc_size]
        !
        !> @date
        !> 19_03_2015 - initial version - J.L. Desmarais
        !
        !>@param gen_coords
        !> coordinates of the grid points expressed in the general
        !> reference frame
        !
        !>@return mainlayer_id
        !> cardinal coordinate identifying the main layer to which
        !> the grid-point belongs
        !--------------------------------------------------------------
        function get_mainlayer_coord(gen_coords)
     $       result(mainlayer_id)
        
          implicit none

          integer(ikind), dimension(2), intent(in) :: gen_coords
          integer                                  :: mainlayer_id
        

          if(gen_coords(2).le.align_S) then
             mainlayer_id = S
          else

             if(gen_coords(2).ge.align_N) then
                mainlayer_id = N
             else

                if(gen_coords(1).le.align_W) then
                   mainlayer_id = W
                else

                   if(gen_coords(1).ge.align_E) then
                      mainlayer_id = E
                   else

                      mainlayer_id = interior
                   end if

                end if
             end if
          end if

        end function get_mainlayer_coord

      end module bf_increase_coords_module

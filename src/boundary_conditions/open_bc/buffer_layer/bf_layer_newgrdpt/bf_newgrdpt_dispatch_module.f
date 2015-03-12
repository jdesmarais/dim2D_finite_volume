      !determine whether it is possible to evaluate
      !the new grdpt using the data from the buffer
      !layer and compute it if it is possible, otherwise
      !ask to compute the new grid-point from the main
      !structure (interior+buffer layer) to gather data
      module bf_newgrdpt_dispatch_module

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use bf_newgrdpt_verification_module, only :
     $       are_grdpts_available,
     $       get_newgrdpt_verification_bounds
        
        use parameters_constant, only :
     $       N,S,E,W

        use parameters_kind, only :
     $       ikind


        private
        public ::
     $       are_grdpts_available_to_get_newgrdpt_data,
     $       are_grdpts_available_to_get_newgrdpt_proc


        contains


        subroutine are_grdpts_available_to_get_newgrdpt_data(
     $     bf_grdpts_id0,
     $     bf_newgrdpt_coords0,
     $     procedure_type,
     $     gradient_type,
     $     data_needed_bounds0,
     $     tmp_array_needed,
     $     grdpts_available)

          implicit none
          
          integer       , dimension(2,2), intent(in)    :: bf_grdpts_id0
          integer(ikind), dimension(2)  , intent(in)    :: bf_newgrdpt_coords0
          integer                       , intent(in)    :: procedure_type
          integer                       , intent(in)    :: gradient_type
          integer       , dimension(2,2), intent(inout) :: data_needed_bounds0
          logical                       , intent(out)   :: tmp_array_needed
          logical                       , intent(out)   :: grdpts_available


          integer(ikind)                              :: size_x
          integer(ikind)                              :: size_y
          integer                                     :: nb_bounds
          integer       , dimension(2,2,2)            :: bounds
          integer                                     :: k
          logical                                     :: all_grdpts_exists

          
          !> estimate the bounds for the grid points to be checked
          call get_newgrdpt_verification_bounds(
     $         procedure_type,
     $         gradient_type,
     $         nb_bounds,
     $         bounds)


          !> update the bounds for the data needed
          do k=1, nb_bounds
             data_needed_bounds0(1,1) = min(data_needed_bounds0(1,1),bounds(1,1,k))
             data_needed_bounds0(1,2) = max(data_needed_bounds0(1,2),bounds(1,2,k))
             data_needed_bounds0(2,1) = min(data_needed_bounds0(2,1),bounds(2,1,k))
             data_needed_bounds0(2,2) = max(data_needed_bounds0(2,2),bounds(2,2,k))
          end do


          !> determine whether an intermediate array is needed to
          !> evaluate the data available
          tmp_array_needed = .false.

          size_x = size(bf_grdpts_id0,1)
          size_y = size(bf_grdpts_id0,2)

          do k=1, nb_bounds
             tmp_array_needed = tmp_array_needed.and.(
     $            ((bf_newgrdpt_coords0(1)+bounds(1,1,k)).lt.1).or.
     $            ((bf_newgrdpt_coords0(1)+bounds(1,2,k)).gt.size_x).or.
     $            ((bf_newgrdpt_coords0(2)+bounds(2,1,k)).lt.1).or.
     $            ((bf_newgrdpt_coords0(2)+bounds(2,2,k)).gt.size_y)
     $            )
          end do


          !> if no tmp_array is needed, the data needed for the new grdpt
          !> are directly checked on grdpts_id
          if(.not.tmp_array_needed) then

             all_grdpts_exists = .true.

             do k=1,nb_bounds

                if(.not.are_grdpts_available(
     $               bf_grdpts_id0,
     $               reshape((/
     $                  bf_newgrdpt_coords0(1)+bounds(1,1,k),
     $                  bf_newgrdpt_coords0(2)+bounds(2,1,k),
     $                  bf_newgrdpt_coords0(1)+bounds(1,2,k),
     $                  bf_newgrdpt_coords0(2)+bounds(2,1,k)/),
     $                  (/2,2/))
     $               )) then
                   
                   all_grdpts_exists = .false.

                end if

             end do

             if(.not.all_grdpts_exists) then

                print '(''bf_newgrdpt_dispatch_module'')'
                print '(''are_grdpts_available_to_get_newgrdpt_data'')'
                print '(''the data needed by the newgrdpt are in the buffer layer'')'
                print '(''but it does not seem they are all available'')'
                stop ''

             else
                grdpts_available = .true.
             end if

          else
             grdpts_available = .false.
          end if

        end subroutine are_grdpts_available_to_get_newgrdpt_data


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine whether there are enough grid points to
        !> identify the procedure computing the new grid point
        !
        !> @date
        !> 13_03_2015 - initial version - J.L. Desmarais
        !
        !>@param bf_localization
        !> cardinal coordinate identifying the mainlayer to which the
        !> buffer layer belongs
        !
        !>@param size_x
        !> size along the x-direction of the buffer layer
        !
        !>@param size_y
        !> size along the y-direction of the buffer layer
        !
        !>@param bf_can_exchange_with_neighbor1
        !> logical identifying whether the buffer layer can exchange
        !> grid points with its neighbor1
        !
        !>@param bf_can_exchange_with_neighbor2
        !> logical identifying whether the buffer layer can exchange
        !> grid points with its neighbor2
        !
        !>@param bf_newgrdpt_coords
        !> local coordinates of the new grid point to be computed
        !
        !>@return grdpts_available
        !> logical determining whether there are enough grid points
        !> to identify the procedure computing the new grid point
        !--------------------------------------------------------------
        function are_grdpts_available_to_get_newgrdpt_proc(
     $       bf_localization,
     $       size_x, size_y,
     $       bf_can_exchange_with_neighbor1,
     $       bf_can_exchange_with_neighbor2,
     $       bf_newgrdpt_coords)
     $       result(grdpts_available)
        
          implicit none

          integer                     , intent(in) :: bf_localization
          integer                     , intent(in) :: size_x
          integer                     , intent(in) :: size_y
          logical                     , intent(in) :: bf_can_exchange_with_neighbor1
          logical                     , intent(in) :: bf_can_exchange_with_neighbor2
          integer(ikind), dimension(2), intent(in) :: bf_newgrdpt_coords
          logical                                  :: grdpts_available


          select case(bf_localization)
            case(N)

               !([i-1] or [i+1]) and (j-1)
               grdpts_available = 
     $              ( ((bf_newgrdpt_coords(1)-1).ge.1).or.
     $                ((bf_newgrdpt_coords(1)+1).le.size_x)).and.
     $              ( ((bf_newgrdpt_coords(2)-1).ge.1) )

            case(S)

               !([i-1] or [i+1]) and (j+1)
               grdpts_available = 
     $              ( ((bf_newgrdpt_coords(1)-1).ge.1).or.
     $                ((bf_newgrdpt_coords(1)+1).le.size_x)).and.
     $              ( ((bf_newgrdpt_coords(2)+1).le.size_y) )
               
            case(E,W)

               if(bf_localization.eq.E) then
                  grdpts_available = (bf_newgrdpt_coords(1)-1).ge.1
               else
                  grdpts_available = (bf_newgrdpt_coords(1)+1).le.size_x
               end if

               if(bf_can_exchange_with_neighbor1) then

               !([i+1]) and (j-1,j+1)
                  if(bf_can_exchange_with_neighbor2) then
                     grdpts_available = 
     $                    grdpts_available.and.
     $                    ( (bf_newgrdpt_coords(2)-1).ge.1 ).and.
     $                    ( (bf_newgrdpt_coords(2)+1).le.size_y )

               !([i+1]) and (j-1)
                  else
                     grdpts_available = 
     $                    grdpts_available.and.
     $                    ( (bf_newgrdpt_coords(2)-1).ge.1 )
                  end if

               else

               !([i+1]) and (j+1)
                  if(bf_can_exchange_with_neighbor2) then
                     grdpts_available = 
     $                    grdpts_available.and.
     $                    ( (bf_newgrdpt_coords(2)+1).le.size_y )

               !([i+1]) and ((j-1) or (j+1))
                  else
                     grdpts_available = 
     $                    grdpts_available.and.(
     $                    ( (bf_newgrdpt_coords(2)-1).ge.1 ).or.
     $                    ( (bf_newgrdpt_coords(2)+1).le.size_y ))

                  end if                     

               end if

            case default
               call error_mainlayer_id(
     $              'bf_compute_newgrdpt_class',
     $              'are_grdpts_available_to_get_newgrdpt_proc',
     $              bf_localization)

          end select

        end function are_grdpts_available_to_get_newgrdpt_proc

      end module bf_newgrdpt_dispatch_module

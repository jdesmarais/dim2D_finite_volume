      program test_sd_scheme_pattern

        use parameters_bf_layer, only :
     $       interior_pt

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use sd_operators_x_oneside_L0_class, only :
     $       sd_operators_x_oneside_L0

        use sd_operators_x_oneside_L1_class, only :
     $       sd_operators_x_oneside_L1

        use sd_operators_x_oneside_R1_class, only :
     $       sd_operators_x_oneside_R1

        use sd_operators_x_oneside_R0_class, only :
     $       sd_operators_x_oneside_R0

        use sd_operators_y_oneside_L0_class, only :
     $       sd_operators_y_oneside_L0

        use sd_operators_y_oneside_L1_class, only :
     $       sd_operators_y_oneside_L1

        use sd_operators_y_oneside_R1_class, only :
     $       sd_operators_y_oneside_R1

        use sd_operators_y_oneside_R0_class, only :
     $       sd_operators_y_oneside_R0

        use sd_operators_n1_oneside_L0_class, only :
     $       sd_operators_n1_oneside_L0

        use sd_operators_n1_oneside_L1_class, only :
     $       sd_operators_n1_oneside_L1

        use sd_operators_n1_oneside_R1_class, only :
     $       sd_operators_n1_oneside_R1

        use sd_operators_n1_oneside_R0_class, only :
     $       sd_operators_n1_oneside_R0

        use sd_operators_n2_oneside_L0_class, only :
     $       sd_operators_n2_oneside_L0

        use sd_operators_n2_oneside_L1_class, only :
     $       sd_operators_n2_oneside_L1

        use sd_operators_n2_oneside_R1_class, only :
     $       sd_operators_n2_oneside_R1

        use sd_operators_n2_oneside_R0_class, only :
     $       sd_operators_n2_oneside_R0


        implicit none        

        integer, dimension(5,5)          :: grdpts_interior
                                         
        integer, dimension(5,5)          :: grdpts_x_L0
        integer, dimension(5,5)          :: grdpts_x_L1
        integer, dimension(5,5)          :: grdpts_x_R1
        integer, dimension(5,5)          :: grdpts_x_R0
                                         
        integer, dimension(5,5)          :: grdpts_y_L0
        integer, dimension(5,5)          :: grdpts_y_L1
        integer, dimension(5,5)          :: grdpts_y_R1
        integer, dimension(5,5)          :: grdpts_y_R0
                                         
        integer, dimension(5,5)          :: grdpts_n1_L0
        integer, dimension(5,5)          :: grdpts_n1_L1
        integer, dimension(5,5)          :: grdpts_n1_R1
        integer, dimension(5,5)          :: grdpts_n1_R0
                                         
        integer, dimension(5,5)          :: grdpts_n2_L0
        integer, dimension(5,5)          :: grdpts_n2_L1
        integer, dimension(5,5)          :: grdpts_n2_R1
        integer, dimension(5,5)          :: grdpts_n2_R0
                                         
                                         
        type(sd_operators)               :: sd_interior
                                         
        type(sd_operators_x_oneside_L0)  :: sd_x_L0
        type(sd_operators_x_oneside_L1)  :: sd_x_L1
        type(sd_operators_x_oneside_R1)  :: sd_x_R1
        type(sd_operators_x_oneside_R0)  :: sd_x_R0
                                         
        type(sd_operators_y_oneside_L0)  :: sd_y_L0
        type(sd_operators_y_oneside_L1)  :: sd_y_L1
        type(sd_operators_y_oneside_R1)  :: sd_y_R1
        type(sd_operators_y_oneside_R0)  :: sd_y_R0

        type(sd_operators_n1_oneside_L0) :: sd_n1_L0
        type(sd_operators_n1_oneside_L1) :: sd_n1_L1
        type(sd_operators_n1_oneside_R1) :: sd_n1_R1
        type(sd_operators_n1_oneside_R0) :: sd_n1_R0

        type(sd_operators_n2_oneside_L0) :: sd_n2_L0
        type(sd_operators_n2_oneside_L1) :: sd_n2_L1
        type(sd_operators_n2_oneside_R1) :: sd_n2_R1
        type(sd_operators_n2_oneside_R0) :: sd_n2_R0


        !for each space discretization scheme, determine
        !which gridpoints are used for the computation of
        !the fluxes
        grdpts_interior   = get_grdpts_int(sd_interior)
                          
        grdpts_x_L0       = get_grdpts_x_oneside(sd_x_L0)
        grdpts_x_L1       = get_grdpts_x_oneside(sd_x_L1)
        grdpts_x_R1       = get_grdpts_x_oneside(sd_x_R1)
        grdpts_x_R0       = get_grdpts_x_oneside(sd_x_R0)
                          
        grdpts_y_L0       = get_grdpts_y_oneside(sd_y_L0)
        grdpts_y_L1       = get_grdpts_y_oneside(sd_y_L1)
        grdpts_y_R1       = get_grdpts_y_oneside(sd_y_R1)
        grdpts_y_R0       = get_grdpts_y_oneside(sd_y_R0)

        grdpts_n1_L0      = get_grdpts_n1_oneside(sd_n1_L0)
        grdpts_n1_L1      = get_grdpts_n1_oneside(sd_n1_L1)
        grdpts_n1_R1      = get_grdpts_n1_oneside(sd_n1_R1)
        grdpts_n1_R0      = get_grdpts_n1_oneside(sd_n1_R0)
                          
        grdpts_n2_L0      = get_grdpts_n2_oneside(sd_n2_L0)
        grdpts_n2_L1      = get_grdpts_n2_oneside(sd_n2_L1)
        grdpts_n2_R1      = get_grdpts_n2_oneside(sd_n2_R1)
        grdpts_n2_R0      = get_grdpts_n2_oneside(sd_n2_R0)

        
        call print_grdpts_used('sd_int' ,grdpts_interior)

        call print_grdpts_used('sd_x_L0',grdpts_x_L0)
        call print_grdpts_used('sd_x_L1',grdpts_x_L1)
        call print_grdpts_used('sd_x_R1',grdpts_x_R1)
        call print_grdpts_used('sd_x_R0',grdpts_x_R0)

        call print_grdpts_used('sd_y_L0',grdpts_y_L0)
        call print_grdpts_used('sd_y_L1',grdpts_y_L1)
        call print_grdpts_used('sd_y_R1',grdpts_y_R1)
        call print_grdpts_used('sd_y_R0',grdpts_y_R0)


        call print_grdpts_used('sd_n1_L0',grdpts_n1_L0)
        call print_grdpts_used('sd_n1_L1',grdpts_n1_L1)
        call print_grdpts_used('sd_n1_R1',grdpts_n1_R1)
        call print_grdpts_used('sd_n1_R0',grdpts_n1_R0)

        call print_grdpts_used('sd_n2_L0',grdpts_n2_L0)
        call print_grdpts_used('sd_n2_L1',grdpts_n2_L1)
        call print_grdpts_used('sd_n2_R1',grdpts_n2_R1)
        call print_grdpts_used('sd_n2_R0',grdpts_n2_R0)


        contains


        function get_grdpts_int(sd_used)
     $       result(grdpts)

          implicit none

          class(sd_operators), intent(in) :: sd_used
          integer, dimension(5,5)         :: grdpts
          integer, dimension(5,5)         :: grdpts_id

          real(rkind)                     :: dx
          real(rkind)                     :: dy
          real(rkind), dimension(5)       :: x_map
          real(rkind), dimension(5)       :: y_map
          real(rkind), dimension(5,5,ne)  :: nodes

          type(pmodel_eq)                 :: p_model

          real(rkind), dimension(6,5,ne)  :: flux_x
          real(rkind), dimension(5,6,ne)  :: flux_y  

          logical :: grdpt_needed

          integer(ikind) :: i,j
          integer        :: k

          call ini_grdpts(grdpts)
          call ini_grdpts_id(grdpts_id)

          do j=1,5
             do i=1,5

                !initialize the x_map and y_map
                call initialize_maps(dx,dy,x_map,y_map)

                
                !initialize the nodes with one wrong grdpt
                call initialize_nodes(nodes)

                do k=1,ne
                   nodes(i,j,k) = 1e20
                end do
                

                !compute the x-fluxes and the y-fluxes needed
                !to compute the time derivatives at (3,3)
                call p_model%compute_flux_x_nopt(
     $               nodes,dx,dy,sd_used,
     $               grdpts_id,
     $               flux_x,
     $               [3,3],[3,3])

                call p_model%compute_flux_y_nopt(
     $               nodes,dx,dy,sd_used,
     $               grdpts_id,
     $               flux_y,
     $               [3,3],[3,3])


                !test if there is one flux which could not be
                !computed such that the nodes(i,j,:) were needed
                !for the computation of the fluxes
                grdpt_needed = check_if_grdpt_is_needed(flux_x,flux_y)
                
                if(grdpt_needed) then
                   grdpts(i,j)=1
                end if
                
             end do
          end do
                

        end function get_grdpts_int


        function get_grdpts_x_oneside(sd_used)
     $       result(grdpts)

          implicit none

          class(sd_operators), intent(in) :: sd_used
          integer, dimension(5,5)         :: grdpts

          real(rkind)                     :: dx
          real(rkind)                     :: dy
          real(rkind), dimension(5)       :: x_map
          real(rkind), dimension(5)       :: y_map
          real(rkind), dimension(5,5,ne)  :: nodes

          type(pmodel_eq)                 :: p_model

          real(rkind), dimension(6,5,ne)  :: flux_x
          real(rkind), dimension(5,6,ne)  :: flux_y 

          logical :: grdpt_needed

          integer(ikind) :: i,j
          integer        :: k

          integer(ikind) :: i_f,j_f


          call ini_grdpts(grdpts)

          do j=1,5
             do i=1,5
                
                !initialize the x_map and y_map
                call initialize_maps(dx,dy,x_map,y_map)

                
                !initialize the nodes with one wrong grdpt
                call initialize_nodes(nodes)

                do k=1,ne
                   nodes(i,j,k) = 1e20
                end do
                
                i_f = 3
                j_f = 3


                !compute the x-fluxes and the y-fluxes needed
                !to compute the time derivatives at (3,3)
                flux_x(i_f,j_f,:)   = [0,0,0,0]
                flux_x(i_f+1,j_f,:) = [0,0,0,0]

                flux_y(i_f,j_f,:)   = p_model%compute_flux_y_oneside(
     $               nodes,dx,dy,
     $               i_f,j_f,
     $               sd_used)

                flux_y(i_f,j_f+1,:) = p_model%compute_flux_y_oneside(
     $               nodes,dx,dy,
     $               i_f,j_f+1,
     $               sd_used)


                !test if there is one flux which could not be
                !computed such that the nodes(i,j,:) were needed
                !for the computation of the fluxes
                grdpt_needed = check_if_grdpt_is_needed(flux_x,flux_y)

                if(grdpt_needed) then
                   grdpts(i,j)=1
                end if
                
             end do
          end do

        end function get_grdpts_x_oneside


        function get_grdpts_y_oneside(sd_used)
     $       result(grdpts)

          implicit none

          class(sd_operators), intent(in) :: sd_used
          integer, dimension(5,5)         :: grdpts

          real(rkind)                     :: dx
          real(rkind)                     :: dy
          real(rkind), dimension(5)       :: x_map
          real(rkind), dimension(5)       :: y_map
          real(rkind), dimension(5,5,ne)  :: nodes

          type(pmodel_eq)                 :: p_model

          real(rkind), dimension(6,5,ne)  :: flux_x
          real(rkind), dimension(5,6,ne)  :: flux_y 

          logical :: grdpt_needed

          integer(ikind) :: i,j
          integer        :: k

          integer(ikind) :: i_f,j_f


          call ini_grdpts(grdpts)

          do j=1,5
             do i=1,5
                
                !initialize the x_map and y_map
                call initialize_maps(dx,dy,x_map,y_map)


                !initialize the nodes with one wrong grdpt
                call initialize_nodes(nodes)

                do k=1,ne
                   nodes(i,j,k) = 1e20
                end do
                
                i_f = 3
                j_f = 3
                

                !compute the x-fluxes and the y-fluxes needed
                !to compute the time derivatives at (3,3)
                flux_x(i_f,j_f,:)   = p_model%compute_flux_x_oneside(
     $               nodes,dx,dy,
     $               i_f,j_f,
     $               sd_used)

                flux_x(i_f+1,j_f,:) = p_model%compute_flux_x_oneside(
     $               nodes,dx,dy,
     $               i_f+1,j_f,
     $               sd_used)

                flux_y(i_f,j_f,:)   = [0,0,0,0]
                flux_y(i_f,j_f+1,:) = [0,0,0,0]


                !test if there is one flux which could not be
                !computed such that the nodes(i,j,:) were needed
                !for the computation of the fluxes
                grdpt_needed = check_if_grdpt_is_needed(flux_x,flux_y)

                if(grdpt_needed) then
                   grdpts(i,j)=1
                end if
                
             end do
          end do

        end function get_grdpts_y_oneside        


        function get_grdpts_n1_oneside(sd_used)
     $       result(grdpts)

          implicit none

          class(sd_operators), intent(in) :: sd_used
          integer, dimension(5,5)         :: grdpts

          real(rkind)                     :: dx
          real(rkind)                     :: dy
          real(rkind), dimension(5)       :: x_map
          real(rkind), dimension(5)       :: y_map
          real(rkind), dimension(5,5,ne)  :: nodes

          type(pmodel_eq)                 :: p_model

          real(rkind), dimension(6,5,ne)  :: flux_x
          real(rkind), dimension(5,6,ne)  :: flux_y 

          logical :: grdpt_needed

          integer(ikind) :: i,j
          integer        :: k

          integer(ikind) :: i_f,j_f


          call ini_grdpts(grdpts)

          do j=1,5
             do i=1,5
                
                !initialize the x_map and y_map
                call initialize_maps(dx,dy,x_map,y_map)

                
                !initialize the nodes with one wrong grdpt
                call initialize_nodes(nodes)

                do k=1,ne
                   nodes(i,j,k) = 1e20
                end do
                
                i_f = 3
                j_f = 3


                !compute the x-fluxes and the y-fluxes needed
                !to compute the time derivatives at (3,3)
                flux_x(i_f,j_f,:)   = [0,0,0,0]
                flux_x(i_f+1,j_f,:) = [0,0,0,0]

                flux_y(i_f,j_f,:)   = p_model%compute_flux_y_oneside(
     $               nodes,dx,dy,
     $               i_f,j_f,
     $               sd_used)

                flux_y(i_f,j_f+1,:) = p_model%compute_flux_y_oneside(
     $               nodes,dx,dy,
     $               i_f+1,j_f+1,
     $               sd_used)


                !test if there is one flux which could not be
                !computed such that the nodes(i,j,:) were needed
                !for the computation of the fluxes
                grdpt_needed = check_if_grdpt_is_needed(flux_x,flux_y)

                if(grdpt_needed) then
                   grdpts(i,j)=1
                end if
                
             end do
          end do

        end function get_grdpts_n1_oneside


        function get_grdpts_n2_oneside(sd_used)
     $       result(grdpts)

          implicit none

          class(sd_operators), intent(in) :: sd_used
          integer, dimension(5,5)         :: grdpts

          real(rkind)                     :: dx
          real(rkind)                     :: dy
          real(rkind), dimension(5)       :: x_map
          real(rkind), dimension(5)       :: y_map
          real(rkind), dimension(5,5,ne)  :: nodes

          type(pmodel_eq)                 :: p_model

          real(rkind), dimension(6,5,ne)  :: flux_x
          real(rkind), dimension(5,6,ne)  :: flux_y 

          logical :: grdpt_needed

          integer(ikind) :: i,j
          integer        :: k

          integer(ikind) :: i_f,j_f


          call ini_grdpts(grdpts)

          do j=1,5
             do i=1,5
                
                !initialize the x_map and y_map
                call initialize_maps(dx,dy,x_map,y_map)


                !initialize the nodes with one wrong grdpt
                call initialize_nodes(nodes)

                do k=1,ne
                   nodes(i,j,k) = 1e20
                end do
                
                i_f = 3
                j_f = 3
                

                !compute the x-fluxes and the y-fluxes needed
                !to compute the time derivatives at (3,3)
                flux_x(i_f,j_f,:)   = p_model%compute_flux_x_oneside(
     $               nodes,dx,dy,
     $               i_f,j_f,
     $               sd_used)

                flux_x(i_f+1,j_f,:) = p_model%compute_flux_x_oneside(
     $               nodes,dx,dy,
     $               i_f+1,j_f-1,
     $               sd_used)

                flux_y(i_f,j_f,:)   = [0,0,0,0]
                flux_y(i_f,j_f+1,:) = [0,0,0,0]


                !test if there is one flux which could not be
                !computed such that the nodes(i,j,:) were needed
                !for the computation of the fluxes
                grdpt_needed = check_if_grdpt_is_needed(flux_x,flux_y)

                if(grdpt_needed) then
                   grdpts(i,j)=1
                end if
                
             end do
          end do

        end function get_grdpts_n2_oneside


        function check_if_grdpt_is_needed(flux_x,flux_y)
     $     result(grdpt_needed)

          implicit none

          real(rkind), dimension(6,5,ne), intent(in) :: flux_x
          real(rkind), dimension(5,6,ne), intent(in) :: flux_y
          logical                                    :: grdpt_needed


          integer(ikind) :: i,j
          integer        :: k


          grdpt_needed = .false.
          i=3
          do k=1,ne
             do j=3,4
                if ( abs(flux_y(i,j,k)).gt.(1e10)) then
                   grdpt_needed = .true.
                end if
             end do
          end do
                
          j=3
          do k=1,ne
             do i=3,4
                if ( abs(flux_x(i,j,k)).gt.(1e10)) then
                   grdpt_needed = .true.
                end if
             end do
          end do

        end function check_if_grdpt_is_needed


        subroutine initialize_maps(dx,dy,x_map,y_map)

          implicit none

          real(rkind)              , intent(out) :: dx
          real(rkind)              , intent(out) :: dy
          real(rkind), dimension(5), intent(out) :: x_map
          real(rkind), dimension(5), intent(out) :: y_map

          
          dx = 1.0
          dy = 1.0

          x_map = [1.0,2.0,3.0,4.0,5.0]
          y_map = [1.0,2.0,3.0,4.0,5.0]

        end subroutine initialize_maps


        subroutine initialize_nodes(nodes)

          implicit none

          real(rkind), dimension(5,5,ne), intent(out) :: nodes


          nodes = reshape((/
     $         9.75, 9.76, 9.74, 9.81, 9.82,
     $         9.75, 9.73, 9.72, 9.74, 9.81,
     $         9.67, 9.76, 9.56, 9.81, 9.72,
     $         9.70, 9.78, 9.85, 9.73, 9.86,
     $         9.73, 9.79, 9.82, 9.81, 9.72,
     $         
     $         0.15,0.18,0.20,0.12,0.20,0.23,
     $         0.15,0.15,0.18,0.13,0.23,0.23,
     $         0.19,0.23,0.21,0.12,0.18,0.23,
     $         0.14,0.16,0.25,0.15,0.16,0.21,
     $         0.13,0.18,0.27,0.18,0.22,0.23,
     $         
     $         0.3 ,0.26,0.21,0.28,0.29,0.31,
     $         0.26,0.21,0.30,0.22,0.25,0.33,
     $         0.21,0.26,0.21,0.28,0.29,0.31,
     $         0.25,0.22,0.24,0.28,0.35,0.31,
     $         0.28,0.23,0.28,0.32,0.29,0.31,
     $         
     $         7.26,7.20,7.26,7.28,7.12,7.26,
     $         7.26,7.20,7.26,7.28,7.12,7.26,
     $         7.26,7.20,7.26,7.28,7.12,7.26,
     $         7.26,7.20,7.26,7.28,7.12,7.26,
     $         7.26,7.20,7.26,7.28,7.12,7.26/),
     $         (/5,5,ne/))


        end subroutine initialize_nodes

      
        subroutine ini_grdpts_id(grdpts_id)

          implicit none

          integer, dimension(5,5), intent(out) :: grdpts_id

          integer :: i,j

          do j=1,5
             do i=1,5
                grdpts_id(i,j) = interior_pt
             end do
          end do

        end subroutine ini_grdpts_id


        subroutine ini_grdpts(grdpts)

          implicit none

          integer, dimension(5,5), intent(out) :: grdpts

          integer :: i,j

          do j=1,5
             do i=1,5
                grdpts(i,j) = 0
             end do
          end do

        end subroutine ini_grdpts


        subroutine print_grdpts_used(operator,grdpts)

          implicit none

          character*(*)          , intent(in) :: operator
          integer, dimension(5,5), intent(in) :: grdpts

          integer :: j

          print *, operator
          print '(''--------------------'')'
          print '()'

          do j=5,1,-1
             print '(5I2)', grdpts(:,j)
          end do

          print '()'

        end subroutine print_grdpts_used

      end program test_sd_scheme_pattern

      !> @file
      !> determine the grid points needed to compute the new grid point
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> subroutines to determine the grid points needed to compute
      !> the new grid point
      !
      !> @date
      !> 27_02_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_newgrdpt_verification_module

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
     $       no_gradient_type,
     $       gradient_I_type,
     $       gradient_L0_type,
     $       gradient_R0_type,
     $       gradient_xLR0_yI_type,
     $       gradient_xI_yLR0_type,
     $       gradient_xLR0_yLR0_type


        implicit none

        private
        public ::
     $       get_newgrdpt_verification_bounds

        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the procedures needed to compute the new grid points
        !> in association with the configuration of the grid points
        !> around it, given by grdpts_id
        !
        !> @date
        !> 27_02_2015 - initial version - J.L. Desmarais
        !
        !>@param procedure_type
        !> type of procedure used to compute the new grid point
        !
        !>@param gradient_type
        !> type of gradient used to compute the transverse terms
        !> when computing the new grid points
        !
        !>@param nb_bounds
        !> number of loops needed to check the grid points around
        !> the new grid point
        !
        !>@param bounds
        !> loop bounds for checking the grid points around the new
        !> grid point
        !--------------------------------------------------------------
        subroutine get_newgrdpt_verification_bounds(
     $       procedure_type,
     $       gradient_type,
     $       nb_bounds,
     $       bounds)

          implicit none

          integer                  , intent(in)  :: procedure_type
          integer                  , intent(in)  :: gradient_type
          integer                  , intent(out) :: nb_bounds
          integer, dimension(2,2,2), intent(out) :: bounds

          
          select case(procedure_type)

            case(N_edge_type)
               
               nb_bounds     = 1

               bounds(2,1,1) = -2
               bounds(2,2,1) = -1

               select case(gradient_type)
                 case(gradient_I_type)

                    bounds(1,1,1) = -1
                    bounds(1,2,1) =  1

                 case(gradient_L0_type)
                    
                    bounds(1,1,1) =  0
                    bounds(1,2,1) =  1
                    
                 case(gradient_R0_type)
                    
                    bounds(1,1,1) = -1
                    bounds(1,2,1) =  0

                 case default
                    print '(''bf_newgrdpt_verification_module'')'
                    print '(''get_newgrdpt_verification_bounds'')'
                    print '(''for N edge'')'
                    print '(''gradient_type not recognized: '',I2)', gradient_type
                    stop ''

               end select

            case(S_edge_type)
               
               nb_bounds     = 1

               bounds(2,1,1) = 1
               bounds(2,2,1) = 2

               select case(gradient_type)
                 case(gradient_I_type)

                    bounds(1,1,1) = -1
                    bounds(1,2,1) =  1

                 case(gradient_L0_type)
                    
                    bounds(1,1,1) =  0
                    bounds(1,2,1) =  1
                    
                 case(gradient_R0_type)
                    
                    bounds(1,1,1) = -1
                    bounds(1,2,1) =  0

                 case default
                    print '(''bf_newgrdpt_verification_module'')'
                    print '(''get_newgrdpt_verification_bounds'')'
                    print '(''for S edge'')'
                    print '(''gradient_type not recognized: '',I2)', gradient_type
                    stop ''

               end select

            case(E_edge_type)
               
               nb_bounds     = 1

               bounds(1,1,1) = -2
               bounds(1,2,1) = -1

               select case(gradient_type)
                 case(gradient_I_type)

                    bounds(2,1,1) = -1
                    bounds(2,2,1) =  1

                 case(gradient_L0_type)
                    
                    bounds(2,1,1) =  0
                    bounds(2,2,1) =  1
                    
                 case(gradient_R0_type)
                    
                    bounds(2,1,1) = -1
                    bounds(2,2,1) =  0

                 case default
                    print '(''bf_newgrdpt_verification_module'')'
                    print '(''get_newgrdpt_verification_bounds'')'
                    print '(''for E edge'')'
                    print '(''gradient_type not recognized: '',I2)', gradient_type
                    stop ''

               end select

            case(W_edge_type)
               
               nb_bounds     = 1

               bounds(1,1,1) = 1
               bounds(1,2,1) = 2

               select case(gradient_type)
                 case(gradient_I_type)

                    bounds(2,1,1) = -1
                    bounds(2,2,1) =  1

                 case(gradient_L0_type)
                    
                    bounds(2,1,1) =  0
                    bounds(2,2,1) =  1
                    
                 case(gradient_R0_type)
                    
                    bounds(2,1,1) = -1
                    bounds(2,2,1) =  0

                 case default
                    print '(''bf_newgrdpt_verification_module'')'
                    print '(''get_newgrdpt_verification_bounds'')'
                    print '(''for W edge'')'
                    print '(''gradient_type not recognized: '',I2)', gradient_type
                    stop ''

               end select

            case(NW_corner_type)
               
               nb_bounds     =  1

               bounds(1,1,1) =  1
               bounds(1,2,1) =  3
               bounds(2,1,1) = -3
               bounds(2,2,1) = -1

            case(NE_corner_type)
               
               nb_bounds     =  1

               bounds(1,1,1) = -3
               bounds(1,2,1) = -1
               bounds(2,1,1) = -3
               bounds(2,2,1) = -1

            case(SW_corner_type)
               
               nb_bounds     =  1

               bounds(1,1,1) =  1
               bounds(1,2,1) =  3
               bounds(2,1,1) =  1
               bounds(2,2,1) =  3

            case(SE_corner_type)
               
               nb_bounds     =  1

               bounds(1,1,1) = -3
               bounds(1,2,1) = -1
               bounds(2,1,1) =  1
               bounds(2,2,1) =  3

            case(SW_edge_type)
               
               nb_bounds     =  2

               bounds(1,1,1) =  1
               bounds(1,2,1) =  2
               bounds(2,2,1) =  0

               bounds(1,2,2) =  2
               bounds(2,1,2) =  1
               bounds(2,2,2) =  2

               select case(gradient_type)
                 case(gradient_I_type)

                    bounds(2,1,1) = -1
                    bounds(1,1,2) = -1

                 case(gradient_xLR0_yI_type)
                    
                    bounds(2,1,1) = -1
                    bounds(1,1,2) =  0
                    
                 case(gradient_xI_yLR0_type)
                    
                    bounds(2,1,1) =  0
                    bounds(1,1,2) = -1

                 case(gradient_xLR0_yLR0_type)
                    
                    bounds(2,1,1) =  0
                    bounds(1,1,2) =  0

                 case default
                    print '(''bf_newgrdpt_verification_module'')'
                    print '(''get_newgrdpt_verification_bounds'')'
                    print '(''for SW_edge_type'')'
                    print '(''gradient_type not recognized: '',I2)', gradient_type
                    stop ''

               end select

            case(SE_edge_type)
               
               nb_bounds     =  2

               bounds(1,1,1) = -2
               bounds(1,2,1) = -1
               bounds(2,2,1) =  0

               bounds(1,1,2) = -2
               bounds(2,1,2) =  1
               bounds(2,2,2) =  2

               select case(gradient_type)
                 case(gradient_I_type)

                    bounds(2,1,1) = -1
                    bounds(1,2,2) =  1

                 case(gradient_xLR0_yI_type)
                    
                    bounds(2,1,1) = -1
                    bounds(1,2,2) =  0
                    
                 case(gradient_xI_yLR0_type)
                    
                    bounds(2,1,1) =  0
                    bounds(1,2,2) =  1

                 case(gradient_xLR0_yLR0_type)
                    
                    bounds(2,1,1) =  0
                    bounds(1,2,2) =  0

                 case default
                    print '(''bf_newgrdpt_verification_module'')'
                    print '(''get_newgrdpt_verification_bounds'')'
                    print '(''for SE_edge_type'')'
                    print '(''gradient_type not recognized: '',I2)', gradient_type
                    stop ''

               end select

            case(NE_edge_type)
               
               nb_bounds     =  2

               bounds(1,1,1) = -2
               bounds(2,1,1) = -2
               bounds(2,2,1) = -1

               bounds(1,1,2) = -2
               bounds(1,2,2) = -1
               bounds(2,1,2) =  0

               select case(gradient_type)
                 case(gradient_I_type)

                    bounds(1,2,1) =  1
                    bounds(2,2,2) =  1

                 case(gradient_xLR0_yI_type)
                    
                    bounds(1,2,1) =  0
                    bounds(2,2,2) =  1
                    
                 case(gradient_xI_yLR0_type)
                    
                    bounds(1,2,1) =  1
                    bounds(2,2,2) =  0

                 case(gradient_xLR0_yLR0_type)
                    
                    bounds(1,2,1) =  0
                    bounds(2,2,2) =  0

                 case default
                    print '(''bf_newgrdpt_verification_module'')'
                    print '(''get_newgrdpt_verification_bounds'')'
                    print '(''for NE_edge_type'')'
                    print '(''gradient_type not recognized: '',I2)', gradient_type
                    stop ''

               end select

            case(NW_edge_type)
               
               nb_bounds     =  2

               bounds(1,2,1) =  2
               bounds(2,1,1) = -2
               bounds(2,2,1) = -1

               bounds(1,1,2) =  1
               bounds(1,2,2) =  2
               bounds(2,1,2) =  0

               select case(gradient_type)
                 case(gradient_I_type)

                    bounds(1,1,1) = -1
                    bounds(2,2,2) =  1

                 case(gradient_xLR0_yI_type)
                    
                    bounds(1,1,1) =  0
                    bounds(2,2,2) =  1
                    
                 case(gradient_xI_yLR0_type)
                    
                    bounds(1,1,1) = -1
                    bounds(2,2,2) =  0

                 case(gradient_xLR0_yLR0_type)
                    
                    bounds(1,1,1) =  0
                    bounds(2,2,2) =  0

                 case default
                    print '(''bf_newgrdpt_verification_module'')'
                    print '(''get_newgrdpt_verification_bounds'')'
                    print '(''for NW_edge_type'')'
                    print '(''gradient_type not recognized: '',I2)', gradient_type
                    stop ''

               end select

            case default
               print '(''bf_newgrdpt_verification_module'')'
               print '(''get_newgrdpt_verification_bounds'')'
               print '(''procedure_type not recognized: '',I2)', procedure_type
               stop ''

          end select

        end subroutine get_newgrdpt_verification_bounds

      end module bf_newgrdpt_verification_module

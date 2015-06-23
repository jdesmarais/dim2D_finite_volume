      !> @file
      !> identification of the procedure needed to compute the new
      !> grid points
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module with the subroutines identifying the procedure
      !> to compute the new grid points
      !
      !> @date
      ! 27_02_2015 - initial version - J.L. Desmarais
      !----------------------------------------------------------------
      module bf_newgrdpt_procedure_module

        use parameters_bf_layer, only :
     $       bc_pt,
     $       no_pt,     
     $       no_gradient_type,
     $       no_gradient_type,
     $       gradient_I_type,
     $       gradient_L0_type,
     $       gradient_R0_type,
     $       gradient_xLR0_yI_type,
     $       gradient_xI_yLR0_type,
     $       gradient_xLR0_yLR0_type,
     $       BF_SUCCESS

        use parameters_constant, only :
     $       no_bc_procedure_type,
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
     $       SW_corner_type

        use parameters_kind, only :
     $       ikind

        implicit none


        private
        public ::
     $       get_newgrdpt_procedure,
     $       get_config_id,
     $       get_config_procedures


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
        !>@param i
        !> x-index identifying the position of the new grid point
        !> in grdpts_id
        !
        !>@param j
        !> y-index identifying the position of the new grid point
        !> in grdpts_id
        !
        !>@param grdpts_id
        !> configuration of the grid points around the new grid point
        !
        !>@param nb_procedures
        !> number of procedures needed to compute the new grid point
        !
        !>@param procedure_type
        !> integer identifying the procedure applied to compute the
        !> newgrdpt
        !
        !>@param gradient_type
        !> integer identifying the procedure applied to compute the
        !> transverse terms when computing the new grid point
        !
        !>@return ierror
        !> integer identifying whether the identificaqtion was
        !> successful
        !--------------------------------------------------------------
        function get_newgrdpt_procedure(
     $       i,j,
     $       grdpts_id,
     $       nb_procedures,
     $       procedure_type,
     $       gradient_type)
     $       result(ierror)

          implicit none

          integer(ikind)                , intent(in)  :: i
          integer(ikind)                , intent(in)  :: j
          integer       , dimension(:,:), intent(in)  :: grdpts_id
          integer                       , intent(out) :: nb_procedures
          integer       , dimension(4)  , intent(out) :: procedure_type
          integer       , dimension(4)  , intent(out) :: gradient_type
          logical                                     :: ierror

          integer :: config_id

          ierror = BF_SUCCESS

          !the configuration of the bc_pt and no_pt
          !arund the new grid point determines an integer
          !identifying the procedure and gradient applied
          !for the computation
          config_id = get_config_id(i,j,grdpts_id)

          !the procedure and gradient are extracted
          call get_config_procedures(
     $         config_id,
     $         nb_procedures,
     $         procedure_type,
     $         gradient_type,
     $         ierr=ierror)

        end function get_newgrdpt_procedure


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the integer config_id identifying the
        !> configuration around the new grid point
        !
        !> @date
        !> 27_02_2015 - initial version - J.L. Desmarais
        !
        !>@param i
        !> x-index identifying the position of the new grid point
        !> in grdpts_id
        !
        !>@param j
        !> y-index identifying the position of the new grid point
        !> in grdpts_id
        !
        !>@param grdpts_id
        !> configuration of the grid points around the new grid point
        !
        !>@return config_id
        !> integer identifying the configuration around the new
        !> grid point
        !--------------------------------------------------------------
        function get_config_id(i,j,grdpts_id) result(config_id)

          implicit none

          integer(ikind)                , intent(in) :: i
          integer(ikind)                , intent(in) :: j
          integer       , dimension(:,:), intent(in) :: grdpts_id
          integer                                    :: config_id

          integer, dimension(3,3) :: grdpts_id_tmp
          integer, dimension(2,2) :: send
          integer, dimension(2,2) :: recv
          integer, dimension(8)   :: analyse_tmp
          integer                 :: k

          ! we extract the configuration of
          ! the grid points around the new
          ! grid point (i,j)
          !   ___________
          !  | 6   7   8 |
          !  | 4 (i,j) 5 |
          !  |_1___2___3_|
          !
          ! non-existing grid points are set
          ! to no_pt by default
          grdpts_id_tmp = reshape((/
     $         no_pt, no_pt, no_pt, 
     $         no_pt, no_pt, no_pt, 
     $         no_pt, no_pt, no_pt/),
     $         (/3,3/))

          send(1,1) = max(i-1,1)
          send(1,2) = min(i+1,size(grdpts_id,1))
          send(2,1) = max(j-1,1)
          send(2,2) = min(j+1,size(grdpts_id,2))

          recv(1,1)  = send(1,1) - i + 2
          recv(1,2)  = send(1,2) - i + 2
          recv(2,1)  = send(2,1) - j + 2
          recv(2,2)  = send(2,2) - j + 2

          grdpts_id_tmp(
     $         recv(1,1):recv(1,2),
     $         recv(2,1):recv(2,2)) = 
     $    grdpts_id(
     $         send(1,1):send(1,2),
     $         send(2,1):send(2,2))

          ! the configuration of the grid points
          ! around the new grid point (i,j) are
          ! gathered in a continuous table
          !  _________________
          ! |_1_2_3_4_5_6_7_8_|
          analyse_tmp(1:3) = grdpts_id_tmp(1:3,1)
          analyse_tmp(4)   = grdpts_id_tmp(1,2)
          analyse_tmp(5)   = grdpts_id_tmp(3,2)
          analyse_tmp(6:8) = grdpts_id_tmp(1:3,3)

          config_id = 0

          do k=1,8
             
             select case(analyse_tmp(k))
               case(no_pt)
                  config_id = config_id
               case(bc_pt)
                  config_id = config_id + 2**(k-1)
               case default
                  print '(''bf_newgrdpt_procedure_module'')'
                  print '(''get_config_id'')'
                  print '(''gridpt_id not recognized'')'
                  print '(''analyse_tmp('',I1,''): '',I2)',
     $                 k,
     $                 analyse_tmp(k)
                  stop ''
             end select

          end do

        end function get_config_id


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the procedures needed to compute the new grid points
        !> in association with the configuration of the grid points
        !> around it, expressed by an integer, config_id
        !
        !> @date
        !> 27_02_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param config_id
        !> identification of the configuration of the grid points
        !> around the new grid point (0-255)
        !
        !>@param nb_procedures
        !> number of procedures needed to compute the new grid point
        !
        !>@param procedure_type
        !> integer identifying the procedure applied to compute the
        !> newgrdpt
        !
        !>@param gradient_type
        !> integer identifying the procedure applied to compute the
        !> transverse terms when computing the new grid point
        !--------------------------------------------------------------
        subroutine get_config_procedures(
     $     config_id,
     $     nb_procedures,
     $     procedure_type,
     $     gradient_type,
     $     ierr)

          implicit none

          integer                , intent(in)  :: config_id
          integer                , intent(out) :: nb_procedures
          integer, dimension(4)  , intent(out) :: procedure_type
          integer, dimension(4)  , intent(out) :: gradient_type
          logical, optional      , intent(out) :: ierr
          
          logical :: ierr_op

          ierr_op = BF_SUCCESS

          !to prevent warnings related to non-initialization
          !of output variables
          nb_procedures     = 0
          procedure_type(1) = no_bc_procedure_type
          gradient_type(1)  = no_gradient_type

          select case(config_id)
            case(1)
               nb_procedures     = 1
               procedure_type(1) = NE_corner_type
            case(3)
               nb_procedures     = 1
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_R0_type
            case(4)
               nb_procedures     = 1
               procedure_type(1) = NW_corner_type
            case(5)
               nb_procedures     = 2
               procedure_type(1) = NE_corner_type
               procedure_type(2) = NW_corner_type
            case(6)
               nb_procedures     = 1
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_L0_type
            case(7)
               nb_procedures     = 1
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_I_type
            case(9)
               nb_procedures     = 1
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_R0_type
            case(11)
               nb_procedures     = 1
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xLR0_yLR0_type
            case(13)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = NW_corner_type
            case(15)
               nb_procedures     = 1
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xI_yLR0_type
            case(20)
               nb_procedures     = 1
               procedure_type(1) = W_edge_type
               gradient_type(1)  = gradient_R0_type
            case(21)
               nb_procedures     = 2
               procedure_type(1) = NE_corner_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_R0_type
            case(22)
               nb_procedures     = 1
               procedure_type(1) = NW_edge_type
               gradient_type(1)  = gradient_xLR0_yLR0_type
            case(23)
               nb_procedures     = 1
               procedure_type(1) = NW_edge_type
               gradient_type(1)  = gradient_xI_yLR0_type
            case(29)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_R0_type
            case(31)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xI_yLR0_type
               procedure_type(2) = NW_edge_type
               gradient_type(2)  = gradient_xI_yLR0_type
            case(32)
               nb_procedures     = 1
               procedure_type(1) = SE_corner_type
            case(33)
               nb_procedures     = 2
               procedure_type(1) = NE_corner_type
               procedure_type(2) = SE_corner_type
            case(35)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = SE_corner_type
            case(36)
               nb_procedures     = 2
               procedure_type(1) = NW_corner_type
               procedure_type(2) = SE_corner_type
            case(37)
               nb_procedures     = 3
               procedure_type(1) = NE_corner_type
               procedure_type(2) = NW_corner_type
               procedure_type(3) = SE_corner_type
            case(38)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = SE_corner_type
            case(39)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = SE_corner_type
            case(40)
               nb_procedures     = 1
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_L0_type
            case(41)
               nb_procedures     = 1
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_I_type
            case(43)
               nb_procedures     = 1
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xLR0_yI_type
            case(44)
               nb_procedures     = 2
               procedure_type(1) = NW_corner_type
               procedure_type(2) = E_edge_type
               gradient_type(2)  = gradient_L0_type
            case(45)
               nb_procedures     = 2
               procedure_type(1) = NW_corner_type
               procedure_type(2) = E_edge_type
               gradient_type(2)  = gradient_I_type
            case(46)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = E_edge_type
               gradient_type(2)  = gradient_L0_type
            case(47)
               nb_procedures     = 1
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_I_type
            case(52)
               nb_procedures     = 2
               procedure_type(1) = W_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = SE_corner_type
            case(53)
               nb_procedures     = 3
               procedure_type(1) = NE_corner_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_R0_type
               procedure_type(3) = SE_corner_type
            case(54)
               nb_procedures     = 2
               procedure_type(1) = NW_edge_type
               gradient_type(1)  = gradient_xLR0_yLR0_type
               procedure_type(2) = SE_corner_type
            case(55)
               nb_procedures     = 2
               procedure_type(1) = NW_edge_type
               gradient_type(1)  = gradient_xI_yLR0_type
               procedure_type(2) = SE_corner_type
            case(60)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_R0_type
            case(61)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_R0_type
            case(62)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = NW_edge_type
               gradient_type(2)  = gradient_xLR0_yLR0_type
            case(63)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = NW_edge_type
               gradient_type(2)  = gradient_xI_yLR0_type
            case(96)
               nb_procedures     = 1
               procedure_type(1) = S_edge_type
               gradient_type(1)  = gradient_R0_type
            case(97)
               nb_procedures     = 2
               procedure_type(1) = NE_corner_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_R0_type
            case(99)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_R0_type
            case(100)
               nb_procedures     = 2
               procedure_type(1) = NW_corner_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_R0_type
            case(101)
               nb_procedures     = 3
               procedure_type(1) = NE_corner_type
               procedure_type(2) = NW_corner_type
               procedure_type(3) = S_edge_type
               gradient_type(3)  = gradient_R0_type
            case(102)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_R0_type
            case(103)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_R0_type
            case(104)
               nb_procedures     = 1
               procedure_type(1) = SE_edge_type
               gradient_type(1)  = gradient_xLR0_yLR0_type
            case(105)
               nb_procedures     = 1
               procedure_type(1) = SE_edge_type
               gradient_type(1)  = gradient_xLR0_yI_type
            case(107)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xLR0_yI_type
               procedure_type(2) = SE_edge_type
               gradient_type(2)  = gradient_xLR0_yI_type
            case(108)
               nb_procedures     = 2
               procedure_type(1) = NW_corner_type
               procedure_type(2) = SE_edge_type
               gradient_type(2)  = gradient_xLR0_yLR0_type
            case(109)
               nb_procedures     = 2
               procedure_type(1) = NW_corner_type
               procedure_type(2) = SE_edge_type
               gradient_type(2)  = gradient_xLR0_yI_type
            case(110)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = SE_edge_type
               gradient_type(2)  = gradient_xLR0_yLR0_type
            case(111)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = SE_edge_type
               gradient_type(2)  = gradient_xLR0_yI_type
            case(116)
               nb_procedures     = 2
               procedure_type(1) = W_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_R0_type
            case(117)
               nb_procedures     = 3
               procedure_type(1) = NE_corner_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_R0_type
               procedure_type(3) = S_edge_type
               gradient_type(3)  = gradient_R0_type
            case(118)
               nb_procedures     = 2
               procedure_type(1) = NW_edge_type
               gradient_type(1)  = gradient_xLR0_yLR0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_R0_type
            case(119)
               nb_procedures     = 2
               procedure_type(1) = NW_edge_type
               gradient_type(1)  = gradient_xI_yLR0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_R0_type
            case(124)
               nb_procedures     = 2
               procedure_type(1) = W_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = SE_edge_type
               gradient_type(2)  = gradient_xLR0_yLR0_type
            case(125)
               nb_procedures     = 2
               procedure_type(1) = W_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = SE_edge_type
               gradient_type(2)  = gradient_xLR0_yI_type
            case(126)
               nb_procedures     = 2
               procedure_type(1) = NW_edge_type
               gradient_type(1)  = gradient_xLR0_yLR0_type
               procedure_type(2) = SE_edge_type
               gradient_type(2)  = gradient_xLR0_yLR0_type
            case(127)
               nb_procedures     = 3
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = NW_edge_type
               gradient_type(2)  = gradient_xI_yLR0_type
               procedure_type(3) = SE_edge_type
               gradient_type(3)  = gradient_xLR0_yI_type
            case(128)
               nb_procedures     = 1
               procedure_type(1) = SW_corner_type
            case(129)
               nb_procedures     = 2
               procedure_type(1) = NE_corner_type
               procedure_type(2) = SW_corner_type
            case(131)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = SW_corner_type
            case(132)
               nb_procedures     = 2
               procedure_type(1) = NW_corner_type
               procedure_type(2) = SW_corner_type
            case(133)
               nb_procedures     = 3
               procedure_type(1) = NE_corner_type
               procedure_type(2) = NW_corner_type
               procedure_type(3) = SW_corner_type
            case(134)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = SW_corner_type
            case(135)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = SW_corner_type
            case(137)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = SW_corner_type
            case(139)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xLR0_yLR0_type
               procedure_type(2) = SW_corner_type
            case(141)
               nb_procedures     = 3
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = NW_corner_type
               procedure_type(3) = SW_corner_type
            case(143)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xI_yLR0_type
               procedure_type(2) = SW_corner_type
            case(144)
               nb_procedures     = 1
               procedure_type(1) = W_edge_type
               gradient_type(1)  = gradient_L0_type
            case(145)
               nb_procedures     = 2
               procedure_type(1) = NE_corner_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_L0_type
            case(147)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_L0_type
            case(148)
               nb_procedures     = 1
               procedure_type(1) = W_edge_type
               gradient_type(1)  = gradient_I_type
            case(149)
               nb_procedures     = 2
               procedure_type(1) = NE_corner_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_I_type
            case(150)
               nb_procedures     = 1
               procedure_type(1) = NW_edge_type
               gradient_type(1)  = gradient_xLR0_yI_type
            case(151)
               nb_procedures     = 1
               procedure_type(1) = NW_edge_type
               gradient_type(1)  = gradient_I_type
            case(153)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_L0_type
            case(155)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xLR0_yLR0_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_L0_type
            case(157)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_I_type
            case(159)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xI_yLR0_type
               procedure_type(2) = NW_edge_type
               gradient_type(2)  = gradient_I_type
            case(160)
               nb_procedures     = 2
               procedure_type(1) = SE_corner_type
               procedure_type(2) = SW_corner_type
            case(161)
               nb_procedures     = 3
               procedure_type(1) = NE_corner_type
               procedure_type(2) = SE_corner_type
               procedure_type(3) = SW_corner_type
            case(163)
               nb_procedures     = 3
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = SE_corner_type
               procedure_type(3) = SW_corner_type
            case(164)
               nb_procedures     = 3
               procedure_type(1) = NW_corner_type
               procedure_type(2) = SE_corner_type
               procedure_type(3) = SW_corner_type
            case(165)
               nb_procedures     = 4
               procedure_type(1) = NE_corner_type
               procedure_type(2) = NW_corner_type
               procedure_type(3) = SE_corner_type
               procedure_type(4) = SW_corner_type
            case(166)
               nb_procedures     = 3
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = SE_corner_type
               procedure_type(3) = SW_corner_type
            case(167)
               nb_procedures     = 3
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = SE_corner_type
               procedure_type(3) = SW_corner_type
            case(168)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = SW_corner_type
            case(169)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = SW_corner_type
            case(171)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xLR0_yI_type
               procedure_type(2) = SW_corner_type
            case(172)
               nb_procedures     = 3
               procedure_type(1) = NW_corner_type
               procedure_type(2) = E_edge_type
               gradient_type(2)  = gradient_L0_type
               procedure_type(3) = SW_corner_type
            case(173)
               nb_procedures     = 3
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = NW_corner_type
               procedure_type(3) = SW_corner_type
            case(174)
               nb_procedures     = 3
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = E_edge_type
               gradient_type(2)  = gradient_L0_type
               procedure_type(3) = SW_corner_type
            case(175)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = SW_corner_type
            case(176)
               nb_procedures     = 2
               procedure_type(1) = SE_corner_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_L0_type
            case(177)
               nb_procedures     = 3
               procedure_type(1) = NE_corner_type
               procedure_type(2) = SE_corner_type
               procedure_type(3) = W_edge_type
               gradient_type(3)  = gradient_L0_type
            case(179)
               nb_procedures     = 3
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = SE_corner_type
               procedure_type(3) = W_edge_type
               gradient_type(3)  = gradient_L0_type
            case(180)
               nb_procedures     = 2
               procedure_type(1) = SE_corner_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_I_type
            case(181)
               nb_procedures     = 3
               procedure_type(1) = NE_corner_type
               procedure_type(2) = SE_corner_type
               procedure_type(3) = W_edge_type
               gradient_type(3)  = gradient_I_type
            case(182)
               nb_procedures     = 2
               procedure_type(1) = SE_corner_type
               procedure_type(2) = NW_edge_type
               gradient_type(2)  = gradient_xLR0_yI_type
            case(183)
               nb_procedures     = 2
               procedure_type(1) = SE_corner_type
               procedure_type(2) = NW_edge_type
               gradient_type(2)  = gradient_I_type
            case(184)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_L0_type
            case(185)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_L0_type
            case(187)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xLR0_yI_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_L0_type
            case(188)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_I_type
            case(189)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = W_edge_type
               gradient_type(2)  = gradient_I_type
            case(190)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = NW_edge_type
               gradient_type(2)  = gradient_xLR0_yI_type
            case(191)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = NW_edge_type
               gradient_type(2)  = gradient_I_type
            case(192)
               nb_procedures     = 1
               procedure_type(1) = S_edge_type
               gradient_type(1)  = gradient_L0_type
            case(193)
               nb_procedures     = 2
               procedure_type(1) = NE_corner_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_L0_type
            case(195)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_L0_type
            case(196)
               nb_procedures     = 2
               procedure_type(1) = NW_corner_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_L0_type
            case(197)
               nb_procedures     = 3
               procedure_type(1) = NE_corner_type
               procedure_type(2) = NW_corner_type
               procedure_type(3) = S_edge_type
               gradient_type(3)  = gradient_L0_type
            case(198)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_L0_type
            case(199)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_L0_type
            case(201)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_L0_type
            case(203)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xLR0_yLR0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_L0_type
            case(205)
               nb_procedures     = 3
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = NW_corner_type
               procedure_type(3) = S_edge_type
               gradient_type(3)  = gradient_L0_type
            case(207)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xI_yLR0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_L0_type
            case(208)
               nb_procedures     = 1
               procedure_type(1) = SW_edge_type
               gradient_type(1)  = gradient_xLR0_yLR0_type
            case(209)
               nb_procedures     = 2
               procedure_type(1) = NE_corner_type
               procedure_type(2) = SW_edge_type
               gradient_type(2)  = gradient_xLR0_yLR0_type
            case(211)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = SW_edge_type
               gradient_type(2)  = gradient_xLR0_yLR0_type
            case(212)
               nb_procedures     = 1
               procedure_type(1) = SW_edge_type
               gradient_type(1)  = gradient_xLR0_yI_type
            case(213)
               nb_procedures     = 2
               procedure_type(1) = NE_corner_type
               procedure_type(2) = SW_edge_type
               gradient_type(2)  = gradient_xLR0_yI_type
            case(214)
               nb_procedures     = 2
               procedure_type(1) = NW_edge_type
               gradient_type(1)  = gradient_xLR0_yI_type
               procedure_type(2) = SW_edge_type
               gradient_type(2)  = gradient_xLR0_yI_type
            case(215)
               nb_procedures     = 2
               procedure_type(1) = NW_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = SW_edge_type
               gradient_type(2)  = gradient_xLR0_yI_type
            case(217)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = SW_edge_type
               gradient_type(2)  = gradient_xLR0_yLR0_type
            case(219)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xLR0_yLR0_type
               procedure_type(2) = SW_edge_type
               gradient_type(2)  = gradient_xLR0_yLR0_type
            case(221)
               nb_procedures     = 2
               procedure_type(1) = E_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = SW_edge_type
               gradient_type(2)  = gradient_xLR0_yI_type
            case(223)
               nb_procedures     = 3
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xI_yLR0_type
               procedure_type(2) = NW_edge_type
               gradient_type(2)  = gradient_I_type
               procedure_type(3) = SW_edge_type
               gradient_type(3)  = gradient_xLR0_yI_type
            case(224)
               nb_procedures     = 1
               procedure_type(1) = S_edge_type
               gradient_type(1)  = gradient_I_type
            case(225)
               nb_procedures     = 2
               procedure_type(1) = NE_corner_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_I_type
            case(227)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_R0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_I_type
            case(228)
               nb_procedures     = 2
               procedure_type(1) = NW_corner_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_I_type
            case(229)
               nb_procedures     = 3
               procedure_type(1) = NE_corner_type
               procedure_type(2) = NW_corner_type               
               procedure_type(3) = S_edge_type
               gradient_type(3)  = gradient_I_type
            case(230)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_L0_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_I_type
            case(231)
               nb_procedures     = 2
               procedure_type(1) = N_edge_type
               gradient_type(1)  = gradient_I_type
               procedure_type(2) = S_edge_type
               gradient_type(2)  = gradient_I_type
            case(232)
               nb_procedures     = 1
               procedure_type(1) = SE_edge_type
               gradient_type(1)  = gradient_xI_yLR0_type
            case(233)
               nb_procedures     = 1
               procedure_type(1) = SE_edge_type
               gradient_type(1)  = gradient_I_type
            case(235)
               nb_procedures     = 2
               procedure_type(1) = NE_edge_type
               gradient_type(1)  = gradient_xLR0_yI_type
               procedure_type(2) = SE_edge_type
               gradient_type(2)  = gradient_I_type
            case(236)
               nb_procedures      = 2
               procedure_type(1)  = NW_corner_type
               procedure_type(2)  = SE_edge_type
               gradient_type(2)   = gradient_xI_yLR0_type
            case(237)
               nb_procedures      = 2
               procedure_type(1)  = NW_corner_type
               procedure_type(2)  = SE_edge_type
               gradient_type(2)   = gradient_I_type
            case(238)
               nb_procedures      = 2
               procedure_type(1)  = N_edge_type
               gradient_type(1)   = gradient_L0_type
               procedure_type(2)  = SE_edge_type
               gradient_type(2)   = gradient_xI_yLR0_type
            case(239)
               nb_procedures      = 2
               procedure_type(1)  = NE_edge_type
               gradient_type(1)   = gradient_I_type
               procedure_type(2)  = SE_edge_type
               gradient_type(2)   = gradient_I_type
            case(240)
               nb_procedures      = 1
               procedure_type(1)  = SW_edge_type
               gradient_type(1)   = gradient_xI_yLR0_type
            case(241)
               nb_procedures      = 2
               procedure_type(1)  = NE_corner_type
               procedure_type(2)  = SW_edge_type
               gradient_type(2)   = gradient_xI_yLR0_type
            case(243)
               nb_procedures      = 2
               procedure_type(1)  = N_edge_type
               gradient_type(1)   = gradient_R0_type
               procedure_type(2)  = SW_edge_type
               gradient_type(2)   = gradient_xI_yLR0_type
            case(244)
               nb_procedures      = 1
               procedure_type(1)  = SW_edge_type
               gradient_type(1)   = gradient_I_type
            case(245)
               nb_procedures      = 2
               procedure_type(1)  = NE_corner_type
               procedure_type(2)  = SW_edge_type
               gradient_type(2)   = gradient_I_type
            case(246)
               nb_procedures      = 2
               procedure_type(1)  = NW_edge_type
               gradient_type(1)   = gradient_xLR0_yI_type
               procedure_type(2)  = SW_edge_type
               gradient_type(2)   = gradient_I_type
            case(247)
               nb_procedures      = 2
               procedure_type(1)  = NW_edge_type
               gradient_type(1)   = gradient_I_type
               procedure_type(2)  = SW_edge_type
               gradient_type(2)   = gradient_I_type
            case(248)
               nb_procedures      = 2
               procedure_type(1)  = SE_edge_type
               gradient_type(1)   = gradient_xI_yLR0_type
               procedure_type(2)  = SW_edge_type
               gradient_type(2)   = gradient_xI_yLR0_type
            case(249)
               nb_procedures      = 2
               procedure_type(1)  = SE_edge_type
               gradient_type(1)   = gradient_I_type
               procedure_type(2)  = SW_edge_type
               gradient_type(2)   = gradient_xI_yLR0_type
            case(251)
               nb_procedures      = 3
               procedure_type(1)  = NE_edge_type
               gradient_type(1)   = gradient_xLR0_yI_type
               procedure_type(2)  = SE_edge_type
               gradient_type(2)   = gradient_I_type
               procedure_type(3)  = SW_edge_type
               gradient_type(3)   = gradient_xI_yLR0_type
            case(252)
               nb_procedures      = 2
               procedure_type(1)  = SE_edge_type
               gradient_type(1)   = gradient_xI_yLR0_type
               procedure_type(2)  = SW_edge_type
               gradient_type(2)   = gradient_I_type
            case(253)
               nb_procedures      = 2
               procedure_type(1)  = SE_edge_type
               gradient_type(1)   = gradient_I_type
               procedure_type(2)  = SW_edge_type
               gradient_type(2)   = gradient_I_type
            case(254)
               nb_procedures      = 3
               procedure_type(1)  = NW_edge_type
               gradient_type(1)   = gradient_xLR0_yI_type
               procedure_type(2)  = SE_edge_type
               gradient_type(2)   = gradient_xI_yLR0_type
               procedure_type(3)  = SW_edge_type
               gradient_type(3)   = gradient_I_type
            case(255)
               nb_procedures      = 4
               procedure_type(1)  = NE_edge_type
               gradient_type(1)   = gradient_I_type
               procedure_type(2)  = NW_edge_type
               gradient_type(2)   = gradient_I_type
               procedure_type(3)  = SE_edge_type
               gradient_type(3)   = gradient_I_type
               procedure_type(4)  = SW_edge_type
               gradient_type(4)   = gradient_I_type
            case default
               ierr_op = .not.BF_SUCCESS
c$$$               print '(''bf_newgrdpt_procedure_module'')'
c$$$               print '(''get_config_procedures'')'
c$$$               print '(''configuration not recognized'')'
c$$$               print '(''config_id: '',I3)', config_id
c$$$               print '()'
          end select

          if(present(ierr)) then
             ierr = ierr_op
          end if

        end subroutine get_config_procedures

      end module bf_newgrdpt_procedure_module

      !> @file
      !> generalized boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> object that encapulates the different types of boundary
      !> conditions implemented as relatives of bc_operators_abstract.
      !> The boundary conditions are then applied by choosing according
      !> to the user's choice
      !>
      !> @date
      !> 07_06_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_gen_class

        use bc_operators_hedstrom_xy_class, only :
     $       bc_operators_hedstrom_xy

        use bc_operators_reflection_xy_class, only :
     $       bc_operators_reflection_xy

        use bc_operators_wall_xy_class, only :
     $       bc_operators_wall_xy

        use ISO_FORTRAN_ENV, only :
     $       ERROR_UNIT

        use parameters_bf_layer, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type

        use parameters_constant, only :
     $       reflection_xy_choice,
     $       wall_xy_choice

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        implicit none


        integer, dimension(4), parameter :: SW_corner_section = [SW_corner_type, 1           , 1           , 0           ]
        integer, dimension(4), parameter :: S_edge_section    = [S_edge_type   , bc_size+1   , 1           , nx-bc_size  ]
        integer, dimension(4), parameter :: SE_corner_section = [SE_corner_type, nx-bc_size+1, 1           , 0           ]
        integer, dimension(4), parameter :: W_edge_section    = [W_edge_type   , 1           , bc_size+1   , ny-bc_size+1]
        integer, dimension(4), parameter :: E_edge_section    = [E_edge_type   , nx-bc_size+1, bc_size+1   , ny-bc_size+1]
        integer, dimension(4), parameter :: NW_corner_section = [NW_corner_type, 1           , ny-bc_size+1, 0           ]
        integer, dimension(4), parameter :: N_edge_section    = [N_edge_type   , bc_size+1   , ny-bc_size+1, nx-bc_size  ]
        integer, dimension(4), parameter :: NE_corner_section = [NE_corner_type, nx-bc_size+1, ny-bc_size+1, 0           ]


        private
        public :: bc_operators_gen


        !> @class bc_operators_gen
        !> object that encapulates the different types of boundary
        !> conditions implemented as relatives of bc_operators_abstract.
        !> The boundary conditions are then applied by choosing according
        !> to the user's choice
        !
        !>@param apply_bc_on_nodes
        !> interface to apply the boundary conditions along
        !> the x and y directions on the nodes at the edge of
        !> the computational domain
        !
        !>@param apply_bc_on_fluxes
        !> interface to apply the boundary conditions along
        !> the x and y directions on the fluxes at the edge of
        !> the computational domain
        !
        !>@param apply_bc_on_timedev
        !> interface to apply the boundary conditions along
        !> the x and y directions on the time derivatives at
        !> the edge of the computational domain
        !
        !>@param apply_bc_on_nodes_nopt
        !> interface to apply the boundary conditions along
        !> the x and y directions on the nodes at the edge of
        !> the computational domain for an adaptive computational
        !> domain
        !
        !>@param apply_bc_on_fluxes_nopt
        !> interface to apply the boundary conditions along
        !> the x and y directions on the fluxes at the edge of
        !> the computational domain for an adaptive computational
        !> domain
        !
        !>@param apply_bc_on_timedev_nopt
        !> interface to apply the boundary conditions along
        !> the x and y directions on the time derivatives at
        !> the edge of the computational domain for an adaptive
        !> computational domain
        !---------------------------------------------------------------
        type :: bc_operators_gen

          contains

          procedure, nopass :: apply_bc_on_nodes
          procedure, pass   :: apply_bc_on_fluxes
          procedure, pass   :: apply_bc_on_timedev

          !procedure(nodes_nopt_proc),   pass, deferred :: apply_bc_on_nodes_nopt
          !procedure(tdev_nopt_proc) ,   pass, deferred :: apply_bc_on_timedev_nopt

        end type bc_operators_gen


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> along the x and y directions at the edge of the
        !> computational domain
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> coordinate map along the x-direction
        !
        !>@param y_map
        !> coordinate map along the y-direction
        !
        !>@param nodes_tmp
        !> governing variables at t-dt
        !
        !>@param nodes
        !> governing variables at t
        !-------------------------------------------------------------
        subroutine apply_bc_on_nodes(t,x_map,y_map,nodes_tmp,nodes)
        
          implicit none
        
          real(rkind)                     , intent(in)    :: t
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

          integer, dimension(8) :: bc_nodes_order
          integer               :: k


          ! choose order in which the boundary conditions are applied
          call choose_bc_nodes_order(bc_nodes_order)


          ! apply the boundary conditions for each layer
          do k=1, size(bc_nodes_order,1)

             select case(bc_nodes_order(k))
             
               ! N layer
               case(N_edge_type)
                  call apply_bc_on_nodes_with_bc_operator(
     $                 N_edge_section,
     $                 bc_N_choice,
     $                 t,x_map,y_map,
     $                 nodes_tmp,nodes)

               ! S layer
               case(S_edge_type)
                  call apply_bc_on_nodes_with_bc_operator(
     $                 S_edge_section,
     $                 bc_S_choice,
     $                 t,x_map,y_map,
     $                 nodes_tmp,nodes)

               ! E layer
               case(E_edge_type)
                  call apply_bc_on_nodes_with_bc_operator(
     $                 E_edge_section,
     $                 bc_E_choice,
     $                 t,x_map,y_map,
     $                 nodes_tmp,nodes)

               ! W layer
               case(W_edge_type)
                  call apply_bc_on_nodes_with_bc_operator(
     $                 W_edge_section,
     $                 bc_W_choice,
     $                 t,x_map,y_map,
     $                 nodes_tmp,nodes)

               ! NW corner
               case(NW_corner_type)
                  call apply_bc_on_nodes_with_bc_operator(
     $                 NW_corner_section,
     $                 bc_N_choice,
     $                 t,x_map,y_map,
     $                 nodes_tmp,nodes)

               ! NE corner
               case(NE_corner_type)
                  call apply_bc_on_nodes_with_bc_operator(
     $                 NE_corner_section,
     $                 bc_N_choice,
     $                 t,x_map,y_map,
     $                 nodes_tmp,nodes)

               ! SW corner
               case(SW_corner_type)
                  call apply_bc_on_nodes_with_bc_operator(
     $                 SW_corner_section,
     $                 bc_S_choice,
     $                 t,x_map,y_map,
     $                 nodes_tmp,nodes)

               ! SE corner
               case(SE_corner_type)
                  call apply_bc_on_nodes_with_bc_operator(
     $                 SE_corner_section,
     $                 bc_S_choice,
     $                 t,x_map,y_map,
     $                 nodes_tmp,nodes)

             end select

          end do

        end subroutine apply_bc_on_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the nodes using the bc_operators
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> coordinate map along the x-direction
        !
        !>@param y_map
        !> coordinate map along the y-direction
        !
        !>@param nodes_tmp
        !> governing variables at t-dt
        !
        !>@param nodes
        !> governing variables at t
        !-------------------------------------------------------------
        subroutine apply_bc_on_nodes_with_bc_operator(
     $     bc_section, bc_choice,
     $     t,x_map,y_map,nodes_tmp,nodes)
        
          implicit none
        
          integer    , dimension(4)       , intent(in)    :: bc_section
          integer                         , intent(in)    :: bc_choice
          real(rkind)                     , intent(in)    :: t
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes


          type(bc_operators_reflection_xy) :: bc_reflection_xy
          type(bc_operators_wall_xy)       :: bc_wall_xy
          

          select case(bc_choice)

            case(reflection_xy_choice)

               call bc_reflection_xy%apply_bc_on_nodes(
     $              bc_section,
     $              t,x_map,y_map,
     $              nodes_tmp,nodes)

            case(wall_xy_choice)

               call bc_wall_xy%apply_bc_on_nodes(
     $              bc_section,
     $              t,x_map,y_map,
     $              nodes_tmp,nodes)
          
            case default
               call error_bc_choice(
     $              'bc_operators_gen_class',
     $              'apply_bc_on_nodes_with_bc_operator',
     $              bc_choice)
               
          end select

        end subroutine apply_bc_on_nodes_with_bc_operator


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes along the x directions at the
        !> edge of the computational domain
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates
        !
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operators
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !-------------------------------------------------------------
        subroutine apply_bc_on_fluxes(
     $     t,x_map,y_map,nodes,s,flux_x,flux_y)
        
          implicit none
        
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(sd_operators)                , intent(in)    :: s
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
        

          ! SW corner
          call apply_bc_on_fluxes_with_bc_operator(
     $         SW_corner_section,
     $         bc_S_choice,
     $         t,x_map,y_map,nodes,
     $         s,flux_x,flux_y)

          ! S layer
          call apply_bc_on_fluxes_with_bc_operator(
     $         S_edge_section,
     $         bc_S_choice,
     $         t,x_map,y_map,nodes,
     $         s,flux_x,flux_y)

          ! SE corner
          call apply_bc_on_fluxes_with_bc_operator(
     $         SE_corner_section,
     $         bc_S_choice,
     $         t,x_map,y_map,nodes,
     $         s,flux_x,flux_y)

          ! W layer
          call apply_bc_on_fluxes_with_bc_operator(
     $         W_edge_section,
     $         bc_W_choice,
     $         t,x_map,y_map,nodes,
     $         s,flux_x,flux_y)

          ! E layer
          call apply_bc_on_fluxes_with_bc_operator(
     $         E_edge_section,
     $         bc_E_choice,
     $         t,x_map,y_map,nodes,
     $         s,flux_x,flux_y)

          ! NW corner
          call apply_bc_on_fluxes_with_bc_operator(
     $         NW_corner_section,
     $         bc_N_choice,
     $         t,x_map,y_map,nodes,
     $         s,flux_x,flux_y)

          ! N layer
          call apply_bc_on_fluxes_with_bc_operator(
     $         N_edge_section,
     $         bc_N_choice,
     $         t,x_map,y_map,nodes,
     $         s,flux_x,flux_y)
          
          ! NE corner
          call apply_bc_on_fluxes_with_bc_operator(
     $         NE_corner_section,
     $         bc_N_choice,
     $         t,x_map,y_map,nodes,
     $         s,flux_x,flux_y)

        end subroutine apply_bc_on_fluxes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes using the bc_operators
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> coordinate map along the x-direction
        !
        !>@param y_map
        !> coordinate map along the y-direction
        !
        !>@param nodes_tmp
        !> governing variables at t-dt
        !
        !>@param nodes
        !> governing variables at t
        !-------------------------------------------------------------
        subroutine apply_bc_on_fluxes_with_bc_operator(
     $     bc_section, bc_choice,
     $     t,x_map,y_map,nodes,
     $     s,flux_x,flux_y)
        
          implicit none
        
          integer    , dimension(4)         , intent(in)    :: bc_section
          integer                           , intent(in)    :: bc_choice
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(sd_operators)                , intent(in)    :: s
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y


          type(bc_operators_wall_xy) :: bc_wall_xy

          
          select case(bc_choice)

            case(wall_xy_choice)

               ! there is here that the fluxes should be splitted
               ! into viscid and invisicd fluxes for the Yoo-Lodato
               ! b.c. using optional arguments in apply_bc_on_fluxes
               ! for the wall b.c.
               call bc_wall_xy%apply_bc_on_fluxes(
     $              bc_section,
     $              t,x_map,y_map,
     $              nodes,s,flux_x,flux_y)
          
            case default
               call error_bc_choice(
     $              'bc_operators_gen_class',
     $              'apply_bc_on_fluxes_with_bc_operator',
     $              bc_choice)
               
          end select

        end subroutine apply_bc_on_fluxes_with_bc_operator


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes along the x directions at the
        !> edge of the computational domain
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates
        !
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operators
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev(
     $     t,x_map,y_map,nodes,
     $     p_model,
     $     flux_x,flux_y,timedev)
        
          implicit none
        
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev
        

          ! SW corner
          call apply_bc_on_timedev_with_bc_operator(
     $         [SW_corner_type,1,1,0],
     $         bc_S_choice,
     $         t,x_map,y_map,nodes,
     $         p_model,
     $         flux_x, flux_y, timedev)

          ! S layer
          call apply_bc_on_timedev_with_bc_operator(
     $         [S_edge_type,bc_size+1,1,nx-bc_size],
     $         bc_S_choice,
     $         t,x_map,y_map,nodes,
     $         p_model,
     $         flux_x, flux_y, timedev)

          ! SE corner
          call apply_bc_on_timedev_with_bc_operator(
     $         [SE_corner_type,nx-bc_size+1,1,0],
     $         bc_S_choice,
     $         t,x_map,y_map,nodes,
     $         p_model,
     $         flux_x, flux_y, timedev)

          ! W layer
          call apply_bc_on_timedev_with_bc_operator(
     $         [W_edge_type,1,bc_size+1,ny-bc_size+1],
     $         bc_W_choice,
     $         t,x_map,y_map,nodes,
     $         p_model,
     $         flux_x, flux_y, timedev)

          ! E layer
          call apply_bc_on_timedev_with_bc_operator(
     $         [E_edge_type,nx-bc_size+1,bc_size+1,ny-bc_size+1],
     $         bc_E_choice,
     $         t,x_map,y_map,nodes,
     $         p_model,
     $         flux_x, flux_y, timedev)

          ! NW corner
          call apply_bc_on_timedev_with_bc_operator(
     $         [NW_corner_type,1,ny-bc_size+1,0],
     $         bc_N_choice,
     $         t,x_map,y_map,nodes,
     $         p_model,
     $         flux_x, flux_y, timedev)

          ! N layer
          call apply_bc_on_timedev_with_bc_operator(
     $         [N_edge_type,bc_size+1,ny-bc_size+1,nx-bc_size],
     $         bc_N_choice,
     $         t,x_map,y_map,nodes,
     $         p_model,
     $         flux_x, flux_y, timedev)
          
          ! NE corner
          call apply_bc_on_timedev_with_bc_operator(
     $         [NE_corner_type,nx-bc_size+1,ny-bc_size+1,0],
     $         bc_N_choice,
     $         t,x_map,y_map,nodes,
     $         p_model,
     $         flux_x, flux_y, timedev)

        end subroutine apply_bc_on_timedev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes using the bc_operators
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> coordinate map along the x-direction
        !
        !>@param y_map
        !> coordinate map along the y-direction
        !
        !>@param nodes_tmp
        !> governing variables at t-dt
        !
        !>@param nodes
        !> governing variables at t
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev_with_bc_operator(
     $     bc_section, bc_choice,
     $     t,x_map,y_map,nodes,
     $     p_model,
     $     flux_x, flux_y, timedev)
        
          implicit none
        
          integer    , dimension(4)         , intent(in)    :: bc_section
          integer                           , intent(in)    :: bc_choice
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev


          type(bc_operators_hedstrom_xy) :: bc_hedstrom_xy

          
          select case(bc_choice)

            case(hedstrom_xy_choice)

               ! there is here that the fluxes should be splitted
               ! into viscid and invisicd fluxes for the Yoo-Lodato
               ! b.c. using optional arguments in apply_bc_on_fluxes
               ! for the wall b.c.
               call bc_hedstrom_xy%apply_bc_on_timedev(
     $              bc_section,
     $              t,x_map,y_map,nodes,
     $              p_model,
     $              flux_x,flux_y,timedev)
          
            case default
               call error_bc_choice(
     $              'bc_operators_gen_class',
     $              'apply_bc_on_timedev_with_bc_operator',
     $              bc_choice)
               
          end select

        end subroutine apply_bc_on_timedev_with_bc_operator


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> error printed if bc_choice nto recognized
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param mainlayer_id
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine choose_bc_nodes_order(bc_nodes_order)

          implicit none

          integer, dimension(8), intent(out) :: bc_nodes_order

          ! the application of the wall b.c. in N and S layers 
          ! requires to have the nodes in the E and W layers
          ! computed first (for the x-gradient of the mass
          ! density in the wall b.c....)
          if((bc_N_choice.eq.wall_xy_choice).or.
     $       (bc_S_choice.eq.wall_xy_choice)) then

             bc_nodes_order = [
     $            W_edge_type,
     $            E_edge_type,
     $            SW_corner_type,
     $            S_edge_type,
     $            SE_corner_type,
     $            NW_corner_type,
     $            N_edge_type,
     $            NE_corner_type]

          ! the application of the wall b.c./reflection/periodic
          ! in the E and W layers requires to have the nodes
          ! in the N and S layers
          else

                bc_nodes_order = [
     $               S_edge_type,
     $               N_edge_type,
     $               W_edge_type,
     $               E_edge_type,
     $               SW_corner_type,
     $               SE_corner_type,
     $               NW_corner_type,
     $               NE_corner_type]             

          end if

        end subroutine choose_bc_nodes_order


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> error printed if bc_choice nto recognized
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param mainlayer_id
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine error_mainlayer_id(
     $     file_name, fct_name, mainlayer_id)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: mainlayer_id

          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'bc_choice not recognized'
          write(ERROR_UNIT, '(''bc_choice: '',I2)') mainlayer_id

          stop 'error_bc_choice'

        end subroutine error_mainlayer_id

      end module bc_operators_gen_class

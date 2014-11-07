       !> @file
       !> module implemeting subroutines to apply open boundary
       !> conditions at the edge of the computational domain using
       !> yoo and lodato like boundary conditions
       !
       !> @author 
       !> Julien L. Desmarais
       !
       !> @brief
       !> module implemeting subroutines to apply open boundary
       !> conditions at the edge of the computational domain using
       !> Yoo and Lodato boundary conditions
       !
       !> @date
       !> 09_09_2014 - initial version - J.L. Desmarais
       !-----------------------------------------------------------------
       module lodi_timedev_xy_module

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        use openbc_operators_module, only :
     $       add_body_forces

        use lodi_edge_inflow_class, only :
     $       lodi_edge_inflow

        use lodi_edge_outflow_class, only :
     $       lodi_edge_outflow

        use lodi_corner_inflow_inflow_class, only :
     $       lodi_corner_inflow_inflow

        use lodi_corner_inflow_outflow_class, only :
     $       lodi_corner_inflow_outflow

        use lodi_corner_outflow_outflow_class, only :
     $       lodi_corner_outflow_outflow

        use parameters_constant, only :
     $       left,
     $       inflow_type,
     $       outflow_type,
     $       ask_flow, always_inflow, always_outflow

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq        
        
        implicit none


        private
        public :: 
     $       compute_timedev_x_edge,
     $       compute_timedev_y_edge,
     $       compute_timedev_corner,
     $       get_flow_config

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives at a boundary
        !> normal to the x-direction using Yoo boundary conditions
        !
        !> @date
        !> 09_09_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> governing equations of the physical model
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array of grid points
        !
        !>@param x_map
        !> coordinate map in the x-direction
        !
        !>@param y_map
        !> coordinate map in the y-direction
        !
        !>@param i
        !> index identifying the grid point position along the x-axis
        !
        !>@param j
        !> index identifying the grid point position along the y-axis
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param transverse_lodi
        !> vector of the transverse lodi terms
        !
        !>@param viscous_lodi
        !> vector of the viscous lodi terms
        !
        !>@param side_x
        !> boolean determining the spatial location of the b.c.
        !
        !>@param gradient_x
        !> procedure computing the gradient along the x-direction
        !
        !>@param inflow_bc
        !> procedure computing the contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives using an inflow
        !> boundary condition
        !
        !>@param outflow_bc
        !> procedure computing the contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives using an outflow
        !> boundary condition
        !
        !>@param flow_x_user
        !> type of boundary determined by the user (force teh boundary
        !> to be of type inflow or outflow or let the flow decides)
        !
        !>@return timedev
        !> time derivatives of the governing variables
        !-------------------------------------------------------------
        subroutine compute_timedev_x_edge(
     $       p_model,
     $       t,nodes,x_map,y_map,i,j,
     $       flux_y,
     $       transverse_lodi, viscous_lodi,
     $       side_x,
     $       gradient_x,
     $       inflow_bc,
     $       outflow_bc,
     $       flow_x_user,
     $       timedev)

          implicit none

          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          integer(ikind)                    , intent(in)    :: i
          integer(ikind)                    , intent(in)    :: j
          real(rkind), dimension(nx,ny+1,ne), intent(in)    :: flux_y
          real(rkind), dimension(ne)        , intent(in)    :: transverse_lodi
          real(rkind), dimension(ne)        , intent(in)    :: viscous_lodi
          logical                           , intent(in)    :: side_x
          procedure(gradient_x_proc)                        :: gradient_x
          type(lodi_edge_inflow)            , intent(in)    :: inflow_bc
          type(lodi_edge_outflow)           , intent(in)    :: outflow_bc
          integer                           , intent(in)    :: flow_x_user
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          
          logical     :: flow_x_type
          real(rkind) :: dy


          dy = y_map(2) - y_map(1)


          !determine whether the flow at the
          !x-boundary is inflow or outflow
          flow_x_type = get_flow_config(nodes(i,j,2),side_x,flow_x_user)


          !if the b.c. is of inlet type, the lodi_edge_inflow
          !is applied
          if(flow_x_type.eqv.inflow_type) then

             timedev(i,j,:) =
     $            inflow_bc%compute_x_timedev(
     $               p_model,
     $               t,nodes,x_map,y_map,i,j,
     $               transverse_lodi, viscous_lodi,
     $               side_x,
     $               gradient_x) +
     $            
     $            1.0d0/dy*(flux_y(i,j,:) - flux_y(i,j+1,:)) +
     $         
     $            add_body_forces(
     $            p_model,
     $            t,x_map(i),y_map(j),nodes(i,j,:))


          !otherwise, if the b.c. is of outlet type,
          !the lodi_edge_outflow is applied
          else

             timedev(i,j,:) =
     $            outflow_bc%compute_x_timedev(
     $               p_model,
     $               t,nodes,x_map,y_map,i,j,
     $               transverse_lodi, viscous_lodi,
     $               side_x,
     $               gradient_x) +
     $            
     $            1.0d0/dy*(flux_y(i,j,:) - flux_y(i,j+1,:)) +
     $         
     $            add_body_forces(
     $            p_model,
     $            t,x_map(i),y_map(j),nodes(i,j,:))

          end if

        end subroutine compute_timedev_x_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives at a boundary
        !> normal to the y-direction using Yoo boundary conditions
        !
        !> @date
        !> 09_09_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> governing equations of the physical model
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array of grid points
        !
        !>@param x_map
        !> coordinate map in the x-direction
        !
        !>@param y_map
        !> coordinate map in the y-direction
        !
        !>@param i
        !> index identifying the grid point position along the x-axis
        !
        !>@param j
        !> index identifying the grid point position along the y-axis
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param transverse_lodi
        !> vector of the transverse lodi terms
        !
        !>@param viscous_lodi
        !> vector of the viscous lodi terms
        !
        !>@param side_y
        !> boolean determining the spatial location of the b.c.
        !
        !>@param gradient_y
        !> procedure computing the gradient along the y-direction
        !
        !>@param inflow_bc
        !> procedure computing the contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives using an inflow
        !> boundary condition
        !
        !>@param outflow_bc
        !> procedure computing the contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives using an outflow
        !> boundary condition
        !
        !>@param flow_y_user
        !> type of boundary determined by the user (force the boundary
        !> to be of type inflow or outflow or let the flow decides)
        !
        !>@return timedev
        !> time derivatives of the governing variables
        !-------------------------------------------------------------
        subroutine compute_timedev_y_edge(
     $     p_model,
     $     t,nodes,x_map,y_map,i,j,
     $     flux_x,
     $     transverse_lodi, viscous_lodi,
     $     side_y,
     $     gradient_y,
     $     inflow_bc,
     $     outflow_bc,
     $     flow_y_user,
     $     timedev)

          implicit none

          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          integer(ikind)                    , intent(in)    :: i
          integer(ikind)                    , intent(in)    :: j
          real(rkind), dimension(nx+1,ny,ne), intent(in)    :: flux_x
          real(rkind), dimension(ne)        , intent(in)    :: transverse_lodi
          real(rkind), dimension(ne)        , intent(in)    :: viscous_lodi
          logical                           , intent(in)    :: side_y
          procedure(gradient_y_proc)                        :: gradient_y
          type(lodi_edge_inflow)            , intent(in)    :: inflow_bc
          type(lodi_edge_outflow)           , intent(in)    :: outflow_bc
          integer                           , intent(in)    :: flow_y_user
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          
          logical     :: flow_y_type
          real(rkind) :: dx


          dx = x_map(2) - x_map(1)


          !determine whether the flow at the
          !x-boundary is inflow or outflow
          flow_y_type = get_flow_config(nodes(i,j,3),side_y,flow_y_user)


          !if the b.c. is of inlet type, the lodi_edge_inflow
          !is applied
          if(flow_y_type.eqv.inflow_type) then

             timedev(i,j,:) =
     $            inflow_bc%compute_y_timedev(
     $               p_model,
     $               t,nodes,x_map,y_map,i,j,
     $               transverse_lodi, viscous_lodi,
     $               side_y,
     $               gradient_y) +
     $            
     $            1.0d0/dx*(flux_x(i,j,:) - flux_x(i+1,j,:)) +
     $         
     $            add_body_forces(
     $            p_model,
     $            t,x_map(i),y_map(j),nodes(i,j,:))


          !otherwise, if the b.c. is of outlet type,
          !the lodi_edge_outflow is applied
          else

             timedev(i,j,:) =
     $            outflow_bc%compute_y_timedev(
     $               p_model,
     $               t,nodes,x_map,y_map,i,j,
     $               transverse_lodi, viscous_lodi,
     $               side_y,
     $               gradient_y) +
     $            
     $            1.0d0/dx*(flux_x(i,j,:) - flux_x(i+1,j,:)) +
     $         
     $            add_body_forces(
     $            p_model,
     $            t,x_map(i),y_map(j),nodes(i,j,:))

          end if

        end subroutine compute_timedev_y_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives at the corner
        !> using Lodato et al. b.c.
        !
        !> @date
        !> 09_09_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> governing equations of the physical model
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array of grid points
        !
        !>@param x_map
        !> coordinate map in the x-direction
        !
        !>@param y_map
        !> coordinate map in the y-direction
        !
        !>@param i
        !> index identifying the grid point position along the x-axis
        !
        !>@param j
        !> index identifying the grid point position along the y-axis
        !
        !>@param side_x
        !> boolean determining the spatial location of the b.c.
        !
        !>@param side_y
        !> boolean determining the spatial location of the b.c.
        !
        !>@param gradient_x
        !> procedure computing the gradient along the x-direction
        !
        !>@param gradient_y
        !> procedure computing the gradient along the y-direction
        !
        !>@param inflow_bc
        !> procedure computing the contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives using an inflow
        !> boundary condition
        !
        !>@param outflow_bc
        !> procedure computing the contribution of the hyperbolic terms
        !> along the x-direction to the time derivatives using an outflow
        !> boundary condition
        !
        !>@param flow_x_user
        !> type of boundary determined by the user (force the boundary
        !> to be of type inflow or outflow or let the flow decides)
        !
        !>@param flow_y_user
        !> type of boundary determined by the user (force the boundary
        !> to be of type inflow or outflow or let the flow decides)
        !
        !>@return timedev
        !> time derivatives of the governing variables
        !-------------------------------------------------------------
        subroutine compute_timedev_corner(
     $     p_model,
     $     t,nodes,x_map,y_map,i,j,
     $     side_x, side_y,
     $     gradient_x, gradient_y,
     $     corner_inflow_inflow,
     $     corner_inflow_outflow,
     $     corner_outflow_outflow,
     $     flow_x_user, flow_y_user,
     $     timedev)

          implicit none

          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          integer(ikind)                    , intent(in)    :: i
          integer(ikind)                    , intent(in)    :: j
          logical                           , intent(in)    :: side_x
          logical                           , intent(in)    :: side_y
          procedure(gradient_y_proc)                        :: gradient_x
          procedure(gradient_y_proc)                        :: gradient_y
          type(lodi_corner_inflow_inflow)   , intent(in)    :: corner_inflow_inflow
          type(lodi_corner_inflow_outflow)  , intent(in)    :: corner_inflow_outflow
          type(lodi_corner_outflow_outflow) , intent(in)    :: corner_outflow_outflow
          integer                           , intent(in)    :: flow_x_user
          integer                           , intent(in)    :: flow_y_user
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev


          logical :: flow_x_type
          logical :: flow_y_type


          !determine whether the flow at the
          !x-boundary is inflow or outflow
          flow_x_type = get_flow_config(nodes(i,j,2),side_x,flow_x_user)
          flow_y_type = get_flow_config(nodes(i,j,3),side_y,flow_y_user)


          !choice of b.c.:
          ! x:inlet/y:inlet  -> corner_inflow_inflow
          ! x:inlet/y:outlet -> corner_inflow_outflow
          ! x:outlet/y:inlet -> corner_inflow_outflow
          ! x:outlet/y:outlet -> corner_outflow_outflow
          if(flow_x_type.eqv.inflow_type) then
             
             if(flow_y_type.eqv.inflow_type) then

                timedev(i,j,:) =
     $               corner_inflow_inflow%compute_x_and_y_timedev(
     $               p_model,
     $               t, nodes, x_map, y_map, i,j,
     $               side_x, side_y,
     $               flow_x_type, flow_y_type,
     $               gradient_x, gradient_y)
                
             else

                timedev(i,j,:) =
     $               corner_inflow_outflow%compute_x_and_y_timedev(
     $               p_model,
     $               t, nodes, x_map, y_map, i,j,
     $               side_x, side_y,
     $               flow_x_type, flow_y_type,
     $               gradient_x, gradient_y)

             end if

          else

             if(flow_y_type.eqv.inflow_type) then

                timedev(i,j,:) =
     $               corner_inflow_outflow%compute_x_and_y_timedev(
     $               p_model,
     $               t, nodes, x_map, y_map, i,j,
     $               side_x, side_y,
     $               flow_x_type, flow_y_type,
     $               gradient_x, gradient_y)
                
             else

                timedev(i,j,:) =
     $               corner_outflow_outflow%compute_x_and_y_timedev(
     $               p_model,
     $               t, nodes, x_map, y_map, i,j,
     $               side_x, side_y,
     $               flow_x_type, flow_y_type,
     $               gradient_x, gradient_y)

             end if

          end if

        end subroutine compute_timedev_corner


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives at the corner
        !> using Lodato et al. b.c.
        !
        !> @date
        !> 09_09_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> governing equations of the physical model
        !
        !>@param vector
        !> vector with the same direction as the velocity component
        !> normal to the boundary
        !
        !>@param side
        !> boolean identifying the spatial location of the boundary
        !
        !>@param usr_config
        !> configuration of the flow from the user (either forced
        !> to be inflow or outflow or ask the flow to decide)
        !-------------------------------------------------------------
        function get_flow_config(vector,side,usr_config)
     $     result(flow_config)

          implicit none

          real(rkind), intent(in) :: vector
          logical    , intent(in) :: side
          integer    , intent(in) :: usr_config
          logical                 :: flow_config

          
          if(usr_config.eq.ask_flow) then
             if(side.eqv.left) then
                if(vector.gt.0) then
                   flow_config = inflow_type
                else
                   flow_config = outflow_type
                end if

             else
                if(vector.gt.0) then
                   flow_config = outflow_type
                else
                   flow_config = inflow_type
                end if

             end if

          else
             select case(usr_config)
               case(always_inflow)
                  flow_config = inflow_type
               case(always_outflow)
                  flow_config = outflow_type
               case default
                  print '(''lodi_timedev_xy_module'')'
                  print '(''get_flow_config'')'
                  print '(''usr_config not recognized'')'
                  print '(''usr_config: '',I2)', usr_config
                  stop
             end select 
             
          end if

        end function get_flow_config

      end module lodi_timedev_xy_module

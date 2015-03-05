      !> @file
      !> class encapsulating subroutines to compute
      !> the time derivatives of the main variables
      !> using the finite volume method
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the time derivatives of the main variables
      !> using the finite volume method
      !
      !> @date
      !> 13_08_2013 - initial version               - J.L. Desmarais
      !> 22_10_2014 - time dev for buffer layers    - J.L. Desmarais
      !-----------------------------------------------------------------
      module td_operators_class

        use bc_operators_class, only :
     $       bc_operators

        use sd_operators_class, only :
     $       sd_operators

        use parameters_constant, only :
     $       bc_fluxes_choice,
     $       bc_timedev_choice

        use parameters_bf_layer, only :
     $       interior_pt,
     $       bc_interior_pt,
     $       BF_SUCCESS

        use parameters_input, only : 
     $       nx,ny,ne,bc_size,
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice

        use parameters_kind, only :
     $       rkind, ikind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use td_operators_abstract_class, only :
     $       td_operators_abstract

        implicit none

        private
        public :: td_operators


        !> @class td_operators
        !> class encapsulating operators to compute
        !> the time derivatives of the main variables
        !> using the finite volume method
        !>
        !> @param compute_time_dev
        !> compute the time derivatives
        !
        !> @param compute_time_dev_nopt
        !> compute the time derivatives without optimizing the
        !> size of the arrays passed as arguments
        !---------------------------------------------------------------
        type, extends(td_operators_abstract) :: td_operators

          contains

          procedure, nopass :: compute_time_dev
          procedure, nopass :: compute_time_dev_nopt

        end type td_operators


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives using the
        !> space discretisation operators and the physical model
        !
        !> @date
        !> 13_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param s
        !> space discretization operators
        !
        !>@param p_model
        !> physical model
        !
        !>@param bc_used
        !> boundary conditions
        !
        !>@param time_dev
        !> time derivatives
        !--------------------------------------------------------------
        function compute_time_dev(
     $       t,nodes,x_map,y_map,
     $       s,p_model,bc_used,
     $       bc_sections)
     $       result(time_dev)

          implicit none

          real(rkind)                                          , intent(in) :: t
          real(rkind), dimension(nx,ny,ne)                     , intent(in) :: nodes
          real(rkind), dimension(nx)                           , intent(in) :: x_map
          real(rkind), dimension(ny)                           , intent(in) :: y_map
          type(sd_operators)                                   , intent(in) :: s
          type(pmodel_eq)                                      , intent(in) :: p_model
          type(bc_operators)                                   , intent(in) :: bc_used
          integer(ikind), dimension(:,:), allocatable, optional, intent(in) :: bc_sections
          real(rkind), dimension(nx,ny,ne)                                  :: time_dev

          real(rkind)                        :: dx
          real(rkind)                        :: dy
          integer                            :: k
          integer(ikind)                     :: i,j
          real(rkind), dimension(nx+1,ny,ne) :: flux_x
          real(rkind), dimension(nx,ny+1,ne) :: flux_y

          integer(ikind), dimension(2,2) :: bf_alignment


          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)
          

          !<compute the fluxes
          !FORCEINLINE RECURSIVE
          flux_x = p_model%compute_flux_x(nodes,dx,dy,s)

          !FORCEINLINE RECURSIVE
          flux_y = p_model%compute_flux_y(nodes,dx,dy,s)


          !< if the boundary conditions influence the computation
          !> of the fluxes, then we need to modify the fluxes
          if((bc_N_type_choice.eq.bc_fluxes_choice).or.
     $       (bc_S_type_choice.eq.bc_fluxes_choice).or.
     $       (bc_E_type_choice.eq.bc_fluxes_choice).or.
     $       (bc_W_type_choice.eq.bc_fluxes_choice)
     $    ) then

             call bc_used%apply_bc_on_fluxes(
     $            nodes,dx,dy,s,flux_x,flux_y)

          end if


          !<compute the time derivatives
          do k=1, ne
             do j=1+bc_size, ny-bc_size
                do i=1+bc_size, nx-bc_size

                   time_dev(i,j,k)=
     $                  (flux_x(i,j,k)-flux_x(i+1,j,k))/dx+
     $                  (flux_y(i,j,k)-flux_y(i,j+1,k))/dy
                   
                   time_dev(i,j,k)=time_dev(i,j,k)+
     $                  p_model%compute_body_forces(
     $                  t,x_map(i),y_map(j),
     $                  nodes(i,j,:),k)

                end do
             end do
          end do

          !< if the boundary conditions influence the computation
          !> of the time derivatives, then we need to compute the
          !> time derivatives at the boundary
          if(
     $         (bc_N_type_choice.eq.bc_timedev_choice).or.
     $         (bc_S_type_choice.eq.bc_timedev_choice).or.
     $         (bc_E_type_choice.eq.bc_timedev_choice).or.
     $         (bc_W_type_choice.eq.bc_timedev_choice)
     $    ) then

             !if not all the time derivatives of the boundary
             !layers should be computed, the selected time
             !derivatives to be computed are given by
             !bc_sections
             if(present(bc_sections)) then

                bf_alignment(1,1) = bc_size+1
                bf_alignment(1,2) = nx-bc_size
                bf_alignment(2,1) = bc_size+1
                bf_alignment(2,2) = ny-bc_size

                call bc_used%apply_bc_on_timedev_nopt(
     $               p_model,t,
     $               nodes,
     $               bf_alignment,
     $               nodes,x_map,y_map,
     $               flux_x,flux_y,
     $               time_dev,
     $               bc_sections)

             !if all the time derivatives of the boundary
             !layers are computed
             else
                call bc_used%apply_bc_on_timedev(
     $               p_model,
     $               t,nodes,x_map,y_map,
     $               flux_x,flux_y,
     $               time_dev)
             end if          
          end if

        end function compute_time_dev



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives using the
        !> space discretisation operators and the physical model
        !> for non-fixed size arrays
        !
        !> @date
        !> 22_10_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param s
        !> space discretization operators
        !
        !>@param p_model
        !> physical model
        !
        !>@param bc_used
        !> boundary conditions
        !
        !>@param time_dev
        !> time derivatives
        !
        !>@param grdpts_id
        !> array containing the role of the grid points (interior_pt,
        !> bc_interior_pt, bc_pt, no_pt)
        !
        !>@param bc_sections
        !> array identifying the boundary layers
        !
        !>@param x_borders
        !> array containing the limits of the computed grid points in
        !> the x-direction
        !
        !>@param y_borders
        !> array containing the limits of the computed grid points in
        !> the y-direction
        !
        !>@param N_bc_sections
        !> determine whether the last two north lines of the array
        !> should be computed by the buffer layer itself or not
        !
        !>@param S_bc_sections
        !> determine whether the last two south lines of the array
        !> should be computed by the buffer layer itself or not
        !--------------------------------------------------------------
        subroutine compute_time_dev_nopt(
     $     t,nodes,x_map,y_map,
     $     s,p_model,bc_used,
     $     time_dev,
     $     bf_alignment,
     $     grdpts_id,
     $     interior_nodes,
     $     bc_sections,
     $     x_borders,
     $     y_borders)

            implicit none

            real(rkind)                                  , intent(in)  :: t
            integer(ikind), dimension(2,2)               , intent(in)  :: bf_alignment
            integer       , dimension(:,:)               , intent(in)  :: grdpts_id
            real(rkind)   , dimension(:)                 , intent(in)  :: x_map
            real(rkind)   , dimension(:)                 , intent(in)  :: y_map
            real(rkind)   , dimension(:,:,:)             , intent(in)  :: nodes
            type(sd_operators)                           , intent(in)  :: s
            type(pmodel_eq)                              , intent(in)  :: p_model
            type(bc_operators)                           , intent(in)  :: bc_used
            real(rkind)   , dimension(:,:,:)             , intent(out) :: time_dev
            real(rkind)   , dimension(nx,ny,ne)          , intent(in)  :: interior_nodes
            integer       , dimension(:,:)  , allocatable, intent(in)  :: bc_sections
            integer(ikind), dimension(2)                 , intent(in)  :: x_borders
            integer(ikind), dimension(2)                 , intent(in)  :: y_borders

            real(rkind), dimension(:,:,:), allocatable :: flux_x
            real(rkind), dimension(:,:,:), allocatable :: flux_y
            real(rkind)                                :: dx
            real(rkind)                                :: dy
            integer(ikind)                             :: i,j
            integer                                    :: k

            real(rkind) :: t_s
            logical     :: ierror

            t_s = t
            ierror = BF_SUCCESS          


            !allocate space for the flux computation
            allocate(flux_x(size(nodes,1)+1,size(nodes,2),ne))
            allocate(flux_y(size(nodes,1),size(nodes,2)+1,ne))


            !space steps
            dx = x_map(2)-x_map(1)
            dy = y_map(2)-y_map(1)


            !compute the fluxes
            !FORCEINLINE RECURSIVE
            call p_model%compute_flux_x_nopt(
     $           nodes,dx,dy,s,
     $           grdpts_id,
     $           flux_x,
     $           [bc_size+1,size(nodes,1)-bc_size],
     $           [bc_size+1,size(nodes,2)-bc_size])


            !FORCEINLINE RECURSIVE
            call p_model%compute_flux_y_nopt(
     $           nodes,dx,dy,s,
     $           grdpts_id,
     $           flux_y,
     $           [bc_size+1,size(nodes,1)-bc_size],
     $           [bc_size+1,size(nodes,2)-bc_size])


            !if the boundary conditions influence the computation
            !of the fluxes, then we need to modify the fluxes
            if(
     $           (bc_N_type_choice.eq.bc_fluxes_choice).or.
     $           (bc_S_type_choice.eq.bc_fluxes_choice).or.
     $           (bc_E_type_choice.eq.bc_fluxes_choice).or.
     $           (bc_W_type_choice.eq.bc_fluxes_choice)
     $      ) then

               call bc_used%apply_bc_on_fluxes(
     $              nodes,dx,dy,s,flux_x,flux_y)

            end if


            !compute the time derivatives
            do k=1, ne
               do j=y_borders(1), y_borders(2)
                  do i=x_borders(1), x_borders(2)
                        
                     if(grdpts_id(i,j).eq.interior_pt) then
                              
                        time_dev(i,j,k)=
     $                       (flux_x(i,j,k)/dx-flux_x(i+1,j,k)/dx)+
     $                       (flux_y(i,j,k)/dy-flux_y(i,j+1,k)/dy)
                        
                        time_dev(i,j,k)=
     $                       time_dev(i,j,k)+
     $                       p_model%compute_body_forces(
     $                       t, x_map(i), y_map(j),
     $                       nodes(i,j,:),k)
                           
                     end if
                        
                  end do
               end do
            end do


            !if the boundary conditions influence the computation
            !of the time derivatives, then we need to compute the
            !time derivatives at the boundary
            if((bc_N_type_choice.eq.bc_timedev_choice).or.
     $         (bc_S_type_choice.eq.bc_timedev_choice).or.
     $         (bc_E_type_choice.eq.bc_timedev_choice).or.
     $         (bc_W_type_choice.eq.bc_timedev_choice)) then

               call bc_used%apply_bc_on_timedev_nopt(
     $              p_model,t,
     $              interior_nodes,
     $              bf_alignment,
     $              nodes,x_map,y_map,
     $              flux_x,flux_y,
     $              time_dev,
     $              bc_sections,
     $              grdpts_id=grdpts_id)

            end if

            deallocate(flux_x)
            deallocate(flux_y)

        end subroutine compute_time_dev_nopt

      end module td_operators_class

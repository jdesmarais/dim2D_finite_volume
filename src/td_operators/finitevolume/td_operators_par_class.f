      !> @file
      !> class encapsulating subroutines to compute
      !> the time derivatives of the main variables
      !> using the finite volume method in a
      !> distributed memory system
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the time derivatives of the main variables
      !> using the finite volume method in a
      !> distributed memory system
      !
      !> @date
      !> 25_09_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module td_operators_par_class

        use bc_operators_par_class, only :
     $     bc_operators_par

        use parameters_constant, only :
     $       bc_fluxes_choice

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

        use sd_operators_class, only :
     $       sd_operators

        use td_operators_abstract_par_class , only :
     $       td_operators_abstract_par

        implicit none

        private
        public :: td_operators_par


        !> @class td_operators_par
        !> class encapsulating operators to compute
        !> the time derivatives of the main variables
        !> using the finite volume method on a
        !> distributed memory system
        !>
        !> @param compute_time_dev
        !> compute the time derivatives
        !---------------------------------------------------------------
        type, extends(td_operators_abstract_par) :: td_operators_par

          contains

          procedure, nopass :: compute_time_dev

        end type td_operators_par


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives using the
        !> space discretisation operators and the physical model
        !> on a distributed memory system
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !
        !>@param s
        !> space discretization operators
        !
        !>@param p_model
        !> physical model
        !
        !>@param bc_par_used
        !> boundary conditions for a distributed memory system
        !
        !>@param time_dev
        !> time derivatives
        !--------------------------------------------------------------
        function compute_time_dev(
     $       comm_2d, usr_rank,
     $       t, nodes, x_map, y_map, 
     $       s, p_model, bc_par_used)
     $       result(time_dev)

           implicit none

           integer                         , intent(in) :: comm_2d
           integer                         , intent(in) :: usr_rank
           real(rkind)                     , intent(in) :: t
           real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
           real(rkind), dimension(nx)      , intent(in) :: x_map
           real(rkind), dimension(ny)      , intent(in) :: y_map
           type(sd_operators)              , intent(in) :: s
           type(pmodel_eq)                 , intent(in) :: p_model
           type(bc_operators_par)          , intent(in) :: bc_par_used
           real(rkind), dimension(nx,ny,ne)             :: time_dev

           real(rkind)                        :: dx
           real(rkind)                        :: dy
           integer                            :: k
           integer(ikind)                     :: i,j
           real(rkind), dimension(nx+1,ny,ne) :: flux_x
           real(rkind), dimension(nx,ny+1,ne) :: flux_y


           !for debuggging: initialize the timedev with
           !a large initial value
           if(debug_initialize_timedev) then
              time_dev = reshape((/
     $             (((debug_real,i=1,nx),j=1,ny),k=1,ne)/),
     $             (/nx,ny,ne/))
           end if

           
           dx = x_map(2) - x_map(1)
           dy = y_map(2) - y_map(1)


           !<compute the fluxes
           !FORCEINLINE RECURSIVE
           flux_x = p_model%compute_flux_x(nodes,dx,dy,s)

           !FORCEINLINE RECURSIVE
           flux_y = p_model%compute_flux_y(nodes,dx,dy,s)


           !<if the boundary conditions influence the computation
           !> of the fluxes, then we need to modify the fluxes
           if(
     $          (bc_N_type_choice.eq.bc_fluxes_choice).or.
     $          (bc_S_type_choice.eq.bc_fluxes_choice).or.
     $          (bc_E_type_choice.eq.bc_fluxes_choice).or.
     $          (bc_W_type_choice.eq.bc_fluxes_choice)
     $     ) then
              call bc_par_used%apply_bc_on_fluxes(
     $             comm_2d, usr_rank,
     $             nodes, dx, dy, s, flux_x,flux_y)
           end if


           !<compute the time derivatives
           do k=1, ne
              do j=1+bc_size, ny-bc_size
                 do i=1+bc_size, nx-bc_size

                    time_dev(i,j,k)=
     $                   (flux_x(i,j,k)/dx-flux_x(i+1,j,k)/dx)+
     $                   (flux_y(i,j,k)/dy-flux_y(i,j+1,k)/dy)

                    time_dev(i,j,k)=time_dev(i,j,k)+
     $                   p_model%compute_body_forces(
     $                   t,x_map(i),y_map(j),
     $                   nodes(i,j,:),k)

                 end do
              end do
           end do


           !< if the boundary conditions influence the computation
           !> of the time derivatives, then we need to compute the
           !> time derivatives at the boundary
           if(
     $          (bc_N_type_choice.eq.bc_timedev_choice).or.
     $          (bc_S_type_choice.eq.bc_timedev_choice).or.
     $          (bc_E_type_choice.eq.bc_timedev_choice).or.
     $          (bc_W_type_choice.eq.bc_timedev_choice)
     $     ) then

             !in parallel not all the time derivatives of the
             !boundary : the selected time derivatives to be
             !computed are given by bc_sections
              bf_alignment(1,1) = bc_size+1
              bf_alignment(1,2) = nx-bc_size
              bf_alignment(2,1) = bc_size+1
              bf_alignment(2,2) = ny-bc_size
              
              call bc_par_used%apply_bc_on_timedev_nopt(
     $             t,
     $             bf_alignment,
     $             bf_grdpts_id,
     $             x_map,y_map,nodes,
     $             nodes,
     $             p_model,
     $             flux_x,flux_y,
     $             bc_sections,
     $             time_dev)

           end if
              
        end function compute_time_dev

      end module td_operators_par_class

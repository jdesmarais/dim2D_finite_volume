      !> @file
      !> class encapsulating the main temporary tables for the time 
      !> integration of the buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the main temporary tables for the time 
      !> integration of the buffer layer
      !
      !> @date
      !> 16_07_2014 - initial version         - J.L. Desmarais
      !> 15_10_2014 - interface modifications - J.L. Desmarais
      !> (to have a unique sd_operators, p_model... shared between the
      !>  field and the buffer layer objects)
      !-----------------------------------------------------------------
      module bf_compute_class

        use bf_layer_bc_procedure_module, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type

        use bc_operators_class, only :
     $       bc_operators

        use interface_integration_step, only :
     $       timeInt_step_nopt

        use parameters_constant, only :
     $       bc_timedev_choice

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne

        use parameters_input, only :
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use td_operators_class, only :
     $       td_operators


        implicit none

        private
        public :: bf_compute


        !>@class bf_compute
        !> class encapsulating the main temporary tables for the time 
        !> integration of the buffer layer
        !
        !>@param bc_sections
        !> identification of the boundary layers
        !
        !>@param grdpts_id_tmp
        !> temporary array saving the grdpts_id at the previous time step
        !
        !>@param nodes_tmp
        !> temporary array whose size is the same as the array containing
        !> the grid points integrated in time
        !
        !>@param time_dev
        !> temporary array containing the time derivatives and whose size
        !> is the same as the array containing the grid points integrated
        !> in time
        !
        !>@param allocate_tables
        !> allocate the grdpts_id_tmp, nodes_tmp and time_dev attributes
        !
        !>@param deallocate_tables
        !> deallocate the grdpts_id_tmp, nodes_tmp and time_dev attributes
        !
        !>@param compute_time_dev
        !> compute the time_dev attribute
        !
        !>@param compute_integration_step
        !> compute the nodes and nodes_tmp using the integration
        !> procedure
        !---------------------------------------------------------------
        type :: bf_compute

          integer    , dimension(:,:)  , allocatable, private :: bc_sections

          integer    , dimension(:,:)  , allocatable, private :: grdpts_id_tmp
          real(rkind), dimension(:,:,:), allocatable, private :: nodes_tmp
          real(rkind), dimension(:,:,:), allocatable, private :: time_dev

          contains

          procedure, pass :: allocate_tables
          procedure, pass :: deallocate_tables

          procedure, pass :: determine_bc_sections

          procedure, pass :: compute_time_dev
          procedure, pass :: compute_integration_step

          procedure, pass :: get_time_dev !only for tests

        end type bf_compute

        contains


        !allocate the nodes_tmp and time_dev tables
        subroutine allocate_tables(this, size_x, size_y)

          implicit none

          class(bf_compute), intent(inout) :: this
          integer(ikind)   , intent(in)    :: size_x
          integer(ikind)   , intent(in)    :: size_y

          allocate(this%grdpts_id_tmp(size_x,size_y))
          allocate(this%nodes_tmp(size_x,size_y,ne))
          allocate(this%time_dev(size_x,size_y,ne))

        end subroutine allocate_tables


        !allocate the nodes_tmp and time_dev tables
        subroutine deallocate_tables(this)

          implicit none

          class(bf_compute), intent(inout) :: this

          deallocate(this%bc_sections)
          deallocate(this%grdpts_id_tmp)
          deallocate(this%nodes_tmp)
          deallocate(this%time_dev)

        end subroutine deallocate_tables        


        !compute the time derivatives
        subroutine compute_time_dev(
     $     this,
     $     nodes,
     $     dx,
     $     dy,
     $     sd_operators_used,
     $     pmodel_eq_used,
     $     bc_operators_used,
     $     td_operators_used,
     $     grdpts_id)

          implicit none

          class(bf_compute)            , intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          real(rkind)                  , intent(in)    :: dx
          real(rkind)                  , intent(in)    :: dy
          type(sd_operators)           , intent(in)    :: sd_operators_used
          type(pmodel_eq)              , intent(in)    :: pmodel_eq_used
          type(bc_operators)           , intent(in)    :: bc_operators_used
          type(td_operators)           , intent(in)    :: td_operators_used
          integer    , dimension(:,:)  , intent(in)    :: grdpts_id
          
          call td_operators_used%compute_time_dev_nopt(
     $         nodes,
     $         dx,
     $         dy,
     $         sd_operators_used,
     $         pmodel_eq_used,
     $         bc_operators_used,
     $         this%time_dev,
     $         grdpts_id,
     $         this%bc_sections)

        end subroutine compute_time_dev


        !apply the boundary conditions on the time derivatives
        subroutine apply_bc_on_timedev(
     $     this,
     $     bc_used, p_model,
     $     t, nodes, x_map, y_map,
     $     flux_x, flux_y,
     $     timedev)

          implicit none

          class(bc_compute)            , intent(in)    :: this
          type(bc_operators)           , intent(in)    :: bc_used
          type(pmodel_eq)              , intent(in)    :: p_model
          real(rkind)                  , intent(in)    :: t
          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map
          real(rkind), dimension(:,:,:), intent(inout) :: flux_x
          real(rkind), dimension(:,:,:), intent(inout) :: flux_y
          real(rkind), dimension(:,:,:), intent(inout) :: timedev


          !if there are effectively boundary layers
          !in the buffer layer computed, the time
          !derivatives corresponding to the boundary
          !grid points are computed
          if(allocated(this%bc_sections)) then


             !compute the fluxes at the grid points
             !corresponding to edge-type boundary
             !layers
             call compute_edge_fluxes(this,flux_x,flux_y)


             !compute the time derivatives
             !corresponding to the buffer layers


          end if          

        end subroutine apply_bc_on_timedev


        subroutine compute_edge_fluxes(
     $     this,
     $     bc_used,
     $     p_model,
     $     nodes, dx, dy,
     $     flux_x, flux_y)

          implicit none
          
          class(bc_compute)            , intent(in)    :: this
          type(bc_operators)           , intent(in)    :: bc_used
          type(pmodel_eq)              , intent(in)    :: p_model
          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          real(rkind)                  , intent(in)    :: dx
          real(rkind)                  , intent(in)    :: dy
          real(rkind), dimension(:,:,:), intent(inout) :: flux_x
          real(rkind), dimension(:,:,:), intent(inout) :: flux_y


          !go through the boundary layers
          !if the boundary actually needs the computation
          !of the fluxes in the direction of the edge, the
          !fluxes are computed
          do k=1, size(this%bc_section,2)

             !identify the type of boundary layer
             select case(this%bc_section(1,k))

               case(N_edge_type)

                  !do not compute the edge fluxes only if
                  !(y.ge.bc_y_max).and.
                  !(bc_N_type_choice.ne.bc_timedev_choice)

                  j = this%bc_section(3,k)

                  compute_edge = .not.(
     $                 (y_map(j).ge.y_max).and.
     $                 (bc_N_type_choice.ne.bc_timedev_choice))

                  !determine the extent of the edge from the
                  !bc_section and compute the fluxes
                  if(compute_edge) then
                  
                     i_min = this%bc_section(2,k)
                     i_max = this%bc_section(4,k)

                     call bc_used%compute_fluxes_for_bc_y_edge(
     $                    p_model,
     $                    nodes,
     $                    s_y_L0, s_y_L1,
     $                    s_y_R1, s_y_R0,
     $                    dx, dy,
     $                    i_min, i_max, j,
     $                    N,
     $                    flux_x)

                  end if
                     
               case(S_edge_type)

                  !do not compute the edge fluxes only if
                  !(y.le.bc_y_min).and.
                  !(bc_S_type_choice.ne.bc_timedev_choice)

                  j = this%bc_section(3,k)

                  compute_edge = .not.(
     $                 (y_map(j+1).le.y_min).and.
     $                 (bc_S_type_choice.ne.bc_timedev_choice))

                  !determine the extent of the edge from the
                  !bc_section and compute the fluxes
                  if(compute_edge) then
                  
                     i_min = this%bc_section(2,k)
                     i_max = this%bc_section(4,k)

                     call bc_used%compute_fluxes_for_bc_y_edge(
     $                    p_model,
     $                    nodes,
     $                    s_y_L0, s_y_L1,
     $                    s_y_R1, s_y_R0,
     $                    dx, dy,
     $                    i_min, i_max, j,
     $                    S,
     $                    flux_x)

                  end if

               case(E_edge_type)

                  !do not compute the edge fluxes only if
                  !(x.ge.bc_x_max).and.
                  !(bc_E_type_choice.ne.bc_timedev_choice)

                  i = this%bc_section(2,k)

                  compute_edge = .not.(
     $                 (x_map(i).ge.y_max).and.
     $                 (bc_E_type_choice.ne.bc_timedev_choice))

                  !determine the extent of the edge from the
                  !bc_section and compute the fluxes
                  if(compute_edge) then
                  
                     j_min = this%bc_section(3,k)
                     j_max = this%bc_section(4,k)

                     call bc_used%compute_fluxes_for_bc_x_edge(
     $                    p_model,
     $                    nodes,
     $                    s_x_L0, s_x_L1,
     $                    s_x_R1, s_x_R0,
     $                    dx, dy,
     $                    j_min, j_max, i,
     $                    E,
     $                    flux_y)

                  end if

               case(W_edge_type)

                  !do not compute the edge fluxes only if
                  !(x.le.bc_x_min).and.
                  !(bc_W_type_choice.ne.bc_timedev_choice)

                  i = this%bc_section(2,k)

                  compute_edge = .not.(
     $                 (x_map(i+1).le.x_min).and.
     $                 (bc_W_type_choice.ne.bc_timedev_choice))

                  !determine the extent of the edge from the
                  !bc_section and compute the fluxes
                  if(compute_edge) then
                  
                     j_min = this%bc_section(3,k)
                     j_max = this%bc_section(4,k)

                     call bc_used%compute_fluxes_for_bc_x_edge(
     $                    p_model,
     $                    nodes,
     $                    s_x_L0, s_x_L1,
     $                    s_x_R1, s_x_R0,
     $                    dx, dy,
     $                    j_min, j_max, i,
     $                    W,
     $                    flux_y)

                  end if

             end select

          end do




        end subroutine compute_edge_fluxes


        !compute the integration step
        !det_bc_sections : determine the boundary sections
        subroutine compute_integration_step(
     $     this,
     $     grdpts_id, nodes, dt,
     $     integration_step_nopt)

          implicit none

          class(bf_compute)            , intent(inout) :: this
          integer    , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind)                  , intent(in)    :: dt
          procedure(timeInt_step_nopt)                 :: integration_step_nopt

          call integration_step_nopt(
     $         nodes, dt, this%nodes_tmp, this%time_dev, grdpts_id,
     $         this%bc_sections)

        end subroutine compute_integration_step


        !get the time_dev attribute
        subroutine get_time_dev(this, time_dev)

          implicit none

          class(bf_compute)                         , intent(in)  :: this
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: time_dev


          if(allocated(this%time_dev)) then
             allocate(time_dev(
     $            size(this%time_dev,1),
     $            size(this%time_dev,2),
     $            size(this%time_dev,3)))

             time_dev = this%time_dev

          else
             print '(''bf_compute_class'')'
             print '(''get_time_dev'')'
             stop 'time dev not allocated'

          end if

        end subroutine get_time_dev        

      end module bf_compute_class

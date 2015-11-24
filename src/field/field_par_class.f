      !> @file
      !> class extending the field class to encapsulate mpi
      !> attributes that will allows the computational domain
      !> simulated on one processor to communicate with the
      !> other processors computing the neighboring regions
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class extending the field class to encapsulate mpi
      !> attributes that will allows the computational domain
      !> simulated on one processor to communicate with the
      !> other processors computing the neighboring regions
      !
      !> @date
      !> 07_04_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module field_par_class
      
        use field_class, only :
     $       field

        use io_operators_par_class, only :
     $       io_operators_par

        use interface_integration_step, only :
     $       timeInt_step

        use mpi

        use mpi_interface_class, only :
     $       mpi_interface

        use mpi_process_class, only :
     $       mpi_process

        use parameters_constant, only :
     $       hedstrom_xy_choice,
     $       hedstrom_xy_corners_choice,
     $       hedstrom_x_reflection_y_choice,
     $       poinsot_xy_choice,
     $       yoolodato_xy_choice

        use parameters_input, only :
     $       nx,ny,ne,
     $       npx,npy,
     $       ntx,nty,
     $       bc_size,
     $       x_min,x_max,
     $       y_min,y_max,
     $       bc_choice,
     $       io_onefile_per_proc

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: field_par


        !> @class field_par
        !> class extending the field class to encapsulate the
        !> attributes needed to exchange data between processors
        !> computing the entire computational domain, one processor=
        !> one tile
        !---------------------------------------------------------------
        type, extends(field) :: field_par

          integer             :: comm2d                                !<@brief identify the mpi main communicator between the tiles
          integer             :: usr_rank                              !<@brief identify the processor ID computing the tile
          type(mpi_interface) :: mpi_interface_used                    !<@brief operator to exchange data between processors

          type(io_operators_par) :: io_operators_par_used              !<@brief operator to collect data from multiple processors and write the outputs

          integer(ikind), dimension(:,:), allocatable :: bc_sections_x !<@brief identify how the boundary points are computed for the E and W boundary regions
          integer(ikind), dimension(:,:), allocatable :: bc_sections_y !<@brief identify how the boundary points are computed for the N and S boundary regions

          contains

          procedure, pass :: ini               !<@brief initialize the computational field and the mpi attributes
          procedure, pass :: ini_coordinates   !<@brief define the coordinate maps depending on the position of the tile computed
          procedure, pass :: apply_bc_on_nodes !<@brief compute the boundary points of the domain (solve PDE, exchange data with neighboring tile...)
          procedure, pass :: write_data        !<@brief write the data related to the tile computed in the file for the entire computational domain

        end type field_par


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the field_par by:
        !> -# check the simulation inputs
        !> -# initialize the cartesian communicator between the tiles
        !> -# initializing the boundary conditions bc_operators_used
        !> -# initializing the coordinates
        !> -# applying the initial conditions
        !> -# initializing the i/o operators io_operators_used
        !
        !> @date
        !> 07_04_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(field_par), intent(inout) :: this

          type(mpi_process) :: mpi_process_used
          
          integer(ikind), dimension(2) :: x_borders
          integer(ikind), dimension(2) :: y_borders


          !1) check the inputs
          call this%check_inputs()

          !2) initialize the cartesian communicator
          call mpi_process_used%ini_cartesian_communicator(
     $         this%comm2d,
     $         this%usr_rank)

          !3) initialize the mpi_interface to identify the
          !   processor computing the neighboring tiles and 
          !   data that will be exchanged b/w the tiles
          call this%mpi_interface_used%ini(
     $         this%comm2d)

          !4) initialize the integration borders
          call this%mpi_interface_used%ini_for_timeInt(
     $         this%comm2d,
     $         x_borders,
     $         y_borders,
     $         this%bc_sections_x,
     $         this%bc_sections_y)

          call this%ini_for_timeInt()

          this%x_borders(1) = max(this%x_borders(1),x_borders(1))
          this%x_borders(2) = min(this%x_borders(2),x_borders(2))
          this%y_borders(1) = max(this%y_borders(1),y_borders(1))
          this%y_borders(2) = min(this%y_borders(2),y_borders(2))

          !5) initialize the boundary conditions
          call this%bc_operators_used%ini(this%pmodel_eq_used)

          !6) initialize the time+x_map,y_map+nodes+io_operators
          this%time = 0.0

          call this%ini_coordinates()
          call this%apply_initial_conditions()

          call this%io_operators_used%ini()
          call this%io_operators_par_used%ini(this%comm2d,this%usr_rank)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to initialize the coordinates along the
        !> x and y directions of the tile
        !
        !> @date
        !> 07_04_20135- initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !--------------------------------------------------------------
        subroutine ini_coordinates(this)

          implicit none

          class(field_par), intent(inout) :: this

          type(mpi_process)     :: mpi_op
          integer               :: dims_nb
          integer, dimension(2) :: cart_coord
          integer(ikind)        :: i,j
          integer               :: ierror
          real(rkind)           :: x_min_tile
          real(rkind)           :: y_min_tile


          ! initialize the space steps along the 
          ! x and y directions
          if(npx.gt.1) then
             this%dx = (x_max-x_min)/(npx*(nx-2*bc_size)-1)
          else
             this%dx = (x_max-x_min)/(ntx-1-2*bc_size)
          end if

          if(npy.gt.1) then
             this%dy = (y_max-y_min)/(npy*(ny-2*bc_size)-1)
          else
             this%dy = (y_max-y_min)/(nty-1-2*bc_size)
          end if


          ! get the cartesian coordinates of the tile
          dims_nb=2
          call MPI_CART_COORDS(
     $         this%comm2d, this%usr_rank,
     $         dims_nb, cart_coord,
     $         ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             print '(''field_par_class'')'
             print '(''ini_coordinates'')'
             stop 'MPI_CART_RANK failed'
          end if


          ! get the x_min corresponding to the tile
          if(cart_coord(1).eq.0) then
             x_min_tile = x_min
          else
             x_min_tile =
     $            x_min + cart_coord(1)*(nx-2*bc_size)*this%dx
          end if


          ! get the y_min corresponding to the tile
          if(cart_coord(2).eq.0) then
             y_min_tile = y_min
          else
             y_min_tile =
     $            y_min + cart_coord(2)*(ny-2*bc_size)*this%dy
          end if


          ! initialize the coordinates along the
          ! x-direction
          do i=1, nx
             this%x_map(i)=x_min_tile + (i-1-bc_size)*this%dx
          end do


          ! initialize the coordinates along the
          ! y-direction
          do j=1, ny
             this%y_map(j)=y_min_tile + (j-1-bc_size)*this%dy
          end do

        end subroutine ini_coordinates


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the boundary conditions on the grid points
        !> : solve PDE, exchange with neighboring tiles...
        !
        !> @date
        !> 07_04_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(this)

          implicit none

          class(field_par), intent(inout) :: this

          
          !1) create requests for exchanging data in the
          !   x-direction
          call this%mpi_interface_used%MPI_ISENDRECV_XDIR(
     $         this%comm2d,
     $         this%nodes)

          !2) overlap the communications of boundary points 
          !   in the x-direction by computation of the other
          !   boundary points in the x-direction
          if((bc_choice.ne.hedstrom_xy_choice).and.
     $       (bc_choice.ne.hedstrom_xy_corners_choice).and.
     $       (bc_choice.ne.poinsot_xy_choice).and.
     $       (bc_choice.ne.yoolodato_xy_choice)) then

             call this%bc_operators_used%apply_bc_on_nodes_nopt(
     $            this%nodes, this%bc_sections_x)

          end if

          !3) wait for the MPI requests in the x-direction
          call this%mpi_interface_used%MPI_WAITALL_XDIR()


          !4) create requests for exchanging data in the
          !   y-direction
          call this%mpi_interface_used%MPI_ISENDRECV_YDIR(
     $         this%comm2d,
     $         this%nodes)

          !5) overlap the communications of boundary points 
          !   in the y-direction by computation of the other
          !   boundary points in the y-direction
          if((bc_choice.ne.hedstrom_xy_choice).and.
     $       (bc_choice.ne.hedstrom_xy_corners_choice).and.
     $       (bc_choice.ne.poinsot_xy_choice).and.
     $       (bc_choice.ne.yoolodato_xy_choice)) then

             call this%bc_operators_used%apply_bc_on_nodes_nopt(
     $            this%nodes, this%bc_sections_y)

          end if

          !6) wait for the MPI requests in the y-direction
          call this%mpi_interface_used%MPI_WAITALL_YDIR()

        end subroutine apply_bc_on_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> write the interior domain data on output files depending
        !> on the i/o operators used
        !
        !> @date
        !> 07_04_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !--------------------------------------------------------------
        subroutine write_data(this)

          implicit none

          class(field_par), intent(inout) :: this

          if(io_onefile_per_proc) then

             call this%io_operators_used%write_data(
     $            this%nodes,
     $            this%x_map,
     $            this%y_map,
     $            this%pmodel_eq_used,
     $            this%time,
     $            rank=this%usr_rank)

          else

             call this%io_operators_par_used%write_data(
     $            this%comm2d,
     $            this%nodes,
     $            this%x_map,
     $            this%y_map,
     $            this%pmodel_eq_used,
     $            this%time)

          end if

        end subroutine write_data      

      end module field_par_class

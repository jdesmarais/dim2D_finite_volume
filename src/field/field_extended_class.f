      !> @file
      !> class encapsulating the main tables for the variables and the
      !> coordinates with a possibility to extend the computational domain
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the main tables for the variables and the
      !> coordinates with a possibility to extend the computational domain
      !
      !> @date
      ! 17_07_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module field_extended_class

        use bf_interface_class, only :
     $       bf_interface

        use cmd_operators_class, only :
     $       cmd_operators

        use field_abstract_class, only :
     $       field_abstract

        use interface_integration_step, only :
     $       timeInt_step,
     $       timeInt_step_nopt

        use parameters_input, only :
     $       nx,ny,ne,
     $       write_domain_extension,
     $       debug_restart_for_geometry,
     $       debug_adapt_computational_domain

        use parameters_kind, only :
     $       ikind,rkind

        use td_integrator_class, only :
     $       td_integrator

        implicit none


        private
        public :: field_extended


        !>@class field_extended
        !> object combining the interior domain and its extension
        !
        !>@param domain_extension
        !> object containing the domain extension and how data are
        !> exchanged between the interior computational domain and
        !> the buffer layers
        !
        !>@param td_integrator_used
        !> time integration method
        !
        !>@param ini
        !> initialize the interior computational domain and its
        !> extension
        !
        !>@param apply_initial_conditions
        !> apply the initial conditions of the physical model on the
        !> interior domain and its extension if any
        !
        !>@param compute_time_dev_ext
        !> compute the time derivatives of the interior computational
        !> domain and the domain extension
        !
        !>@param compute_integration_step_ext
        !> compute the time integration step for the interior
        !> computational domain and its extension
        !
        !>@param integrate
        !> integrate the interior computational domain and its extension
        !> in time
        !
        !>@param write_data
        !> write the data contained in the interior domain and its
        !> extension on external netcdf files
        !
        !>@param adapt_domain
        !> adapt the computational domain by increasing or decreasing
        !> the buffer layers depending on the activation of the grid-points
        !> in the computational domain
        !---------------------------------------------------------------
        type, extends(field_abstract) :: field_extended

          type(bf_interface)                          :: domain_extension
          type(td_integrator)                         :: td_integrator_used
          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections

          contains

          procedure, pass :: ini
          procedure, pass :: apply_initial_conditions
          procedure, pass :: compute_time_dev_ext
          procedure, pass :: compute_integration_step_ext
          procedure, pass :: integrate
          procedure, pass :: write_data
          procedure, pass :: adapt_domain

        end type field_extended


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the interior computational domain and its
        !> extension
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables        
        !--------------------------------------------------------------
        subroutine ini(this)
        
          implicit none

          class(field_extended)  , intent(inout) :: this

          type(cmd_operators)   :: cmd_operators_used
          integer, dimension(4) :: nb_bf_layers


          !initialize the interior domain
          call this%field_abstract%ini()

          !analyse the command line arguments
          call cmd_operators_used%analyse_cmd_line_arg()

          !initialize the domain extension
          if(cmd_operators_used%is_restart_activated()) then

             nb_bf_layers = cmd_operators_used%get_nb_bf_layers()

             print '(''field_extended_class'')'
             stop 'domain extension restart not implemented'

c$$$             call this%domain_extension%restart(
c$$$     $            this%x_map,
c$$$     $            this%y_map,
c$$$     $            this%nodes,
c$$$     $            nb_bf_layers,
c$$$     $            this%pmodel_eq_used,
c$$$     $            this%io_operators_used%get_nb_timesteps_written())

             !if the restart was only used to have the geometry
             !of the previous computational domain, the initial
             !conditions should be re-applied on the complete
             !computational domain
             if(debug_restart_for_geometry) then
                call this%apply_initial_conditions()
             end if
             
          else

             call this%domain_extension%ini(this%x_map,this%y_map)

          end if

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the initial conditions of the physical model on
        !> the entire field (interior+buffer layers)
        !
        !> @date
        !> 23_01_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !--------------------------------------------------------------
        subroutine apply_initial_conditions(this)

          implicit none

          class(field_extended), intent(inout) :: this
          
          !interior domain
          call this%field_abstract%apply_initial_conditions()

          !buffer layers
          call this%domain_extension%apply_initial_conditions(
     $         this%pmodel_eq_used)

        end subroutine apply_initial_conditions


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the interior field
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !
        !>@return time_dev
        !> array containing the time derivatives of the grid points
        !--------------------------------------------------------------
        function compute_time_dev(this) result(time_dev)

          implicit none

          class(field_extended), intent(in) :: this
          real(rkind), dimension(nx,ny,ne)  :: time_dev

          !make use of the time discretization operator
          !to compute the time derivative of the field
          time_dev = this%td_operators_used%compute_time_dev(
     $         this%time,
     $         this%nodes,
     $         this%x_map,
     $         this%y_map,
     $         this%sd_operators_used,
     $         this%pmodel_eq_used,
     $         this%bc_operators_used,
     $         bc_sections=this%interior_bc_sections)

          !the boundary conditions are applied on the time
          !derivatives of the interior by using the local
          !bc_sections for the interior

        end function compute_time_dev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the interior computational
        !> domain and the domain extension
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@return time_dev
        !> time derivative of the interior domain
        !--------------------------------------------------------------
        function compute_time_dev_ext(this) result(time_dev)

          implicit none

          class(field_extended), intent(inout) :: this
          real(rkind), dimension(nx,ny,ne)     :: time_dev

  
          !compute the time derivatives of the interior domain
          !(the boundary conditions on the time derivatives are
          ! applied on the interior domain using the interior
          ! bc_sections)
          time_dev = compute_time_dev(this)


          !compute the time derivatives of the domain extension
          !(the boundary conditions on the time derivatives are
          ! applied on the buffer layers using the bc_sections
          ! for each buffer layer)
          call this%domain_extension%compute_time_dev(
     $         this%td_operators_used,
     $         this%time,
     $         this%sd_operators_used,
     $         this%pmodel_eq_used,
     $         this%bc_operators_used,
     $         this%nodes)

        end function compute_time_dev_ext


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the integration step of the interior domain and its
        !> extension
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@param dt
        !> time step
        !
        !>@param nodes_tmp
        !> array containing the temporary grid points for the
        !> time integration of the interior computational domain
        !
        !>@param time_dev
        !> time derivatives of the interior domain
        !
        !>@param integration_step
        !> procedure for the time integration of the interior domain
        !
        !>@param integration_step_nopt
        !> procedure for the time integration of the domain extension
        !
        !>@param full
        !> logical to enforce whether the full domain is computed
        !> discarding the x_borders and y_borders supplied
        !> (important for the first integration step when the
        !> previous integration step is saved temporary in another
        !> array)
        !--------------------------------------------------------------
        subroutine compute_integration_step_ext(
     $     this, dt, nodes_tmp, time_dev,
     $     integration_step, integration_step_nopt)

          implicit none

          class(field_extended)           , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          procedure(timeInt_step)                         :: integration_step
          procedure(timeInt_step_nopt)                    :: integration_step_nopt

          
          !compute the integration step for the interior domain
          call this%field_abstract%compute_integration_step(
     $         dt, nodes_tmp, time_dev, integration_step)

                    
          !compute the integration step for the domain extension
          call this%domain_extension%compute_integration_step(
     $         dt, integration_step_nopt, this%nodes)

        end subroutine compute_integration_step_ext


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> integrate the interior computational domain and its extension
        !> in time
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@param dt
        !> time step
        !--------------------------------------------------------------
        subroutine integrate(this, dt)

          implicit none

          class(field_extended), intent(inout) :: this
          real(rkind)          , intent(in)    :: dt


          !allocate memory space for the temporary tables
          !used in the time integration of the domain extension
          call this%domain_extension%initialize_before_timeInt(
     $         this%interior_bc_sections)

          !integrate the computational domain: interior + extension
          call this%td_integrator_used%integrate_ext(this,dt)

          !deallocate memory space for the temporary tables
          !used in the time integration of the domain extension
          call this%domain_extension%finalize_after_timeInt()


          !update the time
          this%time = this%time + dt

        end subroutine integrate


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the boundary conditions on the gridpoints
        !
        !> @date
        !> 03_11_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the main variables
        !--------------------------------------------------------------
        subroutine write_data(this)

          implicit none

          class(field_extended), intent(inout) :: this

          character(len=10), dimension(ne) :: name_var
          character(len=33), dimension(ne) :: longname_var
          character(len=23), dimension(ne) :: unit_var

          integer :: nb_timesteps


          !write the data from the domain extension if set by the user
          if(write_domain_extension) then

             name_var     = this%pmodel_eq_used%get_var_name()
             longname_var = this%pmodel_eq_used%get_var_longname()
             unit_var     = this%pmodel_eq_used%get_var_unit()
             
             nb_timesteps = this%io_operators_used%get_nb_timesteps_written()

             call this%domain_extension%print_netcdf(
     $            nb_timesteps,
     $            name_var,
     $            longname_var,
     $            unit_var,
     $            this%time)
             
          end if

          !write the interior data using the function encapsulated
          !in field_abstract
          call this%field_abstract%write_data()

        end subroutine write_data


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> adapt the computational domain
        !
        !> @date
        !> 22_11_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the main variables at t
        !
        !> @param t
        !> time
        !
        !> @param dt
        !> time step
        !
        !> @param nodes0
        !> nodes at the previous time step (t-dt)
        !--------------------------------------------------------------
        subroutine adapt_domain(this, dt, nodes0)

          implicit none

          class(field_extended)           , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes0


          if(debug_adapt_computational_domain) then

             !allocate memory space for the temporary tables
             !used in the time integration of the domain extension
             call this%domain_extension%adapt_domain_extension(
     $            this%x_map,
     $            this%y_map,
     $            nodes0,
     $            this%nodes,
     $            this%pmodel_eq_used,
     $            this%time,
     $            dt,
     $            this%interior_bc_sections)

          end if

        end subroutine adapt_domain

      end module field_extended_class

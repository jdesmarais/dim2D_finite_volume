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

        use bf_interface_dcr_class    , only : bf_interface_dcr
        use field_abstract_class      , only : field_abstract
        use interface_integration_step, only : timeInt_step,
     $                                         timeInt_step_nopt
        use parameters_input          , only : nx,ny,ne
        use parameters_kind           , only : rkind
        use td_integrator_class       , only : td_integrator

        implicit none


        private
        public :: field_extended


        !>@class field
        !> class encapsulating the variables of the governing equations
        !> and the discretisation maps
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
        !---------------------------------------------------------------
        type, extends(field_abstract) :: field_extended

          type(bf_interface_dcr) :: domain_extension
          type(td_integrator)    :: td_integrator_used

          contains

          procedure, pass :: ini
          procedure, pass :: compute_time_dev_ext
          procedure, pass :: compute_integration_step_ext
          procedure, pass :: integrate

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

          class(field_extended), intent(inout) :: this

          !initialize the interior domain
          call this%field_abstract%ini()

          !initialize the domain extension
          call this%domain_extension%ini()

        end subroutine ini


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
          time_dev = this%field_abstract%compute_time_dev()


          !compute the time derivatives of the domain extension
          call this%domain_extension%compute_time_dev()

          
          !WARNING: depending on the way to implement the computation
          !of the boundary conditions on the time derivatives, it may
          !be needed to add the computation of the boundary conditions
          !here and to replace this%field_abstract%compute_time_dev()
          !by its implementation to distinguish what is computed from
          !what is exchanged for the interior domain
          print '()'
          print '(''**********************************************'')'
          print '(''field_extended_class'')'
          print '(''compute_time_dev'')'
          print '(''think on how to implement the b.c. on time dev'')'
          print '(''**********************************************'')'
          print '()'

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
        !--------------------------------------------------------------
        subroutine compute_integration_step_ext(
     $     this, dt, nodes_tmp, time_dev,
     $     integration_step, integration_step_nopt)

          implicit none

          class(field_extended)           , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          procedure(timeInt_step)      :: integration_step
          procedure(timeInt_step_nopt) :: integration_step_nopt

          
          !compute the integration step for the interior domain
          call this%field_abstract%compute_integration_step(
     $         dt, nodes_tmp, time_dev, integration_step)
                    
          !compute the integration step for the domain extension
          call this%domain_extension%compute_integration_step(
     $         dt, integration_step_nopt)

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
          call this%domain_extension%allocate_before_timeInt()

          !integrate the computational domain: interior + extension
          call this%td_integrator_used%integrate_ext(this,dt)

          !deallocate memory space for the temporary tables
          !used in the time integration of the domain extension
          call this%domain_extension%deallocate_after_timeInt()

        end subroutine integrate

      end module field_extended_class

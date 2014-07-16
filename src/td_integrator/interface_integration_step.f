      module interface_integration_step

        use parameters_kind, only : rkind

        implicit none

        private
        public :: timeInt_step_nopt


        abstract interface

          subroutine timeInt_step_nopt(
     $          nodes, dt, nodes_tmp, time_dev, grdpts_id)
           
             import rkind

             real(rkind), dimension(:,:,:), intent(inout) :: nodes
             real(rkind)                  , intent(in)    :: dt
             real(rkind), dimension(:,:,:), intent(inout) :: nodes_tmp
             real(rkind), dimension(:,:,:), intent(in)    :: time_dev
             integer    , dimension(:,:)  , intent(in)    :: grdpts_id

          end subroutine timeInt_step_nopt

        end interface

      end module interface_integration_step

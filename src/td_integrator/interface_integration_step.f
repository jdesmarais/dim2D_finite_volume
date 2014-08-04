      module interface_integration_step

        use parameters_input, only : nx,ny,ne
        use parameters_kind , only : ikind,rkind

        implicit none

        private
        public :: timeInt_step,
     $            timeInt_step_nopt


        abstract interface

          subroutine timeInt_step(
     $       nodes, dt, nodes_tmp, time_dev,
     $       x_borders, y_borders)
           
             import ikind,rkind
             import nx,ny,ne

             real(rkind), dimension(nx,ny,ne)      , intent(inout) :: nodes
             real(rkind)                           , intent(in)    :: dt
             real(rkind), dimension(nx,ny,ne)      , intent(inout) :: nodes_tmp
             real(rkind), dimension(nx,ny,ne)      , intent(in)    :: time_dev
             integer(ikind), dimension(2), optional, intent(in)    :: x_borders
             integer(ikind), dimension(2), optional, intent(in)    :: y_borders

          end subroutine timeInt_step


          subroutine timeInt_step_nopt(
     $     nodes, dt, nodes_tmp, time_dev, grdpts_id,
     $     x_borders, y_borders)
           
             import ikind,rkind

             real(rkind), dimension(:,:,:)         , intent(inout) :: nodes
             real(rkind)                           , intent(in)    :: dt
             real(rkind), dimension(:,:,:)         , intent(inout) :: nodes_tmp
             real(rkind), dimension(:,:,:)         , intent(in)    :: time_dev
             integer    , dimension(:,:)           , intent(in)    :: grdpts_id
             integer(ikind), dimension(2), optional, intent(in)    :: x_borders
             integer(ikind), dimension(2), optional, intent(in)    :: y_borders

          end subroutine timeInt_step_nopt

        end interface

      end module interface_integration_step

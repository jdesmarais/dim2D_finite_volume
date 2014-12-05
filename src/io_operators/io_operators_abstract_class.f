      !> @file
      !> class encapsulating subroutines to write
      !> data on output files (ex:netcdf files)
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to write
      !> data on output files (ex:netcdf files)
      !
      !> @date
      !> 14_08_2013 - initial version              - J.L. Desmarais
      !> 15_07_2014 - composition over inheritance - J.L. Desmarais
      !-----------------------------------------------------------------
      module io_operators_abstract_class

        use parameters_input, only : nx,ny,ne,bc_size
        use parameters_kind , only : rkind
        use pmodel_eq_class , only : pmodel_eq

        implicit none
        
        private
        public  :: io_operators_abstract


        !> @class io_operators_abstract
        !> class encapsulating subroutines to write
        !> data on output files (ex:netcdf files)
        !
        !> @param ini
        !> initialize the writer
        !
        !> @param write_data
        !> write data on output files
        !
        !>@param read_data
        !> read data from output files
        !---------------------------------------------------------------
        type, abstract :: io_operators_abstract

          integer :: nb_timesteps_written

          contains
          
          procedure            , pass           :: ini
          procedure(write_proc), pass, deferred :: write_data
          procedure(read_proc) , pass, deferred :: read_data
          procedure            , pass           :: increment_counter
          procedure            , pass           :: get_nb_timesteps_written
          procedure            , pass           :: set_nb_timesteps_written

        end type io_operators_abstract


        interface

           subroutine write_proc(this,nodes,x_map,y_map,p_model,time)

             import io_operators_abstract
             import pmodel_eq
             import nx,ny,ne
             import rkind

             class(io_operators_abstract)    , intent(inout) :: this
             real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
             real(rkind), dimension(nx)      , intent(in)    :: x_map
             real(rkind), dimension(ny)      , intent(in)    :: y_map
             type(pmodel_eq)                 , intent(in)    :: p_model
             real(rkind)                     , intent(in)    :: time

           end subroutine write_proc

      
           subroutine read_proc(
     $       this,
     $       filename,
     $       nodes,
     $       x_map,
     $       y_map,
     $       p_model,
     $       time)

             import io_operators_abstract
             import pmodel_eq
             import nx,ny,ne
             import rkind

             class(io_operators_abstract)    , intent(inout) :: this
             character*(*)                   , intent(in)    :: filename
             real(rkind), dimension(nx,ny,ne), intent(out)   :: nodes
             real(rkind), dimension(nx)      , intent(out)   :: x_map
             real(rkind), dimension(ny)      , intent(out)   :: y_map
             type(pmodel_eq)                 , intent(in)    :: p_model
             real(rkind)                     , intent(out)   :: time

           end subroutine read_proc

        end interface


        contains


        !initialize the writer
        subroutine ini(this, counter)

          implicit none

          class(io_operators_abstract), intent(inout) :: this
          integer           , optional, intent(in)    :: counter

          if(present(counter)) then
             this%nb_timesteps_written=counter
          else
             this%nb_timesteps_written=0
          end if
          
        end subroutine ini


        !increment the counter for tracking the number of timesteps
        !written
        subroutine increment_counter(this)

          implicit none

          class(io_operators_abstract), intent(inout) :: this

          this%nb_timesteps_written=this%nb_timesteps_written+1

        end subroutine increment_counter

        
        !get the nb_timesteps_written attribute
        function get_nb_timesteps_written(this) result(nb_timesteps_written)
        
          implicit none

          class(io_operators_abstract), intent(in) :: this
          integer                                  :: nb_timesteps_written

          nb_timesteps_written = this%nb_timesteps_written

        end function get_nb_timesteps_written


        !set the nb_timesteps_written attribute
        subroutine set_nb_timesteps_written(this,nb_timesteps_written)
        
          implicit none

          class(io_operators_abstract), intent(inout) :: this
          integer                     , intent(in)    :: nb_timesteps_written

          this%nb_timesteps_written = nb_timesteps_written

        end subroutine set_nb_timesteps_written

      end module io_operators_abstract_class

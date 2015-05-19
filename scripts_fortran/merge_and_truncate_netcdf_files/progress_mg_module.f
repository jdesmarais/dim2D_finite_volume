      module progress_mg_module

        implicit none

        private
        public :: 
     $       ini_progress_mg,
     $       update_progress_mg,
     $       finalize_progress_mg

        contains

        subroutine ini_progress_mg(unit_mg_progress)

          implicit none

          integer, intent(in) :: unit_mg_progress

          open(unit_mg_progress, carriagecontrol ='fortran')
          write(unit_mg_progress,'(1h+"files merged: ...")')

        end subroutine ini_progress_mg


        subroutine update_progress_mg(unit_mg_progress,file_merged,nb_files)

          implicit none

          integer, intent(in) :: unit_mg_progress
          integer, intent(in) :: file_merged
          integer, intent(in) :: nb_files

          write(unit_mg_progress,'(1h+"files merged: ",i4," /",i4)') file_merged, nb_files

        end subroutine update_progress_mg


        subroutine finalize_progress_mg(unit_mg_progress)

          implicit none

          integer, intent(in) :: unit_mg_progress

          write(unit_mg_progress,'(1h+"files merged: done          ")')

        end subroutine finalize_progress_mg


      end module progress_mg_module

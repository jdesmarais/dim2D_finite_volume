      module bf_restart_module

        use parameters_input, only :
     $       nx,ny,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public ::
     $       get_restart_alignment,
     $       get_nb_detectors,
     $       read_detectors_from_file


        contains


        !get the bf_alignment corresponding to the buffer layer
        function get_restart_alignment(
     $       interior_x_map,
     $       interior_y_map,
     $       x_borders,
     $       y_borders)
     $       result(bf_alignment)

          implicit none

          real(rkind), dimension(nx) , intent(in) :: interior_x_map
          real(rkind), dimension(ny) , intent(in) :: interior_y_map
          real(rkind), dimension(2)  , intent(in) :: x_borders
          real(rkind), dimension(2)  , intent(in) :: y_borders
          integer    , dimension(2,2)             :: bf_alignment


          bf_alignment(1,1) = get_align(
     $         x_borders(1),
     $         interior_x_map,
     $         nx)

          bf_alignment(1,2) = get_align(
     $         x_borders(2),
     $         interior_x_map,
     $         nx)

          bf_alignment(2,1) = get_align(
     $         y_borders(1),
     $         interior_y_map,
     $         ny)

          bf_alignment(2,2) = get_align(
     $         y_borders(2),
     $         interior_y_map,
     $         ny)

          bf_alignment(1,1) = bf_alignment(1,1)+bc_size
          bf_alignment(1,2) = bf_alignment(1,2)-bc_size
          bf_alignment(2,1) = bf_alignment(2,1)+bc_size
          bf_alignment(2,2) = bf_alignment(2,2)-bc_size

        end function get_restart_alignment


        function get_align(x_border,interior_x_map,size_x)
     $     result(i_border_x)

          implicit none

          real(rkind)              , intent(in) :: x_border
          real(rkind), dimension(:), intent(in) :: interior_x_map
          integer                  , intent(in) :: size_x
          integer                               :: i_border_x

          integer     :: i
          real(rkind) :: dx
          real(rkind) :: x1,x2

          if(x_border.lt.interior_x_map(1)) then

             i=0
             dx = interior_x_map(2)-interior_x_map(1)

             do while(x_border.lt.(interior_x_map(1)+(i-1)*dx))
                i = i-1
             end do

             x1 = interior_x_map(1)+(i-1)*dx
             x2 = interior_x_map(1)+i*dx

             if((abs(x_border-x1)).lt.(abs(x_border-x2))) then
                i_border_x = i
             else
                i_border_x = i+1
             end if

          else

             if(x_border.le.interior_x_map(size_x)) then

                i=1
                
                do while(x_border.ge.interior_x_map(i))
                   i=i+1
                   if(i.eq.size_x) then
                      exit
                   end if
                end do

                x1 = interior_x_map(i-1)
                x2 = interior_x_map(i)

                if((abs(x_border-x1)).lt.(abs(x_border-x2))) then
                   i_border_x = i-1
                else
                   i_border_x = i
                end if
                      
             else

                dx = interior_x_map(size_x) - interior_x_map(size_x-1)
                i  = size_x

                do while(x_border.gt.(interior_x_map(size_x)+(i-size_x)*dx))
                   i = i+1
                end do

                x1 = interior_x_map(size_x)+(i-1-size_x)*dx
                x2 = interior_x_map(size_x)+(i  -size_x)*dx

                if((abs(x_border-x1)).lt.(abs(x_border-x2))) then
                   i_border_x = i-1
                else
                   i_border_x = i
                end if

             end if

          end if

        end function get_align


        !get the number of detectors in the file
        function get_nb_detectors(filename)
     $     result(nb_detectors)

          implicit none

          character*(*), intent(in)   :: filename
          integer      , dimension(4) :: nb_detectors

          character(1) :: line
          integer      :: ios
          integer      :: k
          integer      :: nb_lines


          !open for the file for reading
          open(unit=1,
     $         file=filename,
     $         iostat=ios)
          if(ios.ne.0) then
             print *, 'error opening detector file'
          end if

          !read first line
          read(1,*,iostat=ios) line
          k=0
          nb_lines=0

          !read lines
          do while(ios.eq.0)

             if(line.eq.'#') then
                if(k>0) then
                   nb_detectors(k) = nb_lines-1
                   nb_lines=0
                end if
                k=k+1                
             end if

             read(1,*,iostat=ios) line

             if(ios.eq.0) then
                nb_lines=nb_lines+1
             else
                nb_detectors(k) = nb_lines-1
             end if
             
          end do
          nb_detectors(4) = nb_detectors(4)+1

          !close file
          close(unit=1)

        end function get_nb_detectors


        !read the detectors
        subroutine read_detectors_from_file(
     $     filename,
     $     N_detectors,
     $     S_detectors,
     $     E_detectors,
     $     W_detectors)

          implicit none

          character*(*)                           , intent(in)  :: filename
          real(rkind), dimension(:,:), allocatable, intent(out) :: N_detectors
          real(rkind), dimension(:,:), allocatable, intent(out) :: S_detectors
          real(rkind), dimension(:,:), allocatable, intent(out) :: E_detectors
          real(rkind), dimension(:,:), allocatable, intent(out) :: W_detectors

          
          integer, dimension(4) :: nb_detectors
          integer               :: ios
          character(1)          :: line


          !1) get the number of detectors saved
          nb_detectors = get_nb_detectors(filename)

          allocate(N_detectors(2,nb_detectors(1)))
          allocate(S_detectors(2,nb_detectors(2)))
          allocate(E_detectors(2,nb_detectors(3)))
          allocate(W_detectors(2,nb_detectors(4)))


          !2) read the detectors
          !open for the file for reading
          open(unit=1,
     $         file=filename,
     $         form='formatted',
     $         access='sequential',
     $         action='read',
     $         status='old',
     $         iostat=ios)
          if(ios.ne.0) then
             print *, 'error opening detector file'
          end if

          !skip first line
          read(1,*,iostat=ios) line

          !read N_detectors
          call read_detectors(N_detectors,nb_detectors(1))
          
          !skip intermediate lines
          read(1,*,iostat=ios) line

          !read S_detectors
          call read_detectors(S_detectors,nb_detectors(2))
          
          !skip intermediate lines
          read(1,*,iostat=ios) line

          !read E_detectors
          call read_detectors(E_detectors,nb_detectors(3))
          
          !skip intermediate lines
          read(1,*,iostat=ios) line

          !read W_detectors
          call read_detectors(W_detectors,nb_detectors(4))
          
          !close file
          close(unit=1)

        end subroutine read_detectors_from_file


        subroutine read_detectors(detectors,nb_detectors)

          implicit none

          real(rkind), dimension(:,:), intent(out) :: detectors
          integer                    , intent(in)  :: nb_detectors

          integer :: i
          integer :: ios
          
          do i=1, nb_detectors
             
             read(1,FMT='(2F8.4)',iostat=ios) detectors(:,i)

             if(ios.ne.0) then
                print '(''error reading detectors'')'
                stop ''
             end if
             
          end do

        end subroutine read_detectors
        

      end module bf_restart_module

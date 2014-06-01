      module bf_interface_icr_class

        use bf_nbc_template_module
        use bf_sublayer_class  , only : bf_sublayer
        use bf_interface_class , only : bf_interface
        use parameters_constant, only : interior
        use parameters_input   , only : nx,ny,ne,bc_size,
     $                                  dx,dy,dt,search_nb_dt


        implicit none

        private
        public :: bf_interface_icr


        type, extends(bf_interface) :: bf_interface_icr

          integer(ikind), dimension(:,:), allocatable :: N_detectors_list
          integer(ikind), dimension(:,:), allocatable :: S_detectors_list
          integer(ikind), dimension(:,:), allocatable :: E_detectors_list
          integer(ikind), dimension(:,:), allocatable :: W_detectors_list

          contains

          procedure,   pass :: ini
          procedure,   pass :: get_modified_grdpts_list

          procedure, nopass, private :: is_i_detector_activated
          procedure, nopass, private :: get_central_pt
          procedure,   pass, private :: check_neighboring_bc_interior_pts

          prcoedure, nopass, private :: create_nbc_interior_pt_template

        end type bf_interface_increase


        contains


        !< initialize the detector lists
        subroutine ini(this)

          implicit none

          class(bf_interface_icr), intent(inout) :: this

          integer(ikind) :: i, dbf_distance

          dbf_distance = bc_size

          allocate(N_detectors_list(2,nx-2*(bc_size+dbf_distance)+2))
          allocate(S_detectors_list(2,nx-2*(bc_size+dbf_distance)+2))
          allocate(E_detectors_list(2,ny-2*(bc_size+dbf_distance)))
          allocate(W_detectors_list(2,nx-2*(bc_size+dbf_distance)))

          do i=bc_size+dbf_distance, nx-(bc_size+dbf_distance)+1
             S_detectors_list(1,i-(bc_size+dbf_distance)+1) = i
             S_detectors_list(2,i-(bc_size+dbf_distance)+1) = bc_size+dbf_distance
          end do

          do i=bc_size+dbf_distance, nx-(bc_size+dbf_distance)+1
             N_detectors_list(1,i-(bc_size+dbf_distance)+1) = i
             N_detectors_list(2,i-(bc_size+dbf_distance)+1) = ny-(bc_size+dbf_distance)+1
          end do

          do i=bc_size+dbf_distance+1, ny-(bc_size+dbf_distance)
             W_detectors_list(1,i-(bc_size+dbf_distance)) = bc_size+dbf_distance
             W_detectors_list(2,i-(bc_size+dbf_distance)) = i
          end do

          do i=bc_size+dbf_distance+1, ny-(bc_size+dbf_distance)
             E_detectors_list(1,i-(bc_size+dbf_distance)) = nx-(bc_size+dbf_distance)+1
             E_detectors_list(2,i-(bc_size+dbf_distance)) = i
          end do          

        end subroutine ini


        function is_i_detector_activated()
        
          implicit none
          
          is_i_detector_activated = .true.

        end function is_i_detector_activated


        !< from a detector position, get the bc_interior_pt activated
        !> d_coords     : detector general coordinates
        !> cpt_coords_p : previous grid point whose neighboring grid points
        !>                were checked to find bc_interior_pt
        !> nb_mgrdpts   : number of grid points to be modified
        !> mgrdpts      : table containing the coordinates of the grid
        !>                points to be modified
        !> ndt_list     : temporary new detector list
        subroutine get_modified_grdpts_list(
     $     this, d_coords,
     $     nb_mgrdpts, mgrdpts, ndt_list)

          implicit none

          class(bf_interface_icr)       , intent(inout) :: this
          integer(ikind), dimension(2)  , intent(in)    :: d_coords
          integer(ikind), dimension(2)  , intent(in)    :: cpt_coords_p
          integer                       , intent(out)   :: nb_mgrdpts
          integer(ikind), dimension(2,9), intent(out)   :: mgrdpts
          type(bf_detector_i_list)      , intent(inout) :: ndt_list


          real(rkind), dimension(ne)   :: node_var
          type(bf_sublayer), pointer   :: sublayer
          integer(ikind), dimension(2) :: l_coords
          real(rkind)   , dimension(2) :: velocity
          integer(ikind), dimension(2) :: cpt_coords
          integer(ikind), dimension(2) :: d_coords_n


          !if the detector is activated, then we check
          !whether grid points need to be modified
          if(is_i_detector_activated()) then

             !extract the velocity at the coordinates of the detector
             node_var = this%get_nodes(d_coords)
             velocity(1) = node_var(2)/node_var(1)
             velocity(2) = node_var(3)/node_var(1)
             

             !get the first point from which we should look for a
             !bc_interior_pt to be activated and teh new coordinates
             !from the detector
             cpt_coords = get_central_pt(d_coords, velocity, d_coords_n)
             
             !add the new coordinates of the detector of the ndt_list
             call ndt_list%add_new_detector(d_coords_n)

             !look for a bc_interior_pt around the point previously
             !computed whose coordinates are: cpt_coords
             !we make use of the previously checked neighboring points
             !whose center was cpt_coords_p
             call this%check_neighboring_bc_interior_pts(
     $            cpt_coords_p,
     $            cpt_coords,
     $            nb_mgrdpts,
     $            mgrdpts)

          !otherwise, the coordinates of the new detector are simply
          !the previous ones, and are saved in the ndt_list
          else
             call ndt_list%add_new_detector(d_coords)
             nb_mgrdpts = 0
          end if

        end function get_modified_grdpts_list


        function get_central_grdpt(d_coords, velocity, d_coords_n)
     $     result(cpt_coords)

          implicit none

          integer(ikind), dimension(2), intent(in)  :: d_coords
          real(rkind)   , dimension(2), intent(in)  :: velocity
          integer(ikind), dimension(2), intent(out) :: d_coords_n
          integer(ikind), dimension(2)              :: cpt_coords

          real(rkind) :: dir_x, dir_y

          !1) get the direction to look for a bc_interior_pt
          dir_x  = velocity(1)*search_nb_dt*dt/dx
          dir_y  = velocity(2)*search_nb_dt*dt/dy
          
          !2) get the point indices in the direction given
          !   by the velocity vector
          cpt_coords(1) = d_coords(1) + nint(dir_x)
          cpt_coords(2) = d_coords(2) + nint(dir_y)
          
          !3) compute the new detector position
          d_coords_n(1) = d_coords(1) + nint(dir_x/SQRT(dir_x**2+dir_y**2))
          d_coords_n(2) = d_coords(2) + nint(dir_y/SQRT(dir_x**2+dir_y**2))

        end function get_central_grdpt


        subroutine check_neighboring_bc_interior_pts(
     $     this,
     $     cpt_coords_p, cpt_coords,
     $     nb_mgrdpts, mgrdpts)

          implicit none

          class(bf_interface)         , intent(in)  :: this
          integer(ikind), dimension(2), intent(in)  :: cpt_coords_p
          integer(ikind), dimension(2), intent(in)  :: cpt_coords
          integer                     , intent(out) :: nb_mgrdpts
          integer, dimension(:,:)     , intent(out) :: mgrdpts


          !to be implemented
          stop 'to be implemented: l196, bc_interface_icr_class.f'

        end subroutine check_neighboring_bc_interior_pts        


        function create_nbc_interior_pt_template(cpt_coords)
     $     result(nbc_template)

          implicit none

          integer(ikind), dimension(2)  , intent(in) :: cpt_coords
          integer       , dimension(3,3)             :: nbc_template


          integer :: i_lim1, i_lim2, i_lim4, i_lim5
          integer :: j_lim1, j_lim2, j_lim4, j_lim5
          

          i_lim1 = bc_size
          i_lim2 = bc_size+1
          i_lim4 = nx-bc_size
          i_lim5 = nx-bc_size+1

          j_lim1 = bc_size
          j_lim2 = bc_size+1
          j_lim4 = ny-bc_size
          j_lim5 = ny-bc_size+1

          
          select case(cpt_coords(1))

            case(i_lim1)
               select case(cpt_coords(2))
                 case(j_lim1)
                    nbc_template = make_nbc_template_11()
                 case(j_lim2)
                    nbc_template = make_nbc_template_12()
                 case(j_lim4)
                    nbc_template = make_nbc_template_14()
                 case(j_lim5)
                    nbc_template = make_nbc_template_15()
                 case default
                    nbc_template = make_nbc_template_13()
               end select

            case(i_lim2)
               select case(cpt_coords(2))
                 case(j_lim1)
                    nbc_template = make_nbc_template_21()
                 case(j_lim2)
                    nbc_template = make_nbc_template_22()
                 case(j_lim4)
                    nbc_template = make_nbc_template_24()
                 case(j_lim5)
                    nbc_template = make_nbc_template_25()
                 case default
                    nbc_template = make_nbc_template_23()
               end select

            case(i_lim4)
               select case(cpt_coords(2))
                 case(j_lim1)
                    nbc_template = make_nbc_template_41()
                 case(j_lim2)
                    nbc_template = make_nbc_template_42()
                 case(j_lim4)
                    nbc_template = make_nbc_template_44()
                 case(j_lim5)
                    nbc_template = make_nbc_template_45()
                 case default
                    nbc_template = make_nbc_template_43()
               end select

            case(i_lim5)
               select case(cpt_coords(2))
                 case(j_lim1)
                    nbc_template = make_nbc_template_51()
                 case(j_lim2)
                    nbc_template = make_nbc_template_52()
                 case(j_lim4)
                    nbc_template = make_nbc_template_54()
                 case(j_lim5)
                    nbc_template = make_nbc_template_55()
                 case default
                    nbc_template = make_nbc_template_53()
               end select

            case default
               select case(cpt_coords(2))
                 case(j_lim1)
                    nbc_template = make_nbc_template_31()
                 case(j_lim2)
                    nbc_template = make_nbc_template_32()
                 case(j_lim4)
                    nbc_template = make_nbc_template_34()
                 case(j_lim5)
                    nbc_template = make_nbc_template_35()
               end select
          end select

        end function create_nbc_interior_pt_template

      end module bf_interface_icr_class

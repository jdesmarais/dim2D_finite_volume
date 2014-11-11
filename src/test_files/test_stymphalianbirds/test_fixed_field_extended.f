      !> @file
      !> simulation test of the wave2d governing equations
      !> on a domain with buffer layers
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> run a simulation using the physical model
      !> (in makefile_header choices) on a computational
      !> domain with buffer layers
      !-----------------------------------------------------------------
      program test_fixed_field_extended

        use bf_interface_class, only :
     $       bf_interface

        use bf_mainlayer_class, only :
     $       bf_mainlayer

        use bf_sublayer_class, only :
     $       bf_sublayer

        use field_extended_class, only :
     $       field_extended

        use parameters_bf_layer, only :
     $       align_N,
     $       align_S,
     $       align_E,
     $       align_W,
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt,
     $       no_pt

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       t_max,
     $       dt,
     $       detail_print

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none


        integer :: bf_config

        
        !configuration of the domain extended:
        !
        ! 1 : simple extension:
        !          _________
        !         |_________|
        !         |  |   |  |
        !         |__|___|__|
        !         |_________|
        !
        ! 2 : cropped buffer layers:
        !               ____
        !             _|____|
        !          __|   |__|
        !         |__|___|_
        !         |    ____|
        !         |____|
        !
        bf_config = 2

        call make_test_fixed_field_extended(bf_config)

        
        contains

        subroutine make_test_fixed_field_extended(bf_config)

          implicit none

          ! configuration for the domain
          ! extension
          ! 1: simple extension
          ! 2: cropped buffer layers
          integer, intent(in) :: bf_config


          ! size of the buffer layer added to the interior domain
          integer(ikind) :: size_bf

          ! computational domain simulated
          type(field_extended) :: f_simulated
        
          ! intermediate variables for the
          ! simulation
          integer(ikind) :: nt, output_print
          integer(ikind) :: t
  
          ! CPU recorded times
          real :: time1, time2, time3  
          
          !choice of the buffer layer size
          if(nx.le.20) then
             size_bf = 10
          else
             size_bf = nint(nx/3.0)
          end if

  
          ! get the initial CPU time
          call CPU_TIME(time1)
  
  
          ! time parameters
          nt           = int(t_max/dt)
          output_print = int(1.0d0/detail_print)
  
  
          ! initialize the field
          call f_simulated%ini()
          call f_simulated%apply_bc_on_nodes()
  
          call add_buffer_layers_for_test(f_simulated,bf_config,size_bf)
  
          call f_simulated%write_data()
  
          ! initialization time
          call CPU_TIME(time2)
          print *, 'time_elapsed: ', time2-time1
  
  
          ! integrate the field until t=t_max
          do t=1, nt
  
             !DEC$ FORCEINLINE RECURSIVE
             call f_simulated%integrate(dt)
  
             !  write the output data
             if((output_print.eq.1).or.
     $            ((output_print.ne.0).and.(mod(t,output_print).eq.0))) then
                call f_simulated%write_data()
             end if
  
          end do
  
  
          ! print the time needed for the simulation
          call CPU_TIME(time3)
          print *, 'time_elapsed: ', time3-time1
  
  
          ! write the last timestep
          if((output_print.eq.0).or.(mod(nt,output_print).ne.0)) then
             call f_simulated%write_data()
          end if

        end subroutine make_test_fixed_field_extended


        subroutine add_buffer_layers_for_test(
     $       field_ext_used,
     $       config,
     $       size_bf)

          implicit none

          class(field_extended), intent(inout) :: field_ext_used
          integer              , intent(in)    :: config
          integer              , intent(in)    :: size_bf

          select case(config)
            case(1)
               call four_bf_layer_config(field_ext_used,size_bf)

            case(2)
               call difficult_bf_layer_config(field_ext_used,size_bf)

            case default
               print '(''test_field_extended_wave2d'')'
               print '(''add_buffer_layers_for_test'')'
               stop 'config not recognized'

          end select

        end subroutine add_buffer_layers_for_test


        subroutine four_bf_layer_config(field_ext_tested,size_bf)

          implicit none

          class(field_extended), intent(inout) :: field_ext_tested
          integer(ikind)       , intent(in)    :: size_bf

          real(rkind)   , dimension(nx)       :: interior_x_map
          real(rkind)   , dimension(ny)       :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
          integer(ikind), dimension(2,2)      :: alignment
          
          interior_x_map = field_ext_tested%get_x_map()
          interior_y_map = field_ext_tested%get_y_map()
          interior_nodes = field_ext_tested%get_nodes()


          !extend the interior domain to have
          !four buffer layers
          !----------------------------------
          !add the north buffer layer
          alignment(1,1) = align_W-(size_bf-1)
          alignment(1,2) = align_E+(size_bf-1)
          alignment(2,1) = align_N
          alignment(2,2) = align_N+(size_bf-1)
          call add_sublayer(
     $         field_ext_tested,
     $         N,
     $         alignment,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         size_bf)
          
          !add the south buffer layer
          alignment(1,1) = align_W-(size_bf-1)
          alignment(1,2) = align_E+(size_bf-1)
          alignment(2,1) = align_S-(size_bf-1)
          alignment(2,2) = align_S
          call add_sublayer(
     $         field_ext_tested,
     $         S,
     $         alignment,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         size_bf)

          !add the east buffer layer
          alignment(1,1) = align_E
          alignment(1,2) = align_E+(size_bf-1)
          alignment(2,1) = align_S+1
          alignment(2,2) = align_N-1
          call add_sublayer(
     $         field_ext_tested,
     $         E,
     $         alignment,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         size_bf)

          !add the west buffer layer
          alignment(1,1) = align_W-(size_bf-1)
          alignment(1,2) = align_W
          alignment(2,1) = align_S+1
          alignment(2,2) = align_N-1
          call add_sublayer(
     $         field_ext_tested,
     $         W,
     $         alignment,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         size_bf)

        end subroutine four_bf_layer_config

      
        subroutine difficult_bf_layer_config(field_ext_tested,size_bf)

          implicit none

          class(field_extended), intent(inout) :: field_ext_tested
          integer(ikind)       , intent(in)    :: size_bf

          real(rkind)   , dimension(nx)       :: interior_x_map
          real(rkind)   , dimension(ny)       :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
          integer(ikind), dimension(2,2)      :: alignment
          
          interior_x_map = field_ext_tested%get_x_map()
          interior_y_map = field_ext_tested%get_y_map()
          interior_nodes = field_ext_tested%get_nodes()


          !extend the interior domain to have
          !four buffer layers
          !----------------------------------
          !add the north buffer layer
          alignment(1,1) = align_E-(size_bf-1)
          alignment(1,2) = align_E+(size_bf-1)
          alignment(2,1) = align_N
          alignment(2,2) = align_N+(size_bf-1)
          call add_sublayer(
     $         field_ext_tested,
     $         N,
     $         alignment,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         size_bf)
          
          !add the south buffer layer
          alignment(1,1) = align_W-(size_bf-1)
          alignment(1,2) = align_E+(size_bf-1)
          alignment(2,1) = align_S-(size_bf-1)
          alignment(2,2) = align_S
          call add_sublayer(
     $         field_ext_tested,
     $         S,
     $         alignment,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         size_bf,
     $         difficult_geometry=.true.)

          !add the east buffer layer
          alignment(1,1) = align_E
          alignment(1,2) = align_E+(size_bf-1)
          alignment(2,1) = align_N-size_bf
          alignment(2,2) = align_N-1
          call add_sublayer(
     $         field_ext_tested,
     $         E,
     $         alignment,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         size_bf,
     $         difficult_geometry=.true.)

          !add the west buffer layer
          alignment(1,1) = align_W-(size_bf-1)
          alignment(1,2) = align_W
          alignment(2,1) = align_S+1
          alignment(2,2) = align_S+size_bf
          call add_sublayer(
     $         field_ext_tested,
     $         W,
     $         alignment,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         size_bf,
     $         difficult_geometry=.true.)

        end subroutine difficult_bf_layer_config


        subroutine add_sublayer(
     $     field_ext_tested,
     $     mainlayer_id,
     $     alignment,
     $     x_map, y_map, nodes,
     $     size_bf,
     $     difficult_geometry)

          implicit none

          type(field_extended)            , intent(inout) :: field_ext_tested
          integer                         , intent(in)    :: mainlayer_id
          integer(ikind), dimension(4)    , intent(inout) :: alignment
          real(rkind)   , dimension(:)    , intent(in)    :: x_map
          real(rkind)   , dimension(:)    , intent(in)    :: y_map
          real(rkind)   , dimension(:,:,:), intent(in)    :: nodes
          integer(ikind)                  , intent(in)    :: size_bf
          logical       , optional        , intent(in)    :: difficult_geometry

          type(bf_sublayer), pointer :: added_sublayer

          real(rkind)   , dimension(:)    , allocatable :: bf_x_map
          real(rkind)   , dimension(:)    , allocatable :: bf_y_map
          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes
          integer(ikind), dimension(:,:)  , allocatable :: bf_grdpts_id

          integer(ikind) :: i,j

          logical :: difficult_geometry_op

          !check whether the geometry simulated is
          !the conmplex case
          if(present(difficult_geometry)) then
             difficult_geometry_op = difficult_geometry
          else
             difficult_geometry_op = .false.
          end if

          !get x_map and y_map
          added_sublayer => field_ext_tested%domain_extension%allocate_sublayer(
     $         mainlayer_id,
     $         x_map,
     $         y_map,
     $         nodes,
     $         alignment)

          call added_sublayer%get_x_map(bf_x_map)
          call added_sublayer%get_y_map(bf_y_map)

          !set nodes
          allocate(bf_nodes(size(bf_x_map,1),size(bf_y_map,1),ne))
          call initialize_nodes(
     $         field_ext_tested,
     $         bf_x_map,
     $         bf_y_map,
     $         bf_nodes)

          call added_sublayer%set_nodes(bf_nodes)

          !set grdpts_id
          allocate(bf_grdpts_id(
     $         size(bf_x_map,1),
     $         size(bf_y_map,1)))

          deallocate(bf_x_map)
          deallocate(bf_y_map)
          
          select case(mainlayer_id)
            case(N)
               
               !----------------------------------
               do j=1, size(bf_grdpts_id,2)-bc_size

                  i=1
                  bf_grdpts_id(i,j) = bc_pt

                  i=bc_size
                  bf_grdpts_id(i,j) = bc_interior_pt

                  do i=bc_size+1, size(bf_grdpts_id,1)-bc_size
                     bf_grdpts_id(i,j) = interior_pt
                  end do

                  i=size(bf_grdpts_id,1)-bc_size+1
                  bf_grdpts_id(i,j) = bc_interior_pt

                  i=size(bf_grdpts_id,1)
                  bf_grdpts_id(i,j) = bc_pt

               end do

               j=size(bf_grdpts_id,2)-1
               !----------------------------------
               i=1
               bf_grdpts_id(i,j)=bc_pt

               do i=bc_size,size(bf_grdpts_id,1)-1
                  bf_grdpts_id(i,j) = bc_interior_pt
               end do

               i=size(bf_grdpts_id,1)
               bf_grdpts_id(i,j) = bc_pt

               j=size(bf_grdpts_id,2)
               !----------------------------------
               do i=1, size(bf_grdpts_id,1)
                  bf_grdpts_id(i,j) = bc_pt
               end do
               
            case(S)

               if(difficult_geometry_op) then

                  j=1
                  !----------------------------------
                  do i=1,2*size_bf
                     bf_grdpts_id(i,j) = bc_pt
                  end do

                  do i=2*size_bf+1,size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = no_pt
                  end do

                  
                  j=2
                  !----------------------------------
                  i=1
                  bf_grdpts_id(i,j) = bc_pt
                     
                  do i=2, 2*size_bf-1
                     bf_grdpts_id(i,j) = bc_interior_pt
                  end do
                  
                  i=2*size_bf
                  bf_grdpts_id(i,j) = bc_pt
                     
                  do i=2*size_bf+1, size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = no_pt
                  end do
                  

                  do j=3,4
                  !----------------------------------

                     i=1
                     bf_grdpts_id(i,j) = bc_pt
                     
                     i=2
                     bf_grdpts_id(i,j) = bc_interior_pt

                     do i=3, 2*size_bf-bc_size
                        bf_grdpts_id(i,j) = interior_pt
                     end do
                     
                     i=2*size_bf-1
                     bf_grdpts_id(i,j) = bc_interior_pt

                     i=2*size_bf
                     bf_grdpts_id(i,j) = bc_pt
                     
                     do i=2*size_bf+1, size(bf_grdpts_id,1)
                        bf_grdpts_id(i,j) = no_pt
                     end do

                  end do

                  
                  j=5
                  !----------------------------------
                  i=1
                  bf_grdpts_id(i,j) = bc_pt
                  
                  i=2
                  bf_grdpts_id(i,j) = bc_interior_pt

                  do i=3, 2*size_bf-bc_size
                     bf_grdpts_id(i,j) = interior_pt
                  end do
                  
                  i=2*size_bf-1
                  bf_grdpts_id(i,j) = bc_interior_pt

                  do i=2*size_bf, size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = bc_pt
                  end do

                  
                  j=6
                  !----------------------------------
                  i=1
                  bf_grdpts_id(i,j) = bc_pt
                  
                  i=2
                  bf_grdpts_id(i,j) = bc_interior_pt

                  do i=3, 2*size_bf-bc_size
                     bf_grdpts_id(i,j) = interior_pt
                  end do

                  do i=2*size_bf-bc_size+1, size(bf_grdpts_id,1)-1
                     bf_grdpts_id(i,j) = bc_interior_pt
                  end do

                  i=size(bf_grdpts_id,1)
                  bf_grdpts_id(i,j) = bc_pt


                  do j=7, size(bf_grdpts_id,2)
                  !----------------------------------

                     i=1
                     bf_grdpts_id(i,j) = bc_pt
                     
                     i=2
                     bf_grdpts_id(i,j) = bc_interior_pt

                     do i=3, size(bf_grdpts_id,1)-2
                        bf_grdpts_id(i,j) = interior_pt
                     end do

                     i=size(bf_grdpts_id,1)-1
                     bf_grdpts_id(i,j) = bc_interior_pt

                     i=size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = bc_pt

                  end do

               else

                  j=1
                  !----------------------------------
                  do i=1, size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = bc_pt
                  end do

                  j=bc_size
                  !----------------------------------
                  i=1
                  bf_grdpts_id(i,j)=bc_pt

                  do i=bc_size,size(bf_grdpts_id,1)-1
                     bf_grdpts_id(i,j) = bc_interior_pt
                  end do

                  i=size(bf_grdpts_id,1)
                  bf_grdpts_id(i,j) = bc_pt
                  !----------------------------------
                  do j=bc_size+1, size(bf_grdpts_id,2)

                     i=1
                     bf_grdpts_id(i,j) = bc_pt

                     i=bc_size
                     bf_grdpts_id(i,j) = bc_interior_pt

                     do i=bc_size+1, size(bf_grdpts_id,1)-bc_size
                        bf_grdpts_id(i,j) = interior_pt
                     end do

                     i=size(bf_grdpts_id,1)-bc_size+1
                     bf_grdpts_id(i,j) = bc_interior_pt

                     i=size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = bc_pt

                  end do

               end if

            case(E)
               
               if(difficult_geometry_op) then

                  j=1
                  !----------------------------------
                  do i=1, size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = bc_pt
                  end do

                  j=2
                  !----------------------------------
                  do i=1, size(bf_grdpts_id,1)-1
                     bf_grdpts_id(i,j) = bc_interior_pt
                  end do

                  i=size(bf_grdpts_id,1)
                  bf_grdpts_id(i,j) = bc_pt

                  
                  do j=3, size(bf_grdpts_id,2)
                  !----------------------------------
                     do i=1, size(bf_grdpts_id,1)-bc_size
                        bf_grdpts_id(i,j) = interior_pt
                     end do
                     
                     i=size(bf_grdpts_id,1)-1
                     bf_grdpts_id(i,j) = bc_interior_pt

                     i=size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = bc_pt

                  end do

               else

                  do j=1, size(bf_grdpts_id,2)
                     do i=1, size(bf_grdpts_id,1) -bc_size
                        bf_grdpts_id(i,j) = interior_pt
                     end do
                     
                     i=size(bf_grdpts_id,1)-1
                     bf_grdpts_id(i,j) = bc_interior_pt

                     i=size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = bc_pt

                  end do

               end if

            case(W)

               if(difficult_geometry_op) then

                  do j=1, size(bf_grdpts_id,2)-bc_size
                  !----------------------------------

                     i=1
                     bf_grdpts_id(i,j) = bc_pt
                     
                     i=bc_size
                     bf_grdpts_id(i,j) = bc_interior_pt

                     do i=bc_size+1, size(bf_grdpts_id,1)
                        bf_grdpts_id(i,j) = interior_pt
                     end do

                  end do


                  j=size(bf_grdpts_id,2)-1
                  !----------------------------------
                  i=1
                  bf_grdpts_id(i,j) = bc_pt
                  
                  do i=bc_size, size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = bc_interior_pt
                  end do

                  
                  j=size(bf_grdpts_id,2)
                  !----------------------------------
                  do i=1, size(bf_grdpts_id,1)
                     bf_grdpts_id(i,j) = bc_pt
                  end do                  

               else

                  do j=1, size(bf_grdpts_id,2)

                     i=1
                     bf_grdpts_id(i,j) = bc_pt
                     
                     i=bc_size
                     bf_grdpts_id(i,j) = bc_interior_pt

                     do i=bc_size+1, size(bf_grdpts_id,1)
                        bf_grdpts_id(i,j) = interior_pt
                     end do

                  end do

               end if

          end select

          call added_sublayer%set_grdpts_id(bf_grdpts_id)
          
          !update integration borders
          call field_ext_tested%domain_extension%update_integration_borders(added_sublayer)

        end subroutine add_sublayer


        subroutine initialize_nodes(
     $     field_ext_used,
     $     x_map,
     $     y_map,
     $     nodes)

          implicit none

          class(field_extended)        , intent(in)  :: field_ext_used
          real(rkind), dimension(:)    , intent(in)  :: x_map
          real(rkind), dimension(:)    , intent(in)  :: y_map
          real(rkind), dimension(:,:,:), intent(out) :: nodes

          call field_ext_used%pmodel_eq_used%apply_ic(
     $         nodes,
     $         x_map,
     $         y_map)

        end subroutine initialize_nodes

      end program test_fixed_field_extended

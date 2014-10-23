      !> @file
      !> test file for the object 'bf_compute'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the computation of the time derivatives
      !> using the finite volume method by comparing 
      !> with expected data
      !
      !> @date
      ! 16_07_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_bf_compute_prog

        use bc_operators_class , only : bc_operators

        use bf_compute_class   , only : bf_compute

        use parameters_bf_layer, only : interior_pt

        use parameters_constant, only : periodic_xy_choice,
     $                                  bc_nodes_choice

        use parameters_input   , only : ne,bc_choice,bc_size,
     $                                  bc_N_type_choice,
     $                                  bc_S_type_choice,
     $                                  bc_E_type_choice,
     $                                  bc_W_type_choice

        use parameters_kind    , only : ikind, rkind

        use pmodel_eq_class    , only : pmodel_eq

        use sd_operators_class , only : sd_operators

        use td_operators_class , only : td_operators

        implicit none

        
        !<operators tested
        integer    , parameter             :: nxt = 10
        integer    , parameter             :: nyt = 6
        real(rkind), dimension(nxt,nyt,ne) :: nodes
        integer    , dimension(nxt,nyt)    :: grdpts_id
        real(rkind)                        :: dx
        real(rkind)                        :: dy
        type(bf_compute)                   :: bf_compute_used
        type(sd_operators)                 :: sd_operators_used
        type(pmodel_eq)                    :: pmodel_eq_used
        type(bc_operators)                 :: bc_operators_used
        type(td_operators)                 :: td_operators_used

        real(rkind), dimension(:,:,:), allocatable :: time_dev


        !<CPU recorded times
        real :: time1, time2

        !<test parameters
        logical, parameter         :: detailled=.true.
        integer(ikind)             :: i,j
        logical                    :: test_global,test_local
        logical                    :: test_parameter


        !<if nx<4, ny<4 then the test cannot be done
        test_parameter=.true.
        test_parameter=test_parameter.and.(nxt.eq.10)
        test_parameter=test_parameter.and.(nyt.eq.6)
        test_parameter=test_parameter.and.(ne.eq.1)
        test_parameter=test_parameter.and.(bc_choice.eq.periodic_xy_choice)
        test_parameter=test_parameter.and.(bc_N_type_choice.eq.bc_nodes_choice)
        test_parameter=test_parameter.and.(bc_S_type_choice.eq.bc_nodes_choice)    
        test_parameter=test_parameter.and.(bc_E_type_choice.eq.bc_nodes_choice)    
        test_parameter=test_parameter.and.(bc_W_type_choice.eq.bc_nodes_choice)    
        if(.not.test_parameter) then
           print *, 'the test requires several parameters'
           print *, 'test designed for simpletest eq'
           print *, 'nx=10'
           print *, 'ny=6'
           print *, 'ne=1'
           print *, 'bc_choice=periodic_xy_choice'
           print *, 'bc_N_type_choice=bc_nodes_choice'
           print *, 'bc_S_type_choice=bc_nodes_choice'
           print *, 'bc_E_type_choice=bc_nodes_choice'
           print *, 'bc_W_type_choice=bc_nodes_choice'
           stop ''
        end if


        !<get the initial CPU time
        call CPU_TIME(time1)


        !<initialize the tables for the field
        dx=1.0
        dy=1.0

        do j=1, nyt
           do i=1, nxt
              nodes(i,j,1) = i + (j-1)*nxt
              grdpts_id(i,j) = interior_pt
           end do
        end do


        !<compute the time derivatives
        call bf_compute_used%allocate_tables(nxt,nyt)
        call bf_compute_used%compute_time_dev(
     $       nodes,
     $       dx,
     $       dy,
     $       sd_operators_used,
     $       pmodel_eq_used,
     $       bc_operators_used,
     $       td_operators_used,
     $       grdpts_id)
        call bf_compute_used%get_time_dev(time_dev)
        call bf_compute_used%deallocate_tables()


        !<print the time derivatives
        test_global=.true.
        do j=1+bc_size, nyt-bc_size
           do i=1+bc_size, nxt-bc_size
              test_local=(time_dev(i,j,1).eq.(-20.0d0))
              test_global=test_global.and.test_local
              if(detailled) then
                 print *, i,j, test_local
              end if
           end do
        end do
        print '(''test_validated: '', L1)', test_global

        call CPU_TIME(time2)
        print '(''time elapsed:'', F6.2)', time2-time1

      end program test_bf_compute_prog

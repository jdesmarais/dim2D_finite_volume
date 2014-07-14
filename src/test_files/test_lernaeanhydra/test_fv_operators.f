      !> @file
      !> test file for the object 'fv_operators'
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
      ! 13_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_fv_operators

        use bc_operators_class , only : bc_operators
        use sd_operators_class , only : sd_operators
        use td_operators_class , only : td_operators
        use parameters_constant, only : periodic_xy_choice,
     $                                  bc_nodes_choice
        use parameters_input   , only : nx,ny,ne,bc_size,bc_choice,
     $                                  bcx_type_choice,bcy_type_choice
        use parameters_kind    , only : ikind, rkind
        use pmodel_eq_class    , only : pmodel_eq

        implicit none

        
        !<operators tested
        real(rkind), dimension(nx,ny,ne) :: nodes
        real(rkind)                      :: dx
        real(rkind)                      :: dy
        type(sd_operators)               :: s
        type(pmodel_eq)                  :: p_model
        type(bc_operators)               :: bc_used
        type(td_operators)               :: t_operator

        real(rkind), dimension(nx,ny,ne) :: time_dev

        !<CPU recorded times
        real :: time1, time2

        !<test parameters
        logical, parameter         :: detailled=.true.
        integer(ikind)             :: i,j
        real(rkind), dimension(12) :: test_data
        logical                    :: global,local
        logical                    :: test_parameter


        !<if nx<4, ny<4 then the test cannot be done
        test_parameter=.true.
        test_parameter=test_parameter.and.(nx.eq.10)
        test_parameter=test_parameter.and.(ny.eq.6)
        test_parameter=test_parameter.and.(ne.eq.1)
        test_parameter=test_parameter.and.(bc_choice.eq.periodic_xy_choice)
        test_parameter=test_parameter.and.(bcx_type_choice.eq.bc_nodes_choice)
        test_parameter=test_parameter.and.(bcy_type_choice.eq.bc_nodes_choice)    
        if(.not.test_parameter) then
           print *, 'the test requires several parameters'
           print *, 'test designed for simpletest eq'
           print *, 'nx=10'
           print *, 'ny=6'
           print *, 'ne=1'
           print *, 'bc_choice=periodic_xy_choice'
           print *, 'bcx_type_choice=bc_nodes_choice'
           print *, 'bcy_type_choice=bc_nodes_choice'
           stop ''
        end if


        !<get the initial CPU time
        call CPU_TIME(time1)


        !<initialize the tables for the field
        dx=1.0
        dy=1.0

        do j=1, ny
           do i=1, nx
              nodes(i,j,1) = i + (j-1)*nx
           end do
        end do


        !<compute the time derivatives
        time_dev = t_operator%compute_time_dev(
     $       nodes,dx,dy,s,p_model,bc_used)


        !<print the time derivatives
        global=.true.
        do j=1+bc_size, ny-bc_size
           do i=1+bc_size, nx-bc_size
              local=(time_dev(i,j,1).eq.(-20.0d0))
              global=global.and.local
              if(detailled) then
                 print *, i,j, local
              end if
           end do
        end do
        print '(''test_validated: '', 1L)', global

        call CPU_TIME(time2)
        print '(''time elapsed:'', F6.2)', time2-time1

      end program test_fv_operators

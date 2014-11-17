      !> @file
      !> test of the computation of gradients in n1- and n2- directions
      !> using Mattsson's operators
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test of the computation of gradients in n1- and n2- directions
      !> using Mattsson's operators
      !
      !> @date
      !> 17_11_2014 - Initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_sd_operators_fd_n_module

        use interface_primary, only :
     $       gradient_n_proc

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use sd_operators_fd_n_module, only :
     $       gradient_n1_xL0_yL0,
     $       gradient_n1_xL0_yL1,
     $       gradient_n1_xL0_yI,
     $       gradient_n1_xL0_yR1,
     $       gradient_n1_xL0_yR0,
     $       
     $       gradient_n1_xL1_yL0,
     $       gradient_n1_xL1_yL1,
     $       gradient_n1_xL1_yI,
     $       gradient_n1_xL1_yR1,
     $       gradient_n1_xL1_yR0,
     $       
     $       gradient_n1_xI_yL0,
     $       gradient_n1_xI_yL1,
     $       gradient_n1_xI_yI,
     $       gradient_n1_xI_yR1,
     $       gradient_n1_xI_yR0,
     $       
     $       gradient_n1_xR1_yL0,
     $       gradient_n1_xR1_yL1,
     $       gradient_n1_xR1_yI,
     $       gradient_n1_xR1_yR1,
     $       gradient_n1_xR1_yR0,
     $       
     $       gradient_n1_xR0_yL0,
     $       gradient_n1_xR0_yL1,
     $       gradient_n1_xR0_yI,
     $       gradient_n1_xR0_yR1,
     $       gradient_n1_xR0_yR0,
     $       
     $       gradient_n2_xL0_yL0,
     $       gradient_n2_xL0_yL1,
     $       gradient_n2_xL0_yI,
     $       gradient_n2_xL0_yR1,
     $       gradient_n2_xL0_yR0,
     $       
     $       gradient_n2_xL1_yL0,
     $       gradient_n2_xL1_yL1,
     $       gradient_n2_xL1_yI,
     $       gradient_n2_xL1_yR1,
     $       gradient_n2_xL1_yR0,
     $       
     $       gradient_n2_xI_yL0,
     $       gradient_n2_xI_yL1,
     $       gradient_n2_xI_yI,
     $       gradient_n2_xI_yR1,
     $       gradient_n2_xI_yR0,
     $       
     $       gradient_n2_xR1_yL0,
     $       gradient_n2_xR1_yL1,
     $       gradient_n2_xR1_yI,
     $       gradient_n2_xR1_yR1,
     $       gradient_n2_xR1_yR0,
     $       
     $       gradient_n2_xR0_yL0,
     $       gradient_n2_xR0_yL1,
     $       gradient_n2_xR0_yI,
     $       gradient_n2_xR0_yR1,
     $       gradient_n2_xR0_yR0


        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        real(rkind), dimension(5,5,1) :: nodes
        real(rkind)                   :: dx
        real(rkind)                   :: dy

        nodes = reshape((/
     $       0.1d0, -2.6d0, 6.86d0, 9.48d0, 7.563d0,
     $       14.3d0, 8.23d0, 2.13d0, 2.536d0, 8.12d0,
     $       -9.26d0, 6.15d0, 7.5d0, -9.15d0, 6.48d0,
     $       5.23d0, 7.123d0, 4.19d0, 9.26d0, -3.25d0,
     $       1.02d0, 6.52d0, 8.15d0, -2.15d0, -6.48d0/),
     $       (/5,5,1/))

        test_validated = .true.
        detailled      = .false.

        dx = 0.1d0
        dy = 0.2d0


        !n1_xL0_yL0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xL0_yL0,
     $       -106.03066d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xL0_yL0: '',L1)', test_loc


        !n1_xL0_yL1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xL0_yL1,
     $       -121.374879d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xL0_yL1: '',L1)', test_loc


        !n1_xL0_yI
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xL0_yI,
     $       -121.374879d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xL0_yI: '',L1)', test_loc


        !n1_xL0_yR1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xL0_yR1,
     $       -121.374879d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xL0_yR1: '',L1)', test_loc


        !n1_xL0_yR0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xL0_yR0,
     $       -136.7190961d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xL0_yR0: '',L1)', test_loc


        !n1_xL1_yL0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xL1_yL0,
     $       -42.39105153d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xL1_yL0: '',L1)', test_loc


        !n1_xL1_yL1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xL1_yL1,
     $       -57.73526868d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xL1_yL1: '',L1)', test_loc


        !n1_xL1_yI
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xL1_yI,
     $       -57.73526868d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xL1_yI: '',L1)', test_loc


        !n1_xL1_yR1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xL1_yR1,
     $       -57.73526868d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xL1_yR1: '',L1)', test_loc


        !n1_xL1_yR0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xL1_yR0,
     $       -73.07948584d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xL1_yR0: '',L1)', test_loc


        !n1_xI_yL0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xI_yL0,
     $       -42.39105153d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xI_yL0: '',L1)', test_loc


        !n1_xI_yL1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xI_yL1,
     $       -57.73526868d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xI_yL1: '',L1)', test_loc


        !n1_xI_yI
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xI_yI,
     $       -57.73526868d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xI_yI: '',L1)', test_loc


        !n1_xI_yR1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xI_yR1,
     $       -57.73526868d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xI_yR1: '',L1)', test_loc


        !n1_xI_yR0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xI_yR0,
     $       -73.07948584d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xI_yR0: '',L1)', test_loc


        !n1_xR1_yL0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xR1_yL0,
     $       -42.39105153d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xR1_yL0: '',L1)', test_loc


        !n1_xR1_yL1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xR1_yL1,
     $       -57.73526868d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xR1_yL1: '',L1)', test_loc


        !n1_xR1_yI
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xR1_yI,
     $       -57.73526868d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xR1_yI: '',L1)', test_loc


        !n1_xR1_yR1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xR1_yR1,
     $       -57.73526868d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xR1_yR1: '',L1)', test_loc


        !n1_xR1_yR0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xR1_yR0,
     $       -73.07948584d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xR1_yR0: '',L1)', test_loc


        !n1_xR0_yL0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xR0_yL0,
     $       21.24855877d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xR0_yL0: '',L1)', test_loc


        !n1_xR0_yL1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xR0_yL1,
     $       5.904341623d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xR0_yL1: '',L1)', test_loc


        !n1_xR0_yI
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xR0_yI,
     $       5.904341623d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xR0_yI: '',L1)', test_loc


        !n1_xR0_yR1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xR0_yR1,
     $       5.904341623d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xR0_yR1: '',L1)', test_loc


        !n1_xR0_yR0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n1_xR0_yR0,
     $       -9.439875529d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n1_xR0_yR0: '',L1)', test_loc


        !n2_xL0_yL0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xL0_yL0,
     $       -129.4358963d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xL0_yL0: '',L1)', test_loc


        !n2_xL0_yL1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xL0_yL1,
     $       -114.0916791d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xL0_yL1: '',L1)', test_loc


        !n2_xL0_yI
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xL0_yI,
     $       -114.0916791d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xL0_yI: '',L1)', test_loc


        !n2_xL0_yR1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xL0_yR1,
     $       -114.0916791d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xL0_yR1: '',L1)', test_loc


        !n2_xL0_yR0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xL0_yR0,
     $       -98.74746199d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xL0_yR0: '',L1)', test_loc


        !n2_xL1_yL0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xL1_yL0,
     $       -65.79628599d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xL1_yL0: '',L1)', test_loc


        !n2_xL1_yL1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xL1_yL1,
     $       -50.45206884d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xL1_yL1: '',L1)', test_loc


        !n2_xL1_yI
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xL1_yI,
     $       -50.45206884d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xL1_yI: '',L1)', test_loc


        !n2_xL1_yR1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xL1_yR1,
     $       -50.45206884d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xL1_yR1: '',L1)', test_loc


        !n2_xL1_yR0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xL1_yR0,
     $       -35.10785169d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xL1_yR0: '',L1)', test_loc


        !n2_xI_yL0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xI_yL0,
     $       -65.79628599d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xI_yL0: '',L1)', test_loc


        !n2_xI_yL1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xI_yL1,
     $       -50.45206884d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xI_yL1: '',L1)', test_loc


        !n2_xI_yI
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xI_yI,
     $       -50.45206884d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xI_yI: '',L1)', test_loc


        !n2_xI_yR1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xI_yR1,
     $       -50.45206884d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xI_yR1: '',L1)', test_loc


        !n2_xI_yR0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xI_yR0,
     $       -35.10785169d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xI_yR0: '',L1)', test_loc


        !n2_xR1_yL0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xR1_yL0,
     $       -65.79628599d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xR1_yL0: '',L1)', test_loc


        !n2_xR1_yL1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xR1_yL1,
     $       -50.45206884d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xR1_yL1: '',L1)', test_loc


        !n2_xR1_yI
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xR1_yI,
     $       -50.45206884d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xR1_yI: '',L1)', test_loc


        !n2_xR1_yR1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xR1_yR1,
     $       -50.45206884d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xR1_yR1: '',L1)', test_loc


        !n2_xR1_yR0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xR1_yR0,
     $       -35.10785169d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xR1_yR0: '',L1)', test_loc


        !n2_xR0_yL0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xR0_yL0,
     $       -2.156675683d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xR0_yL0: '',L1)', test_loc


        !n2_xR0_yL1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xR0_yL1,
     $       13.18754147d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xR0_yL1: '',L1)', test_loc


        !n2_xR0_yI
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xR0_yI,
     $       13.18754147d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xR0_yI: '',L1)', test_loc


        !n2_xR0_yR1
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xR0_yR1,
     $       13.18754147d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xR0_yR1: '',L1)', test_loc


        !n2_xR0_yR0
        test_loc = test_gradient_n(
     $       nodes,dx,dy,
     $       gradient_n2_xR0_yR0,
     $       28.53175862d0,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_n2_xR0_yR0: '',L1)', test_loc


        contains

        
        function test_gradient_n(
     $       nodes,dx,dy,
     $       grad_n_proc,
     $       grad_n_data,
     $       detailled)
     $       result(test_validated)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          procedure(gradient_n_proc)                :: grad_n_proc
          real(rkind)                  , intent(in) :: grad_n_data
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated
          
          real(rkind) :: test_grad

          test_grad = grad_n_proc(nodes,3,3,simple_var_proc,dx,dy)

          test_validated = is_test_validated(
     $         test_grad, grad_n_data, detailled)

        end function test_gradient_n


        function simple_var_proc(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var = nodes(i,j,1)

        end function simple_var_proc


        function is_test_validated(var,cst,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical    , intent(in) :: detailled
          logical                 :: test_validated

          if(detailled) then
             print *, int(var*1e5)
             print *, int(cst*1e5)
          end if
          
          test_validated=abs(
     $         int(var*1e5)-
     $         int(cst*1e5)).le.1
          
        end function is_test_validated

      end program test_sd_operators_fd_n_module

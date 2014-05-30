      module test_nbf_list_module

        use bf_sublayer_class  , only : bf_sublayer
        use nbf_list_class     , only : nbf_list
        use nbf_element_class  , only : nbf_element
        use parameters_bf_layer, only : align_N, align_S,
     $                                  align_E, align_W
        use parameters_constant , only : N,S,E,W
        use parameters_input    , only : bc_size
        use parameters_kind     , only : ikind


        implicit none


        public :: ini_alignment,
     $            print_nbf_list

        contains


        subroutine ini_alignment(
     $     test_alignment)

          implicit none

          integer(ikind), dimension(4,5,2,2), intent(out) :: test_alignment
          
          integer(ikind) :: large_layer
          integer(ikind) :: small_layer

          large_layer = 5
          small_layer = 3

          !north
          test_alignment(N,1,1,1) = align_W+1
          test_alignment(N,1,1,2) = test_alignment(N,1,1,1)+large_layer
          test_alignment(N,1,2,1) = align_N
          test_alignment(N,1,2,2) = align_N

          test_alignment(N,2,1,1) = test_alignment(N,1,1,2)+2*bc_size+small_layer
          test_alignment(N,2,1,2) = test_alignment(N,2,1,1)+large_layer
          test_alignment(N,2,2,1) = align_N
          test_alignment(N,2,2,2) = align_N

          test_alignment(N,3,1,2) = align_E-1
          test_alignment(N,3,1,1) = test_alignment(N,3,1,2)-2*bc_size-large_layer
          test_alignment(N,3,2,1) = align_N
          test_alignment(N,3,2,2) = align_N

          test_alignment(N,4,1,2) = align_W-(3*bc_size)
          test_alignment(N,4,1,1) = test_alignment(N,4,1,2)-small_layer
          test_alignment(N,4,2,1) = align_N
          test_alignment(N,4,2,2) = align_N

          test_alignment(N,5,1,1) = align_E+3*bc_size
          test_alignment(N,5,1,2) = test_alignment(N,5,1,1)+small_layer
          test_alignment(N,5,2,1) = align_N
          test_alignment(N,5,2,2) = align_N

          !south
          test_alignment(S,1,1,1) = align_W+1
          test_alignment(S,1,1,2) = test_alignment(S,1,1,1)+large_layer
          test_alignment(S,1,2,1) = align_S
          test_alignment(S,1,2,2) = align_S

          test_alignment(S,2,1,1) = test_alignment(S,1,1,2)+2*bc_size+small_layer
          test_alignment(S,2,1,2) = test_alignment(S,2,1,1)+large_layer
          test_alignment(S,2,2,1) = align_S
          test_alignment(S,2,2,2) = align_S

          test_alignment(S,3,1,2) = align_E-1
          test_alignment(S,3,1,1) = test_alignment(S,3,1,2)-2*bc_size-large_layer
          test_alignment(S,3,2,1) = align_S
          test_alignment(S,3,2,2) = align_S

          test_alignment(S,4,1,2) = align_W-(3*bc_size)
          test_alignment(S,4,1,1) = test_alignment(S,4,1,2)-small_layer
          test_alignment(S,4,2,1) = align_S
          test_alignment(S,4,2,2) = align_S
                         
          test_alignment(S,5,1,1) = align_E+3*bc_size
          test_alignment(S,5,1,2) = test_alignment(S,5,1,1)+small_layer
          test_alignment(S,5,2,1) = align_S
          test_alignment(S,5,2,2) = align_S
          
          !east
          test_alignment(E,1,1,1) = align_E
          test_alignment(E,1,1,2) = align_E
          test_alignment(E,1,2,1) = align_S+1
          test_alignment(E,1,2,2) = test_alignment(E,1,2,1)+large_layer

          test_alignment(E,2,1,1) = align_E
          test_alignment(E,2,1,2) = align_E
          test_alignment(E,2,2,1) = test_alignment(E,1,2,2)+2*bc_size+small_layer
          test_alignment(E,2,2,2) = test_alignment(E,2,2,1)+large_layer

          test_alignment(E,3,1,1) = align_E
          test_alignment(E,3,1,2) = align_E
          test_alignment(E,3,2,2) = align_N-1
          test_alignment(E,3,2,1) = test_alignment(E,3,2,2)-2*bc_size-large_layer


          !west
          test_alignment(W,1,1,1) = align_W
          test_alignment(W,1,1,2) = align_W
          test_alignment(W,1,2,1) = align_S+1
          test_alignment(W,1,2,2) = test_alignment(W,1,2,1)+large_layer

          test_alignment(W,2,1,1) = align_W
          test_alignment(W,2,1,2) = align_W
          test_alignment(W,2,2,1) = test_alignment(W,1,2,2)+2*bc_size+small_layer
          test_alignment(W,2,2,2) = test_alignment(W,2,2,1)+large_layer

          test_alignment(W,3,1,1) = align_W
          test_alignment(W,3,1,2) = align_W
          test_alignment(W,3,2,2) = align_N-1
          test_alignment(W,3,2,1) = test_alignment(W,3,2,2)-2*bc_size-large_layer

        end subroutine ini_alignment


        subroutine print_nbf_list(nbf_list_printed)
        
          implicit none

          type(nbf_list), intent(in) :: nbf_list_printed

          type(nbf_element), pointer :: current_element
          type(bf_sublayer), pointer :: bf_sublayer_ptr
          integer(ikind), dimension(2,2) :: alignment
          character(2)  , dimension(4)   :: bf_layer_char

          integer :: i

          bf_layer_char = ['N_','S_','E_','W_']

          current_element => nbf_list_printed%get_head()

          do i=1, nbf_list_printed%get_nb_elements()
             bf_sublayer_ptr => current_element%get_ptr()
             alignment       = bf_sublayer_ptr%get_alignment_tab()
             
             print '(A1,'': ('',I3,I3,'')  ('',I3,I3,'')'')',
     $            bf_layer_char(bf_sublayer_ptr%get_localization()),
     $            alignment(1,1), alignment(1,2),
     $            alignment(2,1), alignment(2,2)

             current_element => current_element%get_next()
          end do
          print '()'

        end subroutine print_nbf_list

      end module test_nbf_list_module

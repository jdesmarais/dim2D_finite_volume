      !test: cmd_operators_class
      program test_cmd_operators

        use cmd_operators_class, only :
     $     cmd_operators


        type(cmd_operators)   :: cmd_operators_used
        logical               :: restart_activated
        character(len=1024)   :: restart_filename
        integer, dimension(4) :: nb_bf_layers
        
        call cmd_operators_used%analyse_cmd_line_arg()

        restart_activated = cmd_operators_used%is_restart_activated()
        restart_filename  = cmd_operators_used%get_restart_filename()
        nb_bf_layers      = cmd_operators_used%get_nb_bf_layers()


        print '(''analysis of the command line arguments'')'
        print '(''----------------------------------------'')'
        print '(''restart_activated: '',L1)', restart_activated
        if(restart_activated) then
           print *, '-> restart_filename : ', trim(restart_filename)
        end if
        print '(''nb_bf_layers     : '',4I2)', nb_bf_layers
        print '()'

      end program test_cmd_operators

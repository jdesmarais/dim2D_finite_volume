      program test_cmd_operators_error_max

        use cmd_operators_error_max_class, only :
     $     cmd_operators_error_max

        implicit none

        type(cmd_operators_error_max) :: cmd_operators_used

        logical             :: inputs_provided
        character(len=1024) :: dir
        character(len=1024) :: filename
        integer             :: nb_files

        call cmd_operators_used%analyse_cmd_line_arg()

        inputs_provided = cmd_operators_used%are_inputs_provided()
        call cmd_operators_used%get_directory_error_files(dir)
        call cmd_operators_used%get_filename_error_max(filename)
        call cmd_operators_used%get_nb_files(nb_files)

        print *, 'are_inputs_provided:' , inputs_provided
        print *, 'directory_error_max: ', trim(dir)
        print *, 'filename_error_max : ', trim(filename)
        print *, 'nb_files           : ', nb_files


      end program test_cmd_operators_error_max

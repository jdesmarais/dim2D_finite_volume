#!/usr/bin/python -i

from library_bf_plot_visit import\
    create_grdpts_id_pictures,\
    create_movie,\
    create_field_pictures

import sys
import os
import getopt
sys.path.append(os.environ['VISIT_PYTHON_LIB'])
import visit


#describe the usage for using options
def usage():
    
    print '-h                            : display this help'
    print '-d dir                        : directory where the .nc files are saved'
    print '-g                            : create movie for the configuration of the grid-points'
    print '-f field                      : create movie for the field given (ex: mass)'
    print '--x_min=x_min                 : x_min for the field display'
    print '--x_max=x_max                 : x_max for the field display'
    print '--y_min=y_min                 : y_min for the field display'
    print '--y_max=y_max                 : y_max for the field display'
    print '--min=min                     : min for the field'
    print '--max=max                     : max for the field'
    print '--diff                        : create the movie for the error'
    print '--tsteps=[t_min,t_max,t_incr] : time steps analyzed'
    print '-w                            : with window'
    print '--load_var=path_var.xml       : file describing extra field variables'

    return


#parse the options
def parse_options(argv):
    
    try:
        opts, args = getopt.getopt(argv,
                                   "hgf:d:w",
                                   ["help",
                                    "grdpts_id",
                                    "field=",
                                    "dir=",
                                    "x_min=",
                                    "x_max=",
                                    "y_min=",
                                    "y_max=",
                                    "diff",
                                    "tsteps=",
                                    "min=",
                                    "max=",
                                    "load_var="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
    display_grdpts = False
    display_field  = False
    field_displayed = ''
    dirInputNcFiles = '/home/jdesmarais/Code/augeanstables/src/test_files/test_geometry_updates'
    x_min=-12
    x_max= 30
    y_min=-17
    y_max= 17
    diff = False
    no_window=False
    fieldMin='None'
    fieldMax='None'
    loadVar='None'

    for opt, arg in opts:

        if opt in ("-h","--help"):
            usage()
            os._exit(1) #sys.exit(0)
            
        elif opt in ("-g", "--grdpts_id"):
            display_grdpts=True
            
        elif opt in ("-f", "--field"):
            display_field = True
            field_displayed = arg
            
        elif opt in ("-d", "--dir"):
            dirInputNcFiles = arg
            
            if(not os.path.isdir(dirInputNcFiles)):
                print 'requested directory does not exist'
                print 'dir: ', dirInputNcFiles
                sys.exit(2)

        elif opt == "--x_min":
            x_min = float(arg)

        elif opt == "--x_max":
            x_max = float(arg)

        elif opt == "--y_min":
            y_min = float(arg)

        elif opt == "--y_max":
            y_max = float(arg)

        elif opt == "--min":
            fieldMin = float(arg)

        elif opt == "--max":
            fieldMax = float(arg)

        elif opt == "--diff":
            diff = True

        elif opt == "--tsteps":
            tsteps = arg.split('[')[1]
            tsteps = tsteps.split(']')[0]
            tsteps = tsteps.split(',')
            for i in range(0,len(tsteps)):
                tsteps[i]=int(tsteps[i])
            print 'timesteps: ', tsteps

        elif opt == "-w":
            no_window = True

        elif opt == "--loadVar":
            loadVar=arg
                

    print 'display_grdpts: ', display_grdpts
    print 'display_field: ', display_field
    if(display_field): print 'field_displayed: ', field_displayed
    print 'dir: ', dirInputNcFiles
    print ''


    #return the options in a dictionnary
    options = dict([('display_grdpts',display_grdpts),
                    ('display_field',display_field),
                    ('field_displayed',field_displayed),
                    ('dirInputNcFiles',dirInputNcFiles),
                    ('x_min',x_min),
                    ('x_max',x_max),
                    ('y_min',y_min),
                    ('y_max',y_max),
                    ('field_min',fieldMin),
                    ('field_max',fieldMax),
                    ('diff',diff),
                    ('timesteps',tsteps),
                    ('no_window',no_window),
                    ('')
                   ])
    return options


if __name__=="__main__":


    #main body
    argv = sys.argv[1:]
    options = parse_options(argv)

    if(options['no_window']):
        visit.LaunchNowin()
    else:
        visit.Launch()
    visit.SuppressMessages(1)


    #analyze options
    dirNc     = options['dirInputNcFiles']
    xmin      = options['x_min']
    xmax      = options['x_max']
    ymin      = options['y_min']
    ymax      = options['y_max']
    diff      = options['diff']
    timesteps = options['timesteps']
    

    #grdpts_id
    if(options['display_grdpts']):        
	
	#create the frames for the movie
	create_grdpts_id_pictures(
	    dirNc,
	    xmin,xmax,
	    ymin,ymax,
            diff=diff,
            timesteps=timesteps)
	
	#create the movie from the frames
	filePattern='grdpts_id%4d.png'
	movieFileName='grdpts_id.avi'
	
	create_movie(
	    os.path.join(dirNc,'visit_grdpts_id'),
	    filePattern,
	    movieFileName)

    #field
    if(options['display_field']):
        
	#field position
	field=options['field_displayed']
        fieldMin=options['field_min']
        fieldMax=options['field_max']
	
	create_field_pictures(
	    dirNc,
	    field,
	    fieldMin,
	    fieldMax,
	    xmin,xmax,
	    ymin,ymax,
            timesteps=timesteps,
            legend=False)
	#
	#filePattern=field+'%4d.png'
	#movieFileName=field+'.avi'
	#
	#create_movie(
	#    os.path.join(dirNc,'visit_'+field),
	#    filePattern,
	#    movieFileName)

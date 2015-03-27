#!/usr/bin/python -i

from library_bf_plot_visit import create_grdpts_id_pictures,create_movie,create_field_pictures

import sys
import os
import getopt
sys.path.append(os.environ['VISIT_PYTHON_LIB'])
import visit


#describe the usage for using options
def usage():
    
    print '-h       : display this help'
    print '-g       : create movie for the configuration of the grid-points'
    print '-f=field : create movie for the field'
    print '--x_min  : x_min for the field display'
    print '--x_max  : x_min for the field display'
    print '--y_min  : x_min for the field display'
    print '--y_max  : x_min for the field display'

    return

#parse the options
def parse_options(argv):
    
    try:
        opts, args = getopt.getopt(argv,
                                   "hgf:d:",
                                   ["help", "grdpts_id", "field=","dir=","x_min=","x_max=","y_min=","y_max="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
    display_grdpts = True
    display_field = False
    field_displayed = ''
    dirInputNcFiles = '/home/jdesmarais/Code/augeanstables/src/test_files/test_geometry_updates'
    x_min=-12
    x_max= 30
    y_min=-17
    y_max= 17

    for opt, arg in opts:

        if opt in ("-h","--help"):
            usage()
            sys.exit(2)
            
        elif opt in ("-g", "--grdpts_id"):
            display_grdpts=True
            
        elif opt in ("-f", "--field"):
            display_field = True
            field_displayed = arg
            
        elif opt in ("-d", "--dir"):
            dirInputNcFiles = arg
            
            if(not os.path.isdir(dirInputNcFiles)):
                print 'requested directory does not exists'
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
                    ('y_max',y_max)
                   ])
    return options


#main body
argv = sys.argv[1:]
options = parse_options(argv)
visit.Launch()

#grdpts_id
if(options['display_grdpts']):

    #analyze options
    movieStep=1
    dirNc= options['dirInputNcFiles']
    xmin = options['x_min']
    xmax = options['x_max']
    ymin = options['y_min']
    ymax = options['y_max']

    #create the frames for the movie
    create_grdpts_id_pictures(
        dirNc,
        xmin,xmax,
        ymin,ymax,
        step=movieStep)

    #create the movie from the frames
    filePattern='grdpts_id%4d.png'
    movieFileName='grdpts_id.avi'
    
    create_movie(
        os.path.join(dirNc,'visit_grdpts_id'),
        filePattern,
        movieFileName)

##field position
#field="position"
#fieldMin=-5
#fieldMax=1
#
#create_field_pictures(
#    dirInputNcFiles,
#    field,
#    fieldMin,
#    fieldMax,
#    xmin,
#    xmax,
#    ymin,
#    ymax,
#    step=movieStep)
#
#filePattern=field+'%4d.png'
#movieFileName=field+'.avi'
#
#create_movie(
#    os.path.join(dirInputNcFiles,'visit_'+field),
#    filePattern,
#    movieFileName,
#    movieRate=2)

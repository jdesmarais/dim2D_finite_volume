#!/usr/bin/python -i

from library_bf_plot_visit import create_grdpts_id_pictures,create_movie,create_field_pictures


#grdpts_id
movieStep=1
dirInputNcFiles="/home/jdesmarais/projects/dim2d_0.95_0.1_md0.0001_detailled/sm_domain"
xmin=-0.20
xmax= 0.60
ymin=-0.20
ymax= 0.20

create_grdpts_id_pictures(
    dirInputNcFiles,
    xmin,xmax,
    ymin,ymax,
    step=movieStep)


filePattern='grdpts_id%4d.png'
movieFileName='grdpts_id.avi'

create_movie(
    os.path.join(dirInputNcFiles,'visit_grdpts_id'),
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

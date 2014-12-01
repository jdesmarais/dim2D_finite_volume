#!/usr/bin/python -i

from library_bf_plot_visit import create_grdpts_id_pictures,create_movie,create_field_pictures


#grdpts_id
movieStep=10
dirInputNcFiles="/home/jdesmarais/projects/test_wave2d_hedstrom_xy_bf_apt2/"
xmin=-15
xmax= 15
ymin=-15
ymax= 15

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


#field position
field="position"
fieldMin=-5
fieldMax=1

create_field_pictures(
    dirInputNcFiles,
    field,
    fieldMin,
    fieldMax,
    xmin,
    xmax,
    ymin,
    ymax,
    step=movieStep)

filePattern=field+'%4d.png'
movieFileName=field+'.avi'

create_movie(
    os.path.join(dirInputNcFiles,'visit_'+field),
    filePattern,
    movieFileName,
    movieRate=2)

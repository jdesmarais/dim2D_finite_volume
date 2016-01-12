import sys
import os

sys.path.append(os.environ['VISIT_PYTHON_LIB'])
import visit

import inspect
import subprocess
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],
    "../")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from plot_bubble_with_error import\
    launch_visit,\
    postprocess_DIM_soundwave_with_error


if __name__=='__main__':

    # launch the graphical engine
    launch_visit(window=False)

    # set two windows
    visit.AddWindow()

    # postprocess for T=0.95, v=0.1, dct \in [2,8,12,16]
    for dct in [12,16]: #2,8,12,16]:

        # set the options
        dirNc = os.path.join('/home','jdesmarais',
                             'projects','jcp2015_submission',
                             '20150509_dim2d_bb_trans_cv_r3.5_lin',
                             'dim2d_0.95_0.1_dct'+str(dct))
        
        nbTimesteps = 1036
        
        soundMin =-0.05
        soundMax = 0.05
        errorMax = 5.0e-2
        bubbleValue = 1.022
        borders = [-0.25,0.25,-0.25,0.65]

        # postprocess
        postprocess_DIM_soundwave_with_error(
            dirNc,
            nbTimesteps,
            soundMin,
            soundMax,
            errorMax,
            bubbleValue,
            borders)
    

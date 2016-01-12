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

    # postprocess for T=0.99, v=0.05, md \in [0.05,0.1,0.2,0.3]
    for md in [0.3]: #0.05,0.1,0.2,0.3]:

        # set the options
        dirNc = os.path.join('/home','jdesmarais',
                             'projects','jcp2015_submission',
                             '20150509_dim2d_bb_trans_cv_r3.5_lin',
                             'dim2d_0.99_0.05_md'+str(md))
        
        nbTimesteps = 1017
        
        soundMin =-3.0e-3
        soundMax = 3.0e-3
        errorMax = 1.0e-3
        bubbleValue = 1.012
        borders = [-0.5,0.5,-0.5,1.2]

        # postprocess
        postprocess_DIM_soundwave_with_error(
            dirNc,
            nbTimesteps,
            soundMin,
            soundMax,
            errorMax,
            bubbleValue,
            borders)
    

#!/usr/bin/python

import unittest

from create_bb_tr_input import\
    compute_inputsToBeModified


# tests for methods in library_bb_tr_inputs
class TestCreate_bb_tr_inputs(unittest.TestCase):
    
    def test_compute_inputsToBeModified(self):

        inputs = compute_inputsToBeModified(temperature,
                                            flow_velocity,
                                            nb_pts_in_interface,
                                            ratio_interface_separation,
                                            ratio_bubble_interface,
                                            ratio_interface_influence,
                                            CFL_constant,
                                            total_nb_files,
                                            dct_distance,
                                            md_threshold_ac,
                                            md_threshold,
                                            flow_direction='x')

        # check whether the parameters in dim2d_parameters.f
        # suits well the parameters for the test

        # extract length_c, dim2d_a, dim2d_b, dim2d_M, dim2d_cv, dim2d_R
        # and dim2d_K from the dim2d_parameters.f fortran file
        dim2dParamPath = os.path.join(os.getenv('augeanstables'),
                                      'src',
                                      'physical_models',
                                      'dim2d',
                                      'dim2d_parameters.f')

        if(os.path.isfile(dim2dParamPath)):
            length_c  = float(get_parameter('length_c', dim2dParamPath))
            dim2d_a   = float(get_parameter( 'dim2d_a', dim2dParamPath))
            dim2d_b   = float(get_parameter( 'dim2d_b', dim2dParamPath))
            dim2d_M   = float(get_parameter( 'dim2d_M', dim2dParamPath))
            dim2d_cv  = float(get_parameter('dim2d_cv', dim2dParamPath))
            dim2d_R   = float(get_parameter( 'dim2d_R', dim2dParamPath))
            dim2d_K   = float(get_parameter( 'dim2d_K', dim2dParamPath))

        else:
            raise ValueError('test_create_bb_tr_input',
                             'test_compute_inputsToBeModified',
                             'dim2dParamPath does not exist')

        length_c_default = 

        inputParamFitTest = True
        inputParamFitTest = inputParamFitTest.and.(assertAlmostEqual(length_c,))
        
        if(not inputParamFitTest):
            print 'lenght_c', assertAlmostEqual(length_c,)
            
            raise valueError('parameters in dim2d_parameters do not fit test')
        


        self.assertAlmostEqual()

# main run
if __name__=="__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestLibrary_bb_tr_inputs)
    unittest.TextTestRunner(verbosity=2).run(suite)

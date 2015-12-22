#!/usr/bin/python

import unittest

from library_bb_tr_inputs import\
    get_domain_sizes,\
    get_domain_extent,\
    get_simulation_time


# tests for methods in library_bb_tr_inputs
class TestLibrary_bb_tr_inputs(unittest.TestCase):
    
    def test_get_domain_sizes(self):

        [Lx,Ly] = get_domain_sizes(2.0,4.0,5.0,0.25)

        self.assertAlmostEqual(Lx,45.0)
        self.assertAlmostEqual(Ly,42.0)

    def test_get_domain_extent(self):
        
        domain_extent = get_domain_extent(3.5,5.0,1.0,2.5)
        
        self.assertAlmostEqual(domain_extent[0][0],-2.0)
        self.assertAlmostEqual(domain_extent[1][0], 2.0)
        self.assertAlmostEqual(domain_extent[0][1],-2.5)
        self.assertAlmostEqual(domain_extent[1][1], 2.5)

    def test_get_simulation_time(self):

        self.assertAlmostEqual(get_simulation_time(6.0,1.0,0.5,2.0,2.0),
                               4.125)


# main run
if __name__=="__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestLibrary_bb_tr_inputs)
    unittest.TextTestRunner(verbosity=2).run(suite)

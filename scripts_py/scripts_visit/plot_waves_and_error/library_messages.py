#!/usr/bin/python

import sys

'''
@description: module with functions useful to print messages on screen
   - print_mg_error: print error message
   - print_mg_progress: print a message with return carriage
   - print_mg_final: print a message w/o return carriage
'''

# error message print
def print_mg_error(error_mg):
    sys.stdout.write('****'+error_mg+'**** \n')


# progress message print
def print_mg_progress(mg_progress):
    '''
    @description:
    print a message which is overwritten
    '''
    sys.stdout.write('%s\r' % mg_progress)
    sys.stdout.flush()


# final progress message print
def print_mg_final(mg_progress):
    '''
    @description:
    print a message which is not overwritten
    '''
    sys.stdout.write('%s' % mg_progress+'\n')
    sys.stdout.flush()

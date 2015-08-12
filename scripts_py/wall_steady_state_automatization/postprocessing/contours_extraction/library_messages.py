#!/usr/bin/python -i

import sys

def create_mg_progress(mg_progress):
    '''
    @description:
    print a message which is overwritten
    '''
    sys.stdout.write('%s\r' % mg_progress)
    sys.stdout.flush()


def create_mg_final(mg_progress):
    '''
    @description:
    print a message which is not overwritten
    '''
    sys.stdout.write('%s' % mg_progress)
    sys.stdout.flush()
    print '\n'

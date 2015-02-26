#!/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def iround(x):
    '''
    @description
    round a number to the nearest integer.
    '''
    y = round(x) - .5
    return int(y) + (y > 0)


def get_decomposition(x,base=2):
    '''
    @param
    x: integer
    @param
    base: integer
    @description
    get the decomposition of x in its base
    '''

    if(x==0):
        n=1
        d_array = np.zeros([1])

    else:
        n = int(math.log(x,base))+1
        
        d_array = np.empty([n])

        y=x
        for i in range(0,n):
            
            r = y%(base**(n-i-1))
            d = int((y-r)/(base**(n-i-1)))

            d_array[i] = d
            
            y = r

    return d_array


def get_grdpt_configuration(x):
    '''
    @param
    x: configuration: \f$x \in [1,2^8]\f$
    @return
    a: grdpt configuration
    '''

    debug = False

    d = get_decomposition(x,base=2)
    
    if(len(d)>8):
        print 'get_grdpt_configuration'
        print 'error: 0=<x<=254'
        sys.exit(0)

    c = np.zeros([9])

    if(len(d)>4):
        d1 = d[0:len(d)-4]
        d2 = d[len(d)-4:len(d)]
    else:
        d1 = np.zeros([4])
        d2 = d


    c[4-len(d1):4] = d1
    c[9-len(d2):9] = d2
    
    if(debug):
        print 'd:  ', d
        print 'd1: ', d1
        print 'd2: ', d2
        print c

    c = c[::-1]
    c = c.reshape([3,3])
    c = c[::-1,::]

    return c

if __name__ == "__main__":

    #a = get_grdpt_configuration(input())
    #
    #fig=plt.figure(figsize=(12,6))
    #ax = fig.add_subplot(1,1,1)
    #res = ax.imshow(a[:,:], cmap=cm.summer, interpolation='nearest')
    #fig.colorbar(res)
    #plt.show()

    fig=plt.figure(figsize=(12,6))
    fig.canvas.set_window_title("new grdpt configurations")

    for k in range(0,256):

        a = get_grdpt_configuration(k)

        ax = fig.add_subplot(16,16,k+1)
        res = ax.imshow(a[:,:], cmap=cm.summer, interpolation='nearest')
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_xaxis().set_label(str(k))

    plt.show()
    

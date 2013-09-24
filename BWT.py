#BWT.py
'''
Burrows-Wheeler Transform.
See http://www.homolog.us/blogs/blog/2011/10/03/finding-us-in-homolog-us/
'''

import numpy as np

sequence = "homolog.us$"

def BWT(sequence):
    n = len(sequence)
    seq = list(sequence)
    my_array = np.zeros([n,n],dtype = str)
    my_array[:,0] = seq
    for i in range(1,n):
        seq = seq[1:]+[seq[0]]
        my_array[:,i] = seq
    sortedlist = sorted(my_array.tolist())
    sortedarray = np.array(sortedlist)
    BWT ="".join(sortedarray[:,-1].tolist())
    print 'BWT:',BWT
BWT(sequence)

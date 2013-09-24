#BWT.py
'''
Burrows-Wheeler Transform.
See http://www.homolog.us/blogs/blog/2011/10/03/finding-us-in-homolog-us/
and
http://www.cs.nthu.edu.tw/~wkhon/ds/ds10/tutorial/tutorial7.pdf
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
    return BWT
BWT = BWT(sequence)



def Occ(c,q,BWT):
    return BWT[0:q].count(c)


def ctable(sequence):
    sortedsequence = sorted(sequence)
    sortedunique = sorted(list(set(sortedsequence)))
    charcount = [sortedsequence.count(character) for character in sortedunique]
    C =[]
    C.append((sortedunique[0],0))
    for i in range(1,len(sortedunique)):
        C.append((sortedunique[i],sum(charcount[0:i])))
    return C
seq = 'mississippi$'
C = ctable(seq)
print 'C:',C

def lasttofront(C,BWT,i):
    '''LF(i) = C[BWT[i]]+  Occ(BWT[i,i]'''
    LF = ctable(BWT[i])[i][1] + Occ(BWT[i],i)
    return LF
    




def BWTsearch(query,subject):
    pass

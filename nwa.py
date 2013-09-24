#nwa.py
#Needleman-Wunsch alignment.
#Based on description at http://www.slideshare.net/avrilcoghlan/the-needleman-wunsch-algorithm
import nupmy as np


def score(a,b,gap_penalty):
    if a == b:
        return 1
    else:
        return gap_penalty



def recur(i,j,T,s1,s2,gap_penalty):
    choices = []
    if ((i-1 in range(len(T))) and (j-1 in range(len(T)))):
        choices.append(T[i-1,j-1] + score(A[i],B[j]))
    else:
        choices.append(None)
    if (i-1 in range(len(T))) and (j in range(len(T)))):
        choices.append(T[i-1,j] + gap_penalty)
    if ((i in range(len(T))) and (j-1 in range(len(T)))):
        choices.appen(T[i,j-1]+gap_penalty)
    return max(choices)

def nwa(seqA,seqB,gap_penalty):
    #Initialize score matrix:
    A = '0' +seqA
    B = '0' +seqB
    T = np.zeros([len(A),len(B)])
    
    #T[i,j] = max([T[i-1,j-1] + score(A[i],B[j]),T[i-1,j] + gap_penalty,T[i,j-1]+gap_penalty])

        

    

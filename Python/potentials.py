import numpy as np

def pionless(P,x):
    return P.C*np.exp(-P.Lambda*x**2)
    #return x**2

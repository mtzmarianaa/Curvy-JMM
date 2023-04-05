# Gradient sampling algorithm for optimization on a triangle fan\
# without tops
# From Burke, Lewis, Overton's paper
# A robost gradient sampling algorithm for nonsmooth, nonconvex optimization

import numpy as np
from numpy.linalg import norm
from math import sqrt, pi
import intermediateTests as itt
from scipy.optimize import root_scalar

def sampleEpsBall(r, m, epsilonk):
    '''
    Draw m samples from the ball with radius epsilonk
    in R^r. According to the paper we must have m>r
    '''
    uk = []
    for j in range(r):
        random_vec = np.random.rand(r) - 0.5
        random_vec = random_vec/norm(random_vec)
        uk.append( epsilonk*random_vec )
    return uk







import numpy as np
import scipy.special as sc


def beta_entropy(alpha, beta):
    """See wikipedia"""
    result = np.log(sc.beta(alpha, beta))-(alpha-1)*sc.psi(alpha)-(beta-1)*sc.psi(beta)+\
             (alpha+beta-2)*sc.psi(alpha+beta)
    return result


def normal_entropy(sigma):
    result = 0.5 * np.log(2*np.pi*sigma**2) + 0.5
    return result


def bernoulli_entropy(p):
    result = -(1-p)*np.log(1-p) - p*np.log(p)
    return result


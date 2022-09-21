import typing
import numpy as np


def truncated_stick_breaking(n, k, v_prior):
    """See section 3.2"""
    v_s = np.random.beta(v_prior[0], v_prior[1], size =k)
    pi_s = np.cumprod(v_s)
    z_mat = np.zeros((n,k))
    for col in range(k):
        z_mat[:, col] = np.random.binomial(1, pi_s[col], n)
    return v_s, pi_s, z_mat


def linear_gaussian(
        n: int = 500,
        d: int = 500,
        k: int = 20,
        mu_a: float = 0,
        sigma_a: float = 1,
        sigma_eps: float = 0.1,
        v_prior: typing.Tuple[float, float] = (2,1)
        ):
    """
        n: number of observations
        d: data dimensionality
        k: num features
        v_prior: beta prior for stick breaking
    """
    v, pi, z = truncated_stick_breaking(n, k, v_prior)
    a = np.random.normal(mu_a, sigma_a, (k,d))
    epsilon = np.random.normal(0, sigma_eps)
    x = np.matmul(z,a) + epsilon
    result = {"v": v, "pi": pi, "Z": z, "A": a, "eps": epsilon, "X": x}
    return result
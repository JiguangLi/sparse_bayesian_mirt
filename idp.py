from generate import truncated_stick_breaking
import numpy as np
from scipy.stats import beta, norm, bernoulli


class IDB(object):
    """Indian Buffet Process with Variational Inference"""
    def __init__(self,
                 obs_data: np.array,
                 alpha: float = 2.0,
                 sigma_a: float = 1.0,
                 sigma_eps: float = 0.1,
                 k: int = 30
                 ):
        self.data = obs_data
        self.alpha = alpha
        self.sigma_a = sigma_a
        self.sigma_eps = sigma_eps
        self.n, self.d = obs_data.shape
        self.k = k
        # initialize latent variables of interest
        self.v_vars = np.random.beta(self.alpha, 1, size=self.k)
        self.a_mat = np.random.normal(0, self.sigma_a, (self.k, self.d))
        self.z_mat = truncated_stick_breaking(self.n, self.k, (1, self.alpha))[-1]

    def compute_elbo(self):
        # entropy term H(q)


        # epy_v = -np.mean(beta.logpdf(self.v_vars, self.alpha, 1))
        # epy_a = -np.mean(norm.logpdf(self.a_mat, loc=0, scale=self.sigma_a))
        # pi_terms = np.cumprod(self.v_vars)
        # epy_z = 0
        # for i, pi in enumerate(pi_terms):
        #     epy_z += np.sum(bernoulli.logpmf(self.z_mat[:,i], pi))
        # epy_z = -epy_z/(self.n*self.k)
        # epy = epy_v+epy_a+epy_z
        # full joint probability term














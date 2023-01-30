import numpy as np
import typing
import scipy.special as sc

class IBP(object):
    """Indian Buffet Process with Variational Inference"""
    def __init__(self,
                 obs_data: np.array,
                 alpha: float = 2.0,
                 sigma_a: float = 1.0,
                 sigma_eps: float = 0.1,
                 k: int = 30,
                 ):
        self.x = obs_data
        self.alpha = alpha
        self.sigma_a = sigma_a
        self.sigma_eps = sigma_eps
        self.n, self.d = obs_data.shape
        self.k = k
        # initialize variational parameters
        self.q_phi_mu = np.zeros((self.k, self.d))
        self.q_phi_cov = np.tile(np.identity(self.d), (self.k, 1, 1))
        self.q_vi = np.ones((self.n, self.k)) * 0.5
        self.q_pi_beta1 = np.ones(self.k) * 2
        self.q_pi_beta2 = np.ones(self.k)
        self.sigma_n = np.std(self.x)

    def precompute_q_ks(self):
        result = np.zeros(self.k)
        psi_beta1_cum_sum = np.append(0, np.cumsum(sc.psi(self.q_pi_beta1)))
        psi_sum_cum_sum = np.cumsum(sc.psi(self.q_pi_beta1 + self.q_pi_beta2))
        for i in range(self.k):
            result[i] = np.exp(sc.psi(self.q_pi_beta2[i]) + psi_beta1_cum_sum[i] + psi_sum_cum_sum[i])
        result = result / np.sum(result)
        return result

    def multinomial_approximation(self):
        """ Compute E[log(1-prod_{1}^{k} v_m)] for each k"""
        q_ks = self.precompute_q_ks()
        term1 = np.cumsum(q_ks * sc.psi(self.q_pi_beta2))
        term2 = np.append(0, np.cumsum(np.cumsum(q_ks[1:]) * sc.psi(self.q_pi_beta1[:-1])))
        term3 = np.cumsum(np.cumsum(q_ks) * sc.psi(self.q_pi_beta1 + self.q_pi_beta2))
        term4 = np.cumsum(q_ks * np.log(q_ks))
        return term1 + term2 - term3 - term4

    def compute_elbo(self):
        # Joint probability dist of X and W (see p12 of the technical report for derivation)
        vi_term = self.k * np.log(self.alpha) + \
                  (self.alpha - 1) * np.sum(sc.psi(self.q_pi_beta1) - sc.psi(self.q_pi_beta2 + self.q_pi_beta1))
        mn_approx = self.multinomial_approximation()
        inner_psi_term = np.cumsum(sc.psi(self.q_pi_beta2) - sc.psi(self.q_pi_beta2 + self.q_pi_beta1))
        z_term = np.sum(self.q_vi * inner_psi_term + (1 - self.q_vi) * mn_approx)
        a_term = np.sum([
            -(np.trace(self.q_phi_cov[k]) + np.dot(self.q_phi_mu[k, :], self.q_phi_mu[k, :])) / (2 * self.sigma_a ** 2)
            for k in range(self.k)]) - np.log(2 * np.pi * self.sigma_a ** 2) * self.d * self.k / 2
        # TODO: optimize x term implementation
        x_term = - np.log(2 * np.pi * self.sigma_n ** 2) * self.d / 2
        for n in range(self.n):
            temp = np.dot(self.x[n, :], self.x[n, :])
            temp = temp - 2 * np.dot(np.matmul(self.q_phi_mu, self.x[n, :].T), self.q_vi[n, :])
            for k in range(self.k - 1):
                for k_prime in range(k + 1, self.k):
                    temp = temp + 2 * self.q_vi[n, k] * self.q_vi[n, k_prime] * np.dot(self.q_phi_mu[k, :],
                                                                                       self.q_phi_mu[k_prime, :])
                temp = temp + self.q_vi[n, k] * (np.trace(self.q_phi_cov[k]) + np.dot(self.q_phi_mu[k, :],
                                                                                      self.q_phi_mu[k, :]))
            x_term = x_term - temp / (2 * self.sigma_n ** 2)
        # entropy part H[q]
        epy_vi = np.sum(
            np.log(
                sc.gamma(self.q_pi_beta1) * sc.gamma(self.q_pi_beta2) /
                sc.gamma(self.q_pi_beta1 + self.q_pi_beta2)
            ) - (self.q_pi_beta1 - 1) * sc.psi(self.q_pi_beta1)
            - (self.q_pi_beta2 - 1) * sc.psi(self.q_pi_beta2)
            + (self.q_pi_beta1 + self.q_pi_beta2 - 2) * sc.psi(self.q_pi_beta1 + self.q_pi_beta2)
        )
        epy_a = np.sum([self.d / 2 * (np.log(2 * np.pi * np.linalg.det(self.q_phi_cov[k])) + 1) for k in range(self.k)])
        epy_z = np.sum(-self.q_vi * np.log(self.q_vi)) - np.sum((1 - self.q_vi) * np.log(1 - self.q_vi))
        result = vi_term + z_term + a_term + x_term + epy_vi + epy_a + epy_z
        return result

    def update_a(self):
        for k in range(self.k):
            shared_term = (1/self.sigma_a**2 + np.sum(self.q_vi[:, k])/self.sigma_n**2)**-1
            self.q_phi_cov[k] = shared_term * np.identity(self.d)
            phi_mu = 0
            for n in range(self.n):
                phi_mu += self.q_vi[n, k]*(
                        self.x[n, :] - np.matmul(self.q_vi[n, :], self.q_phi_mu) + self.q_vi[n,k] * self.q_phi_mu[k, :]
                )
            self.q_phi_mu[k, :] = phi_mu * shared_term / self.sigma_n ** 2

    def update_z(self):
        mn_approx = self.multinomial_approximation()
        psi_diff = np.cumsum(sc.psi(self.q_pi_beta1) - sc.psi(self.q_pi_beta1+self.q_pi_beta2))
        trace_terms = [-(np.trace(self.q_phi_cov[k]) +
                         np.dot(self.q_phi_mu[k, :], self.q_phi_mu[k, :]))/(2*self.sigma_n**2)
                       for k in range(self.k)]
        for n in range(self.n):
            n_temp = self.x[n, :] - np.matmul(self.q_vi[n, :], self.q_phi_mu)
            for k in range(self.k):
                k_temp = mn_approx[k] - psi_diff[k] - trace_terms[k]
                total = k_temp + (1/self.sigma_n**2)*np.dot(self.q_phi_mu[k, :],
                                                            n_temp + self.q_vi[n, k] * self.q_phi_mu[k, :])
                self.q_vi[n, k] = 1/(1+np.exp(-total))

    def update_v(self):
        q_ks = self.precompute_q_ks()
        for k in range(self.k):
            if k == self.k-1:
                pi_1_prod_sum = 0
            else:
                pi_1_prod_sum = sum([(self.n - np.sum(self.q_vi[:, m]))*(np.sum(q_ks[k+1:m+1]))
                                     for m in range(k+1, self.k)])
            self.q_pi_beta1[k] = self.alpha + np.sum(self.q_vi[:, k:]) + pi_1_prod_sum
            self.q_pi_beta2[k] = 1 + sum([(self.n - np.sum(self.q_vi[:, m]))*q_ks[m] for m in range(k, self.k)])













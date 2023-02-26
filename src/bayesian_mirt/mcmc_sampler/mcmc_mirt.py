import typing
import numpy as np
from numpy.linalg import inv
import pandas as pd
from polyagamma import random_polyagamma
from tqdm import tqdm
from scipy.stats import norm
np.random.seed(0)
class MCMC_MIRT:
    def __init__(
            self,
            data: typing.Union[np.ndarray, pd.DataFrame],
            dim: int,
            num_samples: int,
            alpha_prior: str = "normal",
            loading_constraints: typing.Dict[typing.Tuple[int, int], float] = {}
    ):
        """
        Create a Bayesian MCMC object for Exploratory MIRT

            Args:
                data: a numpy array, with one row per participant and one column per response.
                      responses should be coded as positive integers. Shape: (n, m)
                dim: number of factors
                num_samples: num of MCMC draws
                alpha_prior: normal, or spike-and-slab normal prior
                loading_constraints: dictionary fixing (j,k) entry of the loading matrix

            Raises:
                ValueError: if given inputs do not meet object expectations.
        """
        if isinstance(data, pd.DataFrame):
            self.data = data.values
        elif isinstance(data, np.ndarray):
            self.data = data
        else:
            raise ValueError(f"Input data has unrecognised type {type(data)}, try a numpy.array instead.")
        if alpha_prior not in ["normal", "ss"]:
            raise ValueError("Invalid prior for alpha - 'normal' or 'ss'")

        self.n = self.data.shape[0]
        self.m = self.data.shape[1]
        self.k = dim
        self.alpha_prior = alpha_prior
        self.n_s = num_samples
        self.params = {}
        self.loading_constraints = loading_constraints
        if self.alpha_prior == "ss":
            self.ss_prior_prob = 0.5
            self.spike_var, self.slab_var = 0.01, 1

    def initialize_thetas(self):
        """Initialize Theta Params"""
        self.params["thetas"] = np.zeros((self.n_s + 1, self.n, self.k), dtype=float)
        theta_mu = np.zeros(self.k)
        theta_cov = np.identity(self.k)
        self.params["thetas"][0] = np.random.multivariate_normal(theta_mu, theta_cov, size=self.n)

    def initialize_alphas(self):
        """Initialize Alpha(loading) Params. Alphas_holder record the most recent alpha updates"""
        self.params["alphas"] = np.zeros((self.n_s + 1, self.m, self.k), dtype=float)
        alpha_mu = np.zeros(self.k)
        alpha_cov = np.identity(self.k)
        self.params["alphas"][0] = np.random.multivariate_normal(alpha_mu, alpha_cov, size=self.m)
        for pos, val in self.loading_constraints.items():
            j_entry, k_entry = pos
            self.params["alphas"][:, j_entry, k_entry] = val
        self.params["alphas_holder"] = np.copy(self.params["alphas"][0])

    def initialize_ss_alphas(self):
        """Initialize Alpha(loading) Params with spike-and-slab Prior"""
        self.params["alphas"] = np.zeros((self.n_s + 1, self.m, self.k), dtype=float)
        alpha_mu = np.zeros(self.k)
        alpha_cov_spike, alpha_cov_slab = np.identity(self.k)*self.spike_var, np.identity(self.k)*self.slab_var
        alphas_spike = np.random.multivariate_normal(alpha_mu, alpha_cov_spike, size=self.m)
        alphas_slab = np.random.multivariate_normal(alpha_mu, alpha_cov_slab, size=self.m)
        self.params["ss"] = np.zeros((self.n_s + 1, self.m, self.k), dtype=int)
        self.params["ss"][0] = np.random.binomial(1, self.ss_prior_prob, size=(self.m, self.k))
        self.params["alphas"][0] = self.params["ss"][0] * alphas_slab + (1-self.params["ss"][0])*alphas_spike
        for pos, val in self.loading_constraints.items():
            j_entry, k_entry = pos
            self.params["alphas"][:, j_entry, k_entry] = val
            if val == 0:
                self.params["ss"][:, j_entry, k_entry] = 0
            else:
                self.params["ss"][:, j_entry, k_entry] = 1
        self.params["alphas_holder"] = np.copy(self.params["alphas"][0])

    def initialize_intercepts(self):
        """Initialize Intercept Params"""
        self.params["intercepts"] = np.zeros((self.n_s + 1, self.m, 1), dtype=float)
        self.params["intercepts"][0] = np.random.normal(0, 1, self.m).reshape(-1, 1)

    def initialize_pg(self):
        """Initialize Polya-Gamma Variables w_{ij} for Data Augmentation"""
        ## TODO: update Polya-Gamma: think whether we need to keep all PG sample
        self.params["pg"] = np.zeros((self.n_s + 1, self.n, self.m), dtype=float)
        self.params["pg"][0] = random_polyagamma(1, 0, size=(self.n, self.m))
        self.update_pg_zmat(self.params["pg"][0])

    def update_pg_zmat(self, pg_mat: np.ndarray):
        """Each entry of Z has form (y_{ij}-0.5)/w_{ij}"""
        self.params["z_mat"] = (self.data-0.5)/pg_mat

    def initialize(self) -> None:
        """Initalize All Parameters If We want to Run Full-scale MCMC"""
        self.initialize_thetas()
        self.initialize_intercepts()
        if self.alpha_prior == "normal":
            self.initialize_alphas()
        elif self.alpha_prior == "ss":
            self.initialize_ss_alphas()
        self.initialize_pg()

    def sample_pg(self, draw_idx: int, pred_mat: np.ndarray) -> None:
        """Sample polya-gamma variables"""
        self.params["pg"][draw_idx] = random_polyagamma(1, pred_mat)

    def sample_thetas(self, draw_idx: int,
                      load_mat: np.ndarray,
                      intercepts: np.ndarray):
        """Gibb Sampling Thetas

        Arg:
            load_mat: m by k loading matrix
            intercepts: m by 1 intercept
        """
        for i in range(self.n):
            omega_i = np.diag(self.params["pg"][draw_idx][i])
            z_i = self.params["z_mat"][i].reshape(-1,1)
            cov_mat = inv(np.matmul(np.matmul(load_mat.T, omega_i), load_mat) + np.identity(self.k))
            mu_theta = np.matmul(cov_mat, np.matmul(load_mat.T, np.matmul(omega_i, z_i-intercepts)))
            self.params["thetas"][draw_idx][i] = np.random.multivariate_normal(mu_theta.flatten(), cov_mat)

    def sample_intercepts(self, draw_idx: int, thetas: np.ndarray):
        """Gibb Sampling Intercepts

        Args:
            thetas: N by K latent traits

        Note: We assume the program sample intercepts before sampling alphas for now.
        """
        v_mat = self.params["z_mat"] - np.matmul(thetas, self.params["alphas"][draw_idx-1].T)
        inv_w = 1/self.params["pg"][draw_idx]
        for i in range(self.m):
            mu_s = v_mat[:, i]
            var_s = inv_w[:, i]
            lik_mu, lik_var = self.compute_gaussian_kernel(mu_s, var_s)
            pos_mu = lik_mu/(1+lik_var)
            pos_var = lik_var/(lik_var+1)
            self.params["intercepts"][draw_idx][i] = np.random.normal(pos_mu, pos_var)

    def sample_alphas(self, draw_idx: int, thetas: np.ndarray):
        """Gibb Sampling Loading Matrix Alphas (with Normal prior)
            Args:
                thetas: N by K latent traits

            Note: We assume the program sample intercepts before sampling alphas.
        """
        for j in range(self.m):
            for k in range(self.k):
                if (j, k) not in self.loading_constraints:
                    boolean_slicer = np.array([True]*self.k)
                    boolean_slicer[k] = False
                    v_inner_prod = np.matmul(thetas[:, boolean_slicer],
                                             self.params["alphas_holder"][j, boolean_slicer].reshape(-1, 1)).flatten()
                    v = self.params["z_mat"][:, j] - v_inner_prod - np.ones(self.n) * self.params["intercepts"][draw_idx][j][0]
                    theta_k = thetas[:, k]
                    mu_s = v/theta_k
                    var_s = (1/self.params["pg"][draw_idx][:, j])*(1/np.square(theta_k))
                    lik_mu, lik_var = self.compute_gaussian_kernel(mu_s, var_s)
                    pos_mu = lik_mu / (1 + lik_var)
                    pos_var = lik_var / (lik_var + 1)
                    new_draw = np.random.normal(pos_mu, pos_var)
                    self.params["alphas"][draw_idx][j, k] = new_draw
                    self.params["alphas_holder"][j, k] = new_draw

    def sample_ss_alphas(self, draw_idx: int, thetas: np.ndarray):
        """Gibb Sampling Loading Matrix Alphas (with spike-and-slab prior)
            Args:
                thetas: N by K latent traits

            Note: We assume the program sample intercepts before sampling alphas.
        """
        for j in range(self.m):
            for k in range(self.k):
                if (j, k) not in self.loading_constraints:
                    boolean_slicer = np.array([True]*self.k)
                    boolean_slicer[k] = False
                    v_inner_prod = np.matmul(thetas[:, boolean_slicer],
                                             self.params["alphas_holder"][j, boolean_slicer].reshape(-1, 1)).flatten()
                    v = self.params["z_mat"][:, j] - v_inner_prod - np.ones(self.n) * self.params["intercepts"][draw_idx][j][0]
                    theta_k = thetas[:, k]
                    mu_s = v/theta_k
                    var_s = (1/self.params["pg"][draw_idx][:, j])*(1/np.square(theta_k))
                    lik_mu, lik_var = self.compute_gaussian_kernel(mu_s, var_s)
                    ss_variable = self.params["ss"][draw_idx-1][j, k]
                    prior_var = self.slab_var*ss_variable + self.spike_var*(1-ss_variable)
                    pos_mu = lik_mu*prior_var/(lik_var + prior_var)
                    pos_var = lik_var*prior_var/(lik_var + prior_var)
                    new_draw = np.random.normal(pos_mu, pos_var)
                    self.params["alphas"][draw_idx][j, k] = new_draw
                    self.params["alphas_holder"][j, k] = new_draw
                    ss_prob = self.compute_ss_prob(new_draw)
                    self.params["ss"][draw_idx][j, k] = np.random.binomial(1, ss_prob)

    def fit(self):
        """Fit a Full-Scaled MCMC Model"""
        self.initialize()
        for i in tqdm(range(1, 1 + self.n_s)):
            pred_mat = np.matmul(self.params["thetas"][i-1],
                                 self.params["alphas"][i-1].T) + self.params["intercepts"][i-1].T
            self.sample_pg(i, pred_mat)
            self.update_pg_zmat(self.params["pg"][i])
            self.sample_thetas(i, self.params["alphas"][i-1], self.params["intercepts"][i-1])
            self.sample_intercepts(i, self.params["thetas"][i])
            if self.alpha_prior == "normal":
                self.sample_alphas(i, self.params["thetas"][i])
            elif self.alpha_prior == "ss":
                self.sample_ss_alphas(i, self.params["thetas"][i])

    def estimate_thetas(self,
                       fixed_alphas: np.ndarray,
                       fixed_intercepts: np.ndarray):
        """Fix factor Loadings, estimate latent traits only.
         Args:
            fixed_alphas: loading matrix, shape: m by k
            fixed_intercepts: intercept term: shape m by 1

        Call this function when we only want to estimate latent trait theta. Each iteration has two steps:
            Step 1: Gibb sampling PG augmentation
            Step 2: Gibb sampling theta
        """
        self.initialize_thetas()
        self.initialize_pg()
        for i in tqdm(range(1, 1 + self.n_s)):
            pred_mat = np.matmul(self.params["thetas"][i-1], fixed_alphas.T) + fixed_intercepts.T
            self.sample_pg(i, pred_mat)
            self.update_pg_zmat(self.params["pg"][i])
            self.sample_thetas(i, fixed_alphas, fixed_intercepts)

    def estimate_loadings(self, fixed_thetas: np.ndarray):
        """Fixed Latent Traits, Estimate Loading matrix and Intercepts
        Args:
            fixed_thetas: N by K numpy array
        """
        self.initialize_pg()
        if self.alpha_prior == "normal":
            self.initialize_alphas()
        elif self.alpha_prior == "ss":
            self.initialize_ss_alphas()
        else:
            raise NotImplementedError("alpha prior type unimplemented")
        self.initialize_intercepts()
        for i in tqdm(range(1, 1 + self.n_s)):
            pred_mat = np.matmul(fixed_thetas, self.params["alphas"][i-1].T) + self.params["intercepts"][i-1].T
            self.sample_pg(i, pred_mat)
            self.update_pg_zmat(self.params["pg"][i])
            self.sample_intercepts(i, fixed_thetas)
            if self.alpha_prior == "normal":
                self.sample_alphas(i, fixed_thetas)
            elif self.alpha_prior == "ss":
                self.sample_ss_alphas(i, fixed_thetas)

    def predict(self,
                burn_nums: int,
                fixed_thetas: typing.Optional[np.array] = None,
                fixed_alphas: typing.Optional[np.array] = None,
                fixed_intercepts: typing.Optional[np.array] = None,
                ):
        """Predict the probability of student i answering question j correctly using posterior means of prediction

           Args:
               data: students item responses to predict
               new_thetas: students abilities estimation for the new dataset. N*K
               new_alphas: user input alphas. (default None, using estimated ones): M*K
               new_intercepts: user input intercepts. (default None, using estimated ones): M*1

           Returns:
               probability_df: a dataframe of probabilities/score predictions with same dimension as input

           Raises:
               ValueError: if given arguments do not fit requirements.
        """
        samples_left = self.n_s + 1 - burn_nums
        if fixed_thetas is None:
            thetas = self.params["thetas"][burn_nums+1:]
        else:
            thetas = np.repeat(fixed_thetas[np.newaxis, :, :], samples_left, axis=0)

        if fixed_alphas is None:
            alphas = self.params["alphas"][burn_nums+1:]
        else:
            alphas = np.repeat(fixed_alphas[np.newaxis, :, :], samples_left, axis=0)

        if fixed_intercepts is None:
            intercepts = self.params["intercepts"][burn_nums+1:]
        else:
            intercepts = np.repeat(fixed_intercepts[np.newaxis, :, :], samples_left, axis=0)
        linear_term = np.matmul(thetas, alphas.transpose(0,2,1)) + intercepts.transpose(0,2,1)
        prob_mat = 1/(1 + np.exp(-linear_term))
        result = np.mean(prob_mat, axis=0)
        return result

    def thetas(self, burn_nums: int):
        """Posterior Mean of Thetas"""
        if "thetas" not in self.params:
            raise ValueError("No thetas parameter available: model hasn't been fit ")
        pos_mean = np.mean(self.params["thetas"][burn_nums+1:], axis=0)
        return pos_mean

    def alphas(self, burn_nums: int):
        """Posterior Mean of Alphas"""
        if "alphas" not in self.params:
            raise ValueError("No alphas parameter available: model hasn't been fit ")
        pos_mean = np.mean(self.params["alphas"][burn_nums + 1:], axis=0)
        return pos_mean

    def intercepts(self, burn_nums: int):
        """Posterior Mean of Intercepts"""
        if "intercepts" not in self.params:
            raise ValueError("No intercepts parameter available: model hasn't been fit ")
        pos_mean = np.mean(self.params["intercepts"][burn_nums + 1:], axis=0).flatten()
        return pos_mean

    def ss(self, burn_nums: int):
        """Posterior for spike-slab-normal Variable Selection Variable"""
        if "ss" not in self.params:
            raise ValueError("No ss  parameter available: model hasn't been fit ")
        pos_mean = np.mean(self.params["ss"][burn_nums + 1:], axis=0)
        return pos_mean

    def compute_ss_prob(self, alpha_draw):
        """Compute P(ss_variable=1 | alpha_{jk})"""
        p_alpha1 = norm.pdf(alpha_draw, 0, self.slab_var**0.5) * self.ss_prior_prob
        p_alpha0 = norm.pdf(alpha_draw, 0, self.spike_var**0.5) * (1-self.ss_prior_prob)
        return p_alpha1/(p_alpha1 + p_alpha0)

    @staticmethod
    def compute_gaussian_kernel(mu_s: np.array, var_s: np.array) -> typing.Tuple[float, float]:
        """Take arrays of gaussian means and variances to compute mean and variance of their density products"""
        mu_result, var_result = mu_s[0], var_s[0]
        for i in range(1, len(mu_s)):
            new_mu, new_var = mu_s[i], var_s[i]
            denominator = var_result + new_var
            mu_result = (var_result*new_mu + new_var*mu_result)/denominator
            var_result = (var_result*new_var) / denominator
        return mu_result, var_result






























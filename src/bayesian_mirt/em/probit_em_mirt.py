import typing
import numpy as np
from numpy.linalg import inv
import pandas as pd
from tqdm import tqdm
from .minimax_tilting import TruncatedMVN
import cvxpy as cp


class EMProbitMIRT:

    def __init__(
            self,
            data: typing.Union[np.ndarray, pd.DataFrame],
            dim: int,
            ibp_alpha: float,
            num_mc_samples: int,
            ssl_lambdas: typing.Tuple[float, float] = (50, 1),
            max_iterations: int = 1000,
            epsilon: float = 1e-4,
            loading_constraints: typing.Optional[typing.Dict[typing.Tuple[int, int], float]] = None,
            start: typing.Optional[typing.Dict[str, np.ndarray]] = None,
            random_state: int = 42,
    ):
        """
        Create an Expectation-Maximization Object for Probit Multidimensional Item Response Theory Model

            Args:
                data: a numpy array, with one row per participant and one column per response.
                      responses should be coded as positive integers. Shape: (n, m)
                dim: truncated dimension for Indian Buffet Process Stick-Breaking Construction
                ibp_alpha: alpha prior for Indian Buffet Process
                ssl_lambdas: lambda0 and lambda1 controlling the variance of laplace distributions
                max_iterations: max number of iterations for EM algorithm
                epsilon: stops the program when maximum absolute difference between two consecutive estimated
                         loadings are smaller than epsilon.
                loading_constraints: dictionary fixing (j,k) entry of the loading matrix
                start: start values for thetas, intercepts, and alphas. Must be a dictionary
                random_state: seed
            Raises:
                ValueError: if given inputs do not meet object expectations.
        """

        if isinstance(data, pd.DataFrame):
            self.data = data.values
        elif isinstance(data, np.ndarray):
            self.data = data
        else:
            raise ValueError(f"Input data has unrecognised type {type(data)}, try a numpy.array instead.")
        self.rng = np.random.default_rng(random_state)
        self.n = self.data.shape[0]
        self.m = self.data.shape[1]
        self.k = dim
        self.ibp_alpha = ibp_alpha
        self.num_samples = num_mc_samples
        self.lambda0, self.lambda1 = ssl_lambdas
        self.epsilon = epsilon
        self.max_iterations = max_iterations
        self.params = {}
        self.loading_constraints = loading_constraints
        self.start = start
        self.random_state = random_state

    def initialize_thetas(self):
        """Initialize Latent Traits/Factors"""
        if self.start:
            if "thetas" in self.start:
                self.params["thetas"] = self.start["thetas"]
        else:
            self.params["thetas"] = np.zeros((self.n, self.num_samples, self.k), dtype=float)

    def initialize_alphas(self):
        """Initialize Factor Loadings"""
        self.params["alphas"] = np.zeros((self.m, self.k), dtype=float)
        alpha_mu = np.zeros(self.k)
        alpha_cov = np.identity(self.k)
        if self.start:
            if "alphas" in self.start:
                self.params["alphas"] = self.start["alphas"]
        else:
            self.params["alphas"] = self.rng.multivariate_normal(alpha_mu, alpha_cov, size=self.m)
        for pos, val in self.loading_constraints.items():
            j_entry, k_entry = pos
            self.params["alphas"][j_entry, k_entry] = val
        self.params["alphas_prev"] = np.copy(self.params["alphas"])

    def initialize_gammas_and_cs(self):
        """Initialize Gammas Variable Using Stick-breaking"""
        v_s = self.rng.beta(self.ibp_alpha, 1, size=self.k)
        c_s = np.cumprod(v_s)
        self.params["c_params"] = c_s
        self.params["gammas"] = np.zeros((self.m, self.k))
        for i in range(self.k):
            self.params["gammas"][:, i] = self.rng.binomial(1, c_s[i], 40)
        if self.start:
            if "c_params" in self.start:
                self.params["c_params"] = self.start["c_params"]
            if "gammas" in self.start:
                self.params["gammas"] = self.start["gammas"]

    def initialize_intercepts(self):
        """Initialize Intercepts"""
        self.params["intercepts"] = np.zeros(self.m, dtype=float)
        if self.start:
            if "intercepts" in self.start:
                self.params["intercepts"] = self.start["intercepts"]
        else:
            self.params["intercepts"] = self.rng.normal(0, 1, self.m)
        self.params["intercepts_prev"] = np.copy(self.params["intercepts"])

    def initialize_parameters(self):
        """Initialize Model Parameters"""
        self.initialize_alphas()
        self.initialize_thetas()
        self.initialize_gammas_and_cs()
        self.initialize_intercepts()

    def sample_thetas(self):
        """Sample from the Posterior Distribution of Thetas for all i's
        The function starts by creating the following variables
            d1_tensor: shape m by k by n, containing D1 matrix for all i's
            d2_array: shape n by m, containing D2 matrix for all i's
            s_array: shape n by m, containing s matrix for all i's
        """
        d1_tensor = self.params["alphas"][:, :, None] * (2*self.data-1).T[:, None, :]
        d2_array = (2*self.data-1)*self.params["intercepts"].reshape(1, -1)
        s_array = np.transpose(np.sqrt(np.sum(np.square(d1_tensor), axis=1)+1))
        for i in range(self.n):
            d1i = d1_tensor[:, :, i]
            cov_v0 = np.identity(self.k) - d1i.T@(inv(d1i@d1i.T+np.identity(self.m)))@d1i
            cov_v1 = np.diag(1/s_array[i])@(d1i@d1i.T+np.identity(self.m))@np.diag(1/s_array[i])
            truc_lev = - np.diag(1/s_array[i])@(d2_array[i].reshape(-1, 1)).flatten()
            v0_s = self.rng.multivariate_normal(np.zeros(self.k), cov_v0, self.num_samples)
            v1_s = TruncatedMVN(np.zeros(self.m), cov_v1, truc_lev,
                                np.ones_like(truc_lev) * np.inf, seed=self.random_state).sample(self.num_samples)
            linear_t = d1i.T@(inv(d1i@d1i.T + np.identity(self.m)))@np.diag(s_array[i])
            self.params["thetas"][i] = v0_s + (linear_t @ v1_s).T

    def update_gammas(self):
        """Update Gamma Variables for Variable Selection"""
        lambda0, lambda1 = self.lambda0, self.lambda1
        lambda1_part = (lambda1/2 * np.exp(-lambda1 * np.abs(self.params["alphas"]))
                        )*self.params["c_params"].reshape(1, -1)
        lambda0_part = (lambda0/2 * np.exp(-lambda0 * np.abs(self.params["alphas"]))
                        )*(1-self.params["c_params"].reshape(1, -1))
        self.params["gammas"] = lambda1_part/(lambda1_part+lambda0_part)

    def e_step(self):
        """Perform E-Step of EM Algorithm"""
        self.sample_thetas()
        self.update_gammas()

    def constrain_loadings(self, item_idx):
        """Check if there is loading constraints for a given item index"""
        if not self.loading_constraints:
            return False
        else:
            items = [x[0] for x in self.loading_constraints.keys()]
            result = item_idx in items
            return result

    def optimize_loadings(self):
        """Maximizing Loglikelihood and Penalty Terms to Obtain Factor Loadings and Intercepts"""
        thetas = self.params["thetas"].reshape(self.n*self.num_samples, -1)
        l1_penalty = self.params["gammas"] * self.lambda1 + (1-self.params["gammas"]) * self.lambda0
        for j in range(self.m):
            alpha_j = cp.Variable((self.k, 1))
            d_j = cp.Variable((1, 1))
            y_j = np.repeat(self.data[:, j], self.num_samples).reshape(-1, 1)
            log_likelihood = cp.sum(cp.multiply(y_j, cp.log_normcdf(thetas@alpha_j+d_j)) +
                                    cp.multiply(1 - y_j, cp.log_normcdf(-(thetas@alpha_j + d_j))))
            penal_term = cp.norm(cp.multiply(alpha_j, l1_penalty[j, :].reshape(-1, 1)), 1)
            if self.constrain_loadings(j):
                constraints = []
                for pos, val in self.loading_constraints.items():
                    if pos[0] == j:
                        constraints += [alpha_j[pos[1]] == val]
                problem = cp.Problem(
                    cp.Maximize(log_likelihood / self.num_samples - penal_term - 0.5 * cp.square(cp.norm(d_j, 2))),
                    constraints
                )
            else:
                problem = cp.Problem(
                    cp.Maximize(log_likelihood/self.num_samples - penal_term - 0.5 * cp.square(cp.norm(d_j, 2)))
                )
            problem.solve()
            self.params["alphas"][j, :] = alpha_j.value.flatten()
            self.params["intercepts"][j] = d_j.value[0][0]

    def optimize_c_params(self):
        """MAP for C-params, the Probabilities of Inclusion"""
        c_params = cp.Variable((1, self.k))
        objective = cp.sum(cp.multiply(self.params["gammas"],
                                       cp.log(c_params)) + cp.multiply(1 - self.params["gammas"],
                                                                       cp.log(1 - c_params))
                           ) + (self.ibp_alpha - 1) * cp.log(c_params[:, -1])
        constraints = []
        for i in range(self.k):
            constraints += [
                c_params[:, i] >= 0,
                c_params[:, i] <= 1,
            ]
            if i < self.k-1:
                constraints += [c_params[:, i + 1] - c_params[:, i] <= 0]
        problem = cp.Problem(cp.Maximize(objective), constraints)
        problem.solve()
        self.params["c_params"] = c_params.value[0]

    def m_step(self):
        """Perform M-Step of EM Algorithm"""
        self.optimize_loadings()
        self.optimize_c_params()

    def fit(self):
        """Fit a Probit EM Algorithm"""
        self.initialize_parameters()
        for i in tqdm(range(self.max_iterations)):
            self.e_step()
            self.m_step()
            max_diff = max(np.max(np.abs((self.params["alphas"] - self.params["alphas_prev"]))),
                           np.max(np.abs((self.params["intercepts"] - self.params["intercepts_prev"]))))
            if max_diff < self.epsilon:
                break
            print(max_diff)
            self.params["alphas_prev"] = np.copy(self.params["alphas"])
            self.params["intercepts_prev"] = np.copy(self.params["intercepts"])

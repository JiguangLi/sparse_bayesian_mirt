import typing
import numpy as np
import pandas as pd

def irt_parameters(
    n: int, m: int, theta_normal_prior: typing.Tuple = (0, 1), alpha_uniform_prior: typing.Tuple = (0.3, 2.5)
) -> typing.Dict[str, typing.Any]:
    """Generate paraemters for an IRT model.

    Args:
        n: number of respondents
        m: number of items
        theta_normal_prior: tuple in form (mean, variance) of normal prior
        alpha_uniform_prior: tuple in form (?, ?) of normal prior

    Returns:
        A dictionary with parameter values.
    """
    np.random.seed(0)
    thetas = np.random.normal(*theta_normal_prior, n)
    alphas = np.random.uniform(*alpha_uniform_prior, m)
    intercepts = np.random.uniform(-3, 3, m)
    true_params = {"alphas": alphas, "thetas": thetas, "intercepts": intercepts}
    return true_params


def mirt_parameters(n: int, m: int, k: int, theta_sigma: float, alpha_sigma: float) -> typing.Dict[str, typing.Any]:
    """Generate paraemters for an IRT model.

    Args:
        n: num of students
        m: num of items
        k: num of dimensions
        theta_sigma: entries on the diagonal of covariance matrix for thetas
        alpha_sigma: entries on the diagonal of covariance matrix for alphas

    Returns:
        Dictionary with generated parameter values, including
            * alphas (k by m array),
            * thetas: n by k array
            * intercepts: 1 by m array
    """
    np.random.seed(0)
    thetas = np.random.multivariate_normal(np.zeros(k), np.identity(k) * theta_sigma ** 2, n)
    alphas = np.random.multivariate_normal(np.ones(k), np.identity(k) * alpha_sigma ** 2, m).T
    intercepts = np.random.uniform(-3, 3, m).reshape(1, -1)
    true_params = {"alphas": alphas, "thetas": thetas, "intercepts": intercepts}
    return true_params

def ss_mirt_parameters(
    n: int, m: int, k: int, zeros_per_dim: int = 0, theta_sigma: float = 1.0, alpha_sigma: float = 1.0,
        compensatory_flag: bool= True, lower_triangularize = False
) -> typing.Dict[str, typing.Any]:
    """Generate paraemters for an SS mIRT model.

    Args:
        n: num of respondents.
        m: num of items.
        k: num of dimensions.
        theta_sigma: entries on the diagonal of covariance matrix for thetas.
        alpha_sigma: entries on the diagonal of covariance matrix for alphas.
        zeros_per_dim: number of zeros per dimension
        compensartory_flag: whether the parameter is for compensatory model
    Returns:
        alphas: k by m array
        thetas: n by k array
        intercepts: 1 by m array for compensatory, k by m for non-compensatory
    """
    np.random.seed(0)
    if zeros_per_dim >= k or zeros_per_dim <0:
        raise ValueError("The number of zero dimensions is invalid, must between 0 and k")
    thetas = np.random.multivariate_normal(np.zeros(k), np.identity(k) * theta_sigma ** 2, n)
    if zeros_per_dim > 0:
        alphas_temp = np.random.multivariate_normal(np.ones(k), np.identity(k) * alpha_sigma ** 2, m).T
        zero_idx_x_coords = np.zeros(m * zeros_per_dim).astype(int)
        for i in range(m):
            zero_idx_x_coords[(i * zeros_per_dim) : ((i + 1) * zeros_per_dim)] = np.random.choice(
                np.arange(0, k, 1).astype(int), zeros_per_dim, replace=False
            ).astype(int)
        zero_idx_y_coords = np.repeat(np.arange(0, m).astype(int), zeros_per_dim).astype(int)
        boolean_mask = np.ones(alphas_temp.shape)
        boolean_mask[zero_idx_x_coords, zero_idx_y_coords] = 0
        alphas = alphas_temp * boolean_mask
    else: # zeros_per_dim ==0
        alphas = np.random.multivariate_normal(np.ones(k), np.identity(k) * alpha_sigma ** 2, m).T
        boolean_mask = np.nan
    if compensatory_flag:
        intercepts = np.random.normal(loc=0, scale=1.5, size=(1, m))
    else:
        intercepts= np.random.normal(loc=0, scale=1.5, size=(k, m))
    if lower_triangularize:
        alphas = np.tril(alphas.T)
        np.fill_diagonal(alphas, 1)
        alphas = alphas.T

    true_params = {"alphas": alphas, "thetas": thetas, "intercepts": intercepts, "boolean_mask": boolean_mask,
                   "compensatory_flag": compensatory_flag}
    return true_params



def threshold_ss_mirt_parameters(
    n: int, m: int, k: int, zeros_per_dim: int = 0, threshold= 0.2, lower_triangularize = False
) -> typing.Dict[str, typing.Any]:
    """Generate paraemters for an SS mIRT model threshold parameters.

    Args:
        n: num of respondents.
        m: num of items.
        k: num of dimensions.
        zeros_per_dim: number of zeros per dimension
    Returns:
        alphas: k by m array
        thetas: n by k array
        intercepts: 1 by m array for compensatory, k by m for non-compensatory
    """
    np.random.seed(0)
    if zeros_per_dim >= k or zeros_per_dim <0:
        raise ValueError("The number of zero dimensions is invalid, must between 0 and k")
    thetas = np.random.multivariate_normal(np.zeros(k), np.identity(k), n)
    if zeros_per_dim > 0:
        alphas_temp = np.random.multivariate_normal(np.ones(k), np.identity(k), m).T
        # make abs(arr) < threshold 0
        alphas_temp = alphas_temp*(np.abs(alphas_temp) > threshold)
        # fill zeros with uniform(threshold, threshold+0.5)
        alphas_temp[alphas_temp == 0] = np.random.uniform(threshold, threshold+0.5,  np.count_nonzero(alphas_temp == 0))
        zero_idx_x_coords = np.zeros(m * zeros_per_dim).astype(int)
        for i in range(m):
            zero_idx_x_coords[(i * zeros_per_dim) : ((i + 1) * zeros_per_dim)] = np.random.choice(
                np.arange(0, k, 1).astype(int), zeros_per_dim, replace=False
            ).astype(int)
        zero_idx_y_coords = np.repeat(np.arange(0, m).astype(int), zeros_per_dim).astype(int)
        boolean_mask = np.ones(alphas_temp.shape)
        boolean_mask[zero_idx_x_coords, zero_idx_y_coords] = 0
        alphas = alphas_temp * boolean_mask
    else: # zeros_per_dim ==0
        alphas = np.random.multivariate_normal(np.ones(k), np.identity(k), m).T
        boolean_mask = np.nan

    # make the first k intercepts positive to make sure the number of correct responses are sufficient
    intercepts1 = np.random.uniform(0.5, 1.5, k).reshape(1, -1)
    intercepts2 = np.random.uniform(-2, 2, m-k).reshape(1, -1)
    intercepts = np.hstack((intercepts1, intercepts2))

    if lower_triangularize:
        alphas = np.tril(alphas.T)
        np.fill_diagonal(alphas, 1)
        alphas = alphas.T

    true_params = {"alphas": alphas, "thetas": thetas, "intercepts": intercepts, "boolean_mask": boolean_mask}
    return true_params


def lower_trig_constraints(num_items:int, dim:int) -> typing.Dict[typing.Tuple[int, int], float]:
    """Return a dictionary where the entries of main diagonal is  1, and all the upper triangular part 0

    Note we assume num_items >= factor dimension dim
    """
    result = {}
    for i in range(dim):
        result[(i, i)] = 1
    for i in range(dim):
        for j in range(num_items):
            if i > j:
                result[(j, i)] = 0
    return result


def lower_trig_constraints_truncated(dim:int, truncated_level:int) -> typing.Dict[typing.Tuple[int, int], float]:
    """Return a dictionary where the entries of main diagonal is  1, and all the upper triangular part 0, and
    right part of triangular is also 0

    dim: true dimension

    Note we assume num_items >= factor dimension dim, truncated_level > dim
    """
    result = {}
    for i in range(dim):
        result[(i, i)] = 1
    for i in range(dim):
        for j in range(i+1, truncated_level):
            result[(i,j)] = 0
    return result
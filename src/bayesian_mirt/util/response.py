import typing
import numpy as np
import pandas as pd


def irt_binary_data(n: int, m: int, params: typing.Dict[str, typing.Any]):
    # generate samples
    np.random.seed(0)
    linear_term = np.matmul(params["thetas"].reshape(-1, 1), params["alphas"].reshape(1, -1)) + params[
        "intercepts"
    ].reshape(1, -1)
    prob_matrix = 1 / (1 + np.exp(-linear_term))
    result_matrix = (np.random.uniform(0, 1, (n, m)) < prob_matrix).astype(int)
    # create dataframe
    col_names = ["item" + str(i) for i in range(1, m + 1)]
    df = pd.DataFrame(data=result_matrix, columns=col_names)
    return df


def mirt_binary_data(n: int, m: int, params: typing.Dict[str, typing.Any]):
    # generate samples
    np.random.seed(0)
    if params["compensatory_flag"]:
        linear_term = np.matmul(params["thetas"], params["alphas"]) + params["intercepts"]
        prob_matrix = 1 / (1 + np.exp(-linear_term))
    else:
        prob_matrix = np.ones(n, m)
        dim = params["thetas"].shape[1]
        for i in range(dim):
            prob_matrix = (params["thetas"][:, i].reshape(-1, 1) - params["intercepts"][i, :].reshape(1, -1)) * (
                params["alphas"][i, :].reshape(1, -1)) * prob_matrix
    result_matrix = (np.random.uniform(0, 1, (n, m)) < prob_matrix).astype(int)
    col_names = ["item" + str(i) for i in range(1, m + 1)]
    df = pd.DataFrame(data=result_matrix, columns=col_names)
    return df
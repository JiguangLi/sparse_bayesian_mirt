from generate import ss_mirt_parameters, mirt_binary_data
import pandas as pd
import numpy as np
import pathlib
import yaml
import typing
import argparse
import pickle

def parse_args(verbose: bool = False) -> typing.Dict[str, typing.Any]:
    parser = argparse.ArgumentParser("Generate fake mirt data")
    parser.add_argument("--config_filepath", default=pathlib.Path("/Users", "jiguangli", "IBP_IRT", "config", "data_config.yaml"))
    parser.add_argument("--verbose", default=True)
    parser.add_argument("--n", help="number of responses", default=10000, type=int)
    parser.add_argument("--m", help="number of questions", default=40, type=int)
    parser.add_argument("--k", help="number of factors", default=3, type=int)
    parser.add_argument("--zeros_per_dim", help="number of zero loadings per dimension", default=2, type=int)

    arguments = vars(parser.parse_args())

    # get remaining arguments from config
    with open(arguments["config_filepath"]) as config_file:
        config = yaml.full_load(config_file)

    arguments.update(config)

    # fin
    if verbose:
        print(arguments)
    return arguments


if __name__ == '__main__':
    arguments = parse_args()
    n = arguments["n"]
    m = arguments["m"]
    k = arguments["k"]
    zeros_per_dim = arguments["zeros_per_dim"]
    params = ss_mirt_parameters(n, m, k, zeros_per_dim)
    data = mirt_binary_data(n, m, params)
    filename = "_".join(["fake", str(n), str(m), str(k), str(zeros_per_dim)]) + ".feather"
    data.to_feather(pathlib.Path(arguments["output_data_dir"]).joinpath(filename))
    params_filename = "_".join(["fake", str(n), str(m), str(k), str(zeros_per_dim), "params"]) + ".pickle"
    with open(pathlib.Path(arguments["output_data_dir"]).joinpath(params_filename), "wb") as handle:
        pickle.dump(params, handle, protocol=pickle.HIGHEST_PROTOCOL)
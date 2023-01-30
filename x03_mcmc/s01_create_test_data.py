import pandas as pd
import numpy as np
import pathlib
import yaml
import bayesian_mirt as bmirt
import argparse
import typing
import pickle


def parse_args(verbose: bool = False) -> typing.Dict[str, typing.Any]:
    """Get command line arguments, with defaults set for local testing."""
    # parse command line arguments first
    parser = argparse.ArgumentParser("Generate a few synthetic dataset for testing")
    parser.add_argument("--config_filepath",
                        default= pathlib.Path.home().joinpath(pathlib.Path("IBP_IRT", "config", "data_config.yaml")))
    parser.add_argument("--n", default=10000, type=int)
    parser.add_argument("--m", default=40, type=int)
    parser.add_argument("--verbose", default=True)
    arguments = vars(parser.parse_args())
    # get remaining arguments from config
    with open(arguments["config_filepath"]) as config_file:
        config = yaml.full_load(config_file)
    arguments.update(config)
    # fin
    if verbose:
        print(arguments)
    return arguments

def save_data(model_df: pd.DataFrame, model_params: typing.Dict[str, typing.Any],
              output_dir:str, file_name:str):
    param_filename= file_name+ "_params.pickle"
    data_filename = file_name+ ".feather"
    model_df.to_feather(pathlib.Path(output_dir).joinpath(data_filename))
    with open(pathlib.Path(output_dir).joinpath(param_filename), "wb") as handle:
        pickle.dump(model_params, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    # data setup
    arguments = parse_args()
    n, m = arguments["n"], arguments["m"]
    # fake dataset 1: irt model
    irt_params = bmirt.util.irt_parameters(n, m)
    irt_data = bmirt.util.irt_binary_data(n, m, params= irt_params)
    save_data(irt_data, irt_params, arguments["fake_data_dir"], "fake1")
    # fake dataset 2: 2 dimensional mirt without sparsity
    mirt_params2 = bmirt.util.ss_mirt_parameters(n, m, 2, zeros_per_dim=0)
    mirt_data2 = bmirt.util.mirt_binary_data(n, m, params=mirt_params2)
    save_data(mirt_data2, mirt_params2, arguments["fake_data_dir"], "fake2")
    # fake dataset 3: 2 dimensional mirt with one dimension 0
    mirt_params3 = bmirt.util.ss_mirt_parameters(n, m, 2, zeros_per_dim=1)
    mirt_data3 = bmirt.util.mirt_binary_data(n, m, params=mirt_params3)
    save_data(mirt_data3, mirt_params3, arguments["fake_data_dir"], "fake3")
    # fake dataset 4: 3 dimensional mirt without sparsity
    mirt_params4 = bmirt.util.ss_mirt_parameters(n, m, 3, zeros_per_dim=0)
    mirt_data4 = bmirt.util.mirt_binary_data(n, m, params=mirt_params4)
    save_data(mirt_data4, mirt_params4, arguments["fake_data_dir"], "fake4")
    # fake dataset 5: 3 dimensional mirt with one dimension 1
    mirt_params5 = bmirt.util.ss_mirt_parameters(n, m, 3, zeros_per_dim=1)
    mirt_data5 = bmirt.util.mirt_binary_data(n, m, params=mirt_params5)
    save_data(mirt_data5, mirt_params5, arguments["fake_data_dir"], "fake5")
    # fake dataset 6: 3 dimensional mirt with one dimension 2
    mirt_params6 = bmirt.util.ss_mirt_parameters(n, m, 3, zeros_per_dim=2)
    mirt_data6 = bmirt.util.mirt_binary_data(n, m, params=mirt_params6)
    save_data(mirt_data6, mirt_params6, arguments["fake_data_dir"], "fake6")









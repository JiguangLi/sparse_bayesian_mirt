import pandas as pd
import numpy as np
import pathlib
import yaml
import argparse
import typing
import pickle


def parse_args(verbose: bool = False) -> typing.Dict[str, typing.Any]:
    """Get command line arguments, with defaults set for local testing."""
    # parse command line arguments first
    parser = argparse.ArgumentParser("Extract Gibbs Inclusion Probability for Thresholding")
    parser.add_argument("--config_filepath",
                        default= pathlib.Path.home().joinpath(pathlib.Path("IBP_IRT", "config", "data_config.yaml")))
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


def read_small_thetas(small_model_output_dir, burn_in):
    with open(small_model_output_dir.joinpath("gibbs_large_k.pkl"), 'rb') as handle:
        model = pickle.load(handle)
        thetas_small.append(model.ss_thetas(burn_in))
    with open(small_model_output_dir.joinpath("gibbs_trig_large_k.pkl"), 'rb') as handle:
        model = pickle.load(handle)
        thetas_small.append(model.ss_thetas(burn_in))
    with open(small_model_output_dir.joinpath("gibbs_true_k.pkl"), 'rb') as handle:
        model = pickle.load(handle)
        thetas_small.append(model.ss_thetas(burn_in))
    with open(small_model_output_dir.joinpath("gibbs_trig_true_k.pkl"), 'rb') as handle:
        model = pickle.load(handle)
        thetas_small.append(model.ss_thetas(burn_in))


def read_ibp_thetas(ibp_model_output_dir, burn_in):
    with open(ibp_model_output_dir.joinpath("gibbs_large_k.pkl"), 'rb') as handle:
        model = pickle.load(handle)
        thetas_ibp.append(model.ss_thetas(burn_in))
    with open(ibp_model_output_dir.joinpath("gibbs_trig_large_k.pkl"), 'rb') as handle:
        model = pickle.load(handle)
        thetas_ibp.append(model.ss_thetas(burn_in))
    with open(ibp_model_output_dir.joinpath("gibbs_true_k.pkl"), 'rb') as handle:
        model = pickle.load(handle)
        thetas_ibp.append(model.ss_thetas(burn_in))
    with open(ibp_model_output_dir.joinpath("gibbs_trig_true_k.pkl"), 'rb') as handle:
        model = pickle.load(handle)
        thetas_ibp.append(model.ss_thetas(burn_in))


if __name__ == "__main__":
    # data setup
    arguments = parse_args()
    small_model_output_dir = pathlib.Path.home().joinpath(pathlib.Path("IBP_IRT", "models", "small_n"))
    ibp_model_output_dir = pathlib.Path.home().joinpath(pathlib.Path("IBP_IRT", "models", "ibp"))
    burn_in = 500
    thetas_small = []
    thetas_ibp = []
    read_small_thetas(small_model_output_dir, burn_in)
    read_ibp_thetas(ibp_model_output_dir, burn_in)
    model_types = ["large_k", "large_k_trig", "true_k", "true_k_trig"]
    col_names = ["Model_type", "small_n", "IBP"]
    thetas_df = pd.DataFrame(np.vstack((model_types, thetas_small, thetas_ibp)).T, columns=col_names)
    output_dir = pathlib.Path.home().joinpath(pathlib.Path("IBP_IRT", "paper_markdown", "synthetic_data"))
    thetas_df.to_feather(output_dir.joinpath("gibbs_inclusion_prob.feather"))








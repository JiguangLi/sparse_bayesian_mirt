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
    parser = argparse.ArgumentParser("Fit Gibbs Samplers on DESE data, 2022 grade 10")
    parser.add_argument("--config_filepath",
                        default= pathlib.Path.home().joinpath(pathlib.Path("IBP_IRT", "config", "data_config.yaml")))
    parser.add_argument("--verbose", default=True)
    parser.add_argument("--model_output_dir", default= pathlib.Path.home().joinpath(pathlib.Path("IBP_IRT", "models", "dese", "2022")))
    parser.add_argument("--year", default="22")
    parser.add_argument("--grade", default="10")
    arguments = vars(parser.parse_args())
    # get remaining arguments from config
    with open(arguments["config_filepath"]) as config_file:
        config = yaml.full_load(config_file)
    arguments.update(config)
    # fin
    if verbose:
        print(arguments)
    return arguments


def save_results(model, model_output_dir, model_name, burn_in):
    # save model
    with open(model_output_dir.joinpath(model_name+".pkl"), 'wb') as handle:
        pickle.dump(model, handle, pickle.HIGHEST_PROTOCOL)
    alphas = model.alphas(burn_in)
    intercepts = model.intercepts(burn_in)
    preds = model.predict(burn_in)
    adapt_ss = model.adapt_ss(burn_in)[0]
    # save loading
    loading_filename = model_name + "_loadings.feather"
    loading_col_names = ["dim" + str(i) for i in range(alphas.shape[1])] + ["intercepts"] + ["gamma" + str(j) for j in range(alphas.shape[1])]
    loading_df = pd.DataFrame(np.hstack((alphas, intercepts.reshape(-1, 1), adapt_ss)), columns= loading_col_names)
    loading_df.to_feather(model_output_dir.joinpath(loading_filename))
    # save Prediction
    pred_file_name = model_name + "_preds.feather"
    pred_col_names = ["item" + str(i) for i in range(1, intercepts.shape[0]+1)]
    pred_df = pd.DataFrame(preds, columns= pred_col_names)
    pred_df.to_feather(model_output_dir.joinpath(pred_file_name))


if __name__ == "__main__":
    # data setup
    arguments = parse_args()
    data_input_dir = pathlib.Path(arguments["dese_data_dir"]).joinpath(pathlib.Path("sample_processed", "2022"))
    data_filename = "grade_" + arguments["grade"] + "_sampled.feather"
    data = pd.read_feather(data_input_dir.joinpath(data_filename))
    model_output_dir = pathlib.Path.home().joinpath(pathlib.Path("IBP_IRT", "models", "dese", "2022", "gibbs"))
    item_names = [col for col in data.columns if "mitem" in col or "eitem" in col]
    large_k = 10
    burn_in = 500
    loading_constraints_large = bmirt.util.lower_trig_constraints(len(item_names), large_k)
    # Gibbs sampler, large k, not triangularize the loading
    model1 = bmirt.MCMC_MIRT(data[item_names], dim=large_k, num_samples=2000, alpha_prior="adaptive-ss",
                             loading_constraints={})
    model1.fit()
    save_results(model1, model_output_dir, "gibbs_grade" + arguments["grade"], burn_in)
    del model1
    # Gibbs sampler, large k, triangularize the loading
    model2 = bmirt.MCMC_MIRT(data[item_names], dim=large_k, num_samples=2000, alpha_prior="adaptive-ss",
                            loading_constraints = loading_constraints_large)
    model2.fit()
    save_results(model2, model_output_dir, "gibbs_trig_grade" + arguments["grade"], burn_in)
    del model2






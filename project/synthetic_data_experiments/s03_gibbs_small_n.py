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
    parser = argparse.ArgumentParser("Fit Gibbs Samplers on Small Overlap Data")
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
    data_input_dir = pathlib.Path.home().joinpath(pathlib.Path("IBP_IRT", "data", "small_n"))
    model_output_dir = pathlib.Path.home().joinpath(pathlib.Path("IBP_IRT", "models", "small_n"))
    data = pd.read_feather(pathlib.Path(data_input_dir).joinpath("small_n_logit.feather"))
    print(data.shape)
    num_items = data.shape[1]
    true_k = 4
    large_k = 10
    burn_in = 500
    loading_constraints_small = bmirt.util.lower_trig_constraints(num_items, true_k)
    loading_constraints_large = bmirt.util.lower_trig_constraints(num_items, large_k)
    # Gibbs sampler, correct k, not triangularize the loading
    model1 = bmirt.MCMC_MIRT(data, dim=true_k, num_samples=2000, alpha_prior="adaptive-ss", loading_constraints={})
    model1.fit()
    save_results(model1, model_output_dir, "gibbs_true_k", burn_in)
    del model1
    # Gibbs sampler, correct k, triangularize the loading
    model2 = bmirt.MCMC_MIRT(data, dim=true_k, num_samples=2000, alpha_prior="adaptive-ss",
                            loading_constraints=loading_constraints_small)
    model2.fit()
    save_results(model2, model_output_dir, "gibbs_trig_true_k", burn_in)
    del model2
    # Gibbs sampler, misspecify k, not triangularize the loading
    model3 = bmirt.MCMC_MIRT(data, dim=large_k, num_samples=2000, alpha_prior="adaptive-ss", loading_constraints={})
    model3.fit()
    save_results(model3, model_output_dir, "gibbs_large_k", burn_in)
    del model3
    # Gibbs sampler, misspecify k, triangularize the loading
    model4 = bmirt.MCMC_MIRT(data, dim=large_k, num_samples=2000, alpha_prior="adaptive-ss",
                            loading_constraints=loading_constraints_large)
    model4.fit()
    save_results(model4, model_output_dir, "gibbs_trig_large_k", burn_in)
    del model4






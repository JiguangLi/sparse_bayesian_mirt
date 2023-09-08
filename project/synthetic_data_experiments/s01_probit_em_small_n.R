setwd("/Users/jiguangli/IBP_IRT/project")
source("probit_em_mirt.R")
source("probit_em_util.R")
library(pacman)
p_load("tidyverse", "argparser", "feather",  "magic", "hash" , "pander", "tictoc", "dplyr",  "plotmo", "caret",  "cvms", "magrittr", "arrow", "yaml", "plot.matrix", "here")


# parse arguments
arguments <- arg_parser("Test Probit on small dataset") %>% 
  add_argument(
    "--model_output_dir", 
    help= "model_output_dir", 
    default= here("models", "small_n")
  ) %>%
  add_argument(
    "--data_output_dir", 
    help= "data_output_dir", 
    default= here("data", "small_n")
  ) %>%
  parse_args()


set.seed(1)

# Create probit dataset and use the same loading of to generate logistic dataset
small_data_probit <- generate_overlap_data(n=250,
                               k = 4,
                               items_per_dim = 120,
                               overlap = 60,
                               alpha_threshold = NA,
                               loading_prior = "boolean",
                               seed= 1)
small_data_logit <- generate_logistic_data(n= 250, 
                                           loading_matrix = small_data_probit$alphas, 
                                           intercepts = small_data_probit$intercepts)
saveRDS(small_data_probit, file.path(arguments[["data_output_dir"]], "small_n_probit.rds"))
write_feather(as.data.frame(small_data_logit), file.path(arguments[["data_output_dir"]], "small_n_logit.feather"))


# Initialize Loadings
large_k <- 10
nitems <- nrow(small_data_probit$alphas)
loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.2), nitems, large_k),
                        "intercepts"= runif(nitems , -0.2,0.2), "c_params" = rep(0.5, large_k))
lambda0_path <- c(0.2, 1, 5, 10, 20, 30, 40, 50)
lambda1 <- 0.2


# PX-EM
tic()
px_em <- dynamic_posterior_exploration(data = small_data_probit$response, k = large_k, ibp_alpha = 2, mc_samples =50,
                                       ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
                                       max_iterations= 100, epsilon = 0.07, PX = TRUE, varimax = FALSE,
                                       loading_constraints= NULL, start = loading_starts_large,
                                       plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
toc()
par(mfrow=c(1,2), mar=c(2, 2, 2, 2))
plot(small_data_probit$alphas, digits=1, main= "True Loading", text.cell=list(cex=0.6), key=NULL) 
plot(px_em$lambda0_50$alphas, digits=1, main= "Estimated Loading",  text.cell=list(cex=0.6), key=NULL) 
saveRDS(px_em, file.path(arguments[["model_output_dir"]], "px_em.rds"))

# # EM
tic()
em <- dynamic_posterior_exploration(data = small_data_probit$response, k = large_k, ibp_alpha = 2, mc_samples =50,
                                       ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
                                       max_iterations= 100, epsilon = 0.07, PX = FALSE, varimax = FALSE,
                                       loading_constraints= NULL, start = loading_starts_large,
                                       plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
toc()
par(mfrow=c(1,2), mar=c(2, 2, 2, 2))
plot(small_data_probit$alphas, digits=1, main= "True Loading", text.cell=list(cex=0.6), key=NULL)
plot(em$lambda0_50$alphas, digits=1, main= "Estimated Loading",  text.cell=list(cex=0.6), key=NULL)
# save models
saveRDS(em, file.path(arguments[["model_output_dir"]], "em.rds"))


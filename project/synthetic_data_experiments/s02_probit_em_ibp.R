setwd("/Users/jiguangli/IBP_IRT/project")
source("probit_em_mirt.R")
source("probit_em_util.R")
library(pacman)
p_load("tidyverse", "argparser", "feather",  "magic", "hash" , "pander", "tictoc", "dplyr",  "plotmo", "caret",  "cvms", "magrittr", "arrow", "yaml", "plot.matrix", "here")


# parse arguments
arguments <- arg_parser("Test Probit on IBP Loading") %>% 
  add_argument(
    "--model_output_dir", 
    help= "model_output_dir", 
    default= here("models", "ibp")
  ) %>%
  add_argument(
    "--data_output_dir", 
    help= "data_output_dir", 
    default= here("data", "ibp")
  ) %>%
  parse_args()


set.seed(1)

# Create probit dataset and use the same loading of to generate logistic dataset
ibp_data_probit <- generate_ibp_data(n=250,
                                       num_items = 350,
                                       dim = 5,
                                       ibp_alpha = 2,
                                       seed= 1)
ibp_data_logit <- generate_logistic_data(n= 250, 
                                        loading_matrix = ibp_data_probit$alphas, 
                                        intercepts = ibp_data_probit$intercepts)


saveRDS(ibp_data_probit, file.path(arguments[["data_output_dir"]], "ibp_probit.rds"))
write_feather(as.data.frame(ibp_data_logit), file.path(arguments[["data_output_dir"]], "ibp_logit.feather"))


par(mfrow=c(1,2), mar=c(2, 2, 2, 2))
plot(ibp_data_probit$alphas_unsorted, digits=1, main= "Generated IBP Loading", text.cell=list(cex=0.6), key=NULL) 
plot(ibp_data_probit$alphas, digits=1, main= "Sorted_IBP_Loading",  text.cell=list(cex=0.6), key=NULL) 

# Initialize Params
large_k <- 10
nitems <- nrow(ibp_data_probit$alphas)
loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.02), nitems, large_k),
                             "intercepts"= runif(nitems , -0.02,0.02), "c_params" = rep(0.5, large_k))


lambda0_path <- c(0.5, 1, 3, 6, 10, 20, 30, 40)
lambda1 <- 0.5
tic()
px_em <- dynamic_posterior_exploration(data = ibp_data_probit$response, k = large_k, ibp_alpha = 2, mc_samples =50,
                                       ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
                                       max_iterations= 100, epsilon = 0.06, PX = TRUE, varimax = FALSE,
                                       loading_constraints= NULL, start = loading_starts_large,
                                       plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
toc()
par(mfrow=c(1,2), mar=c(2, 2, 2, 2))
plot(ibp_data_probit$alphas, digits=1, main= "True Loading", text.cell=list(cex=0.6), key=NULL)
plot(px_em$lambda0_40$alphas, digits=1, main= "Estimated Loading",  text.cell=list(cex=0.6), key=NULL)
saveRDS(px_em, file.path(arguments[["model_output_dir"]], "px_em.rds"))




# # EM
tic()
em <- dynamic_posterior_exploration(data = ibp_data_probit$response, k = large_k, ibp_alpha = 2, mc_samples =50,
                                    ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
                                    max_iterations= 100, epsilon = 0.06, PX = FALSE, varimax = FALSE,
                                    loading_constraints= NULL, start = loading_starts_large,
                                    plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
toc()
par(mfrow=c(1,2), mar=c(2, 2, 2, 2))
plot(ibp_data_probit$alphas, digits=1, main= "True Loading", text.cell=list(cex=0.6), key=NULL)
plot(em$lambda0_40$alphas, digits=1, main= "Estimated Loading" ,  text.cell=list(cex=0.6), key=NULL)
# save models
saveRDS(em, file.path(arguments[["model_output_dir"]], "em.rds"))







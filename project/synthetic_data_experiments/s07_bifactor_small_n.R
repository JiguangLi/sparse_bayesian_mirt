setwd("/Users/jiguangli/IBP_IRT/project")
source("probit_em_util.R")
library(pacman)
p_load("mirt", "tidyverse", "argparser", "feather",  "magic", "hash" , "pander", "tictoc", "dplyr",  "plotmo", "caret",  "cvms", "magrittr", "arrow", "yaml", "plot.matrix", "here")


# parse arguments
arguments <- arg_parser("Test Bifactor on small dataset") %>% 
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

# read small dataset
response_data <- arrow::read_feather(file.path(arguments[["data_output_dir"]], "small_n_logit.feather"))
data <- readRDS(file.path(arguments[["data_output_dir"]], "small_n_probit.rds"))
true_k <- 4
# Approximate Bifactor Loading
loading_list <- approximate_bifactor_loadings(data$alphas)
print(loading_list)
# Bifactor model
tic()
bifactor_model <- bfactor(response_data, loading_list)
toc()
# save results
saveRDS(bifactor_model, file.path(arguments[["model_output_dir"]], "bifactor_approx.rds"))
factor_coefs <- coef(bifactor_model)
num_items <- dim(response_data)[2]
loadings <- matrix(0, num_items, true_k+1)
for(j in 1:num_items){
  loadings[j, ] <- factor_coefs[[j]][1, 1:(true_k+1)]
}
colnames(loadings) <- c(paste0("dim",1:true_k), "intercepts")
arrow::write_feather(loadings%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "bifactor_approx_loading.feather"))


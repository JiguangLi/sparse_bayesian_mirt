setwd("/Users/jiguangli/IBP_IRT/project")
library(pacman)
p_load("tidyverse", "argparser", "feather",  "magic", "hash" , "pander", "tictoc", "dplyr",  "plotmo", "caret",  "cvms", "magrittr", "arrow", "yaml", "plot.matrix", "here", "mirt")


# parse arguments
arguments <- arg_parser("Test MHRM on small dataset") %>% 
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
data <- arrow::read_feather(file.path(arguments[["data_output_dir"]], "small_n_logit.feather"))
true_k <-4
large_k <- 10

# true model
tic()
model <- mirt(data, true_k, itemtype = "2PL",  method='MHRM', SE=TRUE ,NCYCLES=2000)
toc()
# save results
saveRDS(model, file.path(arguments[["model_output_dir"]], "mhrm_true.rds"))
factor_coefs <- coef(model)
num_items <- dim(data)[2]
loadings <- matrix(0, num_items, true_k+1)
for(j in 1:num_items){
  loadings[j, ] <- factor_coefs[[j]][1, 1:(true_k+1)]
}
colnames(loadings) <- c(paste0("dim",1:true_k), "intercepts")
arrow::write_feather(loadings%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "mhrm_true_loading.feather"))



# misspecified model
# tic()
# model_mis <- mirt(data, large_k, itemtype = "2PL",  method='MHRM', SE=TRUE ,NCYCLES=2000)
# toc()
# # save results
# saveRDS(model_mis, file.path(arguments[["model_output_dir"]], "mhrm_misspecified.rds"))
# factor_coefs <- coef(model_mis)
# loadings <- matrix(0, num_items, large_k+1)
# for(j in 1:num_items){
#   loadings[j, ] <- factor_coefs[[j]][1, 1:(large_k+1)]
# }
# colnames(loadings) <- c(paste0("dim",1:large_k), "intercepts")
# arrow::write_feather(loadings%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "mhrm_misspecified_loading.feather"))
# 
# 



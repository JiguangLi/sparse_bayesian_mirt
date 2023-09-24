setwd("/Users/jiguangli/IBP_IRT/project")
library(pacman)
p_load("mirt","tidyverse", "argparser", "feather",  "magic" , "pander", "dplyr",  "plotmo", "caret",  "cvms", "magrittr", "arrow", "here", "hash", "tictoc")

# parse arguments
arguments <- arg_parser("FIT QOL data") %>% 
  add_argument(
    "--data_input_dir", 
    help= "data_input_dir", 
    default= here("data", "qol")
  ) %>%
  add_argument(
    "--model_output_dir", 
    help= "model_output_dir", 
    default= here("models", "qol")
  ) %>%
  add_argument(
    "--dicho_threshold", 
    help= "how to dichotomize data", 
    default= 4
  ) %>%
  parse_args()

# read data and convert to binary response
lines <- readLines(file.path(arguments[["data_input_dir"]],"QOL.txt"), warn=FALSE)
cleaned_lines <- trimws(lines, "both")
data <- data.frame(cleaned_lines, stringsAsFactors = FALSE)
data <- data %>%
  tidyr::separate(col = cleaned_lines, into = paste0("col", 1:35), sep = 1:34, convert = TRUE)
binary_response <- +(data > arguments[["dicho_threshold"]]) %>% as.matrix()

# fit mhrm
loading_list <- c(1, rep(2, 4), rep(3, 4), rep(4, 6), rep(5, 6), rep(6, 5), rep(7, 5), rep(8, 4))
set.seed(1)
tic()
model <- bfactor(binary_response, loading_list)
toc()
# save results
saveRDS(model, file.path(arguments[["model_output_dir"]], "bifactor.rds"))
factor_coefs <- coef(model)
num_items <- dim(binary_response)[2]
alphas_params <- matrix(0, num_items, 10)
for(j in 1:num_items){
  alphas_params[j, ] <- factor_coefs[[j]][1, 1:10]
}
colnames(alphas_params) <- c(paste0("dim",1:9), "intercepts")
arrow::write_feather(alphas_params%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "bifactor_alphas.feather"))
loading_summary <- summary(model)
arrow::write_feather(loading_summary$rotF%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "bifactor_loadings.feather"))

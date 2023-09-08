setwd("/Users/jiguangli/IBP_IRT/project")
source("probit_em_mirt.R")
source("probit_em_util.R")
library(pacman)
p_load("tidyverse", "argparser", "feather",  "magic" , "pander", "dplyr",  "plotmo", "caret",  "cvms", "magrittr", "arrow", "here", "hash", "tictoc")

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
binary_response <- +(data > arguments[["dicho_threshold"]]) %>% as.matrix() %>% unname()

# fit px-em
lambda1 <- 0.1
lambda0_path <- c(0.1, 0.5, 1, 5, 10, 20, 30, 40)
large_k <- 10
nitems <- 35
set.seed(1)
loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.1), nitems, large_k),
                             "intercepts"= runif(nitems , -0.1,0.1), "c_params" = rep(0.5, large_k))
tic()
px_em <- dynamic_posterior_exploration(data =binary_response, k = large_k, ibp_alpha = 2, mc_samples =50,
                                       ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
                                       max_iterations= 100, epsilon = 0.04, PX = TRUE, varimax = FALSE,
                                       loading_constraints= NULL, start = loading_starts_large,
                                       plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
toc()
par(mfrow=c(1,1), mar=c(2, 2, 2, 2))
plot(px_em$lambda0_40$alphas, digits=1, main= "Estimated Loading",  text.cell=list(cex=0.6), key=NULL) 
saveRDS(px_em, file.path(arguments[["model_output_dir"]], "px_em.rds"))




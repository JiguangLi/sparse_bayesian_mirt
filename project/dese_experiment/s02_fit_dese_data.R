setwd("/Users/jiguangli/IBP_IRT/project")
source("probit_em_mirt.R")
source("probit_em_util.R")
library(pacman)
p_load("tidyverse", "argparser", "feather",  "magic" , "pander", "dplyr",  "plotmo", "caret",  "cvms", "magrittr", "arrow", "here", "hash", "tictoc")

# parse arguments
arguments <- arg_parser("FIT DESE data in 2022") %>% 
  add_argument(
    "--data_input_dir", 
    help= "data_input_dir", 
    default= here("DESE_data", "sample_processed", "2022")
  ) %>%
  add_argument(
    "--model_output_dir", 
    help= "model_output_dir", 
    default= here("models", "dese", "2022", "px_em")
  ) %>%
  add_argument(
    "--grades", 
    help= "grades_considered", 
    default= c("03", "04", "05", "06", "07", "08", "10")
  ) %>%
  add_argument(
    "--year", 
    help= "exam_year", 
    default= c("22")
  ) %>%
  add_argument(
    "--sample_size", 
    help= "number of students per grade", 
    default= 2000
  ) %>%
  parse_args()

# read data
data <- paste0(paste0("grade_", arguments[["grades"]], "_sampled"), ".feather") %>%
  set_names(., paste0("grade", arguments[["grades"]])) %>%
  map(function(.x) arrow::read_feather(file.path(arguments[["data_input_dir"]], .x)))
set.seed(1)
large_k <- 10
lambda1 <- 0.1

# grade 10 model
print("grade 10")
ir_grade10 <- data$grade10 %>% 
  select(contains("eitem")|contains("mitem")) 
print(colnames(ir_grade10))
nitems <- ncol(ir_grade10)
loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.2), nitems, large_k),
                             "intercepts"= runif(nitems , -0.2,0.2), "c_params" = rep(0.5, large_k))
lambda0_path <- c(1, 5, 10, 20, 40, 60, 80, 100)
tic()
px_em <- dynamic_posterior_exploration(data = ir_grade10 %>% as.matrix() %>% unname(), k = large_k, ibp_alpha = 2, mc_samples =50,
                                       ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
                                       max_iterations= 100, epsilon = 0.04, PX = TRUE, varimax = FALSE,
                                       loading_constraints= NULL, start = loading_starts_large,
                                       plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
toc()
par(mfrow=c(1,1), mar=c(2, 2, 2, 2))
plot(px_em$lambda0_100$alphas, digits=1, main= "Estimated Loading(Grade 10)",  text.cell=list(cex=0.6), key=NULL) 
saveRDS(px_em, file.path(arguments[["model_output_dir"]], "px_em_grade10.rds"))
rm(px_em)

# # grade 8 model
# print("grade 8")
# ir_grade8 <- data$grade08 %>% 
#   select(contains("eitem")|contains("mitem")) 
# print(colnames(ir_grade8))
# nitems <- ncol(ir_grade8)
# loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.2), nitems, large_k),
#                              "intercepts"= runif(nitems , -0.2,0.2), "c_params" = rep(0.5, large_k))
# lambda0_path <- c(0.01, 1, 5, 10, 20, 40, 60, 80, 100)
# tic()
# px_em <- dynamic_posterior_exploration(data = ir_grade8 %>% as.matrix() %>% unname(), k = large_k, ibp_alpha = 1, mc_samples =50,
#                                        ssl_lambda0_path = lambda0_path, ssl_lambda1 =0.01, pos_init =TRUE,
#                                        max_iterations= 100, epsilon = 0.04, PX = TRUE, varimax = FALSE,
#                                        loading_constraints= NULL, start = loading_starts_large,
#                                        plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
# toc()
# par(mfrow=c(1,1), mar=c(2, 2, 2, 2))
# plot(px_em$lambda0_20$alphas, digits=1, main= "Estimated Loading(Grade 8)",  text.cell=list(cex=0.6), key=NULL) 
# saveRDS(px_em, file.path(arguments[["model_output_dir"]], "px_em_grade8.rds"))
# rm(px_em)
# 
# 
# # grade 7 model
# print("grade 7")
# ir_grade7 <- data$grade07 %>% 
#   select(contains("eitem")|contains("mitem")) 
# print(colnames(ir_grade7))
# nitems <- ncol(ir_grade7)
# loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.2), nitems, large_k),
#                              "intercepts"= runif(nitems , -0.2,0.2), "c_params" = rep(0.5, large_k))
# lambda0_path <- c(1, 5, 10, 20, 40, 60, 80, 100)
# tic()
# px_em <- dynamic_posterior_exploration(data = ir_grade7 %>% as.matrix() %>% unname(), k = large_k, ibp_alpha = 2, mc_samples =50,
#                                        ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
#                                        max_iterations= 100, epsilon = 0.04, PX = TRUE, varimax = FALSE,
#                                        loading_constraints= NULL, start = loading_starts_large,
#                                        plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
# toc()
# par(mfrow=c(1,1), mar=c(2, 2, 2, 2))
# plot(px_em$lambda0_100$alphas, digits=1, main= "Estimated Loading(Grade 7)",  text.cell=list(cex=0.6), key=NULL) 
# saveRDS(px_em, file.path(arguments[["model_output_dir"]], "px_em_grade7.rds"))
# rm(px_em)
# 
# 
# # grade 6 model
# print("grade 6")
# ir_grade6 <- data$grade06 %>% 
#   select(contains("eitem")|contains("mitem")) 
# print(colnames(ir_grade6))
# nitems <- ncol(ir_grade6)
# loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.2), nitems, large_k),
#                              "intercepts"= runif(nitems , -0.2,0.2), "c_params" = rep(0.5, large_k))
# lambda0_path <- c(1, 5, 10, 20, 40, 60, 80, 100)
# tic()
# px_em <- dynamic_posterior_exploration(data = ir_grade6 %>% as.matrix() %>% unname(), k = large_k, ibp_alpha = 2, mc_samples =50,
#                                        ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
#                                        max_iterations= 100, epsilon = 0.04, PX = TRUE, varimax = FALSE,
#                                        loading_constraints= NULL, start = loading_starts_large,
#                                        plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
# toc()
# par(mfrow=c(1,1), mar=c(2, 2, 2, 2))
# plot(px_em$lambda0_100$alphas, digits=1, main= "Estimated Loading(Grade 6)",  text.cell=list(cex=0.6), key=NULL) 
# saveRDS(px_em, file.path(arguments[["model_output_dir"]], "px_em_grade6.rds"))
# rm(px_em)
# 
# 
# # grade 5 model
# print("grade 5")
# ir_grade5 <- data$grade05 %>% 
#   select(contains("eitem")|contains("mitem")) 
# print(colnames(ir_grade5))
# nitems <- ncol(ir_grade5)
# loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.2), nitems, large_k),
#                              "intercepts"= runif(nitems , -0.2,0.2), "c_params" = rep(0.5, large_k))
# lambda0_path <- c(1, 5, 10, 20, 40, 60, 80, 100)
# tic()
# px_em <- dynamic_posterior_exploration(data = ir_grade5 %>% as.matrix() %>% unname(), k = large_k, ibp_alpha = 2, mc_samples =50,
#                                        ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
#                                        max_iterations= 100, epsilon = 0.04, PX = TRUE, varimax = FALSE,
#                                        loading_constraints= NULL, start = loading_starts_large,
#                                        plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
# toc()
# par(mfrow=c(1,1), mar=c(2, 2, 2, 2))
# plot(px_em$lambda0_100$alphas, digits=1, main= "Estimated Loading(Grade 5)",  text.cell=list(cex=0.6), key=NULL) 
# saveRDS(px_em, file.path(arguments[["model_output_dir"]], "px_em_grade5.rds"))
# rm(px_em)
# 
# # grade 4 model
# print("grade 4")
# ir_grade4 <- data$grade04 %>% 
#   select(contains("eitem")|contains("mitem")) 
# print(colnames(ir_grade4))
# nitems <- ncol(ir_grade4)
# loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.2), nitems, large_k),
#                              "intercepts"= runif(nitems , -0.2,0.2), "c_params" = rep(0.5, large_k))
# lambda0_path <- c(1, 5, 10, 20, 40, 60, 80, 100)
# tic()
# px_em <- dynamic_posterior_exploration(data = ir_grade4 %>% as.matrix() %>% unname(), k = large_k, ibp_alpha = 2, mc_samples =50,
#                                        ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
#                                        max_iterations= 100, epsilon = 0.04, PX = TRUE, varimax = FALSE,
#                                        loading_constraints= NULL, start = loading_starts_large,
#                                        plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
# toc()
# par(mfrow=c(1,1), mar=c(2, 2, 2, 2))
# plot(px_em$lambda0_100$alphas, digits=1, main= "Estimated Loading(Grade 4)",  text.cell=list(cex=0.6), key=NULL) 
# saveRDS(px_em, file.path(arguments[["model_output_dir"]], "px_em_grade4.rds"))
# rm(px_em)
# 
# # grade 3 model
# print("grade 3")
# ir_grade3 <- data$grade03 %>% 
#   select(contains("eitem")|contains("mitem")) 
# print(colnames(ir_grade3))
# nitems <- ncol(ir_grade3)
# loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.2), nitems, large_k),
#                              "intercepts"= runif(nitems , -0.2,0.2), "c_params" = rep(0.5, large_k))
# lambda0_path <- c(1, 5, 10, 20, 40, 60, 80, 100)
# tic()
# px_em <- dynamic_posterior_exploration(data = ir_grade3 %>% as.matrix() %>% unname(), k = large_k, ibp_alpha = 2, mc_samples =50,
#                                        ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
#                                        max_iterations= 100, epsilon = 0.04, PX = TRUE, varimax = FALSE,
#                                        loading_constraints= NULL, start = loading_starts_large,
#                                        plot=FALSE, stop_rotation=100, random_state = 1, cores=8)
# toc()
# par(mfrow=c(1,1), mar=c(2, 2, 2, 2))
# plot(px_em$lambda0_100$alphas, digits=1, main= "Estimated Loading(Grade 3)",  text.cell=list(cex=0.6), key=NULL) 
# saveRDS(px_em, file.path(arguments[["model_output_dir"]], "px_em_grade3.rds"))
# rm(px_em)
# 
# 
# 
# 
# 
# 
# 

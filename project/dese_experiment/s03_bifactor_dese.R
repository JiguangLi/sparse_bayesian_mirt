setwd("/Users/jiguangli/IBP_IRT/project")
library(pacman)
p_load("mirt", "tidyverse", "argparser", "feather",  "magic", "hash" , "pander", "tictoc", "dplyr",  "plotmo", "caret",  "cvms", "magrittr", "arrow", "yaml", "plot.matrix", "here")

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
    default= here("models", "dese", "2022", "bifactor")
  ) %>%
  add_argument(
    "--grades", 
    help= "grades_considered", 
    default= c("10")
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

# grade 10, 2-factor model, main component is english, secondary component is math
ir_grade10 <- data$grade10 %>% 
  select(contains("eitem")|contains("mitem")) 
nitems <- ncol(ir_grade10)
n_eitems <- ir_grade10 %>% select(contains("eitem")) %>% ncol()
n_mitems <- nitems - n_eitems
loading_list <- c(rep(NA, n_eitems), rep(1, n_mitems))
tic()
bifactor_model <- bfactor(ir_grade10 %>% as.matrix(), loading_list)
toc()
# save results
saveRDS(bifactor_model, file.path(arguments[["model_output_dir"]], "bifactor_grade10.rds"))
factor_coefs <- coef(bifactor_model)
num_items <- dim(ir_grade10)[2]
alphas_params <- matrix(0, num_items, 3)
for(j in 1:num_items){
  alphas_params[j, ] <- factor_coefs[[j]][1, 1:3]
}
colnames(alphas_params) <- c(paste0("dim",1:2), "intercepts")
arrow::write_feather(alphas_params%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "bifactor_grade10_alphas.feather"))
loading_summary <- summary(bifactor_model)
arrow::write_feather(loading_summary$rotF%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "bifactor_grade10_loadings.feather"))

# grade 10, 3-factor model, main component is general ability, two components: math and english items
loading_list <- c(rep(1, n_eitems), rep(2, n_mitems))
tic()
bifactor_model_3factor <- bfactor(ir_grade10 %>% as.matrix(), loading_list)
toc()
# save results
saveRDS(bifactor_model_3factor, file.path(arguments[["model_output_dir"]], "bifactor_grade10_3factor.rds"))
factor_coefs <- coef(bifactor_model_3factor)
alphas_params <- matrix(0, num_items, 4)
for(j in 1:num_items){
  alphas_params[j, ] <- factor_coefs[[j]][1, 1:4]
}
colnames(alphas_params) <- c(paste0("dim",1:3), "intercepts")
arrow::write_feather(alphas_params%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "bifactor_grade10_alphas_3factor.feather"))
loading_summary <- summary(bifactor_model_3factor)
arrow::write_feather(loading_summary$rotF%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "bifactor_grade10_loadings_3factor.feather"))


# grade 10, 6-factor model, main component is general ability, one secondary english ability, 4 secondary math abilities 
# 2: Algebra; 3: Geometry; 4: Number and Quantity; 5: statistics
loading_list <- c(rep(1, n_eitems), 2,2,3,2,2,3,2,3,2,4, 3, 2, 2, 3, 2,3,3,4,5, 2, 3, 3, 2, 3, 4, 5, 3, 5, 3, 2, 3, 2)
tic()
bifactor_model_6factor <- bfactor(ir_grade10 %>% as.matrix(), loading_list)
toc()
# save results
saveRDS(bifactor_model_6factor, file.path(arguments[["model_output_dir"]], "bifactor_grade10_6factor.rds"))
factor_coefs <- coef(bifactor_model_6factor)
alphas_params <- matrix(0, num_items, 7)
for(j in 1:num_items){
  alphas_params[j, ] <- factor_coefs[[j]][1, 1:7]
}
colnames(alphas_params) <- c(paste0("dim",1:6), "intercepts")
arrow::write_feather(alphas_params%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "bifactor_grade10_alphas_6factor.feather"))
loading_summary <- summary(bifactor_model_6factor)
arrow::write_feather(loading_summary$rotF%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "bifactor_grade10_loadings_6factor.feather"))



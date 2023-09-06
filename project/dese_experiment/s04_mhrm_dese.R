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
    default= here("models", "dese", "2022", "mhrm")
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

# grade 10
ir_grade10 <- data$grade10 %>% 
  select(contains("eitem")|contains("mitem")) 

set.seed(1)

tic()
model <- mirt(ir_grade10, 10, itemtype = "2PL",  method='MHRM', SE=TRUE ,NCYCLES=2000)
toc()
# save results
saveRDS(model, file.path(arguments[["model_output_dir"]], "mhrm_grade10.rds"))
factor_coefs <- coef(model)
num_items <- dim(ir_grade10)[2]
loadings <- matrix(0, num_items, 11)
for(j in 1:num_items){
  loadings[j, ] <- factor_coefs[[j]][1, 1:11]
}
colnames(loadings) <- c(paste0("dim",1:10), "intercepts")
arrow::write_feather(loadings%>% as.data.frame(), file.path(arguments[["model_output_dir"]], "mhrm_grade10_loading.feather"))



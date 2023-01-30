library(pacman)
p_load("tidyverse", "argparser", "feather",  "yaml", "mirt", "sparsepca")


config_filepath <- "/Users/jiguangli/IBP_IRT/config/data_config.yaml"
config <- read_yaml(config_filepath)
# read real data
data_filename_format <- paste(config[["year"]], config[["grade"]], config[["subject"]], "%s.feather", sep= "_")
data_filename_format <- paste(config[["output_data_dir"]], data_filename_format, sep= "/")
train_df <- read_feather(sprintf(data_filename_format, "train"))
test_df <- read_feather(sprintf(data_filename_format, "test"))

# read fake data
#train_df <- arrow::read_feather(config[["fake_data3"]])
items_col <- train_df %>% dplyr::select(contains("item")) %>% names()
train_df <- train_df %>% dplyr::select(items_col)

# PCA on covariance mat (or correlation mat?)
# cov_mat <- cov(train_df %>% select(contains("item")))
# eigen_cov <- eigen(cov_mat)

# Train Baseline IRT model with PCA factor loadings
# train_single_cmirt <- function(df, num_factor, loading){
#   if (num_factor == 1){
#         return(mirt(df, 1, itemtype = "2PL", SE=TRUE))
#     }else{
#         num_items <- dim(df)[2]
#         factor_format <- paste0("F", "%s = 1-", num_items)
#         factor_str <- paste(1:num_factor %>% 
#                               map(function(x) sprintf(factor_format, x)), collapse="\n")
#         start_arg <- vector(length = num_factor)
#         for(i in 1: num_factor){
#             start_format <- paste0("(%s,", "a", i, ",", "%s)")
#             start_arg[[i]]<- paste(map2(1:num_items, loading[, i], function(x,y) sprintf(start_format, x, y)), collapse=",")
#         }
#         start_arg <- paste0("START = ",paste(start_arg, collapse=","))
#         
#         end_arg <- vector(length = num_factor)
#         for(i in 1: num_factor){
#           end_format <- paste0("(", "%s", ",a", i,")")
#           end_arg[[i]]<- paste(map(1:num_items, function(x) sprintf(end_format, x)), collapse=",")
#         }
#         end_arg <- paste0("FIXED = ",paste(end_arg, collapse=","))
#         
#         fs <- paste(factor_str, start_arg, end_arg, sep="\n")
#         return(mirt(df, fs , itemtype= "2PL", SE=TRUE, methods="QMCEM"))
#       }
#     
# }
# 
# train_all_cmirt <- function(df, num_factors, eigenvecs){
#   models <- num_factors %>%
#     set_names (. , paste(., "factor(s)" ,sep="_")) %>%
#     map(function(.x) train_single_cmirt(df, .x, eigenvecs[, 1:{.x}]))
#   return(models)
# 
# }
# num_factors <- 1:4
# cmirt_models <- train_all_cmirt(train_df %>% select(contains("item")), num_factors, eigen_cov$vectors)

# Train Exploratory MIRT models
train_single_cmirt <- function(df, num_factor){
  if(num_factor <= 2){
    return(mirt(df, num_factor, itemtype= "2PL", SE=TRUE, NCYCLES=500))
  }else{
    return(mirt(df, num_factor, itemtype= "2PL", SE=TRUE, methods="MHRM", technical = list(NCYCLES = 1800)))
    
  }
}

train_all_cmirt <- function(df, num_factors){
  models <- num_factors %>% 
    set_names (. , paste(., "factors" ,sep="_")) %>%
    map(function(.x) train_single_cmirt(df, .x))
  return(models)
  
}

num_factors <- 1:2
cmirt_models <- train_all_cmirt(train_df, num_factors)


# Estimate latent traits
estimate_latent_params <- function(df, model){
  num_factor <- model@Model$nfact
  method <- if_else(num_factor > 2, "MAP", "EAP")
  return(fscores(model, method = method, response.pattern= df %>% as.matrix(), full.scores = TRUE, rotate = "oblimin") %>% as.data.frame())
  
}

get_latent_traits <- function(cmirt_models, df){
  latent_traits <- cmirt_models %>%
    map(function(.x) estimate_latent_params(df, .x))
  
  return(latent_traits)
  
}

factor_coefs <- cmirt_models %>%
  map(coef, rotate = "oblimin")

latent_traits_train <- get_latent_traits(cmirt_models, train_df %>% select(contains("item"))) 


# Get Predictions

inverse_logit <- function (x){
  return(1/(1+exp(-x)))
}

predict_single_question <- function(coef_ests, latent_traits, num_factors){
  slopes <- coef_ests[1, 1:num_factors] %>% as.numeric()
  intercept <- coef_ests[1, (num_factors+1)]%>% as.numeric()
  
  col_predictions <- latent_traits %>% 
    as.matrix() %>%
    magrittr::multiply_by_matrix (diag(slopes, num_factors, num_factors)) %>%
    rowSums() %>%
    magrittr::add(intercept) %>%
    inverse_logit()
  
  return(col_predictions)
}


construct_prediction_df <- function(coef_est, latent_traits, q_names){
  
  num_factors <- dim(latent_traits)[2]
  predict_df <- q_names %>%
    set_names(.,.) %>%
    map(function(.x) predict_single_question(coef_est[[.x]], latent_traits, num_factors)) %>%
    as.data.frame()
  
  return(predict_df)
  
}

construct_prediction_df_all <- function(models, factor_coefs,  latent_traits, q_names){
  result_dfs <- names(models) %>%
    set_names(.,.) %>%
    map(function(.x) construct_prediction_df(factor_coefs[[.x]], latent_traits[[.x]] %>% select(starts_with("F")), q_names) )
  return(result_dfs)
}

train_predictions <- construct_prediction_df_all(cmirt_models, factor_coefs, latent_traits_train,items_col)



# Evaluation
Y <- train_df %>% select(contains("item")) %>% as.matrix()
train_size <- prod(dim(train_df %>% select(contains("item"))))
irt_mse_insample <- names(cmirt_models) %>%
  set_names(.,.) %>%
  map(function(.x) sum((as.matrix(train_df %>% select(contains("item"))) - as.matrix(train_predictions[[.x]]))^2)/train_size)

irt_acc <- names(cmirt_models) %>%
  set_names(.,.) %>%
  map(function(.x) 1- sum(abs(as.matrix(train_df %>% select(contains("item"))) - round(train_predictions[[.x]])))/train_size)

irt_loglik <- names(cmirt_models) %>%
  set_names(.,.) %>%
  map(function(.x) cmirt_models[[.x]]@Fit$logLik)

irt_myloglik <- names(cmirt_models) %>%
  set_names(.,.) %>%
  map(function(.x) (sum(Y*log(train_predictions[[.x]])) + sum((1-Y)*log(1-train_predictions[[.x]]))))


# compute_loglikelihood <- function(pred, true_data){
#   term1<- 
#   
# }

# save_cmirt_models <- function(cmirt_models, outpur_dir){
#   names(cmirt_models) %>%
#     paste( . , config[["year"]], config[["grade"]], config[["subject"]], "cmirt.rds", sep= "_") %>%
#     map2(cmirt_models, function(.x, .y) saveRDS(.y, file.path(outpur_dir, .x)))
# }
# 
# save_cmirt_models(cmirt_models, config[["model_output_dir"]])


create_alpha_matrix <-function(coef_est, q_num, f_num){
  alpha_mat <- matrix(0, q_num, f_num)
  for(i in 1:q_num){
    alpha_mat[i, ] <- coef_est[[i]][1:f_num]
  }
  return(alpha_mat)
}

irt_alphas <- names(cmirt_models) %>%
  set_names(.,.) %>%
  map(function(.x) create_alpha_matrix(factor_coefs[[.x]], length(items_col), cmirt_models[[.x]]@Model$nfact ))


#save(cmirt_models, train_predictions, irt_alphas, irt_acc, irt_loglik, irt_mse_insample, irt_myloglik, file = file.path(config[["model_output_dir"]], "fake_dataset3_irt.Rdata"))
#write.table(irt_alphas[[3]],file="/Users/jiguangli/IBP_IRT/markdown/dataset3_3factor_alphas.txt")




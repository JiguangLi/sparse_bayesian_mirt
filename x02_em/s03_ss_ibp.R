source("FACTOR_CODE_update.R")
library(pacman)
p_load("tidyverse", "argparser", "feather",  "yaml", "mirt", "sparsepca", 
       "mvtnorm", "partitions", "nloptr", "glmnet")

# read real data
config_filepath <- "/Users/jiguangli/IBP_IRT/config/data_config.yaml"
config <- read_yaml(config_filepath)
data_filename_format <- paste(config[["year"]], config[["grade"]], config[["subject"]], "%s.feather", sep= "_")
data_filename_format <- paste(config[["output_data_dir"]], data_filename_format, sep= "/")
train_df <- read_feather(sprintf(data_filename_format, "train"))
test_df <- read_feather(sprintf(data_filename_format, "test"))

# read fake data
#train_df <- arrow::read_feather(config[["fake_data"]])


# INITIALIZATIONS and Center
items_col <- train_df %>% dplyr::select(contains("item")) %>% names()
Y <- train_df %>% dplyr::select(all_of(items_col)) %>% as.matrix()
diff_terms <- colMeans(Y)
Y_center <- sweep(Y, 2, diff_terms)

K <- 4
G <-length(items_col)
startB<-matrix(rnorm(G*K),G,K)
alpha<-1/G
lambda1<-0.001
epsilon<-0.05
myImagePlot(abs(startB),F)
title("Initialization")

# Dumb Run (Approximate M-step)
# q1: why do we need to specify theta=0.5, isn't it determined by alpha
# q2: why p=10 != k in the example
start<-list(B=startB,sigma=rep(1,10),theta=rep(0.5,K))

lambda0<-1
result_1<-FACTOR_ROTATE(Y_center,lambda0,lambda1,start,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-2
result_2<-FACTOR_ROTATE(Y_center,lambda0,lambda1,result_1,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-3
result_3<-FACTOR_ROTATE(Y_center,lambda0,lambda1,result_2,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-4
result_4<-FACTOR_ROTATE(Y_center,lambda0,lambda1,result_3,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-5
result_5<-FACTOR_ROTATE(Y_center,lambda0,lambda1,result_4,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-10
result_10<-FACTOR_ROTATE(Y_center,lambda0,lambda1,result_5,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-20
result_20<-FACTOR_ROTATE(Y_center,lambda0,lambda1,result_10,K,epsilon,alpha,TRUE,TRUE,100,TRUE)
 
lambda0<-30
result_30<-FACTOR_ROTATE(Y_center,lambda0,lambda1,result_20,K,epsilon,alpha,TRUE,TRUE,100,TRUE)


lambda0 <- c("1","2" ,"3", "4","5","10","20","30")

predictions <- list(result_1, result_2, result_3, result_4, result_5, result_10, result_20, result_30) %>%
  set_names(., lambda0) %>%
  map(function(.x) sweep(.x$Omega %*% t(.x$B),2, diff_terms, "+" ))


train_size <- prod(dim(Y))

mse_insample <- names(predictions) %>%
  set_names(.,.) %>%
  map(function(.x) sum((Y - predictions[[.x]])^2)/train_size)

acc_insample <- names(predictions) %>%
  set_names(.,.) %>% 
  map(function(.x) 1- sum(abs(Y - round(predictions[[.x]])))/train_size)

loglik_insample <- names(predictions) %>%
  set_names(.,.) %>%
  map(function(.x) sum(Y*log(predictions[[.x]])) + sum((1-Y)*log(1-predictions[[.x]])))


#write.table(result_30$B,file="/Users/jiguangli/IBP_IRT/markdown/dataset1_ibp_B.txt")

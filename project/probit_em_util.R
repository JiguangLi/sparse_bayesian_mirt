library(pacman)
p_load("tidyverse", "argparser", "feather",  "yaml", "here", "hash" , "parallel", "dplyr")

# Make Probablistic Predictions of Item Response Data after Estimation
em_probit_predict <-function(data, params){
  m <-params[["mc_samples"]]
  predict_all <- pnorm(t(params[["alphas"]] %*% t(params[["thetas"]]) + params[["intercepts"]]))
  avg_pred <- rowsum(predict_all, rep(1:(nrow(predict_all)/m), each=m))/m
  return(avg_pred)
}

# Compute log likelihood
loglikhood <-function(data, params){
  m <-params[["mc_samples"]]
  linear_term <- t(params[["alphas"]] %*% t(params[["thetas"]]) + params[["intercepts"]])
  ytilde <- data[rep(seq_len(nrow(data)), each = m), ]
  result <- sum(log(pnorm((2*ytilde-1)*linear_term)))/m
  return(result)
}

# Evaluate Probit Models to compute predictions, prediction MSE, mean absolute error, logliklihood
em_probit_eval <-function(data, params){
  pred <- em_probit_predict(data,params)
  mse <- mean((data-pred)^2)
  mae <- mean(abs(data-pred))
  likhood <- loglikhood(data, params)
  result <- list(pred=pred, mse=mse, mae=mae, loglik= likhood)
  return(result)
}

# Output lower triangularized loading constraints for identification purpose
lower_trig_loading <- function(k, trunc_lev){
  if(k==trunc_lev){
    num_constraints <- (k+1)*k/2
  }else{
    num_constraints <- ((2*trunc_lev) - k +1)*k/2
  }
  result <- matrix(NA, num_constraints, 3)
  row_count <- 1
  for (l in 1:k){
    for (m in l:trunc_lev){
        result[row_count, 1]=l
        result[row_count, 2] = m
        if(l==m){
          result[row_count, 3] <- 1
        }else{
          result[row_count, 3] <- 0  
        }
      row_count <- row_count +1
    }
  }
  return(result)
  
}
  
  
# generate overlap data (Section 2 of the Paper)
# Args:
#  n: number of observations
#  k: true dimension
#  item_per_dim: how many items per k
#  overlap: how many items are overlapping between dimension k and k+1
#  loading_prior is boolean or generated from normal
#  alpha_thereshold: any abs(loading) smaller than alpha_threshold would be set 0. Irrelavant when loading prior is boolean
generate_overlap_data <- function(n, k, items_per_dim, overlap, alpha_threshold= NA, loading_prior="boolean", loading_constraints=NULL, seed=1){
  set.seed(seed)
  m <- items_per_dim + (items_per_dim-overlap)*(k-1) # num of items
  gamma_mat <- matrix(0, m, k)
  for(i in 1:k){
    row_start <- 1+ (i-1)*(items_per_dim-overlap)
    gamma_mat[(row_start:(row_start+items_per_dim-1)), i] <- 1
  }
  thetas <- MASS::mvrnorm(n=n, mu = rep(0, k), diag(k)) # n by k
  intercepts <- runif(m, -1.5, 1.5) # 1 by m
  if(loading_prior=="boolean"){
    alphas <- gamma_mat
  }else{
    alphas_temp <-  MASS::mvrnorm(n=m, mu = rep(1, k), diag(k))
    # make abs(alphas)> thereshold 0
    alphas_temp <- alphas_temp * (abs(alphas_temp) > alpha_threshold)
    # fill these zeros with small values larger than the threshold
    alphas_temp[alphas_temp==0] <- runif(length(alphas_temp[alphas_temp ==0]), alpha_threshold, alpha_threshold+0.5)
    alphas <- alphas_temp * gamma_mat
  }
  if(is.null(loading_constraints)==FALSE){
    for(i in 1:nrow(loading_constraints)){
      alphas[loading_constraints[i, 1], loading_constraints[i, 2]] <- loading_constraints[i, 3]
    }
  }
  prob_mat <- pnorm(sweep(thetas %*% t(alphas), 2, intercepts, "+"))
  response <- +(matrix(runif(m*n, 0, 1), n, m) < prob_mat)
  result <- hash("response"= response, "alphas"= alphas, "intercepts"= intercepts, "thetas" = thetas, "gammas" = gamma_mat, "items_per_dim" = items_per_dim, "overlap"= overlap)
  return(result)
}

# generate Bifactor Data
# Args:
#  n: number of observations
#  k: num of dimensions including primart factor
#  item_per_dim: how many items per k
#  loading_prior is boolean or generated from normal
#  alpha_thereshold: any abs(loading) smaller than alpha_threshold would be set 0. Irrelavant when loading prior is boolean
generate_bifactor_data <- function(n, k, items_per_dim, alpha_threshold=0.2, loading_prior="boolean", loading_constraints=NULL, seed=1){
  set.seed(seed)
  m <- items_per_dim *(k-1) # num of items
  gamma_mat <- matrix(0, m, k)
  print(dim(gamma_mat))
  gamma_mat[, 1] <- 1
  for(i in 2:k){
    row_start <- 1+ (i-2)*(items_per_dim)
    print(row_start)
    gamma_mat[(row_start:(row_start+items_per_dim-1)), i] <- 1
  }
  thetas <- MASS::mvrnorm(n=n, mu = rep(0, k), diag(k)) # n by k
  intercepts <- runif(m, -1.5, 1.5) # 1 by m
  if(loading_prior=="boolean"){
    alphas <- gamma_mat
  }else{
    alphas_temp <-  MASS::mvrnorm(n=m, mu = rep(1, k), diag(k))
    # make abs(alphas)> thereshold 0
    alphas_temp <- alphas_temp * (abs(alphas_temp) > alpha_threshold)
    # fill these zeros with small values larger than the threshold
    alphas_temp[alphas_temp==0] <- runif(length(alphas_temp[alphas_temp ==0]), alpha_threshold, alpha_threshold+0.5)
    alphas <- alphas_temp * gamma_mat
  }
  if(is.null(loading_constraints)==FALSE){
    for(i in 1:nrow(loading_constraints)){
      alphas[loading_constraints[i, 1], loading_constraints[i, 2]] <- loading_constraints[i, 3]
    }
  }
  prob_mat <- pnorm(sweep(thetas %*% t(alphas), 2, intercepts, "+"))
  response <- +(matrix(runif(m*n, 0, 1), n, m) < prob_mat)
  result <- hash("response"= response, "alphas"= alphas, "intercepts"= intercepts, "thetas" = thetas, "gammas" = gamma_mat, "items_per_dim" = items_per_dim)
  return(result)
}
  
  
# generate data from logistic model given loading and intercepts
# Args:
#   loading_matrix: m by k
#   intercepts: m-vector
generate_logistic_data <- function(n, loading_matrix, intercepts){
  set.seed(42)
  m <- dim(loading_matrix)[1]
  k <- dim(loading_matrix)[2]
  factors <- MASS::mvrnorm(n=n, mu = rep(0, k), diag(k)) # n by k
  linear_term <- sweep(factors %*% t(loading_matrix), 2, intercepts, "+")
  prob_mat <- 1/(1+exp(-linear_term))
  response <- +(matrix(runif(m*n, 0, 1), n, m) < prob_mat)
  return(response)
}

# Generate data with IBP loading (Section 6.2 of the paper)
generate_ibp_data <- function(n, num_items, dim, ibp_alpha, seed=1){
  set.seed(seed)
  v_s <- rbeta(dim, ibp_alpha, 1, ncp = 0)
  c_s <- cumprod(v_s)
  alphas <- matrix(0, num_items, dim)
  for(i in 1:dim){
    alphas[, i] <- rbinom(num_items, 1, c_s[dim-i+1])
  }
  # sort alphas by rows
  integer_vec <- numeric(num_items)
  for(j in 1:num_items){
    integer_vec[j] <- as.integer(paste0(rev(alphas[j, ]), collapse = ""))
  }
  ordering_index <- order(integer_vec)
  sorted_alphas <- matrix(0, num_items, dim)
  for(l in 1:num_items){
    sorted_alphas[l , ] <-alphas[ordering_index[l] ,] 
  }
  # generate responses
  thetas <- MASS::mvrnorm(n=n, mu = rep(0, dim), diag(dim)) # n by k
  intercepts <- runif(num_items, -1.5, 1.5) # 1 by m
  prob_mat <- pnorm(sweep(thetas %*% t(sorted_alphas), 2, intercepts, "+"))
  response <- +(matrix(runif(num_items*n, 0, 1), n, num_items) < prob_mat)
  result <- hash("response"= response, "alphas"= sorted_alphas, "intercepts"= intercepts, "thetas" = thetas, "alphas_unsorted" = alphas)
  return(result)
}


# Approximate Bifactor Loadings when Models are misspecified
# Assuming The First Column is the Main Factor
approximate_bifactor_loadings <- function(true_loading){
  num_items <- dim(true_loading)[1]
  num_dim <- dim(true_loading)[2]
  result_vec <- numeric(num_items)
  for(i in 1:num_items){
      if(max(true_loading[i, 2:num_dim]) == 0){
        result_vec[i] <- NA
      }else{
        result_vec[i]<- which.max(abs(true_loading[i, 2:num_dim]))
      }
  }
  return(result_vec)
}

# Reorder Columns of Estimated Loading to Match the Truth based on MSE. We also allow switching the sign for the ENTIRE latent dimension.
reorder_est_loadings <- function(ref_mat, est_mat){
  num_items <- dim(ref_mat)[1]
  num_dim <- dim(ref_mat)[2]
  est_dim <- dim(est_mat)[2]
  reorder_mat <- matrix(0, num_items, num_dim)
  new_order <- numeric(num_dim)
  reorder_whole <- matrix(0 , num_items, est_dim)
  for(i in 1:num_dim){
    set_difference <- setdiff(1:est_dim, new_order)
    mse_temp <- numeric(length(set_difference))
    sign_temp <-  numeric(length(set_difference))
    for(j in 1:length(set_difference)){
      mse_pos <- mean((ref_mat[, i]- est_mat[, set_difference[j]])^2)
      mse_neg <-  mean((ref_mat[, i]+ est_mat[, set_difference[j]])^2)
      if(mse_pos < mse_neg){
        mse_temp[j] <- mse_pos
        sign_temp[j] <- 1
      }else{
        mse_temp[j] <- mse_neg
        sign_temp[j] <- 0
      }
      
    }
    temp_idx <- which.min(mse_temp)
    new_order[i] <- set_difference[temp_idx]
    if(sign_temp[temp_idx] == 1){
      reorder_mat[, i] <- est_mat[, new_order[i]]
      reorder_whole[, i] <- est_mat[, new_order[i]]
    }else{
      reorder_mat[, i] <- -est_mat[, new_order[i]]
      reorder_whole[, i] <- -est_mat[, new_order[i]]
    }
    
  }
  if(est_dim > num_dim){
    reorder_whole[, (num_dim+1):est_dim] <- est_mat[, setdiff(1:est_dim, new_order)]  
  }
  
  return(list("new_order"= new_order, "reorder_mat" = reorder_mat, "reorder_whole" = reorder_whole))
}


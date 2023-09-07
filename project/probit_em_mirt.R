library(pacman)
p_load("tidyverse", "dplyr", "glmnet", "TruncatedNormal", "progress", "hash", "nloptr", "furrr", "purrr", "plot.matrix")


# Initialize Loading Matrix: m by k
initialize_alphas <- function(m, k, loading_constraints, start){
  if (is.null(start)==FALSE){
      if(has.key("alphas",start)){
        alphas <- start[["alphas"]]
      }else{
        alphas <- matrix(runif(m*k, 0, 0.02), m, k)
      }
  }else{
        alphas <- matrix(runif(m*k, 0, 0.02), m, k)
  }
  if(is.null(loading_constraints)==FALSE){
    for(i in 1:nrow(loading_constraints)){
      alphas[loading_constraints[i, 1], loading_constraints[i, 2]] <- loading_constraints[i, 3]
    }
  }
  return(alphas)
}

# Initialize Factors thetas: size (n*mc_samples, k)
initialize_thetas <- function(n, k , mc_samples, start){
  if (is.null(start)==FALSE){
      if(has.key("thetas",start)){
        return(start[["thetas"]]) 
      }
  }
  thetas<- matrix(0, n*mc_samples, k)
  return(thetas)
}

# Initialize Gammas and C-parameters:
# gammas: m by k, c params: k dimensional vector
initialize_gammas_and_cs <- function(m, k, ibp_alpha, start){
  v_s <- rbeta(k, ibp_alpha, 1, ncp = 0)
  c_s <- cumprod(v_s)
  gammas <- matrix(0, m, k)
  for(i in 1:k){
    gammas[, i] <- rbinom(m, 1, c_s[i])
  }
  cg_result <- hash("c_params"= c_s, "gammas" = gammas)
  if (is.null(start)==FALSE){
    if(has.key("gammas",start)){
      cg_result[["gammas"]] <- start[["gammas"]]
    }
    if(has.key("c_params", start)){
      cg_result[["c_params"]] <- start[["c_params"]]
    }
  }
  return(cg_result)
}

# Initialize Intercepts: m-dimensional vector
initialize_intercepts <-function(m, start){
  if (is.null(start)==FALSE){
    if(has.key("intercepts",start)){
      return(start[["intercepts"]])
    }
  }
  return(runif(m , -0.2,0.2))
  
}

# Initialize all model parameters and save them in a hash object
initialize_params <-function(n, m, k, ibp_alpha, mc_samples, loading_constraints, start, 
                            ssl_lambda0, ssl_lambda1){
  result <- hash("n"=n, "m"=m, "k"=k, "ibp_alpha" = ibp_alpha , "mc_samples"= mc_samples,"lambda0"= ssl_lambda0, "lambda1" = ssl_lambda1)
  result[["alphas"]] <- initialize_alphas(m, k, loading_constraints, start)
  result[["thetas"]] <- initialize_thetas(n, k, mc_samples, start)
  result[["alphas_prev"]] <- result[["alphas"]]
  gamma_c <- initialize_gammas_and_cs(m, k, ibp_alpha, start)
  result[["gammas"]] <- gamma_c[["gammas"]]
  result[["c_params"]] <- gamma_c[["c_params"]]
  result[["intercepts"]] <- initialize_intercepts(m, start)
  result[["intercepts_prev"]] <- result[["intercepts"]]
  return(result)
}

# Obtain mc_samples of factors for each observations:
#   Theta_i has size(mc_samples , k) 
sample_single_thetas <- function(m, k, mc_samples, d1, d2, s_arr){
  inv_term <- solve(d1%*%t(d1)+diag(m))
  cov_vo <- diag(k) - t(d1)%*%inv_term %*% d1
  cov_v1 <- diag(1/s_arr) %*% (d1%*%t(d1) + diag(m)) %*% diag(1/s_arr)
  trunc_lev <- -diag(1/s_arr) %*% d2
  v0_s <- MASS::mvrnorm(n=mc_samples, mu = rep(0, k), cov_vo)
  v1_s <- rtmvnorm(n=mc_samples, mu = rep(0, m), sigma= cov_v1, lb = trunc_lev )
  linear_term <- t(d1)%*%inv_term%*%diag(s_arr)
  theta_i <- v0_s + t(linear_term %*% t(v1_s))
  return(theta_i)
}

# E-step to obtain mc_samples of factors for each observation
sample_thetas <- function(data, params){
  n <- params[["n"]]
  m <- params[["m"]]
  k <- params[["k"]]
  mc_samples <-params[["mc_samples"]]
  d1_array <- sweep((rep(1, n) %x% params[["alphas"]]), 1, c(t(2*data-1)), "*")
  d2_array <- sweep(2*data-1, 2, params[["intercepts"]], "*")
  s_array <- (rowSums(params[["alphas"]]^2)+1)^0.5
  thetas <- 1:n %>% 
    future_map(function(.x) sample_single_thetas(m, k, mc_samples, d1_array[((.x-1)*m+1):(.x*m), ], d2_array[.x, ], s_array), .options = furrr_options(seed = TRUE)) %>% 
    do.call(rbind, .)
  return(thetas)
}

# E Step to updata Gamma (variable selection) variables
update_gammas <- function(params){
  lambda0 <- params[["lambda0"]]
  lambda1 <- params[["lambda1"]]
  lambda1_part <- (lambda1/2 * exp(-lambda1* abs(params[["alphas"]])))%>%
    sweep(., 2, params[["c_params"]], "*")
  lambda0_part <- (lambda0/2 * exp(-lambda0* abs(params[["alphas"]])))%>%
    sweep(., 2, 1-params[["c_params"]], "*")
  gammas_temp <- lambda1_part/(lambda1_part+lambda0_part)
  return(gammas_temp)
}

# E-step
e_step <- function(data, params){
  params[["thetas"]] <- sample_thetas(data, params)
  params[["gammas"]] <- update_gammas(params)
  return(params)
}

# M-step to run a penalized probit regression on idx-th item
train_single_probit<- function(idx, y, thetas, gammas, lambda0, lambda1, loading_constraints){
  k <- dim(thetas)[2]
  n <- length(y)
  mc_samples <- dim(thetas)[1]/n
  free_mask <- rep(1, k)
  y_tilde <- rep(y, each=mc_samples)
  # check if all the latent dimensions are free to estimate
  if(is.null(loading_constraints)){
    free_dim <-k
  }else{
    temp_mat <- loading_constraints[loading_constraints[,1] == idx,]
    # no constriant for item idx
    if(all(is.na(temp_mat))){
      free_dim <- k
    }else{ 
      # one constraint
      if(is.vector(temp_mat)){
        free_dim <- k-1 # bug spotted: previously, free_dim <- 1
        fixed_pos <- temp_mat[2]  
        fixed_vals <- temp_mat[3]
      }else{ # more than one constraint
        free_dim <- k-nrow(temp_mat)
        fixed_pos <- temp_mat[, 2]  
        fixed_vals <- temp_mat[, 3]
      }
      free_pos <- setdiff(1:k, fixed_pos)
      free_mask[fixed_pos] <- 0
    }
  }
  # fit probit model
  if(free_dim >0){
    penalty <- lambda0*(1-gammas)+lambda1*gammas  # bugs spotted previously lambda0*(1-gammas)+lambda1*(1+gammas)
    lambda_scale <- 1/n*sum(free_mask*penalty)/free_dim
    # no loading constraints
    if(free_dim ==k){
      result<-glmnet(thetas, y_tilde, family = binomial(link = "probit") ,alpha = 1, lambda=lambda_scale, penalty.factor = penalty,
                     standardize = FALSE, intercept=TRUE, standardize.response=FALSE, parallel = TRUE)
      return(result)
    }else{
      thetas_new <- thetas[, -fixed_pos]
      penalty_filter <-  penalty[free_pos]
      if(length(fixed_pos) > 1){ # more than one constraint
        offset <- rowSums(t(t(thetas[, fixed_pos]) * fixed_vals))
      }else{ # just one constraint
        offset <- thetas[, fixed_pos] * fixed_vals
      }
      if(free_dim == 1){ # glmnet needs two x's
        if(is.vector(thetas_new)){
            aug_num <- length(thetas_new)
        }else{
            aug_num <- nrow(thetas_new)
        }
        thetas_new <- cbind(thetas_new, rep(0, aug_num))
        penalty_filter <- c(penalty_filter, 0)
      }

      result<-glmnet(thetas_new, y_tilde, family = binomial(link = "probit") ,offset= offset, alpha = 1, lambda=lambda_scale, penalty.factor = penalty_filter,
                     standardize = FALSE, intercept=TRUE, standardize.response=FALSE, parallel = TRUE)
      alpha_est <- result$beta[1:free_dim]
      alpha_result <- rep(NA, k)
      alpha_result[fixed_pos] <- fixed_vals
      alpha_result[free_pos] <- alpha_est
      result$beta <- alpha_result
      return(result)
    }
  }else{ # no free dimension, only estimating intercepts
    offset <- rowSums(t(t(thetas) * fixed_vals))
    result <- glm.fit(rep(1,nrow(thetas)), y_tilde, family = binomial(link = "probit"), offset =offset, intercept = FALSE)
    result$beta <- fixed_vals
    result$a0 <- result$coefficients
    return(result)
  }
  return(result)
}

# M-Step to optimize the loading matrix alphas
optimize_loadings <- function(data, params, loading_constraints){
  probit_regs <- 1:params[["m"]] %>% 
    future_map(function(.x) train_single_probit(.x, data[, .x], params[["thetas"]], params[["gammas"]][.x, ], params[["lambda0"]], params[["lambda1"]], loading_constraints), .options = furrr_options(seed = TRUE))
  
  intercepts <- 1:params[["m"]] %>% 
    map(function(.z) probit_regs[[.z]]$a0) %>%
    unlist() %>%
    unname()
  
  betas <- 1:params[["m"]] %>% 
    map(function(.z) probit_regs[[.z]]$beta %>% as.numeric()) %>%
    unlist() %>%
    matrix(., params[["k"]], params[["m"]]) %>%
    t()
  
  new_loadings <- hash("alphas"= betas, "intercepts"= intercepts)
  return(new_loadings)
}


# M-step to optimize the C-paramter
optimize_c_params <- function(alpha, gammas){
  coefs<-apply(gammas,2,sum)
  K<-ncol(gammas)
  N<-nrow(gammas)
  
  eval_jac_g<-function(x){
    Matrix<--1*diag(K)
    for (i in (1:K-1)){
      Matrix[i,i+1]<-1}
    Matrix[K,K]<-0
    Matrix[K,1]<--1
    return(Matrix)
    }

  eval_g_ineq<-function(x){
    text<-paste("-x[",1:(K-1),"]+x[",(2:K),"]",sep="",collapse=",")
    text<-paste("c(",text,",-x[1])",sep="")
    eval(parse(text=text))
    }
  
  eval_grad_f<-function(x){
    paste1<-paste("-coefs[",1:(K-1),"]/","x[",1:(K-1),"]","+(N-coefs[",1:(K-1),"])/","(1-x[",1:(K-1),"])",sep="",collapse=",")
    paste2<-paste("-(alpha-1+coefs[",K,"])/","x[",K,"]","+(N-coefs[",K,"])/","(1-x[",K,"])",sep="",collapse=",")
    text<-paste(c(paste1,paste2),collapse=",")
    text<-paste("c(",text,")",sep="")
    eval(parse(text=text))}
  
  eval_f<-function(x){
    paste1<-paste("-coefs[",1:K,"]*","log(x[",1:K,"])",sep="",collapse="+")
    paste2<-paste("(N-coefs[",1:K,"])*","log(1-x[",1:K,"])",sep="",collapse="-")
    paste3<-paste("(alpha-1)*log(x[",K,"])",sep="")
    text<-paste(c(paste1,paste2,paste3),collapse="-")
    eval(parse(text=text))}
  
  
  opts<-list("algorithm"="NLOPT_LD_MMA","check_derivatives"=F,"xtol_rel"=10^-10)
  x0<-sort(rbeta(K,1,1),decreasing=TRUE)
  res<-nloptr(x0=x0,eval_f=eval_f,eval_grad_f=eval_grad_f,eval_g_ineq=eval_g_ineq,
              eval_jac_g_ineq=eval_jac_g,
              opts=opts,lb=rep(0,K),ub=rep(1,K))
  return(res$solution)
}

# m-step
m_step <- function(data, params, loading_constraints){
  loadings <- optimize_loadings(data, params, loading_constraints)
  params[["alphas"]] <- loadings[["alphas"]]
  params[["intercepts"]] <- loadings[["intercepts"]]
  params[["c_params"]] <- optimize_c_params(params[["ibp_alpha"]], params[["gammas"]])
  return(params)
}

# compute the loglikelihood
loglikhood <-function(data, param0s, normalize=TRUE){
  m <-params[["mc_samples"]]
  linear_term <- t(params[["alphas"]] %*% t(params[["thetas"]]) + params[["intercepts"]])
  ytilde <- data[rep(seq_len(nrow(data)), each = m), ]
  result <- sum(log(pnorm((2*ytilde-1)*linear_term)))/m
  if(normalize ==TRUE){
    result <- result /(nrow(data)*ncol(data))
  }
  return(result)
}



plot_loading_mat<-function(x, iter_idx, rotate){
  
  title <- paste0("iteartion ", iter_idx, ":  ", rotate, " rotation, abs(loadings)" )
  par(mar=c(2.1, 2.1, 1.1, 1.1))
  # topo.colors
  plot(x, digits=2, col= topo.colors, reorder=FALSE, axis.col=list(side=3, cex.axis=1), axis.row=list(cex.axis=1), main=title, key=NULL)
  
}

plot_rotation_mat<-function(x, iter_idx){
  title <- paste0("iteartion ", iter_idx, "PX-Rotation Matrix" )
  par(mar=c(2.1, 2.1, 1.1, 1.1))
  # topo.colors
  plot(x, digits=2, col= topo.colors, axis.col=list(side=3, cex.axis=1), axis.row=list(cex.axis=1), main=title, key=NULL)
  
}


# Main Program: EM Algroithm to Learn Factor Loadings from Item Response Data
# Args:
#   data: n by m item response
#   k: trunctaed dimension for the IBP prior
#   ibp_alpha: IBP intensity parameter
#   mc_samples: number of monte-carlo sampls to approximate E step
#   ssl_lambda0, ssl_lambda1: controlling the variance for SSL prior
#   max_iterations: maximum iterations allowed
#   epsilon: stopping threshold 
#   PX: parameter Expansion
#   varimax: whether performing varimax
#   Loading_constraints: 3-column matrix, first two columns are positions, last column is value
#   start: parameter initialization, use hash object
#   plot: whether to plot loading matrix after each iteration, useful for experimentation
#   stop_rotation: switch from PX-EM to EM after "stop_rotation" iterations
#   random_state: set seed
#   cores: number of cores for parallelization

probit_em_mirt <- function(data, k, ibp_alpha, mc_samples, ssl_lambda0, 
                           ssl_lambda1, max_iterations, epsilon, PX = TRUE, 
                           varimax = FALSE, loading_constraints= NULL, 
                           start = NULL, plot=FALSE, stop_rotation=100, random_state = 1,
                           cores=8){
    #future::plan(multicore, workers = cores)
    future::plan(multisession, workers = cores)
    set.seed(random_state)
    n <- dim(data)[1]
    m <- dim(data)[2]
    params <- initialize_params(n, m, k, ibp_alpha, mc_samples, loading_constraints, start, 
                                ssl_lambda0, ssl_lambda1)
    for (i in 1:max_iterations){
        params <- e_step(data, params)
        params <- m_step(data, params, loading_constraints)
        if(plot== TRUE & (i<=10 || i %% 1 ==0)){
          plot_loading_mat(abs(params[["alphas"]]), i, "before") 
        }
        # parameter expansion
        if(PX==TRUE & i== stop_rotation+1){
          print("switch to EM")
        }
        if(PX == TRUE & i <= stop_rotation ){
          a_mat <- t(params[["thetas"]]) %*% params[["thetas"]] /(n*mc_samples)
          new_alphas <- params[["alphas"]] %*% t(chol(a_mat))
          if(is.null(loading_constraints)==FALSE){
            for(j in 1:nrow(loading_constraints)){
              new_alphas[loading_constraints[j, 1], loading_constraints[j, 2]] <- loading_constraints[j, 3]
            }
          }
          params[["alphas"]] <- new_alphas
          if(plot== TRUE  & (i<=10 || i %% 1 ==0)){
            plot_loading_mat(abs(params[["alphas"]]), i, "PX-EM") 
            plot_rotation_mat(a_mat, i)
          }
        }
        
        # varimax step
        if(varimax==TRUE & i %% 5 ==0){
          new_alphas <-varimax(params[["alphas"]]+0.0000001)$loadings
          if(is.null(loading_constraints)==FALSE){
            for(j in 1:nrow(loading_constraints)){
              new_alphas[loading_constraints[j, 1], loading_constraints[j, 2]] <- loading_constraints[j, 3]
            }
          }
          params[["alphas"]] <- new_alphas
          
          plot_loading_mat(abs(params[["alphas"]]), i, "varimax")
        }
        max_diff <- max(abs(params[["alphas"]] - params[["alphas_prev"]]))
        # loglik <- loglikhood(data, params)
        # print (c(i, max_diff, loglik))
        print (c(i, max_diff))
        if(max_diff < epsilon){
          break
        }
        params[["alphas_prev"]] = params[["alphas"]]
        params[["intercepts_prev"]] = params[["intercepts"]]

    }
    # Extra M-step if rotated
    if(PX==TRUE){
      params <- m_step(data, params, loading_constraints)
    }
    return(params)
}




## Dynamic Posterior Exploration
# Args:
#     same argument meaning as probit_em_mirt
#     ssl_lambda0_path: a vector of ssl_lambda0 values for dynamic posterior exploration
#     pos_init: whether to zero out negative loading estimations after each ssl_lambda0 value. For MIRT model, this should be True.
dynamic_posterior_exploration <- function(data, k, ibp_alpha, mc_samples, ssl_lambda0_path, 
                           ssl_lambda1, pos_init =TRUE, max_iterations, epsilon, PX = FALSE, 
                           varimax = FALSE, loading_constraints= NULL, 
                           start = NULL, plot=FALSE, stop_rotation=100, random_state = 1,
                           cores=8){
  models <- vector("list",length = length(ssl_lambda0_path)) %>% setNames(paste("lambda0",ssl_lambda0_path, sep="_"))
  models[[1]] <- probit_em_mirt(data, k, ibp_alpha, mc_samples, ssl_lambda0_path[[1]], 
                                ssl_lambda1, max_iterations, epsilon, PX, 
                                varimax, loading_constraints,  start, plot=FALSE, stop_rotation, random_state,
                                cores)
  if(plot == TRUE){
    par(mfrow=c(1,2), mar=c(2, 2, 2, 2))
    plot(models[[1]]$alphas, digits=1, main= paste(ssl_lambda0_path[[1]], ssl_lambda1, sep="_"),  text.cell=list(cex=0.6), key=NULL) 
  }
  
  for(i in 2:length(ssl_lambda0_path)){
    if(pos_init == TRUE){
      new_start = hash("alphas"= pmax(models[[i-1]][["alphas"]], 0), "intercepts" = models[[i-1]][["intercepts"]], "c_params"= models[[i-1]][["c_params"]])
    }else{
      new_start = hash("alphas"= models[[i-1]][["alphas"]], "intercepts" = models[[i-1]][["intercepts"]], "c_params"= models[[i-1]][["c_params"]])
    }
    models[[i]] <- probit_em_mirt(data, k, ibp_alpha, mc_samples, ssl_lambda0_path[[i]], 
                                  ssl_lambda1, max_iterations, epsilon, PX, 
                                  varimax, loading_constraints,  new_start, plot=FALSE, stop_rotation, random_state,
                                  cores)
    if(plot==TRUE){
      par(mfrow=c(1,2), mar=c(2, 2, 2, 2))
      plot(models[[i]]$alphas, digits=1, main= paste(ssl_lambda0_path[[i]], ssl_lambda1, sep="_"),  text.cell=list(cex=0.6), key=NULL)
      
    }
    
  }
  return(models)
}
  
  
  
  
  


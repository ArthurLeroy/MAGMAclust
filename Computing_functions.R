library(tidyverse)
library(optimr)
library(Matrix)
library(MASS)
library(maotai)

##### USEFUL FUNCTIONS ########
'%notin%' <- Negate('%in%')

dmvnorm <- function (x, mu, inv_Sigma, log = FALSE, tol = 1e-06) 
{
  if (is.vector(x)) 
    x = t(as.matrix(x))
  n = length(mu)
  if (is.vector(mu)) {
    p <- length(mu)
    if (is.matrix(x)) {
      mu <- matrix(rep(mu, nrow(x)), ncol = p, byrow = TRUE)
    } 
  }
  else {
    p <- ncol(mu)
  }
  if (!all(dim(inv_Sigma) == c(p, p)) || nrow(x) != nrow(mu)) 
    stop("incompatible arguments")
  
  z <- t(x - mu)
  logdetS <- try(- determinant(inv_Sigma, logarithm = TRUE)$modulus,
                 silent=TRUE)
  attributes(logdetS) <- NULL
  
  ssq <- t(z) %*% inv_Sigma %*% z
  loglik <- -(n * (log(2*pi)) +  logdetS + ssq)/2
  if (log) return(loglik) else return(exp(loglik))
}

mat_dist = function(x,y)
{ ## return the matrice of distances between all pairs of vectors in x and y
  outer(x, y, Vectorize(function(p,q) sum((p - q)^2)) ) %>% return()
}

##### KERNELS DEFINITION ######
kernel = function(mat, theta = c(1, 0.5))
{ ## mat : the matrix M[i,j] =  t(i - j) %*% (i - j), for all i,j in the vector of timestamps
  ## theta : list of hyperparamaters of the kernel
  ####
  ## return : the value of dot production <f(t1),f(t2)> computed from the kernel
  
  exp(theta[[1]] - exp(-theta[[2]])/2 * mat) %>% return()
}

kernel_mu = function(mat, theta = c(1,0.5))
{ ## mat : the matrix M[i,j] =  t(i - j) %*% (i - j), for all i,j in the vector of timestamps
  ## theta : list of hyperparamaters of the kernel
  ####
  ## return : the value of dot production <f(t1),f(t2)> computed from the kernel
  
  exp(theta[[1]] - exp(-theta[[2]])/2 * mat) %>% return()
}

kern_to_cov = function(x, kern = kernel, theta = c(1, 0.5), sigma = 0.2)
{ ## x : vector or tibble of all timestamps for each individuals (input vector). 2 columns, 'ID' and 'Timestamp' required
  ## kern : indicates which kernel function to use to compute the covariance matrix
  ## theta : list of the required hyperparameters
  ####
  ## return : list of inverse covariance matrices (1 by individual) of the input vectors according to kernel() function
  
  if(is.vector(x))
  { 
    if(length(x) == 1)
    {
      mat = sigma^2 + exp(theta[[1]])
    } 
    else 
    {
      x = x %>% sort()
      M = dist(x)^2
      mat = as.matrix(kern(M, theta)) + diag(sigma^2 + exp(theta[[1]]) ,length(x))
    }
    
    
    if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
    rownames(mat) = paste0('X', x)
    colnames(mat) = paste0('X', x)
    return(mat)
  }
  
  ## If x is a tibble composed of observations from different individuals, identified by a variable 'ID'
  ## In this case, the list of HP must provide a set of HP for each individuals
  else
  { 
    floop = function(i)
    {
      if(theta %>% is.vector())
      {
        theta = list(c(theta, sigma))
        names(theta) = i
      }
      
      indiv = x %>% filter(ID == i) %>% arrange(Timestamp) %>% pull(Timestamp) %>% sort()
      if(length(indiv) == 1)
      {
        mat = theta[[i]][3]^2 + exp(theta[[i]][1])
      } 
      else 
      {
        M = dist(indiv)^2
        mat = as.matrix(kern(M, theta[[i]][1:2])) + diag(theta[[i]][3]^2 + exp(theta[[i]][1]), length(indiv))
      }
      
      if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
      rownames(mat) = paste0('X', indiv)
      colnames(mat) = paste0('X', indiv)
      return(mat)
    }
    list_mat = sapply(unique(x$ID), floop, simplify = FALSE, USE.NAMES = TRUE)
    if(length(list_mat) == 1){return(list_mat[[1]])}
    
    return(list_mat)
  } 
}

kern_to_inv = function(x, kern = kernel, theta = c(1, 0.5), sigma = 0.2)
{ ## x : vector or tibble of all timestamps for each individuals (input vector). 2 columns, 'ID' and 'Timestamp' required
  ## kern : indicates which kernel function to use to compute the covariance matrix
  ## theta : list of the required hyperparameters
  ####
  ## return : list of inverse covariance matrices (1 by individual) of the input vectors according to kernel() function
  #browser()
  if(is.vector(x))
  { 
    if(length(x) == 1)
    {
      mat = sigma^2 + exp(theta[[1]])
    } 
    else 
    {
      x = x %>% sort()
      M = dist(x)^2
      mat = as.matrix(kern(M, theta)) + diag(sigma^2 + exp(theta[[1]]) ,length(x))
    }
    inv = tryCatch(solve(mat), error = function(e){MASS::ginv(mat)})
    #inv = solve(mat)
    
    if(inv %>% dim() %>% is.null()){inv = as.matrix(inv)}
    rownames(inv) = paste0('X', x)
    colnames(inv) = paste0('X', x)
    return(inv)
  }
  
  ## If x is a tibble composed of observations from different individuals, identified by a variable 'ID'
  ## In this case, the list of HP must provide a set of HP for each individuals
  else
  { 
    floop = function(i)
    { 
      if(theta %>% class() == 'numeric')
      {
        theta = list(c(theta, sigma))
        names(theta) = i
      }
      
      indiv = x %>% filter(ID == i) %>% arrange(Timestamp) %>% pull(Timestamp) %>% sort()
      if(length(indiv) == 1)
      {
        mat = theta[[i]][3]^2 + exp(theta[[i]][1])
      } 
      else 
      {
        M = dist(indiv)^2
        mat = as.matrix(kern(M, theta[[i]][1:2])) + diag(theta[[i]][3]^2 + exp(theta[[i]][1]), length(indiv))
      }
      inv = tryCatch(solve(mat), error = function(e){MASS::ginv(mat)})
      #inv = solve(mat)
      
      if(inv %>% dim() %>% is.null()){inv = as.matrix(inv)}
      rownames(inv) = paste0('X', indiv)
      colnames(inv) = paste0('X', indiv)
      return(inv)
    }
    list_mat = sapply(unique(x$ID), floop, simplify = FALSE, USE.NAMES = TRUE)
    
    if(length(list_mat) == 1){return(list_mat[[1]])}
    
    return(list_mat)
  } 
}

##### LOGLIKELIHOOD FUNCTIONS ####
logL_GP_mod = function(hp, db, mean, kern, new_cov, pen_diag = NULL)
{ ## hp : vector or list of parameters of the kernel with format (a, b, sigma)
  ## db : tibble containing values we want to compute logL on. Required columns : Timestamp, Output
  ## mean : mean of the GP at corresponding timestamps
  ## kern : kernel used to compute the covariance matrix at corresponding timestamps
  ## new_cov : posterior covariance matrix of the mean GP (mu_0). Used to compute correction term (cor_term)
  ####
  ##return : value of the modified Gaussian log-likelihood for one GP as it appears in the model
  
  ## To avoid pathological behaviour of the opm optimization function in rare cases
  #if((hp %>% abs %>% sum) > 20){return(10^10)}
  
  if(length(mean) == 1){mean = rep(mean, nrow(db))} ## mean is equal for all timestamps
  sigma = ifelse((length(hp) == 3), hp[[3]], pen_diag) ## mean GP (mu_0) is noiseless and thus has only 2 hp
  
  inv =  kern_to_inv(db$Timestamp, kern, theta = hp[1:2], sigma) 
  #inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)}) ## fast or slow matrix inversion if singular
  
  LL_norm = - dmvnorm(db$Output, mean, inv, log = T) ## classic gaussian loglikelihood
  cor_term =  0.5 * (inv * new_cov) %>% sum() ## correction term (0.5 * Trace(inv %*% new_cov))
  
  return(LL_norm + cor_term)
}

logL_clust_multi_GP = function(hp, db, mu_k_param, kern)
{ ## hp : vector or list of parameters of the kernel with format (a, b, sigma)
  ## db : tibble containing values we want to compute logL on. Required columns : Timestamp, Output
  ## mu_k_param : list of parameters for the K mean Gaussian processes
  ## kern : kernel used to compute the covariance matrix at corresponding timestamps
  ####
  ##return : value of the modified Gaussian log-likelihood for one GP as it appears in the model
  #browser()
  sigma = ifelse((length(hp) == 3), hp[[3]], 0.1) ## mean GP (mu_0) is noiseless and thus has only 2 hp
  names_k = mu_k_param$mean %>% names()
  t_i = db$Timestamp
  y_i = db$Output
  i = unique(db$ID)
  
  inv =  kern_to_inv(db$Timestamp, kern, theta = hp[1:2], sigma) 
  
  LL_norm = - dmvnorm(y_i, rep(0, length(y_i)), inv, log = T) ## classic gaussian centered loglikelihood
  
  corr1 = 0
  corr2 = 0
  
  for(k in seq_len(length(names_k)))
  {
    tau_i_k = mu_k_param$tau_i_k[[k]][[i]]
    mean_mu_k = mu_k_param$mean[[k]] %>% filter(Timestamp %in% t_i) %>% pull(Output)
    corr1 = corr1 + tau_i_k * mean_mu_k
    corr2 = corr2 + tau_i_k * ( mean_mu_k %*% t(mean_mu_k) + mu_k_param$cov[[k]][paste0('X', t_i), paste0('X', t_i)] )
  }
  
  return( LL_norm - y_i %*% inv %*% corr1 + 0.5 * sum(inv * corr2) )
}

logL_GP_mod_common_hp_k = function(hp, db, mean, kern, new_cov, pen_diag = NULL)
{ ## hp : vector of common hyperparameters for all individuals. Format : c(a, b, sigma)
  ## db : list of the k tibble of data. Format :list( 'ID'  = Timestamp, Output)
  ## mean : list of the k means of the GP at union of observed timestamps
  ## kern : kernel used to compute the covariance matrices at corresponding timestamps
  ## new_cov : list of the k posterior covariance of the mean GP (mu_k). Used to compute correction term (cor_term)
  ####
  ##return : value of the modified Gaussian log-likelihood for the sum of the k mean GPs with same HPs
  
  ## To avoid pathological behaviour of the opm optimization function in rare cases
  #if((hp %>% abs %>% sum) > 20){return(10^10)}
  
  sigma = ifelse((length(hp) == 3), hp[[3]], pen_diag) ## mean GP is noiseless (0.1 is for computational issues) has only 2 hp
  list_ID_k = names(db)
  t_k = db[[1]] %>% pull(Timestamp)
  inv =  kern_to_inv(t_k, kern, theta = hp[1:2], sigma)
  
  LL_norm = 0
  cor_term = 0
  
  for(k in list_ID_k)
  {
    y_k = db[[k]] %>% pull(Output)
    
    LL_norm = LL_norm - dmvnorm(y_k, rep(mean[[k]], length(t_k)), inv, log = T) 
    cor_term = cor_term + 0.5 * (inv * new_cov[[k]]) %>% sum()  ##(0.5 * Trace(inv %*% new_cov))
  }
  return(LL_norm + cor_term)
}

logL_clust_multi_GP_common_hp_i = function(hp, db, mu_k_param, kern)
{ ## hp : vector or list of parameters of the kernel with format (a, b, sigma)
  ## db : tibble containing values we want to compute logL on. Required columns : Timestamp, Output
  ## mu_k_param : list of parameters for the K mean Gaussian processes
  ## kern : kernel used to compute the covariance matrix at corresponding timestamps
  ####
  ##return : value of the modified Gaussian log-likelihood for one GP as it appears in the model
  
  ## To avoid pathological behaviour of the opm optimization function in rare cases
  #if((hp %>% abs %>% sum) > 20){return(10^10)}
  
  sigma = ifelse((length(hp) == 3), hp[[3]], 0.2) ## mean GP (mu_0) is noiseless and thus has only 2 hp
  names_k = mu_k_param$mean %>% names()
  t = unique(db$Timestamp)
  
  sum_i = 0
  t_i_old = NULL
  
  for(i in unique(db$ID))
  {
    t_i = db %>% filter(ID == i) %>% pull(Timestamp)
    input_i = paste0('X', t_i)
    y_i = db %>% filter(ID == i) %>% pull(Output)
    
    corr1 = 0
    corr2 = 0
    
    for(k in seq_len(length(names_k)))
    {
      tau_i_k = mu_k_param$tau_i_k[[k]][[i]]
      mean_mu_k = mu_k_param$mean[[k]] %>% filter(Timestamp %in% t_i) %>% pull(Output)
      corr1 = corr1 + tau_i_k * mean_mu_k
      corr2 = corr2 + tau_i_k * ( mean_mu_k %*% t(mean_mu_k) + mu_k_param$cov[[k]][input_i, input_i] )
    }
    
    if( !identical(t_i, t_i_old) )
    {
      inv = kern_to_inv(t_i, kern, theta = hp[1:2], sigma) 
    }
    LL_norm = - dmvnorm(y_i, rep(0, length(y_i)), inv, log = T) ## classic gaussian centered loglikelihood
    
    sum_i = sum_i + LL_norm - y_i %*% inv %*% corr1 + 0.5 * sum(inv * corr2) 
    
  }
  return(sum_i)
}

logL_monitoring_VEM = function(hp, db, kern_i, kern_0, mu_k_param, m_k)
{ ## hp : list of parameters of the kernel for each individuals. Format : list(list(theta_k)_k, list(theta_i)_i)
  ## db : tibble containing values we want to compute logL on. Required columns : Timestamp, Output
  ## kern_i : kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
  ## kern_0 : kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0)
  ## m_k_param : parameters of the variational distributions of mean GPs (mu_k)
  ## m_k : prior value of the mean parameter of the mean GPs (mu_k). Length = 1 or nrow(db)
  ####
  ##return : value of expectation of joint log-likelihood of the model. The function to be maximised in step M
  #browser()
  
  floop = function(k)
  {
    logL_GP_mod(hp$theta_k[[k]], db = mu_k_param$mean[[k]], mean = m_k[[k]] , kern_0, mu_k_param$cov[[k]]) %>%
      return()
  }
  sum_ll_k = sapply(names(m_k), floop) %>% sum()
  
  floop2 = function(i)
  { 
    t_i = db %>% filter(ID == i) %>% pull(Timestamp)
    logL_clust_multi_GP(hp$theta_i[[i]], db %>% filter(ID == i), mu_k_param, kern_0) %>% return()
  }
  sum_ll_i = sapply(unique(db$ID), floop2) %>% sum()
  return(-sum_ll_k - sum_ll_i)
}

BIC = function(hp, db, kern_0, kern_i, m_k, mu_k_param, clust = T)
{#browser()
  list_clust = mu_k_param$mean %>% names
  nb_indiv = db %>% pull(ID) %>% n_distinct
  
  nb_hp_i = hp$theta_i %>% unlist %>% n_distinct
  nb_hp_k = hp$theta_k %>% unlist %>% n_distinct
  nb_clust = ifelse(clust, hp$pi_k %>% n_distinct, 1)
  size_mu = db %>% pull(Timestamp) %>% n_distinct
  ## Penalty terms for the BIC
  pen =  0.5 * (nb_hp_i + nb_hp_k + nb_clust - 1)*log(nb_indiv) 
  
  floop = function(i)
  {
    t_i = db %>% filter(ID == i) %>% pull(Timestamp)
    input_t_i = paste0('X', t_i)
    y_i = db %>% filter(ID == i) %>% pull(Output)
    inv_i = kern_to_inv(t_i, kern_i, theta = hp$theta_i[[i]][1:2], sigma = hp$theta_i[[i]][[3]])
    
    if(clust)
    {
      sum_LL = 0
      LL_Z = 0
      for(k in list_clust)
      {
        mu_k_ti = mu_k_param$mean[[k]] %>% filter(Timestamp %in% t_i) %>% pull(Output)
        pi_k = hp$pi_k[[k]]
        tau_i_k = mu_k_param$tau_i_k[[k]][[i]]
        cov_k_hat = mu_k_param$cov[[k]][input_t_i, input_t_i]
        
        sum_LL = sum_LL + tau_i_k * (dmvnorm(y_i, mu_k_ti, inv_i, log = T) - 0.5 * sum(inv_i * cov_k_hat)) 
        LL_Z = LL_Z + tau_i_k * ifelse((tau_i_k == 0|(pi_k == 0)), 0, log(pi_k/ tau_i_k)) # To avoid log(0)
      }
    }  
    else
    {
      mu_0 = mu_k_param$mean %>% filter(Timestamp %in% t_i) %>% pull(Output)
      cov_k_hat = mu_k_param$cov[input_t_i, input_t_i]
      
      sum_LL = dmvnorm(y_i, mu_0, inv_i, log = T) - 0.5 * sum(inv_i * cov_k_hat)
      LL_Z = 0
    }
    #paste0('LL = ', sum_LL, '|| LL_Z = ', LL_Z) %>% print()
    return( sum_LL + LL_Z )
  }
  LL_i = sapply(unique(db$ID), floop) %>% sum
  
  t = db %>% pull(Timestamp) %>% unique
  input_t = paste0('X', t)
  if(clust)
  {
    LL_k = 0
    for(k in list_clust)
    {
      mu_k_t  = mu_k_param$mean[[k]] %>% filter(Timestamp %in% t) %>% pull(Output)
      mean_k = rep(m_k[[1]], length(t))
      inv_k = kern_to_inv(t  , kern_0, theta = hp$theta_k[[k]][1:2], sigma = hp$theta_k[[k]][[3]])
      cov_k_hat = mu_k_param$cov[[k]][input_t, input_t] #+ diag(0.01, size_mu) # jitter term 
      diag = diag(cov_k_hat) %>% mean %>% `*`(0) %>% diag(size_mu) # jitter term for numerical problems
      
      LL_k = LL_k + dmvnorm(mu_k_t, mean_k, inv_k, log = T) - 0.5 * sum(inv_k * cov_k_hat) +
        0.5 * (size_mu +  size_mu * log(2*pi) + pdeterminant(cov_k_hat + diag)$modulus )
    }
  }
  else
  {
    mu_0 = mu_k_param$mean %>% filter(Timestamp %in% t) %>% pull(Output)
    mean_0 = rep(m_k[[1]], length(t))
    pen_diag = sapply(hp$theta_i, function(x) x[[3]]) %>% mean
    inv_0 = kern_to_inv(t , kern_0, theta = hp$theta_0[1:2], sigma = pen_diag)
    cov_hat = mu_k_param$cov[input_t, input_t] #+ diag(0.01, size_mu)
    diag = diag(cov_hat) %>% mean %>% `*`(0) %>% diag(size_mu) # jitter term for numerical problems
    
    LL_k = dmvnorm(mu_0, mean_0, inv_0, log = T) - 0.5 * sum(inv_0 * cov_hat) +
      0.5 * (size_mu +  size_mu * log(2*pi) + pdeterminant(cov_hat + diag)$modulus )
  }
  print(c(LL_i, LL_k, - pen))
  return(LL_i + LL_k - pen)
}


##### GRADIENT OF LogL FUNCTION #####
deriv_hp1 = function(mat, theta)
{
  exp(theta[[1]] - exp(-theta[[2]])/2 * mat) %>% return()
}

deriv_hp2 = function(mat, theta)
{
  grad = 0.5 * exp(- theta[[2]]) * mat
  
  return( exp(theta[[1]]) * grad * exp(- grad) )
}

gr_GP_mod = function(hp, db, mean, kern, new_cov, pen_diag = NULL)
{ 
  y = db$Output
  t = db$Timestamp
  
  sigma = ifelse((length(hp) == 3), hp[[3]], pen_diag)
  inv = kern_to_inv(t, kern, theta = hp[1:2], sigma)
  prod_inv = inv %*% (y - mean) 
  cste_term = prod_inv %*% t(prod_inv) + inv %*% ( new_cov %*% inv - diag(1, length(t)) )
  
  
  g_1 = 1/2 * (cste_term %*% kern_to_cov(t, deriv_hp1, theta = hp[1:2], sigma = 0)) %>%  diag() %>% sum()
  if(length(t) == 1)
  { ## Second hp has a 0 diagonal, and dist() return an error for only one observation 
    g_2 = 0
  }
  else
  {
    g_2 = 1/2 * (cste_term %*% as.matrix(deriv_hp2(dist(t)^2, theta = hp[1:2]) ))  %>%  diag() %>% sum()
  }
  
  if(length(hp) == 3)
  {
    g_3 = hp[[3]] * (cste_term %>% diag() %>% sum() )
    (- c(g_1, g_2, g_3)) %>% return()
  }
  else   (- c(g_1, g_2)) %>% return()
} 

gr_clust_multi_GP = function(hp, db, mu_k_param, kern)
{ 
  names_k = mu_k_param$mean %>% names()
  t_i = db$Timestamp
  y_i = db$Output
  i = unique(db$ID)
  
  corr1 = 0
  corr2 = 0
  for(k in seq_len(length(names_k)))
  {
    tau_i_k = mu_k_param$tau_i_k[[k]][[i]]
    mean_mu_k = mu_k_param$mean[[k]] %>% filter(Timestamp %in% t_i) %>% pull(Output)
    corr1 = corr1 + tau_i_k * mean_mu_k
    corr2 = corr2 + tau_i_k * ( mean_mu_k %*% t(mean_mu_k) + mu_k_param$cov[[k]][paste0('X', t_i), paste0('X', t_i)] )
  }
  
  sigma = ifelse((length(hp) == 3), hp[[3]], 0.1)
  inv = kern_to_inv(t_i, kern, theta = hp[1:2], sigma)
  prod_inv = inv %*% y_i 
  cste_term = (prod_inv - 2 * inv %*% corr1) %*% t(prod_inv)  + inv %*% ( corr2 %*% inv - diag(1, length(t_i)) )
  
  g_1 = 1/2 * (cste_term %*% kern_to_cov(t_i, deriv_hp1, theta = hp[1:2], sigma = 0)) %>%  diag() %>% sum()
  g_2 = 1/2 * (cste_term %*% as.matrix(deriv_hp2(dist(t_i)^2, theta = hp[1:2]) ))  %>%  diag() %>% sum()
  
  if(length(hp) == 3)
  {
    g_3 = hp[[3]] * (cste_term %>% diag() %>% sum() )
    (- c(g_1, g_2, g_3)) %>% return()
  }
  else   (- c(g_1, g_2)) %>% return()
}  

gr_GP_mod_common_hp_k = function(hp, db, mean, kern, new_cov, pen_diag = NULL)
{ 
  list_ID_k = names(db)
  t_k = db[[1]] %>% pull(Timestamp)
  sigma = ifelse((length(hp) == 3), hp[[3]], pen_diag)
  inv =  kern_to_inv(t_k, kern, theta = hp[1:2], sigma )
  
  g_1 = 0
  g_2 = 0
  g_3 = 0
  
  for(k in list_ID_k)
  {
    y_k = db[[k]] %>% pull(Output)
    
    prod_inv = inv %*% (y_k - mean[[k]]) 
    cste_term = prod_inv %*% t(prod_inv) + inv %*% ( new_cov[[k]] %*% inv - diag(1, length(t_k)) )
    
    g_1 = g_1 + 1/2 * (cste_term %*% kern_to_cov(t_k, deriv_hp1, theta = hp[1:2], sigma = 0)) %>% diag() %>% sum()
    g_2 = g_2 + 1/2 * (cste_term %*% as.matrix(deriv_hp2(dist(t_k)^2, theta = hp[1:2]) )) %>% diag() %>% sum()
    
    if(length(hp) == 3)
    {
      g_3 = g_3 + hp[[3]] * (cste_term %>% diag() %>% sum() )
    }
  }
  if(length(hp) == 3) return(- c(g_1, g_2, g_3)) else  return(- c(g_1, g_2))
}

gr_clust_multi_GP_common_hp_i = function(hp, db, mu_k_param, kern)
{ 
  sigma = ifelse((length(hp) == 3), hp[[3]], 0.2) ## mean GP (mu_0) is noiseless and thus has only 2 hp
  names_k = mu_k_param$mean %>% names()
  
  g_1 = 0
  g_2 = 0
  g_3 = 0
  t_i_old = NULL
  
  for(i in unique(db$ID))
  {
    
    t_i = db %>% filter(ID == i) %>% pull(Timestamp)
    input_i = paste0('X', t_i)
    y_i = db %>% filter(ID == i) %>% pull(Output)
    
    corr1 = 0
    corr2 = 0
    
    for(k in seq_len(length(names_k)))
    {
      tau_i_k = mu_k_param$tau_i_k[[k]][[i]]
      mean_mu_k = mu_k_param$mean[[k]] %>% filter(Timestamp %in% t_i) %>% pull(Output)
      corr1 = corr1 + tau_i_k * mean_mu_k
      corr2 = corr2 + tau_i_k * ( mean_mu_k %*% t(mean_mu_k) + mu_k_param$cov[[k]][input_i, input_i] )
    }
    
    if( !identical(t_i, t_i_old) )
    { ## We update the inverse cov matrix only if necessary (if different timestamps)
      inv = kern_to_inv(t_i, kern, theta = hp[1:2], sigma)
    }
    
    prod_inv = inv %*% y_i 
    cste_term = (prod_inv - 2 * inv %*% corr1) %*% t(prod_inv)  + 
      inv %*% ( corr2 %*% inv - diag(1, length(t_i)) )
    
    g_1 = g_1 + 1/2 * (cste_term %*% kern_to_cov(t_i, deriv_hp1, theta = hp[1:2], sigma = 0)) %>%  diag() %>% sum()
    g_2 = g_2 + 1/2 * (cste_term %*% as.matrix(deriv_hp2(dist(t_i)^2, theta = hp[1:2]) ))  %>%  diag() %>% sum()
    
    if(length(hp) == 3)
    {
      g_3 = g_3 + hp[[3]] * (cste_term %>% diag() %>% sum() )
    }
    t_i_old = t_i
  }  
  if(length(hp) == 3) return(- c(g_1, g_2, g_3)) else return(- c(g_1, g_2))
} 

##### EM FUNCTIONS #########
e_step_VEM = function(db, m_k, kern_0, kern_i, hp, old_tau_i_k)
{ ## db : full database with all individuals. Columns required : ID, Timestamp, Output
  ## m_k : prior means of the mu_k processes. List of vector of size K
  ## kern_i : kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
  ## kern_0 : kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0) 
  ## hp : set of hyper-parameters theta_k, theta_i, and pi_k (the proportion of individuals in cluster k)
  ## old_tau_i_k : values au tau_i_k from previous iterations. List(list(tau)_i)_k
  ## 
  ####
  ## return : mean and covariance parameters of the mean GP (mu_0)
  #browser()
  pi_k = hp$pi_k
  all_t = unique(db$Timestamp) %>% sort()
  t_clust = tibble('ID' = rep(names(m_k), each = length(all_t)) , 'Timestamp' = rep(all_t, length(m_k)))
  
  inv_k = kern_to_inv(t_clust, kern_0, hp$theta_k, sigma = 0)
  inv_i = kern_to_inv(db, kern_i, hp$theta_i, sigma = 0)
  value_i = base::split(db$Output, list(db$ID))
  
  ## Update each mu_k parameters
  floop = function(k)
  {
    new_inv = update_inv_VEM(prior_inv = inv_k[[k]], list_inv_i = inv_i, old_tau_i_k[[k]])
    tryCatch(solve(new_inv), error = function(e){MASS::ginv(new_inv)}) %>% 
      return()
  }
  cov_k = sapply(names(m_k), floop, simplify = FALSE, USE.NAMES = TRUE)
  
  floop2 = function(k)
  {
    weighted_mean = update_mean_VEM(m_k[[k]], inv_k[[k]], inv_i, value_i, old_tau_i_k[[k]])
    new_mean = cov_k[[k]] %*% weighted_mean %>% as.vector()
    tibble('Timestamp' = all_t, 'Output' = new_mean) %>% return()
  }
  mean_k = sapply(names(m_k), floop2, simplify = FALSE, USE.NAMES = TRUE)
  
  ## Update tau_i_k
  tau_i_k = update_tau_i_k_VEM(db, m_k, mean_k, cov_k, kern_i, hp, pi_k)
  
  list('mean' = mean_k, 'cov' = cov_k, 'tau_i_k' = tau_i_k) %>% return()
}

m_step_VEM = function(db, old_hp, list_mu_param, kern_0, kern_i, m_k, common_hp_k, common_hp_i)
{ ## db : db : full database with all individuals. Columns required : ID, Timestamp, Output
  ## old_hp : the set of hyper-parameters from the previous step of the EM
  ## list_mu_param : List of parameters of the K mean GPs. Format list('mean', 'cov', 'tau_i_k')
  ## kern_0 : kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0)
  ## kern_i : kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
  ## m_k : prior value of the mean parameter of the mean GPs (mu_k). Length = 1 or nrow(unique(db$Timestamp))
  ## common_hp_k : boolean indicating whether hp are common among mean GPs (for each mu_k)
  ## common_hp_i : boolean indicating whether hp are common among individual GPs (for each y_i)
  ####
  ## return : set of optimised hyper parameters for the different kernels of the model, and the pi_k
  list_ID_k = names(m_k)
  list_ID_i = unique(db$ID)
  t1 = Sys.time()
  if(common_hp_i)
  {
    param = opm(old_hp$theta_i[[1]], logL_clust_multi_GP_common_hp_i, gr = gr_clust_multi_GP_common_hp_i, db = db,
                mu_k_param = list_mu_param, kern = kern_i, method = "L-BFGS-B", control = list(kkt = F))[1,1:3]
    new_theta_i = param %>% list() %>% rep(length(list_ID_i))  %>% setNames(nm = list_ID_i) 
  }
  else
  {
    funloop2 = function(i)
    {
      t_i = db %>% filter(ID == i) %>% pull(Timestamp)
      opm(old_hp$theta_i[[i]] %>% unlist(), logL_clust_multi_GP, gr = gr_clust_multi_GP, db = db %>% filter(ID == i),
          mu_k_param = list_mu_param, kern = kern_i, method = "L-BFGS-B", control = list(kkt = F))[1,1:3] %>% return()
    }
    new_theta_i = sapply(list_ID_i, funloop2, simplify = FALSE, USE.NAMES = TRUE)
  }
  t2 = Sys.time()
  pen_diag = sapply(new_theta_i, function(x) x[[3]]) %>% mean
  
  if(common_hp_k)
  {
    param = c(opm(old_hp$theta_k[[1]], logL_GP_mod_common_hp_k, gr = gr_GP_mod_common_hp_k, db = list_mu_param$mean,
                  mean = m_k, kern = kern_0, new_cov = list_mu_param$cov, pen_diag = pen_diag, method = "L-BFGS-B",
                  control = list(kkt = F))[1,1:2] %>% unlist(), pen_diag)
    new_theta_k = param %>% list() %>% rep(length(list_ID_k))  %>% setNames(nm = list_ID_k) 
  }
  else
  {
    funloop = function(k)
    {
      c(opm(old_hp$theta_k[[k]][1:2], logL_GP_mod, gr = gr_GP_mod,  db = list_mu_param$mean[[k]],
            mean = m_k[[k]], kern = kern_0, new_cov = list_mu_param$cov[[k]], pen_diag = pen_diag,
            method = "L-BFGS-B", control = list(kkt = FALSE))[1,1:2] %>% unlist(), pen_diag) %>% return()
    }
    new_theta_k = sapply(list_ID_k, funloop, simplify = FALSE, USE.NAMES = TRUE)
  }
  
  pi_k = sapply( list_mu_param$tau_i_k, function(x) x %>% unlist() %>% mean() ) 
  
  t3 = Sys.time()
  print(c('mu_k:',t3 - t2, 'indiv:', t2 - t1))
  list('theta_k' = new_theta_k, 'theta_i' = new_theta_i, 'pi_k' = pi_k) %>% return()
}

##### UPDATE FUNCTIONS ####
update_inv_VEM = function(prior_inv, list_inv_i, tau_i_k)
{ ## prior_inv : inverse of the covariance matrix of the prior mean GP (mu_0). dim = all timestamps 
  ## list_inv_i : list of inverse of the covariance matrices of each individuals. dim = timestamps of i 
  ## tau_i_k : list of probability for each indiv 'i' to belong to cluster k
  ####
  ## return : inverse of the covariance of the posterior mean GP (mu_0 | (y_i)_i). dim = (all timestamps)^2 
  #browser()
  
  new_inv = prior_inv
  for(x in list_inv_i %>% names())
  {
    inv_i = list_inv_i[[x]]
    common_times = intersect(row.names(inv_i), row.names(new_inv))
    new_inv[common_times, common_times] = new_inv[common_times, common_times] + 
      tau_i_k[[x]] * inv_i[common_times, common_times]
  }
  return(new_inv)
}

update_mean_VEM = function(prior_mean, prior_inv, list_inv_i, list_value_i, tau_i_k)
{ ## prior_mean : mean parameter of the prior mean GP (mu_k)
  ## prior_inv : inverse of the covariance matrix of the prior mean GP (mu_k). dim = (all timestamps)^2 
  ## list_inv_i : list of inverse of the covariance matrices of each individuals. dim = (timestamps of i)^2  
  ## list_value_i : list of outputs (y_i) for each individuals. dim = (timestamps of i) x 1
  ## tau_i_k : list of probability for each indiv 'i' to belong to cluster k 
  ####
  ## return : mean parameter of the posterior mean GP (mu_0 | (y_i)_i). dim = (all timestamps) x 1
  
  if(length(prior_mean) == 1){prior_mean = rep(prior_mean, ncol(prior_inv))}
  weighted_mean = prior_inv %*% prior_mean
  #row.names(weithed_mean) = row.names(prior_inv)
  
  for(i in list_inv_i %>% names()) 
  {
    weighted_i = tau_i_k[[i]] * list_inv_i[[i]] %*% list_value_i[[i]]
    #row.names(weithed_i) = row.names(list_inv_i[[j]])
    
    common_times = intersect(row.names(weighted_i), row.names(weighted_mean))
    weighted_mean[common_times,] = weighted_mean[common_times,] + weighted_i[common_times,]
  }
  return(weighted_mean)
}

update_tau_i_k_VEM = function(db, m_k, mean_k, cov_k, kern_i, hp, pi_k)
{
  #browser()
  c_i = 0
  c_k = 0
  mat_logL = matrix(NA, nrow = length(names(m_k)), ncol = length(unique(db$ID)) )
  
  for(i in unique(db$ID))
  { c_i = c_i + 1
  t_i = db %>% filter(ID == i) %>% pull(Timestamp)
  for(k in names(m_k))
  { c_k = c_k + 1
  
  mat_logL[c_k,c_i] = - logL_GP_mod(hp$theta_i[[i]], db %>% filter(ID == i),
                                    mean_k[[k]] %>% filter(Timestamp %in% t_i) %>% pull(Output), kern_i,
                                    cov_k[[k]][paste0('X',t_i),paste0('X',t_i)])
  if(is.na(mat_logL[c_k,c_i])){print(i)}
  } 
  c_k = 0
  } 
  
  ## We need to use the 'log-sum-exp' trick: exp(x - max(x)) / sum exp(x - max(x)) to remain numerically stable
  mat_L = mat_logL %>% apply(2,function(x) exp(x - max(x)))
  
  (pi_k * mat_L) %>% apply(2,function(x) x / sum(x)) %>%
    `rownames<-`(names(m_k)) %>%
    `colnames<-`(unique(db$ID)) %>% 
    apply(1, as.list) %>%
    return()
}

update_tau_star_k_EM = function(db, mean_k, cov_k, kern, hp, pi_k)
{
  #browser()
  names_k = names(mean_k)
  c_k = 0
  mat_logL = rep(NA, length(names_k))
  t_i = db %>% pull(Timestamp)
  pi = unlist(pi_k)
  
  for(k in names_k)
  { 
    c_k = c_k + 1
    
    mean = mean_k[[k]] %>% filter(Timestamp %in% t_i) %>% pull(Output)
    cov =  (kern_to_cov(db$Timestamp, kern, hp[1:2], hp[3]) + cov_k[[k]][paste0('X',t_i), paste0('X',t_i)] )
    inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)})
    
    mat_logL[c_k] =  dmvnorm(db %>% pull(Output) , mean, inv, log = T) ## classic gaussian loglikelihood
  } 
  ## We need to use the 'log-sum-exp' trick: exp(x - max(x)) / sum exp(x - max(x)) to remain numerically stable
  mat_L = exp(mat_logL - max(mat_logL)) 
  
  ((pi * mat_L)/ sum(pi * mat_L)) %>% as.vector %>% split(names_k) %>% return()
}

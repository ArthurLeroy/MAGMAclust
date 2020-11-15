library(tidyverse)
library(MASS)
library(Matrix)
library(mvtnorm)
library(optimr)
library(fda)
library(plotly)
library(gganimate)
library(transformr)
library(gifski)
library(png)

source('Computing_functions.R')

################ TRAINING FUNCTIONS ####
training_VEM = function(db, prior_mean_k, ini_hp = list('theta_k' = c(1, 1, 0.2), 'theta_i' = c(1, 1, 0.2)),
                        kern_0 = kernel_mu, kern_i = kernel, ini_tau_i_k = NULL, common_hp_k = T, common_hp_i = T)
{ ## db : database with all individuals in training set. Column required : 'ID', Timestamp',  'Output'
  ## prior_mean : prior mean parameter of the K mean GPs (mu_k)
  ## ini_hp : initial values of HP for the kernels
  ## kern_0 : kernel associated to covariance functions of the mean GP
  ## kern_i : kernel associated to common covariance functions of all individuals GPs
  ## ini_tau_i_k : initial values of probabiliy to belong to each cluster for each individuals. 
  ####
  ## return : list of trained HP, boolean to indicate convergence
  #browser()
  n_loop_max = 15
  list_ID = unique(db$ID)
  ID_k = names(prior_mean_k)
  hp = list('theta_k' = ini_hp$theta_k %>% list() %>% rep(length(ID_k))  %>% setNames(nm = ID_k), 
            'theta_i' = ini_hp$theta_i %>% list() %>% rep(length(list_ID))  %>% setNames(nm = list_ID))
  cv = 'FALSE'
  if(is.null(ini_tau_i_k)){ini_tau_i_k = ini_tau_i_k(db, k = length(ID_k), nstart = 50)}
  tau_i_k = ini_tau_i_k
  hp[['pi_k']] = sapply( tau_i_k, function(x) x %>% unlist() %>% mean() ) 
  logLL_monitoring = - Inf
  t1 = Sys.time()
  
  for(i in 1:n_loop_max)
  { 
    print(i)
    ## E-Step
    param = e_step_VEM(db, prior_mean_k, kern_0, kern_i, hp, tau_i_k)  
    
    ## Monitoring of the LL
    new_logLL_monitoring = logL_monitoring_VEM(hp, db, kern_i, kern_0, mu_k_param = param , m_k = prior_mean_k) 
    #0.5 * (length(param$cov) * nrow(param$cov[[1]]) +
    #Reduce('+', lapply(param$cov, function(x) log(det(x)))) )
    print(new_logLL_monitoring)
    diff_moni = new_logLL_monitoring - logLL_monitoring
    
    if(diff_moni < - 0.1){warning('Likelihood descreased')}
    
    ## M-Step
    new_hp = m_step_VEM(db, hp, list_mu_param = param, kern_0, kern_i, prior_mean_k, common_hp_k, common_hp_i)
    
    if(new_hp %>% anyNA(recursive = T))
    {
      print(paste0('The M-step encountered an error at iteration : ', i))
      print('Training has stopped and the function returns values from the last valid iteration')
      break
    }
    
    ## Testing the stoping condition
    logL_new = logL_monitoring_VEM(new_hp, db, kern_i, kern_0, mu_k_param = param, m_k = prior_mean_k)
    eps = (logL_new - logL_monitoring_VEM(hp, db, kern_i, kern_0, mu_k_param = param, m_k = prior_mean_k)) / 
      abs(logL_new)
    
    print(c('eps', eps))
    if(eps < 1e-1)
    {
      if(eps > 0){hp = new_hp}
      cv = TRUE
      break
    }
    hp = new_hp
    tau_i_k = param$tau_i_k
    logLL_monitoring = new_logLL_monitoring
  }
  t2 = Sys.time()
  list('hp' = new_hp, 'convergence' = cv,  'param' = param,
       'Time_train' =  difftime(t2, t1, units = "secs")) %>% 
    return()
}

model_selection = function(db, k_grid = 1:5, ini_hp = list('theta_k' = c(1, 1, 0.2), 'theta_i' = c(1, 1, 0.2)),
                           kern_0 = kernel_mu, kern_i = kernel, ini_tau_i_k = NULL, common_hp_k = T, common_hp_i = T, 
                           plot = T)
{
  floop = function(K)
  {
    print(paste0('K = ', K))
    prior_mean_k = rep(0, K) %>% setNames(paste0('K', seq_len(K))) %>% as.list
    
    if(K == 1)
    {
      model = training(db, 0, list('theta_0'=ini_hp$theta_k[1:2], 'theta_i'=ini_hp$theta_i), kern_0, kern_i, common_hp_i)
    }
    else
    {
      model = training_VEM(db, prior_mean_k, ini_hp, kern_0, kern_i, ini_tau_i_k, common_hp_k, common_hp_i)
    }
    model[['BIC']] = BIC(model$hp, db, kern_0, kern_i, prior_mean_k, model$param, K != 1)
    return(model)
  }
  res = k_grid %>% sapply(floop, simplify = FALSE, USE.NAMES = TRUE) %>% setNames(paste0('K = ', k_grid))
  
  db_plot = tibble('K' = k_grid, 'BIC' = res %>% map_dbl('BIC'))
  
  res$K_max_BIC = db_plot %>% filter(BIC == max(BIC)) %>% pull(K)
  
  if(plot)
  {
    res$plot = ggplot(db_plot, aes(x = K, y = BIC)) + geom_point() +
      geom_line() + theme_classic()
  }
  
  return(res)
}

train_new_gp_EM = function(db, param_mu_k, ini_hp, kern_i, hp_i = NULL)
{
  mean_mu_k = param_mu_k$mean
  cov_mu_k = param_mu_k$cov
  pi_k = lapply(param_mu_k$tau_i_k, function(x) Reduce("+", x)/ length(x))
  if(is.null(hp_i))
  {
    n_loop_max = 25
    hp = ini_hp$theta_i
    
    for(i in 1:n_loop_max)
    { 
      ## E step
      tau_k = update_tau_star_k_EM(db, mean_mu_k, cov_mu_k, kern_i, hp, pi_k)
      
      ## M step
      LL_GP<- function(hp, db, kern_i) 
      {
        floop = function(k)
        {
          mean = mean_mu_k[[k]] %>% filter(Timestamp %in% t) %>% pull(Output)
          cov = (kern_to_cov(db$Timestamp, kern_i, theta = hp[1:2], sigma = hp[3]) +
                   cov_mu_k[[k]][paste0('X',t), paste0('X',t)])
          inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)})
          
          (db$Timestamp - tau_k[[k]] * dmvnorm(db$Output, mean, inv, log = T)) %>%
            return()
        }
        sapply(names(mean_mu_k), floop) %>% sum() %>% return()
      }
      new_hp = opm(hp, LL_GP, db = db, kern = kern_i, method = "L-BFGS-B", control = list(kkt = FALSE))[1,1:3] 
      if(new_hp %>% anyNA(recursive = T))
      {
        print(paste0('The M-step encountered an error at iteration : ', i))
        print('Training has stopped and the function returns values from the last valid iteration')
        break
      }
      
      ## Testing the stoping condition
      eps = (new_hp - hp) %>% abs %>% sum
      print(c('tau_k', tau_k %>% unlist))
      print(c('eps', eps))
      if(eps>0 & eps < 1e-3)
      {
        cv = 'TRUE'
        break
      }
      hp = new_hp
    }
  }
  else
  {
    new_hp = hp_i$hp$theta_i[[1]]
    tau_k = update_tau_star_k_EM(db, mean_mu_k, cov_mu_k, kern_i, new_hp, pi_k)
  }
  list('theta_new' = new_hp , 'tau_k' = tau_k) %>% return()
}

################ PREDICTION FUNCTIONS ################
posterior_mu_k = function(db, timestamps, m_k, kern_0, kern_i, list_hp)
{ ## db : tibble of data columns required ('ID', 'Timestamp',  'Output')
  ## timestamps : timestamps on which we want a prediction
  ## prior_mean : prior mean value of the mean GP (scalar value or vector of same length as 'timestamps')
  ## kern_0 : kernel associated to the covariance function of the mean GP
  ####
  ## return : pamameters of the mean GP at timestamps
  #browser()
  hp_i = list_hp$hp$theta_i
  hp_k = list_hp$hp$theta_k
  for(k in names(hp_k)){hp_k[[k]][3] = 0.1}
  tau_i_k = list_hp$param$tau_i_k
  t_clust = tibble('ID' = rep(names(hp_k), each = length(timestamps)) , 'Timestamp' = rep(timestamps, length(hp_k)))
  inv_k = kern_to_inv(t_clust, kern_0, hp_k, sigma = 0)
  inv_i = kern_to_inv(db, kern_i, hp_i, sigma = 0)
  value_i = base::split(db$Output, list(db$ID))
  
  ## Update each mu_k parameters
  floop = function(k)
  {
    new_inv = update_inv_VEM(prior_inv = inv_k[[k]], list_inv_i = inv_i, tau_i_k[[k]])
    tryCatch(solve(new_inv), error = function(e){MASS::ginv(new_inv)}) %>%
      return()
  }
  cov_k = sapply(names(hp_k), floop, simplify = FALSE, USE.NAMES = TRUE)
  
  floop2 = function(k)
  {
    weighted_mean = update_mean_VEM(m_k[[k]], inv_k[[k]], inv_i, value_i, tau_i_k[[k]])
    new_mean = cov_k[[k]] %*% weighted_mean %>% as.vector()
    tibble('Timestamp' = timestamps, 'Output' = new_mean) %>% return()
  }
  mean_k = sapply(names(hp_k), floop2, simplify = FALSE, USE.NAMES = TRUE)
  
  #names(mean_mu) = paste0('X', t_mu)
  list('mean' = mean_k, 'cov' = cov_k, 'tau_i_k' = tau_i_k) %>% return()
}

pred_gp_clust = function(db, timestamps = NULL, list_mu, kern, hp)
{ ## db : tibble of data columns required ('Timestamp',  'Output')
  ## timestamps : timestamps on which we want a prediction
  ## mean_mu : mean value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
  ## cov_mu : covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
  ## kern : kernel associated to the covariance function of the GP
  ## hp : list of hyperparameters and tau_k for the new individual
  ####
  ## return : pamameters of the gaussian density predicted at timestamps 
  tn = db %>% pull(Timestamp)
  input = paste0('X', db$Timestamp)
  yn = db %>% pull(Output)
  
  if(is.null(timestamps)){timestamps = seq(min(tn), max(tn), length.out = 500)}
  input_t = paste0('X', timestamps)
  
  theta = hp$theta_new[1:2]
  sigma = hp$theta_new[3]
  tau_k = hp$tau_k
  
  floop = function(k)
  {
    mean_mu_obs = list_mu$mean[[k]] %>% filter(Timestamp %in% tn) %>% pull(Output)
    mean_mu_pred = list_mu$mean[[k]] %>% filter(Timestamp %in% timestamps) %>% pull(Output)
    cov_mu = list_mu$cov[[k]]
    
    cov_tn_tn = (kern_to_cov(tn, kern, theta, sigma) + cov_mu[input, input])
    inv_mat = tryCatch(solve(cov_tn_tn), error = function(e){MASS::ginv(cov_tn_tn)}) 
    cov_tn_t = kern(mat_dist(tn, timestamps), theta) + cov_mu[input,input_t]
    cov_t_t = kern_to_cov(timestamps, kern, theta, sigma) + cov_mu[input_t ,input_t] 
    
    tibble('Timestamp' = timestamps, 
           'Mean' = (mean_mu_pred + t(cov_tn_t) %*% inv_mat %*% (yn - mean_mu_obs)) %>% as.vector(),
           'Var' =  (cov_t_t - t(cov_tn_t) %*% inv_mat %*% cov_tn_t) %>% diag %>% as.vector,
           'tau_k' = tau_k[[k]]) %>% return()
  }
  pred = sapply(names(list_mu$mean), floop, simplify = FALSE, USE.NAMES = TRUE)
  
  return(pred)
}

pred_gp_clust_animate = function(db, timestamps = NULL, list_mu, kern, hp)
{
  ## Inputs : same as for a classic GP prediction
  ####
  ## return : tibble of classic GP predictions but with an inscreasing number of data points considered as 'observed'
  db %>% arrange(Timestamp)
  all_pred = tibble()
  
  if(is.null(timestamps)){timestamps = seq(min(db$Timestamp), max(db$Timestamp), length.out = 500)}
  
  for(j in 1:nrow(db))
  {
    pred_j = pred_gp_clust(db[1:j,], timestamps, list_mu, kern, hp)
    
    for(k in list_mu$mean %>% names)
    {
      pred_j[[k]] = pred_j[[k]] %>% mutate(Nb_data = j, Cluster = k)
    }
    pred_j_all = pred_j %>% bind_rows
    
    all_pred = all_pred %>% rbind(pred_j_all)
  }
  return(all_pred)
}

full_algo_clust = function(db, new_db, timestamps, kern_i, ini_tau_i_k = NULL, common_hp_k = T, common_hp_i = T,
                           prior_mean = NULL, kern_0 = NULL, list_hp = NULL, mu_k = NULL, ini_hp = NULL, hp_new_i = NULL)
{ ## db : Database containing all training data from all individuals. Column: ID - Timestamp - Output.
  ## new_db : Database containing data for a new individual we want a prediction on.
  ## timestamps : Timestamps we want to predict at.
  ## kern_i : Kernel associated to individual GPs.
  ## ini_tau_i_k : Initial probability to belong to each cluster for individuals. Optional, if list_hp is given
  ## plot : Boolean indicating whether we want to display a graph at the end of computations.
  ## prior_mean : Prior arbitrary value for the mean processes. Optional, not needed if 'mu' is given.
  ## kern_0 : Kernel associated to the mean GPs. Optional, not needed if 'mu' is given.
  ## list_hp : Hyper-parameters for all individuals in training set. Optional, computed if NULL.
  ## mu_k : list containing parameters of the mean GPs at all prediction timestamps. Optional, computed if NULL.
  ## ini_hp : Initial values of the HP to start the training. Optional, not needed if 'list_hp' is given.
  ## hp_new_i : Hyper-pameters for the new individual to predict. Optional, computed if NULL. 
  ####
  ## return : predicted GP parameters | posterior K mean processes | all trained hyperparameters
  #browser()
  if(is.null(list_hp))
  {
    list_hp = training_VEM(db, prior_mean, ini_hp, kern_0, kern_i, ini_tau_i_k, common_hp_k, common_hp_i)
  }
  
  t_pred = timestamps %>% union(unique(db$Timestamp)) %>% union(unique(new_db$Timestamp)) %>% sort()
  if(is.null(mu_k)){mu_k = posterior_mu_k(db, t_pred, prior_mean, kern_0, kern_i, list_hp)}
  
  if(is.null(hp_new_i))
  {
    if(common_hp_i)
    {
      hp_new_i = train_new_gp_EM(new_db, mu_k, ini_hp, kern_i, hp_i = list_hp) 
    }
    else{hp_new_i = train_new_gp_EM(new_db, mu_k, ini_hp, kern_i)}
  }
  
  pred = pred_gp_clust(new_db, timestamps, mu_k, kern_i, hp_new_i)
  
  # if(plot)
  # {
  #   (pred %>% plot_gp(data = rbind(new_db, db)) + geom_point(mu, aes(Timestamp, Output, col = ID))) %>% print()
  # }
  
  list('Prediction' = pred, 'Mean_processes' = mu_k,
       'Hyperparameters' = list('other_i' = list_hp, 'new_i' = hp_new_i)) %>% 
    return()
}

pred_max_cluster = function(tau_i_k)
{ 
  tau_i_k %>% as_tibble %>% unnest %>% apply(1, which.max) %>% return()
}

################ PLOT FUNCTIONS ######################

plot_db = function(db, cluster = F, legend = F)
{ ## Visualize smoothed raw data. Format : ID, Timestamp, Output
  if(cluster)
  {
    ggplot(db) + geom_smooth(aes(Timestamp, Output, group = ID, color = Cluster)) +
      geom_point(aes(Timestamp, Output, group = ID, color = Cluster)) + guides(col = legend)
  }
  else
  {
    ggplot(db) + geom_smooth(aes(Timestamp, Output, color = ID)) +
      geom_point(aes(Timestamp, Output, color = ID)) + guides(col = legend)
  }
}

plot_gp_clust = function(pred, cluster = 'all', data = NULL, data_train = NULL, mean_k = NULL,
                         col_clust = F, legend = F)
{ ## pred : tibble coming out of the pred_gp() function, columns required : 'Timestamp', 'Mean', 'Var'
  ## cluster : string indicating the cluster to plot from or 'all' for the full GPs mixture.
  ## data : tibble of observational data, columns required : 'Timestamp', 'Output'
  ## data_train : tibble of training dataset, columns required : 'Timestamp', 'Output'
  ## mean : boolean indicating whether we shall plot the mean processes or not. 
  ## col_clust : boolean indicating whether we color according to clusters or individuals. 
  ####
  ## return : plot of the predicted curve of the GP with the 0.95 confidence interval (optional : data points)
  #browser()
  if(cluster == 'all')
  { mean_all = 0
  for(k in names(pred))
  {
    mean_all = mean_all + pred[[k]]$tau_k * pred[[k]]$Mean
  }
  pred_gp = tibble('Timestamp' = pred[[1]]$Timestamp, 'Mean' = mean_all, 'Var' = 0)
  }
  else{pred_gp = pred[[cluster]]}
  
  gg = ggplot() +
    geom_line(data = pred_gp, aes(x = Timestamp, y = Mean), color = 'blue') +
    geom_ribbon(data = pred_gp, aes(x = Timestamp, ymin = Mean - 1.96* sqrt(Var),
                                    ymax = Mean +  1.96* sqrt(Var)), alpha = 0.2) + ylab('Output')
  
  if(!is.null(data_train))
  {
    if(col_clust){gg = gg + geom_point(data = data_train, aes(x = Timestamp, y = Output, col = Cluster), shape = 4)}
    else{gg = gg + geom_point(data = data_train, aes(x = Timestamp, y = Output, col = ID), shape = 4)}
  }
  if(!is.null(data)){gg = gg + geom_point(data = data, aes(x = Timestamp, y = Output), size = 2, shape = 18)}
  if(!is.null(mean_k))
  {
    for(k in names(mean_k$mean))
    {
      gg = gg + geom_line(data = mean_k$mean[[k]],
                          aes(x = Timestamp, y = Output), linetype = 'dashed') 
    }
  }
  return(gg + guides(col = legend))
}

plot_heat_clust =  function(pred, ygrid, data = NULL, data_train = NULL,
                            mean_k = NULL,  col_clust = T, interactive = F, legend = F, 
                            animate = F)
{ ## pred: tibble coming out of the pred_gp_lust() function
  ## cluster : string indicating the cluster to plot from or 'all' for the full GPs mixture.
  ## data : tibble of observational data for the new individual, columns required : 'Timestamp', 'Output' (optional)
  ## data_train : tibble of training data for the model, , columns required : 'Timestamp', 'Output' (optional)
  ## mean : (optional)
  ## ygrid : grid on y-axis to plot the heatmap on (optional, default = range of data +- 2 sd)
  ## interactive : boolean indicating whether the output plot should be interactive (plotly)
  ## CI : boolean. If True, display the heatmap of Credible Intervals. If False, display the heatmap of likelihood
  ####
  ## return : An heatmap illustrating the GP predition results (Optionnal diplay of raw data and mean process)
  #browser()
  combine = function(pred, y, y_decal)
  { proba = 0
  for(k in pred$Cluster %>% unique)
  {
    #pred_gp = pred[[k]]
    pred_gp = pred %>% filter(Cluster == k)
    proba = proba + pred_gp$tau_k * (pnorm(y, mean = pred_gp$Mean, sd =  sqrt(pred_gp$Var)) -
                                       pnorm(y_decal, mean = pred_gp$Mean, sd =  sqrt(pred_gp$Var)))
  }
  tibble('Ygrid' = y, 'Timestamp' = pred %>% filter(Cluster == k) %>% pull(Timestamp),
         'Proba' = proba) %>% return()
  }
  
  if(animate)
  {
    db_heat = c()
    for(a in pred$Nb_data %>% unique)
    {
      y_decal = ygrid[1]
      db_h = c()
      for(i in ygrid[-1])
      {
        db_h = db_h %>% rbind(combine(pred %>% filter(Nb_data == a), i, y_decal))
        y_decal = i 
      }      
      db_h = db_h %>% mutate(Nb_data = a)
      db_heat = db_heat %>% rbind(db_h)
    }
  }
  
  else
  {
    y_decal = ygrid[1]
    db_heat = c()
    for(i in ygrid[-1])
    {
      db_heat = db_heat %>% rbind(combine(pred, i, y_decal))
      y_decal = i 
    }
  }
  
  gg = ggplot() + geom_tile(data = db_heat, aes(Timestamp, Ygrid, fill = Proba)) + 
    #scale_fill_gradientn(low = "white", high = "darkblue") +
    scale_fill_gradientn(colours=c('white', '#AB47BC', '#9C27B0','#8E24AA','#7B1FA2', '#6A1B9A', '#4A148C')) + 
    #scale_fill_distiller(palette = "RdPu", trans = "reverse") + 
    theme_classic() + ylab("Output") + labs(fill = "Proba")
  
  ## Display data/training data/mean process if provided
  if(!is.null(data_train))
  {
    if(col_clust){gg = gg + geom_point(data = data_train, aes(x = Timestamp, y = Output, col = Cluster), shape = 4)}
    else{gg = gg + geom_point(data = data_train, aes(x = Timestamp, y = Output, col = ID), size = 0.4, alpha = 0.5)}
  }
  if(!is.null(data)){gg = gg + geom_point(data = data, aes(x = Timestamp, y = Output), size = 2, shape = 18)}
  if(!is.null(mean))
  {
    for(k in names(mean_k$mean))
    {
      gg = gg + geom_line(data = mean_k$mean[[k]],
                          aes(x = Timestamp, y = Output), linetype = 'dashed') 
    }  
  }
  ## Turn into an interactive plotly (html plot)
  if(interactive){gg = ggplotly(gg)}
  
  return(gg + guides(col = legend))
}

plot_animate_clust = function(pred_gp, ygrid, data = NULL, data_train = NULL, mean_k = NULL, file = "gganim.gif")
{ ## pred_gp : tibble coming out of the pred_gp_animate() function, columns required : 'Timestamp', 'Mean', 'Var'
  ## data : tibble of observational data, columns required : 'Timestamp', 'Output' (Optional)
  ####
  ## return : plot the animated curves of the GP with the 0.95 confidence interval (optional display raw data)
  
  gg = plot_heat_clust(pred_gp, ygrid, data, data_train, mean_k, col_clust = F, animate = T) +
    transition_states(Nb_data, transition_length = 2, state_length = 1)
  animate(gg, renderer = gifski_renderer(file)) %>% return()
}

################ SIMULATED DATA ######################
simu_indiv = function(ID, t, kern = kernel_mu, k, theta, mean, var)
{ ## Return Age and Performance of a simulated individual
  ## nb : Number of timestamps
  ## Id : identifier of the individual
  ## tmin : Value of the first timestamp
  ## tmax : Value of the last timestamp
  ## a : Multiplying coefficient
  ## b : Intercept
  ## var : Variance of the additive noise
  # t = seq(tmin, tmax, length.out = nb)
  # indiv = tibble('ID' = rep(as.character(ID), nb), 'Timestamp' = t,
  #                'Output' = a*t + b + rnorm(nb,0,var),'a' = rep(runif(1,0,5) %>% round(2), nb),
  #                'b' = rep(runif(1,0,0.5) %>% round(2), nb))
  
  db = tibble('ID' = as.character(ID),
              'Timestamp' = t,
              'Output' = rmvnorm(1, rep(mean,length(t)), kern_to_cov(t, kern, theta, sigma = var)) %>% as.vector(),
              'Cluster' = as.character(k))
  return(db)
}

#set.seed(42)

M = 100
N = 20
t = matrix(0, ncol = N, nrow = M)
t_i = sample(seq(10, 20, 0.01),N , replace = F) 
for(i in 1:M){t[i,] = t_i %>% sort()}

db_train = c()#simu_indiv(ID = '1', t[1,], kernel_mu, theta = c(2,1), mean = 0, var = 0.5)
for(i in 1:M)
{
  k = i %% 3
  if(k == 0){db_train = rbind(db_train, simu_indiv(ID = i, t[i,], kernel_mu, k = 1, theta = c(2,1), mean = -20, var = 0.1))}
  if(k == 1){db_train = rbind(db_train, simu_indiv(ID = i, t[i,], kernel_mu, k = 2, theta = c(2,1), mean = 0, var = 0.1))}
  if(k == 2){db_train = rbind(db_train, simu_indiv(ID = i, t[i,], kernel_mu, k = 3, theta = c(2,1), mean = 20, var = 0.1))}
  #if(k == 3){db_train = rbind(db_train, simu_indiv(ID = i, t[i,], kernel_mu, k = 4, theta = c(1,2), mean = 35, var = 0.4))}
  #if(k == 4){db_train = rbind(db_train, simu_indiv(ID = i, t[i,], kernel_mu, k = 5,theta = c(1,1), mean = 55, var = 0.5))}
}

db_obs_test = simu_indiv(ID = M+1, sample(seq(10, 20, 0.01), N, replace = F) %>%
                           sort(), kernel_mu, k = 3, theta = c(2,1), mean = 0, var = 0.5)

################ INITIALISATION ######################
#ini_hp = list('theta_k' = c(1,1), 'theta_i' = c(1, 1, 0.2))

ini_kmeans = function(db, k, nstart = 50, summary = F)
{
  ## db: tibble containing common Timestamps and associated Output values to cluster
  ## k: number of clusters set in the kmeans
  ## nstart: number of initialisations to try in the kmeans
  ###
  ## return: tibble of ID, true clusters, and found clusters with kmeans
  if( !identical(unique(db$Timestamp), db %>% filter(ID == unique(db$ID)[[1]]) %>% pull(Timestamp)) )
  {
    floop = function(i)
    {
      obs_i = db %>% filter(ID == i) %>% pull(Output)
      tibble('ID' = i,
             'Timestamp' = seq_len(3) ,
             'Output' = c(min(obs_i), mean(obs_i) , max(obs_i)) ) %>% 
        return()
    }  
    db_regular = unique(db$ID) %>% lapply(floop) %>% bind_rows
  }
  else{db_regular = db}
  
  res = db_regular %>% dplyr::select(c(ID, Timestamp, Output)) %>% 
    spread(key =  Timestamp, value = Output) %>%
    dplyr::select(-ID) %>% 
    kmeans(centers = k, nstart = nstart)
  
  if(summary){res %>% print}
  
  broom::augment(res, db_regular %>% spread(key =  Timestamp, value = Output)) %>%
    dplyr::select(c(ID, .cluster)) %>% 
    rename(Cluster_ini = .cluster) %>% 
    mutate(Cluster_ini = paste0('K', .$Cluster_ini)) %>% 
    return()
}

ini_tau_i_k = function(db, k, nstart = 50)
{
  ini_kmeans(db, k, nstart) %>% mutate(value = 1) %>%
    spread(key = Cluster_ini, value = value, fill = 0) %>% 
    arrange(as.integer(ID)) %>% 
    column_to_rownames(var = "ID") %>%
    apply(2, as.list) %>% 
    return()
}


## TEST ####
# k = seq_len(3)
# tau_i_k_test = ini_tau_i_k(db = db_train, k = length(k), nstart = 50)
# 
# prior_mean_k = list('K1' = 0, 'K2' = 0, 'K3' = 0)
# ini_hp_test = list('theta_k' = c(2, 0.5, 0.2), 'theta_i' = c(1, 1, 0.2))
# common_hp_k = T
# common_hp_i = T


#  ## Full algo
# pred_test = full_algo_clust(db_train, db_obs_test[1:5,], seq(10, 20, 0.01), kernel, tau_i_k_test,
#                             common_hp_k, common_hp_i, prior_mean = prior_mean_k, kern_0 = kernel_mu, list_hp = NULL,
#                             mu_k = NULL, ini_hp = ini_hp_test, hp_new_i = NULL)
# plot_gp_clust(pred_test$Prediction, cluster = 'K1', data = db_obs_test, data_train = db_train,
#               mean = pred_test$Mean_processes, col_clust = T)
# ##Training
# t1 = Sys.time()
# training_test = training_VEM(db_train, prior_mean_k, ini_hp_test, kernel_mu, kernel,
#                                tau_i_k_test, common_hp_k, common_hp_i)
# t2 = Sys.time()
# c_time = (t2-t1) %>% print()
# ## Posterior mu_k
# timestamps = seq(10, 20, 0.01)
# post_test = posterior_mu_k(db_train, timestamps, prior_mean_k, kernel_mu, kernel, training_test)
# ## Training new HP_*
# new_indiv_hp = train_new_gp_EM(db_obs, post_test, ini_hp_test, kernel)
# ## Pred GP
# pred_gp_test = pred_gp_clust(db_obs[1:5,], timestamps, post_test, kern = kernel,
#              list('theta' = new_indiv_hp$theta_new[1:2], 'sigma' = new_indiv_hp$theta_new[[3]],
#                   'tau_k' = new_indiv_hp$tau_k ))
# # Plot GP
# plot_gp_clust(pred_gp_test, data = db_obs, data_train = db_train, mean = post_test$mean)
# # Plot GIF
# pred_gif = pred_gp_clust_animate(db_obs, timestamps, post_test, kern = kernel,
#                       list('theta' = training_test$theta_i$`1`[1:2], 'sigma' = training_test$theta_i$`1`[[3]],
#                            'tau_k' = pi_k_test ))
# plot_animate(pred_gif, data = db_obs, data_train = db_train, mean = post_test$mean, file = "bla.gif")

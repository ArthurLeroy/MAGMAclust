library(GPFDA)
library(parallel)
library(gridExtra)
library(fda)
library(funHDDC)
library(fossil)

source('MAGMAclust.R')

##### CLUST: COMPETING ALGO IN SIMU ####
train_gpfda = function(db)
{
  t1 = Sys.time()
  db = db %>% dplyr::select(ID, Timestamp, Output)
  M = db$ID %>% unique %>% length
  t_obs = db$Timestamp %>% unique
  N = t_obs %>% length
  
  lx = rep(1, M) %>% as.matrix
  
  fy1 = db %>% spread(key = ID, value = Output) %>% dplyr::select(- Timestamp ) %>% as.matrix %>% t()
  fx1 = db$Timestamp %>% matrix(nrow = N, ncol = M) %>% t()
  
  model = gpfr(response=(fy1), lReg=lx, fReg=NULL, gpReg=list(fx1), fyList=list(nbasis=23,lambda=0.1), fbetaList_l=NULL,
               hyper=NULL, Cov=c('pow.ex', 'linear'), fitting=T, time = t_obs, rPreIdx=T, concurrent=T) 
  t2 = Sys.time()
  model[['Time_train']] =  difftime(t2, t1, units = "secs")
  
  return(model)
}

pred_gpfda = function(model, new_db, timestamps)
{
  a1 = model
  tfx = timestamps %>% as.matrix
  time_new = timestamps
  
  pfy = new_db$Output
  pfx = new_db$Timestamp %>% as.matrix
  ptime = new_db$Timestamp
  
  b1<-gpfrpred(a1, TestData= tfx, NewTime=time_new, lReg = 1, fReg=NULL,
               gpReg=list('response'=(pfy),'input'=(pfx),'time'= ptime))
  
  tibble('Timestamp' = b1$predtime, 'Mean' = as.vector(b1$ypred.mean), 'Var' = as.vector(b1$ypred.sd)^2) %>% return()
}

mean_gpfda =  function(model, timestamps)
{
  a1 = model
  tfx = timestamps %>% as.matrix
  time_new = timestamps
  
  b1<-gpfrpred(a1, TestData= tfx, NewTime=time_new, lReg = 1, fReg=NULL,
               gpReg= NULL)
  
  tibble('Timestamp' = b1$predtime, 'Mean' = as.vector(b1$ypred.mean), 'Var' = as.vector(b1$ypred.sd)^2) %>% return()
}  

clust_kmeans = function(db, K, nstart = 50, summary = F)
{
  ## db: tibble containing common Timestamps and associated Output values to cluster
  ## k: number of clusters set in the kmeans
  ## nstart: number of initialisations to try in the kmeans
  ###
  ## return: tibble of ID, true clusters, and found clusters with kmeans
  #browser()
  basis = create.bspline.basis(c(0,10), nbasis=8)
  timestamps = unique(db$Timestamp)
  obs = db %>% dplyr::select(c(ID, Timestamp, Output)) %>% 
    mutate(ID = as.numeric(ID)) %>% 
    spread(key =  ID, value = Output) %>%
    column_to_rownames('Timestamp') %>% 
    as.matrix
  smoothed_db = smooth.basis(argvals = timestamps , y = obs, fdParobj = basis)$fd$coefs %>% t()
  
  res = smoothed_db %>% kmeans(centers = K, nstart = nstart)
  
  if(summary){ res %>% print }
  
  # broom::augment(res, db %>% spread(key =  Timestamp, value = Output)) %>%
  # dplyr::select(c(ID, Cluster, .cluster)) %>% 
  # rename(Cluster_found = .cluster) %>%
  # mutate_at('Cluster_found', as.integer) %>% 
  # arrange(as.integer(ID)) %>% 
  db %>% distinct(ID, Cluster) %>% 
    mutate(Cluster_found = as.vector(res$cluster)) %>% 
    return()
}

clust_funHDDC = function(db, K)
{
  basis = create.bspline.basis(c(0,10), nbasis=10)
  timestamps = unique(db$Timestamp)
  obs = db %>% dplyr::select(c(ID, Timestamp, Output)) %>% 
    mutate(ID = as.numeric(ID)) %>% 
    spread(key =  ID, value = Output) %>%
    column_to_rownames('Timestamp') %>% 
    as.matrix
  fdata = smooth.basis(argvals = timestamps , y = obs, fdParobj = basis)$fd
  clust = funHDDC(fdata, K = K, model = "all", init = "kmeans", threshold = 0.2)
  
  db %>% dplyr::select(ID, Cluster) %>% 
    unique %>% 
    mutate('Cluster_found' = clust$class) %>% 
    return()
}

##### CLUST: SIMULATION FUNCTIONS #####
simu_indiv = function(ID, t, mean, kern, clust, a, b, sigma)
{ # ID : identification of the individual
  # t : timestamps on which we observe the GP
  # kern : kernel associated to the covariance function of the GP
  # theta : list of hp of the kernel
  # mean : mean parameter of the GP
  # var : variance parameter of the error
  ##
  # return : a simulated individual
  db = tibble('ID' = ID,
              'Timestamp' = t,
              'Output' = rmvnorm(1, mean, kern_to_cov(t, kern, theta = c(a,b), sigma)) %>% as.vector(),
              'Cluster' = clust,
              'a' = a,
              'b' = b,
              'sigma' = sigma)
  return(db)
}

draw = function(int)
{
  return(runif(1,int[1], int[2]) %>% round(2))
}

prior_mean = function(t)
{
  a = draw(c(-2, 2))
  b = draw(c(20,30))
  return(a * t + b)
}

simu_scheme = function(M = 51, N = 30, K= 3, G = seq(0, 10, 0.05), pi = c(0.34,0.33,0.33),
                       common_times = T, common_hp_k = T, common_hp_i = T,
                       kern_0 = kernel_mu, kern_i = kernel,
                       int_mu_a = c(0,5),
                       int_mu_b = c(0,2),
                       int_i_a = c(0,5),
                       int_i_b = c(0,2),
                       int_i_sigma = c(0,0.1))
{
  if(common_hp_k)
  {
    k_a = draw(int_mu_a)
    k_b = draw(int_mu_b)
  }
  floopk = function(k)
  { 
    if(!common_hp_k)
    {
      k_a = draw(int_mu_a)
      k_b = draw(int_mu_b)
    }
    m_k = prior_mean(G)
    simu_indiv(paste0('K',k), G, m_k, kern_i, paste0('K',k), k_a, k_b, sigma = 0) %>% return()
  }
  db_k = lapply(seq_len(K), floopk) %>% bind_rows %>% as_tibble
  
  if(common_times){t_i = sample(G, N, replace = F) %>% sort()}
  if(common_hp_i)
  {
    i_a = draw(int_i_a)
    i_b = draw(int_i_b)
    i_sigma = draw(int_i_sigma)
  }
  
  floopi = function(i)
  { 
    if(!common_times){t_i = sample(G, N, replace = F) %>% sort()}
    if(!common_hp_i)
    {
      i_a = draw(int_i_a)
      i_b = draw(int_i_b)
      i_sigma = draw(int_i_sigma)
    }
    k = sample(seq_len(K), size=1, prob = pi)
    mean_i = db_k %>% filter(Cluster == paste0('K',k)) %>% filter(Timestamp %in% t_i) %>% pull(Output)
    
    simu_indiv(as.character(i), t_i, mean_i, kern_i, paste0('K',k), i_a, i_b, i_sigma) %>% return()
  }
  db_i = lapply(seq_len(M), floopi) %>% bind_rows %>% as_tibble
  return( list('db_i' = db_i, 'db_k' = db_k) )
}

##### CLUST: EVALUATION FUNCTIONS #####

MSE_clust = function(obs, pred)
{
  ## obs : the true observed values
  ## pred : list of probabilities et the parameters of the mixture, coming out from pred_gp_clust()
  ####
  ## return : the RMSE given the vector of errors
  t = obs %>% pull(Timestamp)
  value = obs %>% pull(Output)
  floop = function(k)
  {
    pred[[k]]$tau_k[[1]] * (pred[[k]] %>% filter(Timestamp %in% t) %>% pull(Mean)) %>%  return()
  }
  mix_pred = sapply(names(pred), floop) %>% rowSums
  (value - mix_pred)^2 %>% mean %>% return()
}

loss = function(x, y)
{ ## return : loss function between x and y
  abs(x - y) %>% return()
}

WCIC = function(obs, pred, level)
{ ## obs : the true observed values
  ## pred : list of probabilities et the parameters of the mixture, coming out from pred_gp_clust()
  ## level : confidence level (% of accepted error)
  ####
  ## return : the ratio of observed values lying within the mixture of predicted ICs
  t = obs %>% pull(Timestamp)
  value = obs %>% pull(Output)
  floop = function(k)
  {
    mean = pred[[k]] %>% filter(Timestamp %in% t) %>% pull(Mean)
    sd = pred[[k]] %>% filter(Timestamp %in% t) %>% pull(Var) %>% sqrt
    CI_inf = mean - qnorm(1 - level/2) * sd
    CI_sup = mean + qnorm(1 - level/2) * sd
    100 * pred[[k]]$tau_k[[1]] * ((CI_inf < value) & (value < CI_sup)) %>% mean %>% return()
  }
  sapply(names(pred), floop) %>% sum %>% return()
}

RI = function(group1, group2)
{
  g1 = group1 %>% as.factor %>% as.integer
  g2 = group2 %>% as.factor %>% as.integer
  
  # if( n_distinct(g1) == n_distinct(g2) ){rand.index(g1, g2) %>% return()}
  # else{adj.rand.index(g1, g2) %>% return()}
  adj.rand.index(g1, g2) %>% return()
}

MSE = function(error)
{ ## return : the RMSE given the vector of errors
  mean(error^2) %>% return()
}

ratio_IC = function(obs, IC_inf, IC_sup)
{ ## obs : the true observed values
  ## IC_inf : inferior boundary of the predicted IC_0.95
  ## IC_sup : superior boundary of the predicted IC_0.95
  ####
  ## return : the ratio of observed values lying within the predicted IC_0.95
  nb_between = ((IC_inf < obs) & (obs < IC_sup)) %>% sum()
  nb_between/ length(obs) * 100 %>% return()
}

eval_methods = function(db_results, test)
{ ## db_results : list of results of the different methods. Format : list('algo','one_gp','gpfda')
  ## db_test : vector of observed values on which we test the predictions
  ####
  ## return : Table of results of the evaluation of the different methods through RMSE and ratio_IC
  pred_clust = db_results$MAGMAclust
  
  pred_algo = db_results$MAGMA$Mean 
  sd_algo = db_results$MAGMA$Var %>% sqrt()
  
  pred_one_gp = db_results$GP$Mean
  sd_one_gp = db_results$GP$Var %>% sqrt()
  
  pred_gpfda = db_results$GPFDA$Mean
  sd_gpfda = db_results$GPFDA$Var %>% sqrt()
  
  db_test = test %>% pull(Output)
  
  eval_clust = tibble('MSE' = test %>%  MSE_clust(pred_clust),
                      'WCIC' = test %>% WCIC(pred_clust, 0.05),
                      'Time_train' = db_results$Time_train_magmaclust, 
                      'Time_pred' = db_results$Time_pred_magmaclust)
  
  eval_algo = tibble('MSE' = loss(pred_algo, db_test) %>% MSE(),
                     'WCIC' = ratio_IC(db_test, pred_algo - 1.96 * sd_algo, pred_algo + 1.96 * sd_algo),
                     'Time_train' = db_results$Time_train_algo, 'Time_pred' = db_results$Time_pred_algo)
  eval_one_gp = tibble('MSE' = loss(pred_one_gp, db_test) %>% MSE(),
                       'WCIC' = ratio_IC(db_test, pred_one_gp - 1.96 * sd_one_gp, pred_one_gp + 1.96 * sd_one_gp),
                       'Time_train' = 0, 'Time_pred' = db_results$Time_pred_one_gp)
  
  eval_gpfda = tibble('MSE' = loss(pred_gpfda, db_test) %>% MSE(),
                      'WCIC' = ratio_IC(db_test, pred_gpfda - 1.96 * sd_gpfda, pred_gpfda + 1.96 * sd_gpfda),
                      'Time_train' = db_results$Time_train_gpfda, 'Time_pred' = db_results$Time_pred_gpfda)
  
  rbind(eval_one_gp, eval_gpfda, eval_algo, eval_clust) %>% 
    mutate(Method = c('GP', 'GPFDA', 'MAGMA','MAGMAclust')) %>%
    return()
}

eval_clust = function(train_clust, average = F)
{
  floop = function(i)
  {
    res_clust = train_clust[[i]]
    eval_kmeans = tibble('RI' = RI(res_clust$k_means$Cluster, res_clust$k_means$Cluster_found), 
                         Method = 'kmeans + Bsplines')
    eval_funHDDC = tibble('RI' = RI(res_clust$funHDDC$Cluster, res_clust$funHDDC$Cluster_found),
                          Method = 'funHDDC')
    eval_magmaclust = tibble('RI' = RI(res_clust$MAGMAclust$Cluster, res_clust$MAGMAclust$Cluster_found),
                             Method = 'MAGMAclust')
    
    rbind(eval_kmeans, eval_funHDDC, eval_magmaclust) %>%
      return()
  }
  list_eval = seq_len(100) %>% lapply(floop)
  
  if(average)
  {
    do.call('rbind', list_eval) %>% 
      group_by(Method) %>% 
      summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% 
      return()
  }
  else{do.call('rbind', list_eval) %>% return()}
}

eval_clust_diffk = function(db, train_clust, average = F)
{
  floop = function(i, k)
  {
    res_clust = db %>% filter(ID_dataset == i) %>% 
      filter(ID != 1) %>% 
      dplyr::select(ID, Cluster) %>%
      unique %>% 
      mutate('Cluster_found' = train_clust[[paste0("K=", k)]][[i]]$tau_i_k %>% pred_max_cluster)
    tibble('RI' = RI(res_clust$Cluster, res_clust$Cluster_found), Method = 'MAGMAclust') %>% 
      return()
  }
  
  floop2 = function(k)
  {
    list_eval = seq_len(100) %>% lapply(floop, k = k)
    if(average)
    {
      do.call('rbind', list_eval) %>% 
        group_by(Method) %>% 
        summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% 
        mutate('K' = k) %>% 
        return()
    }
    else{do.call('rbind', list_eval) %>% mutate('K' = k) %>% return()}
  }
  lapply(2:10, floop2) %>% bind_rows %>% return()
}

eval_BIC = function(train_clust, table = T, recompute = F)
{
  floop = function(k)
  {
    k_hat = train_clust[[k]] %>% map_dbl('K_max_BIC')
    
    tibble('K_true' = as.numeric(k), 'K' = k_hat)
  }
  res = train_clust %>%
    names %>%
    lapply(floop) %>%
    bind_rows 
  
  if(table){res %>% table %>% return()}
  else(res %>% return())
}

recompute_BIC = function(db, train_clust, kern_0, kern_i,  plot = T)
{
  train_clust$K_max_BIC = NULL
  train_clust$plot = NULL
  floop = function(K)
  {
    print(K)
    model = train_clust[[K]]
    
    tibble('K' = K, 'BIC' = BIC(model$hp, db, kern_0, kern_i, 0, model$param, !(K == 'K = 1'))) %>% 
      return()
  }
  res = train_clust %>% names %>% lapply(floop) %>% bind_rows
  
  if(plot)
  {
    (res %>%
       mutate(K = factor(K, levels = K)) %>% 
       ggplot(aes(x = K, y = BIC)) + geom_point() + theme_classic()) %>% 
      print()  
  }
  return(res)
}

loop_recompute_BIC = function(db, model_train, kern_0, kern_i)
{  
  model_train$Time_train_tot = NULL
  ## Loop over the different datasets
  floop = function(i, k)
  {
    #browser()
    print(paste0('Dataset n°', i))
    
    ## Select the i-th dataset and remove mean process and testing individual (ID = 0 and 1)
    db %>% 
      filter(Nb_Cluster == k) %>% 
      filter(ID_dataset == i) %>%
      filter(!(ID %in% c(0,1))) %>%
      dplyr::select('ID', 'Timestamp', 'Output', 'Cluster') %>% 
      recompute_BIC(model_train[[k]][[i]], kern_0, kern_i, plot = F) %>% 
      mutate(ID = i) %>% 
      return()
  }
  
  floop2 = function(k)
  {
    print(paste0('True K value:', k))
    
    model_train[[k]] %>% 
      names %>% 
      sapply(floop, k, simplify = FALSE, USE.NAMES = TRUE) %>% 
      bind_rows %>% 
      mutate(Cluster_true = k) %>% 
      return()
  }
  model_train %>% 
    names %>% 
    sapply(floop2, simplify = FALSE, USE.NAMES = TRUE) %>% 
    bind_rows %>% 
    return()
}

##### CLUST: TRAINING AND PRED FUNCTIONS #####

training_diff_k = function(db, kmax, ini_hp, kern_0, kern_i, common_hp_k = T, common_hp_i = T)
{
  #browser()
  ## Loop over the different datasets
  floop = function(i, k = NULL)
  {
    print(paste0('Dataset n°', i, ' || k = ', k))
    
    ## Select the i-th dataset and remove mean process and testing individual (ID = 0 and 1)
    db_train = db %>% filter(ID_dataset == i) %>% 
      filter(ID != 1) %>%
      dplyr::select('ID', 'Timestamp', 'Output')
    
    prior_mean_k = rep(0, k) %>% setNames(paste0('K', seq_len(k))) %>% as.list
    
    train = training_VEM(db_train, prior_mean_k, ini_hp, kern_0, kern_i, ini_tau_i_k = NULL,
                         common_hp_k, common_hp_i)
    list('hp' = train$hp, 'Time_train' = train$Time_train, 'tau_i_k' = train$param$tau_i_k) %>% return()
  }
  
  ## Parallel computing through mclapply not available in windows
  nb_core = ifelse(.Platform$OS.type == 'windows', 1, 1)
  
  list_train <- mclapply(2:kmax, function(j) {
    unique(db$ID_dataset) %>% as.character() %>% 
      sapply(floop, k = j,  simplify = FALSE, USE.NAMES = TRUE) %>% 
      return()
  }, mc.cores= nb_core)
  names(list_train) = paste0('K=', 2:kmax)
  return(list_train)
}

loop_training_for_pred = function(db_loop, k, prior_mean, ini_hp_clust, ini_hp, kern_0, kern_i,
                                  common_hp_k = T, common_hp_i = T, common_times = T)
{  
  ## Loop over the different datasets
  floop = function(i)
  {
    print(paste0('Dataset n°', i))
    
    ## Select the i-th dataset and remove mean process and testing individual (ID = 0 and 1)
    db_train = db_loop %>% filter(ID_dataset == i) %>%
      filter(!(ID %in% c(0,1))) %>%
      dplyr::select('ID', 'Timestamp', 'Output')
    prior_mean_k = rep(prior_mean, k) %>% setNames(paste0('K', seq_len(k))) %>% as.list
    
    if(common_times){model_gpfda = train_gpfda(db_train)}else{model_gpfda = list('Time_train' = 0)}
    model_magma = training(db_train, prior_mean, ini_hp, kern_0, kern_i, common_hp_i)[c('hp', 'Time_train')]
    train = training_VEM(db_train, prior_mean_k, ini_hp_clust, kern_0, kern_i, ini_tau_i_k = NULL,
                         common_hp_k, common_hp_i)
    
    list('GPFDA' = model_gpfda, 'MAGMA' = model_magma,
         'MAGMAclust' = list('hp'=train$hp, 'Time_train'=train$Time_train, 'tau_i_k'=train$param$tau_i_k) ) %>%
      return()
  }
  
  list_train = unique(db_loop$ID_dataset) %>% as.character() %>% sapply(floop, simplify = FALSE, USE.NAMES = TRUE)
  
  list_train %>% c(list('prior_mean' = prior_mean, 'ini_hp_clust' = ini_hp_clust, 'ini_hp' = ini_hp,
                        'common_times' = common_times, 'common_hp_k' = common_hp_k, 'common_hp_i' = common_hp_i)) %>% 
    return()
}

loop_training_for_clust = function(db_loop, k, prior_mean, ini_hp_clust, kern_0, kern_i,
                                   common_hp_k = T, common_hp_i = T)
{  
  ## Loop over the different datasets
  floop = function(i)
  {
    print(paste0('Dataset n°', i))
    
    ## Select the i-th dataset and remove mean process and testing individual (ID = 0 and 1)
    db_train = db_loop %>% filter(ID_dataset == i) %>%
      filter(!(ID %in% c(0,1))) %>%
      dplyr::select('ID', 'Timestamp', 'Output', 'Cluster')
    prior_mean_k = rep(prior_mean, k) %>% setNames(paste0('K', seq_len(k))) %>% as.list
    
    clust_kmeans = clust_kmeans(db_train, k)
    clust_funHDDC = tryCatch(clust_funHDDC(db_train, k), error = function(e){clust_kmeans})
    train = training_VEM(db_train, prior_mean_k, ini_hp_clust, kern_0, kern_i, ini_tau_i_k = NULL,
                         common_hp_k, common_hp_i)
    clust_magmaclust = db_train %>% dplyr::select(ID, Cluster) %>%
      unique %>%
      mutate('Cluster_found' = pred_max_cluster(train$param$tau_i_k))
    list('MAGMAclust' = clust_magmaclust, 'funHDDC' = clust_funHDDC,
         'k_means' = clust_kmeans) %>%
      return()
  }
  
  list_train = unique(db_loop$ID_dataset) %>% as.character() %>% sapply(floop, simplify = FALSE, USE.NAMES = TRUE)
  
  list_train %>% c(list('prior_mean' = prior_mean, 'ini_hp_clust' = ini_hp_clust,
                        'common_hp_k' = common_hp_k, 'common_hp_i' = common_hp_i)) %>%
    return()
}

loop_pred = function(db_loop, train_loop, nb_obs, nb_test, k = 3)
{
  db_loop$ID = db_loop$ID %>% as.character
  ## Get the settings used for training
  prior_mean = train_loop$prior_mean
  prior_mean_k = rep(prior_mean, k) %>% setNames(paste0('K', seq_len(k))) %>% as.list
  ini_hp_clust = train_loop$ini_hp_clust
  ini_hp = train_loop$ini_hp
  kern_0 = kernel
  kern_i = kernel_mu
  common_times = train_loop$common_times
  common_hp_k = train_loop$common_hp_k
  common_hp_i = train_loop$common_hp_i
  common_hp = train_loop$common_hp_i
  
  floop = function(i)
  {
    print(paste('i =' ,i))
    ## Get the trained model for GPFDA and our algo
    model_algo = train_loop[[i]]$MAGMA
    list_hp = model_algo$hp
    model_gpfda = train_loop[[i]]$GPFDA
    model_magmaclust = train_loop[[i]]$MAGMAclust
    list_hp_clust = model_magmaclust
    list_hp_clust[['param']] = list('tau_i_k' = model_magmaclust$tau_i_k)
    
    ## Get the corresponding database
    db_i = db_loop %>% filter(ID_dataset == i)
    db_train_i = db_i %>% filter(ID != 1)
    ## Select the 'nb_obs' first observations of the testing individual to predict with
    db_obs_i = db_i %>% filter(ID == 1) %>% top_n(- nb_obs, Timestamp)
    ## Select the 'nb_test last observations of the testing individual to evaluate predictions on 
    db_pred_i = db_i %>% filter(ID == 1) %>% top_n(nb_test, Timestamp)
    ## Get timestamps to predict on
    t_i_pred = db_pred_i %>% pull(Timestamp)
    
    t1 = Sys.time()
    ## Prediction for our algo (train new indiv + pred)
    res_algo = full_algo(db_train_i, db_obs_i, t_i_pred, kern_i, common_hp, plot = F, prior_mean, kern_0,
                         list_hp, mu = NULL, ini_hp, hp_new_i = NULL)$Prediction
    t2 = Sys.time()
    ## Train new indiv for one GP model
    hp_one_gp = train_new_gp(db_obs_i, rep(prior_mean, nrow(db_obs_i)), cov_mu = 0, ini_hp$theta_i, kern_i)
    ## Prediction for one GP model
    res_one_gp = pred_gp(db_obs_i, t_i_pred, prior_mean, cov_mu = NULL, kern_i, hp_one_gp$theta, hp_one_gp$sigma)
    t3 = Sys.time()
    ## Prediction for GPFDA (unable for uncommon timestamps)
    if(common_times){res_gpfda = pred_gpfda(model_gpfda, db_obs_i, t_i_pred)}
    else{res_gpfda =  tibble('Timestamp' = t_i_pred, 'Mean' = NA ,   'Var' = NA)}
    t4 = Sys.time()
    
    res_magmaclust = full_algo_clust(db_train_i, db_obs_i, t_i_pred, kern_i, tau_i_k, common_hp_k, common_hp_i, 
                                     prior_mean_k, kern_0, list_hp_clust, mu_k = NULL, ini_hp_clust, hp_new_i = NULL)$Prediction
    t5 = Sys.time()
    
    ### Get MSE, RATIO IC95 and computing times on testing points for all methods 
    list('MAGMA' = res_algo, 'Time_train_algo' = model_algo$Time_train, 
         'Time_pred_algo' = difftime(t2, t1, units = "secs"),
         'GP' = res_one_gp, 'Time_pred_one_gp' =  difftime(t3, t2, units = "secs"),
         'GPFDA' = res_gpfda, 'Time_train_gpfda' = model_gpfda$Time_train,
         'Time_pred_gpfda' =  difftime(t4, t3, units = "secs"),
         'MAGMAclust' = res_magmaclust, 'Time_train_magmaclust' = model_magmaclust$Time_train, 
         'Time_pred_magmaclust' = difftime(t5, t4, units = "secs")) %>%
      eval_methods(db_pred_i) %>%
      return()
  }
  
  list_eval = db_loop$ID_dataset %>% unique() %>% lapply(floop)
  
  do.call('rbind', list_eval) %>% 
    mutate(Time_train = as.numeric(Time_train), Time_pred = as.numeric(Time_pred)) %>% 
    return()
}

loop_pred_diffk = function(db_loop, train_loop, nb_obs, nb_test, prior_mean = 0, kern_0 = kernel,
                           kern_i = kernel_mu, common_hp_k = T, common_hp_i = T)
{
  db_loop$ID = db_loop$ID %>% as.character
  ## Get the settings used for training
  floop = function(i, k)
  {
    print(paste('i =' ,i))
    
    model = train_loop[[paste0('K=',k)]][[i]]
    list_hp_clust = model
    list_hp_clust[['param']] = list('tau_i_k' = model$tau_i_k)
    
    prior_mean_k = rep(prior_mean, k) %>% setNames(paste0('K', seq_len(k))) %>% as.list
    
    
    ## Get the corresponding database
    db_i = db_loop %>% filter(ID_dataset == i)
    db_train_i = db_i %>% filter(ID != 1)
    ## Select the 'nb_obs' first observations of the testing individual to predict with
    db_obs_i = db_i %>% filter(ID == 1) %>% top_n(- nb_obs, Timestamp)
    ## Select the 'nb_test last observations of the testing individual to evaluate predictions on 
    db_pred_i = db_i %>% filter(ID == 1) %>% top_n(nb_test, Timestamp)
    ## Get timestamps to predict on
    t_i_pred = db_pred_i %>% pull(Timestamp)
    
    t1 = Sys.time()
    
    pred_clust = full_algo_clust(db_train_i, db_obs_i, t_i_pred, kern_i, ini_tau_i_k= NULL, common_hp_k, common_hp_i, 
                                 prior_mean_k, kern_0, list_hp_clust, mu_k = NULL, ini_hp_clust, hp_new_i = NULL)$Prediction
    t2 = Sys.time()
    
    ### Get MSE, RATIO IC95 and computing times on testing points for all methods 
    tibble('MSE' = db_pred_i %>%  MSE_clust(pred_clust), 'WCIC' = db_pred_i %>% WCIC(pred_clust, 0.05),
           'Time_train' = model$Time_train, 'Time_pred' = difftime(t2, t1, units = "secs"), 
           'K' = k) %>% 
      return()
  }
  
  floop2 = function(k)
  {
    db_loop$ID_dataset %>% unique %>% 
      lapply(floop, k) %>% 
      bind_rows %>% 
      mutate(K = k) %>% 
      mutate(Time_train = as.numeric(Time_train), Time_pred = as.numeric(Time_pred)) %>% 
      return()
  }
  lapply(2:10, floop2) %>% bind_rows %>% return()
}

loop_training_for_BIC = function(db_loop, k_grid, ini_hp_clust, kern_0, kern_i, common_hp_k, common_hp_i)
{  
  ## Loop over the different datasets
  floop = function(i, k)
  {
    #browser()
    print(paste0('Dataset n°', i, ' | True K value:', k))
    
    ## Select the i-th dataset and remove mean process and testing individual (ID = 0 and 1)
    db = db_loop %>% 
      filter(Nb_Cluster == k) %>% 
      filter(ID_dataset == i) %>%
      filter(!(ID %in% c(0,1))) %>%
      dplyr::select('ID', 'Timestamp', 'Output', 'Cluster')
    
    res = tryCatch(model_selection(db, k_grid, ini_hp_clust, kern_0, kern_i,
                                   ini_tau_i_k = NULL, common_hp_k, common_hp_i,
                                   plot = F), 
                   error = function(e){list('K_max_BIC' = NA)})
    return(res)
  }
  
  floop2 = function(k)
  {
    print(paste0('True K value:', k))
    
    db_loop %>%
      filter(Nb_Cluster == k) %>%
      pull(ID_dataset) %>% 
      unique %>% 
      as.character %>% 
      sapply(floop, k, simplify = FALSE, USE.NAMES = TRUE) %>% 
      return()
  }
  db_loop %>% 
    pull(Nb_Cluster) %>%
    unique %>% 
    as.character %>% 
    sapply(floop2, simplify = FALSE, USE.NAMES = TRUE) %>% 
    return()
}

##### CLUST: SIMULATION DATASETS ####

datasets_multi = function(rep, M, N, K, G, common_times, common_hp_i, common_hp_k, kern_0, kern_i,
                          int_mu_a, int_mu_b, int_i_a, int_i_b, int_i_sigma, k_grid = NULL)
{ ## rep : number of dataset to draw
  ## other inputs : same as simu_scheme function
  ####  
  ## return : tibble of rep x length(M) binded datasets. Columns 'nb_M' and 'ID_dataset' are added to distinguish them
  floop = function(j, k)
  {
    pi = rep(1/k, k)
    simu_scheme(M, N, k, G, pi, common_times, common_hp_k, common_hp_i, kern_0, kern_i,
                int_mu_a, int_mu_b, int_i_a, int_i_b, int_i_sigma)$db_i %>% 
      mutate('ID_dataset' = as.character(j)) %>% 
      return()
  }
  floop2 = function(k)
  {
    seq_len(rep) %>%
      lapply(floop, k = k) %>%
      bind_rows %>% 
      mutate('Nb_Cluster' = as.character(k)) %>% 
      return()
  }
  
  if(is.null(k_grid))
  {
    seq_len(rep) %>%
      lapply(floop, k = K) %>%
      bind_rows %>%
      return()
  }
  else
  {
    k_grid %>% 
      lapply(floop2) %>%
      bind_rows %>% 
      return()
  }
}

# set.seed(42)
# db_simu = datasets_multi(rep = 50, M = 101, K=3, N = 30, G = seq(0, 10, 0.05), common_times = T,
#                common_hp_i = T, common_hp_k = T,  kern_0 = kernel_mu, kern_i = kernel, int_mu_a = c(0,3), int_mu_b = c(0,1),
#                int_i_a = c(0,3), int_i_b = c(0,1), int_i_sigma = c(0,0.1), k_grid = 1:5)
# db_simu %>%  write_csv('Simu_clust/Data/db_i_selec_50rep_M100_N30_time_Hoo.csv')

##### CLUST: TABLES OF DATA ####

# # Load the data and ensure IDs are filled as characters
# table_Hoo = read_csv("Simu_clust/Data/db_i_clust_100rep_M50_N30_time_Hoo.csv")
# table_Hoo$ID = as.character(table_Hoo$ID)
# table_Hoo$ID_dataset = as.character(table_Hoo$ID_dataset)

# table_Hoo_diff_t = read_csv("Simu_clust/Data/db_i_clust_100rep_M20_N30_diff_time_Hoo.csv")
# table_Hoo_diff_t$ID = as.character(table_Hoo_diff_t$ID)
# table_Hoo_diff_t$ID_dataset = as.character(table_Hoo_diff_t$ID_dataset)

# table_Hoi = read_csv("Simu_clust/Data/db_i_clust_100rep_M50_N30_time_Hoi.csv")
# table_Hoi$ID = as.character(table_Hoi$ID)
# table_Hoi$ID_dataset = as.character(table_Hoi$ID_dataset)
# 
# table_Hko = read_csv("Simu_clust/Data/db_i_clust_100rep_M50_N30_time_Hko.csv")
# table_Hko$ID = as.character(table_Hko$ID)
# table_Hko$ID_dataset = as.character(table_Hko$ID_dataset)
# 

# table_Hki = read_csv("Simu_clust/Data/db_i_clust_100rep_M50_N30_time_Hki.csv")
# table_Hki$ID = as.character(table_Hki$ID)
# table_Hki$ID_dataset = as.character(table_Hki$ID_dataset)

# table_Hoo_selec = read_csv("Simu_clust/Data/db_i_selec_50rep_M100_N30_time_Hoo.csv")
# table_Hoo_selec$ID = as.character(table_Hoo_selec$ID)
# table_Hoo_selec$ID_dataset = as.character(table_Hoo_selec$ID_dataset)


##### CLUST: TRAIN ALL MODEL ####

# db_to_train = table_Hoo_selec
# t1 = Sys.time()
#train_loop = training_diff_k(db_to_train, kmax = 10, ini_hp =  list('theta_k' = c(1,1,0.2), 'theta_i' = c(1,1,0.2)),
#                             kern_0 = kernel_mu, kern_i = kernel, common_hp_k = T, common_hp_i = T)
# train_loop = loop_training_for_BIC(db_to_train, k_grid = 1:6,
#                                    ini_hp_clust = list('theta_k' = c(1,1,0.2), 'theta_i' = c(1,1,0.2)),
#                                    kern_0 = kernel_mu, kern_i = kernel,
#                                     common_hp_k = T, common_hp_i = T)

# train_loop = loop_training_for_pred(db_to_train, k = 3, prior_mean = 0,
#                                     ini_hp_clust = list('theta_k' = c(1,1,0.2), 'theta_i' = c(1,1,0.2)),
#                                     ini_hp = list('theta_0' = c(1,1), 'theta_i' = c(1,1,0.2)),
#                                     kern_0 = kernel_mu, kern_i = kernel,
#                                     common_hp_k = T, common_hp_i = T, common_times = T)

# train_loop = loop_training_for_clust(db_to_train, k = 3, prior_mean = 0,
#                                     ini_hp_clust = list('theta_k' = c(1,1,0.2), 'theta_i' = c(1,1,0.2)),
#                                     kern_0 = kernel_mu, kern_i = kernel,
#                                     common_hp_k = F, common_hp_i = F)
# t2 = Sys.time()
# train_loop[['Time_train_tot']] = difftime(t2, t1, units = "mins")
# saveRDS(train_loop, 'Simu_clust/Training/train_for_selec_Hoo_M100.rds')

##### CLUST: RESULTS : evaluation of clustering diff K ####
# model_clust = readRDS('Simu_clust/Training/train_diffk_Hoo_M50.rds')
# res_clust = eval_clust_diffk(table_Hoo, model_clust, F)
# write.csv(res_clust, "Simu_clust/Results/res_diffk_clust_M50.csv", row.names=FALSE)

# model_pred = readRDS('Simu_clust/Training/train_diffk_Hoo_M50.rds')
# res_pred = loop_pred_diffk(table_Hoo, model_pred, nb_obs = 20, nb_test = 10, prior_mean = 0, kern_0 = kernel,
#                            kern_i = kernel_mu, common_hp_k = T, common_hp_i = T)
# write.csv(res_pred, "Simu_clust/Results/res_diffk_pred_M50.csv", row.names=FALSE)

# res_pred = read_csv('Simu_clust/Results/res_diffk_pred_M50.csv')
# res_pred %>% group_by(K) %>% summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)
# ggplot(res_clust) + geom_boxplot(aes(x = as.factor(K), y = RI)) #+ scale_y_continuous(limits = c(0,100))


##### CLUST: RESULTS : evaluation of clustering vs alternatives ####

# model_clust = readRDS('Simu_clust/Training/train_for_clust_Hki_M50.rds')
# res_clust = eval_clust(model_clust, F)
# write.csv(res_clust, "Simu_clust/Results/res_clust_Hki_M50.csv", row.names=FALSE)

# res_clust = read_csv('Simu_clust/Results/res_diffk_pred.csv')
# res_clust %>% group_by(K) %>% summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)
# ggplot(res_clust) + geom_boxplot(aes(x = as.factor(K), y = MSE)) #+ scale_y_continuous(limits = c(0,100))

##### CLUST: RESULTS : evaluation of pred vs alternatives ####

# model_pred = readRDS('Simu_clust/Training/train_for_pred_Hoo_M50.rds')
# res_pred = loop_pred(table_Hoo, model_pred, nb_obs = 20, nb_test = 10)
# write.csv(res_pred, "Simu_clust/Results/res_pred_Hoo_M50.csv", row.names=FALSE)

# res_pred = read_csv('Simu_clust/Results/res_pred_Hoo_M50.csv')
# res_pred %>% group_by(Method) %>% summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)
# ggplot(res_clust) + geom_boxplot(aes(x = as.factor(K), y = RI)) #+ scale_y_continuous(limits = c(0,100))


##### CLUST: RESULTS : evaluation of the model selection ####
# model_selec = readRDS('Simu_clust/Training/train_for_selec_Hoo_M100.rds')
# new_BIC_simu = loop_recompute_BIC(table_Hoo_selec, model_selec, kernel, kernel)
# new_BIC = new_BIC_simu %>% mutate(K = str_sub(K, start = 5, end = 5) %>% as.numeric)
# max_new = new_BIC %>% group_by(Cluster_true, ID) %>% summarize(max = which.max(BIC))
# table(max_new$max, max_new$Cluster_true) %>% t()
# model_selec$Time_train_tot = NULL
# res = eval_BIC(model_selec, table = T)

##### CLUST: PLOT OF RESULTS #### 

### Boxplots changing values of K

# res_clust =  read_csv('Simu_clust/Results/res_diffk_clust_M50.csv')
# res_pred = read_csv('Simu_clust/Results/res_diffk_pred_M50.csv')

# res_clust_plot = res_clust %>% mutate(True_clust = ifelse(K == 3, 'Correct K','Incorrect K'))
# plot1 = ggplot(res_clust_plot, aes(x = as.factor(K), y = RI, fill = True_clust)) +
#         geom_boxplot(outlier.shape = NA) + xlab('K') + ylab('ARI') +
#         theme_classic() + scale_x_discrete(labels=c('3'= '3*' )) +
#         scale_fill_manual(name='', values = c("#619CFF", "#F8766D"))

# plot2 = ggplot(res_pred) + geom_boxplot(aes(x = as.factor(K), y = MSE)) + xlab('K') + theme_classic() +
#   scale_x_discrete(labels=c('3'= parse(text = TeX('$\\mathbf{3^*}$'))))

# png("diff_K.png",res=600, height=120, width= 176, units="mm")
# plot1
# dev.off()

### Boxplots comparison in RI vs alternatives

# res_plot_Hoo =  read_csv('Simu_clust/Results/res_clust_Hoo_M50.csv') %>% mutate(Hypothesis = 'Hoo')
# res_plot_Hko =  read_csv('Simu_clust/Results/res_clust_Hko_M50.csv') %>% mutate(Hypothesis = 'Hko')
# res_plot_Hoi =  read_csv('Simu_clust/Results/res_clust_Hoi_M50.csv') %>% mutate(Hypothesis = 'Hoi')
# res_plot_Hki =  read_csv('Simu_clust/Results/res_clust_Hki_M50.csv') %>% mutate(Hypothesis = 'Hki')
# res_plot = res_plot_Hoo %>% bind_rows(res_plot_Hko, res_plot_Hoi, res_plot_Hki)
# 
# res_plot %>% dplyr::select(-Hypothesis) %>% group_by(Method) %>% summarize_all(list(Mean = mean, SD = sd), na.rm = T)
# 
# png("RI_vs_alternatives.png",res=600, height=120, width= 176, units="mm")
# ggplot(res_plot) + geom_boxplot(aes(x = Hypothesis, y = RI, fill = Method), outlier.shape = NA) + ylab('ARI') + theme_classic()
# dev.off()

### Tab prediction vs alternatives

# res_plot = read_csv('Simu_clust/Results/res_pred_Hoo.csv')
# 
# res_plot %>% group_by(Method) %>% summarize_all(list(Mean = mean, SD = sd)) 
# 
# ggplot(res_plot) + geom_boxplot(aes(x = Method, y = MSE))+ theme_classic() +coord_cartesian(ylim = c(0,160))

##### TEST GPFDA ####
# M = 20
# N = 11
# t = matrix(0, ncol = N, nrow = M)
# for(i in 1:M){t[i,] = seq(0,10, length.out = N)}
# 
# db_train = simu_indiv(ID = '1', t[1,], kernel_mu, theta = c(2,1), mean = 45, var = 0.2)
# for(i in 2:M)
# {
#   k = i %% 5
#   if(k == 0){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(2,2), mean = 45, var = 0.6))}
#   if(k == 1){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(2,1), mean = 45, var = 0.2))}
#   if(k == 2){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(3,2), mean = 45, var = 0.3))}
#   if(k == 3){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(1,2), mean = 45, var = 0.4))}
#   if(k == 4){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(1,1), mean = 45, var = 0.5))}
# }
# db_obs = simu_indiv(ID = (M+1) %>% as.character(), seq(0,10, length.out = N),
#                     kernel_mu, theta = c(2,1), mean = 40, var = 0.2)
# 
# 
# plot(b1,type='prediction')
# 
# plot(-1000,col=0,xlim=range(b1$time),ylim=range(b1$ypred),xlab='time',ylab='prediction',
#      main='Prediction by GPFR: type I')
# 
# lines(b1$predtime,b1$ypred[,1])
# lines(b1$predtime,b1$ypred[,2],lty=2,col=2)
# lines(b1$predtime,b1$ypred[,3],lty=2,col=2)
# points(xt,yt)
# 

##### CLUST: ILLUSTRATION EXAMPLE ####
both_db = simu_scheme()
db_i = both_db$db_i
db_k = both_db$db_k
plot_db(db_k)
plot_db(db_i, cluster= T)
db_obs = db_i %>% filter(ID == 1)

### Simple situation and comparision with GP and MAGMA
pred = full_algo(db_i %>% filter(ID != 6), db_obs[1:14,], seq(0,11.5, 0.01), kernel, common_hp = T, plot = F, prior_mean = 0,
                 kernel_mu, list_hp = NULL , mu = NULL, ini_hp = ini_hp, hp_new_i = NULL)

ex_magma = plot_gp(pred$Prediction, data = db_obs[1:14,], data_train = NULL, mean = pred$Mean_process$pred_GP, mean_CI = F) +
  geom_point(data = db_i, aes(Timestamp, Output, color = ID),  size = 0.5, alpha = 0.5) +
  guides(color = FALSE) + geom_point(data = db_obs[15:30,], aes(Timestamp, Output), color ='red') +
  theme_classic()
# scale_y_continuous(limits = c(- 20, 32))

hp_one = train_new_gp(db_obs[1:14,], 0, 0, ini_hp$theta_i, kernel)
pred_one = pred_gp(db_obs[1:14,], timestamps = seq(0,11.5, 0.01), mean_mu = 0, cov_mu = NULL,
                   kern = kernel, theta = c(2,2), sigma = hp_one$sigma )
ex_gp = plot_gp(pred_one, data = db_obs[1:14,]) + geom_point(data = db_obs[15:30,],
                                                             aes(Timestamp, Output), color ='red') + theme_classic()

tau_i_k = replicate(length(k), rep(1,length(unique(db_i$ID[-1])))) %>%
  apply(1,function(x) x / sum(x)) %>%
  `rownames<-`(paste0('K', k)) %>%
  `colnames<-`(unique(db_i$ID[-1])) %>%
  apply(1, as.list)

pred_clust = full_algo_clust(db_i %>% filter(ID != 19), db_obs[1:14,], seq(0,11.5, 0.01), kernel, tau_i_k,
                             common_hp_k = T, common_hp_i = T, prior_mean = list('K1' = -5, 'K2' = 0, 'K3' = 20),
                             kernel_mu, list_hp = NULL , mu = NULL, ini_hp = ini_hp_test, hp_new_i = NULL)

ex_clust_all = plot_gp_clust(pred_clust, cluster= 'all', data = db_obs[1:14,],  mean = T, col_clust = T) +
  geom_point(data = db_i, aes(Timestamp, Output, color = Cluster), size = 0.5, alpha = 0.5) +
  guides(color = FALSE) + geom_point(data = db_obs[15:30,], aes(Timestamp, Output), color ='red')+
  theme_classic()

#  png("illu_2_other_clusters.png",res=600, height=120, width= 352, units="mm")
# # grid.arrange(ex_gp, ex_magma, ex_clust, nrow = 3)
#   grid.arrange(ex_clust1, ex_clust2, ncol = 2)
#  dev.off()

###  Fuzzy situation and heatmap 
pred_clust_2 = full_algo_clust(db_i %>% filter(ID != 1), db_obs[1:14,], seq(0,11.5, 0.01), kernel, tau_i_k,
                               common_hp_k = T, common_hp_i = T, prior_mean = list('K1' = 20 , 'K2' = 30, 'K3' = 60),
                               kernel_mu, list_hp = NULL , mu = NULL, ini_hp = ini_hp_test, hp_new_i = NULL)


heat_plot = plot_heat_clust(pred_clust, data = db_obs[1:14,], mean = T,
                            ygrid = seq(-10, 40, 0.1),  col_clust = T, interactive = F) +
  geom_point(data = db_i, aes(Timestamp, Output, color = Cluster), size = 0.5, alpha = 0.5) +
  guides(color = FALSE) + geom_point(data = db_obs[15:30,], aes(Timestamp, Output), color ='red')+
  theme_classic() + theme(legend.position=c(0.1,0.2))

png("illu_heatmap.png",res=600, height=120, width= 352, units="mm")
# grid.arrange(ex_gp, ex_magma, ex_clust, nrow = 3)
grid.arrange(ex_clust_all, heat_plot, ncol = 2)
dev.off()

### Evaluate clustering

gg2 = ggplot()
for(k in names(pred_clust_2$Prediction))
{
  gg2 = gg2+geom_line(data = pred_clust_2$Mean_processes$mean[[k]], aes(x = Timestamp, y = Output), linetype = 'dashed')
}  
gg2 = gg2 +geom_point(data = db_i, aes(Timestamp, Output, color = Cluster), size = 0.5, alpha = 0.5) +
  guides(color = FALSE) +
  theme_classic() + geom_line(data = db_k, aes(Timestamp, Output, color = Cluster))

png("shape_cluster.png",res=600, height=120, width= 352, units="mm")
grid.arrange(gg, gg2, ncol = 2)
dev.off()

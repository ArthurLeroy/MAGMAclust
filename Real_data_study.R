library(tidyverse)
library(compositions)
library(plotly)
library(png)

source('MAGMAclust.R')

#### SPLITING TRAINING AND TESTING FUNCTIONS ####
split_train = function(db, ratio_train)
{ ## db : the database of all observed data
  ## ratio_train : number between 0 and 1 indicating the proportion of individuals in the training set
  #
  ## return : orginal db with the repartition of individuals between the training and testing set
  n_indiv = db$ID%>% unique()
  n_index = (ratio_train * length(n_indiv)) %>% round()
  index_train = sample(n_indiv, n_index, replace = F)
  
  db %>% mutate(Training = ifelse(ID %in% index_train, 1,0)) %>% return()
}

split_times = function(db, prop_test)
{ ## db : the database of all observed data
  ## prop_test : number between 0 and 1 indicating the proportion of points we predict on
  #
  ## return : orginal db with the repartition of individuals between the training and testing set
  db = db %>% rowid_to_column("Row")
  row_obs = db %>% group_by(ID) %>% top_n(Timestamp, n = - n() * (1 - prop_test)) %>% pull(Row)
  
  db %>% ungroup() %>% mutate(Observed = ifelse(db$Row %in% row_obs, 1, 0)) %>% dplyr::select(- Row) %>% return()
}


#### EVALUATION FUNCTIONS ####
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



#### DATA IMPORT ####
data(juraset)
raw_db_co2 = read_csv('Real_Data_Study/Data/co2-data.csv')
raw_db_weight = read_csv('Real_Data_Study/Data/db_gusto_weight.csv') 
db_j = juraset %>%
  mutate(ID = paste(Rock, Land)) %>% 
  select(ID, X, Y, Zn)

# plot_ly(x=db_jura$X, y=db_jura$Y, z=db_jura$Zn, 
#         type="scatter3d", mode="markers", color=db_jura$ID)
####

#### Weight study ####
set.seed(42)
## Split by gender and draw the training/testing sets
db_w = raw_db_weight %>% 
  rename(Timestamp = Input) %>% 
  dplyr::select(-sex) %>%
  split_train(ratio_train = 0.6)

## Splitting training and testing sets and select observed times
db_w_train = db_w %>%
  filter(Training == 1)
db_w_test = db_w %>%
  filter(Training == 0) %>%
  split_times(prop_test = 0.3)

#mod_select = model_selection(db_w_train, k_grid = 2:10)
#saveRDS(mod_select,'Real_Data_Study/Training/train_weight_Hoo_model_selection.rds')
#mod_select = readRDS('Real_Data_Study/Training/train_weight_Hoo_model_selection.rds')

mod = mod_select$`K = 4`

new_db = db_w_test %>% filter(ID == '010-04050') %>% 
  filter(Observed == 1) %>% 
  dplyr::select(-c(Training, Observed))
new_db_t = db_w_test %>% filter(ID == '010-04050') %>% 
  filter(Observed == 0) %>% 
  dplyr::select(-c(Training, Observed))

timestamps = seq(-1, 100, 0.1)
m_k = list('K1' = 0, 'K2' = 0, 'K3' = 0, 'K4' = 0)
list_mu = posterior_mu_k(db_w_train, timestamps, m_k, kernel, kernel, mod)

new_hp = train_new_gp_EM(new_db, list_mu, NULL, kernel, hp_i = mod)

pred = pred_gp_clust(new_db, timestamps, list_mu, kernel, new_hp)
clust = paste0('K', pred_max_cluster(new_hp$tau_k))

png("weight_plot_pred4.png",res=600, height=120, width= 220, units="mm")
plot_gp_clust(pred, cluster = 'all', data = new_db,
              data_train = db_w_test, mean_k = list_mu) +
  geom_point(data = new_db_t, aes(x = Timestamp, y = Output), col = 'green')+ 
  theme_classic() 
dev.off()

png("weight_plot_pred_bestclust4.png",res=600, height=120, width= 220, units="mm")
plot_gp_clust(pred, cluster = clust, data = new_db,
              data_train = db_w_test, mean_k = list_mu) +
  geom_point(data = new_db_t, aes(x = Timestamp, y = Output), col = 'green')+ 
  theme_classic() 
dev.off()
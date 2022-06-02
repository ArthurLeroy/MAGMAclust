library(tidyverse)
library(png)


#### TRAIN/TEST SPLITING FUNCTIONS ####
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
  row_obs = db %>% group_by(ID) %>% top_n(Input, n = - n() * (1 - prop_test)) %>% pull(Row)
  
  db %>% ungroup() %>% mutate(Observed = ifelse(db$Row %in% row_obs, 1, 0)) %>% dplyr::select(- Row) %>% return()
}

#### EVALUATION FUNCTIONS ####
MSE_clust = function(obs, pred)
{
  ## obs : the true observed values
  ## pred : list of probabilities et the parameters of the mixture, coming out from pred_gp_clust()
  ####
  ## return : the RMSE given the vector of errors
  input = obs %>% pull(Input)
  value = obs %>% pull(Output)
  floop = function(k)
  {
    pred_k = pred$pred[[k]] %>%
      filter(Input %in% input) %>%
      pull(Mean)
    tau_k = pred$pred[[k]]$Proba %>% unique()
    
    (tau_k * pred_k)  %>%
      return()
  }
  mix_pred = sapply(names(pred$pred), floop) %>% rowSums()
  
  (value - mix_pred)^2 %>%
    mean() %>%
    return()
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
  t = obs %>% pull(Input)
  value = obs %>% pull(Output)
  floop = function(k)
  {
    mean = pred$pred[[k]] %>% filter(Input %in% t) %>% pull(Mean)
    sd = pred$pred[[k]] %>% filter(Input %in% t) %>% pull(Var) %>% sqrt
    
    CI_inf = mean - qnorm(1 - level/2) * sd
    CI_sup = mean + qnorm(1 - level/2) * sd
    
    100 * pred$pred[[k]]$Proba[[1]] * ((CI_inf < value) & (value < CI_sup)) %>%
      mean %>%
      return()
  }
  sapply(names(pred$pred), floop) %>% sum %>% return()
}

#### DATA IMPORT ####
raw_db_co2 = read_csv('Real_Data_Study/Data/co2-data.csv')
raw_db_weight = read_csv('Real_Data_Study/Data/db_gusto_weight.csv') 
####

#### Weight study ####
set.seed(42)
## Split by gender and draw the training/testing sets
db_w = raw_db_weight  %>% 
  dplyr::select(-sex) %>%
  split_train(ratio_train = 0.6)

## Splitting training and testing sets and select observed times
db_w_train = db_w %>%
  filter(Training == 1) %>% 
  select(- Training)
db_w_test = db_w %>%
  filter(Training == 0) %>%
  split_times(prop_test = 0.2)

mod_select = select_nb_cluster(data = db_w_train, 
                               fast_approx = F,
                               grid_nb_cluster = 2:6)
#mod = train_magma(db_w_train)
#saveRDS(mod_select,'Real_Data_Study/Training/train_weight_Hoo_mod_select.rds')
mod_select = readRDS('Real_Data_Study/Training/train_weight_Hoo_2clusters.rds')

## Model selection indicates 2 clusters as optimal choice

floop = function(i, db_test)
{ #browser()
  print(i)
  db_obs_i = db_w_test %>%
    filter(ID == i) %>%
    filter(Observed == 1) %>% 
    select(- c(Training, Observed))
  db_pred_i = db_w_test %>%
    filter(ID == i) %>%
    filter(Observed == 0) %>% 
    select(- c(Training, Observed))
  input_i_pred = db_pred_i %>% pull(Input)
  
  pred_clust = pred_magmaclust(db_obs_i, mod,
                               grid_inputs = input_i_pred, plot = F)
  
  eval_clust = tibble('ID' = i,
                      'MSE' = db_pred_i %>%  MSE_clust(pred_clust),
                      'WCIC' = db_pred_i %>% WCIC(pred_clust, 0.05))

  eval_clust %>% return()
}
db_res = db_w_test$ID %>%
  unique %>% 
  lapply(floop, db_w_test) %>%
  bind_rows

db_res %>%
  dplyr::select(-ID) %>%
  summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)

# png("weight_plot_pred4.png",res=600, height=120, width= 220, units="mm")
# plot_magmaclust(pred, cluster = 'all', data = new_db,
#               data_train = db_w_test, mean_k = list_mu) +
#   geom_point(data = new_db_t, aes(x = Input, y = Output), col = 'green')+ 
#   theme_classic() 
# dev.off()
# 
# png("weight_plot_pred_bestclust4.png",res=600, height=120, width= 220, units="mm")
# plot_magmaclust(pred, cluster = clust, data = new_db,
#               data_train = db_w_test, mean_k = list_mu) +
#   geom_point(data = new_db_t, aes(x = Input, y = Output), col = 'green')+ 
#   theme_classic() 
# dev.off()

#### CO2 STUDY ####
list_not_country = c('EU-28', 'Europe', 'Europe (excl. EU-27)', 'World',
                     'North America', 'North America (excl. USA)', 'EU-27',
                     'Europe (excl. EU-28)', 'Asia (excl. China & India)',
                     'Asia', 'South America', 'Oceania', 'Africa',
                     'Sint Maarten (Dutch part)')
db_c = raw_db_co2 %>%
  dplyr::select(country, year, co2_per_capita) %>% 
  rename(ID = country, Input = year, Output = co2_per_capita) %>%
  filter(!(ID %in% list_not_country)) %>% 
  drop_na()

####

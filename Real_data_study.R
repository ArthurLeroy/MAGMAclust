library(tidyverse)
library(png)
library(reticulate)
library(MagmaClustR)
use_condaenv( "r-reticulate")
np <- import("numpy", as = "np")
mogp <- import("mogptk", as = "mogp")



#### TRAIN/TEST SPLITING FUNCTIONS ####
split_train = function(db, ratio_train)
{
  n_indiv = db$ID%>% unique()
  n_index = (ratio_train * length(n_indiv)) %>% round()
  index_train = sample(n_indiv, n_index, replace = F)
  
  db %>% 
    mutate(Training = ifelse(ID %in% index_train, 1,0)) %>%
    return()
}

split_times = function(db, prop_test, last = TRUE)
{ 
  db = db %>% rowid_to_column("Row")
  if(last){
    row_obs = db %>%
      group_by(ID) %>%
      top_n(Input, n = - n() * (1 - prop_test)) %>%
      pull(Row)
  } else {
    row_obs = db %>%
      group_by(ID) %>%
      sample_frac(1- prop_test) %>% 
      pull(Row)
  }
  
  db %>%
    ungroup() %>%
    mutate(Observed = ifelse(db$Row %in% row_obs, 1, 0)) %>%
    dplyr::select(- Row) %>%
    return()
}


#### EVALUATION FUNCTIONS ####
MSE = function(obs, pred)
{
  input = obs %>% pull(Input)
  value = obs %>% pull(Output)
  
  mix_pred = pred %>%
    filter(Input %in% input) %>%
    pull(Mean)
  
  (value - mix_pred)^2 %>%
    mean() %>%
    return()
}

MSE_clust = function(obs, pred)
{
  input = obs %>% pull(Input)
  value = obs %>% pull(Output)
  
  mix_pred = pred$mixture_pred %>%
      filter(Input %in% input) %>%
      pull(Mean)

  (value - mix_pred)^2 %>%
    mean() %>%
    return()
}

loss = function(x, y)
{ ## return : loss function between x and y
  abs(x - y) %>% return()
}

WCIC = function(obs, pred, level)
{
  t = obs %>% pull(Input)
  value = obs %>% pull(Output)

  mean = pred %>% filter(Input %in% t) %>% pull(Mean)
  sd = pred %>% filter(Input %in% t) %>% pull(Var) %>% sqrt
    
  CI_inf = mean - qnorm(1 - level/2) * sd
  CI_sup = mean + qnorm(1 - level/2) * sd
    
  100 * ((CI_inf < value) & (value < CI_sup)) %>%
    mean %>%
    return()
}

WCIC_clust = function(obs, pred, level)
{
  # t = obs %>% pull(Input)
  # value = obs %>% pull(Output)
  # 
  # floop = function(k)
  # { 
  #   mean = pred$pred[[k]] %>% filter(Input %in% t) %>% pull(Mean)
  #   sd = pred$pred[[k]] %>% filter(Input %in% t) %>% pull(Var) %>% sqrt
  #   
  #   CI_inf = mean - qnorm(1 - level/2) * sd
  #   CI_sup = mean + qnorm(1 - level/2) * sd
  #   
  #   100 * pred$pred[[k]]$Proba[[1]] * ((CI_inf < value) & (value < CI_sup)) %>%
  #     mean %>%
  #     return()
  # }
  # sapply(names(pred$pred), floop) %>%
  #   sum %>%
  #   return()
  t = obs %>% pull(Input)
  value = obs %>% pull(Output)
  
  mean = pred$mixture_pred %>% filter(Input %in% t) %>% pull(Mean)
  sd = pred$mixture_pred %>% filter(Input %in% t) %>% pull(Var) %>% sqrt
  
  CI_inf = mean - qnorm(1 - level/2) * sd
  CI_sup = mean + qnorm(1 - level/2) * sd
  
  100 * ((CI_inf < value) & (value < CI_sup)) %>%
    mean %>%
    return()
}

## Prediction loop function over individuals
eval = function(db_test, mod_select)
{ 
  ## Extract trained model for Magma
  mod_magma = mod_select$trained_models[[1]]
  ## Extract list of names for the remaining MagmaClust trained models
  list_clust = mod_select$trained_models[-1]

  floop = function(i){
    print(i)
    db_obs_i = db_test %>%
      filter(ID == i) %>%
      filter(Observed == 1) %>% 
      select(- c(Training, Observed))
    db_pred_i = db_test %>%
      filter(ID == i) %>%
      filter(Observed == 0) %>% 
      select(- c(Training, Observed))
    input_i_pred = db_pred_i %>% pull(Input)

    ## Magma
    pred = pred_magma(db_obs_i,
                      mod_magma, 
                      grid_inputs = input_i_pred,
                      plot = F)
  
    eval = tibble('Method' = 'Magma',
                  'ID' = i,
                  'MSE' = db_pred_i %>% MSE(pred),
                  'WCIC' = db_pred_i %>% WCIC(pred, 0.05))
    ## MagmaClust
    floop_k = function(mod_magma_k){
      
      pred_clust = pred_magmaclust(db_obs_i,
                                   mod_magma_k,
                                   grid_inputs = input_i_pred,
                                   plot = F)
      
      
       tibble('Method' = paste0(
         'MagmaClust - K = ', 
          mod_magma_k$ini_args$nb_cluster),
              'ID' = i,
              'MSE' = db_pred_i %>%  MSE_clust(pred_clust),
              'WCIC' = db_pred_i %>% WCIC_clust(pred_clust, 0.05))
    }
    eval_clust = list_clust %>% 
      lapply(floop_k) %>% 
      bind_rows() %>%
      return()
  
    eval %>% 
      bind_rows(eval_clust) %>%
      return()
  } 
  db_test$ID %>%
    unique %>% 
    lapply(floop) %>%
    bind_rows %>% 
    return()
}


#### DATA IMPORT ####
raw_db_co2 = read_csv('Real_Data_Study/Data/co2-data.csv')
raw_db_weight = read_csv('Real_Data_Study/Data/db_gusto_weight.csv')
raw_db_swimming = read_csv('Real_Data_Study/Data/db_100m_freestyle.csv')

#### WEIGHT STUDY ####
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
  split_times(prop_test = 0.4)

## Model selection and training
# mod_select = select_nb_cluster(data = db_w_train, 
#                                fast_approx = F,
#                                grid_nb_cluster = 1:6, 
#                                cv_threshold = 1e-2)
#saveRDS(mod_select,'Real_Data_Study/Training/train_weight_Hoo_mod_select.rds')
mod_w_select = readRDS('Real_Data_Study/Training/train_weight_Hoo_mod_select.rds')

### MOGP
db_w_train_py = mogp$LoadDataFrame(
  db_w_train,
  x_col='Input',
  y_col='Output',
  name= unique(db_w_train$ID))
  
  py$mogp.DataSet()

## Evaluate the prediction performances
db_res_w = eval(db_w_test, mod_w_select)
# write_csv(db_res_w,'Real_Data_Study/Results/pred_weight_magma_magmaclust.csv')

## Summarise the evaluation results
db_res_w %>%
  dplyr::select(-ID) %>%
  group_by(Method) %>% 
  summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% 
  mutate(MSE_Mean = round(MSE_Mean, 1), WCIC_Mean = round(WCIC_Mean, 1),
         MSE_SD = round(MSE_SD, 1), WCIC_SD = round(WCIC_SD, 1)) %>% 
  mutate('Mean' = paste0(MSE_Mean, ' (', MSE_SD, ')'),
         'WCIC' =  paste0(WCIC_Mean, ' (', WCIC_SD, ')')) %>%
  dplyr::select(c(Method, Mean, WCIC))


# png("weight_plot_pred4.png",res=600, height=120, width= 220, units="mm")
# plot_magmaclust(pred, cluster = 'all', data = new_db,
#               data_train = db_w_test, mean_k = list_mu) +
#   geom_point(data = new_db_t, aes(x = Input, y = Output), col = 'green')+ 
#   theme_classic() 
# dev.off()
# 

#### CO2 STUDY ####
list_not_country = c('EU-28', 'Europe', 'Europe (excl. EU-27)', 'World',
                     'North America', 'North America (excl. USA)', 'EU-27',
                     'Europe (excl. EU-28)', 'Asia (excl. China & India)',
                     'Asia', 'South America', 'Oceania', 'Africa',
                     'Sint Maarten (Dutch part)')
set.seed(42)
db_c = raw_db_co2 %>%
  dplyr::select(country, year, co2_per_capita) %>% 
  rename(ID = country, Input = year, Output = co2_per_capita) %>%
  filter(!(ID %in% list_not_country)) %>% 
  drop_na() %>%
  split_train(ratio_train = 0.6)

## Splitting training and testing sets and select observed times
db_c_train = db_c %>%
  filter(Training == 1) %>% 
  select(- Training)
db_c_test = db_c %>%
  filter(Training == 0) %>%
  split_times(prop_test = 0.4)

## Model selection and training
# mod_c_select = select_nb_cluster(data = db_c_train,
#                                fast_approx = F,
#                                grid_nb_cluster = 1:6,
#                                cv_threshold = 1e-3)
# saveRDS(mod_c_select,'Real_Data_Study/Training/train_co2_Hoo_mod_select.rds')
mod_c_select = readRDS('Real_Data_Study/Training/train_co2_Hoo_mod_select.rds')

## Evaluate the prediction performances
db_res_c = eval(db_c_test, mod_c_select)
# write_csv(db_res_c, 'Real_Data_Study/Results/pred_co2_magma_magmaclust.csv')

## Summarise the evaluation results
db_res_c %>%
  dplyr::select(-ID) %>%
  group_by(Method) %>% 
  summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% 
  mutate(MSE_Mean = round(MSE_Mean, 1), WCIC_Mean = round(WCIC_Mean, 1),
         MSE_SD = round(MSE_SD, 1), WCIC_SD = round(WCIC_SD, 1)) %>% 
  mutate('Mean' = paste0(MSE_Mean, ' (', MSE_SD, ')'),
         'WCIC' =  paste0(WCIC_Mean, ' (', WCIC_SD, ')')) %>%
  dplyr::select(c(Method, Mean, WCIC)) 


#### SWIMMING STUDY ####
set.seed(42)
## Split by gender and draw the training/testing sets
db_m = raw_db_swimming %>%
  filter(Gender == 1) %>% 
  mutate(ID = as.character(ID)) %>% 
  rename(Input = Age, Output = Performance) %>% 
  dplyr::select(- Gender) %>%
  split_train(ratio_train = 0.6)

db_f = raw_db_swimming %>%
  filter(Gender == 2) %>%
  mutate(ID = as.character(ID)) %>% 
  rename(Input = Age, Output = Performance) %>% 
  dplyr::select(- Gender) %>%
  split_train(ratio_train = 0.6)

## Splitting training and testing sets and select observed times
db_m_train = db_m %>%
  filter(Training == 1) %>% 
  select(- Training)
db_m_test = db_m %>%
  filter(Training == 0) %>%
  split_times(prop_test = 0.4)

db_f_train = db_f %>%
  filter(Training == 1) %>% 
  select(- Training)
db_f_test = db_f %>%
  filter(Training == 0) %>%
  split_times(prop_test = 0.4)

## Model selection and training
# mod_m_select = select_nb_cluster(data = db_m_train,
#                                 fast_approx = F,
#                                 grid_nb_cluster = 1:6,
#                                 cv_threshold = 1e-3)
# saveRDS(mod_m_select,'Real_Data_Study/Training/train_male_Hoo_mod_select.rds')
# mod_f_select = select_nb_cluster(data = db_f_train,
#                                 grid_inputs = unique(db_f$Input),
#                                 fast_approx = F,
#                                 grid_nb_cluster = 1:6,
#                                 cv_threshold = 1e-3)
# saveRDS(mod_f_select,'Real_Data_Study/Training/train_female_Hoo_mod_select.rds')

mod_m_select = readRDS('Real_Data_Study/Training/train_male_Hoo_mod_select.rds')
mod_f_select = readRDS('Real_Data_Study/Training/train_female_Hoo_mod_select.rds')

## Evaluate the prediction performances
db_res_m = eval(db_m_test, mod_m_select)
# write_csv(db_res_m, 'Real_Data_Study/Results/pred_male_magma_magmaclust.csv')
# db_res_m = read_csv('Real_Data_Study/Results/pred_male_magma_magmaclust.csv')

## Evaluate the prediction performances
db_res_f = eval(db_f_test, mod_f_select)
 write_csv(db_res_f, 'Real_Data_Study/Results/pred_female_magma_magmaclust.csv')
#db_res_f = read_csv('Real_Data_Study/Results/pred_female_magma_magmaclust.csv')

## Summarise the evaluation results
db_res_m %>%
  dplyr::select(-ID) %>%
  group_by(Method) %>% 
  summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% 
  mutate(MSE_Mean = round(MSE_Mean, 1), WCIC_Mean = round(WCIC_Mean, 1),
         MSE_SD = round(MSE_SD, 1), WCIC_SD = round(WCIC_SD, 1)) %>% 
  mutate('Mean' = paste0(MSE_Mean, ' (', MSE_SD, ')'),
         'WCIC' =  paste0(WCIC_Mean, ' (', WCIC_SD, ')')) %>%
  dplyr::select(c(Method, Mean, WCIC)) 

## Summarise the evaluation results
db_res_f %>%
  dplyr::select(-ID) %>%
  group_by(Method) %>% 
  summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% 
  mutate(MSE_Mean = round(MSE_Mean, 1), WCIC_Mean = round(WCIC_Mean, 1),
         MSE_SD = round(MSE_SD, 1), WCIC_SD = round(WCIC_SD, 1)) %>% 
  mutate('Mean' = paste0(MSE_Mean, ' (', MSE_SD, ')'),
         'WCIC' =  paste0(WCIC_Mean, ' (', WCIC_SD, ')')) %>%
  dplyr::select(c(Method, Mean, WCIC)) 

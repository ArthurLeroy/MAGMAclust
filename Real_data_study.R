library(tidyverse)
library(png)
library(MagmaClustR)

#### TRAIN/TEST SPLITING FUNCTIONS ####
split_train = function(db, ratio_train)
{
  n_indiv = db$ID %>% unique()
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

## Evaluate the prediction performances
# db_res_w = eval(db_w_test, mod_w_select)
# write_csv(db_res_w,'Real_Data_Study/Results/pred_weight_magma_magmaclust.csv')
# db_res_w = read_csv('Real_Data_Study/Results/pred_weight_magma_magmaclust.csv')


## Summarise the evaluation results
# db_res_w %>%
#   dplyr::select(-ID) %>%
#   group_by(Method) %>% 
#   summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% 
#   mutate(MSE_Mean = round(MSE_Mean, 1), WCIC_Mean = round(WCIC_Mean, 1),
#          MSE_SD = round(MSE_SD, 1), WCIC_SD = round(WCIC_SD, 1)) %>% 
#   mutate('Mean' = paste0(MSE_Mean, ' (', MSE_SD, ')'),
#          'WCIC' =  paste0(WCIC_Mean, ' (', WCIC_SD, ')')) %>%
#   dplyr::select(c(Method, Mean, WCIC))
# 
# ## Plot an example
# id_w = '010-21288'
# 
# db_obs_w = db_w_test %>%
#   filter(ID == id_w) %>%
#   filter(Observed == 1) %>%
#   select(- c(Training, Observed))
# db_pred_w = db_w_test %>%
#   filter(ID == id_w) %>%
#   filter(Observed == 0) %>%
#   select(- c(Training, Observed))
# 
# pred_ex_w = pred_magmaclust(
#   db_obs_w,
#   mod_w_select$trained_models[[3]],
#   grid_inputs = seq(0, 72, 0.1),
#   hyperpost = hyp_w,
#   get_hyperpost = FALSE)

# hyp_w = pred_ex_w$hyperpost
# col_db_w = data_allocate_cluster(mod_w_select$trained_models[[3]])

# png("pred_example_weigt_data.png",res=600, height=120, width= 220, units="mm")
# plot_magmaclust(pred_ex_w, cluster = 'all',
#                 data_train = col_db_w, col_clust = TRUE, size_data = 4,
#                 heatmap = TRUE, y_grid = seq(0, 40, 0.2),
#                 prior_mean = hyp_w$mean, size_data_train = 0.5) +
#   geom_point(data = db_pred_w, aes(x = Input, y = Output),
#              col = 'yellow', size = 2) +
#   geom_point(data = db_obs_w, aes(x = Input, y = Output),
#              col = 'black', size = 2) +
#   theme_classic() + ggtitle("") +
#   xlab('Age (in months)') +
#   ylab('Weight (in kg)')
# dev.off()


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
  split_times(prop_test = 0.4, last = F)

## Model selection and training
# mod_c_select = select_nb_cluster(data = db_c_train,
#                                fast_approx = F,
#                                grid_nb_cluster = 1:6,
#                                cv_threshold = 1e-3)
# saveRDS(mod_c_select,'Real_Data_Study/Training/train_co2_Hoo_mod_select.rds')
mod_c_select = readRDS('Real_Data_Study/Training/train_co2_Hoo_mod_select.rds')

## Evaluate the prediction performances
# db_res_c = eval(db_c_test, mod_c_select)
# write_csv(db_res_c, 'Real_Data_Study/Results/pred_co2_magma_magmaclust.csv')
db_res_c = read_csv('Real_Data_Study/Results/pred_co2_magma_magmaclust.csv')

## Summarise the evaluation results
db_res_c %>%
  filter(!(ID %in% c('Brunei', 'Curacao'))) %>%
  dplyr::select(-ID)  %>% 
  group_by(Method) %>% 
  summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% 
  mutate(MSE_Mean = round(MSE_Mean, 1), WCIC_Mean = round(WCIC_Mean, 1),
         MSE_SD = round(MSE_SD, 1), WCIC_SD = round(WCIC_SD, 1)) %>% 
  mutate('Mean' = paste0(MSE_Mean, ' (', MSE_SD, ')'),
         'WCIC' =  paste0(WCIC_Mean, ' (', WCIC_SD, ')')) %>%
  dplyr::select(c(Method, Mean, WCIC)) 

## Plot an example
png("pred_example_co2_data.png",res=600, height=120, width= 220, units="mm")
ggplot(db_c %>% slice(1:1600)) +
  geom_point(aes(x = Input, y = Output, col = ID), size = 0.6) +
  theme_classic() + xlab('Year') + ylab('CO2 emission per capita')
dev.off()

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

## Plot an example
id_f = '296'

db_obs_f = db_f_test %>%
  filter(ID == id_f) %>%
  filter(Observed == 1) %>%
  select(- c(Training, Observed))
db_pred_f = db_f_test %>%
  filter(ID == id_f) %>%
  filter(Observed == 0) %>%
  select(- c(Training, Observed))

pred_ex_f = pred_magmaclust(
  db_obs_f,
  mod_f_select$trained_models[[5]],
  grid_inputs = seq(10, 20, 0.1),
  hyperpost = hyp_f,
  get_hyperpost = FALSE)

#hyp_f = pred_ex_f$hyperpost
col_db_f = data_allocate_cluster(mod_f_select$trained_models[[5]]) %>% 
  slice(1:1000)

png("pred_example_women_swimming_data.png",res=600, height=120, width= 220, units="mm")
plot_magmaclust(pred_ex_f, cluster = 'all', data = db_obs_f,
                data_train = col_db_f, col_clust = TRUE, size_data = 4,
                heatmap = TRUE, y_grid = seq(50, 110, 0.1),
                prior_mean = hyp_f$mean, size_data_train = 0.5) +
  geom_point(data = db_pred_f, aes(x = Input, y = Output),
             col = 'yellow', size = 2) +
  theme_classic() + ggtitle("")
dev.off()

## Plot an example
id_m = "1164" #'778' 

db_obs_m = db_m_test %>%
  filter(ID == id_m) %>%
  filter(Observed == 1) %>%
  select(- c(Training, Observed))
db_pred_m = db_m_test %>%
  filter(ID == id_m) %>%
  filter(Observed == 0) %>%
  select(- c(Training, Observed))

pred_ex_m = pred_magmaclust(
  db_obs_m,
  mod_m_select$trained_models[[5]],
  grid_inputs = seq(10, 20, 0.01),
  hyperpost = hyp_m,
  get_hyperpost = FALSE, plot = F)

#hyp_m = pred_ex_m$hyperpost
col_db_m = data_allocate_cluster(mod_m_select$trained_models[[5]]) %>% 
  slice(1:1000)

png("pred_example_men_swimming_data.png",res=600, height=120, width= 220, units="mm")
plot_magmaclust(pred_ex_m, cluster = 'all', data = db_obs_m,
                data_train = col_db_m, col_clust = TRUE, size_data = 4,
                heatmap = TRUE, y_grid = seq(45, 110, 0.2),
                prior_mean = hyp_m$mean, size_data_train = 0.5) +
  geom_point(data = db_pred_m, aes(x = Input, y = Output),
             col = 'yellow', size = 2) +
  theme_classic() + ggtitle("")
dev.off()

#### Illustration example ####
set.seed(42)
db_illu = simu_db(M = 5, N = 30, K = 3, common_input = FALSE,
                  int_i_sigma = c(3,3), int_mu_v = c(5,5))
db_train_illu = db_illu %>% filter(ID != 'ID5-Clust3')
db_pred_illu = db_illu %>% filter(ID == 'ID5-Clust3') %>% head(6)
db_test_illu = db_illu %>% filter(ID == 'ID5-Clust3') %>% tail(24)

#MagmaClustR:::plot_db(db_illu)

mod_illu = train_magmaclust(db_train_illu)

mod_illu$ini_args$common_hp_i = FALSE
pred_illu = pred_magmaclust(db_pred_illu, mod_illu, get_hyperpost = TRUE,
                            grid_inputs = seq(0,10, 0.1))
  
db_illu_clust = data_allocate_cluster(mod_illu)

plot_magmaclust(pred_illu, data = db_pred_illu, data_train = db_illu_clust, 
                col_clust = TRUE, size_data = 4, y_grid = seq(5, 40, 0.2),
                size_data_train = 0.5, heatmap = TRUE,
                prior_mean = pred_illu$hyperpost$mean) + 
  geom_point(data = db_test_illu, aes(x = Input, y = Output),
             col = 'orange', size = 2) + theme_classic() + ggtitle("")



#### Motivation example ####
db_motiv = db_m %>% 
  filter(ID %in% unique(db_m$ID)[c(5, 14, 87, 73, 1011)])

png("motivation_example.png",res=600, height=90, width= 270, units="mm")
ggplot(db_motiv) +
  geom_point(aes(x = Input, y = Output, col = ID)) +
  theme_classic() + guides(col = 'none') + xlab('Age') + 
  ylab('Performance (in sec)') + scale_x_continuous(breaks = 10:20) + 
  scale_y_continuous(breaks = seq(50, 85, 5))
dev.off()

#### New graphs for presentations ####

## Kernel illustration
library(gganimate)

input = seq(0, 10, 0.05)

kern = kern_to_cov(input, kern = 'SE', hp = tibble(
  se_variance = 1,
  se_lengthscale = 1))
# kern = kern_to_cov(input, kern = 'PERIO', hp = tibble(
#   perio_variance = 1,
#   perio_lengthscale = 1,
#   period= 0.5) )
# kern = kern_to_cov(input, kern = 'LIN', hp = tibble(
#   lin_slope = 1,
#   lin_offset = 0) )
  
draw_samples = function(input = seq(0, 10, 0.05),
                        mean = NULL, 
                        nb_samples = 5,
                        kern){
  
  size = length(input)
  
  floop = function(i){
    if(mean %>% is.null()){
      mean = rep(0, size)
    }
    
      tibble(
    "ID" = rep(1:nb_samples, size) %>%  as.factor(),
    "Input" = rep(input, each = nb_samples),
    "Output" = mvtnorm::rmvnorm(nb_samples, mean, kern)%>%as.vector,
    "Index" = i
    ) %>% return()
  }
  1:nb_samples %>% 
    lapply(floop) %>% 
    bind_rows() %>% 
    return()
}

samples = draw_samples(nb_samples = 5, kern = kern)

gg_anim = ggplot(samples) + geom_line(aes(x = Input, y = Output, col = ID)) + 
  guides(col = 'none') +
  theme_classic() + transition_states(Index)

animate(gg_anim, height = 1600, width = 2000, res = 300)

anim_save("illu_prior_gp.gif")

## Trained GP illustration
set.seed(7)
db = simu_db(M = 1, N = 5)

hp_gp = train_gp(db)
samples = draw_samples(kern = kern_to_cov(input,
                                          kern = 'SE',
                                          hp = hp_gp %>% select(- noise)))

gg_anim = ggplot(samples) + geom_line(aes(x = Input, y = Output, col = ID)) +
  geom_point(data = db, aes(Input, Output)) + 
  guides(col = 'none') + theme_classic() + transition_states(Index)

animate(gg_anim, height = 1600, width = 2000, res = 300)

anim_save("illu_trained_gp.gif")

## Posterior samples GP illustration
sub_db = db[c(1, 2, 4),]
pred_gp = pred_gp(sub_db, grid_inputs = input, hp = hp_gp, get_full_cov = T) 

samples = draw_samples(mean = pred_gp$pred$Mean, kern = pred_gp$cov)

gg_anim = ggplot(samples) + geom_line(aes(x = Input, y = Output, col = ID)) +
  geom_point(data = sub_db, aes(Input, Output)) + 
  guides(col = 'none') + theme_classic() + transition_states(Index)

animate(gg_anim, height = 1600, width = 2000, res = 300)

anim_save("illu_post_gp4.gif")

## Posterior GP illustration
pred_gp = pred_gp(sub_db, grid_inputs = input, hp = hp_gp) 

plot_gp = plot_gp(pred_gp, data = sub_db)

ggsave("illu_post_gp3.png", plot_gp, height=1600, width=2000, dpi=300, units="px")

## Illustration 2-D GIF

db = simu_db(covariate = T)
mod = train_magma(db)
pred_gif = pred_gif(db %>% filter(ID == 1), mod)

plot_gif = plot_gif(pred_gif, data = (db %>% filter(ID == 1)) ) 
animate(plot_gif, height = 1200, width = 2000, res = 300)

anim_save("illu_magma_2D.gif")

## Illustration GP and Magma GIF
set.seed(10)
db = simu_db(common_input = F)
#MagmaClustR:::plot_db(db)

db_i = db %>% filter(ID == 6) %>% filter(Input < 4.5)
test_db_i = db %>% filter(ID == 6) %>% filter(Input > 4.5)


pred_gp = pred_gp(db_i, grid_inputs = seq(0, 10, 0.01))
#mod_magma = train_magma(db)
pred_magma = pred_magma(db_i, mod_magma, grid_inputs = seq(0, 10, 0.01), 
                        get_hyperpost = T)

illu_gp = plot_gp(pred_gp, data = db_i) +
  geom_point(data = test_db_i, aes(Input, Output), col = 'red') + 
  ylim(c(-37, 27))

illu_magma = plot_gp(pred_magma, data = db_i, data_train = db_i, 
        prior_mean = pred_magma$hyperpost$mean) +
  geom_point(data = test_db_i, aes(Input, Output), col = 'red')

ggsave("illu_gp.png", illu_gp, height=1200, width=2000, dpi=300, units="px")
ggsave("illu_magma.png",illu_magma, height=1200, width=2000, dpi=300, units="px")


pred_gp_gif = pred_gif(db %>% filter(ID == 6), grid_inputs = seq(0, 10, 0.01))
pred_magma_gif = pred_gif(db %>% filter(ID == 6), mod_magma,
                          grid_inputs = seq(0, 10, 0.01))

## Use to put add a layer to ggplot backward
# `-.gg` <- function(plot, layer) {
#   if (missing(layer)) {
#     stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
#   }
#   if (!is.ggplot(plot)) {
#     stop('Need a plot on the left side')
#   }
#   plot$layers = c(layer, plot$layers)
#   plot
# }

illu_gp_gif = plot_gif(
  pred_gp_gif %>%filter(Index > 4) %>% mutate(Index = Index - 4),
  data = test_db_i) -
  geom_point(data = test_db_i, aes(Input, Output), col = 'red', size = 2) +
  geom_point(data = db_i, aes(Input, Output), col = 'black', size = 2) 

illu_magma_gif = plot_gif(
  pred_magma_gif %>% filter(Index > 4) %>% mutate(Index = Index - 4),
  data = test_db_i,
  data_train = db_i,
  prior_mean = pred_magma$hyperpost$mean) -
  geom_point(data = test_db_i, aes(Input, Output), col = 'red', size = 2) +
  geom_point(data = db_i, aes(Input, Output), col = 'black', size = 2) 

animate(illu_gp_gif, height = 1200, width = 2000, res = 300)
anim_save("illu_gp_gif.gif")
animate(illu_magma_gif, height = 1200, width = 2000, res = 300)
anim_save("illu_magma_gif.gif")

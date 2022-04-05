library(tidyverse)
library(MagmaClustR)
library(compositions)
library(plotly)

data(juraset)
db_co2 = read_csv('Data/co2-data.csv')
db_weight = read_csv('Data/db_gusto_weight.csv')
db_jura = juraset %>%
  mutate(ID = paste(Rock, Land)) %>% 
  select(ID, X, Y, Zn)


ggplot(db_jura) + geom_tile(aes(x = X, y = Y, fill = Zn))


plot_ly(x=db_jura$X, y=db_jura$Y, z=db_jura$Zn, 
        type="scatter3d", mode="markers", color=db_jura$ID)

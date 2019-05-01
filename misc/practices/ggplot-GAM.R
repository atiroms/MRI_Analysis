# It's possible you don't need all of these
library(tidyverse)
library(magrittr)
library(ggplot2)
library(nlme)
library(MuMIn)
library(mgcv)

# Reading the data from file and making the formats easier for me to work with
land_clim_sa_climate <- read_csv("/users/aperium/Dropbox/Projects/MSES/CLasses/Bios 8700 - Biostats II/WorkingDirectory/LandClim_SA_climate.csv")
land_clim_sa_factor <- read_csv("/users/aperium/Dropbox/Projects/MSES/CLasses/Bios 8700 - Biostats II/WorkingDirectory/LandClim_SA_factor.csv")
sa_climate <- dplyr::select(land_clim_sa_climate, propOpen, climate, gp, diff_Temp2, diff_PPT2) %>% dplyr::rename(open_pasture = propOpen, climate_scenario = climate, grazing_pressure = gp, delta_temp = diff_Temp2, delta_precip = diff_PPT2) %>% mutate(climate_scenario = factor(climate_scenario, levels = unique(climate_scenario)))
sa_factor <- dplyr::select(land_clim_sa_factor, propOpen, climate, gp) %>% dplyr::rename(open_pasture = propOpen, climate_scenario = climate, grazing_pressure = gp) %>% mutate(climate_scenario = factor(climate_scenario, levels = unique(climate_scenario)))

#scaling data so that partial regression coefficients can be used as effect size measures
sa_climate_std <- dplyr::mutate_if(sa_climate, is.numeric, scale)

# A GAM model
m3 <- gam(open_pasture ~ s(grazing_pressure) + climate_scenario, data = sa_factor)

# pulling the prediction and residual data from the model
sa_factor %<>% mutate(resid = resid(m3), predict = predict(m3))

# plotting the data.
plot <- ggplot(sa_factor, aes(x = grazing_pressure, y = open_pasture, color = climate_scenario)) + 
  scale_color_grey(name = "climate scenario") +  
  geom_point(shape = 19, stroke = 1) +
  xlab("grazing pressure") + ylab("openness") +
  theme_bw()

# plotting the predictions lines for each level of the factor
plot %<>% + geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = filter(sa_factor, climate_scenario == unique(sa_factor$climate_scenario)[1]), size = 1, se = FALSE)
plot %<>% + geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = filter(sa_factor, climate_scenario == unique(sa_factor$climate_scenario)[2]), size = 1, se = FALSE)
plot %<>% + geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = filter(sa_factor, climate_scenario == unique(sa_factor$climate_scenario)[3]), size = 1, se = FALSE)
plot %<>% + geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = filter(sa_factor, climate_scenario == unique(sa_factor$climate_scenario)[4]), size = 1, se = FALSE)
plot %<>% + geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = filter(sa_factor, climate_scenario == unique(sa_factor$climate_scenario)[5]), size = 1, se = FALSE)
plot %<>% + geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = filter(sa_factor, climate_scenario == unique(sa_factor$climate_scenario)[6]), size = 1, se = FALSE)
plot %<>% + geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = filter(sa_factor, climate_scenario == unique(sa_factor$climate_scenario)[7]), size = 1, se = FALSE)
plot %<>% + geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = filter(sa_factor, climate_scenario == unique(sa_factor$climate_scenario)[8]), size = 1, se = FALSE)
plot %<>% + geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = filter(sa_factor, climate_scenario == unique(sa_factor$climate_scenario)[9]), size = 1, se = FALSE)
plot %<>% + geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = filter(sa_factor, climate_scenario == unique(sa_factor$climate_scenario)[10]), size = 1, se = FALSE)
plot %<>% + geom_smooth(aes(y = predict),  method = gam, formula = y ~ s(x, bs = "ps"), data = filter(sa_factor, climate_scenario == unique(sa_factor$climate_scenario)[11]), size = 1, se = FALSE)

# printing the final plot.
plot
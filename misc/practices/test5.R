# https://fromthebottomoftheheap.net/2017/12/14/difference-splines-ii/

library('readr')
library('dplyr')
library('ggplot2')
theme_set(theme_bw())
library('mgcv')

####

uri <- 'https://gist.githubusercontent.com/gavinsimpson/eb4ff24fa9924a588e6ee60dfae8746f/raw/geochimica-metals.csv'

metals <- read_csv(uri, skip = 1, col_types = c('ciccd'))
metals <- mutate(metals, SiteCode = factor(SiteCode))
metals <- mutate(metals,oSoilType = ordered(SoilType, levels = c('non-eroded','eroded','thin')))
metals <- mutate(metals,fSoilType = factor(SoilType, levels = c('non-eroded','eroded','thin'),ordered=F))

class(metals$SoilType)
class(metals$fSoilType)
class(metals$oSoilType)

####

m1 <- gam(Hg ~ oSoilType + s(Date) + s(Date, by = oSoilType), data = metals,method = 'REML')
summary(m1)
plot(m1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

####

m2 <- gam(Hg ~ SoilType + s(Date) + s(Date, by = oSoilType), data = metals,method = 'REML')
summary(m2)
plot(m2, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

####

m3 <- gam(Hg ~ fSoilType + s(Date) + s(Date, by = fSoilType), data = metals,method = 'REML')
summary(m3)
plot(m3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

####

m4<- gam(Hg ~ fSoilType +           s(Date, by = fSoilType), data = metals,method = 'REML')
summary(m4)
plot(m4, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

####

m5<- gam(Hg ~                      s(Date)+s(fSoilType,bs='re'), data = metals,method = 'REML')
summary(m5)
plot(m5, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

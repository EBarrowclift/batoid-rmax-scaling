# Load required packages
library(ape) 
library(caper)
library(car)
library(nlme)
library(tidyverse)
library(ggplot2)
library(stringr)

# load data------------------------------------------------------
dat <- read.csv("data/batoid_model_data.csv", stringsAsFactors = FALSE,header = TRUE) # 83 species
phy <- read.tree("data/stein-et-al-single.tree")

# Fit pgls models------------------------------------------------

# combine phylogeny with dataset for use in pgls function
cd <- comparative.data(phy, dat, names.col = "BinomName") 

# Fit 18 pgls models 

# intercept only
zslope_0 <- pgls(log(rmax) ~ 1, data = cd, lambda = "ML") 

# rmax varies with mass only 
zslope_wt <- pgls(log(rmax) ~ log_wt, data = cd, lambda = "ML") 

# depth only
zslope_depth <- pgls(log(rmax) ~ depth_scaled, data = cd, lambda = "ML")

# temp only
zslope_invtemp <- pgls(log(rmax) ~ invtemp_scaled, data = cd, lambda = "ML") 

# mass + depth
zslope_wt_depth <- pgls(log(rmax) ~ log_wt + depth_scaled,
                                  data = cd, lambda = "ML")

# mass + temp
zslope_wt_invtemp <- pgls(log(rmax) ~ log_wt + invtemp_scaled,
                                    data = cd, lambda = "ML")

# mass * depth 
zslope_wt_x_depth <- pgls(log(rmax) ~ log_wt * depth_scaled, 
                                       data = cd, lambda = "ML")

# mass * temp 
zslope_wt_x_invtemp <- pgls(log(rmax) ~ log_wt * invtemp_scaled, 
                                         data = cd, lambda = "ML")

# varying intercepts - rays (coded 0) and skates (coded 1)
zslope_2 <- pgls(log(rmax) ~ 1 + Order, data = cd, lambda = "ML") 

# mass + order
zslope_wt_o <- pgls(log(rmax) ~ log_wt + Order, data = cd, lambda = "ML") 

# depth + order
zslope_depth_o <- pgls(log(rmax) ~ depth_scaled + Order, data = cd, lambda = "ML")

# temp + order
zslope_temp_o <- pgls(log(rmax) ~ invtemp_scaled + Order, data = cd, lambda = "ML")

# mass + depth + order
zslope_wt_depth_o <- pgls(log(rmax) ~ log_wt + depth_scaled + Order,
                                     data = cd, lambda = "ML")

# mass + temp + order
zslope_wt_invtemp_o <- pgls(log(rmax) ~ log_wt + invtemp_scaled + Order,
                                     data = cd, lambda = "ML")

# mass * depth + order
zslope_wt_x_depth_o <- pgls(log(rmax) ~ log_wt * depth_scaled + Order, 
                                       data = cd, lambda = "ML")

# mass * temp + order
zslope_wt_x_invtemp_o <- pgls(log(rmax) ~ log_wt  * invtemp_scaled + Order,
                                         data = cd, lambda = "ML")

# mass + temp + depth
zslope_wt_invtemp_depth <- pgls(log(rmax) ~ log_wt + invtemp_scaled + depth_scaled,
                                  data = cd, lambda = "ML")

# mass + temp * depth
zslope_wt_invtemp_x_depth <- pgls(log(rmax) ~ log_wt  + invtemp_scaled * depth_scaled,
                                  data = cd, lambda = "ML")

#  saving all 18 models in list
capermodels <- list(zslope_0,
                    zslope_wt, 
                    zslope_depth, 
                    zslope_invtemp, 
                    zslope_wt_depth,
                    zslope_wt_invtemp, 
                    zslope_wt_x_depth,
                    zslope_wt_x_invtemp, 
                    zslope_2, 
                    zslope_wt_o, 
                    zslope_depth_o, 
                    zslope_temp_o, 
                    zslope_wt_depth_o,
                    zslope_wt_invtemp_o,
                    zslope_wt_x_depth_o,
                    zslope_wt_x_invtemp_o,
                    zslope_wt_invtemp_depth, 
                    zslope_wt_invtemp_x_depth
) 

# manually extracting information from each model object
aics <- sapply(capermodels, function (x) bbmle::AIC(x)) 
aiccs <- sapply(capermodels, function (x) x$aicc) 
formulas <- sapply(capermodels,   # tidying formulas for easier reading
                   function (x) deparse(formula(x), width.cutoff = 90L)) %>%
  sub("^log\\(\\w+\\)\\s\\~\\s", "", .) %>% 
  gsub("_scaled", "", .) %>%
  gsub("\\_wt", "(M)", .)  %>%
  gsub("O", "o", .) # so order is a lower case
r2s <- sapply(capermodels, function (x) summary(x)$r.squared) 
ar2s <- sapply(capermodels, function (x) summary(x)$adj.r.squared)
LL <- sapply(capermodels, function (x) x$model$log.lik) 
ks <- sapply(capermodels, function (x) x$k) 

models_table <- data.frame(formulas, ks, LL, aics, aiccs, r2s, ar2s) %>%
  rename(Model = formulas, n = ks, AIC = aics,
         AICc = aiccs, R_sq = r2s, adj_R_sq = ar2s) %>%
  mutate(Model = as.character(formulas), 
         LL = round(LL, 1), 
         AIC = round(AIC, 1), AICc = round(AICc, 1), 
         dAIC = round(AIC - min(AIC), 2),
         dAICc = round(AICc - min(AICc), 2),
         R_sq = round(R_sq, 2), adj_R_sq = round(adj_R_sq, 2), 
         Weights = round(exp(-dAICc/2)/sum(exp(-dAICc/2)), 3)) 
         
models_table

### Top-ranked model diagnostic plots and analyses---------------------------------
topmod <- which(models_table$dAICc == 0)
bestmodel <- capermodels[[topmod]]

# Re-running the model using nlme:gls for diagnostics
# parameter estimates are obtained
# corPagel() is used as it's equivalent to `lambda = "ML"` in pgls
bestmodel_gls <- nlme::gls(formula(bestmodel), data = cd$data, 
                           correlation = corPagel(1, cd$phy)) 

# almost identical coefficent estimates
rbind(coef(bestmodel), coef(bestmodel_gls))

# residual plots almost identical as well, homoscedastic
plot(resid(bestmodel) ~ fitted(bestmodel)) # pgls
plot(bestmodel_gls) # gls

# estimate Variance inflation factors (on gls model) to test for collinearity
vif(bestmodel_gls)

# qqplot (on gls model) looks normal
qqnorm(bestmodel_gls)

# profile plot of lambda
plot(pgls.profile(bestmodel)) # lambda closer to 1 indicates strong phylogenetic signal

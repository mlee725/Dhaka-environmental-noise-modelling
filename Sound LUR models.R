##############################
#
## R script for Dhaka sound LUR model
#
##############################

# This script runs forward predictor selection for sound (LAeq24, Lday, Lnight) models (mixed effect)

##################
## Load in packages =====
#################
library(lme4)
library(sf)
library(dplyr)
library(MASS)
library(lubridate)
library(hms)
library(MuMIn)
library(performance)
library(stringr)
library(splines)

##################
## Load in data =====
#################
source("~/Source script for pm and noise for Dhaka.R")
Ldata <- read.csv("~/LUR spatial data.csv")

##################
## Functions =====
#################
## function to get R2 for linear, splines, and log - this function was created with the help of Google Gemini
get_r2_models_wide <- function(y, X) {
  if (!is.vector(y)) stop("y must be a vector")
  X <- as.data.frame(X)
  # Ensure all predictors are numeric 
  X <- X[sapply(X, is.numeric)]
  # Function to safely compute R²
  safe_r2 <- function(formula) {
    model <- try(lm(formula), silent = TRUE)
    if (inherits(model, "try-error")) return(NA)
    summary(model)$r.squared
  }
  # Initialize empty list for long-form results
  r2_list <- list()
  for (var in names(X)) {
    x <- X[[var]]
    df_temp <- data.frame(
      variable = var,
      model = c("linear", "log", "spline"),
      r2 = c(
        safe_r2(y ~ x),
        safe_r2(y ~ log(x+0.001)),
        safe_r2(y ~ splines::ns(x, df = 3))
      )
    )
    r2_list[[var]] <- df_temp
  }
  # Combine and pivot
  r2_long <- bind_rows(r2_list)
  r2_long[(!r2_long$variable %in% spl_vars) & r2_long$model != "linear", "r2"] <- NA # remove values for log and spline values for variables if they should not have them
  r2_wide <- r2_long %>%
    tidyr::drop_na() %>%
    mutate(col_name = paste(variable, model, sep = ".")) %>%
    dplyr::select(col_name, r2) %>%
    tidyr::pivot_wider(names_from = col_name, values_from = r2)
  return(r2_wide)
}

##################
## LUR modelling =====
#################

##################
## Set-up

## List of types of models
List <- c("noise__All", "noise__Day", "noise__Night")

## Vars that are candidates for spline transformation
spl_vars <- c("green_100", "green_200", "green_500", "green_1000",
              "com_ind_100", "com_ind_200", "com_ind_500", "com_ind_1000","TertOM_50", "TertOM_100",
              "MajOM_50", "MajOM_100", "MajOM_200", "MajOM_500", "MajOM_1000",
              "ndvi_50", "ndvi_100", "ndvi_200", "ndvi_500", "ndvi_1000",
              "popden_50", "popden_100", "popden_200", "popden_500", "popden_1000",
              "pow_200", "pow_500", "pow_1000")

## Data with sampling site lat long
Sit <- data.frame(Site = Ssites$Site,
                  X = st_coordinates(Ssites)[,1],
                  Y = st_coordinates(Ssites)[,2])

## emtpy data frames and lists
FullPred <- data.frame()
mod_pred <- data.frame()
mod_stat <- data.frame()
y_predicted_cv <- data.frame()
x20_Cv_dataset_ols <- data.frame()
mod_coef <- data.frame()
models <- list()

##################
## For-loop for predictor selection
for (i in c(1:length(List))) {
  type <- sub("^([^_]*_[^_]*).", "\\2", List[i])
  data <- noise2
  Ldata2 <- Ldata
  if(type == "Day"|type == "Night"){
    data <- noise2[noise2$TimeofDay == type,]
  }
  
  noisehour <- aggregate(data = data, LEQ ~ Site + List + Type + hour + RH + temp, median)
  leq_ave <- aggregate(data = noisehour, LEQ ~ Site + List + Type + RH + temp, median)
  
  fulldata <- merge(Ldata2, leq_ave, by.all = "Site")
  
  fulldata <- fulldata[fulldata$Site != "BARC-F1",] # remove rooftop site
  
  cols <- fulldata[,c("steel_200", "steel_500", "steel_1000", "airport_1000", "airport_3000")] # the binary varaibles
  Site <- fulldata$Site
  
  # -- Step 1 - remove predictors with low variability in values -- #
  fulldata <- fulldata %>% 
    select_if(is.numeric) %>%
    select_if(~0 != quantile(., 0.8, na.rm = T)[[1]]) %>% # remove columns 0 and 80 # remove colums with 10 and 80 quantiles have the same values
    mutate(Site = Site) %>%
    cbind(cols)
  
  # -- Step 2 - rank predictors R2 with LEQ -- #
  Correlation_highest <- get_r2_models_wide(fulldata$LEQ, dplyr::select(fulldata, -c(LEQ, Site))) # remove outcome and random intercept vars
  Correlation_highest <- dplyr::select(Correlation_highest, -dplyr::contains(".spline")) # remove splines since they did not improve model preformance 
  
  ## pull the predictor with the highest assoication with outcome and order from highest to lowest
  Correlation_highest <- Correlation_highest %>%
    reshape2::melt() %>%
    na.omit() %>%
    mutate(var = as.character(variable)) %>%
    mutate(buf = sub(".*\\_+", "", var)) %>%
    mutate(cat = sub("\\_.*", "", var)) %>%
    mutate(cat = sub("\\..*", "", cat)) %>%
    mutate(var2 = case_when(str_detect(var, "log") ~ paste0("log(", sub("\\..*", "", var), " + 0.01)"),
                            str_detect(var, "spline") ~ paste0("ns(", sub("\\..*", "", var), ", df = 3)"),
                            str_detect(var, "linear") ~ paste0(sub("\\..*", "", var)))) %>%
    dplyr::select(c("value","var","buf","cat", "var2" )) %>%
    group_by(cat) %>%
    slice(which.max(abs(value))) %>% # pull out top in group
    dplyr::arrange(desc(abs(value))) # order from highest to lowest

  list_var <- Correlation_highest$var # create list for selection for loop
  
  selected_vars <- c() # list for the selected variables
  
  # -- Step 3 - forward selection -- #
  for (j in 1:length(list_var)) {
    if(str_detect(list_var[j], "log")){
      nvar <- paste0("log(", sub("\\..*", "", list_var[j]), "+0.01)")
    }else if(str_detect(list_var[j], "spline")){
      nvar <- paste0("ns(", sub("\\..*", "", list_var[j]), ", df = 3)")
    }else{
      nvar <- paste0(sub("\\..*", "", list_var[j]))
    }
    
    if(length(selected_vars) == 0){ # if no variables selected yet
      formula <- as.formula(noquote(paste("LEQ ~ ", nvar, " + (1|Site)")))
      old_model <- lmer(LEQ ~ (1|Site), data = fulldata)
    }else{
      formula <- as.formula(noquote(paste("LEQ ~ ", nvar, " + ", paste0(selected_vars,  collapse = " + "), "+ (1|Site)")))
    }
    temp_model <- lmer(formula, data=fulldata)
    summary(temp_model)
    
    temp_adj <- summary(lm(fulldata$LEQ ~ predict(temp_model, fulldata, re.form=NA)))$r.squared
    old_adj <- summary(lm(fulldata$LEQ ~ predict(old_model, fulldata, re.form=NA)))$r.squared
    perchg <- ((temp_adj-old_adj)/old_adj)*100
    perchg
    if (perchg >= 2) {
      selected_vars <- c(selected_vars, nvar)
      old_model <- temp_model
    }
  }
  
  # -- model performance -- #
  formula <- as.formula(noquote(paste("LEQ ~ ", paste0(selected_vars,  collapse = " + "), "+ (1|Site)")))
  mod <- lmer(formula, data = fulldata)
  
  output <- data.frame(sites = fulldata$Site,
                       obs = fulldata$LEQ,
                       pred = predict(mod, fulldata, re.form=NA))
  output$res <- output$obs-output$pred
  
  output_agg <- output %>%
    group_by(sites) %>% 
    summarise(across(c(obs, pred, res), .f = list(mean = mean), na.rm = TRUE))
  
  r2_lm <- summary(lm(fulldata$LEQ ~ predict(mod, fulldata, re.form=NA)))$r.squared
  r2_lm_agg <- summary(lm(output_agg$obs_mean ~ output_agg$pred_mean))$r.squared
  
  models[[i]] <- mod
  
  ## output the coefs with stats
  end <- length(fixef(mod))+2
  coefs <- data.frame(var = names(fixef(mod)),
                      vals = fixef(mod),
                      ci_l = confint(mod)[3:end,1],
                      ci_u = confint(mod)[3:end,2],
                      vif = c(0, car::vif(mod)),
                      r2_marg = MuMIn::r.squaredGLMM(mod)[1],
                      r2_lm = r2_lm,
                      # r2_lm_agg = r2_lm_agg,
                      rmse = performance(mod)[[7]],
                      me_out = mean(output$res),
                      me_out_agg = mean(output_agg$res_mean),
                      mae_out = mean(abs(output$res)),
                      mae_out_agg = mean(abs(output_agg$res_mean)),
                      rmse_out = sqrt(mean(output$res^2)),
                      rmse_out_agg = sqrt(mean(output_agg$res_mean^2)))
  ## check if the direction of the univariate and multivariate analysis match
  compare <- coefs %>%
    merge(., Correlation_highest, by.x = "var", by.y = "var2") %>%
    mutate(dir = ifelse(value > 0, 1, -1)) %>%
    mutate(match = ifelse(dir * vals > 0, "match", "not")) %>%
    mutate(model = paste(List[i]))

  TempDf <- data.frame(res = residuals(mod), 
                       outcome = fulldata$LEQ,
                       Site = Site)
  TempDf <- merge(TempDf, Sit, by = "Site")
  
  colnames(TempDf) <- c("site",'res','outcome', 'Longitude', 'Latitude')
  
  ## Moran's I
  Moran <- tvsews::calc_moransI(TempDf, lon_col = "Longitude", lat_col = "Latitude", meas_col = "res")
  sample.metric.dists <- as.matrix(dist(cbind(TempDf$Longitude, TempDf$Latitude)))
  sample.metric.dists.inv <- 1/sample.metric.dists
  sample.metric.dists.inv[!is.finite(sample.metric.dists.inv)] <- 0
  MI_res_highest <- ape::Moran.I(TempDf$res, sample.metric.dists.inv)
  
  
  Temp <- data.frame(site = fulldata$Site, LEQ = fulldata$LEQ, Pred = predict(mod, fulldata, re.form=NA), Res = residuals(mod), 
                     me = mean(residuals(mod)), mae = mean(abs(residuals(mod))), 
                     MI = Moran, MI2 = MI_res_highest$observed, MI2_pv <- MI_res_highest$p.value, Model = List[i])
  
  FullPred <- rbind(FullPred, cbind(Ldata$Site, predict(mod, Ldata, re.form=NA), List[i])) 
  mod_coef <- rbind(mod_coef, coefs)
  mod_pred <- rbind(mod_pred, Temp)
  mod_stat <- rbind(mod_stat, compare)
  
  ##########################
  ## Run for loop for CV R2 and Adj R2 
  ##########################
  
  k_fold_groupings20 <- read.csv("~/CVGroups.csv")
  
  for (x in 1:20){ 
    Grp <- k_fold_groupings20[,c(1,4+x)]
    colnames(Grp) <- c("Site","fold5_grp")
    y_predicted_cv_temp <- data.frame()
    
    for (k in 1:5){ 
      Temp_dat2 <- merge(Grp, fulldata, by.all = "Site")
      CV_data <- Temp_dat2[Temp_dat2[, "fold5_grp"] != k, ]
      CV_pred_data <- Temp_dat2[Temp_dat2[, "fold5_grp"] == k, ] 
      
      Model_cv <- lmer(formula, data = CV_data) 
      results.val <- predict(Model_cv , CV_pred_data, re.form=NA)
      y_predicted_cv_temp <- rbind(y_predicted_cv_temp ,cbind(CV_pred_data$Site, CV_pred_data$LEQ, results.val, rep(k, times = length(results.val)), rep(x, times = length(results.val)), List[i]))
    }
    y_predicted_cv_temp$V2 <- as.numeric(y_predicted_cv_temp$V2)
    y_predicted_cv_temp$results.val <- as.numeric(y_predicted_cv_temp$results.val)
    y_predicted_cv <- rbind(y_predicted_cv, y_predicted_cv_temp)
    
    # each of the 20 fold 5 folds
    x20_Cv_r2 <- summary(lm(y_predicted_cv_temp$V2 ~ y_predicted_cv_temp$results.val))$r.squared
    x20_Cv_slope <- summary(lm(y_predicted_cv_temp$V2 ~ y_predicted_cv_temp$results.val))$coef[2]
    x20_Cv_Resid <- y_predicted_cv_temp$V2-y_predicted_cv_temp$results.val
    x20_Cv__mean_res <- mean(x20_Cv_Resid) # out of fold mean error
    x20_Cv__mae <- mean(abs(x20_Cv_Resid)) # out of fold mean abs error
    x20_Cv__mse <- mean(x20_Cv_Resid^2) # out of fold mean squared error
    x20_Cv__rmse <- sqrt(mean(x20_Cv_Resid^2)) # out of fold root mean squared error
    x20_Cv__meanSPR <- mean(x20_Cv_Resid)/sd(x20_Cv_Resid) # out of fold mean standardized prediction residuals
    x20_Cv__rms <- sqrt((mean(x20_Cv_Resid)/sd(x20_Cv_Resid))^2) # out of fold root mean square of the standardized prediction residuals
    
    # dataset of the x20_5 CV group results
    x20_Cv_dataset_ols <- rbind(x20_Cv_dataset_ols, cbind(x20_Cv_r2, x20_Cv_slope, x20_Cv__mean_res, x20_Cv__mae, 
                                                          x20_Cv__mse, x20_Cv__rmse, x20_Cv__meanSPR, x20_Cv__rms, x, List[i]))
  }
  
  Coef_highest <- data.frame(fixef(mod)) # variables in model and the coefficients
  Coef_highest$Model_subset <- List[i] # what model this is
  Coef_highest$OOF_Slope <- summary(lm(y_predicted_cv$V2 ~ y_predicted_cv$results.val))$coef[2] # out of fold slope
  Coef_highest$OOF_R2 <- summary(lm(y_predicted_cv$V2 ~ y_predicted_cv$results.val))$r.squared # out of fold R2
  OOF_Resid <- y_predicted_cv$obs-y_predicted_cv$results.val
  Coef_highest$OOF_mean_res <- mean(OOF_Resid) # out of fold mean error
  Coef_highest$OOF_mae <- mean(abs(OOF_Resid)) # out of fold mean abs error
  Coef_highest$OOF_mse <- mean(OOF_Resid^2) # out of fold mean squared error
  Coef_highest$OOF_rmse <- sqrt(mean(OOF_Resid^2)) # out of fold root mean squared error
  Coef_highest$OOF_meanSPR <- mean(OOF_Resid)/sd(OOF_Resid) # out of fold mean standardized prediction residuals
  Coef_highest$OOF_rms <- sqrt((mean(OOF_Resid)/sd(OOF_Resid))^2) # out of fold root mean square of the standardized prediction residuals
} 
colnames(y_predicted_cv) <- c('site', 'obs', 'pred', "fold5_grp5", "fold5_grp20", "model")

clipr::write_clip(mod_coef)
clipr::write_clip(mod_stat)
clipr::write_clip(mod_pred)
clipr::write_clip(x20_Cv_dataset_ols)
clipr::write_clip(y_predicted_cv)

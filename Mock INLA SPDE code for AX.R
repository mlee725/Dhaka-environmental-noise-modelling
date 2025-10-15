#########################################
#  Mock inla modelling script
#########################################

## This is mock code for a mock INLA-SPDE model for noise. This code was written for someone else so they could 
## use it as a reference for their noise LUR model. I pulled parts of this code from older scripts of mine and a number of online sources.

## packages ======
library(tvsews)
library(sf)
library(dplyr)
library(stringr)
library(ape)
library(INLA)

## data ======
# mock data
data <- data.frame(Site = c("A1", "B1", "C1", "A2", "B2", "C3"),
                   y = c(65, 72, 56, 66, 59, 62),
                   UTM_X = c(),
                   UTM_Y = c(),
                   green_50 = c(0,0,0,2,0,0),
                   green_100 = c(0,0,0,2,1,1),
                   ind_50 = c(0,0,0,2,3,6),
                   ind_100 = c(0,1,0,3,5,6),
                   road_50 = c(1,5,3,6,0,1),
                   road_100 = c(1,6,5,6,1,3),
                   ndvi_50 = c(0,0.5,0.3,0,0,1),
                   ndvi_100 = c(0.2,0.5,0.2,0.1,0,0.8))

Pred <- read.csv("~/prediction_surface_points.csv")

## inla model ======
# selected model predictors
Mod_preds <- c("road_50", "ndvi_100")

# -- create mesh -- #
# https://rpubs.com/jafet089/886687
coor <- cbind(data$Longitude, data$Latitude) # get coordinates
coor_sf <- st_as_sf(data, coords = c("UTM_X", "UTM_Y"), crs = 3772) # can use lat/long 

bound.outer = diff(range(st_coordinates(coor_sf)[,1]))/3
max.edge = diff(range(st_coordinates(coor_sf)[,1]))/(3*5)
mesh = inla.mesh.2d(loc=st_coordinates(coor_sf),
                    max.edge = c(1,2)*max.edge, # to avoid boundary effect
                    offset=c(max.edge, bound.outer),
                    cutoff = max.edge/5)

# -- spde -- #
######## this is the default, change as needed ######
spde1 <- inla.spde2.pcmatern(mesh, alpha = 2)

# assuming only one time point
s.index <- inla.spde.make.index("spatial.field", n.spde = spde1$n.spde) 
A_point <- inla.spde.make.A(mesh, loc = coor)  # The projector matrix automatically computes the weight vector for the neighborhood of each point and is calculated by providing the mesh and the locations of the datapoints

# I am assuming you want to output a prediction surface. If you don't, could use the coordinates from the sampling sites
# instead of the predictions grid (Pred) - this is a surface of points.
Pred_sf <- st_as_sf(Pred, coords = c("UTM_X", "UTM_Y"), crs = 3772)
Co <- st_coordinates(Pred_sf)
summary(Co)
Pred <- cbind(Pred, Co)
#Convert it to data frame
Pred_df <- as.data.frame(Pred_sf)
Pred_df <- cbind(Pred_df, Co)

Ap <- inla.spde.make.A(mesh = mesh, loc = Co)

# stack for the estimates
StackEsta <- inla.stack(data = list(y = data$y),               # the response variable
                        A = list(A_point, 1),                              # Then the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
                        effects = list(c(s.index, list(Intercept = 1)),    # The effects are organised in a list of lists. spatial effect and intercept first
                                       data[,Mod_preds]),
                        tag="Est")
# stack for the predictions
stackPred <- inla.stack(data = list(y = NA),  # NAs in the response variable  
                        A = list(A_point, 1),                              # Then the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
                        effects = list(c(s.index, list(Intercept = 1)),
                                       data[,Mod_preds]),
                        tag = "Pred")
# joint stack
StackEst <- inla.stack(StackEsta, stackPred)

Model_Formula <- as.formula(noquote(paste("y ~ -1 + Intercept + ", paste0(Mod_preds,  collapse = " + "), 
                                          " + f(spatial.field, model = spde1)")))
# fit the model
Model <- inla(Model_Formula, data = inla.stack.data(StackEst, spde = spde1), 
              control.predictor=list(A=inla.stack.A(StackEst), compute =TRUE),
              control.compute = list(cpo = T, dic = T, waic = T)) 

# -- output GRF surface -- #
proj <- inla.mesh.projector(mesh)
Mean <- inla.mesh.project(proj,Model$summary.random$spatial.field$mean) # random spatial effect (GRF) - mean
Sd <- inla.mesh.project(proj,Model$summary.random$spatial.field$sd) # random spatial effect (GRF) - sd

# We need to create spatial objects for the mean and variance of the GRF
xmean <- t(Mean)
xsd <- t(Sd)

xmean <- xmean[rev(1:length(xmean[,1])),]
xsd <- xsd[rev(1:length(xsd[,1])),]
xmean_ras <- raster(xmean,
                    xmn = range(proj$x)[1], xmx = range(proj$x)[2],
                    ymn = range(proj$y)[1], ymx = range(proj$y)[2],
                    crs = 3772)
xsd_ras <- raster(xsd,
                  xmn = range(proj$x)[1], xmx = range(proj$x)[2],
                  ymn = range(proj$y)[1], ymx = range(proj$y)[2],
                  crs = 3772)

# -- output prediction surface -- #
index.pred <- inla.stack.index(StackEst, "Pred")$data
Pred_data <- data.frame(X = Co[,1], 
                        Y = Co[,2],
                        Pred_mean = Model$summary.fitted.values[index.pred, "mean"],
                        Pred_ll = Model$summary.fitted.values[index.pred, "0.025quant"],
                        Pred_ul = Model$summary.fitted.values[index.pred, "0.975quant"])

post.mean.pred <- Model$summary.fitted.values[index.pred, "mean"]
post.sd.pred <- Model$summary.fitted.values[index.pred, "sd"]


# -- cv r2 -- #
k_fold_groupings20 <- data.frame(Site = c("A1", "B1", "C1", "A2", "B2", "C3"), # mock data
                                 G1 = c(3, 5, 1, 2, 4, 1),
                                 G2 = c(1, 4, 3, 5, 2, 1),
                                 G3 = c(4, 1, 5, 3, 2, 4),
                                 G4 = c(1, 2, 5, 4, 3, 1),
                                 G5 = c(5, 3, 1, 2, 4, 5),
                                 G6 = c(4, 2, 3, 1, 5, 4),
                                 G7 = c(2, 5, 1, 4, 3, 2),
                                 G8 = c(5, 1, 3, 2, 4, 5),
                                 G9 = c(3, 5, 2, 1, 4, 3),
                                 G10 = c(1, 4, 5, 2, 3, 1),
                                 G11 = c(4, 1, 3, 5, 2, 4),
                                 G12 = c(1, 3, 5, 4, 2, 1),
                                 G13 = c(2, 4, 1, 5, 3, 2),
                                 G14 = c(3, 1, 4, 2, 5, 3),
                                 G15 = c(5, 2, 1, 4, 3, 5),
                                 G16 = c(4, 1, 2, 5, 3, 4),
                                 G17 = c(1, 5, 3, 2, 4, 1),
                                 G18 = c(2, 4, 1, 5, 3, 2),
                                 G19 = c(3, 1, 5, 2, 4, 3),
                                 G20 = c(5, 3, 2, 1, 4, 5))

y_predicted_cv <- data.frame()

for (j in 1:20){
  Group <- k_fold_groupings20[,c(1,1+j)]
  colnames(Group) <- c("Site","fold5_grp")
  y_predicted_cv_temp <- data.frame()
  
  for (k in 1:5){
    Temp_dat2 <- merge(Group, data, by.all = "Site")
    CV_data <- Temp_dat2[Temp_dat2[, "fold5_grp"] != k, ]
    CV_pred_data <- Temp_dat2[Temp_dat2[, "fold5_grp"] == k, ] 
    
    spde1 <- inla.spde2.matern(mesh,
                               alpha = 2)
    # Define the weights
    A.train <- inla.spde.make.A(mesh=mesh, loc=cbind(CV_data$Longitude, CV_data$Latitude))
    A.val <- inla.spde.make.A(mesh=mesh, loc=cbind(CV_pred_data$Longitude, CV_pred_data$Latitude))
    
    s.index <- inla.spde.make.index("spatial.field", n.spde = spde1$n.spde) 
    
    stackTrain <- inla.stack(data = list(y = CV_data$y),               # the response variable
                             A = list(A.train, 1),                              # Then the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
                             effects = list(c(s.index, list(Intercept = 1)),    # The effects are organised in a list of lists. spatial effect and intercept first
                                            CV_data[,Mod_preds]),
                             tag="train")
    stackVal <- inla.stack(data = list(y = NA),  # NAs in the response variable  
                           A = list(A.val, 1),                              # Then the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
                           effects = list(c(s.index, list(Intercept = 1)),
                                          CV_pred_data[,Mod_preds]),
                           tag = "val")
    
    StackJoin <- inla.stack(stackTrain, stackVal)
    Model_Formula <- as.formula(noquote(paste("y ~ -1 + Intercept + ", paste0(Mod_preds,  collapse = " + "), 
                                              " + f(spatial.field, model = spde1)")))
    # fit the model
    
    Model_CV <- inla(Model_Formula, data = inla.stack.data(StackJoin, spde = spde1), 
                     control.predictor=list(A=inla.stack.A(StackJoin), compute =TRUE),
                     control.compute = list(cpo = T, dic = T, waic = T), verbose=TRUE
    )
    
    Coef_cv <- Model_CV$summary.fixed
    
    # Extract the fitted values
    index_inla_val = inla.stack.index(StackJoin,"val")$data
    results.val=Model_CV$summary.fitted$mean[index_inla_val]
    
    y_predicted_cv_temp <- rbind(y_predicted_cv_temp, 
                                 cbind(CV_pred_data$y,
                                       results.val, 
                                       rep(k, times = length(results.val)), 
                                       rep(j, times = length(results.val))))
    
  }
  y_predicted_cv <- rbind(y_predicted_cv, y_predicted_cv_temp)
  
}

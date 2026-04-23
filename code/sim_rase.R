############################################
#### LOADING ELLLIPSE EVOLUTION PACKAGE ####
############################################

library(devtools)
devtools::load_all("./package")

library(spatstat)
library(rase)

#########################
#### SIMULATING DATA ####
#########################

args <- commandArgs(trailingOnly = TRUE)
NUMBER <- args[1]
CONDITION <- args[2]

set.seed(NUMBER)
print(paste0("Number: ",NUMBER))
print(paste0("Condition: ",CONDITION))
name <- paste0(CONDITION,".sim.",NUMBER)
sim_directory <- paste0("./simulations/output/",name,"/")

sim_tree <- ape::read.tree(paste0(sim_directory,"true.tree.txt"))
sim_V <- read.csv(paste0(sim_directory,"true.V.tsv"),sep="\t")
sim_XY <- read.csv(paste0(sim_directory,"true.XY.tsv"),sep="\t")
sim_data <- cbind(sim_XY[1:length(sim_tree$tip.label),2:3],sim_V[1:length(sim_tree$tip.label),2:4])

owin_data <- c()
for (i in 1:nrow(sim_data)) {
  ellipse <- make_ellipse_coords(sim_data$x[i],sim_data$y[i],sim_data$r[i],sim_data$s[i],sim_data$a[i],10)
  owin_shape <- owin(xrange=c(min(ellipse$x),max(ellipse$x)),yrange=c(min(ellipse$y),max(ellipse$y)),poly=ellipse)
  owin_data[[i]] <- owin_shape
}

rase_data <- name.poly(owin_data,sim_tree,poly.names=sim_tree$tip.label)

# Can use polygons instead of ellipses (will need for empirical analyses)
# sp_data <- as_Spatial(spheno_shapes)
# owin_data <- shape.to.rase(sp_data)

##############
#### MCMC ####
##############

ace_x <- as.vector(ace(sim_data$x,sim_tree)$ace)
ace_y <- as.vector(ace(sim_data$y,sim_tree)$ace)
print(ace_x)
print(ace_y)
start_sigma2_x <- var(ace_x) / 40
start_sigma2_y <- var(ace_y) / 40
#start_sigma2_x <- 1
#start_sigma2_y <- 1
print(start_sigma2_x)
print(start_sigma2_y)
#ace_x <- ace(sim_data$x,sim_tree)$sigma2
#ace_y <- ace(sim_data$y,sim_tree)$sigma2
#print(paste0("Starting Values: sigma2_x = ",ace_x,", sigma2_y = ",ace_y))
#sigma2_scale <- 1.5 / length(sim_tree$tip.label)
sigma2_scale <- length(sim_tree$tip.label) / 10
params_init <- c(ace_x,ace_y,start_sigma2_x,start_sigma2_y)
#rase_data <- rase(sim_tree, rase_data, params0=params_init, logevery=10, niter=100000, sigma2_scale = 0.01, nGQ=20)
rase_data <- rase(sim_tree, rase_data, params0=params_init, logevery=10, niter=100000, sigma2_scale = sigma2_scale, nGQ=20)

write.csv(rase_data,paste0(sim_directory,"rase.csv"),row.names=FALSE,quote=FALSE)

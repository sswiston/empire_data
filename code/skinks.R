############################################
#### LOADING ELLLIPSE EVOLUTION PACKAGE ####
############################################

library(devtools)
devtools::load_all("./package")

######################
#### LOADING DATA ####
######################

args <- commandArgs(trailingOnly = TRUE)
NUMBER <- args[1]
CONDITION <- args[2]

print(paste0("Number: ",NUMBER))
print(paste0("Condition: ",CONDITION))

name <- paste0(CONDITION,".skinks.",NUMBER)

if (CONDITION=="under_prior") {
  testing_conditions <- c("under_prior")
  thinning <- 1
  print("Running under prior")
} else {
  testing_conditions <- c("None")
  thinning <- 10
}

print(paste0("Performing MCMC with Seed ",NUMBER))
set.seed(NUMBER)

example_tree <- ape::read.tree("./skinks/data/spheno.tre")
example_tree$tip.label <- gsub("_"," ",example_tree$tip.label)
data <- read.csv("./skinks/data/ellipse_data.csv",header=TRUE)

if (CONDITION=="shuffled") {
	print("Shuffling tips")
	data$taxon <- sample(data$taxon)
}

outpath <- paste0("./skinks/output/",name)

##############
#### MCMC ####
##############

tree <- dataTree$new(tree=example_tree,tip_data=data)
mcmc <- make_MCMC(tree,
                  prior_root_x=distributions3::Normal(mean(data$x),(max(data$x)-min(data$x))),
                  prior_root_y=distributions3::Normal(mean(data$y),(max(data$y)-min(data$y))),
                  prior_root_r=distributions3::Normal(mean(data$r),(max(data$r)-min(data$r))),
                  prior_root_s=distributions3::Normal(mean(data$s),(max(data$s)-min(data$s))),
                  prior_root_a=distributions3::Normal(mean(data$a),(max(data$a)-min(data$a))),
                  prior_sigma_x=distributions3::Uniform(0,10),
                  prior_sigma_y=distributions3::Uniform(0,10),
                  prior_sigma_r=distributions3::Uniform(0,10),
                  prior_sigma_s=distributions3::Uniform(0,10),
                  prior_sigma_a=distributions3::Uniform(0,40),
                  prior_mu=distributions3::Normal(mean(data$a),(max(data$a)-min(data$a))),
                  prior_kappa=distributions3::Uniform(0,10),
                  proposal_weights=list(sigma_x=5,
                                        sigma_y=5,
                                        sigma_r=5,
                                        sigma_s=5,
                                        sigma_a=5,
                                        root_x=3,
                                        root_y=3,
                                        root_r=3,
                                        root_s=3,
                                        root_a=3,
                                        mu=5,
                                        kappa=5,
                                        W_d=1.5,
                                        W_m=1.5,
                                        W_c=1.5,
                                        W_h=1.5,
                                        V_r=4,
                                        V_s=4,
                                        V_a=4,
                                        tip=1.5))

# Initializing parameters to reasonable values
mcmc$dataTree$param <- ellipseParam$new(sigma_x=1,
                                        sigma_y=1,
                                        sigma_r=1,
                                        sigma_s=1,
                                        sigma_a=1,
                                        kappa=1,
                                        mu=mean(data$a),
                                        root_x=mean(data$x),
                                        root_y=mean(data$y),
                                        root_r=mean(data$r),
                                        root_s=mean(data$s),
                                        root_a=mean(data$a))

run_MCMC(mcmc,iterations=2000000,moves_per_iteration=1,burnin=200000,thinning=thinning,testing_conditions=testing_conditions,filepath=outpath)

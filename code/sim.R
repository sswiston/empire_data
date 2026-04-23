############################################
#### LOADING ELLLIPSE EVOLUTION PACKAGE ####
############################################

library(devtools)
devtools::load_all("./package")

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
outpath <- paste0("./simulations/output/",name)

if (CONDITION=="under_prior") {
  testing_conditions <- c("under_prior")
  proposal_weights=list(sigma_x=1,sigma_y=1,sigma_r=1,sigma_s=1,sigma_a=1,root_x=1,root_y=1,root_r=1,root_s=1,root_a=1,mu=1,kappa=1,W_d=1,W_m=1,W_c=1,W_h=1,V_r=1,V_s=1,V_a=1,tip=1)
  tree_limits <- c(20,50)
  print("Running under prior")
}

if (CONDITION=="no_augmentation") {
  testing_conditions <- c("None")
  proposal_weights=list(sigma_x=1,sigma_y=1,sigma_r=3,sigma_s=3,sigma_a=3,root_x=3,root_y=3,root_r=3,root_s=3,root_a=3,mu=3,kappa=3,W_d=0,W_m=0,W_c=0,W_h=0,V_r=0,V_s=0,V_a=0,tip=0)
  tree_limits <- c(20,50)
  print("Running without data augmentation")
}

if (CONDITION=="small") {
  testing_conditions <- c("None")
  proposal_weights=list(sigma_x=5,sigma_y=5,sigma_r=5,sigma_s=5,sigma_a=5,root_x=3,root_y=3,root_r=3,root_s=3,root_a=3,mu=5,kappa=5,W_d=1,W_m=1,W_c=1,W_h=1,V_r=3,V_s=3,V_a=3,tip=1)
  tree_limits <- c(20,50)
  print("Simulating only small trees")
}

if (CONDITION=="full") {
  testing_conditions <- c("None")
  proposal_weights=list(sigma_x=5,sigma_y=5,sigma_r=5,sigma_s=5,sigma_a=5,root_x=3,root_y=3,root_r=3,root_s=3,root_a=3,mu=5,kappa=5,W_d=1.5,W_m=1.5,W_c=1.5,W_h=1.5,V_r=4,V_s=4,V_a=4,tip=1.5)
  tree_limits <- c(20,250)
  print("Simulating potentially large trees")
}

complete <- FALSE
while (complete == FALSE) {
  example_tree <- phytools::pbtree(b=.1,d=0,t=40,type="continuous")
  if (length(example_tree$tip.label) <= tree_limits[2]) {
    if (length(example_tree$tip.label) >= tree_limits[1]) {
      if (min(example_tree$edge.length) >= 0.01) {
        complete <- TRUE
      }
    }
  }
}
print(example_tree)

tree <- dataTree$new(tree=example_tree)
tree$sim_data()
tree$save(filepath=outpath,prefix="true")

##############
#### MCMC ####
##############

print(paste0("Performing MCMC with Seed ",NUMBER))
mcmc <- make_MCMC(tree,proposal_weights=proposal_weights)
run_MCMC(mcmc,iterations=2000000,moves_per_iteration=1,burnin=100000,thinning=10,testing_conditions=testing_conditions,filepath=outpath)

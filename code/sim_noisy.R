############################################
#### LOADING ELLLIPSE EVOLUTION PACKAGE ####
############################################

library(devtools)
devtools::load_all("./package")

library(ape)

######################
#### LOADING DATA ####
######################

args <- commandArgs(trailingOnly = TRUE)
NUMBER <- args[1]
CONDITION <- args[2]

set.seed(NUMBER)
print(paste0("Number: ",NUMBER))
print(paste0("Condition: ",CONDITION))
name <- paste0(CONDITION,".sim.",NUMBER)
sim_directory <- paste0("./simulations/output/",name,"/")
outpath <- paste0(sim_directory,"noisy")

sim_tree <- ape::read.tree(paste0(sim_directory,"true.tree.txt"))
n_taxa <- length(sim_tree$tip.label)
sim_V <- read.csv(paste0(sim_directory,"true.V.tsv"),sep="\t")
sim_XY <- read.csv(paste0(sim_directory,"true.XY.tsv"),sep="\t")
sim_data <- cbind(sim_XY[1:n_taxa,2:3],sim_V[1:n_taxa,2:4])
sim_data$taxon <- sim_tree$tip.label

# Adding noise to the tip data
x_error <- (max(sim_data$x) - min(sim_data$x)) * .1
y_error <- (max(sim_data$y) - min(sim_data$y)) * .1
r_error <- (max(sim_data$r) - min(sim_data$r)) * .1
s_error <- (max(sim_data$s) - min(sim_data$s)) * .1
a_error <- (max(sim_data$a) - min(sim_data$a)) * .1

noisy_data <- sim_data
for (i in 1:nrow(noisy_data)) {
  noisy_data$x[i] <- noisy_data$x[i] + runif(1, -1 * x_error, x_error)
  noisy_data$y[i] <- noisy_data$y[i] + runif(1, -1 * y_error, y_error)
  noisy_data$r[i] <- noisy_data$r[i] + runif(1, -1 * r_error, r_error)
  noisy_data$s[i] <- noisy_data$s[i] + runif(1, -1 * s_error, s_error)
  noisy_data$a[i] <- noisy_data$a[i] + runif(1, -1 * a_error, a_error)
}

tree <- dataTree$new(tree=sim_tree)
tree$tip_data <- noisy_data
tree$save(filepath=outpath,prefix="true")

##############
#### MCMC ####
##############

proposal_weights=list(sigma_x=5,sigma_y=5,sigma_r=5,sigma_s=5,sigma_a=5,root_x=3,root_y=3,root_r=3,root_s=3,root_a=3,mu=5,kappa=5,W_d=1,W_m=1,W_c=1,W_h=1,V_r=3,V_s=3,V_a=3,tip=1)

print(paste0("Performing MCMC with Seed ",NUMBER))
mcmc <- make_MCMC(tree,proposal_weights=proposal_weights)
run_MCMC(mcmc,iterations=2000000,moves_per_iteration=1,burnin=100000,thinning=10,filepath=outpath)

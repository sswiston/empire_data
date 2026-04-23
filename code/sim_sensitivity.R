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
outpath <- paste0(sim_directory,"sensitivity")

sim_tree <- ape::read.tree(paste0(sim_directory,"true.tree.txt"))
sim_V <- read.csv(paste0(sim_directory,"true.V.tsv"),sep="\t")
sim_XY <- read.csv(paste0(sim_directory,"true.XY.tsv"),sep="\t")
sim_data <- cbind(sim_XY[1:length(sim_tree$tip.label),2:3],sim_V[1:length(sim_tree$tip.label),2:4])
sim_data$taxon <- sim_tree$tip.label

# Sampling half (rounded up) of taxa to remove due to extinction
num_total <- length(sim_tree$tip.label)
num_removed <- floor(num_total/2)
removed_taxa <- sample(1:num_total,num_removed)
trimmed_tree <- drop.tip(sim_tree,removed_taxa)
trimmed_data <- sim_data[-c(removed_taxa),]

tree <- dataTree$new(tree=trimmed_tree)
tree$tip_data <- trimmed_data
tree$save(filepath=outpath,prefix="true")

##############
#### MCMC ####
##############

proposal_weights=list(sigma_x=5,sigma_y=5,sigma_r=5,sigma_s=5,sigma_a=5,root_x=3,root_y=3,root_r=3,root_s=3,root_a=3,mu=5,kappa=5,W_d=1,W_m=1,W_c=1,W_h=1,V_r=3,V_s=3,V_a=3,tip=1)

print(paste0("Performing MCMC with Seed ",NUMBER))
mcmc <- make_MCMC(tree,proposal_weights=proposal_weights)
run_MCMC(mcmc,iterations=2000000,moves_per_iteration=1,burnin=100000,thinning=10,filepath=outpath)

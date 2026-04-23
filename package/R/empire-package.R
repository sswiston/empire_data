#' @keywords internal
#' @title empire
#' @author Sarah Swiston, \email{sarah.k.swiston@@gmail.com}
#' @references stuff
#' @details
#'
#' Roxygen created this manual page on `r Sys.Date()` using R version
#' `r getRversion()`.
#'
#' The EMPIRE package is written in an *object-oriented* manner. Data is stored in objects with custom classes, and those objects have "methods" (functions) that can be applied. For example, you can use `dataTree$new()` to create a new dataTree object, or `mcmc$run()` to run your mcmc object.
#'
#' However, there are also a number of "wrapper" functions provided for your convenience. These wrapper functions allow you to interact with EMPIRE in a more conventional "functional" way. The tutorial examples use these wrapper functions.
#'
#' # Important Objects
#'
#' * [ellipseParam()] stores model parameter values
#' * [detailedTree()] stores a phylogenetic tree and associated information
#' * [dataTree()] stores a detailedTree and any simulated / attached data
#' * [ellipseMCMC()] stores a dataTree object and necessary values / methods for performing MCMC
#'
#'# Examples
#'
#' ## SIMULATION EXAMPLE
#'
#' An example which simulates a tree with data and tip states.
#'
#' Functions: [make_dataTree()], [sim_data()], [save_dataTree()]
#'
#' ```
#' # 1. Set directory for output files (modify for your file system)
#' outpath <- "./output"
#'
#' # 2. Simulate an example tree using phytools
#' set.seed(4)
#' example_tree <- phytools::pbtree(b=.05, d=0, t=40, type="continuous")
#'
#' # 3. Make this "phylo" tree into a dataTree (an EMPIRE data type)
#' tree <- make_dataTree(tree=example_tree)
#'
#' # 4. Simulate some data for this tree
#' tree <- sim_data(tree=tree)
#'
#' # 5. Save the tree and tip data
#' save_dataTree(tree=tree, filepath=outpath, prefix="true")
#'
#' # 6. Plot the tree and node labels
#' plot_label_tree(tree=tree, filepath=outpath)
#'
#' # 7. Plot the grid system for the simulated cladogenetic scenarios
#' plot_grid(model=tree)
#' ```
#'
#' ## ANALYSIS EXAMPLE
#'
#' An example which performs MCMC on the previously-simulated dataset. You will need:
#'
#' * A bifurcating, ultrametric phylogeny of type **phylo** (i.e. read in with the "ape" package)
#' * A **data.frame** with tip ellipse data; headers should be *taxon*, *x*, *y*, *r*, *s*, *a*.
#'
#' Make sure your taxon names match the tip labels of the phylogeny (can be in any order).
#'
#' Functions: [make_dataTree()], [make_MCMC()], [set_empirical_priors_MCMC()], [run_MCMC()]
#'
#' ```
#' # 1. Read tree and tip data from directory
#' outpath <- "./output"
#' example_tree <- ape::read.tree(paste0(outpath,"/true.tree.txt"))
#' data <- read.csv(file=paste0(outpath,"/true.tip_data.tsv"), sep="\t", header=TRUE)
#'
#' # 2. Turn tree and data into a dataTree (an EMPIRE data type)
#' tree <- make_dataTree(tree=example_tree, tip_data=data)
#'
#' # 3. Choose proposal weights (optional)
#' proposal_weights=list(sigma_x=5,
#'                       sigma_y=5,
#'                       sigma_r=5,
#'                       sigma_s=5,
#'                       sigma_a=5,
#'                       root_x=3,
#'                       root_y=3,
#'                       root_r=3,
#'                       root_s=3,
#'                       root_a=3,
#'                       mu=5,
#'                       kappa=5,
#'                       W_d=1.5,
#'                       W_m=1.5,
#'                       W_c=1.5,
#'                       W_h=1.5,
#'                       V_r=4,
#'                       V_s=4,
#'                       V_a=4,
#'                       tip=1.5)
#'
#' # 4. Make MCMC object
#' mcmc <- make_MCMC(dataTree=tree, proposal_weights=proposal_weights)
#'
#' # 5. Set priors for MCMC based on the range of the data
#' mcmc <- set_empirical_priors_MCMC(mcmc)
#'
#' # 6. Run MCMC
#' run_MCMC(mcmc=mcmc, iterations=1000, burnin=100, thinning=10, filepath=outpath)
#' ```
#'
#' ## POST-PROCESSING EXAMPLE
#'
#' An example for processing the output of an EMPIRE MCMC analysis.
#'
#' Functions: [reconstruct_ellipses()], [plot_scenario()]
#'
#' ```
#' # 1. Read tree and tip data from directory
#' outpath <- "./output"
#' example_tree <- ape::read.tree(paste0(outpath,"/true.tree.txt"))
#' data <- read.csv(file=paste0(outpath,"/true.tip_data.tsv"), sep="\t", header=TRUE)
#' tree <- make_dataTree(tree=example_tree, tip_data=data)
#'
#' # 2. Generate reconstructions for each internal node
#' reconstruction <- reconstruct_ellipses(outpath,burnin=10)
#' print(reconstruction)
#' # Note that it includes both indices and values for concentric circles (c, cval) and direction lines (h, hval)
#' # The VALUES are used for plotting ellipses (e.g. plot_scenario()) or making daughters (transform_node())
#'
#' # 3. Plot the phylogenetic tree on a map, with nodes located at their reconstructed centroids
#' plot_tree_map(tree,reconstruction)
#'
#' # 4. Plot the scenario at the MRCA -- note that the range of the blue (left) daughter is very small
#' mrca <- reconstruction[1,]
#' plot_scenario_annotated(x=mrca$x,y=mrca$y,r=mrca$r,s=mrca$s,a=mrca$a,d=mrca$d,m=mrca$m,cval=mrca$cval,hval=mrca$hval,z=10,alpha=-5)
#'
#' # 5. Plot the scenario at a different node (e.g. node 28)
#' node <- reconstruction[which(reconstruction$node==28),]
#' plot_scenario_annotated(x=node$x,y=node$y,r=node$r,s=node$s,a=node$a,d=node$d,m=node$m,cval=node$cval,hval=node$hval,z=10,alpha=-5)
#'
#' # 6. Plot a lineage, including cladogenetic scenarios, from the MRCA to a specific tip
#' plot_lineage(tree,reconstruction,1)
#'
#' # 7. Plot the posterior density for the ancestral and daughter ellipses alongside the best reconstruction at a set of nodes
#' plot_posterior(outpath,c(21,23))
#' ```
#'
#' @md
"_PACKAGE"
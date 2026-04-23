#' @title reconstruct_ellipses() function
#' @description Processes output of an ellipse MCMC and generates ancestral ellipses for each node
#' @param filepath An object of class "character" (a string) that gives the filepath where the output log files from your MCMC analysis are saved. Default behavior will look for the files in the current directory.
#' @param burnin The number of burnin *rows* produced during the MCMC analysis (burnin divided by thinning); will be removed before processing
#' @param cval An object of class "vector" containing any number of numeric elements indicating the values associated with different concentric circles relative to the size of the original circle. The default value includes 4 concentric circles: 0, 0.5, 1.0, and 1.5.
#' @param hval An object of class "vector" containing any number of numeric elements indicating the values associated with different direction lines in radians. The default value includes 8 direction lines: 0, pi/4, pi/2, 3pi/4, pi, 5pi/4, 3pi/2, and 7pi/4.
#' @export
reconstruct_ellipses <- function(filepath=".", burnin=0, cval=NULL, hval=NULL) {

  # Default values associated with c and h indices
  if (is.null(cval)) {cval=seq(0,1.5,.5)}
  if (is.null(hval)) {hval=2*pi*seq(0,7/8,1/8)}

  # Starting dataframe for ancestral ellipse reconstructions for each internal node
  reconstruction <- data.frame(matrix(NA,nrow=0,ncol=12))

  # Getting the number of rows in the MCMC output and saving as "length"
  data_d_pre <- read.table(paste0(filepath,"/d_log.tsv"),sep="\t",header=TRUE)
  length <- nrow(data_d_pre)

  # Reading the MCMC output for each ellipse character, removing burnin lines
  data_d <- read.table(paste0(filepath,"/d_log.tsv"),sep="\t",header=TRUE)[burnin:length,]
  data_m <- read.table(paste0(filepath,"/m_log.tsv"),sep="\t",header=TRUE)[burnin:length,]
  data_c <- read.table(paste0(filepath,"/c_log.tsv"),sep="\t",header=TRUE)[burnin:length,]
  data_h <- read.table(paste0(filepath,"/h_log.tsv"),sep="\t",header=TRUE)[burnin:length,]
  data_r <- read.table(paste0(filepath,"/r_log.tsv"),sep="\t",header=TRUE)[burnin:length,]
  data_s <- read.table(paste0(filepath,"/s_log.tsv"),sep="\t",header=TRUE)[burnin:length,]
  data_a <- read.table(paste0(filepath,"/a_log.tsv"),sep="\t",header=TRUE)[burnin:length,]
  data_x <- read.table(paste0(filepath,"/x_log.tsv"),sep="\t",header=TRUE)[burnin:length,]
  data_y <- read.table(paste0(filepath,"/y_log.tsv"),sep="\t",header=TRUE)[burnin:length,]

  # Making lists of all nodes, tip nodes, and internal nodes
  node_list <- c(1:(length(colnames(data_x)) - 1))
  tip_list <- c(1:ceiling(length(node_list) / 2))
  internal_list <- c((length(tip_list)+1):length(node_list))

  # Looping over each internal node to perform a reconstruction
  for (node in internal_list) {
    # Getting the name of the node in the output files
    # Equals the node number preceded by an X
    node_name <- paste0("X",node)
    # Getting the data columns associated with the node
    node_data <- data.frame(cbind(d=data_d[[node_name]],
                                  m=data_m[[node_name]],
                                  c=data_c[[node_name]],
                                  h=data_h[[node_name]],
                                  r=data_r[[node_name]],
                                  s=data_s[[node_name]],
                                  a=data_a[[node_name]],
                                  x=data_x[[node_name]],
                                  y=data_y[[node_name]]))

    # Getting the counts for each combination of discrete values for cladogenetic scenarios
    scenario_counts <- dplyr::summarize(.data=node_data,count=dplyr::n(),.by=c(d,m,c,h))
    # Getting the most common scenario(s)
    ind <- which(scenario_counts$count == max(scenario_counts$count))
    scenario <- scenario_counts[ind,]
    # Randomly choose a scenario if multiple are equally preferred
    if (nrow(scenario) > 1) {scenario <- scenario[sample(1:nrow(scenario),1),]}
    # Getting the individual values associated with the best cladogenetic scenario at the node
    scenario_d <- scenario[[1]]
    scenario_m <- scenario[[2]]
    scenario_c <- scenario[[3]]
    scenario_h <- scenario[[4]]

    # Getting the actual values associated with the c and h indices
    scenario_cval <- cval[scenario_c + 1]
    scenario_hval <- hval[scenario_h + 1]

    # Getting the MCMC iterations consistent with the chosen scenario
    scenario_data <- node_data[which(node_data$d==scenario_d & node_data$m==scenario_m & node_data$c==scenario_c & node_data$h==scenario_h),]
    # Getting the mean values for continuous characters under that scenario
    scenario_r <- mean(scenario_data$r)
    scenario_s <- mean(scenario_data$s)
    scenario_a <- mean(scenario_data$a)
    scenario_x <- mean(scenario_data$x)
    scenario_y <- mean(scenario_data$y)

    # Adding this reconstruction to the dataframe of reconstructions
    row <- c(node,
             scenario_d,
             scenario_m,
             scenario_c,
             scenario_cval,
             scenario_h,
             scenario_hval,
             scenario_r,
             scenario_s,
             scenario_a,
             scenario_x,
             scenario_y)
    reconstruction <- rbind(reconstruction,row)
  }

  # Returning a dataframe of ancestral ellipse reconstructions
  colnames(reconstruction) <- c("node","d","m","c","cval","h","hval","r","s","a","x","y")
  return(reconstruction)
}
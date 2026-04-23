#' @title get_prob_w() function
#' @description Accepts a dataTree object and returns the log likelihood of observing the data-augmented cladogenetic scenario W
#' @param dataTree A dataTree object containing a detailedTree object, an ellipseParam object, and W, a dataframe of data-augmented values for cladogenetic scenarios at nodes
#' @export
get_prob_w <- function(dataTree) {
  # Checking input
  if (!("ellipseParam" %in% class(dataTree$param))) {stop("param object must be of class ellipseParam",call.=FALSE)}
  if (!("data.frame" %in% class(dataTree$W))) {stop("W must be of class data.frame",call.=FALSE)}
  # Calculating probability
  # Since this is a log probability, it is computed through sums rather than products
  # Probabilities for each data-augmented discrete value at each node are obtained, then summed (values of m,c,h are independent)
  p_m <- sum(log(dataTree$param$rho_m[dataTree$W$m+1]))
  p_c <- sum(log(dataTree$param$rho_c[dataTree$W$c+1]))
  p_h <- sum(log(dataTree$param$rho_h[dataTree$W$h+1]))
  prob <- sum(p_m,p_c,p_h)
  # THIS IS A LOG PROBABILITY
  return(prob)
}

#' @title get_prob_v() function
#' @description Accepts a dataTree object and returns the log likelihood of observing the data-augmented values V at nodes based on the history of cladogenetic scenarios W
#' @param dataTree A dataTree object containing a detailedTree object, and ellipseParam object, and V, a dataframe of data-augmented values (r,s,a) at nodes
#' @export
get_prob_v <- function(dataTree) {
  # Checking input
  if (!("ellipseParam" %in% class(dataTree$param))) {stop("param object must be of class ellipseParam",call.=FALSE)}
  if (!("data.frame" %in% class(dataTree$V))) {stop("V must be of class data.frame",call.=FALSE)}
  # If any of the rate parameters for V (r,s,a) are equal to 0, then the log probability is -Inf (0 rates are not allowed)
  find_zero = function(x1,x2,x3,x4) {
    if (x1 <= 0) {return(TRUE)}
    if (x2 <= 0) {return(TRUE)}
    if (x3 <= 0) {return(TRUE)}
    if (x4 <= 0) {return(TRUE)}
    return(FALSE)
  }
  if (find_zero(dataTree$param$sigma_r,dataTree$param$sigma_s,dataTree$param$sigma_a,dataTree$param$kappa)) {prob <- -Inf} # rate parameters are sigma_r, sigma_s, sigma_a, kappa
  else {
    # Calculating probability
    # Getting parameters
    sigma_r <- dataTree$param$sigma_r
    sigma_s <- dataTree$param$sigma_s
    sigma_a <- dataTree$param$sigma_a
    mu <- dataTree$param$mu
    kappa <- dataTree$param$kappa
    # Loop through all edges in the system
    prob <- 0
    for (i in 1:length(dataTree$tree$edges[,1])) {
      # Get edge information
      start_node <- dataTree$tree$edges[[i,1]]
      end_node <- dataTree$tree$edges[[i,2]]
      length <- dataTree$tree$edges[[i,3]]
      children <- dataTree$tree$edges[which(dataTree$tree$edges[,1]==start_node),2]
      daughter <- NULL
      if (end_node==children[1]) {daughter <- "right"} else {daughter <- "left"}
      # Get ancestor values
      anc <- dataTree$V[start_node,2:4]
      # Get cladogenetic scenario
      w <- dataTree$W[which(dataTree$W$node==start_node),]
      # Get starting values from ancestor values and scenario
      start <- anc
      start_a <- transform_node(char="a",daughter=daughter,d=w$d,m=w$m,cval=w$cval,a=anc$a,alpha=dataTree$param$alpha)
      start$a <- start_a
      # Get ending values
      end <- dataTree$V[end_node,2:4]
      # Getting branch probability
      branch_prob_r <- dnorm(end$r, start$r, sigma_r*sqrt(length), log=TRUE)
      branch_prob_s <- dnorm(end$s, start$s, sigma_s*sqrt(length), log=TRUE)
      branch_prob_a <- dnorm(end$a, mu + (start$a - mu) * exp(-1 * kappa * length), sqrt(sigma_a^2 / (2 * kappa) * (1 - exp(-2 * kappa * length))),log=TRUE)
      # Getting total probability (sum of branches)
      prob <- prob+branch_prob_r+branch_prob_s+branch_prob_a
    }
  }
  # THIS IS A LOG PROBABILITY
  return(prob)
}

#' @title get_expectations() function
#' @description Accepts a dataTree object and returns the expected values of x and y values at the tips based on the root values, tree structure, and data-augmented values
#' @param dataTree A dataTree object containing a detailedTree object and dataframes V and W containing values for cladogenetic scenarios at nodes
#' @export
get_expectations <- function(tree,param,XY,W,V) {
  # Checking input
  if (!("detailedTree" %in% class(tree))) {stop("tree object must be of class detailedTree",call.=FALSE)}
  if (!("ellipseParam" %in% class(param))) {stop("param object must be of class ellipseParam",call.=FALSE)}
  if (!("data.frame" %in% class(XY))) {stop("XY must be of class data.frame",call.=FALSE)}
  if (!("data.frame" %in% class(W))) {stop("W must be of class data.frame",call.=FALSE)}
  if (!("data.frame" %in% class(V))) {stop("V must be of class data.frame",call.=FALSE)}
  # Initializing tip expectations
  # These will start at the root value and be updated through tree traversal according to cladogenetic scenarios
  n_tips <- tree$n_tips
  expectations <- matrix(NA,nrow=n_tips,ncol=2)
  for (i in 1:n_tips) {
    expectations[i,1] <- XY[(n_tips+1),2]
    expectations[i,2] <- XY[(n_tips+1),3]
  }
  # Looping over nodes
  for (i in (n_tips+1):(2*n_tips-1)) {
    # Getting children
    children <- tree$edges[which(tree$edges[,1]==i),2]
    right_child <- children[1]
    left_child <- children[2]
    # Getting cladogenetic scenario
    w <- W[which(W$node==i),]
    v <- V[which(V$node==i),]
    # Getting appropriate transformations for each child
    for (j in 1:length(expectations[,1])) {
      if (tree$tip_descendants[left_child,j]==1) {
        # Applying transformation to all tip descendants
        expectations[j,1] <- transform_node(char="x",daughter="left",d=w$d,m=w$m,cval=w$cval,hval=w$hval,x=expectations[j,1],r=v$r,s=v$s,a=v$a,z=param$z)
        expectations[j,2] <- transform_node(char="y",daughter="left",d=w$d,m=w$m,cval=w$cval,hval=w$hval,y=expectations[j,2],r=v$r,s=v$s,a=v$a,z=param$z)
      }
      if (tree$tip_descendants[right_child,j]==1) {
        # Applying transformation to all tip descendants
        expectations[j,1] <- transform_node(char="x",daughter="right",d=w$d,m=w$m,cval=w$cval,hval=w$hval,x=expectations[j,1],r=v$r,s=v$s,a=v$a,z=param$z)
        expectations[j,2] <- transform_node(char="y",daughter="right",d=w$d,m=w$m,cval=w$cval,hval=w$hval,y=expectations[j,2],r=v$r,s=v$s,a=v$a,z=param$z)
      }
    }
  }
  return(expectations)
}

#' @title get_prob_xy() function
#' @description Accepts a dataTree object and returns the log likelihood of observing the x and y values at the tips based on the root values, tree structure, and data-augmented values
#' @param dataTree A dataTree object containing a detailedTree object, an ellipseParam object, and dataframes XY, V, and W containing values for cladogenetic scenarios at nodes
#' @export
get_prob_xy <- function(dataTree) {
  # If any of the rate parameters for XY are equal to 0, then the log probability is -Inf (0 rates are not allowed)
  find_zero = function(x1,x2) {
    if (x1 <= 0) {return(TRUE)}
    if (x2 <= 0) {return(TRUE)}
    return(FALSE)
  }
  if (find_zero(dataTree$param$sigma_x,dataTree$param$sigma_y)) {prob <- -Inf}
  else {
    # Calculating probability
    n_tips <- dataTree$tree$n_tips
    # Getting tip expectations
    expectations <- get_expectations(dataTree$tree,dataTree$param,dataTree$XY,dataTree$W,dataTree$V)
    # Likelihood calculation
    # This uses the phylogenetic variance/covariance matrix and rates of evolution sigma_x, sigma_y
    inv <- dataTree$tree$inv_matrix
    matrix <- dataTree$tree$phy_matrix
    matrix_det <- unlist(determinant(matrix,logarithm=TRUE))[[1]]
    x_data <- t(as.matrix(dataTree$XY$x))[1:n_tips]
    x_expected <- t(as.matrix(expectations[,1]))
    x <- x_data - x_expected
    y_data <- t(as.matrix(dataTree$XY$y))[1:n_tips]
    y_expected <- t(as.matrix(expectations[,2]))
    y <- y_data - y_expected
    sigma_x <- dataTree$param$sigma_x
    sigma_y <- dataTree$param$sigma_y
    # Probabilities for Brownian motion on tree (x and y independently)
    prob_x <- (-1/2) * (length(x) * log(2 * pi) + length(x) * 2 * log(sigma_x) + matrix_det + x %*% ((1 / sigma_x^2) * inv) %*% t(x))
    prob_y <- (-1/2) * (length(y) * log(2 * pi) + length(y) * 2 * log(sigma_y) + matrix_det + y %*% ((1 / sigma_y^2) * inv) %*% t(y))
    # Total probability for x and y
    prob <- prob_x[1,1] + prob_y[1,1]
  }
  # THIS IS A LOG PROBABILITY
  return(prob)
}

#' @title get_prob() function
#' @description Accepts a dataTree object and returns the log likelihood of observing the data augmented history (W,V) and tip values (x and y). Optionally, any of these probabilities can be provided (prob_w, prob_v, prob_xy) and will not be re-computed. If all probabilities are provided, a dataTree is not required.
#' @param dataTree A dataTree object containing a detailedTree object, and ellipseParam object, and dataframes V, W, and XY. If a probability is provided for a particular set of values (W, V, or XY), the dataFrame(s) required to calculate this probability may be missing from the dataTree without consequence.
#' @param prob_w Optionally, the log likelihood of observing the data-augmented cladogenetic scenarios at nodes
#' @param prob_v Optionally, the log likelihood of observing the data-augmented values at nodes
#' @param prob_xy Optionally, the log likelihood of observing the x and y values at the tips based on the root values, tree structure, and data-augmented values
#' @export
get_prob <- function(dataTree=NULL,prob_w=NULL,prob_v=NULL,prob_xy=NULL) {
  if (!("ellipseParam" %in% class(dataTree$param))) {stop("param object must be of class ellipseParam",call.=FALSE)}

  # Little function to check if any of the sigmas are 0
  # This should not happen (proposal should be automatically rejected)
  # But just in case, will return a -Inf log probability
  find_zero = function(x1,x2,x3,x4,x5) {
    if (x1 <= 0) {return(TRUE)}
    if (x2 <= 0) {return(TRUE)}
    if (x3 <= 0) {return(TRUE)}
    if (x4 <= 0) {return(TRUE)}
    if (x5 <= 0) {return(TRUE)}
    return(FALSE)
  }
  if (find_zero(dataTree$param$sigma_x,dataTree$param$sigma_y,dataTree$param$sigma_r,dataTree$param$sigma_s,dataTree$param$sigma_a)) {prob <- -Inf}
  else {
    # Type checking user input probability, if provided
    if(!("numeric" %in% class(prob_w))) {
      # If no user provided probability, type checking the dataTree
      if(!("dataTree" %in% class(dataTree))) {stop("get_prob() requires prob_w or a dataTree",call.=FALSE)}
      # If no user provided probability, calculates it from the dataTree
      prob_w <- get_prob_w(dataTree)
    }
    if(!("numeric" %in% class(prob_v))) {
      if(!("dataTree" %in% class(dataTree))) {stop("get_prob() requires prob_v or a dataTree",call.=FALSE)}
      prob_v <- get_prob_v(dataTree)
    }
    if(!("numeric" %in% class(prob_xy))) {
      if(!("dataTree" %in% class(dataTree))) {stop("get_prob() requires prob_xy or a dataTree",call.=FALSE)}
      prob_xy <- get_prob_xy(dataTree)
    }
    prob <- prob_w + prob_v + prob_xy
  }
  # THIS IS A LOG PROBABILITY
  return(prob)
}
# dataTree object ---------------------------------------------------------

#' @title dataTree object
#' @description This class stores an object of class "detailedTree" alongside an object of class "ellipseParam" containing model parameters and data objects (W, V, XY) containing simulated or data-augmented values at nodes. May take a long time to generate for larger trees.
#' @details This object is used during simulation to contain the simulated history of all cladogenetic events (W) and continuous values (V, XY) at all nodes (including tips). During inference, the data-augmented values for cladogenetic scenarios (W) and continuous parameters (V) will be updated, along with root states for (XY). The (XY) values for internal nodes will not be used during inference, and will therefore be set to "NA". Some elements of the ellipseParam object in the (param) field will also be updated during inference.
#' @field tree An object of class "detailedTree"
#' @field param An object of type ellipseParam containing information about the model parameters associated with a simulation or MCMC iteration
#' @field W An object of class "data.frame" containing the simulated or data-augmented cladogenetic scenario assigned to each non-tip node
#' @field V An object of class "data.frame" containing the simulated or data-augmented values for parameters 'r', 's', and 'a' assigned to each internal node
#' @field XY An object of class "data.frame" containing information about the 'x' and 'y' values. During simulation, this field will contain the full history of 'x' and 'y' values at each node. During inference, non-root internal nodes will be ignored.
#' @field tip_data An object of class "data.frame" containing taxon names "taxon" and the true values of characters 'x', 'y', 'r', 's', and 'a' observed at the tips of the phylogeny (extant taxa). The taxa names must exactly match those provided on the tree, although the order may be different. This data will be used to inform the V and XY dataframes -- if tip_data is not "NULL", simulating or setting V/XY will only apply to internal nodes and will not overwrite the characters at the tips. Due to the structure of the projection equation, extremely large or small values for "a" (log area) will make it impossible to calculate a model likelihood. Limits of [-50,50] will be put on these values, corresponding to true areas [1.93e-22,5.18e21]. All too-small ranges will be resized to have an "a" value of -50. If any too-large ranges exist, the dataTree object will produce an error, and you should consider using different units for your data.
dataTree <- R6::R6Class("dataTree",

private = list(

  # Initializing private fields (correspond to active bindings)
  # These fields (starting with '.') are for accessing information internally
  # Active bindings with corresponding names are for accessing information externally
  .tree = NULL,
  .param = NULL,
  .W = NULL,
  .V = NULL,
  .XY = NULL,
  .tip_data = NULL

),

active = list(

  # Active bindings that return values (stored in private fields) when called without assigning a value, e.g. dataTree$tree
  # Can also be set by the user (after type checking), e.g. dataTree$tree <- input_tree
  tree = function(value) {
    if (missing(value)) {private$.tree}
    else {
      if ("detailedTree" %in% class(value)) {
        private$.tree = value
      } else if ("phylo" %in% class(value)) {
        # If an object of class 'phylo' is provided instead of a 'detailedTree', it will be converted to a 'detailedTree' automatically
        private$.tree = detailedTree$new(value)
      } else {
        stop("tree must be of class detailedTree or phylo",call.=FALSE)
      }
    }
  },

  param=function(value) {
    if (missing(value)) {private$.param}
    else {
      if (!("ellipseParam" %in% class(value))) {
        stop("param object must be of class ellipseParam",call.=FALSE)
      }
      private$.param = value
    }
  },

  W=function(value) {
    if (missing(value)) {private$.W}
    else {
      if (!("data.frame" %in% class(value))) {
        stop("W must be of class data.frame",call.=FALSE)
      }
      if (any(is.na(value))) {
        stop("W must not contain NA or NaN values",call.=FALSE)
      }
      private$.W = value
    }
  },

  V=function(value) {
    if (missing(value)) {private$.V}
    else {
      if (!("data.frame" %in% class(value))) {
        stop("V must be of class data.frame",call.=FALSE)
      }
      if (any(is.na(value))) {
        stop("V must not contain NA or NaN values",call.=FALSE)
      }
      private$.V = value
    }
  },

  XY=function(value) {
    if (missing(value)) {private$.XY}
    else {
      if (!("data.frame" %in% class(value))) {
        stop("XY must be of class data.frame",call.=FALSE)
      }
      if (any(is.na(value))) {
        stop("XY must not contain NA or NaN values",call.=FALSE)
      }
      private$.XY = value
    }
  },

  tip_data=function(value) {
    if (missing(value)) {private$.tip_data}
    else {
      if (!("data.frame" %in% class(value))) {
        stop("tip_data must be of class data.frame",call.=FALSE)
      }
      if (any(is.na(value))) {
        stop("tip_data must not contain NA or NaN values",call.=FALSE)
      }
      if (is.null(self$tree)) {
        # A tree must be input before adding tip data because the tree is used to properly order the tip data
        stop("must input a tree before adding tip data",call.=FALSE)
      }
      # Setting up a dataframe for the tip data
      # This is a new dataframe that will store the input data inside the dataTree object
      temp_data <- data.frame(matrix(nrow=length(self$tree$tree$tip.label),ncol=6))
      colnames(temp_data) <- c("taxon","x","y","r","s","a")
      # Looping through the tip labels on the tree -- this is for reordering the tip data according to the tree
      for (i in 1:length(self$tree$tree$tip.label)) {
        # Getting the taxon name for the tip
        label <- self$tree$tree$tip.label[i]
        # Getting the index of the input data row containing data for a particular tip
        index <- which(value$taxon == label)
        # Making sure that data exists for this taxon
        if (!(length(index)==1)) {stop(paste0("missing data for taxon ",label," , please make sure all tree tips have corresponding ellipse data with matching taxon names"),call.=FALSE)}
        # Getting values for this tip and adding them to the dataframe
        temp_data$taxon[i] <- label
        temp_data$x[i] <- value$x[index]
        temp_data$y[i] <- value$y[index]
        temp_data$r[i] <- value$r[index]
        temp_data$s[i] <- value$s[index]
        temp_data$a[i] <- value$a[index]
      }
      # Checking that the input 'a' values do not have an absolute value >50, this will cause underflow/overflow issues
      # If this problem is encountered, the data should be rescaled
      for (i in 1:nrow(temp_data)) {
        if (as.numeric(temp_data$a[i]) <= -50) {
          print("Some values for 'a' (log area) are less than -50, which will cause problems for the projection equation -- these are being set to -50")
          temp_data$a[i] <- -50
        }
        if (as.numeric(temp_data$a[i]) >= 50) {
          stop("Some values for 'a' (log area) are greater than 50, which will cause problems for the projection equation -- please consider using different units for your data",call.=FALSE)
        }
      }
      private$.tip_data = temp_data
      # Adding to corresponding containers
      # V stores r,s,a values at nodes, including tips
      V <- data.frame(matrix(NA,nrow=self$tree$n_nodes,ncol=4))
      colnames(V) <- c("node","r","s","a")
      for (i in 1:length(self$tree$tree$tip.label)) {
        V$node[i] <- i
        V$r[i] <- self$tip_data$r[i]
        V$s[i] <- self$tip_data$s[i]
        V$a[i] <- self$tip_data$a[i]
      }
      private$.V <- V
      # XY stores x,y values at nodes, including tips
      XY <- data.frame(matrix(NA,nrow=self$tree$n_nodes,ncol=3))
      colnames(XY) <- c("node","x","y")
      for (i in 1:length(self$tree$tree$tip.label)) {
        XY$node[i] <- i
        XY$x[i] <- self$tip_data$x[i]
        XY$y[i] <- self$tip_data$y[i]
      }
      private$.XY <- XY
    }
  }
),

public = list(

  #' @description Creates a new instance of the dataTree class
  #' @param tree An object of class "detailedTree"
  #' @param param An object of type ellipseParam containing information about the model parameters associated with a simulation or MCMC iteration
  #' @param W An object of class "data.frame" containing the simulated or data-augmented cladogenetic scenario assigned to each non-tip node
  #' @param V An object of class "data.frame" containing the simulated or data-augmented values for parameters 'r', 's', and 'a' assigned to each internal node
  #' @param XY An object of class "data.frame" containing information about the 'x' and 'y' values. During simulation, this field will contain the full history of 'x' and 'y' values at each node. During inference, non-root internal nodes will be ignored.
  #' @param tip_data An object of class "data.frame" containing the true values of characters 'x', 'y', 'r', 's', and 'a' observed at the tips of the phylogeny (extant taxa). The taxa names must exactly match those provided on the tree, although the order may be different. This data will be used to inform the V and XY dataframes -- if tip_data is not "NULL", simulating or setting V/XY will only apply to internal nodes and will not overwrite the characters at the tips. When tip-specific moves (inverting r & s values) are performed, the values in tip_data will not be changed, but the values in V will be changed.
  initialize = function(tree=NULL,
                        param = NULL,
                        W = NULL,
                        V = NULL,
                        XY = NULL,
                        tip_data = NULL) {

    # Making sure there is a tree and it is the correct type (detailedTree) or a type that can be converted (phylo)
    # A tree is REQUIRED in order to create the dataTree object
    if ("detailedTree" %in% class(tree)) {
      self$tree <- tree
    } else if ("phylo" %in% class(tree)) {
      tree = detailedTree$new(tree)
      self$tree <- tree
    }

    # Object can be created without setting these fields (default=NULL)
    # This makes it possible to simulate parameters and data for a particular tree
    # This is important for initializing MCMC on real data
    if (!(is.null(param))) {self$param=param}
    if (!(is.null(W))) {self$W=W}
    if (!(is.null(V))) {self$V=V}
    if (!(is.null(XY))) {self$XY=XY}
    if (!(is.null(tip_data))) {self$tip_data=tip_data}
  },

  #' @description Simulates values for the W field of the dataTree object
  #' @details The W object contains information about the simulated or data-augmented cladogenetic scenarios at each internal node. This includes the daughter assignment 'd', the cladogenetic mode 'm', the concentric circle 'c', and the direction line 'h'.
  sim_W = function() {

    # checking if tree exists and is of class "detailedTree"
    if (!("detailedTree" %in% class(self$tree))) {stop("tree must be of class detailedTree",call.=FALSE)}
    # checking if param exists and is of class "ellipseParam"
    if (!("ellipseParam" %in% class(self$param))) {stop("param object must be of class ellipseParam",call.=FALSE)}

    # starting W
    W <- data.frame(matrix(NA,nrow=(self$tree$n_tips-1),ncol=7))
    colnames(W) <- c("node","d","m","c","h","cval","hval")
    W$node <- (self$tree$root_idx):self$tree$n_nodes

    # simulating cladogenetic scenarios
    # Selecting discrete values of d,m,c,h
    W$d <- sample(c(0,1),length(W$node),replace=TRUE)
    W$m <- sample(c(0,1),length(W$node),prob=self$param$rho_m,replace=TRUE)
    W$c <- sample(0:(length(self$param$cval)-1),length(W$node),prob=self$param$rho_c,replace=TRUE)
    W$h <- sample(0:(length(self$param$hval)-1),length(W$node),prob=self$param$rho_h,replace=TRUE)
    # Getting associated quantitative values for c and h
    for (i in 1:(self$tree$n_tips-1)) {
      W$cval[i] <- self$param$cval[W$c[i]+1]
      W$hval[i] <- self$param$hval[W$h[i]+1]
    }

    # setting W for dataTree object
    self$W = W
  },

  #' @description Simulates values for the V field of the dataTree object
  #' @details The V object contains information about the simulated or data-augmented values of continuous characters 'r', 's', and 'a' at each internal node.
  sim_V = function() {

    # ensuring that tree, param, and W exist and are of correct classes
    if (!("detailedTree" %in% class(self$tree))) {stop("Tree must be of class detailedTree",call.=FALSE)}
    if (!("ellipseParam" %in% class(self$param))) {stop("Param must be of class ellipseParam",call.=FALSE)}
    if (!("data.frame" %in% class(self$W))) {stop("W must be of class data.frame",call.=FALSE)}

    # starting V
    V <- data.frame(matrix(NA,nrow=self$tree$n_nodes,ncol=4))
    colnames(V) <- c("node","r","s","a")
    V$node <- 1:self$tree$n_nodes
    # Making sure the root r,s,a values match the parameters in the ellipseParam object
    V$r[self$tree$root_idx] <- self$param$root_r
    V$s[self$tree$root_idx] <- self$param$root_s
    V$a[self$tree$root_idx] <- self$param$root_a

    # function to simulate a node's daughters
    sim_iter_rsa <- function(V,node,edges,W,sigma_r,sigma_s,sigma_a,mu,kappa,alpha) {
      # Checking if the node is a tip (and we don't need to do any more work on this node)
      n_tips <- length(V[,1])-length(W[,1])
      if (node %in% 1:n_tips) {
        return(V)
      }
      # Getting cladogenetic scenario
      w <- W[which(W$node==node),]
      # Getting the values for the current node
      current_vals <- V[node,2:4]
      left_start <- current_vals
      right_start <- current_vals
      current_a <- current_vals[[3]]
      # Getting appropriate starts for each branch
      left_start[3] <- transform_node(char="a",daughter="left",d=w$d,m=w$m,cval=w$cval,a=current_a,alpha=alpha)
      right_start[3] <- transform_node(char="a",daughter="right",d=w$d,m=w$m,cval=w$cval,a=current_a,alpha=alpha)
      # Getting children and branch lengths
      children <- edges[which(edges[,1]==node),2]
      right_child <- children[1]
      left_child <- children[2]
      lengths <- edges[which(edges[,1]==node),3]
      right_length <- lengths[1]
      left_length <- lengths[2]
      # Simulating a Brownian motion for x and y over each branch (left and right)
      left_vals <- c(rnorm(1, left_start[[1]], (sigma_r * sqrt(left_length))),
                     rnorm(1, left_start[[2]], (sigma_s * sqrt(left_length))),
                     rnorm(1, mu + (left_start[[3]] - mu) * exp(-1 * kappa * left_length), sqrt(sigma_a^2 / (2 * kappa) * (1 - exp(-2 * kappa * left_length)))))
      right_vals <- c(rnorm(1, right_start[[1]], sigma_r * sqrt(right_length)),
                      rnorm(1, right_start[[2]], sigma_s * sqrt(right_length)),
                      rnorm(1, mu + (right_start[[3]] - mu) * exp(-1 * kappa * right_length), sqrt(sigma_a^2 / (2 * kappa) * (1 - exp(-2 * kappa * right_length)))))
      # Assigning new node values to descendant nodes (left and right)
      V[left_child,2] <- left_vals[1]
      V[left_child,3] <- left_vals[2]
      V[left_child,4] <- left_vals[3]
      V[right_child,2] <- right_vals[1]
      V[right_child,3] <- right_vals[2]
      V[right_child,4] <- right_vals[3]
      return(V)
    }

    # function to iterate over nodes
    sim_rsa <- function(V,node,nodes,edges,W,sigma_r,sigma_s,sigma_a,mu,kappa,alpha) {
      # calling sim_iter_rsa on a node, which performs a simulation for one node and the branches that connect to its right and left daughters
      V <- sim_iter_rsa(V,node,edges,W,sigma_r,sigma_s,sigma_a,mu,kappa,alpha)
      # if a left daughter exists, recurring on that daughter
      if (!is.na(nodes[node,3])) {
        V <- sim_rsa(V,nodes[node,3],nodes,edges,W,sigma_r,sigma_s,sigma_a,mu,kappa,alpha)
      }
      # if a right daughter exists, recurring on that daughter
      if (!is.na(nodes[node,4])) {
        V <- sim_rsa(V,nodes[node,4],nodes,edges,W,sigma_r,sigma_s,sigma_a,mu,kappa,alpha)
      }
      return(V)
    }
    # The whole process begins with the root node
    V <- sim_rsa(V,node=(self$tree$root_idx),self$tree$nodes,self$tree$edges,self$W,self$param$sigma_r,self$param$sigma_s,self$param$sigma_a,self$param$mu,self$param$kappa,self$param$alpha)

    # setting V for dataTree object
    self$V = V
  },

  #' @description Initializes values for the V field of the dataTree object to max likelihood under simple univariate Brownian motion
  #' @details The V object contains information about the simulated or data-augmented values of continuous characters 'r', 's', and 'a' at each internal node.
  init_V = function() {
    # ensuring that tree and V exist and are of correct classes
    if (!("detailedTree" %in% class(self$tree))) {stop("Tree must be of class detailedTree",call.=FALSE)}
    if (!("data.frame" %in% class(self$V))) {stop("V must be of class data.frame",call.=FALSE)}

    # To do a postorder tree traversal to get maximum likelihood estimates under Brownian motion, we first sort the tree into postorder
    order <- ape::reorder.phylo(self$tree$tree,"postorder")
    # Looping through the edges in the phylogeny
    for (i in 1:nrow(order$edge)) {
      edge <- order$edge[i,]
      node <- edge[1]
      # check if we do not already have a value assigned to the parent of this edge
      if (!(node %in% self$V$node)) {
        # find children & edge lengths
        children <- self$tree$edges[which(self$tree$edges[,1]==node),2]
        lengths <- self$tree$edges[which(self$tree$edges[,1]==node),3]
        values_1 <- unlist(self$V[children[1],2:4])
        values_2 <- unlist(self$V[children[2],2:4])
        # assign weighted average to V
        weighted_average <- (values_1*lengths[1] + values_2*lengths[2]) / (lengths[1]+lengths[2])
        private$.V$node[node] <- node
        private$.V$r[node] <- weighted_average[[1]]
        private$.V$s[node] <- weighted_average[[2]]
        private$.V$a[node] <- weighted_average[[3]]
      }
    }
  },

  #' @description Initializes values for the XY field of the dataTree object to max likelihood under simple univariate Brownian motion
  #' @details The XY object contains information about the simulated or data-augmented values of continuous characters 'x' and 'y' at each internal node, including the root.
  init_XY = function() {
    # ensuring that tree and XY exist and are of correct classes
    if (!("detailedTree" %in% class(self$tree))) {stop("Tree must be of class detailedTree",call.=FALSE)}
    if (!("data.frame" %in% class(self$XY))) {stop("XY must be of class data.frame",call.=FALSE)}

    # To do a postorder tree traversal to get maximum likelihood estimates under Brownian motion, we first sort the tree into postorder
    order <- ape::reorder.phylo(self$tree$tree,"postorder")
    for (i in 1:nrow(order$edge)) {
      # Looping through the edges in the phylogeny
      edge <- order$edge[i,]
      node <- edge[1]
      # check if we do not already have a value assigned to the parent of this edge
      if (!(node %in% self$XY$node)) {
        # find children & edge lengths
        children <- self$tree$edges[which(self$tree$edges[,1]==node),2]
        lengths <- self$tree$edges[which(self$tree$edges[,1]==node),3]
        values_1 <- unlist(self$XY[children[1],2:3])
        values_2 <- unlist(self$XY[children[2],2:3])
        # assign weighted average to V
        weighted_average <- (values_1*lengths[1] + values_2*lengths[2]) / (lengths[1]+lengths[2])
        private$.XY$node[node] <- node
        private$.XY$x[node] <- weighted_average[[1]]
        private$.XY$y[node] <- weighted_average[[2]]
      }
    }
  },

  #' @description Samples values for the XY field of the dataTree object according to multivariate normal distribution using full phylogenetic variance-covariance matrix (including internal nodes), conditioned on observed tips and expected values
  #' @details The XY object contains information about the reconstructed values of continuous characters 'x' and 'y' at each internal node, including the root.
  reconstruct = function() {

    new_XY <- self$XY
    # beginning by getting the root values of x and y, then assigning those values to the remaining internal nodes in the phylogeny
    root <- self$XY[self$tree$n_tips+1,2:3]
    for (i in (self$tree$n_tips+1):(2*self$tree$n_tips-1)) {
      new_XY[i,2] <- root[1]
      new_XY[i,3] <- root[2]
    }

    # getting expected values at the tips given the root
    # This is done by traversing the tree and performing cladogenetic scenarios at the nodes
    expectations <- get_expectations(self$tree,self$param,self$XY,self$W,self$V)

    # subtracting out expected values so that process is centered at (0,0)
    # This will give us the Brownian component of the evolution of x and y, without the deterministic part from the cladogenetic scenarios
    x_data <- t(as.matrix(self$XY$x))[1:self$tree$n_tips]
    x_expected <- t(as.matrix(expectations[,1]))
    x <- x_data - x_expected
    y_data <- t(as.matrix(self$XY$y))[1:self$tree$n_tips]
    y_expected <- t(as.matrix(expectations[,2]))
    y <- y_data - y_expected

    # sampling from multivariate normal conditioned on root and tips -- this gives the Brownian part of the process
    # Note that we are using a root of (0,0) and the new tip values with the expected values REMOVED (the Brownian component ONLY)
    n <- self$tree$n_nodes
    # These are the nodes we want to reconstruct
    dependent.ind <- c((self$tree$n_tips+2):(2*self$tree$n_tips-1))
    # These are the nodes we condition on, or "know" (tips and root)
    given.ind <- c(1:(self$tree$n_tips+1))
    given_x <- c(x,0)
    given_y <- c(y,0)
    # This is our phylogenetic variance/covariance matrix based on branch lengths
    # Includes all nodes in the tree, not just tips
    full_matrix <- self$tree$full_matrix
    full_matrix[self$tree$n_tips+1,self$tree$n_tips+1] <- 1e-10 # Having any 0 values makes it impossible to sample using this method, so we set any 0 values very small
    # The variances and covariances depend on the phylogenetic variance/covariance matrix and the rates of evolution
    matrix_x <- full_matrix * (self$param$sigma_x)^2
    matrix_y <- full_matrix * (self$param$sigma_y)^2
    # This function samples the values we want (internal nodes) conditioned on the values we know (tips and root)
    x_brownian <- condMVNorm::rcmvnorm(n=1,mean=rep(0,n),sigma=matrix_x,dependent.ind=dependent.ind,given.ind=given.ind,X.given=given_x,check.sigma=FALSE)[1,]
    y_brownian <- condMVNorm::rcmvnorm(n=1,mean=rep(0,n),sigma=matrix_y,dependent.ind=dependent.ind,given.ind=given.ind,X.given=given_y,check.sigma=FALSE)[1,]
    # This is the full Brownian component of our reconstruction, now we just have to add our cladogenetic scenarios back in
    brownian <- cbind(node=c((self$tree$n_tips+2):(2*self$tree$n_tips-1)),x=x_brownian,y=y_brownian)

    # Looping over nodes
    # We are starting from the dataframe constructed previously, where all internal node values are equal to the root values
    # We want to get the cladogenetic scenario at each node, and then apply the appropriate transformation to all non-tip descendant nodes
    for (i in (self$tree$n_tips+1):(2*self$tree$n_tips-1)) {
      # Getting children
      children <- self$tree$edges[which(self$tree$edges[,1]==i),2]
      right_child <- children[1]
      left_child <- children[2]
      # Getting non-tip left and right descendants
      right_descendants <- c(right_child,which(self$tree$descendants[right_child,]==1))
      right_descendants <- right_descendants[right_descendants > self$tree$n_tips]
      left_descendants <- c(left_child,which(self$tree$descendants[left_child,]==1))
      left_descendants <- left_descendants[left_descendants > self$tree$n_tips]
      # Getting cladogenetic scenario
      w <- self$W[which(self$W$node==i),]
      v <- self$V[which(self$V$node==i),]
      # Getting appropriate transformations for each set of non-tip descendants
      # Loops through all possible non-tip descendants
      for (j in 1:self$tree$n_nodes) {
        if (j %in% left_descendants) {
          # Applying appropriate transformations due to the cladogenetic scenario being investigated
          new_XY[j,2] <- transform_node(char="x",daughter="left",d=w$d,m=w$m,cval=w$cval,hval=w$hval,x=new_XY[j,2],r=v$r,s=v$s,a=v$a,z=self$param$z)
          new_XY[j,3] <- transform_node(char="y",daughter="left",d=w$d,m=w$m,cval=w$cval,hval=w$hval,y=new_XY[j,3],r=v$r,s=v$s,a=v$a,z=self$param$z)
        }
        if (j %in% right_descendants) {
          new_XY[j,2] <- transform_node(char="x",daughter="right",d=w$d,m=w$m,cval=w$cval,hval=w$hval,x=new_XY[j,2],r=v$r,s=v$s,a=v$a,z=self$param$z)
          new_XY[j,3] <- transform_node(char="y",daughter="right",d=w$d,m=w$m,cval=w$cval,hval=w$hval,y=new_XY[j,3],r=v$r,s=v$s,a=v$a,z=self$param$z)
        }
      }
    }
    # Now we have the deterministic component of the process, with the root values appropriately transformed at each internal node
    # Adding Brownian component
    new_x <- new_XY[,2] + c(rep(0,self$tree$n_tips+1),brownian[,2])
    new_y <- new_XY[,3] + c(rep(0,self$tree$n_tips+1),brownian[,3])
    new_XY[,2] <- new_x
    new_XY[,3] <- new_y
    # Saving values for this reconstruction
    self$XY <- new_XY
  },

  #' @description Simulates values for the XY field of the dataTree object
  #' @details The XY object contains information about the simulated or data-augmented values of continuous characters 'x', and 'y' at each internal node. If the values are input (for inference) rather than simulated, the internal non-root nodes will have NA values.
  sim_XY = function() {

    # making sure tree, param, W, and V exist and are of proper classes
    if (!("detailedTree" %in% class(self$tree))) {stop("Tree must be of class detailedTree",call.=FALSE)}
    if (!("ellipseParam" %in% class(self$param))) {stop("Param must be of class ellipseParam",call.=FALSE)}
    if (!("data.frame" %in% class(self$W))) {stop("W must be of class data.frame",call.=FALSE)}
    if (!("data.frame" %in% class(self$V))) {stop("V must be of class data.frame",call.=FALSE)}

    # starting XY
    XY <- data.frame(matrix(NA,nrow=self$tree$n_nodes,ncol=3))
    colnames(XY) <- c("node","x","y")
    XY$node <- 1:self$tree$n_nodes
    # Making sure the root x,y values match the parameters in the ellipseParam object
    XY$x[self$tree$root_idx] <- self$param$root_x
    XY$y[self$tree$root_idx] <- self$param$root_y

    # Function for simulating descendant nodes from a parent node
    sim_iter_xy <- function(XY,node,edges,W,V,z,sigma_x,sigma_y) {
      # Checking if the node is a tip (and we don't need to do any more work on this node)
      n_tips <- length(V[,1])-length(W[,1])
      if (node %in% 1:n_tips) {
        return(XY)
      }
      # Getting cladogenetic scenario
      w <- W[which(W$node==node),]
      v <- V[which(V$node==node),]
      current_x <- XY[node,2]
      current_y <- XY[node,3]
      # Getting appropriate transformations for each child
      left_start_x <- transform_node(char="x",daughter="left",d=w$d,m=w$m,cval=w$cval,hval=w$hval,x=current_x,r=v$r,s=v$s,a=v$a,z=z)
      left_start_y <- transform_node(char="y",daughter="left",d=w$d,m=w$m,cval=w$cval,hval=w$hval,y=current_y,r=v$r,s=v$s,a=v$a,z=z)
      right_start_x <- transform_node(char="x",daughter="right",d=w$d,m=w$m,cval=w$cval,hval=w$hval,x=current_x,r=v$r,s=v$s,a=v$a,z=z)
      right_start_y <- transform_node(char="y",daughter="right",d=w$d,m=w$m,cval=w$cval,hval=w$hval,y=current_y,r=v$r,s=v$s,a=v$a,z=z)
      # Getting children and branch lengths
      children <- edges[which(edges[,1]==node),2]
      right_child <- children[1]
      left_child <- children[2]
      lengths <- edges[which(edges[,1]==node),3]
      right_length <- lengths[1]
      left_length <- lengths[2]
      # Simulating a Brownian motion for x and y over each branch (left and right)
      left_vals <- c(rnorm(1,left_start_x,(sigma_x*sqrt(left_length))),rnorm(1,left_start_y,(sigma_y*sqrt(left_length))))
      right_vals <- c(rnorm(1,right_start_x,(sigma_x*sqrt(right_length))),rnorm(1,right_start_y,(sigma_y*sqrt(right_length))))
      # Assigning new node values to descendant nodes (left and right)
      XY[left_child,2] <- left_vals[1]
      XY[left_child,3] <- left_vals[2]
      XY[right_child,2] <- right_vals[1]
      XY[right_child,3] <- right_vals[2]
      return(XY)
    }

    # function to iterate over nodes
    sim_xy <- function(XY,node,nodes,edges,W,V,z,sigma_x,sigma_y) {
      # A simulation is performed on this node and its two daughter branches
      XY <- sim_iter_xy(XY,node,edges,W,V,z,sigma_x,sigma_y)
      # If a left daughter exists, the function recurs on the left daughter
      if (!is.na(nodes[node,3])) {
        XY <- sim_xy(XY,nodes[node,3],nodes,edges,W,V,z,sigma_x,sigma_y)
      }
      # If a right daughter exists, the function recurs on the right daughter
      if (!is.na(nodes[node,4])) {
        XY <- sim_xy(XY,nodes[node,4],nodes,edges,W,V,z,sigma_x,sigma_y)
      }
      return(XY)
    }
    # The whole process begins with the root node
    XY <- sim_xy(XY,node=(self$tree$root_idx),self$tree$nodes,self$tree$edges,self$W,self$V,self$param$z,self$param$sigma_x,self$param$sigma_y)

    # setting XY for dataTree object
    self$XY = XY
  },

  #' @description Simulates values for the W, V, and XY fields of the dataTree object
  #' @details This method assumes that the dataTree object already contains a tree of class "detailedTree" in the tree field. If you want to simulate a tree instead, use the sim() method.
  #' @param param An object of type ellipseParam containing information about the model parameters associated with a simulation or MCMC iteration. If no value is provided by the user or by the dataTree object, parameters will be simulated according to default priors (see ?ellipseParam for information about default priors).
  sim_data = function(param=self$param) {
    # make sure a detailedTree is present
    if (!("detailedTree" %in% class(self$tree))) {stop("Tree must be of class detailedTree",call.=FALSE)}
    # get new param object if none is provided
    if (is.null(param)) {param = ellipseParam$new()}
    self$param = param
    # simulate 3 data types
    self$sim_W()
    self$sim_V()
    self$sim_XY()
    # possibly update tip data
    if (is.null(self$tip_data)) {
      tip_df <- data.frame(cbind(taxon=self$tree$tree$tip.label,
                                 x=self$XY$x[1:self$tree$n_tips],
                                 y=self$XY$y[1:self$tree$n_tips],
                                 r=self$V$r[1:self$tree$n_tips],
                                 s=self$V$s[1:self$tree$n_tips],
                                 a=self$V$a[1:self$tree$n_tips]))
      self$tip_data = tip_df
    }
  },

  #' @description Simulates a tree and values for the W, V, and XY fields of the dataTree object
  #' @details This method simulates a tree according to the sim_tree() function, then simulates data associated with that tree.
  #' @param param An object of type ellipseParam containing information about the model parameters associated with a simulation or MCMC iteration. If no value is provided by the user or by the dataTree object, parameters will be simulated according to default priors (see ?ellipseParam for information about default priors).
  sim = function(param=self$param) {
    # get new param object if none is provided
    if (is.null(param)) {param = ellipseParam$new()}
    self$param = param
    # simulate tree
    self$tree <- sim_tree()
    # simulate 3 data types
    self$sim_data(param)
  },

  #' @description This method saves each element of a dataTree in a specified directory with a specified prefix
  #' @details The output files are: <prefix>.tree.txt, <prefix>.W.tsv, <prefix>.V.tsv, <prefix>.XY.tsv, and <prefix>.param.tsv
  #' @param filepath An object of class "character" giving the filepath to the directory where you want the output to be saved. The default behavior will save files to the current directory.
  #' @param prefix An object of class "character" giving the prefix you want to assign to output files. The default prefix is "sim".
  save = function(filepath=".", prefix="sim") {
    # making sure filepath and sim are of class "character"
    if (!("character" %in% class(filepath))) {stop("Filepath must be of class character",call.=FALSE)}
    if (!("character" %in% class(prefix))) {stop("Prefix must be of class character",call.=FALSE)}
    # writing files
    # creating the directory to save files if it does not exist
    dir.create(file.path(filepath), showWarnings = FALSE, recursive=TRUE)
    # writing the tree to file
    ape::write.tree(self$tree$tree,file=paste(filepath,"/",prefix,".tree.txt",sep=""))
    # writing the data to files
    write.table(self$W,file=paste(filepath,"/",prefix,".W.tsv",sep=""),row.names=FALSE,quote=FALSE,sep="\t")
    write.table(self$V,file=paste(filepath,"/",prefix,".V.tsv",sep=""),row.names=FALSE,quote=FALSE,sep="\t")
    write.table(self$XY,file=paste(filepath,"/",prefix,".XY.tsv",sep=""),row.names=FALSE,quote=FALSE,sep="\t")
    write.table(self$tip_data,file=paste(filepath,"/",prefix,".tip_data.tsv",sep=""),row.names=FALSE,quote=FALSE,sep="\t")
    # writing the parameters to file
    param_str <- paste(
      "cval","\t",paste(self$param$cval,collapse="\t"),"\n",
      "hval","\t",paste(self$param$hval,collapse="\t"),"\n",
      "rho_m","\t",paste(self$param$rho_m,collapse="\t"),"\n",
      "rho_c","\t",paste(self$param$rho_c,collapse="\t"),"\n",
      "rho_h","\t",paste(self$param$rho_h,collapse="\t"),"\n",
      "sigma_x","\t",paste(self$param$sigma_x),"\n",
      "sigma_y","\t",paste(self$param$sigma_y),"\n",
      "sigma_r","\t",paste(self$param$sigma_r),"\n",
      "sigma_s","\t",paste(self$param$sigma_s),"\n",
      "sigma_a","\t",paste(self$param$sigma_a),"\n",
      "mu","\t",paste(self$param$mu),"\n",
      "kappa","\t",paste(self$param$kappa),"\n",
      "root_x","\t",paste(self$param$root_x),"\n",
      "root_y","\t",paste(self$param$root_y),"\n",
      "root_r","\t",paste(self$param$root_r),"\n",
      "root_s","\t",paste(self$param$root_s),"\n",
      "root_a","\t",paste(self$param$root_a),"\n",
      "alpha","\t",paste(self$param$alpha),"\n",
      "z","\t",paste(self$param$z),
      sep="")
    dir.create(file.path(filepath), showWarnings = FALSE, recursive=TRUE)
    write(param_str,file=paste(filepath,"/",prefix,".param.tsv",sep=""))
  }
))

# dataTree functions ------------------------------------------------------

#' @title make_dataTree function
#' @description Takes an object of class "detailedTree" or "phylo" and creates an object of class "dataTree". May take a long time to generate for larger trees.
#' @details The dataTree object is used during simulation to contain the simulated history of all cladogenetic events (W) and continuous values (V, XY) at all nodes (including tips). When initialized, the data fields W, V, and XY will be NULL, and the tree field will be an object of class "detailedTree".
#' @param tree An object of class "detailedTree", or an object of class "phylo" (typically from the "ape" package)
#' @export
make_dataTree = function(tree, param=NULL, tip_data=NULL) {
  # This function is a wrapper for dataTree$new(tree)
  data_tree = dataTree$new(tree, param=param, tip_data=tip_data)
  return(data_tree)
}

#' @title sim_data function
#' @description Takes an object of class "dataTree", "detailedTree", or "phylo" and simulates values for the W, V, and XY fields, returning an object of class "dataTree"
#' @details This method assumes that a tree of class "detailedTree" or "phylo is provided, or the dataTree object already contains a tree of class "detailedTree" in the tree field. If you want to simulate a tree instead, use the sim_dataTree() function.
#' @param tree An object of class "dataTree" containing a "detailedTree" in the tree field, an object of class "detailedTree", or an object of class "phylo" (typically from the "ape" package)
#' @param param An object of type ellipseParam containing information about the model parameters associated with a simulation or MCMC iteration. If no value is provided by the user or by the dataTree object, parameters will be simulated according to default priors (see ?ellipseParam for information about default priors).
#' @export
sim_data = function(tree, param=NULL) {
  # making any "detailedTree" or "phylo" objects into "dataTree" objects
  if ("detailedTree" %in% class(tree) | "phylo" %in% class(tree)) {tree = make_dataTree(tree)}
  # making sure tree is of class "dataTree"
  if (!("dataTree" %in% class(tree))) {stop("Tree must be of class dataTree, detailedTree, or phylo.",call.=FALSE)}
  # getting param from dataTree if no param is provided by the user
  if (is.null(param)) {param = tree$param}
  # get new param object if none is provided
  if (is.null(param)) (param = ellipseParam$new())
  tree$param = param
  # simulating data
  tree$sim_data()
  # returning dataTree with simulated data
  return(tree)
}

#' @title sim_dataTree function
#' @description Simulates a dataTree including the tree and all data fields
#' @details This method initializes a new dataTree and simulates using the sim() method
#' @param param An object of type ellipseParam containing information about the model parameters associated with a simulation or MCMC iteration. If no value is provided by the user, parameters will be simulated according to default priors (see ?ellipseParam for information about default priors).
#' @export
sim_dataTree = function(param=NULL) {
  # simulating tree
  tree = dataTree$new()
  # get new param object if none is provided
  if (is.null(param)) (param = ellipseParam$new())
  tree$param <- param
  # simulating tree and data
  tree$sim()
  # returning tree
  return(tree)
}

#' @title save_dataTree function
#' @description This functions saves each element of a dataTree in a speficied directory with a specified prefix
#' @details The output files are: <prefix>.tree.txt, <prefix>.W.tsv, <prefix>.V.tsv, <prefix>.XY.tsv, and <prefix>.param.tsv
#' @param tree An object of class "dataTree" to be saved
#' @param filepath An object of class "character" giving the filepath to the directory where you want the output to be saved. The default behavior will save files to the current directory.
#' @param prefix An object of class "character" giving the prefix you want to assign to output files. The default prefix is "sim".
#' @export
save_dataTree = function(tree, filepath=".", prefix="sim") {
  # making sure tree is of class "dataTree"
  if (!("dataTree" %in% class(tree))) {stop("tree must be of class dataTree",call.=FALSE)}
  # saving dataTree
  dir.create(file.path(filepath), showWarnings = FALSE, recursive=TRUE)
  tree$save(filepath,prefix)
}
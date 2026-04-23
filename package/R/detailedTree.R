# detailedTree object -----------------------------------------------------

#' @title detailedTree object
#' @description This class stores an object of class "phylo" and generates other data relevant to the phylogeny. May take a long time to generate for larger trees.
#' @details The "detailedTree" object stores a tree object of class "phylo", typically from the "ape" package, alongside other important information about its nodes, branches, and structure. This information is useful when implementing phylogenetic analyses, like the ellipse evolution analysis.
#' @field tree An object of class "phylo", typically from the "ape" package. This must be an ultrametric, binary tree.
#' @field n_tips An integer giving the number of tips in the tree
#' @field n_nodes An integer giving the number of nodes in the tree, including root and tips
#' @field root_idx An integer giving the index of the root node, typically 1 greater than the number of tips
#' @field edges An object of class "data.frame" listing all edges with start, end, and length
#' @field nodes An object of class "data.frame" listing all nodes with parent, right child, and left child
#' @field descendants An object of class "matrix" indicating whether a node (column) is a descendant of another node (row) using binary representation (0 = non-descendant, 1 = descendant)
#' @field tip_descendants An object of class "matrix" indicating whether a tip (column) is a descendant of another node (row) using binary representation (0 = non-descendant, 1 = descendant), where all tips are descendants of themselves
#' @field phy_matrix An object of class "matrix" giving the phylogenetic variance-covariance matrix indicated by the tree structure (not scaled by evolutionary rates)
#' @field inv_matrix An object of class "matrix" giving the inverse of the phylogenetic variance-covariance matrix indicated by the tree structure (not scaled by evolutionary rates)
#' @field full_matrix An object of class "matrix" giving the full phylogenetic variance-covariance matrix indicated by the tree structure, including internal nodes (not scaled by evolutionary rates)
detailedTree <- R6::R6Class("detailedTree",

private = list(

  # Initializing private fields (correspond to active bindings)
  # These fields (starting with '.') are for accessing information internally
  # Active bindings with corresponding names are for accessing information externally
  .tree = NULL,
  .n_tips = NULL,
  .n_nodes = NULL,
  .root_idx = NULL,
  .edges = NULL,
  .nodes = NULL,
  .descendants = NULL,
  .tip_descendants = NULL,
  .phy_matrix = NULL,
  .inv_matrix = NULL,
  .full_matrix = NULL,

  # Assigns a node and its descendants to the parent node's list of descendants
  # Information is read from the "nodes" field, which contains parents and children for all nodes
  # Operates on the descendants matrix
  get_descendants = function(nodes,descendants,node) {
    # Checking if there is a parent to operate on -- if not, do nothing
    if (!is.na(nodes[node,2])) {
      # Get the index of the parent
      parent <- nodes[node,2]
      # Get the descendants of the target node (not the parent)
      descendants_list <- which(descendants[node,]==1)
      # Assign the target node as a descendant of the parent
      descendants[parent,node] <- 1
      # Assign each descendant of the target node as a descendant of the parent
      for (i in 1:self$n_nodes) {
        if (i %in% descendants_list) {
          descendants[parent,i] <- 1
        }
      }
    }
    return(descendants)
  },

  # Uses recursion to traverse the tree tips-to-root and generates a pairwise matrix
  # Indicating whether each node (column) is a descendant of another node (row)
  # 1 = descendant, 0 = not descendant
  get_all_descendants = function(nodes,descendants,node) {
    # Assign a node and all of its descendants as descendants of the parent node (if there is no parent, this will do nothing)
    descendants <- private$get_descendants(nodes,descendants,node)
    # If there is a parent node, recur on that parent node
    if (!is.na(nodes[node,2])) {
      descendants <- private$get_all_descendants(nodes,descendants,nodes[node,2])
    }
    return(descendants)
  },

  # Generates the phylogenetic variance/covariance matrix
  # Based on the tree structure and branch lengths
  # THIS IS SAVED --> only need to calculate & invert once
  get_phy_matrix = function(matrix,edges,tip_descendants) {
    # Getting the number of edges in the tree
    n_edges <- length(edges[,1])
    # Looping over all the edges in the tree
    for (edge in 1:n_edges) {
      # Finding the start node of the branch
      node <- edges[edge,2]
      # Finding the length of the branch
      length <- edges[edge,3]
      # Finding which tips share this edge in their evolutionary history
      tip_list <- which(tip_descendants[node,]==1)
      # Looping through each pair of tips and adding the shared evolutionary history from this branch (branch length)
      # Diagonals are included (a tip shares its whole history with itself) to generate variances on the diagonal of the VCM
      for (i in tip_list) {
        for (j in tip_list) {
          matrix[i,j] <- matrix[i,j] + length
        }
      }
    }
    return(matrix)
  },

  # Generates the phylogenetic variance/covariance matrix including internal nodes
  # Based on the tree structure and branch lengths
  # THIS IS SAVED --> only need to calculate & invert once
  get_full_matrix = function(matrix,edges,all_descendants) {
    # Similar to get_phy_matrix, looping over edges
    n_edges <- length(edges[,1])
    for (edge in 1:n_edges) {
      # Getting the starting node of the edge
      node <- edges[edge,2]
      # Getting the length of the edge
      length <- edges[edge,3]
      # Getting all descendants that share this branch, not just tips
      descendant_list <- c(node,which(all_descendants[node,]==1))
      # Looping through each pair of nodes and adding the shared evolutionary history from this branch (branch length)
      # Diagonals are included (a node shares its whole history with itself)
      for (i in descendant_list) {
        for (j in descendant_list) {
          matrix[i,j] <- matrix[i,j] + length
        }
      }
    }
    return(matrix)
  }
),

active = list(

  # Active binding for the "tree" field of the detailedTree object
  # This tree can be set using detailedTree$tree <- input_tree
  # This tree can be accessed using detailedTree$tree
  # When a new tree is set, will check for type and update other informational fields
  tree = function(value) {

    # When called with no value, returns the tree itself (of class "phylo"), assuming user is trying to get the current value
    # If there is a value, the function will assume the user is trying to set that value
    if (missing(value)) {private$.tree} else {

      # phylogenetic tree
      # Making sure the user is setting this field using a tree of class 'phylo' that is ultrametric and binary
      if (!(class(value)=="phylo")) {stop("Tree should be of class phylo",call.=FALSE)}
      if (!(ape::is.ultrametric(value))) {stop("Tree should be ultrametric",call.=FALSE)}
      if (!(ape::is.binary(value))) {stop("Tree should be binary",call.=FALSE)}
      private$.tree <- value

      # Checking tree size and length
      if (length(self$tree$tip.label) < 3) {stop("Tree must have at least 3 tips.",call.=FALSE)}
      if (length(self$tree$tip.label) > 1000) {warning("Tree is very large; this may take a while.",immediate.=TRUE)}
      if (max(ape::branching.times(self$tree)) < 1) {warning("Tree is very young (a tree height >=1 is expected). Behavior may be unusual; consider rescaling tree.",immediate.=TRUE)}
      if (max(ape::branching.times(self$tree)) > 500) {warning("Tree is very old (a tree height <=500 is expected). Behavior may be unusual; consider rescaling tree.",immediate.=TRUE)}

      # Performing a series of calculations on the new tree (these cannot be set independently)
      # number of tips
      private$.n_tips <- length(self$tree$tip.label)

      # number of nodes
      private$.n_nodes <- self$tree$Nnode + self$n_tips

      # index of root
      private$.root_idx <- self$n_tips + 1

      # dataframe of edges with lengths
      df <- self$tree$edge
      df <- cbind(df,self$tree$edge.length)
      colnames(df) <- c("start","end","length")
      private$.edges <- df

      # dataframe of nodes with parents and children
      df <- c()
      # Looping through all nodes
      for (node in 1:self$n_nodes) {
        nprl <- c(node,NA,NA,NA) # node, parent, right daughter, left daughter
        # Checking if the node has a parent, and setting the parent if so
        has_parent <- node %in% self$edges[,2]
        if (has_parent) {nprl[2]=self$edges[which(self$edges[,2]==node),1]}
        # Checking if the node has children, and setting the children if so (there will always be 0 or 2)
        has_children <- !(node %in% seq(1:self$n_tips))
        if (has_children) {nprl[3]=self$edges[which(self$edges[,1]==node),2][1]}
        if (has_children) {nprl[4]=self$edges[which(self$edges[,1]==node),2][2]}
        df <- c(df,nprl)
      }
      df <- data.frame(matrix(df,ncol=4,byrow=TRUE))
      colnames(df) <- c("node","parent","right_child","left_child")
      private$.nodes <- df

      # matrix of all node descendants
      df <- matrix(0,nrow=self$n_nodes,ncol=self$n_nodes)
      # Using get_all_descendants on every tip to assign all descendants to parent nodes
      # This recursive algorithm must be performed on each tip to make sure those tips (and all upstream sections of the tree) are properly assigned to their parents
      for (i in 1:self$n_tips) {df <- private$get_all_descendants(self$nodes,df,i)}
      private$.descendants <- df

      # matrix of tip descendants
      # this is a subset of the matrix with all descendants, including only the tips
      # each tip is also set to be a descendant of itself, which is important for performing operations on the terminal branches
      # (we use tip descendants of left and right daughters to set expected values after cladogenetic events)
      df <- self$descendants[,1:self$n_tips]
      for (i in 1:self$n_tips) {
        for (j in 1:self$n_tips) {
          if (i==j) {
            df[i,j] <- 1
          }
        }
      }
      private$.tip_descendants <- df

      # phylogenetic variance-covariance matrix
      df <- matrix(0,nrow=self$n_tips,ncol=self$n_tips)
      private$.phy_matrix <- private$get_phy_matrix(df,self$edges,self$tip_descendants)

      # inverse of phylogenetic variance-covariance matrix
      private$.inv_matrix <- solve(self$phy_matrix)

      # full phylogenetic variance-covariance matrix (including internal nodes)
      df <- matrix(0,nrow=self$n_nodes,ncol=self$n_nodes)
      private$.full_matrix <- private$get_full_matrix(df,self$edges,self$descendants)
    }
  },

  # Active bindings for accessing informational fields about the provided tree
  # Will return associated values when called with detailedTree${field name}, but cannot be directly set by the user (will be based on the tree)
  n_tips = function(value) {if (missing(value)) {private$.n_tips} else {stop("n_tips is read only",call.=FALSE)}},
  n_nodes = function(value) {if (missing(value)) {private$.n_nodes} else {stop("n_nodes is read only",call.=FALSE)}},
  root_idx = function(value) {if (missing(value)) {private$.root_idx} else {stop("root_idx is read only",call.=FALSE)}},
  edges = function(value) {if (missing(value)) {private$.edges} else {stop("edges is read only",call.=FALSE)}},
  nodes = function(value) {if (missing(value)) {private$.nodes} else {stop("nodes is read only",call.=FALSE)}},
  descendants = function(value) {if (missing(value)) {private$.descendants} else {stop("descendants is read only",call.=FALSE)}},
  tip_descendants = function(value) {if (missing(value)) {private$.tip_descendants} else {stop("tip_descendants is read only",call.=FALSE)}},
  phy_matrix = function(value) {if (missing(value)) {private$.phy_matrix} else {stop("matrix is read only",call.=FALSE)}},
  inv_matrix = function(value) {if (missing(value)) {private$.inv_matrix} else {stop("inv_matrix is read only",call.=FALSE)}},
  full_matrix = function(value) {if (missing(value)) {private$.full_matrix} else {stop("full_matrix is read only",call.=FALSE)}}
),

public = list(

  #' @description Creates a new instance of the detailedTree class
  #' @param tree An object of class "phylo", typically from the "ape" package
  initialize = function(tree) {
    self$tree <- tree
  }

))

# detailedTree functions --------------------------------------------------

#' @title sim_tree() function
#' @description Creates a new object of the "detailedTree" class using a simulated tree. May take a long time to generate for larger trees.
#' @details The simulated tree has height 1 and is generated using the package "phytools" through a pure-birth process with a rate of 4. Trees with less than 25 tips or more than 250 tips are rejected.
#' @export
sim_tree <- function() {
  # ensuring a valid tree is returned
  valid <- FALSE
  while (!valid) {
    # simulating a tree
    tree <- phytools::pbtree(b=4,d=0,t=1,type="continuous")
    # checking the number of tips
    n_tips <- length(tree$tip.label)
    if (n_tips >= 25 & n_tips <= 250) {valid <- TRUE}
  }
  # returning tree
  return(detailedTree$new(tree))
}

#' @title make_detailedTree() function
#' @description Turns a tree of class "phylo" into a "detailedTree" object -- this function is a wrapper for detailedTree$new(). May take a long time to generate for larger trees.
#' @details The "detailedTree" object stores a tree object of class "phylo", typically from the "ape" package, alongside other important information about its nodes, branches, and structure. This information is useful when implementing phylogenetic analyses, like the ellipse evolution analysis.
#' @param tree An object of class "phylo". This must be an ultrametric, binary tree.
#' @export
make_detailedTree = function(tree) {
  # This is just a wrapper for detailedTree$new(tree)
  detailed_tree = detailedTree$new(tree)
  return(detailed_tree)
}
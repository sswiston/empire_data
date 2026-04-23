#' @title plot_ellipse() function
#' @description Accepts parameters describing an ellipse and returns a simple plot
#' @param x The x coordinate of the ellipse's center point
#' @param y The y coordinate of the ellipse's center point
#' @param r The r value of the ellipse, which helps define its oblongness (x coordinate of the tilt vector)
#' @param s The s value of the ellipse, which helps define its oblongness (y coordinate of the tilt vector)
#' @param a The log area of the ellipse
#' @param z The height of the tilt vector, default value is 10
#' @export
plot_ellipse <- function(x,y,r,s,a,z=10) {

  # Making plot coordinates
  ellipse <- make_ellipse_coords(x,y,r,s,a,z)

  # Plotting the ellipse
  ellipse_plot <- plot(ellipse$x,ellipse$y,xlab="",ylab="",asp=1,type="l")
  return(ellipse_plot)
}

#' @title plot_scenario() function
#' @description Accepts parameters describing an ellipse and cladogenetic scenario and returns a simple plot
#' @param x The x coordinate of the ancestral ellipse's center point
#' @param y The y coordinate of the ancestral ellipse's center point
#' @param r The r value of the ancestral ellipse, which helps define its oblongness (x coordinate of the tilt vector)
#' @param s The s value of the ancestral ellipse, which helps define its oblongness (y coordinate of the tilt vector)
#' @param a The log area of the ancestral ellipse
#' @param d The daughter configuration during cladogenesis (0,1), where 0 makes the left daughter D1 and 1 makes the right daughter D1
#' @param m The cladogenetic mode (0,1), where 0 is budding and 1 is splitting
#' @param cval The c value (not index) of the cladogenetic scenario's concentric circle --> ex. 0.5
#' @param hval The h value (not index) of the cladogenetic scenario's direction line --> ex. pi/4, if required
#' @param z The height of the tilt vector, default value is 10
#' @param alpha The alpha parameter, default value is -5
#' @export
plot_scenario <- function(x,y,r,s,a,d,m,cval,hval,z=10,alpha=-5) {

  # Getting daughter info
  # We need x, y, and a for both the left and right daughters after cladogenesis
  # Based on the cladogenetic event input values
  x_left  = transform_node(char="x",daughter="left", d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  y_left  = transform_node(char="y",daughter="left", d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  a_left  = transform_node(char="a",daughter="left", d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  x_right = transform_node(char="x",daughter="right",d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  y_right = transform_node(char="y",daughter="right",d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  a_right = transform_node(char="a",daughter="right",d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)

  # Making plot
  # First, generating ellipse points for ancestor and both daughters
  ancestor <- make_ellipse_coords(x,y,r,s,a,z)
  left_daughter <- make_ellipse_coords(x_left,y_left,r,s,a_left,z)
  right_daughter <- make_ellipse_coords(x_right,y_right,r,s,a_right,z)

  # Determining graph boundaries based on the ellipses being plotted
  xlow <- min(min(ancestor$x),min(left_daughter$x),min(right_daughter$x))
  xhigh <- max(max(ancestor$x),max(left_daughter$x),max(right_daughter$x))
  ylow <- min(min(ancestor$y),min(left_daughter$y),min(right_daughter$y))
  yhigh <- max(max(ancestor$y),max(left_daughter$y),max(right_daughter$y))

  # Plotting all 3 ellipses, with the ancestor in black, the left daughter in blue, and the right daughter in red
  plot(ancestor$x,ancestor$y,xlab="",ylab="",asp=1,type="l",xlim=c(xlow,xhigh),ylim=c(ylow,yhigh))
  lines(left_daughter$x,left_daughter$y,col="blue")
  lines(right_daughter$x,right_daughter$y,col="red")
  ellipse_plot <- recordPlot()
  return(ellipse_plot)
}

#' @title plot_scenario_annotated() function
#' @description Accepts parameters describing an ellipse and cladogenetic scenario and returns an annotated plot using ggplot2
#' @param x The x coordinate of the ancestral ellipse's center point
#' @param y The y coordinate of the ancestral ellipse's center point
#' @param r The r value of the ancestral ellipse, which helps define its oblongness (x coordinate of the tilt vector)
#' @param s The s value of the ancestral ellipse, which helps define its oblongness (y coordinate of the tilt vector)
#' @param a The log area of the ancestral ellipse
#' @param d The daughter configuration during cladogenesis (0,1), where 0 makes the left daughter D1 and 1 makes the right daughter D1
#' @param m The cladogenetic mode (0,1), where 0 is budding and 1 is splitting
#' @param cval The c value (not index) of the cladogenetic scenario's concentric circle --> ex. 0.5
#' @param hval The h value (not index) of the cladogenetic scenario's direction line --> ex. pi/4, if required
#' @param z The height of the tilt vector, default value is 10
#' @param alpha The alpha parameter, default value is -5
#' @export
plot_scenario_annotated <- function(x,y,r,s,a,d,m,cval,hval,z=10,alpha=-5) {

  # Getting daughter info -- see plot_scenario() for details
  x_left  = transform_node(char="x",daughter="left", d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  y_left  = transform_node(char="y",daughter="left", d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  a_left  = transform_node(char="a",daughter="left", d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  x_right = transform_node(char="x",daughter="right",d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  y_right = transform_node(char="y",daughter="right",d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)
  a_right = transform_node(char="a",daughter="right",d=d,m=m,cval=cval,hval=hval,x=x,y=y,r=r,s=s,a=a,alpha=alpha,z=z)

  # Making plot
  ancestor <- make_ellipse_coords(x,y,r,s,a,z)
  left_daughter <- make_ellipse_coords(x_left,y_left,r,s,a_left,z)
  right_daughter <- make_ellipse_coords(x_right,y_right,r,s,a_right,z)
  # We also want an ellipse to represent the concentric circle used for the cladogenetic scenario
  c_ellipse <- make_ellipse_coords(x,y,r,s,a+log(cval^2),z)

  # Determining graph boundaries based on the ellipses and the tilt point
  xlow <- min(min(ancestor$x),min(left_daughter$x),min(right_daughter$x),min(c_ellipse$x),r)
  xhigh <- max(max(ancestor$x),max(left_daughter$x),max(right_daughter$x),max(c_ellipse$x),r)
  ylow <- min(min(ancestor$y),min(left_daughter$y),min(right_daughter$y),min(c_ellipse$y),s)
  yhigh <- max(max(ancestor$y),max(left_daughter$y),max(right_daughter$y),max(c_ellipse$y),s)

  # Expanding the graph slightly, so annotations fit
  x_low_lim <- xlow - (xhigh - xlow) * 0.1
  x_high_lim <- xhigh + (xhigh - xlow) * 0.1
  y_low_lim <- ylow - (yhigh - ylow) * 0.1
  y_high_lim <- yhigh + (yhigh - ylow) * 0.1

  # d is the daughter configuration during cladogenesis (0,1), where 0 makes the left daughter D1 and 1 makes the right daughter D1
  # m is the cladogenetic mode (0,1), where 0 is budding and 1 is splititng
  # c is the c value (not index) of the cladogenetic scenario's concentric circle --> ex. 0.5
  # h is the h value (not index) of the cladogenetic scenario's direction line --> ex. pi/4, if required

  # Gathering information for annotations
  if (d == 0) {left <- "d1"} else {left <- "d2"}
  if (d == 0) {right <- "d2"} else {right <- "d1"}
  if (m == 0) {mode <- paste0("budding, alpha=",alpha)} else {mode <- "splitting"}

  # Annotation describing mode and daughter configuration
  top_text <- paste0(
    "Left Daughter: ",left," (blue)","\n",
    "Right Daughter: ",right," (red)","\n",
    "Mode: ",mode
  )

  # Annotation describing concentric circle, direction line, and log area of ancestor
  bottom_text <- paste0(
    "C*: ",round(cval,2),"\n",
    "H*: ",round(hval/pi,2)," x pi","\n",
    "Log Area (a): ",round(a,2)
  )

  # Annotations for center point and tilt point, adjusted to fit in graph
  center_label <- paste0("Center: (",round(x,2),",",round(y,2),")")
  tilt_label <- paste0("Tilt Point: (",round(r,2),",",round(s,2),",",round(z,2),")")
  if (s > y) {
    vjust_center <- 1
    vjust_tilt <- 0}
  else {
    vjust_center <- 0
    vjust_tilt <- 1}
  # Direction line
  slope <- (y_right-y_left)/(x_right-x_left)
  intercept <- y_left - (x_left * slope)

  # Combining ellipse plots and annotations
  ellipse_plot <- ggplot2::ggplot() +
    ggplot2::lims(x=c(x_low_lim,x_high_lim),y=c(y_low_lim,y_high_lim)) +
    ggplot2::coord_fixed() +
    # Ancestor ellipse
    #ggplot2::geom_polygon(data=ancestor,ggplot2::aes(x=x,y=y),color="black",fill=NA,linewidth=2) +
    ggplot2::geom_polygon(data=ancestor,ggplot2::aes(x=x,y=y),fill="purple",alpha=0.5,color="purple",lwd=2) +
    # Direction line
    ggplot2::geom_abline(slope=slope,intercept=intercept,linetype="dashed") +
    # Left daughter ellipse
    #ggplot2::geom_polygon(data=left_daughter,ggplot2::aes(x=x,y=y),color="blue",fill=NA) +
    ggplot2::geom_polygon(data=left_daughter,ggplot2::aes(x=x,y=y),fill=NA,color="blue",lwd=2,lty="11") +
    ggplot2::geom_polygon(data=left_daughter,ggplot2::aes(x=x,y=y),fill=NA,color="blue",lwd=.75) +
    # Right daughter ellipse
    #ggplot2::geom_polygon(data=right_daughter,ggplot2::aes(x=x,y=y),color="red",fill=NA) +
    ggplot2::geom_polygon(data=right_daughter,ggplot2::aes(x=x,y=y),fill=NA,color="red",lwd=2,lty="11") +
    ggplot2::geom_polygon(data=right_daughter,ggplot2::aes(x=x,y=y),fill=NA,color="red",lwd=.75) +
    # Concentric circle
    ggplot2::geom_polygon(data=c_ellipse,ggplot2::aes(x=x,y=y),color="black",fill=NA,linetype="dashed") +
    # Tilt vector
    ggplot2::geom_segment(ggplot2::aes(x=x,xend=r,y=y,yend=s),linetype="dashed",color="black") +
    # Left daughter centroid
    ggplot2::geom_point(ggplot2::aes(x=x_left,y=y_left),color="blue") +
    # Right daughter centroid
    ggplot2::geom_point(ggplot2::aes(x=x_right,y=y_right),color="red") +
    # Labels
    ggplot2::annotate("label",x=x_low_lim,y=y_high_lim,label=top_text,hjust="inward",vjust="inward") +
    ggplot2::annotate("label",x=x_high_lim,y=y_low_lim,label=bottom_text,hjust="inward",vjust="inward") +
    ggplot2::annotate("label",x=x,y=y,label=center_label,vjust=vjust_center,hjust=0.5) +
    ggplot2::annotate("label",x=r,y=s,label=tilt_label,vjust=vjust_tilt,hjust=0.5) +
    # Ancestor centroid
    ggplot2::geom_point(ggplot2::aes(x=x,y=y),color="purple") +
    # Tilt point
    ggplot2::geom_point(ggplot2::aes(x=r,y=s)) +
    ggplot2::theme_bw()
  return(ellipse_plot)
}

#' @title plot_label_tree() function
#' @description Plots a tree including taxon names and node numbers for easy reference. The function returns a ggplot object, but also saves the tree to a file, which is more legible.
#' @param tree The object of class phylo (ape package), detailedTree, or dataTree to plot
#' @param type The type of tree to plot. Options are "dated" or "cladogram", default is "dated"
#' @param filepath An object of class "character" giving the filepath to the directory where you want the output to be saved. The default behavior will save files to the current directory.
#' @export
plot_label_tree <- function(tree, type="dated", filepath=".") {
  # Type checking
  if ("phylo" %in% class(tree)) {phy <- tree}
  else if ("detailedTree" %in% class(tree)) {phy <- tree$tree}
  else if ("dataTree" %in% class(tree)) {phy <- tree$tree$tree}
  else {stop("Tree must be of class phylo, detailedTree, or dataTree")}
  # Setting up branch lengths for either a cladogram or dated phylogeny
  # And determining the appropriate size for the plot
  branch.length=NA
  if (type=="cladogram") {
    branch.length="none"
    lim <- max(ape::node.depth(phy,method=2))*1.2
  }
  if (type=="dated") {
    branch.length=1
    lim <- max(ape::branching.times(phy))*1.2
  }
  # Getting tip labels and including node numbers
  labels <- phy$tip.label
  numbers <- c(1:length(labels))
  for (i in 1:length(labels)) {
    labels[i] <- paste0("[",numbers[i],"] ",labels[i])
  }
  phy$tip.label <- labels
  # Plotting the tree
  tree_plot <- ggtree::ggtree(phy,branch.length=branch.length) +
    # Including tip labels
    ggtree::geom_tiplab() +
    # Including node numbers
    ggtree::geom_nodelab(geom="label",ggtree::aes(label=node),angle=0,vjust=0.5) +
    ggplot2::xlim(0,lim) +
    ggplot2::scale_y_continuous()
  # Saving tree file
  ggplot2::ggsave(paste0(filepath,"/tree.pdf"),tree_plot,dpi=600,height=length(phy$tip.label)/4,width=log(length(phy$tip.label))*3,limitsize=FALSE)
  return(tree_plot)
}

#' @title plot_grid() function
#' @description Plots the grid used for discrete cladogenetic scenarios, returning a ggplot object.
#' @param model An object of class ellipseParam, dataTree, or ellipseMCMC with a grid system to plot.
#' @export
plot_grid <- function(model) {
  # Type checking
  if ("ellipseParam" %in% class(model)) {params <- model}
  else if ("dataTree" %in% class(model)) {params <- model$param}
  else if ("ellipseMCMC" %in% class(model)) {params <- model$dataTree$param}
  else {stop("Model must be of class ellipseParam, dataTree, or ellipseMCMC")}

  # Getting the values associated with the c and h indices
  cval <- params$cval
  hval <- params$hval

  # Getting the c and h indices
  c_ind <- c(0:(length(cval)-1))
  h_ind <- c(0:(length(hval)-1))

  # Creating concentric circles for the resulting grid
  # And labeling those circles with labels located at c_lab_points
  circles <- data.frame(matrix(NA,nrow=0,ncol=3))
  c_lab_points <- data.frame(matrix(NA,nrow=0,ncol=3))
  for (i in 1:length(c_ind)) {
    a <- log(pi * cval[i]^2)
    circle <- make_ellipse_coords(0,0,0,0,a)
    circle$c <- c_ind[i]
    circles <- rbind(circles,circle)
    row <- c(0,-1*cval[i],paste0("c",c_ind[i],"=",round(cval[i],2)))
    c_lab_points <- rbind(c_lab_points,row)
  }

  circles$c <- as.factor(circles$c)
  colnames(c_lab_points) <- c("x","y","c")
  c_lab_points$x <- as.numeric(c_lab_points$x)
  c_lab_points$y <- as.numeric(c_lab_points$y)

  # Creating the ancestral circle for the resulting grid
  ancestor <- make_ellipse_coords(0,0,0,0,log(pi))

  # Creating direction lines for the resulting grid
  # And labeling those direction lines with labesl located at h_lab_points
  outer_radius <- max(cval) * 1.2
  h_lab_points <- data.frame(matrix(NA,nrow=0,ncol=3))
  for (i in 1:length(h_ind)) {
    row <- c(outer_radius * cos(hval[i]), outer_radius * sin(hval[i]), paste0("h",h_ind[i],"=",round(hval[i],2)))
    h_lab_points <- rbind(h_lab_points,row)
  }

  colnames(h_lab_points) <- c("x","y","h")
  h_lab_points$x <- as.numeric(h_lab_points$x)
  h_lab_points$y <- as.numeric(h_lab_points$y)

  # Including information about daughter configuration and mode
  text <- paste0(
    "d0: left daughter D1\n",
    "d1: right daughter D1\n",
    "m0: budding\n",
    "m1: splitting"
  )

  # Creating grid plot
  grid_plot <- ggplot2::ggplot(circles, ggplot2::aes(x=x,y=y,group=c)) +
    # Concentric circles
    ggplot2::geom_polygon(fill=NA,color="lightgray") +
    # Ancestral circle
    ggplot2::geom_polygon(data=ancestor,ggplot2::aes(x=x,y=y,group=NULL),fill=NA,color="lightgray",lwd=2.5) +
    # Direction lines
    ggplot2::geom_segment(data=h_lab_points,ggplot2::aes(x=0,y=0,xend=x,yend=y,group=NULL),color="lightgray") +
    # Concentric circle labesl
    ggplot2::geom_text(data=c_lab_points,ggplot2::aes(x=x,y=y,label=c,group=NULL),size=3,family="mono") +
    # Direction line labels
    ggplot2::geom_text(data=h_lab_points,ggplot2::aes(x=x,y=y,label=h,group=NULL),size=3,family="mono") +
    ggplot2::lims(x=c(-1.4*max(cval),1.4*max(cval)),y=c(-1.4*max(cval),1.4*max(cval))) +
    # Information about d/m
    ggplot2::geom_text(ggplot2::aes(x=-1.4*max(cval),y=1.4*max(cval),group=NULL),label=text,hjust="inward",vjust="inward",size=2.5,family="mono") +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio=1,
                   panel.grid=ggplot2::element_blank(),
                   axis.text=ggplot2::element_blank(),
                   axis.ticks=ggplot2::element_blank(),
                   axis.title=ggplot2::element_blank())

  return(grid_plot)
}

#' @title plot_tree_map() function
#' @description Plots the phylogenetic tree on a map, with nodes located at their reconstructed centroids
#' @param tree a **dataTree** object containing your phylogeny and tip data
#' @param reconstruction A **data.frame** of ancestral ellipse reconstructions generated by `reconstruct_ellipses()`
#' @param map An **sf** polygon to plot in the background (optional)
#' @param label Whether to include node and tip numbers, default = FALSE
#' @export
plot_tree_map <- function(tree, reconstruction, map=NULL, label=FALSE) {

  # Type checking
  if (!("dataTree" %in% class(tree))) {stop("Tree must be of class dataTree")}

  node_data <- reconstruction
  ellipse_data <- tree$tip_data
  tree_data <- data.frame(tree$tree$tree$edge)
  colnames(tree_data) <- c("N1", "N2")

  node_locations <- c()

  # Tips
  for (i in 1:length(tree$tree$tree$tip.label)) {
    node_locations$node[i] <- i
    taxon <- tree$tree$tree$tip.label[i]
    node_locations$x[i] <- ellipse_data$x[which(ellipse_data$taxon==tree$tree$tree$tip.label[i])]
    node_locations$y[i] <- ellipse_data$y[which(ellipse_data$taxon==tree$tree$tree$tip.label[i])]
  }

  # Non-tips
  node_locations$node <- c(node_locations$node, node_data$node)
  node_locations$x <- c(node_locations$x, node_data$x)
  node_locations$y <- c(node_locations$y, node_data$y)

  node_locations <- data.frame(node_locations)

  # Categorizing node types
  node_locations$highlight <- rep("none",nrow(node_locations))
  node_locations$highlight[tree$tree$n_tips + 1] <- "root"
  node_locations$highlight[1:tree$tree$n_tips] <- "tip"

  edge_points <- c()

  for (i in 1:nrow(tree_data)) {
    node_1 <- tree_data$N1[i]
    edge_points$X1[i] <- node_locations$x[which(node_locations$node==node_1)]
    edge_points$Y1[i] <- node_locations$y[which(node_locations$node==node_1)]
    node_2 <- tree_data$N2[i]
    edge_points$X2[i] <- node_locations$x[which(node_locations$node==node_2)]
    edge_points$Y2[i] <- node_locations$y[which(node_locations$node==node_2)]
  }

  edge_points <- data.frame(edge_points)

  xmin <- min(node_locations$x)
  xmax <- max(node_locations$x)
  ymin <- min(node_locations$y)
  ymax <- max(node_locations$y)
  xscale <- (xmax - xmin) * 0.1
  yscale <- (ymax - ymin) * 0.1

  if (is.null(map)) {
    coords <- matrix(c(xmin - xscale, ymin - yscale,
                       xmax + xscale, ymin - yscale,
                       xmax + xscale, ymax + yscale,
                       xmin - xscale, ymax + yscale,
                       xmin - xscale, ymin - yscale),
                     ncol=2,
                     byrow=TRUE)
    map = sf::st_as_sf(sf::st_sfc(sf::st_polygon(list(cbind(coords[,1],coords[,2])))))
  }

  node_tree <- ggplot2::ggplot() +
    ggplot2::geom_sf(data=map, fill=NA, color="black") +
    ggplot2::geom_segment(data=edge_points,ggplot2::aes(x=X1,y=Y1,xend=X2,yend=Y2),lty="dotted",lwd=.25) +
    ggplot2::geom_point(data=node_locations,ggplot2::aes(x=x,y=y,pch=highlight),fill="white") +
    {if(label)ggplot2::geom_text(data=node_locations,ggplot2::aes(x=x,y=y,label=node),vjust=1.2,position="nudge")}+
    ggplot2::scale_shape_manual(values=c(20,13,21),labels=c("Internal Node", "Root", "Tip")) +
    ggplot2::labs(pch="Node Type") +
    ggplot2::guides(shape=ggplot2::guide_legend(override.aes=list(color="black",fill=NA),keyheight=0.1)) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill="white",color=NA),
                   legend.position="right",
                   legend.title = ggplot2::element_text(color="black",size=12),
                   legend.text = ggplot2::element_text(color="black",size=8),
                   legend.box.background = ggplot2::element_rect(fill="white", color=NA),
                   legend.margin = ggplot2::margin(4,4,4,4),
                   legend.spacing.y = ggplot2::unit(.1,"cm"))

  return(node_tree)
}

#' @title plot_lineage() function
#' @description Plots a single evolving lineage on a map, starting with the root of the clade and ending in a particular tip, including reconstructed ancestral ellipses and cladogenetic scenarios.
#' @param tree a **dataTree** object containing your phylogeny and tip data
#' @param reconstruction A **data.frame** of ancestral ellipse reconstructions generated by `reconstruct_ellipses()`
#' @param tip A **character** (string) taxon name or **numeric** tip number of the lineage to plot
#' @param map An **sf** polygon to plot in the background (optional)
#' @param label Whether to include node and tip numbers, default = FALSE
#' @param z The height of the tilt vector, default value is 10
#' @param alpha The alpha parameter, default value is -5
#' @export
plot_lineage <- function(tree, reconstruction, tip, map=NULL, label=FALSE, z=10, alpha=-5) {

  # Type checking
  if (!("dataTree" %in% class(tree))) {stop("Tree must be of class dataTree")}

  # Getting the reconstructed ancestors
  node_data <- reconstruction
  # Getting the tip ellipses
  ellipse_data <- tree$tip_data
  # Getting the edges from the tree
  tree_data <- data.frame(tree$tree$tree$edge)
  colnames(tree_data) <- c("N1", "N2")

  # Get taxon
  if (is.character(tip)) {
    tip <- which(tree$tree$tree$tip.label==tip)
  }
  taxon <- tree$tree$tree$tip.label[tip]

  # Tip ellipse
  ellipse_row <- which(ellipse_data$taxon==taxon)

  # Get the list of nodes leading to our desired species
  node_history <- ape::nodepath(tree$tree$tree,from=tip,to=(length(tree$tree$tree$tip.label)+1))
  node_history <- rev(node_history)

  # Getting the subset of node data for this history
  rows <- c()
  for (i in 1:(length(node_history)-1)) {
     row_ind <- which(node_data$node == node_history[i])
     rows <- c(rows,row_ind)
  }
  lineage_data <- node_data[rows,]

  # Determine left and right daughter assignments
  daughter_assignments <- rep(NA,length(node_history)-1)
  for (i in 1:(length(node_history)-1)) {
    parent <- node_history[i]
    right_child <- tree$tree$nodes$right_child[parent]
    left_child <- tree$tree$nodes$left_child[parent]
    if (node_history[i+1] == right_child) {daughter_assignments[i] <- "right"}
    if (node_history[i+1] == left_child) {daughter_assignments[i] <- "left"}
  }
  lineage_data$daughter <- daughter_assignments

  # Put together data for all nodes
  labels_data <- data.frame(matrix(data=NA,nrow=0,ncol=3))
  segments_data <- data.frame(matrix(data=NA,nrow=0,ncol=5))
  polygons_data <- data.frame(matrix(data=NA,nrow=0,ncol=4))
  for (i in 1:nrow(lineage_data)) {
    # Adding node to labels
    label_one <- data.frame(matrix(data=c(lineage_data$x[i],lineage_data$y[i],lineage_data$node[i]),nrow=1,ncol=3))
    labels_data <- rbind(labels_data,label_one)
    # Adding anagenetic segment if applicable
    if (exists("daughter_x")) {
      segment <- data.frame(matrix(data=c(daughter_x,daughter_y,lineage_data$x[i],lineage_data$y[i],"anagenetic"),nrow=1,ncol=5))
      segments_data <- rbind(segments_data,segment)
    }
    # Making ancestral ellipse
    ancestor_ellipse <- data.frame(make_ellipse_coords(lineage_data$x[i],
                                                       lineage_data$y[i],
                                                       lineage_data$r[i],
                                                       lineage_data$s[i],
                                                       lineage_data$a[i]))
    n_points <- nrow(ancestor_ellipse)
    ancestor_ellipse$node <- rep(paste0(lineage_data$node[i],"A"),n_points)
    ancestor_ellipse$type <- rep("ancestor",n_points)
    polygons_data <- rbind(polygons_data,ancestor_ellipse)
    # Making daughter ellipse
    daughter_x <- transform_node(char="x",
                                 daughter=lineage_data$daughter[i],
                                 d=lineage_data$d[i],
                                 m=lineage_data$m[i],
                                 cval=lineage_data$cval[i],
                                 hval=lineage_data$hval[i],
                                 x=lineage_data$x[i],
                                 y=lineage_data$y[i],
                                 r=lineage_data$r[i],
                                 s=lineage_data$s[i],
                                 a=lineage_data$a[i],
                                 alpha=alpha,
                                 z=z)
    daughter_y <- transform_node(char="y",
                                 daughter=lineage_data$daughter[i],
                                 d=lineage_data$d[i],
                                 m=lineage_data$m[i],
                                 cval=lineage_data$cval[i],
                                 hval=lineage_data$hval[i],
                                 x=lineage_data$x[i],
                                 y=lineage_data$y[i],
                                 r=lineage_data$r[i],
                                 s=lineage_data$s[i],
                                 a=lineage_data$a[i],
                                 alpha=alpha,
                                 z=z)
    daughter_a <- transform_node(char="a",
                                 daughter=lineage_data$daughter[i],
                                 d=lineage_data$d[i],
                                 m=lineage_data$m[i],
                                 cval=lineage_data$cval[i],
                                 hval=lineage_data$hval[i],
                                 x=lineage_data$x[i],
                                 y=lineage_data$y[i],
                                 r=lineage_data$r[i],
                                 s=lineage_data$s[i],
                                 a=lineage_data$a[i],
                                 alpha=alpha,
                                 z=z)
    daughter_ellipse <- data.frame(make_ellipse_coords(daughter_x,
                                                       daughter_y,
                                                       lineage_data$r[i],
                                                       lineage_data$s[i],
                                                       daughter_a))
    daughter_ellipse$node <- rep(paste0(lineage_data$node[i],"D"),n_points)
    daughter_ellipse$type <- rep("daughter",n_points)
    lineage_data$daughter_x[i] <- daughter_x
    lineage_data$daughter_y[i] <- daughter_y
    polygons_data <- rbind(polygons_data,daughter_ellipse)
    # Adding cladogenetic segment
    segment <- data.frame(matrix(data=c(lineage_data$x[i],lineage_data$y[i],daughter_x,daughter_y,"cladogenetic"),nrow=1,ncol=5))
    segments_data <- rbind(segments_data,segment)
  }

  # Adding ellipse for the tip
  tip_ellipse <- data.frame(make_ellipse_coords(ellipse_data$x[ellipse_row],
                                                ellipse_data$y[ellipse_row],
                                                ellipse_data$r[ellipse_row],
                                                ellipse_data$s[ellipse_row],
                                                ellipse_data$a[ellipse_row]))
  tip_ellipse$node <- rep(paste0(tip,"A"),n_points)
  tip_ellipse$type <- rep("ancestor",n_points)
  polygons_data <- rbind(polygons_data,tip_ellipse)
  tip_label <- data.frame(matrix(data=c(ellipse_data$x[ellipse_row],ellipse_data$y[ellipse_row],tip),nrow=1,ncol=3))
  labels_data <- rbind(labels_data,tip_label)
  tip_segment <- data.frame(matrix(data=c(daughter_x,daughter_y,ellipse_data$x[ellipse_row],ellipse_data$y[ellipse_row],"anagenetic"),nrow=1,ncol=5))
  segments_data <- rbind(segments_data,tip_segment)
  # Formatting output
  polygons_data$node <- as.factor(polygons_data$node)
  polygons_data$type <- as.factor(polygons_data$type)
  colnames(labels_data) <- c("x","y","label")
  colnames(segments_data) <- c("x","y","xend","yend","type")

  map_data <- polygons_data
  map_data$node <- as.factor(map_data$node)
  map_data$type <- as.factor(map_data$type)
  labels_data$label <- as.character(labels_data$label)
  segments_data$type <- as.factor(segments_data$type)
  segments_data$x <- as.numeric(segments_data$x)
  segments_data$y <- as.numeric(segments_data$y)
  segments_data$xend <- as.numeric(segments_data$xend)
  segments_data$yend <- as.numeric(segments_data$yend)

  # Getting extent of polygons for ancestor & daughters
  xmin <- min(map_data$x)
  xmax <- max(map_data$x)
  ymin <- min(map_data$y)
  ymax <- max(map_data$y)
  xscale <- (xmax - xmin) * 0.1
  yscale <- (ymax - ymin) * 0.1

  # If there is no map provided, create one based on extent of polygon data
  if (is.null(map)) {
    coords <- matrix(c(xmin - xscale, ymin - yscale,
                       xmax + xscale, ymin - yscale,
                       xmax + xscale, ymax + yscale,
                       xmin - xscale, ymax + yscale,
                       xmin - xscale, ymin - yscale),
                     ncol=2,
                     byrow=TRUE)
    map = sf::st_as_sf(sf::st_sfc(sf::st_polygon(list(cbind(coords[,1],coords[,2])))))
  }

  # Plot lineage on a map
  lineage_map <- ggplot2::ggplot() +
    # Plot the map
    ggplot2::geom_sf(data=map, fill=NA, color="black") +
    # Plot the polygons for ancestors and daughters
    ggplot2::geom_polygon(data=map_data,ggplot2::aes(x=x,y=y,group=node,color=type,lwd=type),fill=NA) +
    # Color the polygons based on ancestor or daughter
    ggplot2::scale_color_manual(values=c("black","gray"),
                                guide=ggplot2::guide_legend(override.aes=list(color=NA,fill=c("black","gray")))) +
    ggplot2::scale_linewidth_manual(values=c(2,1),guide=NULL) +
    # Plot the connection lines between evolving ellipses
    ggplot2::geom_segment(data=segments_data,ggplot2::aes(x=x,y=y,xend=xend,yend=yend,lty=type),color="black",lwd=1,arrow=ggplot2::arrow(length=ggplot2::unit(0.25,"cm")),lineend="round") +
    # Line type based on cladogenetic or anagenetic movement type
    ggplot2::scale_linetype_manual(values=c("solid","12")) +
    # Potentially label ellipses (gets messy)
    {if(label)ggplot2::geom_text(data=labels_data,ggplot2::aes(x=x,y=y,label=label),vjust=1.2,position="nudge")} +
    ggplot2::labs(color=NULL,lty=NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill="white",color=NA),
                   legend.position="right",
                   legend.title = ggplot2::element_text(color="black",size=12),
                   legend.text = ggplot2::element_text(color="black",size=8),
                   legend.box.background = ggplot2::element_rect(fill="white", color=NA),
                   legend.margin = ggplot2::margin(4,4,4,4),
                   legend.spacing.y = ggplot2::unit(.1,"cm"))

  return(lineage_map)
}

#' @title plot_posterior() function
#' @description Processes output of an ellipse MCMC and generates a posterior distribution for the ancestral ellipse and both daughter ellipses at a single node
#' @param filepath An object of class "character" (a string) that gives the filepath where the output log files from your MCMC analysis are saved. Default behavior will look for the files in the current directory.
#' @param nodes A **numeric** node number or vector node numbers representing the phylogenetic node(s) to plot
#' @param map An **sf** polygon to plot in the background (optional)
#' @param burnin The number of burnin *rows* produced during the MCMC analysis (burnin divided by thinning); will be removed before processing
#' @param thinning Retain every N iterations for plotting the posterior. The default value is 1, retaining all iterations; this may be prohibitively slow with long MCMC chains
#' @param cval An object of class "vector" containing any number of numeric elements indicating the values associated with different concentric circles relative to the size of the original circle. The default value includes 4 concentric circles: 0, 0.5, 1.0, and 1.5.
#' @param hval An object of class "vector" containing any number of numeric elements indicating the values associated with different direction lines in radians. The default value includes 8 direction lines: 0, pi/4, pi/2, 3pi/4, pi, 5pi/4, 3pi/2, and 7pi/4.
#' @param z The height of the tilt vector, default value is 10
#' @param alpha The alpha parameter, default value is -5
#' @export
plot_posterior <- function(filepath=".", nodes, map=NULL, burnin=0, thinning=1, cval=NULL, hval=NULL, z=10, alpha=-5) {

  # Default values associated with c and h indices
  if (is.null(cval)) {cval=seq(0,1.5,.5)}
  if (is.null(hval)) {hval=2*pi*seq(0,7/8,1/8)}

  # Getting the number of rows in the MCMC output and saving as "length"
  data_d_pre <- read.table(paste0(filepath,"/d_log.tsv"),sep="\t",header=TRUE)
  length <- nrow(data_d_pre)

  # Reading the MCMC output for each ellipse character, removing burnin lines & thinning
  data_d <- read.table(paste0(filepath,"/d_log.tsv"),sep="\t",header=TRUE)[seq(burnin,length,thinning),]
  data_m <- read.table(paste0(filepath,"/m_log.tsv"),sep="\t",header=TRUE)[seq(burnin,length,thinning),]
  data_c <- read.table(paste0(filepath,"/c_log.tsv"),sep="\t",header=TRUE)[seq(burnin,length,thinning),]
  data_h <- read.table(paste0(filepath,"/h_log.tsv"),sep="\t",header=TRUE)[seq(burnin,length,thinning),]
  data_r <- read.table(paste0(filepath,"/r_log.tsv"),sep="\t",header=TRUE)[seq(burnin,length,thinning),]
  data_s <- read.table(paste0(filepath,"/s_log.tsv"),sep="\t",header=TRUE)[seq(burnin,length,thinning),]
  data_a <- read.table(paste0(filepath,"/a_log.tsv"),sep="\t",header=TRUE)[seq(burnin,length,thinning),]
  data_x <- read.table(paste0(filepath,"/x_log.tsv"),sep="\t",header=TRUE)[seq(burnin,length,thinning),]
  data_y <- read.table(paste0(filepath,"/y_log.tsv"),sep="\t",header=TRUE)[seq(burnin,length,thinning),]

  # Starting the list of plots for each node
  plot_list <- list()

  # Looping over each node to perform a reconstruction
  for (node in nodes) {

    # Getting the name of the node in the output files
    # Equals the node number preceded by an X
    # And also the numbers of both children
    node_name <- paste0("X",node)
    left_name <- tree$tree$nodes$left_child[which(tree$tree$nodes$node==node)]
    right_name <- tree$tree$nodes$right_child[which(tree$tree$nodes$node==node)]
    print(paste0("Node: ",node),quote=F)

    # Getting the data columns associated with the node
    node_dataframe <- data.frame(cbind(d=data_d[[node_name]],
                                       m=data_m[[node_name]],
                                       c=data_c[[node_name]],
                                       h=data_h[[node_name]],
                                       r=data_r[[node_name]],
                                       s=data_s[[node_name]],
                                       a=data_a[[node_name]],
                                       x=data_x[[node_name]],
                                       y=data_y[[node_name]]))

    # Getting the actual c and h values associated with the indices in the MCMC output
    for (i in 1:nrow(node_dataframe)) {node_dataframe$cval[i] <- cval[node_dataframe$c[i]+1]}
    for (i in 1:nrow(node_dataframe)) {node_dataframe$hval[i] <- hval[node_dataframe$h[i]+1]}

    # Function for obtaining the polygon associated with a daughter based on ancestral values
    get_daughter_polygon <- function(daughter,d,m,cval,hval,x,y,r,s,a) {
      new_x <- transform_node(char="x", daughter=daughter, d=d, m=m, cval=cval, hval=hval, x=x, y=y, r=r, s=s, a=a, alpha=-5, z=10)
      new_y <- transform_node(char="y", daughter=daughter, d=d, m=m, cval=cval, hval=hval, x=x, y=y, r=r, s=s, a=a, alpha=-5, z=10)
      new_a <- transform_node(char="a", daughter=daughter, d=d, m=m, cval=cval, hval=hval, x=x, y=y, r=r, s=s, a=a, alpha=-5, z=10)
      new_polygon <- data.frame(make_ellipse_coords(new_x,new_y,r,s,new_a))
      return(new_polygon)
    }

    # Starting dataframes for ancestor and daughter polygons for each MCMC iteration
    ancestor_data <- data.frame(matrix(data=NA,nrow=0,ncol=3))
    left_data <- data.frame(matrix(data=NA,nrow=0,ncol=4))
    right_data <- data.frame(matrix(data=NA,nrow=0,ncol=4))

    # Looping over the MCMC iterations (possibly thinned)
    for (i in 1:nrow(node_dataframe)) {
      # Getting the ancestral ellipse
      ancestor_ellipse <- data.frame(make_ellipse_coords(node_dataframe$x[i],node_dataframe$y[i],node_dataframe$r[i],node_dataframe$s[i],node_dataframe$a[i]))
      n_points <- nrow(ancestor_ellipse)
      ancestor_ellipse$index <- rep(i,n_points)

      # Getting the left daughter ellipse
      left_ellipse <- get_daughter_polygon("left",node_dataframe$d[i],node_dataframe$m[i],node_dataframe$cval[i],node_dataframe$hval[i],node_dataframe$x[i],node_dataframe$y[i],node_dataframe$r[i],node_dataframe$s[i],node_dataframe$a[i])
      n_points <- nrow(left_ellipse)
      left_ellipse$index <- rep(i,n_points)

      # Getting the right daughter ellipse
      right_ellipse <- get_daughter_polygon("right",node_dataframe$d[i],node_dataframe$m[i],node_dataframe$cval[i],node_dataframe$hval[i],node_dataframe$x[i],node_dataframe$y[i],node_dataframe$r[i],node_dataframe$s[i],node_dataframe$a[i])
      n_points <- nrow(right_ellipse)
      right_ellipse$index <- rep(i,n_points)

      # Adding ancestor, left daughter, and right daughter to dataframes
      ancestor_data <- rbind(ancestor_data,ancestor_ellipse)
      left_data <- rbind(left_data,left_ellipse)
      right_data <- rbind(right_data,right_ellipse)
    }

    # Creating multipolygons out of ancestor, left daughter, and right daughter data
    ancestor_shapes <- sf::as_Spatial(sfheaders::sf_multipolygon(ancestor_data,x="x",y="y",polygon_id="index",multipolygon_id="index"))
    left_shapes <- sf::as_Spatial(sfheaders::sf_multipolygon(left_data,x="x",y="y",polygon_id="index",multipolygon_id="index"))
    right_shapes <- sf::as_Spatial(sfheaders::sf_multipolygon(right_data,x="x",y="y",polygon_id="index",multipolygon_id="index"))

    # Getting the extent of the ancestor/daughter datasets and getting limits
    anc_box <- sp::bbox(ancestor_shapes)
    left_box <- sp::bbox(left_shapes)
    right_box <- sp::bbox(right_shapes)

    xmin <- min(anc_box[1],left_box[1],right_box[1])
    ymin <- min(anc_box[2],left_box[2],right_box[2])
    xmax <- max(anc_box[3],left_box[3],right_box[3])
    ymax <- max(anc_box[4],left_box[4],right_box[4])
    xscale <- (xmax - xmin) * 0.1
    yscale <- (ymax - ymin) * 0.1

    # If there is no map provided, create one based on extent of polygon data
    if (is.null(map)) {
      coords <- matrix(c(xmin - xscale, ymin - yscale,
                         xmax + xscale, ymin - yscale,
                         xmax + xscale, ymax + yscale,
                         xmin - xscale, ymax + yscale,
                         xmin - xscale, ymin - yscale),
                       ncol=2,
                       byrow=TRUE)
      new_map = sf::st_as_sf(sf::st_sfc(sf::st_polygon(list(cbind(coords[,1],coords[,2])))))
    } else {
      new_map <- map
    }

    # Creates a grid to use for mapping the density of polygons
    # There will be 20 grid cells in the longest map direction
    size <- max(c(xmax-xmin),(ymax-ymin))
    scale <- size / 20
    extent <- raster::extent(c(xmin,xmax,ymin,ymax))
    grid <- raster::raster(extent,resolution=scale)
    raster_data <- raster::stack(grid)

    # Rasterizing the data for the ancestor, left daughter, and right daughter
    ancestor_raster <- raster::rasterize(ancestor_shapes,raster_data,"index",fun="count",background=0)
    left_raster <- raster::rasterize(left_shapes,raster_data,"index",fun="count",background=0)
    right_raster <- raster::rasterize(right_shapes,raster_data,"index",fun="count",background=0)
    names(ancestor_raster) <- "ancestor"
    names(left_raster) <- "left"
    names(right_raster) <- "right"
    raster_with_counts <- raster::addLayer(raster_data,ancestor_raster)
    raster_with_counts <- raster::addLayer(raster_with_counts,left_raster)
    raster_with_counts <- raster::addLayer(raster_with_counts,right_raster)

    # Turning the raster data into a dataframe
    count_df <- raster::as.data.frame(raster_with_counts,xy=TRUE)
    count_df$ancestor <- count_df$ancestor / max(count_df$ancestor)
    count_df$left <- count_df$left / max(count_df$left)
    count_df$right <- count_df$right / max(count_df$right)
    count_df$ancestor[count_df$ancestor < 0.05] <- 0
    count_df$left[count_df$left < 0.05] <- 0
    count_df$right[count_df$right < 0.05] <- 0

    # Getting the counts for each combination of discrete values for cladogenetic scenarios
    scenario_counts <- dplyr::summarize(.data=node_dataframe,count=dplyr::n(),.by=c(d,m,c,h))
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
    scenario_data <- node_dataframe[which(node_dataframe$d==scenario_d &
                                            node_dataframe$m==scenario_m &
                                            node_dataframe$c==scenario_c &
                                            node_dataframe$h==scenario_h),]
    # Getting the mean values for continuous characters under that scenario
    scenario_r <- mean(scenario_data$r)
    scenario_s <- mean(scenario_data$s)
    scenario_a <- mean(scenario_data$a)
    scenario_x <- mean(scenario_data$x)
    scenario_y <- mean(scenario_data$y)
    # Creating polygons for reconstructed ancestor & daughters
    best_ancestor <- data.frame(make_ellipse_coords(scenario_x,
                                                    scenario_y,
                                                    scenario_r,
                                                    scenario_s,
                                                    scenario_a))
    best_left <- get_daughter_polygon("left",
                                      scenario_d,
                                      scenario_m,
                                      scenario_cval,
                                      scenario_hval,
                                      scenario_x,
                                      scenario_y,
                                      scenario_r,
                                      scenario_s,
                                      scenario_a)
    best_right <- get_daughter_polygon("right",
                                       scenario_d,
                                       scenario_m,
                                       scenario_cval,
                                       scenario_hval,
                                       scenario_x,
                                       scenario_y,
                                       scenario_r,
                                       scenario_s,
                                       scenario_a)

    # Plotting the posterior distribution for the ancestor
    # Color = purple; best ellipse = black
    plot_ancestor <- ggplot2::ggplot() +
      ggplot2::geom_sf(data=new_map, fill=NA, color="black") +
      ggplot2::geom_tile(data=count_df,ggplot2::aes_string(x="x",y="y",alpha="ancestor"),fill="purple") +
      ggplot2::scale_alpha_continuous(range=c(0,1)) +
      ggplot2::labs(alpha="Probability") +
      ggplot2::geom_sf(fill=NA,color="black",lwd=1) +
      ggplot2::geom_polygon(data=best_ancestor,ggplot2::aes(x=x,y=y),fill=NA,color="black",lwd=1.25) +
      ggplot2::ggtitle(label=paste0("Ancestor\n(",node,")")) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill="white",color=NA),
                     panel.background = ggplot2::element_rect(fill="white",color=NA),
                     legend.position="none",
                     plot.title = ggplot2::element_text(color="black",size=12,hjust=.5,vjust=1))

    # Plotting the posterior distribution for the left daughter
    # Color = blue; best ellipse = black
    plot_left <- ggplot2::ggplot() +
      ggplot2::geom_sf(data=new_map, fill=NA, color="black") +
      ggplot2::geom_tile(data=count_df,ggplot2::aes_string(x="x",y="y",alpha="left"),fill="blue") +
      ggplot2::scale_alpha_continuous(range=c(0,1)) +
      ggplot2::labs(alpha="Probability") +
      ggplot2::geom_sf(fill=NA,color="black",lwd=1) +
      ggplot2::geom_polygon(data=best_left,ggplot2::aes(x=x,y=y),fill=NA,color="black",lwd=1.25) +
      ggplot2::ggtitle(label=paste0("Left Daughter\n(Toward ",left_name,")")) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill="white",color=NA),
                     panel.background = ggplot2::element_rect(fill="white",color=NA),
                     legend.position="none",
                     plot.title = ggplot2::element_text(color="black",size=12,hjust=.5,vjust=1))

    # Plotting the posterior distribution for the right daughter
    # Color = red; best ellipse = black
    plot_right <- ggplot2::ggplot() +
      ggplot2::geom_sf(data=new_map, fill=NA, color="black") +
      ggplot2::geom_tile(data=count_df,ggplot2::aes_string(x="x",y="y",alpha="right"),fill="red") +
      ggplot2::scale_alpha_continuous(range=c(0,1)) +
      ggplot2::labs(alpha="Probability") +
      ggplot2::geom_sf(fill=NA,color="black",lwd=1) +
      ggplot2::geom_polygon(data=best_right,ggplot2::aes(x=x,y=y),fill=NA,color="black",lwd=1.25) +
      ggplot2::ggtitle(label=paste0("Right Daughter\n(Toward ",right_name,")")) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill="white",color=NA),
                     panel.background = ggplot2::element_rect(fill="white",color=NA),
                     legend.position="none",
                     plot.title = ggplot2::element_text(color="black",size=12,hjust=.5,vjust=1))

    # Plotting the reconstruction for the best cladogenetic scenario at the node
    # The ancestor is in purple, left daughter in blue, and right daughter in red
    plot_best <- ggplot2::ggplot() +
      ggplot2::geom_sf(data=new_map, fill=NA, color="black") +
      ggplot2::geom_polygon(data=best_ancestor,ggplot2::aes(x=x,y=y),fill="purple",alpha=0.5,color="purple",lwd=2) +
      ggplot2::geom_polygon(data=best_left,ggplot2::aes(x=x,y=y),fill=NA,color="blue",lwd=2,lty="11") +
      ggplot2::geom_polygon(data=best_left,ggplot2::aes(x=x,y=y),fill=NA,color="blue",lwd=.75) +
      ggplot2::geom_polygon(data=best_right,ggplot2::aes(x=x,y=y),fill=NA,color="red",lwd=2,lty="11") +
      ggplot2::geom_polygon(data=best_right,ggplot2::aes(x=x,y=y),fill=NA,color="red",lwd=.75) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill="white",color=NA),
                     panel.background = ggplot2::element_rect(fill="white",color=NA),
                     legend.position="none")

    node_posterior_plot <- (plot_ancestor + plot_left + plot_right) / plot_best + patchwork::plot_layout(heights=c(1,3.2)) & ggplot2::theme(plot.background=ggplot2::element_rect(fill="white",color=NA))
    ggplot2::ggsave(paste0(filepath,"/","node_",node,".pdf"),node_posterior_plot,dpi=600,limitsize=FALSE)

    plot_list <- append(plot_list,list(node_posterior_plot))
  }

  return(plot_list)
}
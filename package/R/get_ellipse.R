#' @title get_ellipse() function
#' @description Takes a spatial POLYGON or MULTIPOLYGON (sf package) and generates ellipse parameters x, y, r, s, and a. Returns a vector (x,y,r,s,a,z). This is done by sampling points from within the spatial polygon/multipolygon, then creating a confidence ellipse around those points using the stats function cov.wt(). Alternatively, can provide a set of occurrence points in a data.frame or MULTIPOINT (sf) object.
#' @param data The spatial POLYGON or MULTIPOLYGON to turn into an ellipse. Alternatively, the occurrence points to turn into an ellipse, as a data.frame or MULTIPOINT object.
#' @param z The height of the tilt vector for generating the ellipses, default value is 10
#' @param confidence The level of confidence used to generate the confidence ellipse, default value is 0.95
#' @export
get_ellipse <- function(data=NULL,z=10,confidence=0.95) {

  # Checking the type of the input
  # If polygons are provided, sampling points uniformly at random from interior (10000 points)
  # Creating coordinates dataframe from points
  if ("MULTIPOINT" %in% class(data)) {
    coords <- sf::st_coordinates(data)
  } else if ("POLYGON" %in% class(data)) {
    sample <- sf::st_sample(data,10000)
    coords <- sf::st_coordinates(sample)
  } else if ("MULTIPOLYGON" %in% class(data)) {
    sample <- sf::st_sample(data,10000)
    coords <- sf::st_coordinates(sample)
  } else if ("data.frame" %in% class(data)) {
    coords <- data
  } else {
    stop("Input data must be of class MULTIPOINT, POLYGON, MULTIPOLYGON, or data.frame",call.=FALSE)
  }

  # Type checking input for 'z' (height of tilt vector used to project an ellipse, default is 10) <-- make sure this matches during analysis!
  # Type checking input for 'confidence' (the confidence level used to generate the confidence ellipse, default is 0.95)
  if (!(class(z)=="numeric")) {stop("z must be numeric",call.=FALSE)}
  if (!(class(confidence)=="numeric")) {stop("confidence must be numeric",call.=FALSE)}
  if (!(confidence <= 1)) {stop("confidence must be less than or equal to 1",call.=FALSE)}
  if (!(confidence > 0)) {stop("confidence must be greater than 0",call.=FALSE)}

  # If a single occurrence point is provided, returning a very small circle [a=-50]
  if (nrow(coords)==1) {
    # Returning (x, y, r, s, a, z) for an ellipse
    return(list(x=coords[1],y=coords[2],r=0,s=0,a=-50,z=z))
  }

  # Getting the weighted covariance matrix for the points
  # This information is used for a PCA
  info <- cov.wt(coords)

  # Computing eigenvectors and eigenvalues from the weighted covariance matrix
  # Also computing the ellipse radii (major and minor)
  # Ellipse radii are based on eigenvalues (variance of the data in the direction of the eigenvector) and ellipse confidence level
  # Confidence intervals are calculated using chi square distribution with df = 2 (variation in two dimensions)
  # See ConfidenceEllipse package for similar procedure: https://cran.r-project.org/web/packages/ConfidenceEllipse/index.html
  eig <- eigen(info$cov)
  lengths <- sqrt(eig$values * qchisq(confidence, 2))

  # Getting point and angle to use for solving for r,s
  # We are interested in the major axis (eigenvector 1), which has two points on the ellipse
  # We will select the point with the higher x value (both points are explored during MCMC when using the ellipse evolution model)
  pt1 <- lengths[1] * eig$vectors[,1]
  pt2 <- -1 * lengths[1] * eig$vectors[,1]
  point <- if (pt1[1] > pt2[1]) {pt1} else {pt2}
  angle <- atan(point[2]/point[1])

  # Stuff for finding r,s
  # Some intermediates for the projection equation:
  A <- lengths[2] * cos(angle)
  B <- lengths[2] * sin(angle)
  X <- point[1]
  Y <- point[2]

  # Values for r and s are the x and y coordinates of the tilt vector that produces the projected ellipse, at height z
  # x and y are the coordinates of the ellipse's centroid, found from the cov.wt() information
  # 'a' is log area of the projecting circle that forms the ellipse
  # The size (radius) of the projecting circle is equal to the minor radius of the ellipse
  s <- sqrt( (z^2 * (Y - B))  / (A * (X - A) / (Y - B) + B))
  r <- s * (X - A) / (Y - B)
  x <- info$center[[1]]
  y <- info$center[[2]]
  a <- log(pi * lengths[2]^2)

  # Returning (x, y, r, s, a, z) for an ellipse
  return(list(x=x,y=y,r=r,s=s,a=a,z=z))
}

#' @title get_ellipses() function
#' @description Takes a spatial MULTIPOLYGON (sf package) with entries for each species and generates ellipse parameters x, y, r, s, and a for each species. Returns an object of class "data.frame" (taxon,x,y,r,s,a,z). This is done by sampling points from within each spatial polygon/multipolygon, then creating a confidence ellipse around those points using the stats function cov.wt(). Alternatively, can provide a set of occurrence points in a MULTIPOINT (sf) object.
#' @param data The spatial MULTIPOLYGON to turn into ellipses. Alternatively, the occurrence points to turn into ellipses, as a MULTIPOINT object. The first column of this sf object must contain the taxon names. The column containing data (polygons or points) must be called "geometry".
#' @param z The height of the tilt vector for generating the ellipses, default value is 10
#' @param confidence The level of confidence used to generate the confidence ellipse, default value is 0.95
#' @export
get_ellipses <- function(data=NULL,z=10,confidence=0.95) {

  # Getting taxa names from spatial multipolygon or multipoint, assuming that the first column of the sf object has taxon names
  taxa_list <- data[[1]]

  # Starting the dataframe to record values for each ellipse
  out <- data.frame(matrix(ncol=6,nrow=length(taxa_list)))
  colnames(out) <- c("taxon", "x", "y", "r", "s", "a")

  # Looping through the taxa in the multipolygon/multipoint
  for (i in 1:length(taxa_list)) {
    # Getting the data (this must be called 'geometry')
    shape <- data$geometry[[i]]
    # Getting the ellipse
    ellipse <- get_ellipse(shape,z,confidence)
    # Recording the taxon name
    out[i,1] <- taxa_list[i]
    # Recording the ellipse parameters
    out[i,2:6] <- c(ellipse$x,ellipse$y,ellipse$r,ellipse$s,ellipse$a)
  }

  # Returning dataframe with (taxon, x, y, r, s, a, z) for all taxa
  return(out)
}

#' @title make_ellipse_coords() function
#' @description Takes a set of ellipse parameters (x,y,r,s,a,z) and returns a dataframe of points defining the perimeter of the ellipse anticlockwise.
#' @param x The x coordinate of the ellipse centroid
#' @param y The y coordinate of the ellipse centroid
#' @param r The x coordinate of the ellipse tilt point (when combined with s and z, defines oblongness)
#' @param r The y coordinate of the ellipse tilt point (when combined with r and z, defines oblongness)
#' @param z The height of the tilt vector (when combined with r and s, defines oblongness), default value is 10
#' @export
make_ellipse_coords <- function(x,y,r,s,a,z=10) {

  # Type checking
  if (!("numeric" %in% class(x))) {stop("x must be numeric",call.=FALSE)}
  if (!("numeric" %in% class(y))) {stop("y must be numeric",call.=FALSE)}
  if (!("numeric" %in% class(r))) {stop("r must be numeric",call.=FALSE)}
  if (!("numeric" %in% class(s))) {stop("s must be numeric",call.=FALSE)}
  if (!("numeric" %in% class(a))) {stop("a must be numeric",call.=FALSE)}
  if (!("numeric" %in% class(z))) {stop("z must be numeric",call.=FALSE)}

  # The radius of the projecting circle is equal to the minor radius of the resulting ellipse
  radius <- sqrt(exp(a)/pi)
  # To draw the ellipse, we will create a series of points at equal angles around the ellipse and connect them
  # Theta will be our vector of angles to use
  theta <- seq(0, 2*pi, pi/64)
  # fx, fy, and fz will define the locations of the points on the tilted plane (not yet flattened into real space)
  # note that fx and fy are the x and y coordinates of the projecting circle
  fx <- function(theta){radius*cos(theta)}
  fy <- function(theta){radius*sin(theta)}
  fz <- function(theta){(r*fx(theta)+s*fy(theta))}
  # proj_x and proj_y project the points from the tilted plane into real space
  # These are the x and y coordinates of the ellipse points in real space
  proj_x <- lapply(theta,function(theta){x+fx(theta)+r*fz(theta)/z^2})
  proj_y <- lapply(theta,function(theta){y+fy(theta)+s*fz(theta)/z^2})
  return(data.frame(cbind(x=as.numeric(proj_x),y=as.numeric(proj_y))))
}

#' @title transform_node() function
#' @description Accepts model parameters (alpha, z), information about a cladogenetic scenario (d, m, c, h), data-augmented values (r, s), and values for a character prior to cladogenesis (x, y, a) and deterministically calculates the new value of the character immediately after cladogenesis
#' @param char One of three continuous characters: "x", "y", or "a"
#' @param daughter One of two daughters: "left", "right"
#' @param d The daughter configuration during cladogenesis (0,1), where 0 makes the left daughter D1 and 1 makes the right daughter D1
#' @param m The cladogenetic mode (0,1), where 0 is budding and 1 is splitting
#' @param cval The c value (not index) of the cladogenetic scenario's concentric circle --> ex. 0.5
#' @param hval The h value (not index) of the cladogenetic scenario's direction line --> ex. pi/4, if required
#' @param x The x value to be transformed, if required
#' @param y The y value to be transformed, if required
#' @param r The data-augmented 'r' value associated with the cladogenesis event, if required
#' @param s The data-augmented 's' value associated with the cladogenesis event, if required
#' @param a The a value to be transformed, if required
#' @param alpha The 'alpha' model parameter, if required
#' @param z The 'z' model parameter, if required
#' @export
transform_node = function(char=NULL, daughter=NULL, d=NULL, m=NULL, cval=NULL, hval=NULL, x=NULL, y=NULL, r=NULL, s=NULL, a=NULL, alpha=NULL, z=NULL) {

  # make sure the character we are transforming is x, y, or a
  if (!(char %in% c("x", "y", "a"))) {stop("char must be x, y, or a",call.=FALSE)}
  # make sure the daughter we are transforming is left or right
  if (!(daughter %in% c("left","right"))) {stop("daughter must be left or right",call.=FALSE)}
  # make sure we have daughters assigned correctly
  if (!(d %in% c(0,1))) {stop("d must be 0 or 1",call.=FALSE)}
  # make sure we have the cladogenetic mode assigned correctly
  if (!(m %in% c(0,1))) {stop("m must be 0 or 1",call.=FALSE)}
  # make sure a is numeric
  if (!("numeric" %in% class(a))) {stop("a must be numeric",call.=FALSE)}
  # make sure c is numeric
  if (!("numeric" %in% class(cval))) {stop("cval must be numeric",call.=FALSE)}

  # Function for determining the change in x based on the relative radius of the new daughter circle, r, s, a, z, and theta (h or h+pi)
  # Only used if we are looking at the x character
  delta_x <- function(rad, r, s, a, z, theta) {
    # Getting the magnitude of change in unprojected "circle space"
    radius <- exp((1/2) * a + log(rad) - (1/2) * log(pi))
    # Getting the change in x in projected "ellipse space"
    delta <- radius * cos(theta) - r * (-(r * radius * cos(theta) + s * radius * sin(theta)) / z) / z
    return(delta)
  }

  # Function for determining the change in y based on the relative radius of the new daughter circle, r, s, a, z, and theta (h or h+pi)
  # Only used if we are looking at the y character
  delta_y <- function(rad, r, s, a, z, theta) {
    # Getting the magnitude of change in unprojected "circle space"
    radius <- exp((1/2) * a + log(rad) - (1/2) * log(pi))
    # Getting the change in y in projected "ellipse space"
    delta <- radius * sin(theta) - s * (-(r * radius * cos(theta) + s * radius * sin(theta)) / z) / z
    return(delta)
  }

  # If we are examining character 'x', we want to give a new value for 'x' after cladogenesis
  if (char=="x") {
    # make sure x is assigned and numeric
    if (!("numeric" %in% class(x))) {stop("x must be numeric",call.=FALSE)}
    # make sure other required params are assigned and numeric
    if (!("numeric" %in% class(hval))) {stop("hval must be numeric",call.=FALSE)}
    if (!("numeric" %in% class(r))) {stop("r must be numeric",call.=FALSE)}
    if (!("numeric" %in% class(s))) {stop("s must be numeric",call.=FALSE)}
    if (!("numeric" %in% class(z))) {stop("z must be numeric",call.=FALSE)}
    if (m==0) { # m = 0 (budding)
      if (d==0) { # d = 0 (D1 = left, D2 = right)
        if (daughter=="left") {return(x)}
        else if (daughter=="right") {return(x + delta_x(cval,r,s,a,z,hval))}}
      else if (d==1) { # d = 1 (D2 = left, D1 = right)
        if (daughter=="left") {return(x + delta_x(cval,r,s,a,z,hval))}
        else if (daughter=="right") {return(x)}}}
    else if (m==1) { # m = 1 (splitting)
      if (d==0) { # d = 0 (D1 = left, D2 = right)
        if (daughter=="left") {return(x + delta_x(sqrt(1/(1+(1+cval)^2)),r,s,a,z,(hval+pi)))} # Relative radius of new daughter circle is sqrt(1/(1+(1+c)^2))
        else if (daughter=="right") {return(x + delta_x(sqrt((1+cval)^2/(1+(1+cval)^2)),r,s,a,z,hval))}} # Relative radius of new daughter circle is sqrt((1+c)^2/(1+(1+c)^2))
      else if (d==1) { # d = 1 (D2 = left, D1 = right)
        if (daughter=="left") {return(x + delta_x(sqrt((1+cval)^2/(1+(1+cval)^2)),r,s,a,z,hval))} # Relative radius of new daughter circle is sqrt((1+c)^2/(1+(1+c)^2))
        else if (daughter=="right") {return(x + delta_x(sqrt(1/(1+(1+cval)^2)),r,s,a,z,(hval+pi)))}}}} # Relative radius of new daughter circle is sqrt(1/(1+(1+c)^2))

  # If we are examining character 'y', we want to give a new value for 'y' after cladogenesis
  else if (char=="y") {
    # make sure y is assigned and numeric
    if (!("numeric" %in% class(y))) {stop("y must be numeric",call.=FALSE)}
    # make sure other params are assigned and numeric
    if (!("numeric" %in% class(hval))) {stop("h must be numeric",call.=FALSE)}
    if (!("numeric" %in% class(r))) {stop("r must be numeric",call.=FALSE)}
    if (!("numeric" %in% class(s))) {stop("s must be numeric",call.=FALSE)}
    if (!("numeric" %in% class(z))) {stop("z must be numeric",call.=FALSE)}
    if (m==0) { # m = 0 (budding)
      if (d==0) { # d = 0 (D1 = left, D2 = right)
        if (daughter=="left") {return(y)}
        else if (daughter=="right") {return(y + delta_y(cval,r,s,a,z,hval))}}
      else if (d==1) { # d = 1 (D2 = left, D1 = right)
        if (daughter=="left") {return(y + delta_y(cval,r,s,a,z,hval))}
        else if (daughter=="right") {return(y)}}}
    else if (m==1) { # m = 1 (splitting)
      if (d==0) { # d = 0 (D1 = left, D2 = right)
        if (daughter=="left") {return(y + delta_y(sqrt(1/(1+(1+cval)^2)),r,s,a,z,(hval+pi)))} # Relative radius of new daughter circle is sqrt(1/(1+(1+c)^2))
        else if (daughter=="right") {return(y + delta_y(sqrt((1+cval)^2/(1+(1+cval)^2)),r,s,a,z,hval))}} # Relative radius of new daughter circle is sqrt((1+c)^2/(1+(1+c)^2))
      else if (d==1) { # d = 1 (D2 = left, D1 = right)
        if (daughter=="left") {return(y + delta_y(sqrt((1+cval)^2/(1+(1+cval)^2)),r,s,a,z,hval))} # Relative radius of new daughter circle is sqrt((1+c)^2/(1+(1+c)^2))
        else if (daughter=="right") {return(y + delta_y(sqrt(1/(1+(1+cval)^2)),r,s,a,z,(hval+pi)))}}}} # Relative radius of new daughter circle is sqrt(1/(1+(1+c)^2))

  # If we are examining character 'a', we want to give a new value for 'a' after cladogenesis (this is a log area)
  else if (char=="a") {
    # make sure alpha is assigned and numeric
    if (!("numeric" %in% class(alpha))) {stop("alpha must be numeric",call.=FALSE)}
    if (m==0) { # m = 0 (budding)
      if (d==0) { # d = 0 (D1 = left, D2 = right)
        if (daughter=="left") {return(a)}
        else if (daughter=="right") {return(alpha)}} # Alpha is small fixed value for budding scenario
      else if (d==1) { # d = 1 (D2 = left, D1 = right)
        if (daughter=="left") {return(alpha)} # Alpha is small fixed value for budding scenario
        else if (daughter=="right") {return(a)}}}
    else if (m==1) { # m = 1 (splitting)
      if (d==0) { # d = 0 (D1 = left, D2 = right)
        if (daughter=="left") {return(a + log((1+cval)^2/(1+(1+cval)^2)))} # Relative radius of new daughter circle is sqrt((1+c)^2/(1+(1+c)^2))
        else if (daughter=="right") {return(a + log(1/(1+(1+cval)^2)))}} # Relative radius of new daughter circle is sqrt(1/(1+(1+c)^2))
      else if (d==1) { # d = 1 (D2 = left, D1 = right)
        if (daughter=="left") {return(a + log(1/(1+(1+cval)^2)))} # Relative radius of new daughter circle is sqrt(1/(1+(1+c)^2))
        else if (daughter=="right") {return(a + log((1+cval)^2/(1+(1+cval)^2)))}}}} # Relative radius of new daughter circle is sqrt((1+c)^2/(1+(1+c)^2))
}
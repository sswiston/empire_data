# ellipseParam object -----------------------------------------------------

#' @title ellipseParam object
#' @description This class stores model parameters for the ellipse evolution model
#' @details This object is used during simulation to contain the model parameters used for generating data (W, V, XY). It is also used during analysis to store current parameter values.
  #' @field cval An object of class "vector" containing any number of numeric elements indicating the values associated with different concentric circles relative to the size of the original circle. The default value includes 4 concentric circles: 0, 0.5, 1.0, and 1.5.
  #' @field hval An object of class "vector" containing any number of numeric elements indicating the values associated with different direction lines in radians. The default value includes 8 direction lines: 0, pi/4, pi/2, 3pi/4, pi, 5pi/4, 3pi/2, and 7pi/4.
  #' @field rho_m An object of class "vector" containing two numeric elements: the prior probability of cladogenetic mode 1 (budding) and the prior probability of cladogenetic mode 2 (splitting). The default value assigns probability 0.5 to both modes.
  #' @field rho_c An object of class "vector" containing any number of numeric elements (equal to the length of cval) corresponding to the prior probabilities of different concentric circle assignments. The default value assigns probability 0.25 to each of 4 concentric circles.
  #' @field rho_h An object of class "vector" containing any number of numeric elements (equal to the length of hval) corresponding to the prior probabilities of different direction line assignments. The default value assigns probability 0.125 to each of 8 direction lines.
  #' @field sigma_x An object of class "numeric" containing the evolutionary rate of the continuous character 'x', which helps define an ellipse's center location. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
  #' @field sigma_y An object of class "numeric" containing the evolutionary rate of the continuous character 'y', which helps define an ellipse's center location. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
  #' @field sigma_r An object of class "numeric" containing the evolutionary rate of the continuous character 'r', which helps define an ellipse's "oblongness". If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
  #' @field sigma_s An object of class "numeric" containing the evolutionary rate of the continuous character 's', which helps define an ellipse's "oblongness". If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
  #' @field sigma_a An object of class "numeric" containing the evolutionary rate of the continuous character 'a', which helps define an ellipse's size. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 5.
  #' @field root_x An object of class "numeric" containing the root value for the 'x' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and variance 1.
  #' @field root_y An object of class "numeric" containing the root value for the 'y' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and variance 1.
  #' @field root_r An object of class "numeric" containing the root value for the 'r' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and variance 5.
  #' @field root_s An object of class "numeric" containing the root value for the 's' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and variance 5.
  #' @field root_a An object of class "numeric" containing the root value for the 'a' character. If left blank, the default behavior will draw a value from a uniform distribution from -5 to 10.
  #' @field mu An object of class "numeric" containing the central tendency of the continuous character 'a', which helps define an ellipse's size. If left blank, the default behavior will draw a value from a uniform distribution from -5 to 10.
  #' @field kappa An object of class "numeric" containing the rate of mean reversion of the continuous character 'a', which helps define an ellipse's size. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 5.
  #' @field alpha An object of class "numeric" containing the log size of the new ellipse (daughter 2) under the sympatric cladogenetic scenario. Default value is -5.
  #' @field z An object of class "numeric" containing the height of the tilt vector for the ellipses. This model parameter does not evolve over time, and is only used in the interpretation of rates of 'r' and 's'. The default value will suffice for most analyses. The default value is 10.
ellipseParam <- R6::R6Class("ellipseParam",

private = list(

  # Initializing private fields (correspond to active bindings)
  # These fields (starting with '.') are for accessing information internally
  # Active bindings with corresponding names are for accessing information externally
  .cval=NULL,
  .hval=NULL,
  .rho_m=NULL,
  .rho_c=NULL,
  .rho_h=NULL,
  .sigma_r=NULL,
  .sigma_s=NULL,
  .sigma_a=NULL,
  .alpha=NULL,
  .root_r=NULL,
  .root_s=NULL,
  .root_a=NULL,
  .sigma_x=NULL,
  .sigma_y=NULL,
  .root_x=NULL,
  .root_y=NULL,
  .z=NULL,
  .mu=NULL,
  .kappa=NULL
),

active = list(

  # Active bindings that return values (stored in private fields) when called without assigning a value, e.g. ellipseParam$sigma_x
  # Can also be set by the user (after type checking), e.g. ellipseParam$sigma_x <- 1
  cval=function(value) {
    if (missing(value)) {private$.cval}
    else {
      if (!("numeric" %in% class(value))) {
        stop("cval must be numeric",call.=FALSE)
      }
      private$.cval = value
    }
  },

  hval=function(value) {
    if (missing(value)) {private$.hval}
    else {
      if (!("numeric" %in% class(value))) {
        stop("hval must be numeric",call.=FALSE)
      }
      private$.hval = value
    }
  },

  rho_m=function(value) {
    if (missing(value)) {private$.rho_m}
    else {
      if (!("numeric" %in% class(value))) {stop("rho_m must be numeric",call.=FALSE)}
      if (!(length(value)==2)) {stop("rho_m must be of length 2",call.=FALSE)}
      private$.rho_m = value
    }
  },

  rho_c=function(value) {
    if (missing(value)) {private$.rho_c}
    else {
      if (!("numeric" %in% class(value))) {stop("rho_c must be numeric",call.=FALSE)}
      if (!(length(value)==length(self$cval))) {stop("rho_c must be of the same length as cval",call.=FALSE)}
      private$.rho_c = value
    }
  },

  rho_h=function(value) {
    if (missing(value)) {private$.rho_h}
    else {
      if (!("numeric" %in% class(value))) {stop("rho_h must be numeric",call.=FALSE)}
      if (!(length(value)==length(self$hval))) {stop("rho_h must be of the same length as hval",call.=FALSE)}
      private$.rho_h = value
    }
  },

  sigma_x=function(value) {
    if (missing(value)) {private$.sigma_x}
    else {
      if (!("numeric" %in% class(value))) {stop("sigma_x must be numeric",call.=FALSE)}
      if (!(value > 0)) {stop("sigma_x must be positive",call.=FALSE)}
      private$.sigma_x = value
    }
  },

  sigma_y=function(value) {
    if (missing(value)) {private$.sigma_y}
    else {
      if (!("numeric" %in% class(value))) {stop("sigma_y must be numeric",call.=FALSE)}
      if (!(value > 0)) {stop("sigma_y must be positive",call.=FALSE)}
      private$.sigma_y = value
    }
  },

  sigma_r=function(value) {
    if (missing(value)) {private$.sigma_r}
    else {
      if (!("numeric" %in% class(value))) {stop("sigma_r must be numeric",call.=FALSE)}
      if (!(value > 0)) {stop("sigma_r must be positive",call.=FALSE)}
      private$.sigma_r = value
    }
  },

  sigma_s=function(value) {
    if (missing(value)) {private$.sigma_s}
    else {
      if (!("numeric" %in% class(value))) {stop("sigma_s must be numeric",call.=FALSE)}
      if (!(value > 0)) {stop("sigma_s must be positive",call.=FALSE)}
      private$.sigma_s = value
    }
  },

  sigma_a=function(value) {
    if (missing(value)) {private$.sigma_a}
    else {
      if (!("numeric" %in% class(value))) {stop("sigma_a must be numeric",call.=FALSE)}
      if (!(value > 0)) {stop("sigma_a must be positive",call.=FALSE)}
      private$.sigma_a = value
    }
  },

  root_x=function(value) {
    if (missing(value)) {private$.root_x}
    else {
      if (!("numeric" %in% class(value))) {
        stop("root_x must be numeric",call.=FALSE)
      }
      private$.root_x = value
    }
  },

  root_y=function(value) {
    if (missing(value)) {private$.root_y}
    else {
      if (!("numeric" %in% class(value))) {
        stop("root_y must be numeric",call.=FALSE)
      }
      private$.root_y = value
    }
  },

  root_r=function(value) {
    if (missing(value)) {private$.root_r}
    else {
      if (!("numeric" %in% class(value))) {
        stop("root_r must be numeric",call.=FALSE)
      }
      private$.root_r = value
    }
  },

  root_s=function(value) {
    if (missing(value)) {private$.root_s}
    else {
      if (!("numeric" %in% class(value))) {
        stop("root_s must be numeric",call.=FALSE)
      }
      private$.root_s = value
    }
  },

  root_a=function(value) {
    if (missing(value)) {private$.root_a}
    else {
      if (!("numeric" %in% class(value))) {
        stop("root_a must be numeric",call.=FALSE)
      }
      private$.root_a = value
    }
  },

  mu=function(value) {
    if (missing(value)) {private$.mu}
    else {
      if (!("numeric" %in% class(value))) {
        stop("mu must be numeric",call.=FALSE)
      }
      private$.mu = value
    }
  },

  kappa=function(value) {
    if (missing(value)) {private$.kappa}
    else {
      if (!("numeric" %in% class(value))) {stop("kappa must be numeric",call.=FALSE)}
      if (!(value > 0)) {stop("kappa must be positive",call.=FALSE)}
      private$.kappa = value
    }
  },

  alpha=function(value) {
    if (missing(value)) {private$.alpha}
    else {
      if (!("numeric" %in% class(value))) {
        stop("alpha must be numeric",call.=FALSE)
      }
      private$.alpha = value
    }
  },

  z=function(value) {
    if (missing(value)) {private$.z}
    else {
      if (!("numeric" %in% class(value))) {stop("z must be numeric",call.=FALSE)}
      if (!(value > 0)) {stop("z must be positive",call.=FALSE)}
      private$.z = value
    }
  }
),

public = list(

  #' @description Creates a new instance of the ellipseParam class
  #' @param cval An object of class "vector" containing any number of numeric elements indicating the values associated with different concentric circles relative to the size of the original circle. The default value includes 4 concentric circles: 0, 0.5, 1.0, and 1.5.
  #' @param hval An object of class "vector" containing any number of numeric elements indicating the values associated with different direction lines in radians. The default value includes 8 direction lines: 0, pi/4, pi/2, 3pi/4, pi, 5pi/4, 3pi/2, and 7pi/4.
  #' @param rho_m An object of class "vector" containing two numeric elements: the prior probability of cladogenetic mode 1 (budding) and the prior probability of cladogenetic mode 2 (splitting). The default value assigns probability 0.5 to both modes.
  #' @param rho_c An object of class "vector" containing any number of numeric elements (equal to the length of cval) corresponding to the prior probabilities of different concentric circle assignments. The default value assigns probability 0.25 to each of 4 concentric circles.
  #' @param rho_h An object of class "vector" containing any number of numeric elements (equal to the length of hval) corresponding to the prior probabilities of different direction line assignments. The default value assigns probability 0.125 to each of 8 direction lines.
  #' @param sigma_x An object of class "numeric" containing the evolutionary rate of the continuous character 'x', which helps define an ellipse's center location. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
  #' @param sigma_y An object of class "numeric" containing the evolutionary rate of the continuous character 'y', which helps define an ellipse's center location. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
  #' @param sigma_r An object of class "numeric" containing the evolutionary rate of the continuous character 'r', which helps define an ellipse's "oblongness". If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
  #' @param sigma_s An object of class "numeric" containing the evolutionary rate of the continuous character 's', which helps define an ellipse's "oblongness". If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
  #' @param sigma_a An object of class "numeric" containing the evolutionary rate of the continuous character 'a', which helps define an ellipse's size. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 5.
  #' @param root_x An object of class "numeric" containing the root value for the 'x' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and variance 1.
  #' @param root_y An object of class "numeric" containing the root value for the 'y' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and variance 1.
  #' @param root_r An object of class "numeric" containing the root value for the 'r' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and variance 5.
  #' @param root_s An object of class "numeric" containing the root value for the 's' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and variance 5.
  #' @param root_a An object of class "numeric" containing the root value for the 'a' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and variance 5.
  #' @param mu An object of class "numeric" containing the central tendency of the continuous character 'a', which helps define an ellipse's size. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and variance 5.
  #' @param kappa An object of class "numeric" containing the rate of mean reversion of the continuous character 'a', which helps define an ellipse's size. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 5.
  #' @param alpha An object of class "numeric" containing the log size of the new ellipse (daughter 2) under the sympatric cladogenetic scenario. Default value is -5.
  #' @param z An object of class "numeric" containing the height of the tilt vector for the ellipses. This model parameter does not evolve over time, and is only used in the interpretation of rates of 'r' and 's'. The default value will suffice for most analyses. The default value is 10.
  initialize = function(cval=NULL,
                        hval=NULL,
                        rho_m=NULL,
                        rho_c=NULL,
                        rho_h=NULL,
                        sigma_r=NULL,
                        sigma_s=NULL,
                        sigma_a=NULL,
                        alpha=NULL,
                        root_r=NULL,
                        root_s=NULL,
                        root_a=NULL,
                        sigma_x=NULL,
                        sigma_y=NULL,
                        root_x=NULL,
                        root_y=NULL,
                        z=NULL,
                        mu=NULL,
                        kappa=NULL) {

    # default values
    # some are fixed (like the vectors of values and base probabilities)
    # some are sampled (root and stochastic process parameters)
    # the distributions for sampled values are also the default prior distributions for analyses
    if (is.null(cval)) {cval=seq(0,1.5,.5)}
    if (is.null(hval)) {hval=2*pi*seq(0,7/8,1/8)}
    if (is.null(rho_m)) {rho_m=c(.5,.5)}
    if (is.null(rho_c)) {rho_c=rep(1/4,4)}
    if (is.null(rho_h)) {rho_h=rep(1/8,8)}
    if (is.null(sigma_x)) {sigma_x=runif(1,0,1)}
    if (is.null(sigma_y)) {sigma_y=runif(1,0,1)}
    if (is.null(sigma_r)) {sigma_r=runif(1,0,1)}
    if (is.null(sigma_s)) {sigma_s=runif(1,0,1)}
    if (is.null(sigma_a)) {sigma_a=runif(1,0,5)}
    if (is.null(root_r)) {root_r=rnorm(1,0,5)}
    if (is.null(root_s)) {root_s=rnorm(1,0,5)}
    if (is.null(root_a)) {root_a=runif(1,-5,10)}
    if (is.null(root_x)) {root_x=rnorm(1,0,1)}
    if (is.null(root_y)) {root_y=rnorm(1,0,1)}
    if (is.null(alpha)) {alpha=-5}
    if (is.null(z)) {z=10}
    if (is.null(mu)) {mu=runif(1,-5,10)}
    if (is.null(kappa)) {kappa=runif(1,0,5)}

    # setting fields
    # will use user-input values if provided
    # otherwise, will use default values
    self$cval = cval
    self$hval = hval
    self$rho_m = rho_m
    self$rho_c = rho_c
    self$rho_h = rho_h
    self$sigma_r = sigma_r
    self$sigma_s = sigma_s
    self$sigma_a = sigma_a
    self$root_r = root_r
    self$root_s = root_s
    self$root_a = root_a
    self$sigma_x = sigma_x
    self$sigma_y = sigma_y
    self$root_x = root_x
    self$root_y = root_y
    self$alpha = alpha
    self$z = z
    self$mu = mu
    self$kappa = kappa
  }
))

# ellipseParam functions --------------------------------------------------

#' @title sim_params() function
#' @description Creates a new instance of the ellipseParam class -- this function is a wrapper for ellipseParam$new()
#' @param cval An object of class "vector" containing any number of numeric elements indicating the values associated with different concentric circles relative to the size of the original circle. The default value includes 4 concentric circles: 0, 0.5, 1.0, and 1.5.
#' @param hval An object of class "vector" containing any number of numeric elements indicating the values associated with different direction lines in radians. The default value includes 8 direction lines: 0, pi/4, pi/2, 3pi/4, pi, 5pi/4, 3pi/2, and 7pi/4.
#' @param rho_m An object of class "vector" containing two numeric elements: the prior probability of cladogenetic mode 1 (budding) and the prior probability of cladogenetic mode 2 (splitting). The default value assigns probability 0.5 to both modes.
#' @param rho_c An object of class "vector" containing any number of numeric elements (equal to the length of cval) corresponding to the prior probabilities of different concentric circle assignments. The default value assigns probability 0.25 to each of 4 concentric circles.
#' @param rho_h An object of class "vector" containing any number of numeric elements (equal to the length of hval) corresponding to the prior probabilities of different direction line assignments. The default value assigns probability 0.125 to each of 8 direction lines.
#' @param sigma_x An object of class "numeric" containing the evolutionary rate of the continuous character 'x', which helps define an ellipse's center location. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
#' @param sigma_y An object of class "numeric" containing the evolutionary rate of the continuous character 'y', which helps define an ellipse's center location. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
#' @param sigma_r An object of class "numeric" containing the evolutionary rate of the continuous character 'r', which helps define an ellipse's "oblongness". If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
#' @param sigma_s An object of class "numeric" containing the evolutionary rate of the continuous character 's', which helps define an ellipse's "oblongness". If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
#' @param sigma_a An object of class "numeric" containing the evolutionary rate of the continuous character 'a', which helps define an ellipse's size. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
#' @param root_x An object of class "numeric" containing the root value for the 'x' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and standard deviation 1.
#' @param root_y An object of class "numeric" containing the root value for the 'y' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and standard deviation 1.
#' @param root_r An object of class "numeric" containing the root value for the 'r' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and standard deviation 1.
#' @param root_s An object of class "numeric" containing the root value for the 's' character. If left blank, the default behavior will draw a value from a normal distribution with mean 0 and standard deviation 1.
#' @param root_a An object of class "numeric" containing the root value for the 'a' character. If left blank, the default behavior will draw a value from a uniform distribution from -5 to 10.
#' @param mu An object of class "numeric" containing the central tendancy of the continuous character 'a', which helps define an ellipse's size. If left blank, the default behavior will draw a value from a uniform distribution from -5 to 10.
#' @param kappa An object of class "numeric" containing the rate of mean reversion of the continuous character 'a', which helps define an ellipse's size. If left blank, the default behavior will draw a value from a uniform distribution from 0 to 1.
#' @param alpha An object of class "numeric" containing the log size of the new ellipse (daughter 2) under the sympatric cladogenetic scenario. Default value is -5.
#' @param z An object of class "numeric" containing the height of the tilt vector for the ellipses. This model parameter does not evolve over time, and is only used in the interpretation of rates of 'r' and 's'. The default value will suffice for most analyses. The default value is 10.
#' @export
sim_params = function(cval=NULL,
                      hval=NULL,
                      rho_m=NULL,
                      rho_c=NULL,
                      rho_h=NULL,
                      sigma_r=NULL,
                      sigma_s=NULL,
                      sigma_a=NULL,
                      alpha=NULL,
                      root_r=NULL,
                      root_s=NULL,
                      root_a=NULL,
                      sigma_x=NULL,
                      sigma_y=NULL,
                      root_x=NULL,
                      root_y=NULL,
                      z=NULL,
                      mu=NULL,
                      kappa=NULL) {
  # This function is just a wrapper for ellipseParam$new()
  params = ellipseParam$new(cval, hval, rho_m, rho_c, rho_h, sigma_r, sigma_s, sigma_a, alpha, root_r, root_s, root_a, sigma_x, sigma_y, root_x, root_y, z, mu, kappa)
  return(params)
}

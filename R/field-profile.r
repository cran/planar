##' field profile outside of a ML stack
##'
##' returns the field outside of the structure using Fresnel coefficients and transfer matrix for comparison
##' @title field.outside
##' @export
##' @param d distance outside in nm
##' @param theta angle
##' @param epsilon list of permittivities
##' @param thickness vector of layer thickness in nm
##' @param wavelength wavelength in nm
##' @param polarisation polarisation
##' @param ... further args passed to multilayer 
##' @return long format data.frame with field on both sides of the multilayer as a function of angle
##' @author baptiste Auguie
field.outside <- function(d=1, theta,
                          epsilon,
                          thickness,
                          wavelength, polarisation="p", ...){

 res <- multilayer(lambda=wavelength, theta=theta, epsilon=epsilon, d=rep(d,length(thickness)),
                   thickness=thickness, polarisation=polarisation, ...)

 n <- length(epsilon)
 sinl <- sin(theta)
 sinr <- sqrt(epsilon[[1]]) / sqrt(epsilon[[n]]) * sinl

 k0 <- 2*pi/wavelength
 kx <- k0*sqrt(epsilon[[1]]) * sinl
 kzl <- sqrt(epsilon[[1]] * k0^2 - kx^2+0i)
 kzr <- sqrt(epsilon[[n]] * k0^2 - kx^2+0i)
 
 pol.fac1 <- if(polarisation == "p") epsilon[[1]] / epsilon[[n]]  else 1
 pol.fac2 <- if(polarisation == "p") sinl^2  else 1
 pol.fac3 <- if(polarisation == "p") sinr^2  else 1
 
 left <- pol.fac2 * Mod(exp(-1i*d* kzl) + res$reflection * exp(1i*d* kzl))^2
 right <- pol.fac3 * Mod(res$transmission * exp(1i*(d)* kzr))^2
 
 left2 <- if(polarisation == "p") res$Ml.perp[[1]][,1] else res$Ml.par[[1]][,1]
 right2 <- if(polarisation == "p") res$Mr.perp[[2]][,2]  else res$Mr.par[[2]][,2]
 
 d <- data.frame(theta=theta*180/pi,
                 left=left, right=right*pol.fac1,
                 left2 = left2, right2 = right2)

 classify(d, id="theta", vars=list(side = rep(c("left", "right"), 2),
                           model = rep(c("matrix","fresnel"), each=2)))
 
}


##' field profile in a ML stack
##'
##' runs multilayer and returns the LFEF as a function of distance inside and outside of the structure
##' @title field.profile
##' @export
##' @param lambda wavelength
##' @param theta angle
##' @param dmax maximum distance to interface, if > layer thickness
##' @param thickness vector of layer thickness
##' @param res resolution of sampling points
##' @param epsilon list of permittivities
##' @param polarisation polarisation
##' @param ... further args passed to multilayer 
##' @return long format data.frame with positions and LFEF (para and perp)
##' @author baptiste Auguie
##' @references
##' Principles of surface-enhanced Raman spectroscopy and related plasmonic effects
##' Eric C. Le Ru and Pablo G. Etchegoin, published by Elsevier, Amsterdam (2009).
field.profile <- function(lambda=500, theta=0, polarisation='p',
                          thickness = c(0, 20, 140, 20, 0), dmax=200,  res=1e3,
                          epsilon=list(1^2, -12 + (0 + (0+1i)), 1.38^2, -12 + (0 + (0+1i)), 1.46^2), ...){

  d <- seq(0, max(c(dmax,thickness)), length=res)
  
  res <- multilayer(lambda=lambda, theta=theta,
                    epsilon=epsilon,
                    thickness = thickness, d=d,
                    polarisation=polarisation, ...)
  
all <- lapply(seq(1,length(res$dist) -1), function(lay){
  data.frame(x = res$dist[[lay+1]], M.par=res$Mr.par[[lay]],
             M.perp=res$Mr.perp[[lay]])
})

all <- c(list(data.frame(x = res$dist[[1]],
                              M.par=res$Ml.par[[1]],
                              M.perp=res$Ml.perp[[1]])), all)

names(all) <- paste("layer", seq_along(res$dist))
melt(all, id=1)
}

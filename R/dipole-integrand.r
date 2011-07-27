
##' Dipole decay rates near a multilayer interface
##'
##' Integrand of the dipole decay rates near a multilayer interface
##' @title dipole.integrand
##' @export
##' @param d distance in nm
##' @param q normalised in-plane wavevector in [0, infty) OR sqrt(1-q^2) if change.variable
##' @param lambda wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @author baptiste Auguie
dipole.integrand <- function(d=10, q, lambda,
                             epsilon = list(incident=1.5^2, 1.0^2),
                             thickness = c(0, 0)){
  
  ## define constants
  k0 <- 2*pi/lambda
  k1 <- sqrt(epsilon[[1]])*k0

  Nlambda <- length(k0)
  Nq <- length(q)
  
  u <- sqrt(1 - q^2 + 0i)
  
  rp <- recursive.fresnel2(lambda=lambda,
                           q = q,
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="p")$reflection
  
  rs <- recursive.fresnel2(lambda=lambda,
                           q = q,
                           epsilon=epsilon,
                           thickness=thickness,
                           polarisation="s")$reflection
  
  phase <- exp(2i*d*outer(k1,u))
  
  integrand.p <- Re(matrix(q^3 / u, Nlambda, Nq, byrow=TRUE) * rp * phase)
  integrand.s <- Re( (rs / matrix(u, Nlambda, Nq, byrow=TRUE) -
                      rp * matrix(u, Nlambda, Nq, byrow=TRUE)) *
                    matrix(q, Nlambda, Nq, byrow=TRUE) * phase)
   
  list(integrand.p = integrand.p, integrand.s = integrand.s)
}


##' Dipole decay rates near a multilayer interface
##'
##' dipole decay rates near a multilayer interface
##' @title dipole.direct
##' @export
##' @param d distance in nm
##' @param lambda wavelength in nm
##' @param epsilon list of dielectric functions
##' @param thickness list of layer thicknesses
##' @param Nquadrature1 quadrature points in radiative region
##' @param Nquadrature2 quadrature points in SPPs region
##' @param Nquadrature3 quadrature points in dipole image region
##' @param qcut transition between regions 2 and 3
##' @param qmax maximum q of region 3
##' @author baptiste Auguie
dipole.direct <- function(d=1,
                   lambda ,
                   epsilon = list(incident=1.0^2),
                   thickness = c(0, 0),
                   Nquadrature1 = 50, Nquadrature2 = 200, Nquadrature3 = 50,
                   qcut = NULL, qmax = Inf){
   
  require(statmod) # quadrature points in (-1, 1)

  GL1 <- gauss.quad(Nquadrature1)
  GL2 <- gauss.quad(Nquadrature2)
  GL3 <- gauss.quad(Nquadrature3)

  Nq1 <- length(GL1$nodes)
  Nq2 <- length(GL2$nodes)
  Nq3 <- length(GL3$nodes)
  
  Nlambda <- length(lambda)
  
  ## if no qcut provided, estimate one from max of
  ## all possible SPP dispersions
  if(is.null(qcut)){
    qcut <- 1.1

    epsilon_norm <- do.call(cbind, epsilon)
    
    for(ii in seq(1, length(epsilon) - 1)){
      qspp <- sqrt(epsilon_norm[,ii] / epsilon_norm[,1])*
        sqrt(epsilon_norm[,ii+1] / (epsilon_norm[,ii] + epsilon_norm[,ii+1]))
      
      qcut <- max(qcut, max(Re(qspp)))
    }
    
    print(paste("using qcut=",round(qcut,2)))
    
  }

  ## integration from 0 to 1
  qmax1 <- 1; qmin1 <- 0;
  C1 <- (qmax1 - qmin1)/2 ; D1 <- (qmax1+qmin1)/2
  qnodes1 <- C1 * GL1$nodes + D1
  qweights1 <- GL1$weights * C1
  
  in1 <- dipole.integrand(q=qnodes1,
                          d=d, lambda=lambda,
                          epsilon=epsilon, thickness=thickness)
      
  weights1 <- matrix(qweights1, nrow=Nlambda, ncol=Nq1, byrow=TRUE)

  integral1.perp <- rowSums(in1$integrand.p*weights1)
  integral1.par <- rowSums(in1$integrand.s*weights1)
  
  
  ## integration from 1 to qcut
  qmax2 <- qcut; qmin2 <- 1;
  C2 <- (qmax2 - qmin2)/2 ; D2 <- (qmax2+qmin2)/2
  qnodes2 <- C2 * GL2$nodes + D2
  qweights2 <- GL2$weights * C2
  
  in2 <- dipole.integrand(q=qnodes2,
                                  d=d, lambda=lambda,
                                  epsilon=epsilon, thickness=thickness)
  
  weights2 <- matrix(qweights2, nrow=Nlambda, ncol=Nq2, byrow=TRUE)
  
  integral2.perp <- rowSums(in2$integrand.p*weights2)
  integral2.par <- rowSums(in2$integrand.s*weights2)
  

  ## integration from qcut to qmax
  if(is.finite(qmax)){
    ## straight integration from qcut to qmax
    qmax3 <- qmax; qmin3 <- qcut;
    C3 <- (qmax3 - qmin3)/2 ; D3 <- (qmax3+qmin3)/2
    
    qnodes3 <- C3 * GL3$nodes + D3
    qweights3 <- GL3$weights * C3
  
  } else {
    print("performing a change of variable mapping [qcut, infty) -> [0,1]")
    ## change of variables
    ## \int_a^\infty f(x)dx = \int_0^1 f(a + t/(1-t)). 1 / (1-t)^2 dt
    ## as suggested on http://ab-initio.mit.edu/wiki/index.php/Cubature
    qmax3 <- 1; qmin3 <- 0;
    C3 <- (qmax3 - qmin3)/2 ; D3 <- (qmax3+qmin3)/2
    
    qnodes3 <- C3 * GL3$nodes + D3
    qweights3 <- GL3$weights * C3 * 1 / (1 - qnodes3)^2
    qnodes3 <- qcut + qnodes3 / (1 - qnodes3)
    
  }
  
  in3 <- dipole.integrand(q=qnodes3,
                          d=d, lambda=lambda,
                          epsilon=epsilon, thickness=thickness)
      
  weights3 <- matrix(qweights3, nrow=Nlambda, ncol=Nq3, byrow=TRUE)
  
  integral3.perp <- rowSums(in3$integrand.p*weights3)
  integral3.par <- rowSums(in3$integrand.s*weights3)
  
  ## data.frame(wavelength=lambda, Mtot.par= 1 + 3/2*(integral1 + integral2 + integral3))
  
  data.frame(wavelength=lambda,
             Mtot.perp = 1 + 3/2*(integral1.perp + integral2.perp + integral3.perp),
             Mtot.par = 1 + 3/4*(integral1.par + integral2.par + integral3.par) )
  
}

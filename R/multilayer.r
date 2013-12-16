##' Multilayer Fresnel coefficients
##'
##' solves the EM problem of a multilayered interface
##' @title multilayer
##' @export
##' @param wavelength [vector] wavelength in nm
##' @param k0 [vector] wavevector in nm^-1
##' @param angle [vector] incident angles in radians
##' @param q [vector] normalised incident in-plane wavevector
##' @param epsilon list of N+2 dielectric functions, each of length 1 or length(wavelength)
##' @param thickness vector of N+2 layer thicknesses, first and last are dummy
##' @param d vector of distances where LFIEF are evaluated from each interface
##' @param dout vector of distances where LFIEF are evaluated outside the stack
##' @param polarisation [character] switch between p- and s- polarisation
##' @return fresnel coefficients and field profiles
##' @param ... unused
##' @author baptiste Auguie
##' @references
##' Principles of surface-enhanced Raman spectroscopy and related plasmonic effects.
##' Eric C. Le Ru and Pablo G. Etchegoin, published by Elsevier, Amsterdam (2009).
multilayer <- function(wavelength = 2*pi/k0, k0 = 2*pi/wavelength,
                       angle = asin(q), q = sin(angle),
                       epsilon = list(incident=1.5^2, 1.33),
                       thickness = c(0, 0), 
                       polarisation = c('p', 's'),
                       d = 1, dout=d, ...){

  ## checks
  stopifnot(thickness[1]==0L, thickness[length(thickness)]==0L)
  polarisation <- match.arg(polarisation)
  
  epsilon = do.call(cbind, epsilon)
  ## case pure scalars
  if(nrow(epsilon) == 1L)
    epsilon <- matrix(epsilon, nrow=length(k0), ncol=length(thickness), byrow=TRUE)
  
  ## define constants
  Nlambda <- length(k0)
  Nq <- length(q)
  Nlayer <- length(thickness)
  k02 <- k0^2
  kx <- outer(k0*sqrt(epsilon[,1] + 0i), q) # kx = q*k0
  kx2 <- kx^2
  
  ## loop to calculate kiz 
  kiz <- array(0 + 0i, dim=c(Nlambda, Nq, Nlayer))
  
  for (ii in seq(1, Nlayer)){
    kiz[ , , ii] <- sqrt(matrix(epsilon[,ii]*k02, nrow=Nlambda, ncol=Nq) - kx2 + 0i)
  }
  
  ## calculate the transition matrix M
  M11 <- M22 <- 1 + 0i
  M21 <- M12 <- 0 + 0i
  
  Mi11 <- Mi12 <- Mi21 <- Mi22 <- array(1 + 0i, dim=c(Nlambda, Nq, Nlayer-1))
  
  for (ii in seq(1, Nlayer-1)){
    
   if(polarisation == 'p'){
     Ki <- matrix(epsilon[,ii] / epsilon[,ii+1] + 0i , nrow=Nlambda, ncol=Nq) *
       kiz[,,ii+1] / kiz[,,ii]
   } else { # s-polarisation
      Ki <- kiz[,,ii+1] / kiz[,,ii]
   }
   
   phasei <- exp(1i*thickness[ii]*kiz[,,ii])
   
   Mi11[,,ii] <- 0.5*(1+Ki) / phasei
   Mi21[,,ii] <- 0.5*(1-Ki) * phasei
   Mi12[,,ii] <- 0.5*(1-Ki) / phasei
   Mi22[,,ii] <- 0.5*(1+Ki) * phasei
   
   M11new <- M11*Mi11[,,ii] + M12*Mi21[,,ii]
   M21new <- M21*Mi11[,,ii] + M22*Mi21[,,ii]
   M12new <- M11*Mi12[,,ii] + M12*Mi22[,,ii]
   M22new <- M21*Mi12[,,ii] + M22*Mi22[,,ii]
   
   M11 <- M11new
   M12 <- M12new
   M21 <- M21new
   M22 <- M22new
   
  }

  ## calculate the Fresnel coefficients
  transmission <- 1 / M11
  reflection <- M21 * transmission
  
  ## make a list of physical locations to sample the fields
  
  Nd <- length(d)
  sampling <- lapply(thickness, function(t) {
    d[d<=t]
  })
  sampling[[1]] <- sampling[[length(sampling)]] <- dout # outside of the stack
    
  ## enhancement factors, default to 0
  sizes <- lapply(sampling, length)
  
  Ml.perp <- Ml.par <- 
    lapply(sizes[-length(sizes)], function(s) array(0, dim=c(Nlambda, Nq, s)))
  
  Mr.perp <-  Mr.par <-
    lapply(sizes[-1], function(s) array(0, dim=c(Nlambda, Nq, s)))
    
  ## absolute positions, first interface is at 0
  interfaces <- cumsum(thickness)[-length(thickness)]
    
  distance <- list()
  distance[[1]] <- -dout
  
  
  #####################
  ## TM polarisation ##
  #####################
  
  if(polarisation == 'p'){
    
    ## field amplitudes p is for prime, reference is field in region 1
    Hiy.H1y <- Hpiy.H1y <- Eix.E1 <- Epix.E1 <- Eiz.E1 <- Epiz.E1 <-
      array(0+0i, dim=c(Nlambda, Nq, Nlayer))
    
    Hiy.H1y[,,Nlayer] <- transmission
    Hpiy.H1y[,,Nlayer] <- 0i
    
    AuxE1 <- matrix(sqrt(epsilon[,1] + 0i) / k0 / epsilon[,Nlayer], nrow=Nlambda, ncol=Nq)
    Eix.E1[,,Nlayer] <- Hiy.H1y[,,Nlayer] * kiz[,,Nlayer] * AuxE1
    Epix.E1[,,Nlayer] <- 0i
    AuxE2 <- outer(epsilon[,1] / epsilon[,Nlayer], Re(q))
    Eiz.E1[,,Nlayer] <- - Hiy.H1y[,,Nlayer] * AuxE2
    Epiz.E1[,,Nlayer] <- 0i

    ## loop downwards to compute all field amplitudes
    for (ii in seq(Nlayer-1, 1, by=-1)){
       Hiy.H1y[,,ii] <- Mi11[,,ii]*Hiy.H1y[,,ii+1] + Mi12[,,ii]*Hpiy.H1y[,,ii+1]
       Hpiy.H1y[,,ii] <- Mi21[,,ii]*Hiy.H1y[,,ii+1] + Mi22[,,ii]*Hpiy.H1y[,,ii+1]
       
       AuxE1 <- matrix(sqrt(epsilon[,1] + 0i) / k0 / epsilon[,ii]  , nrow=Nlambda, ncol=Nq)
       Eix.E1[,,ii] <- Hiy.H1y[,,ii] * kiz[,,ii] * AuxE1
       Epix.E1[,,ii] <- - Hpiy.H1y[,,ii] * kiz[,,ii] * AuxE1
       AuxE2 <- outer(epsilon[,1] / epsilon[,ii], Re(q))
       Eiz.E1[,,ii] <- - Hiy.H1y[,,ii] * AuxE2
       Epiz.E1[,,ii] <- - Hpiy.H1y[,,ii] * AuxE2
    }
    
    ## loop to compute the local field EFs
    for (ii in seq(1, Nlayer-1, by=1)){

      ## left of interface ii
      ## the relative coordinate is sampling[[ii]] - thickness[ii]
      d1 <- thickness[ii] - sampling[[ii]]
      
      Ml.perp[[ii]] <- sapply(d1, function(.d)
                                Mod(Eiz.E1[,,ii]  * exp(1i*.d*kiz[,,ii]) +
                                    Epiz.E1[,,ii] * exp(-1i*.d*kiz[,,ii]))^2,
                                  simplify="array")
      
      Ml.par[[ii]] <- sapply(d1, function(.d)
                               Mod(Eix.E1[,,ii]  * exp(1i*.d*kiz[,,ii]) +
                                   Epix.E1[,,ii] * exp(-1i*.d*kiz[,,ii]))^2,
                               simplify="array")
      
      ## right of interface ii
      ## the relative coordinate is sampling[[ii+1]]
      
      Mr.perp[[ii]] <- sapply(sampling[[ii+1]], function(.d)
                               Mod(Eiz.E1[,,ii+1]  * exp( 1i*.d*kiz[,,ii+1]) +
                                   Epiz.E1[,,ii+1] * exp(-1i*.d*kiz[,,ii+1]))^2,
                               simplify="array")
      Mr.par[[ii]] <- sapply(sampling[[ii+1]], function(.d)
                              Mod(Eix.E1[,,ii+1]  * exp( 1i*.d*kiz[,,ii+1]) +
                                  Epix.E1[,,ii+1] * exp(-1i*.d*kiz[,,ii+1]))^2,
                              simplify="array")
      
      ## absolute positions on the right of interface ii
      distance[[ii+1]] <-    sampling[[ii+1]] + interfaces[ii]

    } # end loop 
    fields <- list(Eix.E1=Eix.E1, Epix.E1=Epix.E1,
                   Eiz.E1=Eiz.E1, Epiz.E1=Epiz.E1)
  }

  if(polarisation =="s"){
    
    ## ############### ##
    ## TE polarisation ##
    ## ############### ##
    
    ## field amplitudes p is for prime, reference is field in region 1
    Eiy.E1y <- Epiy.E1y <- Hix.H1 <- Hpix.H1 <- Hiz.H1 <- Hpiz.H1 <-
      array(0+0i, dim=c(Nlambda, Nq, Nlayer))
    
    Eiy.E1y[,,Nlayer] <- transmission
    Epiy.E1y[,,Nlayer] <- 0i
    
    AuxH1 <- matrix(1 / (k0 * sqrt(epsilon[,1] + 0i)), nrow=Nlambda, ncol=Nq)
    Hix.H1[,,Nlayer] <- -Eiy.E1y[,,Nlayer] * kiz[,,Nlayer] * AuxH1
    Hpix.H1[,,Nlayer] <- 0i
    AuxH2 <- matrix(Re(q), nrow=Nlambda, ncol=Nq)
    Hiz.H1[,,Nlayer] <- Eiy.E1y[,,Nlayer] * AuxH2
    Hpiz.H1[,,Nlayer] <- 0i

    ## loop downwards to compute all field amplitudes
    for (ii in seq(Nlayer-1, 1, by=-1)){
       Eiy.E1y[,,ii] <- Mi11[,,ii]*Eiy.E1y[,,ii+1] + Mi12[,,ii]*Epiy.E1y[,,ii+1]
       Epiy.E1y[,,ii] <- Mi21[,,ii]*Eiy.E1y[,,ii+1] + Mi22[,,ii]*Epiy.E1y[,,ii+1]
       
       Hix.H1[,,ii]  <- - Eiy.E1y[,,ii] * kiz[,,ii] * AuxH1
       Hpix.H1[,,ii] <-   Epiy.E1y[,,ii] * kiz[,,ii] * AuxH1
       Hiz.H1[,,ii]  <-   Eiy.E1y[,,ii] * AuxH2
       Hpiz.H1[,,ii] <-   Epiy.E1y[,,ii] * AuxH2
    }

    ## loop to compute the local field EFs
    for (ii in seq(1, Nlayer-1, by=1)){
      
      ## left of interface ii
      ## the relative coordinate is sampling[[ii]] - thickness[ii]
      d1 <- thickness[ii] - sampling[[ii]]

      Ml.par[[ii]] <- sapply(d1, function(.d)
                               Mod(Eiy.E1y[,,ii]  * exp(1i*.d*kiz[,,ii]) +
                                   Epiy.E1y[,,ii] * exp(-1i*.d*kiz[,,ii]))^2,
                               simplify="array")
      
      ## right of interface ii
      ## the relative coordinate is sampling[[ii+1]]
     
      Mr.par[[ii]] <- sapply(sampling[[ii+1]], function(.d)
                             Mod(Eiy.E1y[,,ii+1]  * exp( 1i*.d*kiz[,,ii+1]) +
                                 Epiy.E1y[,,ii+1] * exp(-1i*.d*kiz[,,ii+1]))^2,
                             simplify="array")

      ## Ml.perp[[ii]] <- 0 * Ml.par[[ii]] # default value
      ## Mr.perp[[ii]] <- 0 * Mr.par[[ii]] # default value
      
      ## absolute positions on the right of interface ii
      distance[[ii+1]] <- sampling[[ii+1]] + interfaces[ii]

    } # end loop 
 
    fields <- list(Eiy.E1y=Eiy.E1y, Epiy.E1y=Epiy.E1y)
  }  # end swich polarisation
  
  ## T is nt*cos(Ot)*|Et|^2 / ni*cos(Oi)*|Ei|^2
  ## for s-pol, |Et|^2 / |Ei|^2 = |ts|^2, hence T = nt/ni * cos(Ot)/cos(Oi) * |ts|^2
  ## for p-pol, |Et|^2 / |Ei|^2 = (ni/nt)^2 * |tp|^2, hence T = ni/nt * cos(Ot)/cos(Oi) * |tp|^2
  
  # ratio of refractive indices
  index.ratio <- matrix(Re(sqrt(epsilon[,1])/sqrt(epsilon[,Nlayer])), nrow=Nlambda, ncol=Nq)
  # ratio of cosines
  qq <- matrix(q, nrow=Nlambda, ncol=Nq, byrow=TRUE)
  m <- Re(sqrt(1 - (index.ratio * qq)^2 + 0i)/sqrt(1 - qq^2 + 0i))
  
  if(polarisation == "p"){
    rho <- index.ratio
  } else {
    rho <- 1 / index.ratio
  }
  dim(transmission) <- dim(m) # case 1-dims were dropped
  R <- Mod(reflection)^2
  T <- rho * m * Mod(transmission)^2
  T <- drop(T)
  transmission <- drop(transmission)
  
  ## results
  list(wavelength=wavelength, k0 = k0, 
       angle=angle, q=q, 
       reflection=reflection, transmission=transmission,
       R=R, T=T, A = 1 - R - T,
       dist=distance, fields = fields,
       Ml.perp=lapply(Ml.perp, drop), Ml.par=lapply(Ml.par, drop),
       Mr.perp=lapply(Mr.perp, drop), Mr.par=lapply(Mr.par, drop))
}

##' field profile in a ML stack
##'
##' runs multilayer and returns the LFIEF as a function of distance inside and outside of the structure
##' @title field_profile
##' @export
##' @param wavelength wavelength
##' @param angle angle
##' @param dmax maximum distance to interface, if > layer thickness
##' @param thickness vector of layer thickness
##' @param res resolution of sampling points
##' @param res2 resolution of sampling points outside stack
##' @param epsilon list of permittivities
##' @param polarisation polarisation
##' @param displacement logical, Mperp corresponds to displacement squared (D=epsilon x E)
##' @param ... further args passed to multilayer 
##' @return long format data.frame with positions and LFEF (para and perp)
##' @author baptiste Auguie
##' @family helping_functions
##' @references
##' Principles of surface-enhanced Raman spectroscopy and related plasmonic effects
##' 
##' Eric C. Le Ru and Pablo G. Etchegoin, published by Elsevier, Amsterdam (2009).
field_profile <- function(wavelength=500, angle=0, polarisation='p',
                          thickness = c(0, 20, 140, 20, 0), 
                          dmax=200,  res=1e3, res2=res/10,
                          epsilon=list(1^2, -12 , 1.38^2, -12 , 1.46^2), 
                          displacement=FALSE, ...){
  
  d <- seq(0, max(thickness), length=res)
  dout <- seq(0, dmax, length=res2)
  res <- multilayer(wavelength=wavelength, angle=angle,
                    epsilon=epsilon,
                    thickness = thickness, d=d, dout=dout,
                    polarisation=polarisation, ...)
  
  #   loop over 2nd to last interface
  all <- lapply(seq(1,length(res$dist) - 1), function(lay){
    Mperp <- if(displacement) (Mod(epsilon[[lay+1]]))^2 * res$Mr.perp[[lay]] else res$Mr.perp[[lay]]
    data.frame(x = res$dist[[lay+1]], M.par=res$Mr.par[[lay]],
               M.perp=Mperp)
  })
  #   combine with first interface
  Mperp <- if(displacement) (Mod(epsilon[[1]]))^2 * res$Ml.perp[[1]] else res$Ml.perp[[1]]
  all <- c(list(data.frame(x = res$dist[[1]],
                           M.par=res$Ml.par[[1]],
                           M.perp=Mperp)), all)
  
  names(all) <- paste("layer", seq_along(res$dist))
  melt(all, id=1)
}

##' Multilayer Fresnel coefficients
##'
##' solves the EM problem of a multilayered interface
##' @title multilayercpp
##' @export
##' @param wavelength [vector] wavelength in nm
##' @param k0 [vector] wavevector in nm^-1
##' @param angle [vector] incident angles in radians
##' @param q [vector] normalised incident in-plane wavevector
##' @param epsilon list of N+2 dielectric functions, each of length 1 or length(wavelength)
##' @param thickness vector of N+2 layer thicknesses, first and last are dummy
##' @param polarisation [character] switch between p- and s- polarisation
##' @param ... unused
##' @return fresnel coefficients and field profiles
##' @author baptiste Auguie
##' @examples
##' library(planar)
##' demo(package="planar")
multilayercpp <- function(wavelength = 2*pi/k0, k0 = 2*pi/wavelength,
                       angle = asin(q), q = sin(angle),
                       epsilon = list(incident=1.5^2, 1.33),
                       thickness = c(0, 0),
                       polarisation = c('p', 's'), ...){

  
  kx <- outer(k0*sqrt(epsilon[[1]]), q) # kx = q*k0
  epsilon = do.call(cbind, epsilon)
  polarisation = if(polarisation == "p") 0L else 1L
  
  Nlayer <- length(thickness)
  Nlambda <- length(k0)
  Nq <- length(q)
  
  ## checks
  stopifnot(thickness[1]==0L,
            thickness[Nlayer]==0L)
  
  stopifnot(Nlayer == ncol(epsilon),
            nrow(epsilon) == length(k0),
            nrow(kx) == length(k0),
            ncol(kx) == length(q))

  ## call the C++ function
  res <- planar$multilayer(as.vector(k0), as.matrix(kx), as.matrix(epsilon),
                           as.vector(thickness), as.integer(polarisation))
  
  transmission <- drop(res$transmission)
  reflection <- drop(res$reflection)
  
  ## T is nt*cos(Ot)*|Et|^2 / ni*cos(Oi)*|Ei|^2
  ## for s-pol, |Et|^2 / |Ei|^2 = |ts|^2, hence T = nt/ni * cos(Ot)/cos(Oi) * |ts|^2
  ## for p-pol, |Et|^2 / |Ei|^2 = (ni/nt)^2 * |tp|^2, hence T = ni/nt * cos(Ot)/cos(Oi) * |tp|^2
  
  
  # ratio of refractive indices
  index.ratio <- matrix(Re(sqrt(epsilon[,1])/sqrt(epsilon[,Nlayer])), nrow=Nlambda, ncol=Nq)
  # ratio of cosines
  qq <- matrix(q, nrow=Nlambda, ncol=Nq, byrow=TRUE)
  m <- Re(sqrt(1 - (index.ratio * qq)^2 + 0i)/sqrt(1 - qq^2 + 0i))
  
  if(polarisation == 0L){
    rho <- index.ratio
  } else {
    rho <- 1 / index.ratio
  }
  
  dim(transmission) <- dim(m) # case 1-dims were dropped
  R <- Mod(reflection)^2
  T <- rho * m * Mod(transmission)^2
  T <- drop(T)
  transmission <- drop(transmission)
  
  A <- 1 - R - T
  
  list(wavelength=wavelength, k0 = k0, 
       angle=angle, q=q, 
       reflection=reflection, transmission=transmission,
       R=R, T=T, A=A)
}


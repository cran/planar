\name{field.outside}
\alias{field.outside}
\title{field.outside}
\usage{
  field.outside(d = 1, theta, epsilon, thickness, wavelength,
  polarisation = "p", ...)
}
\arguments{
  \item{d}{distance outside in nm}

  \item{theta}{angle}

  \item{epsilon}{list of permittivities}

  \item{thickness}{vector of layer thickness in nm}

  \item{wavelength}{wavelength in nm}

  \item{polarisation}{polarisation}

  \item{...}{further args passed to multilayer}
}
\value{
  long format data.frame with field on both sides of the
  multilayer as a function of angle
}
\description{
  field profile outside of a ML stack
}
\details{
  returns the field outside of the structure using Fresnel
  coefficients and transfer matrix for comparison
}
\author{
  baptiste Auguie
}


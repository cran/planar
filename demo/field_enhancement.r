## comparison of the calculation of near field enhancement outside of a thin metal film
## with Fresnel reflection and transmission coefficients

library(planar)

gold <- epsAu(633)


test <- field.outside(10, theta=seq(0, pi/2, length=1e3),
                      thickness=c(0, 50, 0),
                      polarisation="s", wavelength=633,
                      epsilon=list(incident = 1.45^2, 633, 1.0^2))
p <- 
ggplot(test) + 
  geom_path(aes(theta, value, colour=side,linetype=model), size=1.2, position=position_dodge(width=0.5))

## lines slightly offset for clarity
p

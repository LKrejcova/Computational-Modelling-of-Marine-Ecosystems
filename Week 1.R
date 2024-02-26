#' ---
#' title: Numerical solution of advection-diffusion models
#' author: Léa Krejci
#' date: 5 February 2024
#' ---

#### Packages ####
library(deSolve)

##### Parameters #####
n <- 20
deltaZ <- 5
depth <- n*deltaZ

u <- 2
D <- 0.1

##### Set up grid ####
z <- seq(deltaZ/2,depth-deltaZ/2,deltaZ)

##### Initial conditions #####
time <- seq(0,10,0.1)
phi <- rep(0,n) # set up phi0
phi[7] <- 10


##### Derivative function #####
derivative <- function(t, phi, params) {
  
  Ja <- rep(0, n+1)
  Jd <- rep(0, n+1)
  
  for(i in 2:n){
    Ja[i] <- u*phi[i-1]
    Jd[i] <- -D*(phi[i]-phi[i-1])/deltaZ
  }
  
  Ja[1] <- 0
  Ja[n+1] <- 0
  Jd[1] <- 0
  Jd[n+1] <- 0
  
  J <- Ja + Jd 
  
  dphidt <- seq(1,n,1)
  
  for (i in 1:n){
    dphidt[i] <- - (J[i+1]- J[i])/deltaZ
  }
  
  return(list(dphidt))
}

times <- seq(0,100,1)  

derivative.out <- ode(phi, times, derivative)
derivative.out


#levelplot(derivative.out)
par(mfrow = c(1,1))
image(x = times, y = z, derivative.out[,ncol(derivative.out):2], 
      ylab = "Depth [m]", xlab = "Days", col = hcl.colors(50, "viridis"))


image(t(derivative.out)[ncol(derivative.out):1,])


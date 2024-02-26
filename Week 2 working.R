#### Packages ####
library(deSolve)

##### Parameters ####
D <- 43.2 #[m^2/day]
Iin <- 350 * 60 * 60 * 24 #[µmol photons /m^2 /day]
kw <- 0.2 #[/m]
kp <- 15e-7 #[m^2 / cell]
m <- 0.01 * 24 #[/day]
u <- 0.96 #[/day]
depth <- 100 #[m]
deltaZ <- 2.5 #[m]
n <- depth/deltaZ 
gmax <- 0.04 * 24
H <- 30 * 60 * 60 *24

params <- c(D = D , Iin = Iin, kw = 0.2, kp = kp, m = m,
            u = u, depth = depth, deltaZ = deltaZ, n = depth/deltaZ, 
            gmax = gmax, H = H)

##### Set up grid ####
z = seq(deltaZ/2, by = deltaZ, to = depth)

##### Initial conditions #####
phi <- rep(0,n) # set up phi0
phi[7] <- 10

# set up light intensity
lightloss <- function(P){
  integral <- cumsum((kw + kp * P) * deltaZ) - 0.5 * deltaZ * (kw + kp*P)
  intensity <- Iin * exp(-integral)
  return(intensity)
}  

plot(lightloss(phi))
##### Derivative function #####
derivative <- function(t, phi, params) {
  
  # set up advective and diffusive flow
  Ja <- rep(0, n+1)
  Jd <- rep(0, n+1)
  
  # calculate advective and diffusive flow
  for(i in 2:n){
    Ja[i] <- u*phi[i-1]
    Jd[i] <- -D*(phi[i]-phi[i-1])/deltaZ
  }
  
  # boundary conditions
  Ja[1] <- 0
  Ja[n+1] <- 0
  Jd[1] <- 0
  Jd[n+1] <- 0
  
  # total flow
  J <- Ja + Jd 
  
  # set up time derivative of phi
  dphi_dt <- seq(1,n,1)
  
  # calculate growth rate
  I = lightloss(phi)
  growth <- gmax * I / (H + I) - m
  
  for (i in 1:n){
    dphi_dt[i] <- growth[i]*phi[i] - (J[i+1]- J[i])/deltaZ
    
  }
  
  return(list(dphi_dt))
}

times <- seq(0,200,1)  

derivative.out <- ode(phi, times, derivative, parms = params)

# time <- derivative.out[,1]
# ph <- derivative.out[,ncol(derivative.out:2)]
# image(time, z, ph, col = hcl.colors(20, "viridis"))

image(x = times, y = z, derivative.out[,ncol(derivative.out):2], 
      ylab = "Distance from seabed [m]", xlab = "Days", col = hcl.colors(50, "viridis"))


# for (n in seq(1, length(times), by = 5)) {
#   plot(derivative.out[n, -1], -z, type = "l", col = "blue", xlim = c(0,500),
#        xlab = "Concentration", ylab = "Depth", #ylim =c(0,100), #, xlim =c(-1,100),
#        main = paste("Phytoplankton at Time =", times[n]))
#   if (n < length(times)) Sys.sleep(0.5)  # Pause between frames for visibility
# }

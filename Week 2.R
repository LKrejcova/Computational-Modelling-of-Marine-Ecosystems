#' ---
#' title: Light intensity in the water column
#' author: Léa Krejci
#' date: 6 February 2024
#' ---

#### Packages ####
library(deSolve)

##### Parameters ####
D <- 43.2 #[m^2/day]
Iin <- 350 * 60 * 60 * 24 #[µmol photons /m^2 /day]
kw <- 0.2 #[/m]
kp <- 15e-12 #[m^2 / cell]
m <- 0.01 * 24 #[/day]
u <- 0.96 #[/day]
depth <- 100 #[m]
deltaZ <- 5 #[m]
n <- depth/deltaZ 
gmax <- 0.04 * 24
H <- 30 * 60 * 60 *24

##### Set up grid ####
z <- seq(deltaZ/2,depth-deltaZ/2,deltaZ)

params <- c(D = D , Iin = Iin-8, kw = 0.2, kp = kp, m = m,
            u = u, depth = depth, deltaZ = deltaZ, n = depth/deltaZ, 
            gmax = gmax, H = H)

## Dummy phytoplankton
P <- rep(0,n)
P[5] <- 2

##### Function to define the light #####
# light <- function(kw, kp, deltaZ, Iin, P){
#   lightloss <- cumsum((kw + kp*P) * deltaZ) - 0.5 * deltaZ * (kw + kp*P)
#   intensity <- Iin * exp(-lightloss)
#   return(intensity)
# }

light <- function(params,phi){
  with(as.list(c(params)),{
  lightloss <- cumsum((kw + kp*P) * deltaZ) - 0.5 * deltaZ * (kw + kp*P)
  intensity <- Iin * exp(-lightloss)
  return(intensity)})
}
light(params,P)

###
plot(y = light(params,P), x = z,
     main = "Light intensity",
     ylab = "Light intensity",
     xlab = "Depth", type = "l")

##### growth rate #####
# growth.rate <- function(params){
#   with(as.list(params), {
#     growth <- gmax * light(params) / (H + light(params)) - m
#     return(growth)})
# }
# growth.rate(params)

#growth <- gmax * light(params) / (H + light(params)) - m

##### Initial conditions #####
time <- seq(0,10,0.1)
phi <- rep(0,n) # set up phi0
phi[7] <- 10


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
  
  # set up depth derivative of fluxes
  dJ_dz <- seq(1,n,1)
  
  growth <- gmax * light(params) / (H + light(params)) - m
  
  for (i in 1:n){
    #dJ_dz[i] <- - (J[i+1]- J[i])/deltaZ
    dphi_dt[i] <- growth*phi - (J[i+1]- J[i])/deltaZ
    
  }
  
  return(list(dphi_dt))
}

times <- seq(0,100,1)  

derivative.out <- ode(phi, times, derivative, parms = params)
for (n in seq(1, length(times), by = 5)) {
  plot(derivative.out[n, -1], -z, type = "l", col = "blue", xlim = c(0,5000),
       xlab = "Concentration", ylab = "Depth", #ylim =c(0,100), #, xlim =c(-1,100),
       main = paste("Phytoplankton at Time =", times[n]))
  if (n < length(times)) Sys.sleep(0.5)  # Pause between frames for visibility
}
image(t(derivative.out)[ncol(derivative.out):2,], col = hcl.colors(50, "viridis"))


image(derivative.out[ncol(derivative.out):2,], 
      ylab = "Depth [m]", xlab = "Days", col = hcl.colors(50, "viridis"))

image(x = times, y = z, derivative.out[,ncol(derivative.out):2], 
      ylab = "Depth [m]", xlab = "Days", col = hcl.colors(50, "viridis"))

image(derivative.out[ncol(derivative.out):2,], 
      ylab = "Depth [m]", xlab = "Days", col = hcl.colors(50, "viridis"))

# Calculate color range
color_range <- seq(min(derivative.out[, 2:nrow(derivative.out)]),
                   max(derivative.out[, 2:nrow(derivative.out)]), 10)

# Format color range values to display only two decimal places
color_range_formatted <- format(color_range, nsmall = 0)

image(t(derivative.out)[ncol(derivative.out):2,],col = hcl.colors(50, "viridis"))

derivative.out
volcano

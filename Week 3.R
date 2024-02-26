#' ---
#' title: Nutrient in the water column
#' author: Léa Krejci
#' date: 12 February 2024
#' ---

#### Packages ####
library(deSolve)
library(fields)

##### Parameters ####
D <- 5 #[m^2/day]
Iin <- 350 * 60 * 60 * 24 #[µmol photons /m^2 /day]
kw <- 0.2 #[/m]
kp <- 15e-12 #[m^2 / cell]
m <- 0.01 * 24 #[/day]
u <- 0.96 #[/day]
depth <- 100 #[m]
deltaZ <- 2.5 #[m]
param.n <- depth/deltaZ
gmax <- 0.04 * 24
H.L <- 30 * 60 * 60 *24 # [µmol photons / m^2 / day]
H.N <- 0.0425 # [mmol nutrient / m^3]
Nb <- 5 # [mmol nutrient / m^3]
yield <- 1e-9 # [mmol nutrient / cell]

# D <- 5 #[m^2/day]
# Iin <- 350#[µmol photons /m^2 /day]
# kw <- 0.2 #[/m]
# kp <- 15e-12 #[m^2 / cell]
# m <- 0.01 #[/day]
# u <- 0.04 #[/day]
# depth <- 100 #[m]
# deltaZ <- 1 #[m]
# param.n <- depth/deltaZ 
# n <- param.n
# gmax <- 0.04
# H.L <- 30
# H.N <- 0.0425 # [mmol nutrient / m^3]
# Nb <- 47.5 # [mmol nutrient / m^3]
# yield <- 1e-9 # [mmol nutrient / cell]

# params <- c(D = D , Iin = Iin, kw = 0.2, kp = kp, m = m,
#             u = u, depth = depth, deltaZ = deltaZ, param.n = depth/deltaZ, 
#             gmax = gmax, H.L = H.L, H.N = H.N)

##### Set up grid ####
z = seq(deltaZ/2, by = deltaZ, to = depth)

##### Initial conditions #####
phi <- rep(10,param.n) # set up phi0
#phi[10] <- 1000

N <- rep(1,param.n) # set up N0

P <- rep(0,param.n)
P[5] <- 20000

# set up vector that contains the inital conditions of phi and N
Y.init <- c(phi,N)

# set up light intensity
lightloss <- function(param, P){
  integral <- cumsum(kw + kp * P) * deltaZ - 0.5 * deltaZ * (kw + kp*P)
  intensity <- Iin * exp(-integral)
  
  return(intensity)
}  

I.out <- lightloss(param, P)
plot(lightloss(param, P))

##### Derivative function #####
derivative <- function(t, Y, params) {
  
  # set up phi
  phi <- Y.init[1:param.n]
  
  # set up nutrients
  N <- Y.init[(param.n+1):(2*param.n)]
  
  # set up advective and diffusive flow
  Ja_phi <- rep(0, param.n+1)
  Jd_phi <- rep(0, param.n+1)
  
  Jd_N <- rep(0, param.n+1)
  
  # calculate advective and diffusive flow
  for(i in 2:param.n){
    Ja_phi[i] <- u*phi[i-1]
    Jd_phi[i] <- -D*(phi[i]-phi[i-1])/deltaZ
    Jd_N[i] <- -D*(N[i]-N[i-1])/deltaZ
  }
  
  # boundary conditions
  Ja_phi[1] <- 0
  Ja_phi[param.n+1] <- 0
  Jd_phi[1] <- 0
  Jd_phi[param.n+1] <- 0
  
  Jd_N[1] <- 0
  Jd_N[param.n+1] <- - D * (Nb - N[param.n])/deltaZ
  
  # total flow
  J <- Ja_phi + Jd_phi 
  
  # set up time derivative of phi and nutrients
  dphi_dt <- seq(1,param.n,1)
  dN_dt <- seq(1,param.n,1)
  
  # calculate growth rate
  I = lightloss(param, phi)
  growth.L <- I / (H.L + I)
  growth.N <- N / (H.N + N)
  
  growth <- gmax * pmin(growth.L, growth.N) - m
  
  for (i in 1:param.n){
    dphi_dt[i] <- growth[i]*phi[i] - (J[i+1] - J[i])/deltaZ
    dN_dt[i] <-  -yield*phi[i]*growth[i] - (Jd_N[i+1] - Jd_N[i])/deltaZ
  }
  
  return(list(c(dphi_dt, dN_dt)))
}

times <- seq(0,1000,1)  

derivative.out <- ode(Y.init, times, derivative)

phi.out <- derivative.out[,2:(param.n+1)]
N.out <- derivative.out[,(param.n+2):(2*param.n+1)]

# plots
image.plot(x = times, y = z, phi.out[,ncol(phi.out):1],
      ylab = "Distance from seabed [m]", xlab = "Days", col = hcl.colors(50, "viridis"),
      main = "Phytoplankton")
#
image.plot(x = times, y = z, N.out[,ncol(N.out):2],
      ylab = "Distance from seabed [m]", xlab = "Days", col = hcl.colors(50, "viridis"),
      main = "Nutrients")

# # plot light as a function of depth
# light <- apply(phi.out, MARGIN = c(1,2), FUN = lightloss)
# image.plot(x = times, y = z, light[,ncol(light):1],
#            ylab = "Distance from seabed [m]", xlab = "Days", col = hcl.colors(50, "viridis"),
#            main = "Light intensity [µmol photons/m2/day]")
# # plot 

#' ---
#' title: Nutrient in the water column
#' author: Léa Krejci
#' date: 12 February 2024
#' ---

#### Packages ####
library(deSolve)
library(fields)

##### Parameters ####
# D <- 5 #[m^2/day]
# Iin <- 350 * 60 * 60 * 24 #[µmol photons /m^2 /day]
# kw <- 0.2 #[/m]
# kp <- 15e-7 #[m^2 / cell]
# m <- 0.01 * 24 #[/day]
# u <- 0.96 #[/day]
# depth <- 100 #[m]
# deltaZ <- 2.5 #[m]
# param.n <- depth/deltaZ
# gmax <- 0.04 * 24
# H.L <- 30 * 60 * 60 *24 # [µmol photons / m^2 / day]
# H.N <- 0.0425 # [mmol nutrient / m^3]
# Nb <- 5 # [mmol nutrient / m^3]
# yield <- 1e-9 # [mmol nutrient / cell]

# leonor's parameters
# n <- 100
# d <- 100
# param.deltaZ <- d/n
# t <- 100
# param.u <- 0.04
# param.D <- 5
# z <- c()
# param.kw <- 0.2
# param.kp <- 15*10^(-7)
# param.Iin <- 350
# param.HI <- 30
# param.gmax <-0.04
# param.m <- 0.01
# param.HN <- 0.0425
# param.Nb <- 47.5
# param.yield <- 10^(-9)


# these parameters work
n <- 75
d <- 100 #[m]
param.deltaZ <- d/n #[m]
t <- 300
param.u <- 0.042 * 24 #[m/day]
param.D <- 43.2 #[m^2/day]
z <- c()
param.kw <- 0.2 #[/m]
param.kp <- 15e-12 #[m^2 / cell]
param.Iin <- 350 * 60 * 60 * 24 #[µmol photons /m^2 /day]
param.HI <- 30 * 60 * 60 *24 # [µmol photons / m^2 / day]
param.gmax <- 0.04 * 24 #[/day]
param.m <- 0.01 * 24 #[/day]
param.HN <- 0.0425 # [mmol nutrient / m^3]
param.Nb <- 5 # [mmol nutrient / m^3]
param.yield <- 1e-9 # [mmol nutrient / cell]

params <- list(param.D = param.D , param.Iin = param.Iin, param.kw = 0.2, param.kp = param.kp,
                param.m = param.m, param.u = param.u, d = d, param.deltaZ = param.deltaZ, n = n,
               param.gmax = param.gmax, param.HI = param.HI, param.HN = param.HN)

no.shading <- list(param.D = param.D , param.Iin = param.Iin, param.kw = 0.2, param.kp = 0,
               param.m = param.m, param.u = param.u, d = d, param.deltaZ = param.deltaZ, n = n,
               param.gmax = param.gmax, param.HI = param.HI, param.HN = param.HN)


##### Set up grid ####
z = seq(param.deltaZ/2, by = param.deltaZ, to = d)

##### Initial conditions #####
phi <- rep(10,n) # set up phi0
#phi[10] <- 1000

N <- rep(1,n) # set up N0

P <- phi[1]

Y.init <- c(phi,N)

# set up light intensity
#leonor
lightloss = function (P, param){
  integral = cumsum((param.kw + param.kp*P)*param.deltaZ) - 0.5*param.deltaZ*(param.kw + param.kp*P)
  I <- param.Iin*exp(- integral)
  
  return(I)
}

lightloss(P, params)

##### Derivative function #####
# leonor
derivative = function(t, Y, params) {
  
  phi <- Y[1:n]
  N <- Y[(n+1):(2*n)]
  
  Ja <- rep(0, n+1)
  Jd <- rep(0, n+1)
  JN <- rep(0,n+1)
  
  for(i in 2:n){
    Ja[i] <- param.u*phi[i-1]
    Jd[i] <- -param.D*(phi[i]-phi[i-1])/param.deltaZ
    JN[i] <- -param.D*(N[i]-N[i-1])/param.deltaZ 
  }
  
  #Boundary conditions  
  Ja[1] <- 0
  Ja[n+1] <- 0
  Jd[1] <- 0
  Jd[n+1] <- 0
  JN[n+1] <- -param.D*(param.Nb-N[n])/param.deltaZ
  JN[1] <- 0
  
  J <- Ja + Jd 
  
  dphi_dt <- seq (1,n,1)
  dN_dt <- seq (1,n,1)
  
  # calculate growth rate
  I = lightloss(phi, params)
  growth.L <- I / (param.HI + I)
  growth.N <- N / (param.HN + N)
  
  growth <- param.gmax * pmin(growth.L, growth.N) - param.m
  
  
  for (i in 1:n){
    dphi_dt[i] <- growth[i]*phi[i] - (J[i+1]- J[i])/param.deltaZ
    dN_dt[i] <- -param.yield*growth[i]*phi[i] - (JN[i+1]- JN[i])/param.deltaZ
  }
  
  return( list(c(dphi_dt, dN_dt)) )
}

times <- seq(0,500,1)  

derivative.out <- ode(Y.init, times, derivative, parms = params)

phi.out <- derivative.out[,2:(n+1)]
N.out <- derivative.out[,(n+2):(2*n+1)]

# plots 
image.plot(x = times, y = z, phi.out[,ncol(phi.out):1],
           ylab = "Distance from seabed [m]", xlab = "Days", col = hcl.colors(50, "viridis"),
           main = "Phytoplankton")
#
image.plot(x = times, y = z, N.out[,ncol(N.out):2],
           ylab = "Distance from seabed [m]", xlab = "Days", col = hcl.colors(50, "viridis"),
           main = "Nutrients")

# plot light as a function of depth
light <- sapply(seq_len(ncol(phi.out)), function(i) lightloss(phi.out[, i], param))

image.plot(x = times, y = z, light[,ncol(light):1],
           ylab = "Distance from seabed [m]", xlab = "Days", col = hcl.colors(50, "viridis"),
           main = "Light intensity [µmol photons/m2/day]")

plot(x = lightloss(phi.out[nrow(phi.out),], params), y = n:1, type = "l", main = "Light intensity",
     ylab = "Distance from seabed [m]", xlab = "Light intensity [mmol photon/ m^2 / day]",
     pch = 3, lwd = 2)
param.kp <- 0
lines(x = lightloss(phi.out[nrow(phi.out),],no.shading), y = n:1, type = "l", lty = 2, lwd = 2)
legend("bottomright", legend = c("Self-shading", "No self-shading"), lty = c(1,2), lwd = 2)
param.kp <- 15e-12 #[m^2 / cell]
grid()

# plot light and nutrient limitation as a function of depth
growth.lim.L <- lightloss(phi.out[nrow(phi.out),]) / (param.HI + lightloss(phi.out[nrow(phi.out),]))
growth.lim.N <- N.out[nrow(N.out),] / (param.HN + N.out[nrow(N.out),])


par(mar=c(5, 4, 4, 8), xpd=TRUE) # adds space next to the plot so we can put a legend there
plot(x = growth.lim.L, y = n:1, type = "l", xlim = c(0, max(growth.lim.N)),
     main = "Limitation by light or nutrients", ylab = "Distance from seabed [m]",
     xlab = "Growth limitation factor", lwd = 3)
lines(x = growth.lim.N, y = n:1, type = "l", lwd = 3, lty = 2)
legend("topright", inset=c(-0.5, 0), legend = c("Light", "Nutrients"), lty = c(1,2), lwd = 2)

plot(x = phi.out[nrow(phi.out),], y = n:1, type = "l",
     main = "Limitation by light or nutrients", ylab = "Distance from seabed [m]",
     xlab = "Cell concentration [cell/m^3]", lwd = 3)

# alternative approach
# light <- lightloss(phi.out[,1])
# for (i in 2:ncol(phi.out)){
#   new_light <- lightloss(phi.out[,i])
#   light <- cbind(light, new_light)
# }

#' ---
#' title: Nutrient in the water column
#' author: Léa Krejci
#' date: 26 February 2024
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


# # these parameters work
# n <- 75
# d <- 100 #[m]
# param.deltaZ <- d/n #[m]
# t <- 300
# param.u <- 0.042 * 24 #[m/day]
# param.D <- 43.2 #[m^2/day]
# z <- c()
# param.kw <- 0.2 #[/m]
# param.kp <- 15e-12 #[m^2 / cell]
# param.Iin <- 350 * 60 * 60 * 24 #[µmol photons /m^2 /day]
# param.HI <- 30 * 60 * 60 *24 # [µmol photons / m^2 / day]
# param.gmax <- 0.04 * 24 #[/day]
# param.m <- 0.03 #[/day]
# param.HN <- 0.0425 # [mmol nutrient / m^3]
# param.Nb <- 5 # [mmol nutrient / m^3]
# param.yield <- 1e-9 # [mmol nutrient / cell]
# param.gamma <- 0.25 # [m^3 /mmol N / day]
# param.w <- 5 # [m/day]
# param.rem <- 0.1 # [/day]

# new parameters in mmol.N
n <- 75
d <- 400 #[m]
param.deltaZ <- d/n #[m]
t <- 300
param.u <- 0.1 #[m/day]
param.D <- 5 #[m^2/day]
z <- c()
param.kw <-  0.0375 #[/m]
param.kp <-  0.05 #[m^2 / mmol.N]
param.Iin <- 200 #[W /m^2]
param.m <-  0.03 #[/day]
param.HN <-  0.3 # [mmol.N / m^3]
param.Nb <- 20 # [mmol.N / m^3]
param.gamma <- 1.5 # [m^3 /mmol N / day]
param.w <- 15 # [m/day]
param.rem <- 0.1 # [/day]
param.alpha <- 0.1 # [m^2 / W / day]
param.mu <- 0.5 # [/ day]

# reset graphical parameters
#dev.off()

params <- list(param.D = param.D , param.Iin = param.Iin, param.kw = 0.2, param.kp = param.kp,
               param.m = param.m, param.u = param.u, d = d, param.deltaZ = param.deltaZ, n = n,
               param.mu = param.mu, param.HN = param.HN)

no.shading <- list(param.D = param.D , param.Iin = param.Iin, param.kw = 0.2, param.kp = 0,
                   param.m = param.m, param.u = param.u, d = d, param.deltaZ = param.deltaZ, n = n,
                   param.mu = param.mu, param.HN = param.HN)


##### Set up grid ####
z = seq(param.deltaZ/2, by = param.deltaZ, to = d)

##### Initial conditions #####
phi <- rep(1,n) # set up phi0
#phi[10] <- 1000

N <- rep(20,n) # set up N0

P <- phi[1]

detr <- rep(0,n)

D <- detr[1]

Y.init <- c(phi,N,detr)

# set up light intensity
lightloss = function (P, D, param){
  integral = cumsum(((param.kw + param.kp*(P+D))*param.deltaZ)) - 0.5*param.deltaZ*(param.kw + param.kp*(P+D))
  I <- param.Iin*exp(- integral)
  
  return(I)
}

plot(x = lightloss(phi,detr, params), y = n:1, col = "blue", type = "l")

##### Derivative function #####
derivative = function(t, Y, params) {
  
  phi <- Y[1:n]
  N <- Y[(n+1):(2*n)]
  detr <- Y[(2*n+1):(3*n)]
  
  Ja <- rep(0, n+1)
  Jd <- rep(0, n+1)
  JN <- rep(0, n+1)
  Ja_D <- rep(0, n+1)
  Jd_D <- rep(0, n+1)
  
  for(i in 2:n){
    Ja[i] <- param.u*phi[i-1]
    Jd[i] <- -param.D*(phi[i]-phi[i-1])/param.deltaZ
    JN[i] <- -param.D*(N[i]-N[i-1])/param.deltaZ 
    Ja_D[i] <- param.w*detr[i-1]
    Jd_D[i] <- -param.D*(detr[i]-detr[i-1])/param.deltaZ
  }
  
  #Boundary conditions  
  Ja[1] <- 0
  Ja[n+1] <- 0
  
  Jd[1] <- 0
  Jd[n+1] <- 0
  
  JN[1] <- 0
  JN[n+1] <- -param.D*(param.Nb-N[n])/param.deltaZ
  
  Ja_D[1] <- 0
  Ja_D[n+1] <- param.w*detr[n]
  
  Jd_D[1] <- 0
  Jd_D[n+1] <- 0
  
  J <- Ja + Jd 
  
  JDet <- Ja_D + Jd_D
  
  dphi_dt <- seq(1,n,1)
  dN_dt <- seq(1,n,1)
  ddetr_dt <- seq(1,n,1)
  
  # calculate growth rate
  I = lightloss(phi, detr, params)
  #growth.L <- I / (param.HI + I)
  growth.L <- param.alpha*I / sqrt(param.mu^2 + (param.alpha*I)^2)
  growth.N <- N / (param.HN + N)
  
  growth <- param.mu * pmin(growth.L,growth.N)
  
  
  for (i in 1:n){
    dphi_dt[i] <- growth[i]*phi[i] - param.m*phi[i] - param.gamma*phi[i]^2 - (J[i+1] - J[i])/param.deltaZ
    dN_dt[i] <- -growth[i]*phi[i] + param.rem*detr[i] - (JN[i+1] - JN[i])/param.deltaZ
    ddetr_dt[i] <- param.m*phi[i] + param.gamma*phi[i]^2 - param.rem*detr[i] - (JDet[i+1] - JDet[i])/param.deltaZ
  }
  
  return( list(c(dphi_dt, dN_dt, ddetr_dt)) )
}

times <- seq(0,3000,.1)  

derivative.out <- ode(Y.init, times, derivative, parms = params)

phi.out <- derivative.out[,2:(n+1)]
N.out <- derivative.out[,(n+2):(2*n+1)]
Det.out <- derivative.out[, (2*n+2):(3*n+1)]

#check if converged
max((phi.out[nrow(phi.out),]-phi.out[nrow(phi.out)-1,])/phi.out[nrow(phi.out),])*100
max((N.out[nrow(N.out),]-N.out[nrow(N.out)-1,])/N.out[nrow(N.out),])*100
max((Det.out[nrow(Det.out),]-Det.out[nrow(Det.out)-1,])/Det.out[nrow(Det.out),])*100

# plots 
# image.plot(x = times, y = z, phi.out[,ncol(phi.out):1],
#            ylab = "Distance from seabed [m]", xlab = "Days", col = hcl.colors(50, "viridis"),
#            main = "Phytoplankton", legend.lab="mmol Nitrogen / m^3", legend.line = - 2.5,
#            cex.lab = 1.5, cex.axis = 1.5, legend.cex = 1.5)
# #
# image.plot(x = times, y = z, N.out[,ncol(N.out):2],
#            ylab = "Distance from seabed [m]", xlab = "Days", col = hcl.colors(50, "viridis"),
#            main = "Nutrients", legend.lab="mmol Nitrogen / m^3", legend.line = -2.5,
#            cex.lab = 1.5, cex.axis = 1.5, legend.cex = 1.5)
# 
# image.plot(x = times, y = z, Det.out[,ncol(Det.out):2],
#            ylab = "Distance from seabed [m]", xlab = "Days", col = hcl.colors(50, "viridis"),
#            main = "Detritus", legend.lab="mmol Nitrogen / m^3", legend.line = - 2.5,
#            cex.lab = 1.5, cex.axis = 1.5, legend.cex = 1.5)

# plot light as a function of depth
#light <- sapply(seq_len(ncol(phi.out)), function(i) lightloss(phi.out[, i],  param))

# image.plot(x = times, y = z, light[,ncol(light):1],
#            ylab = "Distance from seabed [m]", xlab = "Days", col = hcl.colors(50, "viridis"),
#            main = "Light intensity [µmol photons/m2/day]")

plot(x = lightloss(phi.out[nrow(phi.out),], Det.out[nrow(Det.out),],  params), y = seq(d,1,-param.deltaZ), type = "l", main = "Light intensity",
     ylab = "Distance from seabed [m]", xlab = "Light intensity [W/m^2]",
     pch = 3, lwd = 2.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
param.kp <- 0
lines(x = lightloss(phi.out[nrow(phi.out),], Det.out[nrow(Det.out),], params), y = seq(d,1,-param.deltaZ), 
      type = "l", lty = 2.5, lwd = 2, col = "red")
grid(lwd = 1.5, col = "dimgrey")
legend("bottomright", legend = c("Self-shading", "No self-shading"), lty = c(1,2), lwd = 2, col = c("black", "red"))
param.kp <- 0.05 #[m^2 / mmolN]


# plot light and nutrient limitation as a function of depth
growth.lim.L <- param.alpha*lightloss(tail(phi.out,1), tail(Det.out,1)) / sqrt(param.mu^2 + (param.alpha*lightloss(tail(phi.out,1), tail(Det.out,1)))^2) 
growth.lim.N <- tail(N.out,1) / (param.HN + tail(N.out,1))


plot(x = growth.lim.L, y = seq(d,1,-param.deltaZ), type = "l", xlim = c(0, max(growth.lim.L)),
     main = "Limitation by light or nutrients", ylab = "Distance from seabed [m]",
     xlab = "Growth limitation factor", lwd = 3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, lty = 2)
lines(x = growth.lim.N, y = seq(d,1,-param.deltaZ), type = "l", lwd = 3, lty = 1)
legend("bottomleft", inset=c(0.1, 0), legend = c("Light", "Nutrients"), lty = c(2,1), lwd = 3, cex = 1.5, bty = "n")
grid(lwd = 1.5, col = "dimgrey")
#legend("topright", inset=c(-0.3, 0), legend = c("Light", "Nutrients"), lty = c(1,2), lwd = 2, cex = 1.1, bty = "n")
lines(x = phi.out[nrow(phi.out),], y = seq(d,1,-param.deltaZ), type = "l",
      main = "Limitation by light or nutrients", ylab = "Distance from seabed [m]",
      xlab = "Cell concentration [cell/m^3]", lwd = 3, col = "red")

# plot steady state solution
# normalize function
normalize_min_max <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

phi.norm <- normalize_min_max(tail(phi.out,1))
N.norm <- normalize_min_max(tail(N.out,1))
Det.norm <- normalize_min_max(tail(Det.out,1))
light.norm <- normalize_min_max(lightloss(tail(phi.out,1), tail(Det.out,1)))

plot(x = phi.norm, y = seq(d,1,-param.deltaZ), type = "l", lwd = 3, col = "#009e73",
     xlab = "Normalized concentrations of P, N, and D, and light intensity",
     ylab = "Distance from seabed [m]", cex.lab = 1.5, cex.axis = 1.5)
grid(lwd = 1.5, col = "dimgrey")
lines(x = N.norm, y = seq(d,1,-param.deltaZ), type = "l", lwd = 3, col = "#cc79a7")
lines(x = Det.norm, y = seq(d,1,-param.deltaZ), type = "l", lwd = 3, col = "#d55e00")
lines(x = light.norm, y = seq(d,1,-param.deltaZ), type = "l", lwd = 3, lty = 2)
legend("bottomleft", inset=c(0.5, 0), legend = c("Phytoplankton", "Nutrients", "Detritus", "Light"),
       lty = c(1,1,1,2), lwd = 3, cex = 1.5, bty = "n", col = c("#009e73", "#cc79a7", "#d55e00",1))


# alternative approach
# light <- lightloss(phi.out[,1])
# for (i in 2:ncol(phi.out)){
#   new_light <- lightloss(phi.out[,i])
#   light <- cbind(light, new_light)
# }
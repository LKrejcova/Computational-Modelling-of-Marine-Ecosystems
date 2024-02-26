library(deSolve)

#1) setting the parameters 

u <- 0.04*24  # Vertical water velocity (advection velocity) settling velocity
D <- 5*(60*60*24)/10000 # Diffusion coefficient
#D[1:20] <- D #When including thermocline in stratified water
#D[21:100] <- 0
#d <- 100 # Depth of the water column
n <- 100 # Number of grid cells along the water column
DeltaZ <- 1 # Vertical spacing between grid cells
kw <- 0.2 
kp <- 15 * (10^-12)
I0 <- 350*(60*60*24)/1000
z_values <- seq(0.5 * DeltaZ, n, by = DeltaZ)
gmax <- 0.04*24
H <- (30*(60*60*24))/1000
m <- 0.01*24

# Initial conditions (setting a patch in the middle)
P <- rep(0,n)
P[1] <- 10

#-----------------------------------------------------
#Creating the lights distribution through the watercolumn

Light <- function(kw,kp,DeltaZ,P,I0){
  
  #Integral calculations - cumulative sums of the light in depth
  integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw*kp*P)
  
  #calculation the light values 
  I = I0*exp(-integral)
  
  return(I)
}

#-------------------------------------------------------------
#Creating the plots for the light attenuation 


I = Light(kw,kp,DeltaZ,P,I0)

plot(I,-z_values, type ="l", xlab = "Light", ylab = "Depth")


#--------------------------------------------------------------
#Calculating the division rate

g <- gmax * (I/(H+I))

#-------------------------------------------------------------
#Calculating the division growth rate 

r <- g-m 

plot(r, -z_values, type = "l")

plot(P, -z_values,type = "l")

#------------------------------------------------
#Doing the diffusion and advection 
#------------------------------------------------

derivative <- function(times, P, params) {
  J_a <- numeric(n + 1)
  J_d <- numeric(n + 1)
  
  ix <- 2:n
  
  # Calculate advective fluxes (Ja)
  J_a[ix] <- u * P[ix - 1]
  
  # Calculate diffusive fluxes (Jd)
  J_d[ix] <- -D * ((P[ix] - P[ix - 1]) / DeltaZ)
  
  # Set boundary fluxes to 0
  J_a[1] <- 0 # Boundary condition
  J_d[1] <- 0
  J_a[n + 1] <- 0
  J_d[n + 1] <- 0
  
  # Calculate dPhi/dt
  dP_dt <- numeric(n)
  dPdz <- numeric(n)
  #for (i in 1:n) {
    dJ_dt <- (-(J_a[ix+1] + J_d[ix+1] - J_a[1:n]-J_d[1:n]) / DeltaZ)
    
    dPdz <-  dJ_dt + (r*P) 
  #}
  
  return(list(dPdz))
}


params <- list(n = n, DeltaZ = DeltaZ, u = u, D = D, r = r, m = m, gmax = gmax, I = I, H=H )
# Time
times <- seq(0, 50, by = 1) # Time vector

# ODE solving
result <- ode(y = P, times = times, func = derivative, parms = params)


#---------------------------------------------------------------
#Plotting
#---------------------------------------------------------------


par(mfrow=c(1,1))

for (n in seq(1, length(times), by = 5)) {
  plot(result[n, -1], -z_values, type = "l", col = "blue", xlim = c(0,5000),
       xlab = "Concentration", ylab = "Depth", #ylim =c(0,100), #, xlim =c(-1,100),
       main = paste("Phytoplankton at Time =", times[n]))
  if (n < length(times)) Sys.sleep(0.5)  # Pause between frames for visibility
}



# Calculate color range
color_range <- seq(min(result[, nrow(result):2]), max(result[, nrow(result):2]), length.out = 10)

# Format color range values to display only two decimal places
color_range_formatted <- format(color_range, nsmall = 0)

# Generate plot
image(result[, nrow(result):2], col = terrain.colors(10))

# Create legend
legend("topleft", legend = color_range_formatted, fill = terrain.colors(10), title = "Values")

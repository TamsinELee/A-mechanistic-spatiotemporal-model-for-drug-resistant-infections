
# Supplementary information 1
# Generate the simulated data. 
# Adapted from Appendix S1 from Hefley et al. 2017, Ecology Letters.
# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fele.12763&file=ele12763-sup-0001-AppendixS1.pdf
# Alterations for this paper denoted with triple hash: ###

# Load packages
library(raster) 
library(fields) 
library(RColorBrewer) 
library(rasterVis) 
library(truncnorm) 
library(mvtnorm) 
library(gridExtra) 
library(fBasics) 
library(coda)
library(matrixcalc)
library(epsiwal)
library(wesanderson)

### Wes Anderson colour palette
pal  <- wes_palette("Zissou1", 100, type = "continuous")

# Link function g(.) 
link <- function(nu){ptruncnorm(nu, a=0, b=Inf, mean = 0, sd = 1)}

# First-order neighbourhood matrix from a RasterLayer object 
neighbourhood <- function(raster){ 
  nn         <- matrix(,length(raster[]),4) 
  for(i in 1:dim(nn)[1]){ 
    loc      <- adjacent(raster,i)[,2] 
    ln       <- loc[which((loc+1)==i)] 
    rn       <- loc[which((loc-1)==i)] 
    bn       <- loc[which((loc-dim(raster)[2])==i)]
    tn       <- loc[which((loc+dim(raster)[2])==i)] 
    nn[i,1]  <- if(length(ln)>0){ln}else{0} 
    nn[i,2]  <- if(length(rn)>0){rn}else{0} 
    nn[i,3]  <- if(length(bn)>0){bn}else{0} 
    nn[i,4]  <- if(length(tn)>0){tn}else{0} 
  } 
  nn 
}

# Propagator matrix for plain diffusion PDE 
propagator.plain <- function(NN, mu, lambda, dx, dy, dt){ 
  H                <- matrix(0, dim(NN)[1], dim(NN)[1]) 
  for(i in 1:dim(H)[1]){ 
    if(length(which(NN[i,] > 0)) == 4){ 
      H[i,i]       <- 1-2*mu[i]*(dt/dx^2 + dt/dy^2) + dt*lambda[i] 
      H[i,NN[i,1]] <- dt/dx^2 * mu[i] 
      H[i,NN[i,2]] <- dt/dx^2 * mu[i] 
      H[i,NN[i,3]] <- dt/dy^2 * mu[i] 
      H[i,NN[i,4]] <- dt/dy^2 * mu[i]} 
  }
  H
}

# RasterStack of c(s,t) 
calc.c <- function(H, c0, t.steps, t.keep){ 
  c.all   <- c0 
  c.all[] <- H %*% c0[] + c0[]
  c.all   <- stack(mget(rep("c.all",t.steps))) 
  for(t in 2:t.steps){ 
    c.all[[t]][] <- H %*% c.all[[t-1]][] + c0[]
  } 
  c.all[[t.keep]] 
}

### Model parameters
theta       <- c(80, 70, 65, 60, 60)
phi         <- c(0.08, 0.09, 0.1, 0.15, 0.12)
alpha       <- c(-8, 1)
gamma       <- c(0.2, 0.1)
beta        <- c(-10)

### Location of resistance hotspots
d           <- cbind(c(0.1, 0.1, 0.5,  0.6, 0.9), c(0.9, 0.1, 0.5, 0.2, 0.5))
Nd          <- nrow(d)

# Number of spatial grid cells 
q           <- 100^2

# set seed 
set.seed(6982)

### Simulate covariates 
w           <- sort(runif(q, -0.5, 0.5))
z           <- sort(runif(q, -0.5, 0.5))

mu          <- raster(, nrows = q^0.5, ncols = q^0.5, xmn = 0, xmx = 1, ymn = 0, ymx = 1, crs = NA) 
mu[]        <- exp(model.matrix(~z) %*% alpha) 
# plot(mu)

### Plot mu 
### Figures/mu.png
levelplot(mu, xlab=list("latitude", cex = 3), ylab=list("longitude", cex = 3),  
          margin = FALSE, scales = list(draw = FALSE), 
          pretty = T, col.regions = pal, colorkey = list(at=seq(-0.6, 0.6, 0.1), 
                                                         labels=list(at=c(-0.6, 0, 0.6), 
                                                                     labels=c("low", "medium", "high"), cex = 3)))#
### Figures/mu_withNos.png
levelplot(mu, xlab=list("latitude", cex = 3), ylab=list("longitude", cex = 3),  
          margin = FALSE, scales = list(draw = FALSE), 
          pretty = T, col.regions = pal, colorkey = list(at=seq(-0.6, 0.6, 0.1), 
          labels = list(at=c(-0.6, 0, 0.6), cex = 3)))# + 
       #   layer(panel.dotplot(x = d[,1], y = d[,2], col="black", cex=1))

lambda      <- raster(, nrows = q^0.5, ncols = q^0.5, xmn = 0, xmx = 1, ymn = 0, ymx = 1, crs = NA) 
lambda[]    <- model.matrix(~w) %*% gamma 

### Plot lambda
### Figures/lambda.png
levelplot(lambda, xlab=list("latitude", cex = 3), ylab=list("longitude", cex = 3),  
          margin = FALSE, scales = list(draw = FALSE), 
          pretty = T, col.regions = pal, colorkey = list(at=seq(-0.6, 0.6, 0.1), 
          labels=list(at=c(-0.6, 0, 0.6), 
          labels=c("low", "medium", "high"), cex = 3)))

### Figures/lambda_withNos.png
levelplot(lambda, xlab=list("latitude", cex = 3), ylab=list("longitude", cex = 3),  
          margin = FALSE, scales = list(draw = FALSE), 
          pretty = T, col.regions = pal, colorkey = list(at=seq(-0.6, 0.6, 0.1), 
                                                         labels=list(at=c(-0.6, 0, 0.6), 
                                                          cex = 3)))

# Scaling factor 
us.fact     <- 5

# Diffusion rate for homogenized pde 
mu.bar      <- aggregate(mu, fact = us.fact, na.rm = TRUE, 
                         fun = function(x, na.rm) { 1/mean(1/x) 
                         })

# Growth rate for homogenized pde 
lambda.bar  <- mu.bar * aggregate(lambda/mu, fact = us.fact, na.rm = TRUE, FUN = mean)

# First-order neighbourhood matrix 
NN          <- neighbourhood(mu.bar)

# Propagator matrix 
H           <- propagator.plain(NN = NN, mu = mu.bar[], lambda = lambda.bar[], dx = us.fact/q^0.5, dy = us.fact/q^0.5, dt = 1)

u0          <- raster(, nrows = q^0.5, ncols = q^0.5, 
                      xmn = 0, xmx = 1, 
                      ymn = 0, ymx = 1, crs = NA) 

### Distance from hotspot
D          <- matrix(NA, q, Nd)
for (kH in 1:Nd){
  dist.d   <- raster(, nrows = q^0.5, ncols = q^0.5, 
                      xmn = 0, xmx = 1, 
                      ymn = 0, ymx = 1, crs = NA) 
  D[, kH] <- rdist(data.frame(SpatialPoints(dist.d)), matrix(d[kH,], 1, 2)) 
}

### Initial conditions based on distance from hotspots
u0[]       <- (exp(-D[,1]^2 / phi[1]^2) / sum(exp(-D[,1]^2 / phi[1]^2))) * theta[1] 
for (kH in 2:Nd){
  u0[]     <- u0[] + (exp(-D[,kH]^2/phi[kH]^2)/sum(exp(-D[,kH]^2/phi[kH]^2))) * theta[kH] 
}
### Plotting initial conditions
u0_plot    <- levelplot(u0, xlab=list("latitude", cex = 2), ylab=list("longitude", cex = 2),  margin = FALSE, scales = list(draw = FALSE), 
              main = "truth", pretty=T, col.regions = pal, colorkey = TRUE) 
u0_plot 

c0         <- raster(, nrows = q^0.5/us.fact, ncols = q^0.5/us.fact, 
                     xmn = 0, xmx = 1, ymn = 0, ymx = 1, crs = NA) 
c0[]       <- extract(mu * u0, SpatialPoints(mu.bar))

# Make raster layers 
c1  <- c0;  c2 <- c0; c3  <- c0; c4  <- c0; c5  <- c0; c6  <- c0; c7  <- c0 
c8  <- c0;  c9 <- c0; c10 <- c0; c11 <- c0; c12 <- c0; c13 <- c0; c14 <- c0 
c15 <- c0; c16 <- c0; c17 <- c0; c18 <- c0; c19 <- c0; c20 <- c0

### Fill raster layers with term to account for distance from hotspot.
c1[]       <- H%*%c0[]  + c0[]
c2[]       <- H%*%c1[]  + c0[]
c3[]       <- H%*%c2[]  + c0[]
c4[]       <- H%*%c3[]  + c0[]
c5[]       <- H%*%c4[]  + c0[]
c6[]       <- H%*%c5[]  + c0[]
c7[]       <- H%*%c6[]  + c0[]
c8[]       <- H%*%c7[]  + c0[]
c9[]       <- H%*%c8[]  + c0[]
c10[]      <- H%*%c9[]  + c0[]
c11[]      <- H%*%c10[] + c0[]
c12[]      <- H%*%c11[] + c0[]
c13[]      <- H%*%c12[] + c0[]
c14[]      <- H%*%c13[] + c0[]
c15[]      <- H%*%c14[] + c0[]
c16[]      <- H%*%c15[] + c0[]
c17[]      <- H%*%c16[] + c0[]
c18[]      <- H%*%c17[] + c0[]
c19[]      <- H%*%c18[] + c0[]
c20[]      <- H%*%c19[] + c0[]

# Stack raster layers 
c.all      <- stack(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20)

c.all      <- calc.c(H,c0,20,1:20)

color      <- colorRampPalette(rev(brewer.pal(11, "Spectral")), bias = 1) 
levelplot(c.all, cuts = 254, margin = FALSE, scales = list(draw = FALSE), 
          names.attr = paste("t =", 1:20), col.regions = colorRampPalette(rev(brewer.pal(11, "Spectral")), bias = 1))

# Calculate u(s,t) 
u.all      <- disaggregate(c.all, us.fact)/mu
#levelplot(u.all, main = list("u: density"), cuts = 254, margin = FALSE, scales = list(draw = FALSE), 
#          names.attr = paste("t =", 1:20), col.regions = colorRampPalette(rev(brewer.pal(11, "Spectral")), bias = 1))

### Plotting the actual (simulated) density
levelplot(u.all, cuts = 254, margin = FALSE, scales = list(draw = FALSE), 
          names.attr = paste("t =", 1:20),  par.strip.text = list(cex = 1.5), 
          col.regions = colorRampPalette(rev(brewer.pal(11, "Spectral"))), 
          colorkey=list(at=seq(0, 20, 0.01), 
                        labels=list(at=c(0, 20), labels=c("none", "high"), cex = 2))) 

### Generating and plotting individual covariate, such as age. 
x           <- rbeta(length(vec(u.all[])), 2, 5)
xr          <- raster(, nrows = q^0.5, ncols = q^0.5, xmn = 0, xmx = 1, ymn = 0, ymx = 1, crs = NA) 
xr[]        <- x

levelplot(xr, xlab=list("latitude", cex = 2), ylab=list("longitude", cex = 2),  margin = FALSE, scales = list(draw = FALSE), 
          main = "x", pretty = T, col.regions = pal, colorkey = list(at=seq(0, 1, 0.1), 
                  labels=list(at=c(0, 0.5, 1), 
                              labels=c("high", "medium", "low")))) 

# Probability of being sampled.
p          <- link(exp(x * beta) * vec(u.all[]))
p.all      <- u.all 
p.all[]    <- matrix(p, dim(p.all[])[1],dim(p.all[])[2])
levelplot(p.all, main = list("p: prob sampled"), cuts = 254, margin = FALSE, scales = list(draw = FALSE), 
          names.attr = paste("t =", 1:20), col.regions = colorRampPalette(rev(brewer.pal(11, "Spectral")), bias = 1))

y          <- rbinom(length(p), 1, p) 
y.all      <- u.all 
y.all[]    <- y 
#levelplot(y.all, main = list("y: data"), margin = FALSE, scales = list(draw = FALSE), 
#          names.attr = paste("t =", 1:20), 
#          col.regions = c("black", "green"), cuts = 1, , colorkey = FALSE)

### Plotting the data
levelplot(y.all, xlab="", ylab="",  margin = FALSE, scales = list(draw = FALSE),
                      names.attr = paste("t =", 1:20), 
                      col.regions = c("lightgray", "#F21A00"), cuts = 1,
                      pretty = T, colorkey = FALSE) 

# Saving the data. 
df1     <- data.frame(y = vec(y.all[]), x = x, z = rep(z, 20), w = rep(w, 20), 
                  t = rep(1:20, each = q), s1 = data.frame(SpatialPoints(y.all))[, 1], 
                  s2 = data.frame(SpatialPoints(y.all))[, 2], cell = rep(1:q, 20))

###
### Exploring different sampling - different levels of bias in the collection.
### df1 is sampling uniformly across space and time. 
### To add bias collection, need to stratify it by time. 
###

set.seed(4223)
NSampled                   <- 20000 # 30000 20000 10000 2000 (corresponds to 15, 10, 5, 1 %)

PresentIdx                 <- which(df1$y > 0) # Bias sampling will collect from here more. 
AbsentIdx                  <- which(df1$y == 0)

PercSampled                <- 0.1 # 0.3, 0.l5, 0.1, 0.05, 0.01

timeidx                    <- seq(1, 200000+1, 10000)
K                          <- c()
PresentIdx                 <- list()
AbsentIdx                  <- list()

for (tidx in 1:20){
  thistimeidx              <- c(timeidx[tidx]:(timeidx[tidx + 1]))
  #print(head(thistimeidx))
  thisPresentIdx           <- thistimeidx[which(df1$y[thistimeidx] > 0)]
  thisAbsentIdx            <- thistimeidx[which(df1$y[thistimeidx] == 0)]
  PresentIdx[[tidx]]       <- thisPresentIdx
  AbsentIdx[[tidx]]        <- thisAbsentIdx
  
  NSampledYr               <- round(PercSampled * 10000, 0)
  if (NSampledYr < length(thisPresentIdx)){
    thisK                  <- sample(thisPresentIdx, replace = FALSE, size = NSampledYr)
  } else {
    thisK1                 <- thisPresentIdx
    thisK2                 <- sample(thisAbsentIdx, replace = FALSE, size = (NSampledYr - length(thisPresentIdx)))
    thisK                  <- c(thisK1, thisK2) 
  }
  K                        <- c(K, thisK)
}

df2                        <- df1[K, ]

# Keep sampled and no resistance
y2          <- y
SampledVec  <- ifelse(y2[K] == 0, 1, 2)
y2[K]       <- SampledVec
y2[-K]      <- 0
y2.all      <- y.all 
y2.all[]    <- y2

# Very biased sampling
levelplot(y2.all, xlab="", ylab="",  margin = FALSE, scales = list(draw = FALSE),
          names.attr = paste("t =", 1:20),  par.strip.text = list(cex = 1.5), 
          col.regions = c("white",  "#F21A00", "#78B7C5"), 
          pretty=T, colorkey = FALSE)

# Unbiased sampling
set.seed(4223) 
K1                           	<- sample(seq(1:200000), replace = FALSE, size = NSampled)
K                            	<- K1
df2                           <- df1[K, ]

# Keep sampled and no resistance
y2          <- y
SampledVec  <- ifelse(y2[K] == 0, 1, 2)
y2[K]       <- SampledVec
y2[-K]      <- 0
y2.all      <- y.all 
y2.all[]    <- y2

# Unbiased sampling
levelplot(y2.all, xlab="", ylab="",  margin = FALSE, scales = list(draw = FALSE),
          names.attr = paste("t =", 1:20),  par.strip.text = list(cex = 1.5), 
          col.regions = c("white",  "#F21A00", "#78B7C5"), 
          pretty=T, colorkey = FALSE) 

# Keep sampled and no resistance
y2          <- y
SampledVec  <- ifelse(y2[K] == 0, 1, 2)
y2[K]       <- SampledVec
y2[-K]      <- 0
y2.all      <- y.all 
y2.all[]    <- y2

levelplot(y2.all, xlab="", ylab="",  margin = FALSE, scales = list(draw = FALSE),
          names.attr = paste("t =", 1:20), 
          col.regions = c("white",  "#F21A00", "#78B7C5"), 
          pretty=T, colorkey = FALSE) 

y2_1         <- raster(, nrows = q^0.5, ncols = q^0.5, 
                       xmn = 0, xmx = 1, ymn = 0, ymx = 1, crs = NA) 
y2_1[]       <- y2[1:q]
levelplot(y2_1, xlab=list("latitude", cex = 2), ylab=list("longitude", cex = 2),  margin = FALSE, scales = list(draw = FALSE), 
          pretty=T, col.regions = c("white",  "#F21A00", "#78B7C5"), colorkey = FALSE)  +  
         layer(panel.dotplot(x = d[,1], y = d[,2], col="black", cex=3))
y2_20         <- raster(, nrows = q^0.5, ncols = q^0.5, 
                        xmn = 0, xmx = 1, ymn = 0, ymx = 1, crs = NA) 
y2_20[]       <- y2[(length(y3) - q + 1):length(y3)]
levelplot(y2_20, xlab=list("latitude", cex = 2), ylab=list("longitude", cex = 2),  margin = FALSE, scales = list(draw = FALSE), 
          pretty=T, col.regions = c("white",  "#F21A00", "#78B7C5"), colorkey = FALSE)  

#y3          <- y2 - y2
#y3[which(y2 == 2)] <- 100 * rbeta(length(which(y2 == 2)), 2, 5) # to show age is different for each data points
#y3.all      <- y.all 
#y3.all[]    <- y3 

#yearidx <- seq(1, 200000, q)
#for (kk in 1:19){
#  print(sum(y3[yearidx[kk + 1]:yearidx[kk]]))
#}

#levelplot(y3.all, xlab="", ylab="",  margin = FALSE, scales = list(draw = FALSE),
#          names.attr = paste("t =", 1:20), 
#          col.regions = c("white", "#F21A00"), cuts = 1,
#          pretty=T, colorkey = FALSE)+ 
#  layer(panel.dotplot(x = d[,1], y = d[,2], col="black", cex=1))


### Plot to show age is different for each data point
y4            <- y
y4[K]         <- 100 * rbeta(length(which(y2 == 2)), 2, 5) 
y4[-K]        <- 0
y4.all        <- y.all 
y4.all[]      <- y4
y4_20         <- raster(, nrows = q^0.5, ncols = q^0.5, 
                        xmn = 0, xmx = 1, ymn = 0, ymx = 1, crs = NA) 
y4_20[]       <- y4[(length(y4) - q + 1):length(y4)]
colfunc       <-  colorRampPalette(c("white", "#F21A00"))
#y4_t20_x.pdf
levelplot(log(y4_20 + 1), xlab=list("latitude", cex = 2), ylab=list("longitude", cex = 2),  margin = FALSE, scales = list(draw = FALSE), 
          pretty=T, col.regions = colfunc,  colorkey = list(at = seq(0, 5, 0.5),
                                                            labels = list(at=c(0, 0.75, 4.75), 
                                                            labels =  c("no data","young", "old"), cex = 2))) +
          layer(panel.dotplot(x = d[,1], y = d[,2], col="black", cex=3))


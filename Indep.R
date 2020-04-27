#The independence example

#install packages
#install.packages("RandomFields")
#install.packages("spatstat")
#install.packages("fields")
#install.packages("rgeos")
# INLA has to be downloaded from r-inla.org

#load packages
library(RandomFields)
library(spatstat)
library(fields)
library(rgeos)
library(INLA)

#observation window
win <- owin(c(0,1), c(0,1))

#resolution of the field
spatstat.options(npixel=100)

## Process 1
#generate spatially structured effect
sigma2x1 <- 0.7
rangenom1 <- 0.5
kappa1 <- sqrt(8)/rangenom1

set.seed(1)
struc1 <- rLGCP('matern', mu=0, 
              var=sigma2x1, scale=1/kappa1, nu=1, win=win)

#get intensity
Lam11 <- attr(struc1, 'Lambda')
summary(as.vector(ll11 <- log(Lam11$v)))

#plotting
par(mfrow=c(1,1))
image.plot(list(x=Lam11$yrow, y=Lam11$xcol, z=ll11), asp=1)

#generate latent field and point pattern
set.seed(1)
mwn1 <- 6 + rnorm(10000, mean = 0, sd = 0.1) #mean plus white noise
mat1 <- matrix(mwn1, nrow=100, ncol=100)
image1 <- im(mat1, xrange=c(0,1), yrange=c(0,1))
max(ll11) - min(ll11)
max(ll21) - min(ll21)

#simulation
set.seed(1)
pp1 <- rLGCP('matern', image1, 
             var=sigma2x1, scale=1/kappa1, nu=1, win=win)

#get point locations
xy1 <- cbind(pp1$x, pp1$y)[,2:1]

#get intensity
Lam12 <- attr(pp1, 'Lambda')
summary(as.vector(ll12 <- log(Lam12$v)))

#plotting
par(mfrow=c(1,1))
image.plot(list(x=Lam12$yrow, y=Lam12$xcol, z=ll12), asp=1)
points(xy1, pch=19) 



## Process 2
#generate spatially structured effect
sigma2x2 <- 0.5
rangenom2 <- 0.3
kappa2 <- sqrt(8)/rangenom2

set.seed(2)
struc2 <- rLGCP('matern', mu=0, 
                var=sigma2x2, scale=1/kappa2, nu=1, win=win)

#get intensity
Lam21 <- attr(struc2, 'Lambda')
summary(as.vector(ll21 <- log(Lam21$v)))

#plotting
image.plot(list(x=Lam21$yrow, y=Lam21$xcol, z=ll21), asp=1)

#generate latent field and point pattern
set.seed(2)
mwn2 <- 5 + rnorm(10000, mean = 0, sd = 0.1) #mean plus white noise
mat2 <- matrix(mwn2, nrow=100, ncol=100)
image2 <- im(mat2, xrange=c(0,1), yrange=c(0,1))

#simulation
set.seed(2)
pp2 <- rLGCP('matern', image2, 
             var=sigma2x2, scale=1/kappa2, nu=1, win=win)

#get point locations
xy2 <- cbind(pp2$x, pp2$y)[,2:1]

#get intensity
Lam22 <- attr(pp2, 'Lambda')
summary(as.vector(ll22 <- log(Lam22$v)))

#plotting
par(mfrow=c(1,1))
image.plot(list(x=Lam22$yrow, y=Lam22$xcol, z=ll22), asp=1)
points(xy2, pch=19) 






## Inference
#create a mesh
loc.d <- t(matrix(c(0,0,1,0,1,1,0,1,0,0), 2))
(nv <- (mesh <- inla.mesh.2d(
  loc.domain=loc.d, offset=c(.1, .5), 
  max.edge=c(.06, 0.5), cutoff=.02))$n) 
 
#SPDE model definition
spde1 <- inla.spde2.pcmatern(
  mesh=mesh, alpha=2,
  prior.range=c(0.05, 0.05), 
  prior.sigma=c(1, 0.05)) 

spde2 <- inla.spde2.pcmatern(
  mesh=mesh, alpha=2,
  prior.range=c(0.05, 0.05), 
  prior.sigma=c(2, 0.05))  

#dual mesh function
inla.mesh.dual <- function(mesh) {
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0) 
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}
dmesh <- inla.mesh.dual(mesh)

#spatial polygons
domainSP <- SpatialPolygons(list(Polygons(
  list(Polygon(loc.d)), '0')))

#weights
sum(w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i,], domainSP))
    return(gArea(gIntersection(dmesh[i,], domainSP)))
  else return(0)
}))

#observations
(n1 <- nrow(xy1))
y1 <- cbind(rep(0:1, c(nv, n1)))
(n2 <- nrow(xy2))
y2 <- rep(0:1, c(nv, n2))

#expected observations
e1 <- c(w, rep(0, n1)) 
e2 <- c(w, rep(0, n2))

#create projector matrix
lmat1 <- inla.spde.make.A(mesh, loc=xy1)
lmat2 <- inla.spde.make.A(mesh, loc=xy2)

imat <- Diagonal(nv, rep(1, nv))

A1 <- rBind(imat, lmat1)
A2 <- rBind(imat, lmat2)

#data stack
ny1 <- nrow(as.matrix(y1))
ny2 <- nrow(as.matrix(y2))
y = matrix(NA, (ny1+ny2), 2)
y[1:ny1, 1] = y1
y[1:ny2 + ny1, 2] = y2

#joint model
stk1 <- inla.stack(data=list(y=cbind(y1,NA), e=e1), 
                      A=list(1,A1), tag='pp1',
                      effects=list(list(b1=rep(1,nv+n1)), list(i=1:nv))) 
stk2 <- inla.stack(data=list(y=cbind(NA,y2), e=e2),
                      A=list(1,A2), tag='pp2',
                      effects=list(list(b2=rep(1,nv+n2)), list(j=1:nv)))

stk <- inla.stack(stk1, stk2)



sys.time.indepm <- system.time(res.ind.m <- inla(y ~ 0 + b1 + b2 + f(i, model=spde1) + f(j, model=spde2), 
                                              family=c('poisson', 'poisson'), data=inla.stack.data(stk),
                                              control.predictor=list(A=inla.stack.A(stk)),
                                              #control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE),
                                              E=inla.stack.data(stk)$e))

#hpd intervals

ms1=res.ind.m$marginals.hyperpar$`Stdev for i`
hpds1 <- inla.hpdmarginal(0.95, ms1)

mr1=res.ind.m$marginals.hyperpar$`Range for i`
hpdr1 <- inla.hpdmarginal(0.95, mr1)

ms2=res.ind.m$marginals.hyperpar$`Stdev for j`
hpds2 <- inla.hpdmarginal(0.95, ms2)

mr2=res.ind.m$marginals.hyperpar$`Range for j`
hpdr2 <- inla.hpdmarginal(0.95, mr2)

#plots for Process 1
#par(mfrow=c(1,3), mar=c(3,3,1,0.3), mgp=c(2,1,0)) 
par(mfrow=c(1,1), mai=c(1,1,1,1))
plot(res.ind.m$marginals.fix[[1]], type='l', mgp=c(2.2,0.6,0), xlim=c(4,8),
     xlab=expression(beta[1]), ylab='Density')
abline(v=mean(mwn1), col=2)
abline(v=res.ind.m$summary.fixed$`0.025quant`[1], lty=2)
abline(v=res.ind.m$summary.fixed$`0.975quant`[1], lty=2)

plot(res.ind.m$marginals.hyperpar[[2]], type='l', xlim = c(0.2,2.5), mgp=c(2.1,0.6,0), 
     xlab=expression(sigma[f[1]]^2), ylab='Density')
abline(v=sigma2x1, col=2)
abline(v=hpds1[,1], lty=2)
abline(v=hpds1[,2], lty=2)

plot(res.ind.m$marginals.hyperpar[[1]], type='l', xlim = c(0.2,2), mgp=c(1.8,0.6,0),
     xlab=expression(rho[1]), ylab='Density')
abline(v=sqrt(8*1)/kappa1, col=2)
abline(v=hpdr1[,1], lty=2)
abline(v=hpdr1[,2], lty=2)


#plots for Process 2
#par(mfrow=c(1,3), mar=c(3,3,1,0.3), mgp=c(2,1,0)) 
plot(res.ind.m$marginals.fix[[2]], type='l',  mgp=c(2,0.6,0), xlim=c(4,6),
     xlab=expression(beta[2]), ylab='Density')
abline(v=mean(mwn2), col=2)
abline(v=res.ind.m$summary.fixed$`0.025quant`[2], lty=2)
abline(v=res.ind.m$summary.fixed$`0.975quant`[2], lty=2)

plot(res.ind.m$marginals.hyperpar[[4]], type='l', xlim = c(0.2,2.5), mgp=c(2.1,0.6,0), 
     xlab=expression(sigma[f[2]]^2), ylab='Density')
abline(v=sigma2x2, col=2)
abline(v=hpds2[,1], lty=2)
abline(v=hpds2[,2], lty=2)

plot(res.ind.m$marginals.hyperpar[[3]], type='l', xlim = c(0,2), mgp=c(1.8,0.6,0),
     xlab=expression(rho[2]), ylab='Density')
abline(v=sqrt(8*1)/kappa2, col=2)
abline(v=hpdr2[,1], lty=2)
abline(v=hpdr2[,2], lty=2)


pgrid0 <- inla.mesh.projector(mesh, xlim=0:1, ylim=0:1, dims=c(101,101))
prd0.m <- inla.mesh.project(pgrid0, (mean(mwn1) + res.ind.m$summary.ran$i$mean))
par(mfrow=c(1,1))
image.plot(prd0.m, asp=1)
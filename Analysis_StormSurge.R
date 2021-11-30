rm(list=ls())
#setwd("~/Projects/MultiFidelity/")

# install ARCokrig package from CRAN
# install.package("ARCokrig")

library(ARCokrig)

###
nump <- function(y.training.hf, y.testing, y.pred, Lower95, Upper95, y.predSE){
  RMSPE=sqrt(mean((y.testing-y.pred)^2))
  PCI95 = mean( (Lower95<y.testing) & (Upper95>y.testing) )
  LCI95 = mean(Upper95 - Lower95)
  crps = mean(CRPS(y.testing, y.pred, y.predSE))
  
  y.testing.mean = apply(y.testing, 2, mean)
  
  y.training.mean = apply(y.training.hf, 2, mean)
  res = y.pred - matrix(rep(y.training.mean, 
                            times=dim(y.testing)[1]), byrow=T, nrow=dim(y.testing)[1])
  NSME = 1 - sum((y.testing - y.pred)^2) / sum((res)^2)
  
  measure = c(RMSPE, PCI95, LCI95, crps, NSME)
  names(measure) <- c("RMSPE", "P95CI", "L95CI", "CRPS", "NSME")
  return(measure)
}

#####################################################################
#####################################################################
##### load data

load("./data/StormSurge.RData")

#### prepare data for cokriging
n_runs = dim(inputs)[1]
k_outputs = dim(lonlat)[1]
set.seed(123)
id.input = 1:n_runs
id.training = sample(id.input, 60)
id.output = 1:k_outputs 
z.hf = PSE.hf[id.training, id.output]
set.seed(234)
id.remain = id.input[-id.training]
id.lf = c(sample(id.remain, 150), id.training)
z.lf = PSE.lf[id.lf, id.output]

# convert long/lat of lanfall to distance
library(fields)
lon0 = min(inputs[ ,6])
lat0 = max(inputs[, 7])
reference = matrix(c(lon0, lat0), nrow=1)
dist = c(rdist.earth(reference, inputs[ ,c(6,7)]))
inputs.6Dim = cbind(inputs[ ,1:5], dist)

input.hf = inputs.6Dim[id.training, ]
input.lf = inputs.6Dim[id.lf, ]

input.max = apply(inputs.6Dim, 2, max)
input.min = apply(inputs.6Dim, 2, min)
delta = input.max - input.min

#####################################################################
############## Non-Nested Design
#####################################################################

set.seed(1234)
id.h.l = c(1:150, 150 + sample(1:60, 50))

zc = z.lf[id.h.l, ]  # only keep 50 designs from high-fidelity level
zf = z.hf
Dc = input.lf[id.h.l, ]
Df = input.hf
output = list(zc, zf)
input = list(Dc, Df)

inputdata.new = inputdata[-id.training, ]
input.new = inputs.6Dim[-id.training, ]
output.new = PSE.hf[-id.training, id.output]
output.new.lf = PSE.lf[-id.training, id.output]
n.new = dim(output.new.lf)[1]

input.new.orig = input.new

y.testing = output.new
z.testing = output.new.lf

##########################################################################
########### PPCokriging 
##########################################################################

##### parameter estimation
cov.model = "matern_5_2"
formula = list(~1, ~1)

opt = list(maxit=800)
prior = list(name="JR")
NestDesign = FALSE

obj = mvcokm(formula=formula, output=output, input=input, 
             cov.model="matern_5_2",prior=prior,opt=opt,
             NestDesign=NestDesign)

ts = proc.time() 
obj.fitted = mvcokm.fit(obj)
te = proc.time() - ts

param.est = mvcokm.param(obj.fitted)
b = param.est$coeff
sigma2 = param.est$var
summary(c(b[[1]]))
summary(c(b[[2]][2,]))
summary(sigma2[[2]])

set.seed(12345)
t1 = proc.time()
pred = mvcokm.condsim(obj=obj.fitted, input.new=input.new)
t.pred = proc.time() - t1

y.testing = output.new
z.testing = output.new.lf
n.new = dim(y.testing)[1]

## overall result
y.testing.overall = y.testing 
y.pred.overall = pred$mu[[2]]
Lower95.overall = pred$lower95[[2]]
Upper95.overall = pred$upper95[[2]]
y.predSE.overall = pred$SE[[2]]

measure.overall = nump(y.training.hf = zf,
                       y.testing=y.testing.overall, 
                       y.pred=y.pred.overall,
                       Lower95=Lower95.overall, 
                       Upper95=Upper95.overall,
                       y.predSE=y.predSE.overall)
measure.overall

##########################################################################
########### PPGaSP for high-fidelity data only 
##########################################################################

library(RobustGaSP)

#### only with 60 high-fidelity runs

input.pp = Df
p.x = dim(Df)[2]
formula = ~1
colnames(input.pp) = paste0("x", 1:p.x)
trend = model.matrix(formula, data.frame(input.pp))
colnames(input.new) = paste0("x", 1:p.x)
trend.testing = model.matrix(formula, data.frame(input.new))

### no nugget
fit.PP = ppgasp(design=input.pp, response=zf, trend=trend)
show(fit.PP)

pred.pp = predict(fit.PP, input.new, testing_trend=trend.testing)

y.pred = pred.pp$mean

measure.overall = nump(y.training.hf=zf,
                      y.testing=y.testing, 
                       y.pred=pred.pp$mean,
                       Lower95=pred.pp$lower95, 
                       Upper95=pred.pp$upper95,
                       y.predSE=pred.pp$sd)
measure.overall


#### only with  low-fidelity runs

input.pp = Dc
p.x = dim(Dc)[2]
formula = ~1
colnames(input.pp) = paste0("x", 1:p.x)
trend = model.matrix(formula, data.frame(input.pp))
colnames(input.new) = paste0("x", 1:p.x)
trend.testing = model.matrix(formula, data.frame(input.new))

### no nugget
fit.PP = ppgasp(design=input.pp, response=zc, trend=trend)
show(fit.PP)

pred.pp = predict(fit.PP, input.new, testing_trend=trend.testing)

y.pred = pred.pp$mean

measure.overall = nump(y.training.hf=zc,
                       y.testing=y.testing, 
                       y.pred=pred.pp$mean,
                       Lower95=pred.pp$lower95, 
                       Upper95=pred.pp$upper95,
                       y.predSE=pred.pp$sd+1e-10*length(pred.pp$sd))
measure.overall


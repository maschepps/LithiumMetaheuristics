
#Load in 
#Source files
setwd('C:/Users/Admin/Desktop/Github/JRSS_Lithium2/Publish/src')
# setwd('~/Desktop/JRSS_Lithium/src')
file_list = list.files()
sapply(file_list, source)


## Genetic covariate optimization
lower = rep(0, 6)
upper = rep(8, 6)
ptm <- proc.time()
g1 = ecr(fitness.fun = D1cov, representation = "float", n.objectives = 1,
         n.dim = 6, survival.strategy = "plus",
         lower = lower, upper = upper,
         mu = 50, lambda = 5,
         mutator = setup(mutGauss, sdev = 2, lower = lower, upper = upper),
         terminators = list(stopOnMaxTime(60)))
eg1=proc.time() - ptm

ptm <- proc.time()
g2 = ecr(fitness.fun = Ds1cov, representation = "float", n.objectives = 1,
         n.dim = 6, survival.strategy = "plus",
         lower = lower, upper = upper,
         mu = 50, lambda = 5,
         mutator = setup(mutGauss, sdev = 2, lower = lower, upper = upper),
         terminators = list(stopOnMaxTime(30)))
eg2=proc.time() - ptm

lower = rep(0, 10)
upper = rep(8, 10)
ptm <- proc.time()
g3 = ecr(fitness.fun = D2cov, representation = "float", n.objectives = 1,
         n.dim = 10, survival.strategy = "plus",
         lower = lower, upper = upper,
         mu = 50, lambda = 5,
         mutator = setup(mutGauss, sdev = 2, lower = lower, upper = upper),
         terminators = list(stopOnMaxTime(60)))
eg3=proc.time() - ptm

ptm <- proc.time()
g4 = ecr(fitness.fun = Ds2cov, representation = "float", n.objectives = 1,
         n.dim = 10, survival.strategy = "plus",
         lower = lower, upper = upper,
         mu = 50, lambda = 5,
         mutator = setup(mutGauss, sdev = 2, lower = lower, upper = upper),
         terminators = list(stopOnMaxTime(60)))
eg4=proc.time() - ptm
## No genetic covariate
## Table 3

lower = rep(0, 6)
upper = rep(8, 6)
ptm <- proc.time()
ng1 = ecr(fitness.fun = D1, representation = "float", n.objectives = 1,
         n.dim = 6, survival.strategy = "plus",
         lower = lower, upper = upper,
         mu = 50, lambda = 5,
         mutator = setup(mutGauss, sdev = 2, lower = lower, upper = upper),
         terminators = list(stopOnMaxTime(30)))
eng1=proc.time() - ptm

ptm <- proc.time()
ng2 = ecr(fitness.fun = Ds1, representation = "float", n.objectives = 1,
         n.dim = 6, survival.strategy = "plus",
         lower = lower, upper = upper,
         mu = 50, lambda = 5,
         mutator = setup(mutGauss, sdev = 2, lower = lower, upper = upper),
         terminators = list(stopOnMaxTime(30)))
eng2=proc.time() - ptm

lower = rep(0, 10)
upper = rep(8, 10)
ptm <- proc.time()
gn3 = ecr(fitness.fun = D2, representation = "float", n.objectives = 1,
         n.dim = 10, survival.strategy = "plus",
         lower = lower, upper = upper,
         mu = 50, lambda = 5,
         mutator = setup(mutGauss, sdev = 2, lower = lower, upper = upper),
         terminators = list(stopOnMaxTime(60)))
eng3=proc.time() - ptm

ptm <- proc.time()
ng4 = ecr(fitness.fun = Ds2, representation = "float", n.objectives = 1,
         n.dim = 10, survival.strategy = "plus",
         lower = lower, upper = upper,
         mu = 50, lambda = 5,
         mutator = setup(mutGauss, sdev = 2, lower = lower, upper = upper),
         terminators = list(stopOnMaxTime(60)))
eng4=proc.time() - ptm
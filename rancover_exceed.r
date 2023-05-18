# Script Title: rancovr_exceed.r
# Author: Massimo Cavallaro
# Description: Main script used to generate the warning scores with rancovr,
# given the baseline (`baseline.matrix.p` and `baseline.matrix.f`) and the
# observation matrices (`Positives` and `Fails`). It can be also applied to
# the downsampled baseline and observation matrices.
# Dependencies: rancovr package, create_matrices.r

set.seed(1)

# We can skeep the initial months are they have zero cases.
t = 350:782
ncols = length(t)

LTLA[,c('y','x')] = rancovr::vlatlong2km(LTLA[, c("latitude", "longitude")])
rho = rancovr::distance.between.points(LTLA[,'x'], LTLA[,'y']) 

ws.p = matrix(NA, nrow = nrow(LTLA), ncol = ncols)
ws.f = matrix(NA, nrow = nrow(LTLA), ncol = ncols)
#
N.p = matrix(NA, nrow = nrow(LTLA), ncol = ncols)
N.f = matrix(NA, nrow = nrow(LTLA), ncol = ncols)

n.cylinders = 500000
mean.p = c()
mean.f = c()

for (TT in t){ # For all observation weeks, recursively.
  # check.point = Sys.time()
  # Generate cylinders to cover observations up until week TT.
  cylinders.p = rancovr::CreateCylinders(
    Positives,
    baseline.matrix.p,
    week.range = c(TT-200, TT),
    GT=16,
    rho = rho,
    n.cylinders = n.cylinders,
    coord.df = LTLA,
    only.last = T)
  cylinders.f = rancovr::CreateCylinders(
    Fails,
    baseline.matrix.f,
    week.range = c(TT-200, TT),
    GT=16,
    rho = rho,
    n.cylinders = n.cylinders,
    coord.df = LTLA,
    only.last = T)

  # we need to check that GT is longer the "weight" vector used to
  # exponentially weight the nomeber of cases in the baseline matrix

  I = TT - t[1] + 1
  mean.p = c(mean.p, mean(cylinders.p$n_obs))
  mean.f = c(mean.f, mean(cylinders.f$n_obs))

  # Find the warning scores applying the function `warning.score2` to all
  # locations. The function `warning.score2` also returns the number of
  # covering cylinders (N.p and N.f) which can be used to estimate
  # uncertainty.
  tmp = apply(LTLA, 1, FUN = rancovr::warning.score2, TT, cylinders.f)
  ws.f[,I] = tmp[1,]
  N.f[,I] = tmp[2,]

  tmp = apply(LTLA, 1, FUN = rancovr::warning.score2, TT, cylinders.p)
  ws.p[,I] = tmp[1,]
  N.p[,I] = tmp[2,]

  print(c(TT, t[length(t)]))
  # print(Sys.time() - check.point)
}

colnames(ws.p) = t[1]:(t[1]+NCOL(ws.p)-1)
colnames(ws.f) = t[1]:(t[1]+NCOL(ws.p)-1)

colnames(N.p) = t[1]:(t[1]+NCOL(ws.p)-1)
colnames(N.f) = t[1]:(t[1]+NCOL(ws.p)-1)

dump(c('ws.p', 'ws.f', 'N.p', 'N.f'),
  file=paste0('ws_', as.character(t[1]),'-', as.character(t[length(t)]), '_',Sys.Date() ,'create_matrices_matt2b.R'))

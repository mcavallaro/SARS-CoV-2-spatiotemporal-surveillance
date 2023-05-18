# Script Title: p_val.r
# Author: Massimo Cavallaro
# Description: Compute the exceedance probabilities for each individual
# geographical locations, withouth considering neighbourhoods.
# Dependencies: extraDistr package, create_matrices.r

# Compute exceedance probability under the Poisson model ...
p.val.f = ppois(Fails, lambda = baseline.matrix.f, lower.tail = T)
p.val.p = ppois(Positives, lambda = baseline.matrix.p, lower.tail = T)

# and under the negative binomial model:
p.val.bb.f = extraDistr::pbbinom(Fails,
                     size=Total,
                     alpha = baseline.matrix.f + 0.5,
                     beta = (1 - baseline.matrix.f/Total) * Total + 0.5, lower.tail = T) |> matrix(ncol = ncol(p.val.f))
p.val.bb.p = extraDistr::pbbinom(Positives,
                     size=Total,
                     alpha = baseline.matrix.p + 0.5,
                     beta = (1 - baseline.matrix.p/Total) * Total + 0.5, lower.tail = T) |> matrix(ncol = ncol(p.val.f))

#############
# Try with quick model fitting

loglikelihood_poisson_model<-function(x, n, lambda){
  idx = !(c(n) == 0)
  sum(dpois(c(x)[idx], c(lambda)[idx], log=T))
}

loglikelihood_bbinom_model<-function(x, n, lambda){
  idx = !(c(n) == 0)
  alpha = c(lambda) + 0.5
  beta = (1 - c(lambda)/c(n)) * c(n) + 0.5
  sum(dbbinom(c(x)[idx], c(n)[idx], alpha[idx], beta[idx], log=T), na.rm=T)
}

c(
  loglikelihood_poisson_model(Fails, Total, baseline.matrix.f),
  loglikelihood_poisson_model(Positives, Total, baseline.matrix.p))

c(loglikelihood_bbinom_model(Fails, Total, baseline.matrix.f),
  loglikelihood_bbinom_model(Positives, Total, baseline.matrix.p))


IC<-function(loglik_function, x, n, lambda, n.parameters.per.ltla = 1){
  loglik = loglik_function(x, n, lambda)
  k = log(prod(dim(x)) )
  n.params = prod(dim(x)) * n.parameters.per.ltla
  ic= c(-2 * loglik + k * n.params,  -2 * loglik + 2 * n.params)
  names(ic) = c('BIC', 'AIC')
  return(ic)
}



c(
  IC(loglikelihood_poisson_model, Fails, Total, baseline.matrix.f, 2),
  IC(loglikelihood_poisson_model, Positives, Total, baseline.matrix.p, 2))
c(
  IC(loglikelihood_bbinom_model, Fails, Total, baseline.matrix.f, 2),
  IC(loglikelihood_bbinom_model, Positives, Total, baseline.matrix.p, 2))

# 
# loglikelihood_bbinom_model_1param<-function(x, n, p, param){
#   idx = !(c(n) == 0)
#   P =   c(p) * param
#   alpha = c(n * P) + 0.5
#   beta = (1 - c(P)) * c(n) + 0.5
#   r = dbbinom(c(x)[idx], c(n)[idx], alpha[idx], beta[idx], log=T)
#   return(sum(r, na.rm = T))
# }

loglikelihood_bbinom_model_2param<-function(x, n, p, param){
  if ((param[1] < 0) | (param[2] < 0)){
    return(-Inf)
  }
  idx = !(c(n) == 0)
  C = param[1] + c(n) * param[2]
  alpha = C * c(p) + 0.5
  beta = (1 - c(p)) * C + 0.5
  r = dbbinom(c(x)[idx], c(n)[idx], alpha[idx], beta[idx], log=T)
  return(sum(r, na.rm = T))
}



values = c()
optimal_params = matrix(NA, nrow=0, ncol=2)
for(i in 1:20){
  o = optim(par=runif(2, 0, 100),
            fn = function(param){-loglikelihood_bbinom_model_2param(Positives, Total, baseline.matrix.p / Total, param)})
  values = c(values, o$value)
  optimal_params = rbind(optimal_params, o$par)
}
optimal_params_means = colMeans(optimal_params)
optimal_params_top = optimal_params[which.min(values),]


Positives_optimal_params_means = optimal_params_means
Positives_optimal_params_top = optimal_params_top

writeLines(sprintf("Positives optimal par %.2f %.2f", optimal_params_means[1], optimal_params_means[2]))
writeLines(sprintf("Positives optimal par %.2f %.2f", optimal_params_top[1], optimal_params_top[2]))
# Positives optimal par 27.03 0.00
# Positives optimal par 23.83 0.00
C1 = optimal_params_top[1]

loglikelihood_bbinom_model_2param(Positives, Total, baseline.matrix.p / Total, Positives_optimal_params_means)
loglikelihood_bbinom_model_2param(Positives, Total, baseline.matrix.p / Total, Positives_optimal_params_top)
loglikelihood_bbinom_model_2param(Positives, Total, baseline.matrix.p / Total, c(0,1))
loglikelihood_poisson_model(Positives, Total, baseline.matrix.p)


values = c()
optimal_params = matrix(NA, nrow=0, ncol=2)
for(i in 1:20){
  o = optim(par=runif(2, 0, 100),
            fn = function(param){-loglikelihood_bbinom_model_2param(Fails, Total, baseline.matrix.f / Total, param)})
  values = c(values, o$value)
  optimal_params = rbind(optimal_params, o$par)
}
optimal_params_means = colMeans(optimal_params)
optimal_params_top = optimal_params[which.min(values),]


Fails_optimal_params_means = optimal_params_means
Fails_optimal_params_top = optimal_params_top

loglikelihood_bbinom_model_2param(Fails, Total, baseline.matrix.f / Total, Fails_optimal_params_top)
loglikelihood_bbinom_model_2param(Fails, Total, baseline.matrix.f / Total, Fails_optimal_params_means)
loglikelihood_bbinom_model_2param(Fails, Total, baseline.matrix.f / Total, c(0,1))
loglikelihood_poisson_model(Fails, Total, baseline.matrix.f)

writeLines(sprintf("Fails optimal par %.2f %.2f", optimal_params_means[1], optimal_params_means[2]))
writeLines(sprintf("Fails optimal par %.2f %.2f", optimal_params_top[1], optimal_params_top[2]))
#Fails optimal par 30.66 0.00
#Fails optimal par 25.23 0.00
C2 = optimal_params_top[1]

# We can compute exceedance probabilities with these fitted parameters
p.val.bb.f = pbbinom(Fails,
                     size=Total,
                     alpha = baseline.matrix.f + 0.5,
                     beta = (1 - baseline.matrix.f/Total) * Total + 0.5, lower.tail = T) |> matrix(ncol = ncol(p.val.f))
p.val.bb.p = pbbinom(Positives,
                     size=Total,
                     alpha = baseline.matrix.p + 0.5,
                     beta = (1 - baseline.matrix.p/Total) * Total + 0.5, lower.tail = T) |> matrix(ncol = ncol(p.val.f))
C1 = 23.83
p.val.fitted_bb.f = pbbinom(Fails,
                     size=Total,
                     alpha = baseline.matrix.f/Total * C1 + 0.5,
                     beta = (1 - baseline.matrix.f/Total) * C1 + 0.5, lower.tail = T) |> matrix(ncol = ncol(p.val.f))
C2 = 25.23
p.val.fitted_bb.p = pbbinom(Positives,
                     size=Total,
                     alpha = baseline.matrix.p/Total * C2 + 0.5,
                     beta = (1 - baseline.matrix.p/Total) * C2 + 0.5, lower.tail = T) |> matrix(ncol = ncol(p.val.f))

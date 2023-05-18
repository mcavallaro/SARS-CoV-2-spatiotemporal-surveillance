# Script name: create_matrices.r
# Author: Massimo Cavallaro
# Description: Generate the baseline matrices, which represent the null model
# for the observation matrices.
# Dependencies: read_data.r

# Proportion of s-gene positives and s-gene negative tests:
PP = Positives / Total
PF = Fails / Total

# p_{L,D} = \sum_{d=1}^T exp( - A d ) P_{L,D-d}
# Q_{L,D}  ~  Poisson ( N_{L,D} * p_{L,D} )

# Create matrix objects to store the baselines for PP and PF:
p = matrix(0,  nrow = nrow(PP), ncol = ncol(PP))
f = matrix(0,  nrow = nrow(PF), ncol = ncol(PF))

GT = 5.5 # generation time
pp = 1 - exp(-1/GT)
window_width = 15
weight = exp(-1/GT * (window_width-1):0 )
weight = weight/sum(weight)
# or
# weight = pp * (1 - pp) ** ((window_width-1):0)

for (i in (window_width+1):NCOL(PP)){
  p[,i] = rowSums( t( t( PP[,(i-window_width):(i-1)]) * weight), na.rm = TRUE)
}

for (i in (window_width+1):NCOL(PP)){
  f[,i] = rowSums( t( t( PF[,(i-window_width):(i-1)]) * weight), na.rm = TRUE)
}

# Baseline matrices for the number of positive and negative tests
baseline.matrix.p = p * Total
baseline.matrix.f = f * Total

# Apply threshold to remove warnings from one event when zero is expected
min_thresh = 0.3553468
baseline.matrix.p[baseline.matrix.p < min_thresh] = min_thresh
baseline.matrix.f[baseline.matrix.f < min_thresh] = min_thresh

# When Total is zero, the baseline matrices should be zero, since there are no
# observations (either Fails or Positives). To locate zero division entries in
# the baseline matrices, search if there are entries such that baseline matrix
# is zero but observations are not.
baseline.index.f = which(
	(baseline.matrix.f == 0) & !(Fails==0),
	arr.ind = T
	)
baseline.index.p = which(
	(baseline.matrix.p == 0) & !(Positives==0),
	arr.ind = T
	)

dump(
	c('Positives', 'Fails', 'Total','baseline.matrix.p', 'baseline.matrix.f'),
    file = 'data.R'
    )

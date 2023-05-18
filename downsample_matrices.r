# Script name: downsample_matrices.r
# Author: Massimo Cavallaro
# Description: Generate the observation and baseline matrices assuming that we had half of the available tests.
# Dependencies:

source('read_data.r')

# Downsample observation matrices
tmp = matrix(
	sample(c(-0.001,0.001), nrow(Positives) * ncol(Positives), replace = T),
	nrow = nrow(Positives)
	)
Positives_DS = round(0.5 * (Positives + tmp))
Fails_DS = round(0.5 * (Fails - tmp))
Total_DS = Positives_DS + Fails_DS
rm(tmp)

# Proportion of s-gene positives and s-gene negative tests:
PP_DS = Positives_DS / Total_DS
PF_DS = Fails_DS / Total_DS

# Create matrix objects to store the baselines:
p_DS = matrix(0, nrow = nrow(PP), ncol = ncol(PP))
f_DS = matrix(0, nrow = nrow(PF), ncol = ncol(PF))

GT = 5.5 # generation time
pp = 1 - exp(-1/GT)
window_width = 15
weight = exp(-1/GT * (window_width-1):0 )
weight = weight/sum(weight)
# or
# weight = pp * (1 - pp) ** ((window_width-1):0)

for (i in (window_width+1):NCOL(PP)){
  p_DS[,i] = rowSums( t( t( PP_DS[,(i-window_width):(i-1)]) * weight), na.rm = TRUE)
}

for (i in (window_width+1):NCOL(PP)){
  f_DS[,i] = rowSums( t( t( PF_DS[,(i-window_width):(i-1)]) * weight), na.rm = TRUE)
}

# Baseline matrices for the number of positive and negative tests
baseline.matrix.p_DS = p_DS * Total_DS
baseline.matrix.f_DS = f_DS * Total_DS

# Apply threshold to remove warnings from one event when zero is expected
min_thresh = 0.3553468
baseline.matrix.p[baseline.matrix.p < min_thresh] = min_thresh
baseline.matrix.f[baseline.matrix.f < min_thresh] = min_thresh

baseline.matrix.p_DS[baseline.matrix.p_DS < min_thresh] = min_thresh
baseline.matrix.f_DS[baseline.matrix.f_DS < min_thresh] = min_thresh

dump(
	c('Positives_DS', 'Fails_DS', 'Total_DS','baseline.matrix.p_DS', 'baseline.matrix.f_DS'),
    file = 'downsampled_data.R'
    )


# Script Title: read_data.r
# Author: Massimo Cavallaro
# Description: Just read some data...
# Dependencies:

LTLA = read.table(
	"Data/LTLA.csv",
	sep = ',',
	header = T
	)
Positives = as.matrix(
	read.table(
		"Data/Number_Sgene_Positive.csv",
		header = F,
        sep = ','
        )
	)
Fails = as.matrix(
	read.table(
		"Data/Number_Sgene_Fails.csv",
        header = F,
        sep = ','
        )
	)
Total = as.matrix(
	read.table(
		"Data/Number_Sgene_Tests.csv",
        header = F,
        sep = ','
        )
	)

colnames(Positives) = as.character(0:(NCOL(Positives)-1))
colnames(Fails) = as.character(0:(NCOL(Fails)-1))
colnames(Total) = as.character(0:(NCOL(Total)-1))

rownames(Positives) = as.character(1:NROW(Positives))
rownames(Fails) = as.character(1:NROW(Fails))
rownames(Total) = as.character(1:NROW(Total))

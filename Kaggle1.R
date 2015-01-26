trainHIV <- read.csv("~/Desktop/Kaggle/trainHIV.txt")
head(trainHIV$PR.Seq)
head(trainHIV$V1)
head(trainHIV$V.1)
head(trainHIV)
length(trainHIV$PR.Seq)
which(trainHIV$Resp==0)
length(which(trainHIV$Resp==0))
head(trainHIV$RT.Seq)
library(Biostrings) # error
library(Bioconductor) # error
length(which(trainHIV$PR.Seq))

# divide data 70/15/15

# meeting 11/1/14
# Adam: run trees
# you might wanna find distances between patient's genes, run trees, rpart
# Yuri: try several different classification approaches, parametric, nonpara.
# brandon: encode sequences into a vector using aa1 = 1, etc.
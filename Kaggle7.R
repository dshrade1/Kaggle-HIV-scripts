cdata <- read.csv("~/Desktop/current projects/Kaggle/datasets/UChicago-train.csv")
names(cdata)
head(cdata)
cdons = cdata[,-c(1,2,4,5)]
dim(cdons)

#rename the codons:
newNames = c('Resp','VL.t0','CD4.t0')
for(i in 1:99){
  codNum = paste('CODON_',as.character(i),sep="")
  newNames = c(newNames,codNum)
}
names(cdons)= newNames

# How many total cases are there?
win = nrow(subset(cdons, cdons$Resp == 1) ) # 175 win
lose = nrow(subset(cdons,cdons$Resp ==0)) # 674 win
#rm(win); rm(lose)

# fisher test function: 
# try1 = 175 & try2 = 674 can be changed in case you want to use this function for a different gene.
fisherRow = function(inp,try1 = 175,try2=674){ 
  dta = matrix(c(inp[1],inp[2],try2-inp[1],try1-inp[2]),ncol=2)
  use = fisher.test(dta)
  pval = matrix(use$p.value)
  inst = matrix(inp[1]+inp[2])
  probs = cbind(pval,inst) # p value and the total number of cases
  return(probs)
}
codonPval = matrix(0,1,3)
for(i in 1:99){
  cdat = cdons[,c(1,i+3)] # get the i-th codon data + the responses
  td = table(cdat); # get the responses by the distinct codon expressions
  td = matrix(td,nrow=2) # turn it into a matrix
  expr = ncol(td) 
  for(j in 1:expr){
    inp = td[,j]
    inp  = matrix(inp, nrow = 2)
    fr = fisherRow(inp)
       codonPval = rbind(codonPval,cbind(fr,matrix(i)))
  }
}
colnames(codonPval) = c('pvalue','Num_inst', 'codon_num')
codonPval = codonPval[-1,]


# sophie's choice email:  characterize the data set
# Task 1:  identify codons of interest
# These are codons that both have SIGNIFICANT influence on determining outcomes 
# but also occur frequently in the data set (more than 40 occurences)

use <- subset(codonPval, codonPval[,2]>40 & codonPval[,1]<.05)
use

# Task 2:  our job is to look at codons 1 through 7.  First, here's an example for just codon 1:
# Note that ce = codon expression, mc = major codon, pz = proportion of outcome that are 0, pc = proportion of outcomes that are 1
table(cdons[,c(1,4)])
ce <- 'cct'
mc <- 'cct'
pz <- 613/674
po <- 150/175
fred1 <- data.frame(cod_expr =ce, maj_cod = mc, pc_zero=pz, pc_one=po)
fred1
fredF <- cbind(data.frame(t(use)),fred)
fredF

# The following code loops through to do all the codons in "use"
win <- nrow(subset(cdons, cdons$Resp == 1) ) # 175 win
lose <- nrow(subset(cdons,cdons$Resp ==0)) # 674 lose
codon.tables <- list()
fred <- data.frame()
for (i in 1:25) {
	codon.tables[[i]] <- table(cdons[,c(1,3+use[i,3])])
	ce[i] <- dimnames(codon.tables[[i]])[2][[1]][which(colSums(codon.tables[[i]])==use[i,2])]
	mc[i] <- dimnames(codon.tables[[i]])[2][[1]][which(colSums(codon.tables[[i]])==max(colSums(codon.tables[[i]])))]
	pz[i] <- codon.tables[[i]][1,which(colSums(codon.tables[[i]])==use[i,2])]/lose
	po[i] <- codon.tables[[i]][2,which(colSums(codon.tables[[i]])==use[i,2])]/win
	fred <- data.frame(cod_expr =ce, maj_cod = mc, pc_zero=pz, pc_one=po) #data frame containing only 
}
fred_all <- cbind(as.data.frame(use),fred) # data frame containing column-bound "use" and "fred" data frames

# a bit of explanation
# where [[codon1]][codon_name_dimension/not Resp][[first set]][2nd entry]
#^^iteratively add [[1]] [2] etc. to demonstrate the point
# id which dimension (colNum) we are dealing with when we find the max & ce

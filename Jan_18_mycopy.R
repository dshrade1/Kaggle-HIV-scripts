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
rm(win); rm(lose)

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

use40 = subset(codonPval, codonPval[,2]>40);
nrow(use40)
use50 = subset(codonPval, codonPval[,2]>50);
nrow(use50)
use80 = subset(codonPval, codonPval[,2]>80);
nrow(use80)

use <- subset(codonPval, codonPval[,2]>40 & codonPval[,1]<.05)

# remember: if the p value returned is small, this is a candidate codon that may distinguish the responses
# for each codon, which have small p-values?
# question: is the pvalue low enough, and is the mutation prevalent enough?  We do need to be careful about overfitting the data.
sorted_codonPval <- codonPval[order(codonPval[,1]),]
hist(log10(sorted_codonPval[,1]))
lowPval <- subset(sorted_codonPval, sorted_codonPval[,1]<0.05)
sorted_lowPval <- subset(codonPval, codonPval[,1]<0.05)

# now we have to associate these low pvalues with the codons to which they're tied.

# we could limit to HIGH num_inst and LOW pval
use40lowPval = subset(codonPval, codonPval[,2]>40 & codonPval[,1]<0.05);
nrow(use40lowPval)

x <- c(613, 150); n<-c(674,175)
prop.test(x,n,conf.level=.95)

#Look at table:  for codon1, use table(cdons[,c(1,4)])

#Example of this week's task for codon1
foo <- use[1,]
ce <- 'cct' # or we could use the names
mc <- 'cct' #
pz <- 613/674
po <- 150/175
# now create a data frame of all these things.
fred <- data.frame(cod_expr =ce, maj_cod = mc, pc_zero=pz, pc_one=po)
fred
fredF <- cbind(data.frame(use[1,]),fred)
use[1,]
u <- data.frame(use[1,])
rbind(u,fred)
cbind(u,fred)
dim(u)
u <- data.frame(use[1,])^t
u <- matrix(u,1,3)
cbind(u,fred)
fredF <- dim(cbind(u,fred))
fredF
use[2,]
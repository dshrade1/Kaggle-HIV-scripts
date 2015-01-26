
library('seqinr')
trainHIV <- read.csv("~/Desktop/Kaggle/trainHIV.txt", header=T,stringsAsFactors=T)
###dim(trainHIV)

#DNA sequences are given in 5'->3', so that's good!
#convert DNA sequences to character vectors

#create PR Sequence vector of characters using s2c
trainHIV$PR.Seq.VecStr <- lapply(trainHIV$PR.Seq, s2c)

#make a loop that converts each element of trainHIV$PR.Seq to character

# I'll just do it with 1 now
try <- as.character(trainHIV$PR.Seq[1])	# works: we get a character string
tryS2C <- s2c(try) # works: we get a vector of individual characters
translate_tryS2C <- translate(tryS2C)
#store in a matrix
aamatrix <- matrix(nrow=length(trainHIV$PR.Seq),ncol=99)
aamatrix[1,]<- translate_tryS2C

# Now I'll do it with another row
try2 <- as.character(trainHIV$PR.Seq[2])	
tryS2C2 <- s2c(try2)
translate_tryS2C2 <- translate(tryS2C2)
aamatrix[2,]<- translate_tryS2C2

# now convert all other things
# convert all PR sequences to character strings
trainHIV$PR.Char <- as.character(trainHIV$PR.Seq)
# convert all PR character strings to character vectors, and paste into a list
charList_PR <- list()
for (i in 1:length(trainHIV$PR.Char)) {
charList_PR[[i]] <- s2c(trainHIV$PR.Char[i])
}

# now we have a 1000-vector-long list, where each item on the list is a PR character vector
# these character vectors can now be handled by the translate function

PR_translated <- list()
PR_translated <- translate(charList_PR) # fails

# translate one or 2...
PR_translated[[1]] <- translate(charList_PR[[1]]) # works
PR_translated[[2]] <- translate(charList_PR[[2]]) # works
PR_translated[[3]] <- translate(charList_PR[[3]])
PR_translated[[4]] <- translate(charList_PR[[4]])
PR_translated[[5]] <- translate(charList_PR[[5]])
PR_translated[[6]] <- translate(charList_PR[[6]])
PR_translated[[7]] <- translate(charList_PR[[7]])
PR_translated[[8]] <- translate(charList_PR[[8]])
PR_translated[[9]] <- translate(charList_PR[[9]])
PR_translated[[10]] <- translate(charList_PR[[10]])

# now translate all... 
PR_translated <- list()
for (i in 1:length(charList_PR)) {
	PR_translated[[i]] <- translate(charList_PR[[i]])
} # fails
# Error in s2n(seq, levels = s2c("tcag")) : 
  # sequence is not a vector of chars
  
# try putting it in a matrix, don't even try to translate yet
PR_matrix <- matrix(nrow=1000,ncol=297)  
for (i in 1:length(charList_PR)) {
	PR_matrix[i,] <- t(as.matrix(charList_PR[[i]]))
} # fails
# Error in PR_matrix[i, ] <- t(as.matrix(charList_PR[[i]])) : 
#  number of items to replace is not a multiple of replacement length

for (i in 1:length(charList_PR)) {
	for (j in 1:297) {
		PR_matrix[i,j] <- t(as.matrix(charList_PR[[i]]))
	}
	}
	

# put it all into a matrix so we can apply the translate functions
PR_matrixdim1 <- lapply(charList_PR, as.matrix) # works
PR_vectors <- lapply(charList_PR, as.vector) # works



charMatrix_PR <- 
# doesn't work!

PR_translated <- sapply(charList_PR,translate) #fails
PR_translated <- apply(charList_PR,translate)
PR_vectors <- lapply(charList_PR, as.matrix) # works
PR_vectors <- lapply(charList_PR, as.vector) # works
PR_translated <- lapply(PR_vectors, translate) #fails

prMatrix <- matrix(nrow=length(charList_PR), ncol=length(charList_PR[[1]]))

charListM <- matrix(nrow=297,ncol=1000)

# first convert all sequences to matrices, so we have a list of character vectors
charList2_PR <- list()
charList2_PR[[1]] <- t(as.matrix(charList_PR[[1]])) # this works

for (k in 1:length(charList_PR)) {
	charList2_PR[[i]] <- t(as.matrix(charList_PR[[i]]))
} # doesn't work

# next, store all character vectors in a matrix.

for (i in 1:length(charList_PR[[1]])) {
	for (j in 1:length(charList_PR)) {
		charListM[i,j] <- as.matrix(charList_PR[[i]][j])
	}
	charListM[i,j] <- as.matrix(charList_PR[[]][j])
}
# doesn't work


for (j in 1:length(charList_PR[[1]][1])) {
	for (i in 1:length(charList_PR[[1]])) {
		prMatrix[i,j] <- charList_PR[[i]][j]
	}
}  # doesn't work

aamatrix <- matrix(nrow=length(trainHIV$PR.Seq),ncol=99)
PR_translated <- list()
for (i in 1:length(charList_PR)) {
	PR_translated[i] <- translate(as.matrix(charList_PR[i], nrow=),ambiguous=T)
}

for (i in 1:length(charList_PR))  {
for (j in 1:length(charList_PR[[1]])) {
	PR_translated[[i]] <- translate(charList_PR[[i]])
	aamatrix[i,j] <- translate(charList_PR[[i]][j])
}	
}




aalist <- list()
aalist[[1]] <- 

}

# have to loop through PRchar to use s2c command
# charMatrix_PR <- matrix(nrow=length(trainHIV$PR.Seq),ncol=nchar(trainHIV$PRchar[1]))
# for (i in 1:length(trainHIV$PRchar)) {
	# charMatrix_PR[i,] <- matrix(s2c(trainHIV$PRchar[i]),nrow=1,ncol=length(trainHIV$PRchar[1])
# }



trainHIV$PRs2c <- s2c(trainHIV$PRchar)

# Cool.  Now, what can I do with these?
which(aamatrix=="P") # this gives me indices as answers.  R calculates indices column-wise.
which(aamatrix[,1]=="P") #I ask: which of the columns = P?  I get the columns that = P, (1 & 2).
# what is the most frequently-occurring amino acid for that locus, and how often does it occur?


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
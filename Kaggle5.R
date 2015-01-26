# get the other stuff aligned using XinQi's program?
require(Biostrings)
train <- read.csv("~/Desktop/Kaggle/Full_Training_Aligned_PR_Codons.csv", header=T, stringsAsFactors=F ,sep=",")
train2 <- read.csv("~/Desktop/Kaggle/Full_Training_Aligned_PR_Codons2.csv", header=T, stringsAsFactors=F ,sep=",")

rt <- train$RT.Seq
RT <- DNAStringSet(rt)
RT_translated <- translate(RT, if.fuzzy.codon="solve")
#Warning messages:
# 1: In .Call2("DNAStringSet_translate", x, skip_code, dna_codes[codon_alphabet],  :
  # in 'x[[405]]': last 2 bases were ignored
# 2: In .Call2("DNAStringSet_translate", x, skip_code, dna_codes[codon_alphabet],  :
  # in 'x[[665]]': last 2 bases were ignored
aaRT <- RT_translated
# excellent!!  Now, how can I plot this?  try a few things...

RT_matrix <- matrix(data=NA,nrow=length(aaRT),ncol=max(nchar(aaRT)))  
RT_char <- as.character(aaRT)
RT_char_split <- strsplit(RT_char,"")
RT_char_split[[1]][1] # = P

for (i in 1:1000) {
	for (j in 1:494) {
		RT_matrix[i,j] <- RT_char_split[[i]][j]
	}
}
# max nchar for aaRT is 494
# min nchar for aaRT is 193

# OK, but what about for the ALIGNED data??

#do it once
# aPR = "aligned Protease"

train2 <- gsub("---","nnn",train)

aPR <- list()
aPR_translated <- list()
aPR[[1]] <- DNAStringSet(train2[,7])

aPR_translated[[1]] <- translate(aPR[[1]], if.fuzzy.codon="solve")
# works!

# now for all of them
for (i in 1:99) {
	aPR[[i]] <- DNAStringSet(train2[,i+6])
	aPR_translated[[i]] <- translate(aPR[[i]], if.fuzzy.codon="solve")
	
}

# works

aPR_matrix <- matrix(data=NA,nrow=1000,ncol=99)

for (i in 1:1000) {
	for (j in 1:99) {
		aPR_matrix[i,j] <- as.character(aPR_translated[[j]][i])
	}
}
# took a long time, but works

typeof(aPR_translated)
class(aPR_translated)
typeof(aPR_translated[[1]])
class(aPR_translated[[1]])
typeof(aPR_translated[[1]][1])
class(aPR_translated[[1]][1])

# My question is:  what's the biggest entity on which I can perform "as.character"?

#now, I want to put this into a file

write.table(aPR_matrix, file = paste("~/Desktop/Kaggle","aPR_matrix",sep="/"),sep=",")

# Conclusion:  I literally did "find and replace" in the CSV file to replace "---" with "nnn"
# I created a matrix of the ALIGNED & Translated values for the protease matrix (thanks XinQi!)
# I saved this matrix to a csv-like file
train_cols1256 <- train[,c(1,2,5,6)]
train_cols1256m <-as.matrix(train_cols1256)
write.table(train_cols1256m, file = paste("~/Desktop/Kaggle","noDNA_matrix",sep="/"),sep=",")
data_aligned_PR_only_matrix <- cbind(train_cols1256m,aPR_matrix)
write.table(data_aligned_PR_only_matrix,file = paste("~/Desktop/Kaggle", "data_aligned_PR_matrix", sep="/"),sep=",")
RT_matrix[1:100,472:494]
data_PR_aligned_with_RT_matrix <- cbind(data_aligned_PR_only_matrix,RT_matrix)
write.table(data_PR_aligned_with_RT_matrix, file = paste("~/Desktop/Kaggle","data_PR_aligned_with_RT_matrix",sep="/"),sep=",")

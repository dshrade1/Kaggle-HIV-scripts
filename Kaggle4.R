#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
hiv[5:dim(hiv)[2]]# works!!!

length(which(PR_matrix[,1]=="P"))
length(which(PR_matrix[,2]=="Q"))

# http://svitsrv25.epfl.ch/R-doc/library/Biostrings/html/translate.html
# https://web.stanford.edu/class/bios221/labs/biostrings/lab_2_biostrings.html

# read the PR sequence as a string set
# read individual codons and translate them to proteins using seqinr?

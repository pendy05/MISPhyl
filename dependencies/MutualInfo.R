#Load seqinr library to access read.alignment function for reading alignment file in fasta format
library(seqinr)


#receive cpu num from python script
options(echo=TRUE)
args<-commandArgs(trailingOnly=TRUE)
folderpath<-args[1]
cpuNum<- strtoi(args[2])

gene_list <-list.files(path=folderpath)
gene_list<-paste(folderpath,gene_list,sep="")


############################################
#Function for computing mutual information #
############################################

minf <- function(dat) {
#Convert dat into a list
X <- lapply(dat$seq, function(j)
factor (unlist(strsplit(j,"")), levels=c("a","g","c","t","-") )
)
#Joint probability distribution for the bases
pxy <- table (unlist(X)) / (length(X) * length(X[[1]]))
#Row marginal probability distribution
px <- lapply(X, function(j)
table(j) / length(j) )
#Column marginal probability distribution
V <- vector("list",length(X[[1]]))
coly <- lapply(X, function(j) {
mapply(function(k)
table(j[k]), 1:length(X[[1]]) )}
)
H <- coly[[1]]
for(k in 2:length(coly)){
	H <- H + coly[[k]]
	}
H <- H / colSums(H)
#Putting everything together
base <-c("a","g","c","t","-") 
mut_inf <- 0
for(y in 1:length(X)){
mut_inf <- mut_inf + mapply(function(s) {
nt <- which(base == X[[y]][s])
pxy[nt] * log(pxy[nt] / H[nt,s] / px[[y]][nt]) }
, 1:length(X[[1]]))
}
#normalised mut_inf
return ( sum(mut_inf) / length(X[[1]]) )
}

####################################
#Computes gap percentage for MSA   #
####################################

gap_pct <- function(dat){

X <- lapply(dat$seq, function(j)
factor (unlist(strsplit(j,"")), levels=c("a","g","c","t","-") )
)

sum ( mapply(function(j)
lapply(X,table)[[j]][5],
1:length(X)) ) / (length(X)*length(X[[1]])) * 100

}




#Assuming that the MSA for all the genes of interest have been obtained

gene_names <- unlist(strsplit(gene_list, ".fasta"))

#Running minf on all these MSA is simple

date()
T <- vector("list",length(gene_list))

for(i in 1:length(T)){
T[[i]] <- read.alignment(file=gene_list[i], "fasta")
}

if (cpuNum <= 1){#serial
gappiness <- mapply(function(j)
gap_pct(T[[j]]), 1:length(T))

MI_genes <- mapply(function(j)
minf(T[[j]]), 1:length(T)
)

} else { #parallelization
gappiness <- parallel::mcmapply(function(j) gap_pct(T[[j]]), 1:length(T), mc.cores=cpuNum)

MI_genes <- parallel::mcmapply(function(j) minf(T[[j]]), 1:length(T), mc.cores=cpuNum)
}

U <- cbind(MI_genes, gappiness)
U <- as.data.frame(U)
rownames(U) <- gene_names
write.csv(U, file="MI_genes.csv")
date()

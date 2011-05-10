############################################################################################################
#Program:	ParallLogicReg
#Programer:	Unitsa Sangket
#Date:		June 18, 2011
#Objective:	To execute logic regression analyses using parallel computing
############################################################################################################

######################################### We're in the parent. #############################################

"ParallLogicReg" <- function(infile, resp, begin=1, end, nperm=20, niter=20){ 

library(LogicReg)

###check the important parameters, which are infile and resp
if (missing(infile)) 
    stop("Please enter an input file")

if (missing(resp)) 
    stop("Please enter resp")

###read data from "infile" file
gene01=read.table(infile,sep="\t",header=T)

###find chromosome number
Chr = gene01[2,2]  

###if "end" is null, it will be equaled number of genes in the "infile" 
if (missing(end)) 
	end=gene01[length(gene01[,1]),5]

glist=as.numeric(levels(as.factor(gene01[,5])))
freq=cbind(glist,table(gene01[,5])/2)   ###each SNP has 2 rows
summary(freq[,2])   

###genes with 2+SNP only
gidlist=freq[freq[,2]>=2,1]

gid_start=gidlist[begin]	###the gid  to start
gid_end=gidlist[end]			###the gid to end

gidlist_here=gidlist[begin:end]	###the gid in this part

gene01=gene01[gene01[,5] %in% gidlist_here,]

len_gid = length(gidlist_here)   ###find length of gidlist_here 

###I have the permuted response outside, and a seed is specified
###so in for each gene, we have the same permuted Y values ################

perYY=matrix(NA,nperm,length(resp))
set.seed(2010)
for(i in 1:nperm){
perYY[i,]=sample(resp)   ###each column is a permutated response variable
}

###### do this way, we always know the permuted y for the nperm


perm_dev=list()  ###each lsit here is a vecotr of size nperm+1--original minimum dev and the permuted min dev among the niter

ngenes_skip = begin - 1  ###ngenes_skip = the number of genes will be skipped

############################ permutation function for some gene ##########################
"genes_nperm" <- function(start=0,stop=0,taskX=0){
	### start = the first gene that will be run
	### stop = the last gene that will be run
	### taskX = the task number that is running

	library(LogicReg)

	load("variables.Rdata")

	ii = 0 ### number of gene have run

	### write the progress to a progress file 
	ngenes_all = stop - start + 1 # number of gene to run
	progress=paste("task", taskX, ": gene id  = ", start, ", ngenes = ", ii, "/", ngenes_all, ", progress = ", (ii/ngenes_all)*100, "%",sep="")
	progress_file_name = paste("progress_task", taskX,".txt", sep="")
	write(progress, file = progress_file_name)
	
	################# loop for each gene ###########	
	for(gg in (start-ngenes_skip):(stop-ngenes_skip)){
		
		gid=gidlist_here[gg]
		path=gene01[gene01[,5]==gid,]
		print(dim(path))

		nsnp=table(path[,5])/2 ###number of snps for the gene
		leaf=min(nsnp,5)   ###number of leaves is the minimum of these two

		set.seed(gg)   ###each gene starts with a unique seed, so we can get he model for the best dev

		ori=Inf
		Dev_per=rep(Inf,nperm)   ###the minimum dev in each permutation , vecotr of size nperm

		
		################# start permutaitons ###########		
		for(j in 1:nperm){ ###j is the permutation number
		
		###for each permutation, run  iterations, save the minimum dev among the  iterations.
			
			################# start iterations ###########	
			for(i in 1:niter){
			### the original Y for the original, do iterations, get the best dev, and the Phat corresponding to it 
				if(j==1){  
					res=logreg(resp=resp,bin=t(path[,-c(1:5)]),type=3,select=1,ntrees=2,nleaves=leaf)
					if(res$model$score<ori) {
						ori=res$model$score 
						}
				}

			###for the permuted, do iterations, get the best dev, and the Phat corresponding to it 
				YY=perYY[j,]    ###choose the jth permuted response variable 
				res=logreg(resp=YY,bin=t(path[,-c(1:5)]),type=3,select=1,ntrees=2,nleaves=leaf)
	
				if(Dev_per[j]>res$model$score){
				Dev_per[j]=res$model$score
				}
		  }################# end of iteration ############# 
			
		 ### write the progress to a progress file 
		 ngenes_all = stop - start + 1 # number of gene to run
		 progress=paste("task",taskX , ": gene id  = ", gg+ngenes_skip, ", ngenes = ", ii, "/", ngenes_all, " nperm = ", j, ", progress = ",(ii+(j/nperm))/ngenes_all*100, "%",sep="")
		 progress_file_name = paste("progress_task", taskX,".txt", sep="")
		 write(progress, file = progress_file_name)
		 
		}################# end of permutation ############# 
		
		ii = ii + 1  ### number of gene have already been run; ii = gg-start+1

		perm_dev[[ii]]=c(ori,Dev_per) ### example: the dev vecore of size 21, 1 original +20 permuted
		### edit index of the list to support combining these lists from compute-nodes
		
	} ################# end of loop for each gene ###########	

	### write progress to a progress file 
	ngenes_all = stop - start + 1 # number of gene to run
	progress=paste("task",taskX ,": gene id  = ", gg+ngenes_skip,", ngenes = ", ii, "/", ngenes_all, " nperm = ", j, ", progress = ",(((ii-1) + (j/nperm))/ngenes_all)*100, "%",sep="")
	progress_file_name = paste("progress_task", taskX,".txt", sep="")
	write(progress, file = progress_file_name)
	
	out = list()
	out$perm_dev = perm_dev

	file.remove(progress_file_name)### remove a progress file
	return(out) ### return output from this function to frontend-node

}
############################ end of permutation function ##########################


############################ parallel imputation start here ########################

library(ParallLogicReg)
#source("Job_management.R")

save(list = ls(),file="variables.Rdata")
output.p <- Job_management()  ### output from imputation

############################ parallel imputation stop #######################

perm_dev = output.p[[1]]
temp_dev = output.p[[1]]

save(perm_dev,file = "out.Rdata")


###write all the results 
###write the deviance 
length(temp_dev)  ###this tells how many gene have run
dev_to_write=matrix(NA,nperm+1, length(temp_dev))
for(aa in 1: length(perm_dev)){	
	dev_to_write[,aa]=temp_dev[[aa]]
}

devout=paste("dev",Chr,".txt",sep="")
write.table(dev_to_write,devout,col.names=F,sep="\t")

Chr 
length(perm_dev)

return("Please see the dev output file")

} ### end of ParallLogicReg function
############################################################################################################
#Program:	SeqLogicReg
#Programer:	Unitsa Sangket
#Date:		June 18, 2011
#Objective:	To execute logic regression analyses using only one CPU.
############################################################################################################

"SeqLogicReg" <- function(){
	load("variables.Rdata")

	start=begin  ### start = the first gene that will be run
	stop=end		### stop = the last gene that will be run
	taskX=1		### taskX = the task number that is running

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

		nsnp=table(path[,5])/2 ##number of snps for the gene
		leaf=min(nsnp,5)   ###number of leaves is the minimum of these two

		set.seed(gg)   ###each gene starts with a unique seed, so we can get he model for the best dev


		ori=Inf
		Dev_per=rep(Inf,nperm)   ###the minimum dev in each permutation , vecotr of size nperm

		
		#################start permutaitons ###########		
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
				res=logreg(resp=resp,bin=t(path[,-c(1:5)]),type=3,select=1,ntrees=2,nleaves=leaf)
	
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
	save(perm_dev,file = "out.Rdata")

	### write all the results 
	###write the deviance 
	temp_dev = perm_dev

	length(temp_dev)  ###this tells how many gene had run
	dev_to_write=matrix(NA,nperm+1, length(temp_dev))
	for(aa in 1: length(perm_dev)){	
		dev_to_write[,aa]=temp_dev[[aa]]
	}

	devout=paste("dev",Chr,sep="")
	write.table(dev_to_write,devout,col.names=F,sep="\t")

	file.remove(progress_file_name)### remove a progress file

	return("Please see the dev output file")

} # end of SeqLogicReg function

############################################################################################################
#Program:	Job_management
#Programer:	Unitsa Sangket
#Date:		June 18, 2011
#Objective:	To control tasks and compute-nodes
############################################################################################################

"Job_management" <- function(){ 

#Initialize MPI
library("Rmpi")

#Notice we just say "give us all the slaves you've got."
mpi.spawn.Rslaves()

if (mpi.comm.size() < 2) {
    # print("More slave processes are required.")
	 # use logic regression sequential version
	 library(ParallLogicReg)
	 # source("SeqLogicReg.R")
	 output=SeqLogicReg()
    mpi.quit()
    }

.Last <- function(){
    if (is.loaded("mpi_initialize")){
        if (mpi.comm.size(1) > 0){
            print("Please use mpi.close.Rslaves() to close slaves.")
            mpi.close.Rslaves()
        }
        print("Please use mpi.quit() to quit R")
        .Call("mpi_finalize")
    }
}


########## Function the slaves will call to perform a validation on the 
# fold equal to their slave number.
# Assumes: fold,foldNumber
foldslave <- function(){
    # Note the use of the tag for sent messages: 
    #     1=ready_for_task, 2=done_task, 3=exiting 
    # Note the use of the tag for received messages: 
    #     1=task, 2=done_tasks 
    junk <- 0 
    done <- 0 

    while (done != 1) {
        # Signal being ready to receive a new task 
        mpi.send.Robj(junk,0,1) 

        # Receive a task 
        task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
        task_info <- mpi.get.sourcetag() 
        tag <- task_info[2] 

        if (tag == 1) {
			
		#************** 3.task for compute-nodes  ******************************
		foldNumber=task$foldNumber

		start = task$start    ###set start gene id 
		stop = task$stop    ###set stop gene id
		taskX = foldNumber
		
		output=genes_nperm(start=start,stop=stop,taskX=taskX)
		
		results <- list(foldNumber=foldNumber,output=output) 

		#***************** end task for compute-nodes ***************************

		mpi.send.Robj(results,0,2)
        }
        else if (tag == 2) {
            done <- 1
        }
        # We'll just ignore any unknown messages
    }

    mpi.send.Robj(junk,0,3)
}


#*********************************** 1.separate task   ******************************

#check number of available processors
npro = mpi.universe.size() 

#load all variables from main program
load("variables.Rdata")

##### separate data 
nrounds = len_gid     #total of genes, let's assign nround is number of genes

npro = nrounds

nrounds_task = floor(nrounds/npro)  #number of genes to be run on a task
pointer = 0

#Create task list for each processor
tasks <- vector('list')

#first process
i = 1
start = begin
stop = begin + (i*nrounds_task)-1
tasks[[i]] <- list(foldNumber=i, start = start, stop = stop) #"stop" is stop gene id

#next processes
for (i in 2:(npro-1)) {
	 start = stop+1
	 stop = start+nrounds_task-1 
	 tasks[[i]] <- list(foldNumber=i, start = start, stop = stop) 

}
#last process
if (nrounds_task > 0){
	i=i+1
	start = stop+1
	stop = end
	tasks[[i]] <- list(foldNumber=i, start = start, stop = stop) 
}

# initial results
results <- vector('list')
for (i in 1:npro) {
    results[[i]] <- list(output=i)
}

#********************************** finish task separation  *****************************

# Now, send the data to the slaves

# Send the function to the slaves
mpi.bcast.Robj2slave(foldslave)

#******************************** 2. distribute task   ******************************
# Call the function in all the slaves to get them ready to
# undertake tasks


# send these variables to compute-node
mpi.bcast.cmd(foldslave())

#********************************* finish task distribution ****************************

junk <- 0 
closed_slaves <- 0 
n_slaves <- mpi.comm.size()-1 


while (closed_slaves < n_slaves) { 
    # Receive a message from a slave 
    message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
    message_info <- mpi.get.sourcetag() 
    slave_id <- message_info[1] 
    tag <- message_info[2] 

    if (tag == 1) { 
        # slave is ready for a task. Give it the next task, or tell it tasks 
        # are done if there are none. 
        if (length(tasks) > 0) { 
            # Send a task, and then remove it from the task list 
            mpi.send.Robj(tasks[[1]], slave_id, 1); 
            tasks[[1]] <- NULL 
            } 
        else{ 
            mpi.send.Robj(junk, slave_id, 2) 
            } 
    } 
    else if (tag == 2) {

		#************************ 4.store result *************************
		# The message contains results. Do something with the results. 
		# Store them in the data structure
		results[[message$foldNumber]] = message$output
		#******************************************************************

    } 
    else if (tag == 3) { 
        # A slave has closed down. 
        closed_slaves <- closed_slaves + 1 
    } 

} 

mpi.close.Rslaves()

#*************************** 5.combine result ************************

r_sep = results	# r_sep = results form every compute-nodes that not yet combine
r_com  =	r_sep[[1]]	# r_com = results that are combined   

n_r_sep = names(r_sep[[1]])

#check structure of result
if (is.list(r_sep)){
	# combine results
	for (i in 2:npro) {
		for (j in 1:length(n_r_sep)){
				r_com[[n_r_sep[j]]] = c(r_com[[n_r_sep[j]]], r_sep[[i]][[n_r_sep[j]]])
		} 
	}
	
}else {
	message_error = paste("Error: structure of result from ", "genes_nperm", " doesn't list or data.frame", sep = "")
	stop(message_error)
}

# remove "variables.Rdata" 
file.remove("variables.Rdata")

return(r_com)

#############

}# End of Job_management function

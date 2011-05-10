\name{ParallLogicReg}
\alias{ParallLogiReg}
\title{Logic regression analyses using parallel computing}
\description{
	This function can be used in many situations, especially identification of SNP interactions
}
\usage{
output = ParallLogicReg(infile,resp,begin,end,nperm,niter) 
}
\arguments{
  \item{infile}{A file name that contains data, which will be analyzed.}
  \item{resp}{vector with the response variables. See more detail in LogicReg manual. }
  \item{begin}{The first gene id will be run.}
  \item{end}{The last gene id will be run.}
  \item{nperm}{number of permutations}
  \item{niter}{number of iterations}
}

\details{
An example of infile:
\tabular{llllllllll}{
snp \tab chr \tab pos \tab gene-symbol \tab gene-id \tab V1 \tab V2 \tab V3 \tab V4 \tab V5\cr
rs1 \tab 1 \tab 12691937 \tab AADACL3 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1\cr
rs1 \tab 1 \tab 12691937 \tab AADACL3 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1\cr
rs31 \tab 1 \tab 12717433 \tab AADACL3 \tab 1 \tab 1 \tab 0 \tab 1 \tab 1 \tab 1\cr
rs31 \tab 1 \tab 12717433 \tab AADACL3 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1\cr
rs2 \tab 1 \tab 12668821 \tab AADACL4 \tab 2 \tab 0 \tab 1 \tab 1 \tab 1 \tab 1\cr
rs2 \tab 1 \tab 12668821 \tab AADACL4 \tab 2 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1\cr
rs58 \tab 1 \tab 12633603 \tab AADACL4 \tab 2 \tab 0 \tab 0 \tab 0 \tab 0 \tab 0\cr
rs58 \tab 1 \tab 12633603 \tab AADACL4 \tab 2 \tab 0 \tab 0 \tab 1 \tab 0 \tab 0\cr
rs59 \tab 1 \tab 12651828 \tab AADACL4 \tab 2 \tab 1 \tab 1 \tab 1 \tab 1 \tab 1
}


where,

snp = snp id\cr
chr = chromosome number\cr
pos = position of the snp on chromosome\cr 
gene-symbol = gene symbol or gene name\cr
gene-id = gene id\cr
V1,V2,V3,V4,V5,... = snp of control or case


Note: One snp has two data rows, and each field must be separated by tab. 
}

\author{Unitsa Sangket and Qi Liu
}

\examples{
#--------------------------------------------------------------------------
#An example run on Sun Grid Engine.
#You have to create two files, which are a ".sh" file, and a ".R" file.
#for more information please visit http://math.acadiau.ca/ACMMaC/Rmpi/submitting.html
#--------------------------------------------------------------------------

#File 1 (R_script.sh):
###########################################################################
#!/bin/bash

# Run in the current directory
#$ -cwd

#$ -j y

#$ -V

# Run using bash
#$ -S /bin/bash

# Set the number of processors
# For example, 2 means one processor is master and slave1, 
# and the rest of processors are slaves.
#$ -pe lam 2

# Run the job. 
# lamrun -np 1 R  --slave CMD BATCH R_script.R R_script_sh.Rout

###########################################################################
#!!! remove "#" out from "lamrun" command when run the file


#File 2 (R_script.R):
###########################################################################
library(ParallLogicReg)

resp=c(rep(0,2935),rep(1,1745))	# number of controls = 2,935; number of cases = 1,745
 
nperm=20			# number of permutations
niter=20				# number of iterations
begin=1				# the first gene id that will be run
end=10				# the last gene id that will be run
in=("input.txt")   	# data will be run

output = ParallLogicReg(infile=in,resp=resp,begin=begin,end=end,nperm=nperm,niter=niter) 

mpi.quit(save="no")

##########################################################################

#--------------------------------------------------------------------------
#You can summit the job using this command:
#qsub R_script.sh

#You can see the progress on R_script_sh.Rout and progress files.
#--------------------------------------------------------------------------
}


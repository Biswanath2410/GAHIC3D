#Publications that use results obtained from this software please include a citation of the paper:
#Zou, C., Zhang, Y., Ouyang, Z. (2016) HSA: integrating multi-track Hi-C data for genome-scale reconstruction of 3D chromatin structure. Submitted.

args <- commandArgs()

args = args[-(1:match("--args", args))]

M=length(args)

if (M<2){
	print('Input should include at least a contact map and outputfilename!')
	q()
}

if(args[M] %in% c("0","1")){
	if (M<3){
	print('Input should include at least a contact map and outputfilename!')
	q()
    }
	outfile=args[M-1]
	mak=as.numeric(args[M])
	Iscovfile=args[M-2]
	K=M-2
    IniS=NULL
}else{
    
    if(args[M-1] %in% c("0","1")){
        
    if (M<4){
    print('Input should include at least a contact map and outputfilename!')
    q()
    }
    
    outfile=args[M-2]
    mak=as.numeric(args[M-1])
    Iscovfile=args[M-3]
    IniS=args[M]
    K=M-3
    }else{
    print('Missing or incorrect value of Markov_Indicator! Markov_Indicator will be set to 0. Error might occur due to wrongly formatted command. Please check your command line format.')
    IniS=args[M]
    outfile=args[M-1]
    mak=0
    Iscovfile=args[M-2]
    K=M-2
    }
}

Iscovfile= !(Iscovfile=="0")

if(Iscovfile){	
	K=K/2
	lscov0=vector("list",K)
    for(i in 1:K){
	covfile=read.table(file=args[K+i],header=F)
	lscov0[[i]]=lapply(3:dim(covfile)[2], function(x) covfile[,x]%*%t(covfile[,x]))
    }
}else{
K=K-1
lscov0=0
}

print(is.numeric(lscov0))

lsmap0=vector("list",K)
for(i in 1:K){
	mat=read.table(file=args[i],header=F)
	lsmap0[[i]]=as.matrix(mat)
}

library(MASS)
source("cstruct1.R")
if(is.null(IniS)){
out=fmain(lsmap0,lscov0,outfile,300,150,50,50,0.001,0,0,mak)
}else{
S=read.table(file=IniS,header=F)
S=as.matrix(S)
out=fmain(lsmap0,lscov0,outfile,300,150,50,50,0.001,0,0,mak,S)
}
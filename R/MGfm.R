#' @title Multi-Group Single-Trait Fine-Mapping credible sets
#' @param stfm.list List of data.frame objects giving model PP from single-trait fine-mapping in each group
#' @param isnps Vector of variants that are present in at least one group's GWAS
#' @param N Vector of sample sizes in same order as groups in stfm.list
#' @param cred Level for credible set, default 0.99 
#' @param cpp cumulative posterior probability threshold for selecting top models; default 0.99
#' @return List of three objects: summary: data.frame of variants that belong to the credible set and their multi-group MPP (mgMPP), details: list of mgMPP for all variants and mgPP for the top multi-group models, cred:credible set level
#' @author Jenn Asimit
#' @export
MGfm <- function(stfm.list,isnps,N,cred=.99,cpp=0.99) {

A <- length(stfm.list)
nsnps <- length(isnps) 	
	
 STR <- PP <- vector("list",A)   
 for(i in 1:A) {
	STR[[i]] <- rownames(stfm.list[[i]])
	PP[[i]] <- stfm.list[[i]][,1]
#	ind0 <- which(STR[[i]]=="1")
#	if(length(ind0)>0) {STR[[i]] <- STR[[i]][-ind0]; PP[[i]] <- PP[[i]][-ind0]}
 	}

 bestmod.thr <- vector("list", A)
    for (i in 1:A) {
    	tmp <- PP2snpmod(data.frame(PP=PP[[i]],str=STR[[i]]))
        bm <- best.models.cpp(tmp, cpp.thr = cpp)
        bestmod.thr[[i]] <- bm$models
    }
    STR <- lapply(bestmod.thr, "[[", "str")
    PP <- lapply(bestmod.thr, "[[", "PP")
   
    if (A < 2) 
        stop("Need at least 2 groups")
    if (is.null(names(STR))) 
        names(STR) <- paste0("group", seq_along(STR))
    qt <- names(STR)

indT <- 1:A
rmT <- c()
for(i in 1:A) {
 check <- 0
 if("1" %in% STR[[i]]) check <- PP[[i]][which(STR[[i]]=="1")] >= 0.9
 if(check) rmT <- c(rmT,i)
}
if(length(check) > 0) indT <- setdiff(indT,rmT)
A <- length(indT)

if(A>1) {
STR <- STR[indT]
PP <- PP[indT]
N <- N[indT]

if(A==2) { 	osnps <- intersect(unlist(strsplit(STR[[1]],"%")),unlist(strsplit(STR[[2]],"%")))
	osnps <- setdiff(osnps,c("0","1"))
  for(j in 1:A) {
  check <- sapply(strsplit(STR[[j]],"%"),function(x) any(osnps %in% x))
  STR[[j]] <- STR[[j]][which(check)]
  PP[[j]] <- PP[[j]][which(check)]
  }
}

  osnp <- c() 	
  for(i in 2:A) {
   for(j in 1:(i-1)) {
    check <- intersect(unlist(strsplit(STR[[i]],"%")),unlist(strsplit(STR[[j]],"%")))
    check <- setdiff(check,c("0","1"))
    osnp <- union(osnp,check)
   }
  }    
  
  if(length(osnp)>0) {
  ## calculate model sizes and adjust each input log PP (PP' instead of PP)
    SS <- lapply(STR,strsplit,"%")
    usnps <- sort(unique(unlist(SS)))
    nsnpspermodel <- lapply(SS,function(x) sapply(x,length))
    for(i in seq_along(STR)) {
        wh <- which(STR[[i]] %in% c("0","1"))
        nsnpspermodel[[i]][wh] <- 0
        names(nsnpspermodel[[i]]) <- STR[[i]]
    }
   
 	Ntot <- sum(N) 
    PP.orig <- PP
    logPP <- vector("list",A)		# logPP adjusted by eta
    for(i in seq_along(STR)) {       
        eta <- 0.5 * nsnpspermodel[[i]] * log(N[i]/Ntot) 
        logPP[[i]] <- log(PP.orig[[i]]) + eta        
    }
PPn <- logPP

	
	message(A," groups")
 
    if(A==2) out <- TA2flashfm1(PPn,nsnpspermodel,nsnps,cred)
    if(A==3) out <- TA3flashfm1(PPn,nsnpspermodel,SS,nsnps,isnps,cred)  
    if(A==4) out <- TA4flashfm1(PPn,nsnpspermodel,SS,nsnps,isnps,cred)  
    if(A==5) out <- TA5flashfm1(PPn,nsnpspermodel,SS,nsnps,isnps,cred)  
    if(A==6) out <- TA6flashfm1(PPn,nsnpspermodel,SS,nsnps,isnps,cred)  
    
   
 	tampp <- sort(out$taMPP,decreasing=T)
    tapp <- out$taPP[order(out$taPP[,1],decreasing=T),,drop=F]
    sg <- data.frame(mgMPP=tampp[out$cs])
	sg <- sg[order(sg[,1],decreasing=T),,drop=FALSE]
	} else{
	sg <- NULL; tampp <- NULL; tapp <- NULL
	} 
	} else{
	sg <- NULL; tampp <- NULL; tapp <- NULL	
	}
    
return(list(summary=sg,details=list(mgMPP=tampp,mgPP=tapp),cred=cred))	

}


#' @title Wrapper for credible sets from Multi-Group Single-Trait Fine-Mapping with JAM 
#' @param gwas.list List of A data.frame objects, where A is the number of groups; gwas.list\[\[i\]\] is a data.frame for group i with 3 columns named: rsID, beta, EAF
#' @param corX.list List of A data.frame objects, where A is the number of groups; corX.list\[\[i\]\] is the SNP correlation matrix for group i 
#' @param Nall Vector of length A; Nall\[i\] is the (effective) sample size for group i
#' @param save.path Path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1").
#' @param cpp cumulative posterior probability threshold for selecting top models; default 0.99
#' @param cred Level for credible set, default 0.99 
#' @param maxcv starting value for maximum number of causal variants
#' @param maxcv_stop maximum value to consider for maximum number of causal variants; maxcv_stop >= maxcv
#' @param NCORES number of cores for parallel computing; recommend NCORES=A, but if on Windows, use NCORES=1
#' @param jam.nM.iter in millions, number of iterations to use in JAM; defailt 1 (1 million)
#' @param extra.java.arguments A character string to be passed through to the java command line. E.g. to specify a
#' different temporary directory by passing "-Djava.io.tmpdir=/Temp".
#' @return List consisting of two objects: CSsummary = List of one data.frame for each trait; each trait data.frame gives the variants in the multi-group  credible set for the trait, the mgMPP, pooled MAF, proportion of studies 
#' that contain the variant, names of studies that contain the variant
#' CSdetail = \[\[1\]\] a list of multi-group  credible sets (variants and their mgMPP), one for each trait; \[\[2\]\] details = for each trait, list of all multi-group  models and variants and their mgPP, mgMPP
#' @author Jenn Asimit
#' @import R2BGLiMS
#' @export
MGFMwithJAM <- function(gwas.list, corX.list,  Nall, save.path,cpp=0.99,cred=0.99, maxcv=1, maxcv_stop = 20,NCORES=1,jam.nM.iter=1,extra.java.arguments=NULL) {


maxcv_autocheck = TRUE
A <- length(gwas.list)
beta.list <- lapply(gwas.list,function(x) {b <- x[,"beta"]; names(b) <- x[,"rsID"]; b})
raf.list <- lapply(gwas.list,function(x) {b <- x[,"EAF"]; names(b) <- x[,"rsID"]; b})

for(i in 1:A) {
 ksnp <- intersect(names(raf.list[[i]]),colnames(corX.list[[i]]))
 raf.list[[i]] <- raf.list[[i]][ksnp]
 beta.list[[i]] <- beta.list[[i]][ksnp]
 corX.list[[i]] <- corX.list[[i]][ksnp,ksnp]
}

 if(!is.null(names(gwas.list) )) {
   studies <- names(gwas.list)
  } else { studies <- paste0("group",1:A)
  		names(gwas.list) <- studies
  }
  
ivec <- vector("list",A)
for(i in 1:A) {
 ivec[[i]] <- i
 dir.create(paste0(save.path,"/a",i))
} 
jamfn <- function(i,beta.list,corX.list,raf.list,ybar,Vy,Nall,save.path,maxcv,maxcv_stop,maxcv_autocheck,extra.java.arguments) {
 JAMcor.tries.maxcv(beta.list[[i]], corX.list[[i]], raf.list[[i]], ybar=0, Vy=1, 
  Nall[i], paste0(save.path,"/a",i), maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)}

stfm <- parallel::mclapply(ivec,jamfn,beta.list,corX.list,raf.list,ybar=0, Vy=1,Nall,save.path, mc.cores =NCORES,maxcv,maxcv_stop,maxcv_autocheck,extra.java.arguments)	

for(i in 1:A) unlink(paste0(save.path,"/a",i,"/*"))

asnps <- lapply(beta.list, names)
usnps <- Reduce(union,asnps)	
csta <- MGfm(stfm,usnps,Nall,cred,cpp)
S <- A 
   if(!is.null(nrow(csta$summary))) {
    mafdf <- NULL
    tsnps <- rownames(csta$summary)
    for(s in 1:S) {
     rownames(gwas.list[[s]]) <- gwas.list[[s]]$rsID
     eaf <- gwas.list[[s]][tsnps,"EAF"]
     maf <- eaf*(eaf<.5) + (1-eaf)*(eaf>=0.5)
     mafdf <- cbind(mafdf,maf)
     }
     names(mafdf) <- studies
    cstaMAF <- data.frame(csta$summary, mafdf)
  } else{
    cstaMAF <- csta$summary
  }

anc <- studies

  if(!is.null(nrow(csta$summary))) {
   Ntot <- sum(Nall)
   Nw <- Nall/Ntot
   pmaf <- apply(cstaMAF[,-1],1,mean,weights=Nw,na.rm=T)
   Aprop <- apply(cstaMAF[,-1],1, function(x) mean(!is.na(x)))
   Aanc <- apply(cstaMAF[,-1],1, function(x) paste(anc[which(!is.na(x))],collapse=","))
   cstaS <- data.frame(csta$summary,maf_pooled=pmaf,prop_groups=Aprop,groups=Aanc)
  } else{
    cstaS <- csta$summary
  }

out <- list(CSsummary=cstaS, CSdetail=csta)
 
return(out)
}

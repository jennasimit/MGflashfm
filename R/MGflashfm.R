
#' @title Wrapper for Multi-group multi-trait fine-mapping, using FLASHFMwithJAM
#' @param gwas.list List of A lists objects, where A is the number of groups; gwas.list\[\[i\]\] is a list for group  i with M data.frames (one for each trait) 
#' with 3 columns named: rsID, beta, EAF; 
#' if trait names are provided for the M data.frames (same names across studies), these trait names are given in output
#' if group  names are provided for the A data.frames, these group  names are given in the output
#' @param LD.list List of A data.frame objects, where A is the number of studies;  groups must be in same order as in gwas.list; LD.list\[\[i\]\] is the SNP correlation matrix for group  i and must have snp names in row.names and col.names
#' @param covY.list List of A trait covariance matrices; trait columns should be in same order as traits are listed in gwas.list\[\[i\]\]
#' @param Nall List of components with same length as number of Studies: Nall\[\[i\]\] is the M-vector of trait sample sizes for group  i, where M is the number of traits
#' @param multi TRUE for multi-group  multi-trait fine-mapping; FALSE for multi-group  single-trait fine-mapping; default TRUE
#' @param TOdds target odds of no sharing to sharing; default is 1
#' @param maxcv starting value for maximum number of causal variants
#' @param maxcv_stop maximum value to consider for maximum number of causal variants; maxcv_stop >= maxcv.
#' @param save.path Path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1").
#' @param cpp cumulative posterior probability threshold for selecting top models; default 0.99
#' @param cred Level for multi-group  credible set, default 0.99 for 99% credible sets
#' @param NCORES number of cores for parallel computing; recommend NCORES=max(A,M), but if on Windows, use NCORES=1
#' @param jam.nM.iter in millions, number of iterations to use in JAM; defailt 1 (1 million)
#' @param flashfmRET TRUE to return single-group  flashfm output; default FALSE 
#' @param extra.java.arguments A character string to be passed through to the java command line. E.g. to specify a
#' different temporary directory by passing "-Djava.io.tmpdir=/Temp".
#' @return List consisting of two objects: CSsummary = List of one data.frame for each trait; each trait data.frame gives the variants in the multi-group  credible set for the trait, the msMPP, pooled MAF, proportion of studies 
#' that contain the variant, names of studies that contain the variant
#' CSdetail = \[\[1\]\] a list of multi-group  credible sets (variants and their msMPP), one for each trait; \[\[2\]\] details = for each trait, list of all multi-group  models and variants and their msPP, msMPP
#' if flashfmRET=TRUE, also returns a list of flashfm output for each group  
#' @export
#' @import R2BGLiMS
#' @author Jenn Asimit
MGFLASHFMwithJAM <- function(gwas.list, LD.list, covY.list, Nall, multi=TRUE, TOdds=1, maxcv=1, maxcv_stop = 20, save.path, cpp=0.99, cred=0.99,NCORES=1,jam.nM.iter=1,flashfmRET=FALSE,extra.java.arguments=NULL) {
  
  maxcv_autocheck = TRUE
  S <- length(gwas.list)
  M <- length(gwas.list[[1]])
  
  if(!is.null(names(gwas.list[[1]]) )) {
   traits <- names(gwas.list[[1]])
  } else { traits <- paste0("trait",1:M)
  		for(j in 1:S) names(gwas.list[[j]]) <- traits
  }
  
  if(!is.null(names(gwas.list) )) {
   studies <- names(gwas.list)
  } else { studies <- paste0("group",1:S)
  		names(gwas.list) <- studies
  }
  
  ybar <- rep(0,M)
  names(ybar) <- traits
 
  flashfm.out <- ssnps <- allflashfmPP <- vector("list",S)
  for(i in 1:S) {
   flashfm.out[[i]] <- FLASHFMwithJAMd(gwas.list[[i]], LD.list[[i]], ybar, Nall[[i]], save.path=save.path, TOdds = TOdds, covY=covY.list[[i]], 
                        cpp = cpp, NCORES=NCORES, maxcv=maxcv, maxcv_stop = maxcv_stop,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)
	ssnps[[i]] <- Reduce(intersect,lapply(gwas.list[[i]],function(x) x$rsID))
	allflashfmPP[[i]] <- flashfm.out[[i]]$mpp.pp
  }
  names(flashfm.out) <- studies
  
  allsnps <- Reduce(union,ssnps)
  
  csta <- MGflashfm(allflashfmPP,Nall,allsnps,cred=cred,multi=multi,cpp=cpp,NCORES=NCORES)

  cstaMAF <- vector("list",length(csta$summary))
  for(i in 1:length(csta$summary)) {
   if(!is.null(nrow(csta$summary[[i]]))) {
    mafdf <- NULL
    tsnps <- rownames(csta$summary[[i]])
    for(s in 1:S) {
     rownames(gwas.list[[s]][[i]]) <- gwas.list[[s]][[i]]$rsID
     eaf <- gwas.list[[s]][[i]][tsnps,"EAF"]
     maf <- eaf*(eaf<.5) + (1-eaf)*(eaf>=0.5)
     mafdf <- cbind(mafdf,maf)
     }
     names(mafdf) <- studies
    cstaMAF[[i]] <- data.frame(csta$summary[[i]], mafdf)
  } else{
    cstaMAF[[i]] <- csta$summary[[i]]
  }
}
names(cstaMAF) <- names(csta$summary)

cstaS <-  vector("list",length(csta$summary))
Ndf <- data.frame(Nall)
anc <- studies

for(i in 1:length(csta$summary)) {
  if(!is.null(nrow(csta$summary[[i]]))) {
   Ntot <- sum(Ndf[,i])
   Nw <- Ndf[,i]/Ntot
   pmaf <- apply(cstaMAF[[i]][,-1],1,mean,weights=Nw,na.rm=T)
   Aprop <- apply(cstaMAF[[i]][,-1],1, function(x) mean(!is.na(x)))
   Aanc <- apply(cstaMAF[[i]][,-1],1, function(x) paste(anc[which(!is.na(x))],collapse=","))
   cstaS[[i]] <- data.frame(csta$summary[[i]],maf_pooled=pmaf,prop_groups=Aprop,groups=Aanc)
  } else{
    cstaS[[i]] <- csta$summary[[i]]
  }
}
names(cstaS) <- names(csta$summary)

if(flashfmRET) {
 out <- list(CSsummary=cstaS, CSdetail=csta, flashfm.out=flashfm.out)
} else {
	out <- list(CSsummary=cstaS, CSdetail=csta)
}
  
return(out)
  
}


#' @title Multi-group  multi-trait fine-mapping, using output from FLASHFMwithJAM or FLASHFMwithFINEMAP; gives both multi- and single-trait results from flashfm object
#' @param gwas.list List of A lists objects, where A is the number of groups; gwas.list\[\[i\]\] is a list for group  i with M data.frames (one for each trait) 
#' with 3 columns named: rsID, beta, EAF; 
#' if trait names are provided for the M data.frames (same names across groups), these trait names are given in output
#' if group  names are provided for the A data.frames, these group  names are given in the output
#' @param flashfm.list List of flashfm output from each of the groups
#' @param Nall List of components with same length as number of groups: Nall\[\[i\]\] is the M-vector of trait sample sizes for group  i, where M is the number of traits
#' @param cred Level for credible set; default 0.99
#' @param multi TRUE for multi-group  multi-trait fine-mapping; FALSE for multi-group  single-trait fine-mapping; default TRUE
#' @param cpp cumulative posterior probability threshold for selecting top models; this is ignored when maxmod is spespecified
#' @param NCORES number of cores for parallel computing; recommend NCORES=M, but if on Windows, use NCORES=1
#' @return List consisting of two objects: CSsummary = List of one data.frame for each trait; each trait data.frame gives the variants in the multi-group  credible set for the trait, the msMPP, pooled MAF, proportion of groups 
#' that contain the variant, names of groups that contain the variant
#' CSdetail = \[\[1\]\] a list of multi-group  credible sets (variants and their msMPP), one for each trait; \[\[2\]\] details = for each trait, list of all multi-group  models and variants and their msPP, msMPP
#' @export
#' @author Jenn Asimit
MGflashfmRET <- function(gwas.list,flashfm.list,Nall,cred=0.99,multi=FALSE,cpp=0.99,NCORES=1) {
 
   S <- length(gwas.list)
  M <- length(gwas.list[[1]])
  
  if(!is.null(names(gwas.list[[1]]) )) {
   traits <- names(gwas.list[[1]])
  } else { traits <- paste0("trait",1:M)
  		for(j in 1:S) names(gwas.list[[j]]) <- traits
  }
  
  if(!is.null(names(gwas.list) )) {
   studies <- names(gwas.list)
  } else { studies <- paste0("group",1:S)
  		names(gwas.list) <- studies
  }

 
  allflashfmPP <- ssnps <- vector("list",S)
  for(i in 1:S) {
     allflashfmPP[[i]] <- flashfm.list[[i]]$mpp.pp
 ssnps[[i]] <- Reduce(intersect,lapply(gwas.list[[i]],function(x) x$rsID))
   }
   names(allflashfmPP) <- studies  
   allsnps <- Reduce(union,ssnps)
 
   csta <- MGflashfm(allflashfmPP,Nall,allsnps,cred=cred,multi=multi,cpp=cpp,NCORES=NCORES)
 
  cstaMAF <- vector("list",length(csta$summary))
   for(i in 1:length(csta$summary)) {
    if(!is.null(nrow(csta$summary[[i]]))) {
     mafdf <- NULL
     tsnps <- rownames(csta$summary[[i]])
     for(s in 1:S) {
      rownames(gwas.list[[s]][[i]]) <- gwas.list[[s]][[i]]$rsID
      eaf <- gwas.list[[s]][[i]][tsnps,"EAF"]
      maf <- eaf*(eaf<.5) + (1-eaf)*(eaf>=0.5)
      mafdf <- cbind(mafdf,maf)
      }
      names(mafdf) <- studies
     cstaMAF[[i]] <- data.frame(csta$summary[[i]], mafdf)
   } else{
     cstaMAF[[i]] <- csta$summary[[i]]
   }
 }
 names(cstaMAF) <- names(csta$summary)
 
 cstaS <-  vector("list",length(csta$summary))
 Ndf <- data.frame(Nall)
 anc <- studies
 
 for(i in 1:length(csta$summary)) {
   if(!is.null(nrow(csta$summary[[i]]))) {
    Ntot <- sum(Ndf[,i])
    Nw <- Ndf[,i]/Ntot
    pmaf <- apply(cstaMAF[[i]][,-1],1,mean,weights=Nw,na.rm=T)
    Aprop <- apply(cstaMAF[[i]][,-1],1, function(x) mean(!is.na(x)))
    Aanc <- apply(cstaMAF[[i]][,-1],1, function(x) paste(anc[which(!is.na(x))],collapse=","))
    cstaS[[i]] <- data.frame(csta$summary[[i]],maf_pooled=pmaf,prop_studies=Aprop,groups=Aanc)
   } else{
     cstaS[[i]] <- csta$summary[[i]]
   }
 }
 names(cstaS) <- names(csta$summary)
 
 
  out <- list(CSsummary=cstaS, CSdetail=csta)
   
 return(out)
   
 }



#' @title Multi-group  multi-trait fine-mapping, using output from FLASHFMwithJAM or FLASHFMwithFINEMAP; gives both multi- and single-trait results from flashfm object
#' @param allflashfmPP List of components with same length as number of Studies: allflashfmPP\[\[i\]\] = flashfmOUT\[\[i\]\]$mpp.pp for group  i
#' @param Nall List of components with same length as number of Studies: Nall\[\[i\]\] is the M-vector of trait sample sizes for group  i, where M is the number of traits
#' @param snps vector of all variants that exist in at least one group 
#' @param cred Level for credible set; default 0.99
#' @param multi TRUE for multi-group  multi-trait fine-mapping; FALSE for multi-group  single-trait fine-mapping; default TRUE
#' @param cpp cumulative posterior probability threshold for selecting top models; this is ignored when maxmod is spespecified
#' @param NCORES number of cores for parallel computing; recommend NCORES=M, but if on Windows, use NCORES=1
#' @return List consisting of two objects: summary = a list of multi-group  credible sets (variants and their msMPP), one for each trait; details = for each trait, list of all multi-group  models and variants and their msPP, msMPP
#' @export
#' @author Jenn Asimit
MGflashfm <- function(allflashfmPP,Nall,snps,cred=.99,multi=TRUE,cpp=0.99,NCORES=1) {

A <- length(allflashfmPP) 
M <- length(allflashfmPP[[1]]$PP)
traits <- names(allflashfmPP[[1]]$PP)

cpp.thr=cpp

fmta <-vector("list",M)
ivec <- vector("list",M)
for(i in 1:M) ivec[[i]] <- i
#for(i in 1:M) {print(i); fmta[[i]] <- fmta1(i,allflashfmPP,Nall,snps,cred,multi)}
fmta <- parallel::mclapply(ivec,fmta1,allflashfmPP,Nall,snps,cred,multi,cpp.thr,mc.cores =NCORES)	
names(fmta) <- traits
 
 if(is.null(names(allflashfmPP))) {
  studies <- paste0("group ",1:A)
  } else studies <- names(allflashfmPP)


sg <- vector("list",M)
for(i in 1:M) { 
 if(!is.null(fmta[[i]]$credset)){
  cs <- fmta[[i]]$credset
  cs <- cs[which(!is.na(cs))] # if single snp in cs, sometimes have NA variant come up b/c of way cs constructed 
  sg[[i]] <- data.frame(fmta[[i]]$mgMPP[cs])
  sg[[i]] <- sg[[i]][order(sg[[i]][,1],decreasing=T),,drop=FALSE]
  colnames(sg[[i]]) <- "mgMPP"
 } else {
  sg[[i]] <- "Insufficient evidence of shared causal variants between groups."
 }
}
names(sg) <- names(fmta)


return(list(summary=sg,details=fmta)) 
}



fmta1 <- function(i,allflashfmPP,Nall,snps,cred,multi,cpp.thr=0.99) {
  A <- length(allflashfmPP)
  flashfm.out <- vector("list",A)
  N <- c()
  for(j in 1:A) {
  	flashfm.out[[j]] <- allflashfmPP[[j]]$PP[[i]] # for trait i, ancestry results
  	N <- c(N,Nall[[j]][i])
  	}

  fmta <- flashfmTA1(flashfm.out, N, snps, cred=cred,multi,cpp.thr)  # trans-ancestry flashfm trait i
  return(fmta)
}  



flashfmTA1 <- function (flashfm.out1, N,snps,cred=0.99,multi,cpp.thr=0.99) 
{ 

	A <- length(N)   
	nsnps <- length(snps) 	
 
 STR <- PP <- vector("list",A)   
 for(i in 1:A) {
	STR[[i]] <- rownames(flashfm.out1[[i]])
	if(multi) PP[[i]] <- flashfm.out1[[i]][,2]
	if(!multi) PP[[i]] <- flashfm.out1[[i]][,1]
 	}
   
  n <- length(STR) # number of studies
    if (n < 2) 
        stop("Need at least 2 groups")
    if (length(STR) != n || length(PP) != n) 
        stop("STR and PP need to have the same lengths")
    if (is.null(names(STR))) 
        names(STR) <- paste0("popn", seq_along(STR))
    qt <- names(STR)

 bestmod.thr <- vector("list", A)
    for (i in 1:A) {
    	tmp <- PP2snpmod(data.frame(PP=PP[[i]],str=STR[[i]]))
        bm <- best.models.cpp(tmp, cpp.thr = cpp.thr)
        bestmod.thr[[i]] <- bm$models
    }
    STR <- lapply(bestmod.thr, "[[", "str")
    PP <- lapply(bestmod.thr, "[[", "PP")

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

if(A==2) { 
	osnps <- intersect(unlist(strsplit(STR[[1]],"%")),unlist(strsplit(STR[[2]],"%")))
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
    if(A==3) out <- TA3flashfm1(PPn,nsnpspermodel,SS,nsnps,snps,cred)  
    if(A==4) out <- TA4flashfm1(PPn,nsnpspermodel,SS,nsnps,snps,cred)  
    if(A==5) out <- TA5flashfm1(PPn,nsnpspermodel,SS,nsnps,snps,cred)  
    if(A==6) out <- TA6flashfm1(PPn,nsnpspermodel,SS,nsnps,snps,cred)  
    
    tampp <- sort(out$taMPP,decreasing=T)
    tapp <- out$taPP[order(out$taPP[,1],decreasing=T),,drop=F]
    
    Fout <- list(mgMPP=tampp,mgPP=tapp,credset=out$cs,cred=cred)
    } else{
     Fout <- list(mgMPP=NULL,mgPP=NULL,credset=NULL,cred=cred)
    }
    } else {
     Fout <- list(mgMPP=NULL,mgPP=NULL,credset=NULL,cred=cred)
    }
    return(Fout)
}


    
### 2 group  functions

calctauTA2 <- function(n1,n2,nsnps) {
    num <- choose(nsnps,n1)
    denom <- choose(nsnps,n1) - choose(nsnps-n2,n1) 
    out <- log(num)-log(denom)
    if(denom==0) out <- 0
    return(out)
}

joint <- function(mod1,mod2,PPn1,PPn2) {	
	x <- unlist(strsplit(mod1,"%",fixed=TRUE))
	y <- unlist(strsplit(mod2,"%",fixed=TRUE))
	check <- any(x %in% y)
 	out <- NA
 	if(check) {
 		ss <- union(x,y)
 		name.out <- paste(sort(ss),collapse="%")
 		jpp <- PPn1[mod1]+PPn2[mod2]
 		names(jpp) <- name.out
# 		out <- list(name.out,jpp)
		out <- jpp
 		}
 	return(out)
}

vjoint <- Vectorize(joint,c("mod1","mod2"),USE.NAMES=FALSE)

tauTA2fn <- function(i,j,tau.mat) tau.mat[i,j]

tauTA2vec <- Vectorize(tauTA2fn,c("i","j")) 

MPPcalc <- function(PP1) {
    mnames <- as.list(rownames(PP1))
    pp <- PP1[,1]
#    msep <- apply(matrix(1:length(mnames), ncol = 1), 1, sep.fn, mnames)
    msep <- lapply(mnames, function(x) unlist(strsplit(as.character(x), "%"))) 
    gnames <- unique(unlist(msep))
    ind.mat <- sapply(msep, function(x) gnames %in% x )
    ind.mat <- t(ind.mat)
    ind.sparse <- Matrix::Matrix(ind.mat, sparse=TRUE)
    mpp1 <- pp*ind.mat
    mpp <- Matrix::colSums(mpp1)
    names(mpp) <- gnames     
    return(mpp)
}


sep.fn <- function (k, mnames) 
{
    msep <- unlist(strsplit(as.character(mnames[k]), "%"))
    return(msep)
}

### 

###### 2 groups ######
TA2flashfm1 <- function(PPn,nsnpspermodel,nsnps,cred=0.99) {
# find joint PP for all pair-wise models between all ancestries , sort, and keep those for cpp=0.99
# go across these model pairs and keep intersecting snps, common to both bits of joint model

  
 namesPPn <- lapply(PPn,names)
 names.all <- expand.grid(namesPPn,stringsAsFactors = FALSE) # all model combinations across ancestries
 modj <- vjoint(names.all[,1],names.all[,2],PPn1=PPn[[1]],PPn2=PPn[[2]])
 ind.keep <- which(!is.na(modj))
 PPj <- modj[ind.keep]

 rm(names.all); gc(verbose = FALSE)

 maxsnps <- max(unlist(nsnpspermodel))
 tau.mat <- matrix(0,nrow=maxsnps,ncol=maxsnps)
 for(i in 1:maxsnps){
  for(j in 1:maxsnps) {
   tau.mat[i,j] <- calctauTA2(i,j,nsnps)
  }
 }

 tau.all <- expand.grid(nsnpspermodel)
 tau.keep <- tau.all[ind.keep,]
 PPtau <- tauTA2vec(tau.keep[,1],tau.keep[,2],tau.mat)
 
 PPjoint <- PPj+PPtau

  tt <- as.matrix(exp(PPjoint),ncol=1)
  ppj <- rowsum(tt, rownames(tt))
  ppj <- ppj/sum(ppj)
 taMPP <- MPPcalc(ppj) 

 tmp <- ppj[order(ppj[,1], decreasing = TRUE),]
 cpp <- cumsum(tmp)
 wh <- which(cpp <= cred)
 if (!length(wh)) wh <- 1
 if(cpp[max(wh)] < cred) wh <- c(wh, max(wh) + 1)
 keepmodPP <- tmp[wh]
 mods <- names(keepmodPP)
 cs <- unique(unlist(strsplit(mods,"%")))
 
 return(list(cs=cs,taMPP=taMPP, taPP=ppj))
 }

### 3 ancestry functions

calctauTA3 <- function(n1,n2,n3,nsnps) {
    ns <- sort(c(n1,n2,n3),decreasing=FALSE)
    num <- choose(nsnps,ns[1])*choose(nsnps,ns[2])
    den <- choose(nsnps,ns[1])*choose(nsnps,ns[2]) - choose(nsnps-ns[3],ns[2])*choose(nsnps-ns[3]-ns[2],ns[1])
    out <- log(num)-log(den)
    if(den==0) out <- 0
    return(out)
}



joint2 <- function(mod1,mod2) {	
	check <- any(mod1 %in% mod2)
 	out <- NA
 	if(check) {
 		out <- sort(union(mod1,mod2))
 		}
 	return(out)
}

vjoint2 <- Vectorize(joint2,c("mod1","mod2"),USE.NAMES=FALSE,SIMPLIFY=FALSE)

joinmod <- function(x,y) paste(sort(union(x,y)),collapse="%")

sortunion <- function(mod1,mod2) {	
 		out <- sort(union(mod1,mod2))
 	return(out)
}

vsortunion <- Vectorize(sortunion,c("mod1","mod2"),USE.NAMES=FALSE,SIMPLIFY=FALSE)



MPPcalcstr <- function(PP1,snps) 
{
    mnames <- as.character(rownames(PP1))
    pp <- PP1[,1]
    msep <- sapply(mnames,function(x) strsplit(x,"%"))
    msnps <- unique(unlist(msep))
    ind.mat <- sapply(msep, function(x) msnps %in% x )
    ind.mat <- t(ind.mat)
    mpp1 <- pp*ind.mat
    mpp <- Matrix::colSums(mpp1)
    names(mpp) <- msnps
    gnames <- as.integer(factor(snps,levels = snps)) 
    dsnps <- setdiff(gnames,msnps)
    mpp0 <- rep(0,length(dsnps))
    names(mpp0) <- dsnps
    mpp <- c(mpp,mpp0)
    names(mpp) <- snps[as.numeric(names(mpp))]     
    return(mpp)
}





int2name <- function(imod,snps) {
 	msnps <- as.integer(unlist(strsplit(imod,"%")))
 	mod <- paste0(snps[msnps],collapse="%")
 	return(mod)
 }
 
Vintname <- Vectorize(int2name,"imod",USE.NAMES=FALSE)	

###    

######## 3 groups  ###########
TA3flashfm1 <- function(PPn,nsnpspermodel,SS,nsnps,snps,cred=0.99) {
 

	nmax <- max(sapply(nsnpspermodel,max))
    tau3 <- vector("list",nmax)
    for(i in 0:nmax) tau3[[i+1]] <- matrix(0,nrow=nmax+1,ncol=nmax+1)
    for(i in 0:nmax){ 
    	for(j in 0:nmax) { 
    		for(k in 0:nmax) { 
    			tau3[[i+1]][j+1,k+1] <- calctauTA3(i,j,k,nsnps) # shift indices by 1 to allow for 0
    			} } }
				
###### PART 1 - overlap T1-T2 models with all T3 models

 namesPPn <- SS
 ## numeric version of namesPPn for speed
   mstr <- lapply(namesPPn, function(ss) {
        lapply(ss, function(x) as.integer(factor(x, levels = snps)))
    })
    names(mstr) <- NULL
 namesPPn <- mstr
 rm(mstr)
 gc()
 
 names12 <- expand.grid(namesPPn[[1]],namesPPn[[2]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) # all model combinations across ancestries 1-2
 modj12 <- vjoint2(names12[,1],names12[,2])
 ind12.keep <- which(!is.na(modj12))
 pp12c <- expand.grid(PPn[[1]],PPn[[2]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) 
 pp12a <- apply(pp12c,1,sum)
 ppj123 <- c()
 nsize12all <- expand.grid(nsnpspermodel[[1]],nsnpspermodel[[2]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)

 if(length(ind12.keep)>0) {
 	MODj12 <- modj12[ind12.keep]  # T1-T2 models with overlap 
 	ppj12 <- pp12a[ind12.keep] # T1-T2 models with overlap lPP1+lPP2
 
#tau 
 	nsize12a <- nsize12all[ind12.keep,]
 
 	tau123a <- matrix(0,nrow=nrow(nsize12a),ncol=length(nsnpspermodel[[3]]))
 	uniqt12 <- unique(nsize12a,MARGIN=1) 
	 uniqt3 <- unique(nsnpspermodel[[3]])
	 for(i in 1:nrow(uniqt12)) {
	 	mrow <- uniqt12[i,]
 		tmp <- apply(nsize12a,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		for(j in uniqt3) {
 			ind2 <- which(nsnpspermodel[[3]] == j)
 			tau123a[ind1,ind2] <- tau3[[mrow[1,1]+1]][mrow[1,2]+1,j+1] # shift indices by 1 to allow for 0
 		}
	 }
# pp123 
	 ppj12a3 <- outer(ppj12,PPn[[3]],function(x,y) x+y)
	 ppj123a <-  ppj12a3 + tau123a

# model names
	 npp3 <- namesPPn[[3]]
	 mod12a3 <- expand.grid(MODj12,npp3,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE)
 
	 mod123a <- mapply(joinmod,mod12a3[,1],mod12a3[,2])
	 
 
	 PPj123a <- as.vector(ppj123a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
	 names(PPj123a) <- mod123a
	 tt <- as.matrix(exp(PPj123a),ncol=1)
	 ppj123 <- rowsum(tt, rownames(tt)) # unique T1-T2-T3 models with T1-T2 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj123a,mod123a,tmp,mod12a3,npp3,ppj123a,ppj12a3,tau123a,nsize12a)
	 gc()
	 }
###### PART 2 - non-overlap T1-T2 models with all T3 models, check if T3 model overlaps T1 or T2 model and include those that do
	

 
 if(length(ind12.keep) < length(modj12) ) {
	
	if(length(ind12.keep)>0) {
	tmp <- names12[-ind12.keep,] # T1-T2 models with NO overlap - check if either overlaps with T3 models - check if T1 U T2 mod overlaps T3 mod
 	modj12x <- vsortunion(tmp[,1],tmp[,2])
 	ppj12x <- pp12a[-ind12.keep] 
  	} else {
  		modj12x <- vsortunion(names12[,1],names12[,2]) # T1-T2 models with NO overlap - check if either overlaps with T3 models - check if T1 U T2 mod overlaps T3 mod
 		ppj12x <- pp12a 
  	}


	
 	names123 <- expand.grid(modj12x,namesPPn[[3]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) # all model combinations across ancestries 1-2(non-overlap) -3
 
 	modj123 <- vjoint2(names123[,1],names123[,2])
 	ind123.keep <- which(!is.na(modj123))
 	
 	if(length(ind123.keep) > 0) {
 		MODj123 <- modj123[ind123.keep]  # T1-T2-T3 models with overlap 
  	 
		pp12c <- expand.grid(ppj12x,PPn[[3]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) 
 		pp12a <- apply(pp12c,1,sum)
 		ppj123b <- pp12a[ind123.keep] # T1-T2 models with non-overlap, but overlap T3 lPP1+lPP2+lPP3 
#		names.ppj123 <- sapply(MODj123,function(x) paste(x,collapse="%"))

		if(length(ind12.keep)>0) {
			nsize12a <- nsize12all[-ind12.keep,] # 1-2
			} else { nsize12a <- nsize12all }
		ns123a.ind <- expand.grid(1:nrow(nsize12a),1:length(nsnpspermodel[[3]]),KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) # indices of 1-2 -3
		ns123indkeep <- ns123a.ind[ind123.keep,] # kept indices - col1 = 1-2, col2 = 3
#	ns1 <- nsize12a[[1]][ns123indkeep[[1]]] # kept model sizes, T1
#	ns2 <- nsize12a[[2]][ns123indkeep[[1]]]
#	ns3 <- nsnpspermodel[[3]][ns123indkeep[[2]]]
	
	# over all non-overlap T1-T2 with all T3, then extract indices ns1,ns2,ns3
 		tau123a <- matrix(0,nrow=nrow(nsize12a),ncol=length(nsnpspermodel[[3]]))
 		uniqt12 <- unique(nsize12a,MARGIN=1) 
		 uniqt3 <- unique(nsnpspermodel[[3]])
	 	for(i in 1:nrow(uniqt12)) {
 			mrow <- uniqt12[i,]
 			tmp <- apply(nsize12a,1,function(x) all(x==mrow))
 			ind1 <- which(tmp)
 			for(j in uniqt3) {
 				ind2 <- which(nsnpspermodel[[3]] == j)
 				tau123a[ind1,ind2] <- tau3[[mrow[1,1]+1]][mrow[1,2]+1,j+1]
 			}
	 	} 
 	
	 	tau123 <- tauTA2vec(ns123indkeep[[1]],ns123indkeep[[2]],tau123a)
 		
# pp123 
		ppj123c <- ppj123b + tau123
	
# model names
		mod123a <- sapply(MODj123, function(x) paste(x,collapse="%"))
	

#	PPj123a <- as.vector(ppj123a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
		names(ppj123c) <- mod123a
		tt <- as.matrix(exp(ppj123c),ncol=1)
		ppj123b <- rowsum(tt, rownames(tt)) # unique T1-T2-T3 models with T1-T2 overlap lPP1+lPP2+lPP3+tau
	
	 	if(!is.null(ppj123)) {
	 		name1 <- rownames(ppj123)
	 		name2 <- rownames(ppj123b)
	 		ind <- match(name1,name2)
	 		in1 <-  which(!is.na(ind))
	 		in2 <- ind[in1]
	 		ppj2 <- ppj123[in1,] + ppj123b[in2,]
	 		ppj1a <- ppj123[which(is.na(ind)),]
	 		in2a <- setdiff((1:nrow(ppj123b)),in2)
	 		ppj2a <- ppj123b[in2a,]
	 		ppj <- c(ppj2,ppj1a,ppj2a)
	 		ppj <- data.frame(ppj)
#	 		tmp <- gtools::smartbind(t(ppj123),t(ppj123b),fill=0) # combine parts 1,2
#	 		ppj <- colSums(tmp)
 			} else {ppj <- ppj123b}
	 	rm(tmp,ppj123b,tt,ppj123c,tau123,pp12c,names123)
	 	gc()
	 	
 	
 		} else { ppj <- ppj123 }
 
} else { ppj <- ppj123 } 
 
 
 
 nameppj <- rownames(ppj)
 ppj <- ppj[,1]/sum(ppj[,1])
 names(ppj) <- nameppj
 
 ppj <- data.frame(ppj)
 taMPP <- MPPcalcstr(ppj,snps) 

 nppj <- rownames(ppj) 
 name.ppj <- Vintname(nppj,snps)
 rownames(ppj) <- name.ppj

 
 indsort <- order(ppj[,1], decreasing = TRUE)
 tmp <- ppj[indsort,]
 ntmp <- rownames(ppj)[indsort]
 cpp <- cumsum(tmp)
 wh <- which(cpp <= cred)
 if (!length(wh)) wh <- 1
 if(cpp[max(wh)] < cred) wh <- c(wh, max(wh) + 1)
 mods <- ntmp[wh]
 cs <- unique(unlist(strsplit(mods,"%")))
 
 return(list(cs=cs,taMPP=taMPP, taPP=ppj))

 }



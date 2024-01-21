
### 4 study functions

calctauTA4 <- function(n1,n2,n3,n4,nsnps) {
    ns <- sort(c(n1,n2,n3,n4),decreasing=FALSE) # ns[1] is min, ns[4] is max
    num <- choose(nsnps,ns[1])*choose(nsnps,ns[2])*choose(nsnps,ns[3])
    den <- choose(nsnps,ns[1])*choose(nsnps,ns[2])*choose(nsnps,ns[3]) - choose(nsnps-ns[4],ns[3])*choose(nsnps-ns[4]-ns[3],ns[2])*choose(nsnps-ns[4]-ns[3]-ns[2],ns[1])
    out <- log(num)-log(den)
    if(den==0) out <- 0
    return(out)
}

joint2F <- function(mod1,mod2) {	
	check <- any(mod1 %in% mod2)
 	out <- name.out <- NA
 	if(check) {
 		out <- sort(union(mod1,mod2))
 		name.out <- paste(out,collapse="%")
 		}
 	return(name.out)
}

vjoint2F <- Vectorize(joint2F,c("mod1","mod2"),USE.NAMES=FALSE,SIMPLIFY=FALSE)

joinmod <- function(mod1,mod2) {	
 		out <- sort(union(mod1,mod2))
 		name.out <- paste(out,collapse="%")
 	return(name.out)
}

vjoinmod <- Vectorize(joinmod,c("mod1","mod2"),USE.NAMES=FALSE,SIMPLIFY=FALSE)


######## 4 studies  ###########
TA4flashfm1 <- function(PPn,nsnpspermodel,SS,nsnps,snps,cred=0.99) {
 

	nmax <- max(sapply(nsnpspermodel,max))
    tau4 <- vector("list",nmax+1)
    for(i in 0:nmax) tau4[[i+1]] <- array(0,dim=rep(nmax+1,3))
    for(i in 0:nmax){ 
    	for(j in 0:nmax) { 
    		for(k in 0:nmax) { 
    			for(l in 0:nmax) {
    				tau4[[i+1]][j+1,k+1,l+1] <- calctauTA4(i,j,k,l,nsnps)
    			} } }}
				
namesPPn <- SS
 ## numeric version of namesPPn for speed
   mstr <- lapply(namesPPn, function(ss) {
        lapply(ss, function(x) as.integer(factor(x, levels = snps)))
    })
    names(mstr) <- NULL
 namesPPn <- mstr
 rm(mstr)
 gc()
 
 ###### PART 1 - find overlap/non-overlap T1-T2 models AND overlap/non-overlap T3-T4 models 
 
 names12 <- expand.grid(namesPPn[[1]],namesPPn[[2]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) # all model combinations across ancestries 1-2
 modj12 <- vjoint2(names12[,1],names12[,2])
 ind12.keep <- which(!is.na(modj12)) # overlap models; those with NA are non-overlap
 pp12c <- expand.grid(PPn[[1]],PPn[[2]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) 
 pp12a <- apply(pp12c,1,sum)
 nsize12all <- expand.grid(nsnpspermodel[[1]],nsnpspermodel[[2]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)

 names34 <- expand.grid(namesPPn[[3]],namesPPn[[4]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) # all model combinations across ancestries 3-4
 modj34 <- vjoint2(names34[,1],names34[,2])
 ind34.keep <- which(!is.na(modj34))
 pp34c <- expand.grid(PPn[[3]],PPn[[4]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) 
 pp34a <- apply(pp34c,1,sum)
 nsize34all <- expand.grid(nsnpspermodel[[3]],nsnpspermodel[[4]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)

ppj1234no <- ppj1234oo <- ppj1234nn <- ppj1234on <- c()

##### PART 2 - if have overlap T1-T2 models
if(length(ind12.keep)>0) {
 	MODj12 <- modj12[ind12.keep]  # T1-T2 models with overlap 
 	ppj12 <- pp12a[ind12.keep] # T1-T2 models with overlap lPP1+lPP2
 	nsize12a <- nsize12all[ind12.keep,]
 
	### and have overlap T3-T4 models 
	if(length(ind34.keep)>0) {
 		MODj34 <- modj34[ind34.keep]  # T3-T4 models with overlap 
 		ppj34 <- pp34a[ind34.keep] # T3-T4 models with overlap lPP1+lPP2

	## JOIN OVERLAP T1-T2 WITH OVERLAP T3-T4
#tau 
 	
 	nsize34a <- nsize34all[ind34.keep,]
 
 	tau1234a <- matrix(0,nrow=nrow(nsize12a),ncol=nrow(nsize34a))
 	uniqt12 <- unique(nsize12a,MARGIN=1) 
	uniqt34 <- unique(nsize34a,MARGIN=1)
	 for(i in 1:nrow(uniqt12)) {
	 	mrow <- uniqt12[i,]
 		tmp <- apply(nsize12a,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		for(j in 1:nrow(uniqt34)) {
 			nrow <- uniqt34[j,]
 			tmp <- apply(nsize34a,1,function(x) all(x==nrow))
 			ind2 <- which(tmp)
 			tau1234a[ind1,ind2] <- tau4[[mrow[1,1]+1]][mrow[1,2]+1,nrow[1,1]+1,nrow[1,2]+1]
 		}
	 }
# pp1234 
	 ppj12a34 <- outer(ppj12,ppj34,function(x,y) x+y)
	 ppj1234a <-  ppj12a34 + tau1234a

# model names
	 mod12a34 <- expand.grid(MODj12,MODj34,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 mod1234a <- vjoinmod(mod12a34[,1],mod12a34[,2])
	 
	 PPj1234a <- as.vector(ppj1234a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
	 names(PPj1234a) <- unlist(mod1234a)
	 tt <- as.matrix(exp(PPj1234a),ncol=1)
	 ppj1234oo <- rowsum(tt, rownames(tt)) # unique T1-T2-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj1234a,mod1234a,tmp,mod12a34,ppj1234a,ppj12a34,tau1234a)
	 gc()
	
	} # length(ind34.keep)>0 
	 
	## JOIN OVERLAP T1-T2 WITH NON-OVERLAP T3-T4 (if exists)
	ppj1234on <- c()
	if(length(ind34.keep)<length(modj34)) {
	
	if(length(ind34.keep) >0){
 		modj34a <- vsortunion(names34[,1],names34[,2])
		modj34x <- modj34a[-ind34.keep]  # T3-T4 models with non-overlap 
 		ppj34x <- pp34a[-ind34.keep] # T3-T4 models with non-overlap lPP1+lPP2
		nsize34ax <- nsize34all[-ind34.keep,]
	} else {
		modj34x <- vsortunion(names34[,1],names34[,2])
 		ppj34x <- pp34a
		nsize34ax <- nsize34all
	}
#tau 
 	
 	
 
 	tau1234a <- matrix(0,nrow=nrow(nsize12a),ncol=nrow(nsize34ax))
 	uniqt12 <- unique(nsize12a,MARGIN=1) 
	uniqt34 <- unique(nsize34ax,MARGIN=1)
	 for(i in 1:nrow(uniqt12)) {
	 	mrow <- uniqt12[i,]
 		tmp <- apply(nsize12a,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		for(j in 1:nrow(uniqt34)) {
 			nrow <- uniqt34[j,]
 			tmp <- apply(nsize34ax,1,function(x) all(x==nrow))
 			ind2 <- which(tmp)
 			tau1234a[ind1,ind2] <- tau4[[mrow[1,1]+1]][mrow[1,2]+1,nrow[1,1]+1,nrow[1,2]+1]
 		}
	 }
# pp1234 
	 ppj12a34 <- outer(ppj12,ppj34x,function(x,y) x+y)
	 ppj1234a <-  ppj12a34 + tau1234a

# model names
	 mod12a34 <- expand.grid(MODj12,modj34x,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 mod1234a <- vjoinmod(mod12a34[,1],mod12a34[,2])
	 
	 PPj1234a <- as.vector(ppj1234a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
	 names(PPj1234a) <- unlist(mod1234a)
	 tt <- as.matrix(exp(PPj1234a),ncol=1)
	 ppj1234on <- rowsum(tt, rownames(tt)) # unique T1-T2-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj1234a,mod1234a,tmp,mod12a34,ppj1234a,ppj12a34,tau1234a)
	 gc()
		} ## if(length(ind34.keep)<length(modj34)) 

	 }  # length(ind12.keep)>0
	 	
	
	 
## JOIN NON-OVERLAP T1-T2 (IF EXISTS) WITH OVERLAP T3-T4 (if exists)
#tau 
 	if(length(ind12.keep)<length(modj12) & length(ind34.keep)>0 ) {
 	
 	if(length(ind12.keep) >0){
		modj12a <- vsortunion(names12[,1],names12[,2])
		modj12x <- modj12a[-ind12.keep]  # T1-T2 models with non-overlap 
 		ppj12x <- pp12a[-ind12.keep] # T1-T2 models with non-overlap lPP1+lPP2
 		nsize12ax <- nsize12all[-ind12.keep,]
 	} else {
 		modj12x <- vsortunion(names12[,1],names12[,2])
 		ppj12x <- pp12a
 		nsize12ax <- nsize12all
 	}
 	
 	nsize34a <- nsize34all[ind34.keep,] 
 	MODj34 <- modj34[ind34.keep]  # T3-T4 models with overlap 
 	ppj34 <- pp34a[ind34.keep] # T3-T4 models with overlap lPP1+lPP2
 	 
 	tau1234a <- matrix(0,nrow=nrow(nsize12ax),ncol=nrow(nsize34a))
 	uniqt12 <- unique(nsize12ax,MARGIN=1) 
	uniqt34 <- unique(nsize34a,MARGIN=1)
	 for(i in 1:nrow(uniqt12)) {
	 	mrow <- uniqt12[i,]
 		tmp <- apply(nsize12ax,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		for(j in 1:nrow(uniqt34)) {
 			nrow <- uniqt34[j,]
 			tmp <- apply(nsize34a,1,function(x) all(x==nrow))
 			ind2 <- which(tmp)
 			tau1234a[ind1,ind2] <- tau4[[mrow[1,1]+1]][mrow[1,2]+1,nrow[1,1]+1,nrow[1,2]+1]
 		}
	 }
# pp1234 
	 ppj12a34 <- outer(ppj12x,ppj34,function(x,y) x+y)
	 ppj1234a <-  ppj12a34 + tau1234a

# model names
	 mod12a34 <- expand.grid(modj12x,MODj34,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 mod1234a <- vjoinmod(mod12a34[,1],mod12a34[,2])
	 
	 PPj1234a <- as.vector(ppj1234a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
	 names(PPj1234a) <- unlist(mod1234a)
	 tt <- as.matrix(exp(PPj1234a),ncol=1)
	 ppj1234no <- rowsum(tt, rownames(tt)) # unique T1-T2-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj1234a,mod1234a,tmp,mod12a34,ppj1234a,ppj12a34,tau1234a)
	 gc()
		} ### if(length(ind12.keep)<length(MODj12)) & length(ind34.keep)>0


## CHECK IF ANY OVERLAP BETWEEN NON-OVERLAP T1-T2 AND NON-OVERLAP T3-T4
if(length(ind12.keep)<length(modj12) & length(ind34.keep)<length(modj34) ) {	 
	
	if(length(ind12.keep) >0){
		modj12a <- vsortunion(names12[,1],names12[,2])
		modj12x <- modj12a[-ind12.keep]  # T1-T2 models with non-overlap 
 		ppj12x <- pp12a[-ind12.keep] # T1-T2 models with non-overlap lPP1+lPP2
 		nsize12ax <- nsize12all[-ind12.keep,]
 	} else {
 		modj12x <- vsortunion(names12[,1],names12[,2])
 		ppj12x <- pp12a
 		nsize12ax <- nsize12all
 	}
 	
 	if(length(ind34.keep) >0){
 		modj34a <- vsortunion(names34[,1],names34[,2])
		modj34x <- modj34a[-ind34.keep]  # T3-T4 models with non-overlap 
 		ppj34x <- pp34a[-ind34.keep] # T3-T4 models with non-overlap lPP1+lPP2
		nsize34ax <- nsize34all[-ind34.keep,]
	} else {
		modj34x <- vsortunion(names34[,1],names34[,2])
 		ppj34x <- pp34a
		nsize34ax <- nsize34all
	}

	
# model names
	 mod12a34 <- expand.grid(modj12x,modj34x,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 modj1234 <- vjoint2F(mod12a34[,1],mod12a34[,2])
	 ind1234keep <- which(!is.na(modj1234))
	 
	 if(length(ind1234keep) > 0) {
		 tau1234a <- matrix(0,nrow=nrow(nsize12ax),ncol=nrow(nsize34ax))
 		uniqt12 <- unique(nsize12ax,MARGIN=1) 
		uniqt34 <- unique(nsize34ax,MARGIN=1)
		 for(i in 1:nrow(uniqt12)) {
		 	mrow <- uniqt12[i,]
 			tmp <- apply(nsize12ax,1,function(x) all(x==mrow))
 			ind1 <- which(tmp)
 			for(j in 1:nrow(uniqt34)) {
 				nrow <- uniqt34[i,]
 				tmp <- apply(nsize34ax,1,function(x) all(x==nrow))
 				ind2 <- which(tmp)
 				tau1234a[ind1,ind2] <- tau4[[mrow[1,1]+1]][mrow[1,2]+1,nrow[1,1]+1,nrow[1,2]+1]
 			}
		 }
# pp1234 
	 ppj12a34 <- outer(ppj12x,ppj34x,function(x,y) x+y)
	 ppj1234a <-  ppj12a34 + tau1234a
	 
	 PPj1234a <- as.vector(ppj1234a[ind1234keep]) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
	 names(PPj1234a) <- modj1234[ind1234keep]
	 tt <- as.matrix(exp(PPj1234a),ncol=1)
	 ppj1234nn <- rowsum(tt, rownames(tt)) # unique T1-T2-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj1234a,tmp,mod12a34,ppj1234a,ppj12a34,tau1234a)
	 gc()
	 }

} # length(ind12.keep)<length(MODj12) & length(ind34.keep)<length(MODj34)


#ppj1234no <- ppj1234oo <- ppj1234nn <- ppj1234on <- c()

all.list <- vector("list",4)
all.list[[1]] <- ppj1234oo
all.list[[2]] <- ppj1234on
all.list[[3]] <- ppj1234no
all.list[[4]] <- ppj1234nn

have <- which(sapply(all.list,function(x) !is.null(x)))


	if(length(have)==1)  ppj <- all.list[[have]]
	if(length(have)>1) {
			ppj123 <- all.list[[have[1]]]
			ppj123b <- all.list[[have[2]]]
	 		name1 <- rownames(ppj123)
	 		name2 <- rownames(ppj123b)
	 		ppj2 <- ppj1a <- ppj2a <- c()
	 		ind <- match(name1,name2)
	 		in1 <-  which(!is.na(ind))
	 		in2a <- c()
	 				if(length(in1)>0) {
	 					in2a <- ind[in1]
	 					ppj2 <- ppj123[in1,] + ppj123b[in2a,]
	 					}
	 				in2 <- which(is.na(ind))
	 				if(length(in2) > 0) {	
	 					ppj1a <- ppj123[in2,]
	 					in2b <- setdiff((1:nrow(ppj123b)),in2a)
	 					ppj2a <- ppj123b[in2b,]
	 					}
	 		ppj <- matrix(c(ppj2,ppj1a,ppj2a),ncol=1,dimnames=list(names(c(ppj2,ppj1a,ppj2a)),NULL))
	 		if(length(have) > 2) {
	 			ppj123 <- ppj
				ppj123b <- all.list[[have[3]]]
	 			name1 <- rownames(ppj123)
		 		name2 <- rownames(ppj123b)
		 		ppj2 <- ppj1a <- ppj2a <- c()
	 			ind <- match(name1,name2)
	 			in1 <-  which(!is.na(ind))
	 			in2a <- c()
	 				if(length(in1)>0) {
	 					in2a <- ind[in1]
	 					ppj2 <- ppj123[in1,] + ppj123b[in2a,]
	 					}
	 				in2 <- which(is.na(ind))
	 				if(length(in2) > 0) {	
	 					ppj1a <- ppj123[in2,]
	 					in2b <- setdiff((1:nrow(ppj123b)),in2a)
	 					ppj2a <- ppj123b[in2b,]
	 					}
	 			ppj <- matrix(c(ppj2,ppj1a,ppj2a),ncol=1,dimnames=list(names(c(ppj2,ppj1a,ppj2a)),NULL))
	 			if(length(have) == 4) {
	 				ppj123 <- ppj
					ppj123b <- all.list[[have[4]]]
	 				name1 <- rownames(ppj123)
		 			name2 <- rownames(ppj123b)
		 			ppj2 <- ppj1a <- ppj2a <- c()
	 				ind <- match(name1,name2)
	 				in1 <-  which(!is.na(ind))
	 				in2a <- c()
	 				if(length(in1)>0) {
	 					in2a <- ind[in1]
	 					ppj2 <- ppj123[in1,] + ppj123b[in2a,]
	 					}
	 				in2 <- which(is.na(ind))
	 				if(length(in2) > 0) {	
	 					ppj1a <- ppj123[in2,]
	 					in2b <- setdiff((1:nrow(ppj123b)),in2a)
	 					ppj2a <- ppj123b[in2b,]
	 					}
	 				ppj <- matrix(c(ppj2,ppj1a,ppj2a),ncol=1,dimnames=list(names(c(ppj2,ppj1a,ppj2a)),NULL))
 					} 
 				}
 			} 				
 		
 		
 		nameppj <- rownames(ppj)
 		ppj <- data.frame(ppj,row.names=nameppj)
 		
 		
	 	gc()
	 	
 
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
 
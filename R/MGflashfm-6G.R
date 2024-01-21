
### 6 study functions

calctauTA6 <- function(n1,n2,n3,n4,n5,n6,nsnps) {
    ns <- sort(c(n1,n2,n3,n4,n5,n6),decreasing=FALSE) # ns[1] is min, ns[6] is max
    num <- choose(nsnps,ns[1])*choose(nsnps,ns[2])*choose(nsnps,ns[3])*choose(nsnps,ns[4])*choose(nsnps,ns[5])
    den <- choose(nsnps,ns[1])*choose(nsnps,ns[2])*choose(nsnps,ns[3])*choose(nsnps,ns[4])*choose(nsnps,ns[5]) - 
    choose(nsnps-ns[6],ns[5])*choose(nsnps-ns[6]-ns[5],ns[4])*choose(nsnps-ns[6]-ns[5]-ns[4],ns[3])*choose(nsnps-ns[6]-ns[5]-ns[4]-ns[3],ns[2])*
    choose(nsnps-ns[6]-ns[5]-ns[4]-ns[3]-ns[2],ns[1])
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


######## 6 studies  ###########
TA6flashfm1 <- function(PPn,nsnpspermodel,SS,nsnps,snps,cred=0.99) {

	nmax <- max(sapply(nsnpspermodel,max))
    tau6 <- vector("list",nmax+1)
    for(i in 0:nmax) tau6[[i+1]] <- array(0,dim=rep(nmax+1,5))
    for(i in 0:nmax){ 
    	for(j in 0:nmax) { 
    		for(k in 0:nmax) { 
    			for(l in 0:nmax) {
    				for(m in 0:nmax) {
    					for(n in 0:nmax) {
    				tau6[[i+1]][j+1,k+1,l+1,m+1,n+1] <- calctauTA6(i,j,k,l,m,n,nsnps)
    			} } }}}}
				
namesPPn <- SS
 ## numeric version of namesPPn for speed
   mstr <- lapply(namesPPn, function(ss) {
        lapply(ss, function(x) as.integer(factor(x, levels = snps)))
    })
    names(mstr) <- NULL
 namesPPn <- mstr
 rm(mstr)
 gc()


###### PART 1 - find overlap/non-overlap T1-T2 models AND overlap/non-overlap T3-T4 models AND overlap/non-overlap T5-T6 models
 
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

 names56 <- expand.grid(namesPPn[[5]],namesPPn[[6]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) # all model combinations across ancestries 3-4
 modj56 <- vjoint2(names56[,1],names56[,2])
 ind56.keep <- which(!is.na(modj56))
 pp56c <- expand.grid(PPn[[5]],PPn[[6]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) 
 pp56a <- apply(pp56c,1,sum)
 nsize56all <- expand.grid(nsnpspermodel[[5]],nsnpspermodel[[6]],KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)



#ppj123456oaa <- ppj123456noa <- ppj123456nno <- ppj123456nnn1 <- ppj123456nnn2  <- NULL
ppj123456ooo <- ppj123456oon <- ppj123456ono <- ppj123456onn <- ppj123456noo <- ppj123456non <- ppj123456nno <- ppj123456nnn1 <- ppj123456nnn2  <- NULL


allMODj34 <- vsortunion(names34[,1],names34[,2])
allMODj56 <- vsortunion(names56[,1],names56[,2])

### NON-OVERLAP MODELS T1-T2
	if(length(ind12.keep)>0) {
 		modj12a <- vsortunion(names12[,1],names12[,2])
		modj12x <- modj12a[-ind12.keep]  # T1-T2 models with non-overlap 
 		ppj12x <- pp12a[-ind12.keep] # T1-T2 models with non-overlap lPP1+lPP2
 		nsize12ax <- nsize12all[-ind12.keep,]
 	} else {
 		modj12x <- vsortunion(names12[,1],names12[,2])
 		ppj12x <- pp12a # T1-T2 models with non-overlap lPP1+lPP2
 		nsize12ax <- nsize12all
 	}

### NON-OVERLAP MODELS T3-T4 	
 	if(length(ind34.keep)>0) {
 		modj34a <- vsortunion(names34[,1],names34[,2])
		modj34x <- modj34a[-ind34.keep]  # models with non-overlap 
 		ppj34x <- pp34a[-ind34.keep] #  models with non-overlap lPP1+lPP2
 		nsize34ax <- nsize34all[-ind34.keep,]
 	} else {
 		modj34x <- vsortunion(names34[,1],names34[,2])
 		ppj34x <- pp34a # models with non-overlap lPP1+lPP2
 		nsize34ax <- nsize34all
 	}

### NON-OVERLAP MODELS T5-T6
 	if(length(ind56.keep)>0) {
 		modj56a <- vsortunion(names56[,1],names56[,2])
		modj56x <- modj56a[-ind56.keep]  # models with non-overlap 
 		ppj56x <- pp56a[-ind56.keep] #  models with non-overlap lPP1+lPP2
 		nsize56ax <- nsize56all[-ind56.keep,]
 	} else {
 		modj56x <- vsortunion(names56[,1],names56[,2])
 		ppj56x <- pp56a # models with non-overlap lPP1+lPP2
 		nsize56ax <- nsize56all
 	}

#print("overlap, non-overlap")
#print(c(length(ind12.keep),length(modj12x)))
#print(c(length(ind34.keep),length(modj34x)))
#print(c(length(ind56.keep),length(modj56x)))

##### PART 2 - if have overlap T1-T2 models - ppj123456ooo ppj123456oon ppj123456ono ppj123456onn
if(length(ind12.keep)>0) {
 	MODj12 <- modj12[ind12.keep]  # T1-T2 models with overlap 
 	ppj12 <- pp12a[ind12.keep] # T1-T2 models with overlap lPP1+lPP2
 	nsize12a <- nsize12all[ind12.keep,]
 
 	## COMBINE OVERLAP T1-T2  WITH OVERLAP T3-T4 (if exists)
 	if(length(ind34.keep)>0) {
 		MODj34 <- modj34[ind34.keep]  # T3-T4 models with overlap 
 		ppj34 <- pp34a[ind34.keep] # T3-T4 models with overlap lPP3+lPP4
 		nsize34a <- nsize34all[ind34.keep,]
 		ppj12a34 <- as.vector(outer(ppj12,ppj34,function(x,y) x+y))
 		indsize12 <- 1:nrow(nsize12a)
 		indsize34 <- 1:nrow(nsize34a)
 		indsize12a34 <- expand.grid(indsize12,indsize34,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
		nsize12a34 <- cbind(nsize12a[indsize12a34[,1],], nsize34a[indsize12a34[,2],])
		mod12a34 <- expand.grid(MODj12,MODj34,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 	
		mod1234a <- vsortunion(mod12a34[,1],mod12a34[,2])
 	 	
	 	rm(indsize12a34,indsize12,indsize34,mod12a34)
 		gc()
 	
 	## COMBINE OVERLAP T1-T2 + OVERLAP T3-T4 WITH OVERLAP T5-T6 (if exists)
 		if(length(ind56.keep)>0) { # ppj123456ooo
 			MODj56 <- modj56[ind56.keep]  # T5-T6 models with overlap 
 			ppj56 <- pp56a[ind56.keep] # T5-T6 models with overlap 
 			nsize56a <- nsize56all[ind56.keep,]

	 	ppnames <- ppj123456a <- matrix(0,nrow=nrow(nsize12a34),ncol=nrow(nsize56a))
	 	uniqt1234 <- unique(nsize12a34,MARGIN=1) 
		uniqt56 <- unique(nsize56a,MARGIN=1)
		 for(i in 1:nrow(uniqt1234)) {
		 	mrow <- uniqt1234[i,]
	 		tmp <- apply(nsize12a34,1,function(x) all(x==mrow))
	 		ind1 <- which(tmp)
	 		for(j in 1:nrow(uniqt56)) {
	 			nrow <- uniqt56[j,]
	 			tmp <- apply(nsize56a,1,function(x) all(x==nrow))
 				ind2 <- which(tmp)
 				tauij <- tau6[[mrow[1,1]+1]][mrow[1,2]+1,mrow[1,3]+1,mrow[1,4]+1,nrow[1,1]+1,nrow[1,2]+1] # shift indices by 1 to allow for 0
 				ppj123456a[ind1,ind2] <- outer(ppj12a34[ind1],ppj56[ind2],function(x,y) x+y) + tauij
 				ppn0 <- expand.grid(mod1234a[ind1],MODj56[ind2],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
 				ppn <- vjoinmod(ppn0[,1],ppn0[,2])
 				ppn <- unlist(ppn)
 				ppnames[ind1,ind2] <- matrix(ppn,nrow=length(ind1),ncol=length(ind2))
	 		}
		 }



# model names
	 
	 PPj123456a <- as.vector(ppj123456a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
	names(PPj123456a) <- as.vector(ppnames)
	 tt <- as.matrix(exp(PPj123456a),ncol=1)
	 ppj123456ooo <- rowsum(tt, rownames(tt)) # unique T1-T2-T5-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj123456a,ppj123456a)
	 gc()
	} # length(ind56.keep)>0
	
	if(length(modj56x)>0) { # ppj123456oon
		ppnames <-  ppj123456a <- matrix(0,nrow=nrow(nsize12a34),ncol=nrow(nsize56ax))
	 	uniqt1234 <- unique(nsize12a34,MARGIN=1) 
		uniqt56 <- unique(nsize56ax,MARGIN=1)
		 for(i in 1:nrow(uniqt1234)) {
		 	mrow <- uniqt1234[i,]
	 		tmp <- apply(nsize12a34,1,function(x) all(x==mrow))
	 		ind1 <- which(tmp)
	 		for(j in 1:nrow(uniqt56)) {
	 			nrow <- uniqt56[j,]
	 			tmp <- apply(nsize56ax,1,function(x) all(x==nrow))
 				ind2 <- which(tmp)
 				tauij <- tau6[[mrow[1,1]+1]][mrow[1,2]+1,mrow[1,3]+1,mrow[1,4]+1,nrow[1,1]+1,nrow[1,2]+1]
 				ppj123456a[ind1,ind2] <- outer(ppj12a34[ind1],ppj56x[ind2],function(x,y) x+y) + tauij
 				ppn0 <- expand.grid(mod1234a[ind1],modj56x[ind2],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
 				ppn <- vjoinmod(ppn0[,1],ppn0[,2])
 				ppn <- unlist(ppn)
 				ppnames[ind1,ind2] <- matrix(ppn,nrow=length(ind1),ncol=length(ind2))
#	 			tau123456a[ind1,ind2] <- tau6[[mrow[1,1]]][mrow[1,2],mrow[1,3],mrow[1,4],nrow[1,1],nrow[1,2]]
	 		}
		 }


# model names
	 
	 PPj123456a <- as.vector(ppj123456a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
	 names(PPj123456a) <- as.vector(ppnames)
	 tt <- as.matrix(exp(PPj123456a),ncol=1)
	 ppj123456oon <- rowsum(tt, rownames(tt)) # unique T1-T2-T5-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj123456a)
	 gc()
	
		} # length(modj56x)>0 (nested in length(ind34.keep)>0)
	
	} # length(ind34.keep)>0

	if(length(modj34x)>0) { # ono, onn
		
 		ppj12a34 <- as.vector(outer(ppj12,ppj34x,function(x,y) x+y))
 		indsize12 <- 1:nrow(nsize12a)
 		indsize34 <- 1:nrow(nsize34ax)
 		indsize12a34 <- expand.grid(indsize12,indsize34,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
		nsize12a34 <- cbind(nsize12a[indsize12a34[,1],], nsize34ax[indsize12a34[,2],])
		mod12a34 <- expand.grid(MODj12,modj34x,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 	
		mod1234a <- vsortunion(mod12a34[,1],mod12a34[,2])
 	 	
	 	rm(indsize12a34,indsize12,indsize34,mod12a34)
 		gc()
		
		## COMBINE OVERLAP T1-T2 + NON-OVERLAP T3-T4 WITH OVERLAP T5-T6 (if exists)
 		if(length(ind56.keep)>0) { # ppj123456ooo
 			MODj56 <- modj56[ind56.keep]  # T5-T6 models with overlap 
 			ppj56 <- pp56a[ind56.keep] # T5-T6 models with overlap 
 			nsize56a <- nsize56all[ind56.keep,]

	 	ppnames <-  ppj123456a <- matrix(0,nrow=nrow(nsize12a34),ncol=nrow(nsize56a))
	 	uniqt1234 <- unique(nsize12a34,MARGIN=1) 
		uniqt56 <- unique(nsize56a,MARGIN=1)
		 for(i in 1:nrow(uniqt1234)) {
		 	mrow <- uniqt1234[i,]
	 		tmp <- apply(nsize12a34,1,function(x) all(x==mrow))
	 		ind1 <- which(tmp)
	 		for(j in 1:nrow(uniqt56)) {
	 			nrow <- uniqt56[j,]
	 			tmp <- apply(nsize56a,1,function(x) all(x==nrow))
 				ind2 <- which(tmp)
 				tauij <- tau6[[mrow[1,1]+1]][mrow[1,2]+1,mrow[1,3]+1,mrow[1,4]+1,nrow[1,1]+1,nrow[1,2]+1]
 				ppj123456a[ind1,ind2] <- outer(ppj12a34[ind1],ppj56[ind2],function(x,y) x+y) + tauij
 				ppn0 <- expand.grid(mod1234a[ind1],modj56[ind2],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
 				ppn <- vjoinmod(ppn0[,1],ppn0[,2])
 				ppn <- unlist(ppn)
 				ppnames[ind1,ind2] <- matrix(ppn,nrow=length(ind1),ncol=length(ind2)) 				
#	 			tau123456a[ind1,ind2] <- tau6[[mrow[1,1]]][mrow[1,2],mrow[1,3],mrow[1,4],nrow[1,1],nrow[1,2]]
	 		}
		 }

# pp123456 
#	 ppj1234a56 <- outer(ppj12a34,ppj56,function(x,y) x+y)
#	 ppj123456a <-  ppj1234a56 + tau123456a

# model names
#	 mod1234a56 <- expand.grid(mod1234a,MODj56,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 #print(head(mod1234a56))
#	 mod123456a <- vjoinmod(mod1234a56[,1],mod1234a56[,2])
	 
	 PPj123456a <- as.vector(ppj123456a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
#	 names(PPj123456a) <- unlist(mod123456a)
	 names(PPj123456a) <- as.vector(ppnames)
	 tt <- as.matrix(exp(PPj123456a),ncol=1)
	 ppj123456ono <- rowsum(tt, rownames(tt)) # unique T1-T2-T5-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj123456a)
	 gc()
	} # length(ind56.keep)>0

	if(length(modj56x)>0) { # ppj123456onn
		ppnames <-  ppj123456a <- matrix(0,nrow=nrow(nsize12a34),ncol=nrow(nsize56ax))
	 	uniqt1234 <- unique(nsize12a34,MARGIN=1) 
		uniqt56 <- unique(nsize56ax,MARGIN=1)
		 for(i in 1:nrow(uniqt1234)) {
		 	mrow <- uniqt1234[i,]
	 		tmp <- apply(nsize12a34,1,function(x) all(x==mrow))
	 		ind1 <- which(tmp)
	 		for(j in 1:nrow(uniqt56)) {
	 			nrow <- uniqt56[j,]
	 			tmp <- apply(nsize56ax,1,function(x) all(x==nrow))
 				ind2 <- which(tmp)
 				tauij <- tau6[[mrow[1,1]+1]][mrow[1,2]+1,mrow[1,3]+1,mrow[1,4]+1,nrow[1,1]+1,nrow[1,2]+1]
 				ppj123456a[ind1,ind2] <- outer(ppj12a34[ind1],ppj56x[ind2],function(x,y) x+y) + tauij
 				ppn0 <- expand.grid(mod1234a[ind1],modj56x[ind2],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
 				ppn <- vjoinmod(ppn0[,1],ppn0[,2])
 				ppn <- unlist(ppn)
 				ppnames[ind1,ind2] <- matrix(ppn,nrow=length(ind1),ncol=length(ind2)) 				

#	 			tau123456a[ind1,ind2] <- tau6[[mrow[1,1]]][mrow[1,2],mrow[1,3],mrow[1,4],nrow[1,1],nrow[1,2]]
	 		}
		 }

# pp123456 
#	 ppj1234a56 <- outer(ppj12a34,ppj56x,function(x,y) x+y)
#	 ppj123456a <-  ppj1234a56 + tau123456a

# model names
#	 mod1234a56 <- expand.grid(mod1234a,modj56x,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 #print(head(mod1234a56))
#	 mod123456a <- vjoinmod(mod1234a56[,1],mod1234a56[,2])
	 
	 PPj123456a <- as.vector(ppj123456a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
#	 names(PPj123456a) <- unlist(mod123456a)
	 names(PPj123456a) <- as.vector(ppnames)
	 tt <- as.matrix(exp(PPj123456a),ncol=1)
	 ppj123456onn <- rowsum(tt, rownames(tt)) # unique T1-T2-T5-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj123456a)
	 gc()
	
		} # length(modj56x)>0 (nested in length(ind34.keep)>0)
	
	
	
	} # length(modj34x)>0



	} # length(ind12.keep)>0
	

##### if have overlap T3-T4 models - combine with non-overlap T1-T2 and with all T5-T6 - ppj123456noo, non

if(length(ind34.keep)>0 & length(modj12x) > 0) {
#print("YES: length(ind34.keep)>0 & length(modj12x) > 0")
 	MODj34 <- modj34[ind34.keep]  # T1-T2 models with overlap 
 	ppj34 <- pp34a[ind34.keep] # T1-T2 models with overlap lPP1+lPP2
 	nsize34a <- nsize34all[ind34.keep,]
 	
 	## COMBINE NON-OVERLAP T1-T2  WITH OVERLAP T3-T4 
 	ppj12a34 <- as.vector(outer(ppj12x,ppj34,function(x,y) x+y))
 	indsize12 <- 1:nrow(nsize12ax)
 	indsize34 <- 1:nrow(nsize34a)
 	indsize12a34 <- expand.grid(indsize12,indsize34,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
	nsize12a34 <- cbind(nsize12ax[indsize12a34[,1],], nsize34a[indsize12a34[,2],])
	mod12a34 <- expand.grid(modj12x,MODj34,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	mod1234a <- vsortunion(mod12a34[,1],mod12a34[,2])
 	
 	## COMBINE NON-OVERLAP T1-T2 + OVERLAP T3-T4 WITH OVERLAP T5-T6
 	if(length(ind56.keep)>0) {
 			MODj56 <- modj56[ind56.keep]  # T5-T6 models with overlap 
 			ppj56 <- pp56a[ind56.keep] # T5-T6 models with overlap 
 			nsize56a <- nsize56all[ind56.keep,] 	
 	
 	ppnames <-  ppj123456a <- matrix(0,nrow=nrow(nsize12a34),ncol=nrow(nsize56a))
 	uniqt1234 <- unique(nsize12a34,MARGIN=1) 
	uniqt56 <- unique(nsize56a,MARGIN=1)
	 for(i in 1:nrow(uniqt1234)) {
	 	mrow <- uniqt1234[i,]
 		tmp <- apply(nsize12a34,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		for(j in 1:nrow(uniqt56)) {
 			nrow <- uniqt56[j,]
 			tmp <- apply(nsize56a,1,function(x) all(x==nrow))
 			ind2 <- which(tmp)
 			tauij <- tau6[[mrow[1,1]+1]][mrow[1,2]+1,mrow[1,3]+1,mrow[1,4]+1,nrow[1,1]+1,nrow[1,2]+1]
 			ppj123456a[ind1,ind2] <- outer(ppj12a34[ind1],ppj56[ind2],function(x,y) x+y) + tauij
 			ppn0 <- expand.grid(mod1234a[ind1],MODj56[ind2],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
 			ppn <- vjoinmod(ppn0[,1],ppn0[,2])
 			ppn <- unlist(ppn)
 			ppnames[ind1,ind2] <- matrix(ppn,nrow=length(ind1),ncol=length(ind2)) 				 			
# 			tau123456a[ind1,ind2] <- tau6[[mrow[1,1]]][mrow[1,2],mrow[1,3],mrow[1,4],nrow[1,1],nrow[1,2]]
 		}
	 }

# pp123456 
#	 ppj1234a56 <- outer(ppj12a34,ppj56,function(x,y) x+y)
#	 ppj123456a <-  ppj1234a56 + tau123456a

# model names
#	 mod1234a56 <- expand.grid(mod1234a,MODj56,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 #print(head( mod1234a56))
#	 mod123456a <- vjoinmod(mod1234a56[,1],mod1234a56[,2])
	 
	 PPj123456a <- as.vector(ppj123456a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
	 names(PPj123456a) <- as.vector(ppnames)
#	 names(PPj123456a) <- unlist(mod123456a)
	 tt <- as.matrix(exp(PPj123456a),ncol=1)
	 ppj123456noo <- rowsum(tt, rownames(tt)) # unique T1-T2-T5-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj123456a)
	 gc()
	 
	 } # length(ind56.keep)>0
	 
	 if(length(modj56x)>0) { 		
 	
 	ppnames <-  ppj123456a <- matrix(0,nrow=nrow(nsize12a34),ncol=nrow(nsize56ax))
 	uniqt1234 <- unique(nsize12a34,MARGIN=1) 
	uniqt56 <- unique(nsize56ax,MARGIN=1)
	 for(i in 1:nrow(uniqt1234)) {
	 	mrow <- uniqt1234[i,]
 		tmp <- apply(nsize12a34,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		for(j in 1:nrow(uniqt56)) {
 			nrow <- uniqt56[j,]
 			tmp <- apply(nsize56ax,1,function(x) all(x==nrow))
 			ind2 <- which(tmp)
 			tauij <- tau6[[mrow[1,1]+1]][mrow[1,2]+1,mrow[1,3]+1,mrow[1,4]+1,nrow[1,1]+1,nrow[1,2]+1]
 			ppj123456a[ind1,ind2] <- outer(ppj12a34[ind1],ppj56x[ind2],function(x,y) x+y) + tauij
 			ppn0 <- expand.grid(mod1234a[ind1],modj56x[ind2],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
 			ppn <- vjoinmod(ppn0[,1],ppn0[,2])
 			ppn <- unlist(ppn)
 			ppnames[ind1,ind2] <- matrix(ppn,nrow=length(ind1),ncol=length(ind2)) 				 			
# 			tau123456a[ind1,ind2] <- tau6[[mrow[1,1]]][mrow[1,2],mrow[1,3],mrow[1,4],nrow[1,1],nrow[1,2]]
 		}
	 }

# pp123456 
#	 ppj1234a56 <- outer(ppj12a34,ppj56x,function(x,y) x+y)
#	 ppj123456a <-  ppj1234a56 + tau123456a

# model names
#	 mod1234a56 <- expand.grid(mod1234a,modj56x,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 #print(head( mod1234a56))
#	 mod123456a <- vjoinmod(mod1234a56[,1],mod1234a56[,2])
	 
	 PPj123456a <- as.vector(ppj123456a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
#	 names(PPj123456a) <- unlist(mod123456a)
	 names(PPj123456a) <- as.vector(ppnames)
	 tt <- as.matrix(exp(PPj123456a),ncol=1)
	 ppj123456non <- rowsum(tt, rownames(tt)) # unique T1-T2-T5-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj123456a)
	 gc()
	 
	 } # length(modj56x)>0
	 
	 
	} # length(ind34.keep)>0 & length(modj12x) > 0
	

	
#### if have overlap T5-T6 models - combine with non-overlap T1-T2 and non-overlap T3-T4  - ppj123456nno 

if(length(ind56.keep)>0 & length(modj12x) > 0 & length(modj34x) > 0) {
 	#print("YES: length(ind56.keep)>0 & length(modj12x) > 0 & length(modj34x) > 0")
 	MODj56 <- modj56[ind56.keep]  # T1-T2 models with overlap 
 	ppj56 <- pp56a[ind56.keep] # T1-T2 models with overlap lPP1+lPP2
 	nsize56a <- nsize56all[ind56.keep,]
 	
 	
	## COMBINE NON-OVERLAP T1-T2  WITH NON-OVERLAP T3-T4 
 	ppj12a34 <- as.vector(outer(ppj12x,ppj34x,function(x,y) x+y))
 	indsize12 <- 1:nrow(nsize12ax)
 	indsize34 <- 1:nrow(nsize34ax)
 	indsize12a34 <- expand.grid(indsize12,indsize34,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
	nsize12a34 <- cbind(nsize12ax[indsize12a34[,1],], nsize34ax[indsize12a34[,2],])
	mod12a34 <- expand.grid(modj12x,modj34x,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	mod1234a <- vsortunion(mod12a34[,1],mod12a34[,2])
 	
 	## COMBINE NON-OVERLAP T1-T2 + NON-OVERLAP T3-T4 WITH OVERLAP T5-T6
 	ppnames <-  ppj123456a <- matrix(0,nrow=nrow(nsize12a34),ncol=nrow(nsize56a))
 	uniqt1234 <- unique(nsize12a34,MARGIN=1) 
	uniqt56 <- unique(nsize56a,MARGIN=1)
	 for(i in 1:nrow(uniqt1234)) {
	 	mrow <- uniqt1234[i,]
 		tmp <- apply(nsize12a34,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		for(j in 1:nrow(uniqt56)) {
 			nrow <- uniqt56[j,]
 			tmp <- apply(nsize56a,1,function(x) all(x==nrow))
 			ind2 <- which(tmp)
 			tauij <- tau6[[mrow[1,1]+1]][mrow[1,2]+1,mrow[1,3]+1,mrow[1,4]+1,nrow[1,1]+1,nrow[1,2]+1]
 			ppj123456a[ind1,ind2] <- outer(ppj12a34[ind1],ppj56[ind2],function(x,y) x+y) + tauij
 			ppn0 <- expand.grid(mod1234a[ind1],MODj56[ind2],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
 			ppn <- vjoinmod(ppn0[,1],ppn0[,2])
 			ppn <- unlist(ppn)
 			ppnames[ind1,ind2] <- matrix(ppn,nrow=length(ind1),ncol=length(ind2)) 				
# 			tau123456a[ind1,ind2] <- tau6[[mrow[1,1]]][mrow[1,2],mrow[1,3],mrow[1,4],nrow[1,1],nrow[1,2]]
 		}
	 }

# pp123456 
#	 ppj1234a56 <- outer(ppj12a34,ppj56,function(x,y) x+y)
#	 ppj123456a <-  ppj1234a56 + tau123456a

# model names
#	 mod1234a56 <- expand.grid(mod1234a,MODj56,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 #print(head(mod1234a56))
#	 mod123456a <- vjoinmod(mod1234a56[,1],mod1234a56[,2])
	 
	 PPj123456a <- as.vector(ppj123456a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
#	 names(PPj123456a) <- unlist(mod123456a)
	 names(PPj123456a) <- as.vector(ppnames)
	 tt <- as.matrix(exp(PPj123456a),ncol=1)
	 ppj123456nno <- rowsum(tt, rownames(tt)) # unique T1-T2-T5-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj123456a)
	 gc()
	} # length(ind56.keep)>0 & length(modj12x) > 0 & length(modj34x) > 0
	

#### ppj123456nnn - ppj123456nnn1 + ppj123456nnn2

# check for overlap between non-overlap T1-T2 and non-overlap T3-T4 IF non-overlap T5-T6 exists
if(length(modj12x)>0 & length(modj34x)>0 & length(modj56x)>0) {
 #print("YES: length(modj12x)>0 & length(modj34x)>0 & length(modj56x)>0")
 names1234 <- expand.grid(modj12x,modj34x,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) # all model combinations across ancestries 1-2
 modj1234 <- vjoint2(names1234[,1],names1234[,2])
 ind1234.keep <- which(!is.na(modj1234)) # overlap models; those with NA are non-overlap
 pp1234c <- expand.grid(ppj12x,ppj34x,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) 
 pp1234a <- apply(pp1234c,1,sum)
 indsize12 <- 1:nrow(nsize12ax)
 indsize34 <- 1:nrow(nsize34ax)
 indsize12a34 <- expand.grid(indsize12,indsize34,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
 nsize1234all <- cbind(nsize12ax[indsize12a34[,1],], nsize34ax[indsize12a34[,2],])
 
 if(length(ind1234.keep) > 0) {
 #print("YES: ind1234.keep > 0")
 	MODj1234 <- modj1234[ind1234.keep]  #  models with overlap 
 	ppj1234 <- pp1234a[ind1234.keep] # models with overlap lPP1+lPP2
 	nsize1234a <- nsize1234all[ind1234.keep,]
 	
 	## COMBINE OVERLAP no(T1-T2) - no(T3-T4) with no(T5-T6)
 	
 	ppnames <-  ppj123456a <- matrix(0,nrow=nrow(nsize1234a),ncol=nrow(nsize56ax))
 	uniqt1234 <- unique(nsize1234a,MARGIN=1) 
	uniqt56 <- unique(nsize56ax,MARGIN=1)
	 for(i in 1:nrow(uniqt1234)) {
	 	mrow <- uniqt1234[i,]
 		tmp <- apply(nsize1234a,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		for(j in 1:nrow(uniqt56)) {
 			nrow <- uniqt56[j,]
 			tmp <- apply(nsize56ax,1,function(x) all(x==nrow))
 			ind2 <- which(tmp)
			tauij <- tau6[[mrow[1,1]+1]][mrow[1,2]+1,mrow[1,3]+1,mrow[1,4]+1,nrow[1,1]+1,nrow[1,2]+1]
 			ppj123456a[ind1,ind2] <- outer(ppj1234[ind1],ppj56x[ind2],function(x,y) x+y) + tauij
 			ppn0 <- expand.grid(MODj1234[ind1],modj56x[ind2],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
 			ppn <- vjoinmod(ppn0[,1],ppn0[,2])
 			ppn <- unlist(ppn)
 			ppnames[ind1,ind2] <- matrix(ppn,nrow=length(ind1),ncol=length(ind2)) 							
# 			tau123456a[ind1,ind2] <- tau6[[mrow[1,1]]][mrow[1,2],mrow[1,3],mrow[1,4],nrow[1,1],nrow[1,2]]
 		}
	 }

# pp123456 
#	 ppj1234a56 <- outer(ppj1234,ppj56x,function(x,y) x+y)
#	 ppj123456a <-  ppj1234a56 + tau123456a

# model names
#	 mod1234a56 <- expand.grid(MODj1234,modj56x,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 #print(head(mod1234a56))
#	 mod123456a <- vjoinmod(mod1234a56[,1],mod1234a56[,2])
	 
	 PPj123456a <- as.vector(ppj123456a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
#	 names(PPj123456a) <- unlist(mod123456a)
	 names(PPj123456a) <- as.vector(ppnames)
	 tt <- as.matrix(exp(PPj123456a),ncol=1)
	 ppj123456nnn1 <- rowsum(tt, rownames(tt)) # unique T1-T2-T5-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
	} # ind1234.keep > 0
 
# check for overlap between non-overlap (no(T1-T2) - no(T3-T4)) and non-overlap T5-T6 
 if(length(ind1234.keep) < nrow(names1234) & length(modj56x)>0 ) {
 		if(length(ind1234.keep)>0) {
 		modj1234a <- vsortunion(names1234[,1],names1234[,2])
		modj1234x <- modj1234a[-ind1234.keep]  # T1-T2 models with non-overlap 
 		ppj1234x <- pp1234a[-ind1234.keep] # T1-T2 models with non-overlap lPP1+lPP2
 		nsize1234ax <- nsize1234all[-ind1234.keep,]
 		} else {
 		modj1234x <- vsortunion(names1234[,1],names1234[,2])
 		ppj1234x <- pp1234a # T1-T2 models with non-overlap lPP1+lPP2
 		nsize1234ax <- nsize1234all
 		}
 
 names123456 <- expand.grid(modj1234x,modj56x,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) # all model combinations across ancestries 1-2
 modj123456 <- vjoint2F(names123456[,1],names123456[,2])
 ind123456.keep <- which(!is.na(modj123456)) # overlap models; those with NA are non-overlap
 pp123456c <- expand.grid(ppj1234x,ppj56x,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) 
 pp123456a <- apply(pp123456c,1,sum)
 indsize1234 <- 1:nrow(nsize1234ax)
 indsize56 <- 1:nrow(nsize56ax)
 indsize1234a56 <- expand.grid(indsize1234,indsize56,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
 nsize123456all <- cbind(nsize1234ax[indsize1234a56[,1],], nsize56ax[indsize1234a56[,2],])

 
 if(length(ind123456.keep) > 0) {
 #print("YES: ind123456.keep > 0")
 	MODj123456 <- modj123456[ind123456.keep]  #  models with overlap 
 	ppj123456 <- pp123456a[ind123456.keep] # models with overlap lPP1+lPP2
 	nsize123456a <- nsize123456all[ind123456.keep,]
 	 	
 	ppj123456a <- numeric(nrow(nsize123456a))
 	uniqt123456 <- unique(nsize123456a,MARGIN=1) 
	for(i in 1:nrow(uniqt123456)) {
	 	mrow <- uniqt123456[i,]
 		tmp <- apply(nsize123456a,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		tauij <- tau6[[mrow[1,1]+1]][mrow[1,2]+1,mrow[1,3]+1,mrow[1,4]+1,mrow[1,5]+1,mrow[1,6]+1]
 		ppj123456a[ind1] <- ppj123456[ind1] + tauij
# 		tau123456a[ind1] <- tau6[[mrow[1,1]]][mrow[1,2],mrow[1,3],mrow[1,4],mrow[1,5],mrow[1,6]]		
	 }

# pp123456 
#	 ppj123456a <-  ppj123456 + tau123456a

# model names
	 
	 PPj123456a <- as.vector(ppj123456a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
	 names(PPj123456a) <- MODj123456

	 tt <- as.matrix(exp(PPj123456a),ncol=1)
	 ppj123456nnn2 <- rowsum(tt, rownames(tt)) # unique T1-T2-T5-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
	} # ind123456.keep > 0
 
} # length(ind1234.keep) < nrow(names1234) & length(modj56x)>0

} # length(modj12x)>0 & length(modj34x)>0 & length(modj56x)>0

 
all.list <- vector("list",5)
all.list[[1]] <- ppj123456ooo
all.list[[2]] <- ppj123456oon
all.list[[3]] <- ppj123456ono
all.list[[4]] <- ppj123456onn
all.list[[5]] <- ppj123456noo 
all.list[[6]] <- ppj123456non 
all.list[[7]] <- ppj123456nno 
all.list[[8]] <- ppj123456nnn1 
all.list[[9]] <- ppj123456nnn2 

#print(lapply(all.list,head))
  
have <- which(sapply(all.list,function(x) !is.null(x)))

#print(have)

	if(length(have)==1) ppj <- all.list[[have]]
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
	 			if(length(have) > 3) {
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
	 				if(length(have) == 5) {
	 					ppj123 <- ppj
						ppj123b <- all.list[[have[5]]]
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
 			} 				
 		nameppj <- rownames(ppj)
 		ppj <- data.frame(ppj,row.names=nameppj)
 		rm(all.list)
 		rm(ppj123456ooo,ppj123456oon,ppj123456ono,ppj123456onn,ppj123456noo,ppj123456non,ppj123456nno,ppj123456nnn1,ppj123456nnn2)
 		
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
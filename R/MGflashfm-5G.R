
### 5 study functions

calctauTA5 <- function(n1,n2,n3,n4,n5,nsnps) {
    ns <- sort(c(n1,n2,n3,n4,n5),decreasing=FALSE) # ns[1] is min, ns[5] is max
    num <- choose(nsnps,ns[1])*choose(nsnps,ns[2])*choose(nsnps,ns[3])*choose(nsnps,ns[4])
    den <- choose(nsnps,ns[1])*choose(nsnps,ns[2])*choose(nsnps,ns[3])*choose(nsnps,ns[4]) - 
    choose(nsnps-ns[5],ns[4])*choose(nsnps-ns[5]-ns[4],ns[3])*choose(nsnps-ns[5]-ns[4]-ns[3],ns[2])*choose(nsnps-ns[5]-ns[4]-ns[3]-ns[2],ns[1])
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


######## 5 studies ###########
TA5flashfm1 <- function(PPn,nsnpspermodel,SS,nsnps,snps,cred=0.99) {
 

	nmax <- max(sapply(nsnpspermodel,max))
    tau5 <- vector("list",nmax+1)
    for(i in 0:nmax) tau5[[i+1]] <- array(0,dim=rep(nmax+1,4))
    for(i in 0:nmax){ 
    	for(j in 0:nmax) { 
    		for(k in 0:nmax) { 
    			for(l in 0:nmax) {
    				for(m in 0:nmax) {
    				tau5[[i+1]][j+1,k+1,l+1,m+1] <- calctauTA5(i,j,k,l,m,nsnps)
    			} } }}}
				
namesPPn <- SS
 ## numeric version of namesPPn for speed
   mstr <- lapply(namesPPn, function(ss) {
        lapply(ss, function(x) as.integer(factor(x, levels = snps)))
    })
    names(mstr) <- NULL
 namesPPn <- mstr
 rm(mstr)
 gc()
 
 #names(nsnpspermodel[[5]]) <- NULL
 
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

 
 ppj12534oo <- ppj12534on <- ppj12534no  <- ppj12345nn <- ppj12345nno <- c()


##### PART 2 - if have overlap T1-T2 models
if(length(ind12.keep)>0) {
 	MODj12 <- modj12[ind12.keep]  # T1-T2 models with overlap 
 	ppj12 <- pp12a[ind12.keep] # T1-T2 models with overlap lPP1+lPP2
 	nsize12a <- nsize12all[ind12.keep,]
  # COMBINE WITH ALL T5 MODELS	
    ppj12a5 <- as.vector(outer(ppj12,PPn[[5]],function(x,y) x+y))
 	indsize12 <- 1:nrow(nsize12a)
 	indsize5 <- 1:length(nsnpspermodel[[5]])
 	indsize12a5 <- expand.grid(indsize12,indsize5,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
	nsize12a5 <- cbind(nsize12a[indsize12a5[,1],], nsnpspermodel[[5]][indsize12a5[,2]])
	mod12a5 <- expand.grid(MODj12,namesPPn[[5]],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	mod125a <- vsortunion(mod12a5[,1],mod12a5[,2])
 
 
	### and have overlap T3-T4 models 
	if(length(ind34.keep)>0) {
 		MODj34 <- modj34[ind34.keep]  # T3-T4 models with overlap 
 		ppj34 <- pp34a[ind34.keep] # T3-T4 models with overlap lPP1+lPP2

	## JOIN OVERLAP T1-T2 + T5 WITH OVERLAP T3-T4 
#tau 
 	nsize34a <- nsize34all[ind34.keep,]
 
 	ppnames <- ppj12534a <- matrix(0,nrow=nrow(nsize12a5),ncol=nrow(nsize34a))
 	uniqt125 <- unique(nsize12a5,MARGIN=1) 
	uniqt34 <- unique(nsize34a,MARGIN=1) 
	 for(i in 1:nrow(uniqt125)) {
	 	mrow <- uniqt125[i,]
 		tmp <- apply(nsize12a5,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		for(j in 1:nrow(uniqt34)) {
 			nrow <- uniqt34[j,]
 			tmp <- apply(nsize34a,1,function(x) all(x==nrow))
 			ind2 <- which(tmp)
 			tauij <- tau5[[mrow[1,1]+1]][mrow[1,2]+1,nrow[1,1]+1,nrow[1,2]+1,mrow[1,3]+1] # shift indices by 1 to allow for 0
 			ppj12534a [ind1,ind2] <- outer(ppj12a5[ind1],ppj34[ind2],function(x,y) x+y) + tauij
 			ppn0 <- expand.grid(mod125a[ind1],MODj34[ind2],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
 			ppn <- vjoinmod(ppn0[,1],ppn0[,2])
 			ppn <- unlist(ppn)
 			ppnames[ind1,ind2] <- matrix(ppn,nrow=length(ind1),ncol=length(ind2)) 	
# 			tau12534a[ind1,ind2] <- tau5[[mrow[1,1]]][mrow[1,2],nrow[1,1],nrow[1,2],mrow[1,3]]
 		}
	 }
# pp1234 
#	 ppj125a34 <- outer(ppj12a5,ppj34,function(x,y) x+y)
#	 ppj12534a <-  ppj125a34 + tau12534a

# model names
#	 mod125a34 <- expand.grid(mod125a,MODj34,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
#	 mod12534a <- vjoinmod(mod125a34[,1],mod125a34[,2])
	 
	 PPj12534a <- as.vector(ppj12534a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
#	 names(PPj12534a) <- unlist(mod12534a)
	 names(PPj12534a) <- as.vector(ppnames)
	 tt <- as.matrix(exp(PPj12534a),ncol=1)
	 ppj12534oo <- rowsum(tt, rownames(tt)) # unique T1-T2-T5-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj12534a,tmp,ppj12534a)
	 gc()
	
	} # length(ind34.keep)>0 
	 
	## JOIN OVERLAP T1-T2 +T5 WITH NON-OVERLAP T3-T4 (if exists) 
	ppj1234on5 <- c()
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
 		

 	ppnames <- ppj12534a <- matrix(0,nrow=nrow(nsize12a5),ncol=nrow(nsize34ax))
 	uniqt125 <- unique(nsize12a5,MARGIN=1) 
	uniqt34 <- unique(nsize34ax,MARGIN=1)
	 for(i in 1:nrow(uniqt125)) {
	 	mrow <- uniqt125[i,]
 		tmp <- apply(nsize12a5,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		for(j in 1:nrow(uniqt34)) {
 			nrow <- uniqt34[j,]
 			tmp <- apply(nsize34ax,1,function(x) all(x==nrow))
 			ind2 <- which(tmp)
 			tauij <- tau5[[mrow[1,1]+1]][mrow[1,2]+1,nrow[1,1]+1,nrow[1,2]+1,mrow[1,3]+1]
 			ppj12534a[ind1,ind2] <- outer(ppj12a5[ind1],ppj34x[ind2],function(x,y) x+y) + tauij
 			ppn0 <- expand.grid(mod125a[ind1],modj34x[ind2],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE)
 			ppn <- vjoinmod(ppn0[,1],ppn0[,2])
 			ppn <- unlist(ppn)
 			ppnames[ind1,ind2] <- matrix(ppn,nrow=length(ind1),ncol=length(ind2)) 
# 			tau12534a[ind1,ind2] <- tau5[[mrow[1,1]]][mrow[1,2],nrow[1,1],nrow[1,2],mrow[1,3]]
 		}
	 }
# pp1234 
#	 ppj125a34 <- outer(ppj12a5,ppj34x,function(x,y) x+y)
#	 ppj12534a <-  ppj125a34 + tau12534a

# model names
#	 mod125a34 <- expand.grid(mod125a,modj34x,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
#	 mod12534a <- vjoinmod(mod125a34[,1],mod125a34[,2])
	 
	 PPj12534a <- as.vector(ppj12534a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
#	 names(PPj12534a) <- unlist(mod12534a)
	 names(PPj12534a) <- as.vector(ppnames)
	 tt <- as.matrix(exp(PPj12534a),ncol=1)
	 ppj12534on <- rowsum(tt, rownames(tt)) # unique T1-T2-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj12534a,tmp,ppj12534a)
	 gc()
		} ## if(length(ind34.keep)<length(modj34)) 

	 }  # length(ind12.keep)>0
	 	
	
	 
## JOIN NON-OVERLAP T1-T2 (IF EXISTS) WITH OVERLAP T3-T4 (if exists)
#tau 
 	if(length(ind12.keep)<length(modj12) & length(ind34.keep)>0 ) {
 	
 	MODj34 <- modj34[ind34.keep]  # T3-T4 models with overlap 
 	ppj34 <- pp34a[ind34.keep] # T3-T4 models with overlap lPP1+lPP2
 	nsize34a <- nsize34all[ind34.keep,]
 	
 	# combine T1+T2 with T5
 	if(length(ind12.keep) >0){ 
		modj12a <- vsortunion(names12[,1],names12[,2])
		MODj12 <- modj12a[-ind12.keep]  # T1-T2 models with non-overlap 
 		ppj12 <- pp12a[-ind12.keep] # T1-T2 models with non-overlap lPP1+lPP2
 		nsize12a <- nsize12all[-ind12.keep,]
 		# COMBINE WITH ALL T5 MODELS	
	    ppj12a5 <- as.vector(outer(ppj12,PPn[[5]],function(x,y) x+y))
    	indsize12 <- 1:nrow(nsize12a)
	 	indsize5 <- 1:length(nsnpspermodel[[5]])
	 	indsize12a5 <- expand.grid(indsize12,indsize5,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
		nsize12a5 <- cbind(nsize12a[indsize12a5[,1],], nsnpspermodel[[5]][indsize12a5[,2]])
		mod12a5 <- expand.grid(MODj12,namesPPn[[5]],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
		mod125a <-vsortunion(mod12a5[,1],mod12a5[,2])
 	} else {
 		MODj12 <- vsortunion(names12[,1],names12[,2])
 		ppj12 <- pp12a
 		nsize12a <- nsize12all
 		# COMBINE WITH ALL T5 MODELS	
	    ppj12a5 <- as.vector(outer(ppj12,PPn[[5]],function(x,y) x+y))
	    indsize12 <- 1:nrow(nsize12a)
	 	indsize5 <- 1:length(nsnpspermodel[[5]])
	 	indsize12a5 <- expand.grid(indsize12,indsize5,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
		nsize12a5 <- cbind(nsize12a[indsize12a5[,1],], nsnpspermodel[[5]][indsize12a5[,2]])
		mod12a5 <- expand.grid(MODj12,namesPPn[[5]],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
		mod125a <- vsortunion(mod12a5[,1],mod12a5[,2]) 
 	}	
 	
 	ppnames <- ppj12534a <- matrix(0,nrow=nrow(nsize12a5),ncol=nrow(nsize34a))
 	uniqt125 <- unique(nsize12a5,MARGIN=1) 
	uniqt34 <- unique(nsize34a,MARGIN=1)
	 for(i in 1:nrow(uniqt125)) {
	 	mrow <- uniqt125[i,]
 		tmp <- apply(nsize12a5,1,function(x) all(x==mrow))
 		ind1 <- which(tmp)
 		for(j in 1:nrow(uniqt34)) {
 			nrow <- uniqt34[j,]
 			tmp <- apply(nsize34a,1,function(x) all(x==nrow))
 			ind2 <- which(tmp)
 			tauij <- tau5[[mrow[1,1]+1]][mrow[1,2]+1,nrow[1,1]+1,nrow[1,2]+1,mrow[1,3]+1]
 			ppj12534a[ind1,ind2] <- outer(ppj12a5[ind1],ppj34[ind2],function(x,y) x+y) + tauij
 			ppn0 <- expand.grid(mod125a[ind1],MODj34[ind2],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
 			ppn <- vjoinmod(ppn0[,1],ppn0[,2])
 			ppn <- unlist(ppn)
 			ppnames[ind1,ind2] <- matrix(ppn,nrow=length(ind1),ncol=length(ind2)) 	
# 			tau12534a[ind1,ind2] <- tau5[[mrow[1,1]]][mrow[1,2],nrow[1,1],nrow[1,2],mrow[1,3]]
 		}
	 }
# pp1234 
#	 ppj125a34 <- outer(ppj12a5,ppj34,function(x,y) x+y)
#	 ppj12534a <-  ppj125a34 + tau12534a

# model names
#	 mod125a34 <- expand.grid(mod125a,MODj34,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
#	 mod12534a <- vjoinmod(mod125a34[,1],mod125a34[,2])
	 
	 PPj12534a <- as.vector(ppj12534a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
#	 names(PPj12534a) <- unlist(mod12534a)
	names(PPj12534a) <- as.vector(ppnames)
	 tt <- as.matrix(exp(PPj12534a),ncol=1)
	 ppj12534no <- rowsum(tt, rownames(tt)) # unique T1-T2-T5-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj12534a,tmp,ppj12534a)
	 gc()
		} ### if(length(ind12.keep)<length(MODj12)) & length(ind34.keep)>0






## CHECK IF ANY OVERLAP BETWEEN NON-OVERLAP T1-T2 AND NON-OVERLAP T3-T4 and if not then check overlap with T5
if(length(ind12.keep)<length(modj12) & length(ind34.keep)<length(modj34) ) {	 
	
	if(length(ind12.keep) >0){
		modj12a <- vsortunion(names12[,1],names12[,2])
		modj12x <- modj12a[-ind12.keep]  # T1-T2 models with non-overlap 
 		ppj12 <- pp12a[-ind12.keep] # T1-T2 models with non-overlap lPP1+lPP2
 		nsize12a <- nsize12all[-ind12.keep,]
 	} else {
 		modj12x <- vsortunion(names12[,1],names12[,2])
 		ppj12 <- pp12a
 		nsize12a <- nsize12all
 	}
 	
 	if(length(ind34.keep) >0){
 		modj34a <- vsortunion(names34[,1],names34[,2])
		modj34x <- modj34a[-ind34.keep]  # T3-T4 models with non-overlap 
 		ppj34 <- pp34a[-ind34.keep] # T3-T4 models with non-overlap lPP1+lPP2
		nsize34a <- nsize34all[-ind34.keep,]
	} else {
		modj34x <- vsortunion(names34[,1],names34[,2])
 		ppj34 <- pp34a
		nsize34a <- nsize34all
	}
 
	
# model names
	 mod12a34 <- expand.grid(modj12x,modj34x,stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 modj1234 <- vjoint2(mod12a34[,1],mod12a34[,2])
	 ind1234keep <- which(!is.na(modj1234)) # OVERLAP BETWEEN NON-OVERLAP T1-T2 AND NON-OVERLAP T3-T4
	
	 pp1234c <- expand.grid(ppj12,ppj34,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE) 
	 pp1234a <- apply(pp1234c,1,sum)
	 ind12 <- 1:nrow(nsize12a)
	 ind34 <- 1:nrow(nsize34a)
	 ind1234all <- expand.grid(ind12,ind34,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
	 nsize1234all <- cbind(nsize12a[ind1234all[,1],], nsize34a[ind1234all[,2],])
	  if(length(ind1234keep) > 0) {
	    MODj1234 <- modj1234[ind1234keep]  #overlapping T1-T2 models with non-overlap  and T3-T4 non-overlap
 		ppj1234 <- pp1234a[ind1234keep] 
 		nsize1234a <- nsize1234all[ind1234keep,]
	   # COMBINE WITH ALL T5 MODELS	
	    ppj1234a5 <- outer(ppj1234,PPn[[5]],function(x,y) x+y)
	    ind1234 <- 1:nrow(nsize1234a)
		ind5 <- 1:length(nsnpspermodel[[5]])
		ind1234a5 <- expand.grid(ind1234,ind5,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
		nsize1234a5 <- cbind(nsize1234a[ind1234a5[,1],], nsnpspermodel[[5]][ind1234a5[,2]])
		mod1234a5 <- expand.grid(MODj1234,namesPPn[[5]],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
		mod12345a <- vjoinmod(mod1234a5[,1],mod1234a5[,2])
	  

		ppj12345a <- matrix(0,nrow=nrow(nsize1234a),ncol=length(nsnpspermodel[[5]]))
 		uniqt1234 <- unique(nsize1234a,MARGIN=1) 
		uniqt5 <- unique(nsnpspermodel[[5]])
		 for(i in 1:nrow(uniqt1234)) {
		 	mrow <- uniqt1234[i,]
 			tmp <- apply(nsize1234a,1,function(x) all(x==mrow))
 			ind1 <- which(tmp)
 			for(j in uniqt5) {
 				ind2 <- which(nsnpspermodel[[5]] == j)
 				tauij <- tau5[[mrow[1,1]+1]][mrow[1,2]+1,mrow[1,3]+1,mrow[1,4]+1,j+1]
 				ppj12345a[ind1,ind2] <- outer(ppj1234[ind1],PPn[[5]][ind2],function(x,y) x+y) + tauij
# 				tau12345a[ind1,ind2] <- tau5[[mrow[1,1]]][mrow[1,2],mrow[1,3],mrow[1,4],j]
 			}
		 }
# pp12345 
#	 ppj1234a5 <- outer(ppj1234,PPn[[5]],function(x,y) x+y)
#	 ppj12345a <-  ppj1234a5 + tau12345a
	 
	 PPj12345a <- as.vector(ppj12345a) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
	 names(PPj12345a) <- unlist(mod12345a)
	 tt <- as.matrix(exp(PPj12345a),ncol=1)
	 ppj12345nn <- rowsum(tt, rownames(tt)) # unique T1-T2-T3-T4 models with T1-T2 overlap and T3-T4 overlap lPP1+lPP2+lPP3+tau
 
	 rm(tt,PPj12345a,mod12345a,tmp,mod1234a5,ppj1234a5,ppj1234)
	 gc()
	 } # length(ind1234keep) > 0
	 
	 
	 if(length(ind1234keep)< length(modj1234)) {
	  
	  if(length(ind1234keep) >0){
	    modj1234 <- vsortunion(mod12a34[,1],mod12a34[,2])
		MODj1234 <- modj1234[-ind1234keep]  #overlapping T1-T2 models with non-overlap  and T3-T4 non-overlap
 		ppj1234 <- pp1234a[-ind1234keep] 
 		nsize1234a <- nsize1234all[-ind1234keep,]
 	} else {
 		MODj1234 <- vsortunion(mod12a34[,1],mod12a34[,2]) #overlapping T1-T2 models with non-overlap  and T3-T4 non-overlap
 		ppj1234 <- pp1234a 
 		nsize1234a <- nsize1234all
 	}
 	
 	 #### CHECK OVERLAP of nn1234 models WITH ALL T5 MODELS	
	 mod1234a5 <- expand.grid(MODj1234,namesPPn[[5]],stringsAsFactors = FALSE,KEEP.OUT.ATTRS = FALSE) 
	 modj1234a5 <- vjoint2F(mod1234a5[,1],mod1234a5[,2])
	 ind12345keep <- which(!is.na(modj1234a5)) # OVERLAP BETWEEN NON-OVERLAP T1-T2 with NON-OVERLAP T3-T4 and T5
	 
	 ind1234 <- 1:nrow(nsize1234a)
	 ind5 <- 1:length(nsnpspermodel[[5]])
	 ind12345 <- expand.grid(ind1234,ind5,KEEP.OUT.ATTRS = FALSE,stringsAsFactors = FALSE)
	 nsize12345all <- cbind(nsize1234a[ind12345[,1],] ,nsnpspermodel[[5]][ind12345[,2]])
	 
 	 
	  if(length(ind12345keep) > 0) {
	    MODj12345 <- modj1234a5[ind12345keep]  #non-overlapping T1-T2 models with non-overlap  T3-T4  that overlap T5
 		nsize12345a <- nsize12345all[ind12345keep,]
	   
		tau12345a <- matrix(0,nrow=nrow(nsize1234a),ncol=length(nsnpspermodel[[5]]))
 		uniqt1234 <- unique(nsize1234a,MARGIN=1) 
		uniqt5 <- unique(nsnpspermodel[[5]])
		 for(i in 1:nrow(uniqt1234)) {
		 	mrow <- uniqt1234[i,]
 			tmp <- apply(nsize1234a,1,function(x) all(x==mrow))
 			ind1 <- which(tmp)
 			for(j in uniqt5) {
 				ind2 <- which(nsnpspermodel[[5]] == j)
 				tau12345a[ind1,ind2] <- tau5[[mrow[1,1]+1]][mrow[1,2]+1,mrow[1,3]+1,mrow[1,4]+1,j+1]
 			}
		 }
# pp12345 
	 ppj1234a5 <- outer(ppj1234,PPn[[5]],function(x,y) x+y)
	 ppj12345a <-  ppj1234a5 + tau12345a
	 
	 PPj12345a <- as.vector(ppj12345a[ind12345keep]) # c(ppj12a3[,1],ppj12a3[,2],...,ppj12a3[,ncol]), so matches order in expand.grid
	 names(PPj12345a) <-  unlist(MODj12345)
	 tt <- as.matrix(exp(PPj12345a),ncol=1)
	 ppj12345nno <- rowsum(tt, rownames(tt)) 
 
	 rm(tt,PPj12345a,tmp,mod1234a5,ppj1234a5,ppj1234,tau12345a)
	 gc()
	 } # length(ind12345keep) > 0

	 } # length(ind1234keep)< nrow(modj1234)
	 

} # length(ind12.keep)<length(MODj12) & length(ind34.keep)<length(MODj34)


#ppj12534oo <- ppj12534on <- ppj12534no  <- ppj12345nn <- ppj12345nno <- c()

all.list <- vector("list",5)
all.list[[1]] <- ppj12534oo
all.list[[2]] <- ppj12534on
all.list[[3]] <- ppj12534no
all.list[[4]] <- ppj12345nn
all.list[[5]] <- ppj12345nno
 
have <- which(sapply(all.list,function(x) !is.null(x)))

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
 		rm(ppj12534oo,ppj12534on,ppj12534no,ppj12345nn,ppj12345nno)
 		
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
 
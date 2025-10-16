tagexpand.mod <- function (snpmod, taglist) 
{
    if (snpmod == "1") {
        df <- data.frame(str = snpmod, snps = snpmod, size = 0, 
            tag = FALSE, stringsAsFactors = FALSE)
    }
    else {
        snps <- unlist(strsplit(snpmod, "%", fixed = TRUE))
        ns <- length(snps)
        tsnps <- vector("list", ns)
        for (i in 1:ns) {
            # exact match instead of grep
            idx <- which(vapply(taglist, function(el) any(el$snps == snps[i]), logical(1)))
            if (length(idx) == 0L) {
                stop(sprintf("No exact match for '%s' found in taglist", snps[i]))
            }
            tsnps[[i]] <- unique(unlist(taglist[idx]))
        }
        emods <- expand.grid(tsnps, stringsAsFactors = FALSE)
        out <- apply(emods, 1, function(x) {
            paste0(x, collapse = "%")
        })
        Imod <- which(out == snpmod)
        istag <- rep(TRUE, length(out))
        istag[Imod] <- FALSE
        df <- data.frame(str = snpmod, snps = out, size = ns, 
            tag = istag, stringsAsFactors = FALSE)
    }
    return(df)
}



cor.refdata2_mod <- function(corX, beta, MAF = NULL, r2 = 0.99) {
  # corX: precomputed correlation matrix (SNP x SNP)
  # beta: named vector of effect sizes for SNPs
  # MAF: optional named vector of minor allele frequencies
  # r2: r-squared threshold for tagging
  
  # call original tagSNP function
  gmat2t <- tagSNP(corX, threshold = sqrt(r2))
  
  # update tag SNP based on |beta| and optional MAF
  for (i in seq_along(gmat2t)) {
    snps_in_bin <- gmat2t[[i]]$snps
    if (length(snps_in_bin) > 1) {
      # pick SNP(s) with largest |beta|
      abs_beta <- abs(beta[snps_in_bin])
      tag_candidates <- snps_in_bin[abs_beta == max(abs_beta)]
      
      # tie-break with smallest MAF if available
      if (!is.null(MAF) && length(tag_candidates) > 1) {
        tag_idx <- tag_candidates[which.min(MAF[tag_candidates])]
      } else {
        tag_idx <- tag_candidates[1]
      }
      
      gmat2t[[i]]$tagsnp <- tag_idx
    } else {
      gmat2t[[i]]$tagsnp <- snps_in_bin
    }
  }
  
  # extract tag SNPs
  tg <- unlist(gmat2t)
  tg <- tg[names(tg) == "tagsnp"]
  
  # create reference correlation matrix for tag SNPs
  refG <- corX[tg, tg, drop = FALSE]
  
  return(list(refG = refG, taglist = gmat2t))
}
###

#' @title Approximate effective sample size from GWAS of related samples (copied from flashfm)
#' @param raf vector of reference (or minor) allele frequencies for all SNPs in the GWAS
#' @param seB vector of standard errors of SNP effect estimates, in same order as SNPs in raf
#' @param Vy trait variance; default is 1 under assumption trait is transformed to a standard Normal distribution
#' @return Estimate of effective sample size 
#' @author Jenn Asimit
#' @export
 Neff <- function(raf,seB,Vy=1) {
  keep <- which(!is.na(seB))
  raf <- raf[keep]
  seB <- seB[keep]
  Vb <- seB^2
  Vx <- 2*raf*(1-raf)
  Nhat <- Vy/(Vx*Vb)
  Ne <- round(median(Nhat))
  return(Ne)
  }




JAMmaxcv <- function(BETA,Vy,refG,mafs.ref,N,save.path,maxcv=2,jam.nM.iter=1,extra.java.arguments=NULL) {
 jam.results <- R2BGLiMS::JAM(marginal.betas = BETA, trait.variance = Vy, 
            cor.ref = refG, mafs.ref = mafs.ref, model.space.priors = list(a = 1, 
                b = length(BETA), Variables = names(BETA)), max.model.dim = maxcv, 
            n = N, xtx.ridge.term = 0.01, save.path = save.path, n.mil.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments) 
 topmods <- R2BGLiMS::TopModels(jam.results, n.top.models = 1000)
        if (is.null(ncol(topmods))  ) {
            stop("A single model was selected with PP=1. This may mean no convergence because a causal variant  is missing from the data and it has no tags in your data.")
        }
 return(topmods)
 }       

JAM.tries.maxcv <- function(BETA,Vy,refG,mafs.ref,N,save.path,maxcv=2, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1,extra.java.arguments=NULL){
  
 
  JAM_output <- NULL
  while (maxcv_autocheck == TRUE){
    
    tryCatch({
      
      JAM_output <-  JAMmaxcv(BETA=BETA,Vy=Vy,refG=refG,mafs.ref=mafs.ref,N=N,save.path=save.path,maxcv=maxcv,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments) 
        
      maxcv_autocheck = FALSE
      print(paste0('Completed the process when ..., maxcv == ', maxcv))
      
    }, error = function(e){
      
      print(paste0('Keep trying maxcv by adding 1 ..., maxcv == ', maxcv+1))
      
    })
    
    maxcv = maxcv+1
    
    if (maxcv > maxcv_stop & is.null(JAM_output)){
      stop("The maxcv_stop reached, no result! Try increasing maxcv_stop.")
      maxcv_autocheck = FALSE
    }
  
  }
  return(JAM_output)
}

            
JAMmulti <- function (gwas.list, corX, ybar, Vy, N, r2 = 0.99, save.path, maxcv=2, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1,extra.java.arguments=NULL) 
{
	 M <- length(gwas.list)
	if (is.null(names(gwas.list))) {
        ts <- paste0("T", 1:M)
        names(ybar) <- ts
    } else {
        ts <- names(gwas.list)
    }
	beta1 <- lapply(gwas.list, function(x) {b <- x[,"beta"]; names(b) <- x[,"rsID"]; b})
	names(beta1) <- ts
	raf1 <- lapply(gwas.list, function(x) {b <- x[,"EAF"]; names(b) <- x[,"rsID"]; b})
	names(raf1) <- ts

    Nlist <- makeNlist.rel(Ne = N)
    N <- diag(Nlist$Nqq)
    
    snps <- Reduce(intersect,lapply(beta1,names))
    snps <- intersect(snps, colnames(corX))
    for (i in 1:M) {
    	beta1[[i]] <- beta1[[i]][snps]
    	raf1[[i]] <- raf1[[i]][snps]
    	}
    corX <- corX[snps, snps]    
    nsnps <- length(snps)
    refG <- cor.refdata.fn(corX, r2)
    corX2 <- corX^2
    taglist <- tagSNP(corX2, threshold = r2)
    out <- list(SM = NULL, mbeta = NULL, Nlist = Nlist)
    out$SM <- vector("list", M)
    names(out$SM) <- ts
    out$mbeta <- vector("list", M)
    names(out$mbeta) <- ts
    dd <- vector("list", M)
    SSy <- vector("list", M)
    Sxy <- vector("list", M)
    for (j in 1:M) {
    	raf <- raf1[[j]]
    	BETA <- beta1[[j]][colnames(refG)]
    	maf <- raf * (raf <= 0.5) + (1 - raf) * (raf > 0.5)
        mafs.ref <- maf[colnames(refG)]
        covX <- cor2cov(corX, sd = sqrt(2 * raf * (1 - raf)))
        covX <- covX[snps, snps]  
		 topmods <- JAM.tries.maxcv(BETA=BETA,Vy=Vy[j],refG=refG,mafs.ref=mafs.ref,N=N[j],save.path=save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)
#        jam.results <- R2BGLiMS::JAM(marginal.betas = BETA, trait.variance = Vy[j], 
#            cor.ref = refG, mafs.ref = mafs.ref, model.space.priors = list(a = 1, 
#                b = length(BETA), Variables = names(BETA)), max.model.dim = maxcv, 
#            n = N[j], xtx.ridge.term = 0.01, save.path = save.path)
#        topmods <- R2BGLiMS::TopModels(jam.results, n.top.models = 1000)
#        if (is.null(ncol(topmods))) {
#            stop("A single model was selected with PP=1. This may mean no convergence because a causal variant is missing from the data and it has no tags in your data.")
#        }
        binout <- as.matrix(topmods[, -ncol(topmods)])
        colnames(binout) <- colnames(topmods)[-ncol(topmods)]
        snpmods <- apply(binout, 1, mod.fn)
        nmod <- apply(binout, 1, sum)
        PP <- topmods[, ncol(topmods)]
        snpPP <- data.frame(rank = 1:length(nmod), size = nmod, 
            logPP = log(PP), PP = PP, str = snpmods, snps = snpmods, 
            stringsAsFactors = FALSE)
        snpPP <- snpPP[order(snpPP$PP, decreasing = TRUE), ]
        expmods <- rlist::list.stack(lapply(snpPP$snps, tagexpand.mod, 
            taglist = taglist))
        wh <- which(duplicated(expmods$snps))
        if (length(wh) > 0) {
            expmods <- expmods[-wh, ]
        }
        row.names(expmods) <- expmods$snps
        check <- sapply(strsplit(expmods[, 2], "%"), function(x) length(x) > 
            length(unique(x)))
        if (sum(check) > 0) 
            expmods <- expmods[-which(check), ]
        mbeta <- lapply(expmods[, 2], multibeta, beta1[[j]], 
            covX, N = N[j], ybar = ybar[j], is.snpmat = FALSE, 
            raf = raf)
        names(mbeta) <- expmods[, 2]
        SSy[[j]] <- Vy[j] * (N[j] - 1) + N[j] * ybar[j]^2
        Vx <- diag(covX)
        Mx <- 2 * raf
        Sxy[[j]] <- c(Sxy.hat(beta1 = beta1[[j]], Mx = Mx, N = N[j], 
            Vx = Vx, muY = ybar[j]), `1` = ybar[j] * N[j])
        names(Sxy[[j]])[length(Sxy[[j]])] <- "one"
        lABF <- sapply(expmods$snps, calcABF, mbeta, SSy = SSy[[j]], 
            Sxy = Sxy[[j]], Vy = Vy[j], N = N[j])
        names(lABF) <- expmods$snps
        wh <- which(expmods$snps == "1")
        if (!length(wh)) {
            dd[[j]] <- data.frame(model = c("1", expmods$snps), 
                tag = c(FALSE, expmods$tag), lBF = c(0, lABF), 
                stringsAsFactors = FALSE)
            l1 <- multibeta("1", beta1[[j]], covX, N = N[j], 
                ybar = ybar[j], is.snpmat = FALSE, raf = raf)
            mbeta <- rlist::list.append(mbeta, `1` = l1)
        }
        else {
            dd[[j]] <- data.frame(model = expmods$snps, tag = expmods$tag, 
                lBF = lABF, stringsAsFactors = FALSE)
        }
        SM <- makesnpmod(dd[[j]], expected = 2, nsnps = nsnps)
        out$SM[[j]] <- SM
        out$mbeta[[j]] <- mbeta
        names(out$SM) <- ts
        names(out$mbeta) <- ts
        out$SM[[j]] <- best.models.cpp(out$SM[[j]], maxmod = 1000)[[1]]
        out$SM[[j]] <- PP2snpmod(out$SM[[j]])
    }
    out$Nlist <- Nlist
    out$nsnps <- nsnps
    out$Gmat = covX
    out$beta1.list = beta1
    out$raf = raf1[[which.max(N)]]
    return(out)
}


JAMmulti2 <- function (gwas.list, corX, ybar, Vy, N, r2 = 0.99, save.path, maxcv=2, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1,extra.java.arguments=NULL) 
{
	 M <- length(gwas.list)
	if (is.null(names(gwas.list))) {
        ts <- paste0("T", 1:M)
        names(ybar) <- ts
    } else {
        ts <- names(gwas.list)
    }
	beta1 <- lapply(gwas.list, function(x) {b <- x[,"beta"]; names(b) <- x[,"rsID"]; b})
	names(beta1) <- ts
	raf1 <- lapply(gwas.list, function(x) {b <- x[,"EAF"]; names(b) <- x[,"rsID"]; b})
	names(raf1) <- ts

    Nlist <- makeNlist.rel(Ne = N)
    N <- diag(Nlist$Nqq)
    
    snps <- Reduce(intersect,lapply(beta1,names))
    snps <- intersect(snps, colnames(corX))
    for (i in 1:M) {
    	beta1[[i]] <- beta1[[i]][snps]
    	raf1[[i]] <- raf1[[i]][snps]
    	}
    corX <- as.matrix(corX)	
    corX <- corX[snps, snps]    
    nsnps <- length(snps)
    corX <- as.matrix(corX)	
    
    maf_vec <- raf1[[1]] * (raf1[[1]] <= 0.5) + (1 - raf1[[1]]) * (raf1[[1]] > 0.5)
    beta_vec <- do.call(pmax, lapply(beta1, abs)) # find largest in magnitude beta across traits
     names(beta_vec) <- names(beta1[[1]])
     names(maf_vec) <- names(raf1[[1]])
#    reftags <- cor.refdata2(corX, r2)
	reftags <- cor.refdata2_mod(corX, beta=beta_vec, MAF=maf_vec, r2 = r2)
    refGt <- as.matrix(reftags$refG)
    taglist <- reftags$taglist
	rm(reftags)
	
 	refG <- lqmm::make.positive.definite(refGt)
 	refG <- as.matrix(refG)   
    out <- list(SM = NULL, mbeta = NULL, Nlist = Nlist)
    out$SM <- vector("list", M)
    names(out$SM) <- ts
    out$mbeta <- vector("list", M)
    names(out$mbeta) <- ts
    dd <- vector("list", M)
    SSy <- vector("list", M)
    Sxy <- vector("list", M)
    for (j in 1:M) {
    	raf <- raf1[[j]]
    	BETA <- beta1[[j]][colnames(refG)]
    	maf <- raf * (raf <= 0.5) + (1 - raf) * (raf > 0.5)
        mafs.ref <- maf[colnames(refG)]
        covX <- cor2cov(corX, sd = sqrt(2 * raf * (1 - raf)))
        covX <- covX[snps, snps]  
		 topmods <- JAM.tries.maxcv(BETA=BETA,Vy=Vy[j],refG=refG,mafs.ref=mafs.ref,N=N[j],save.path=save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)
#        jam.results <- R2BGLiMS::JAM(marginal.betas = BETA, trait.variance = Vy[j], 
#            cor.ref = refG, mafs.ref = mafs.ref, model.space.priors = list(a = 1, 
#                b = length(BETA), Variables = names(BETA)), max.model.dim = maxcv, 
#            n = N[j], xtx.ridge.term = 0.01, save.path = save.path)
#        topmods <- R2BGLiMS::TopModels(jam.results, n.top.models = 1000)
#        if (is.null(ncol(topmods))) {
#            stop("A single model was selected with PP=1. This may mean no convergence because a causal variant is missing from the data and it has no tags in your data.")
#        }
        binout <- as.matrix(topmods[, -ncol(topmods)])
        colnames(binout) <- colnames(topmods)[-ncol(topmods)]
        snpmods <- apply(binout, 1, mod.fn)
        
         # for any models that have a pair of snps with abs(cor) > 0.8, remove one from the pair with largest mag beta until no pairs have high cor       
        snpmods <- prune_models_by_corr(snpmods=snpmods, corX=refG, beta=BETA, thr = 0.8)
  		nmod <- lengths(strsplit(snpmods, "%", fixed = TRUE))
        
 #       nmod <- apply(binout, 1, sum)
 #       PP <- topmods[, ncol(topmods)]
         snpPP <- data.frame(rank = seq_along(nmod), size = nmod, 
             str = snpmods, snps = snpmods, 
            stringsAsFactors = FALSE)
 #       snpPP <- snpPP[order(snpPP$PP, decreasing = TRUE), ]
        expmods <- rlist::list.stack(lapply(snpPP$snps, tagexpand.mod, 
            taglist = taglist))
        wh <- which(duplicated(expmods$snps))
        if (length(wh) > 0) {
            expmods <- expmods[-wh, ]
        }
        row.names(expmods) <- expmods$snps
        check <- sapply(strsplit(expmods[, 2], "%", fixed = TRUE), function(x) length(x) > 
            length(unique(x)))
        if (sum(check) > 0) 
            expmods <- expmods[-which(check), ]
        
        Vy_j <- Vy[j]
		N_j  <- N[j] 
		beta1_j <- beta1[[j]]   
		ybar_j <- ybar[j]

            
#        mbeta <- lapply(expmods[, 2], multibeta, beta1[[j]], 
#            covX, N = N[j], ybar = ybar[j], is.snpmat = FALSE, 
#            raf = raf)
#        names(mbeta) <- expmods[, 2]

		mbeta_for_mod <- function(mod) {
  			out <- multibeta(mod, beta1_j, covX, N = N_j, ybar = ybar_j, is.snpmat = FALSE, raf = raf)
  			if(mod != "1") modname <- paste(rownames(out)[-1], collapse = "%") # not including intercept in name
  			if(mod == "1") modname <- "1"
  			setNames(list(out), modname)
		}
		
		mbeta <- do.call(c, lapply(expmods[, 2], mbeta_for_mod))

        SSy[[j]] <- Vy_j * (N_j - 1) + N_j * ybar_j^2
        SSy_j <- SSy[[j]]
         Vx <- diag(covX)
        Mx <- 2 * raf
        
       Sxy[[j]] <- c(Sxy.hat(beta1 = beta1_j, Mx = Mx, N = N_j, 
            Vx = Vx, muY = ybar_j), `1` = ybar_j * N_j)
        names(Sxy[[j]])[length(Sxy[[j]])] <- "one"
        Sxy_j <- Sxy[[j]]
        
#        lABF <- sapply(expmods$snps, calcABF, mbeta, SSy = SSy[[j]], 
#            Sxy = Sxy[[j]], Vy = Vy[j], N = N[j])
#        names(lABF) <- expmods$snps

		abf_for_mod <- function(mod) {
  			calcABF(mod, mbeta, SSy_j, Sxy_j, Vy_j, N_j)
			}
		expmodnames <- names(mbeta)
		lABF <- vapply(expmodnames, abf_for_mod, numeric(1))
		names(lABF) <- expmodnames        

        wh <- which(expmodnames == "1")
        if (!length(wh)) {
            dd[[j]] <- data.frame(model = c("1", expmodnames), 
                tag = c(FALSE, expmods$tag), lBF = c(0, lABF), 
                stringsAsFactors = FALSE)
            l1 <- multibeta("1", beta1[[j]], covX, N = N[j], 
                ybar = ybar[j], is.snpmat = FALSE, raf = raf)
            mbeta <- rlist::list.append(mbeta, `1` = l1)
        }
        else {
            dd[[j]] <- data.frame(model = expmodnames, tag = expmods$tag, 
                lBF = lABF, stringsAsFactors = FALSE)
        }
        SM <- makesnpmod(dd[[j]], expected = 2, nsnps = nsnps)
        out$SM[[j]] <- SM
        out$mbeta[[j]] <- mbeta
        names(out$SM) <- ts
        names(out$mbeta) <- ts
        out$SM[[j]] <- best.models.cpp(out$SM[[j]], maxmod = 1000)[[1]]
        out$SM[[j]] <- PP2snpmod(out$SM[[j]])
    }
    out$Nlist <- Nlist
    out$nsnps <- nsnps
    out$Gmat = covX
    out$beta1.list = beta1
    out$raf = raf1[[which.max(N)]]
    return(out)
}



JAMmulti.tries <- function (gwas.list, corX, ybar, Vy, N, save.path,maxcv=2,maxcv_stop =20,maxcv_autocheck =TRUE,jam.nM.iter=1,extra.java.arguments=NULL) 
{
    tryCatch({
        JAMmulti(gwas.list=gwas.list, corX=corX, ybar=ybar, Vy=Vy, N=N, r2 = 0.99, 
            save.path=save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)
    }, error = function(e) {
        JAMmulti2(gwas.list=gwas.list, corX=corX, ybar=ybar, Vy=Vy, N=N, 
            r2 = 0.99, save.path=save.path,maxcv=maxcv,maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)
    })
}

#' @title Wrapper for flashfm Multi-Trait Fine-Mapping with JAM  - this is the dynamic number of max causal variant version
#' @param gwas.list List of M data.frame objects, where M is the number of traits (at most 6); gwas.list\[\[i\]\] is a data.frame for  trait i with 3 columns named: rsID, beta, EAF
#' @param corX SNP correlation matrix  
#' @param ybar trait mean vector; if not specified, the default is a vector of zeros
#' @param N Vector of length M; Nall\[i\] is the (effective) sample size for  trait i
#' @param save.path Path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1").
#' @param TOdds target odds of no sharing to sharing; default is 1
#' @param covY trait covariance matrix (for at most 5 traits and all traits should have a signal in the region, e.g. min p < 1E-6)
#' @param cpp cumulative posterior probability threshold for selecting top models; default 0.99
#' @param NCORES number of cores for parallel computing; recommend NCORES=M, but if on Windows, use NCORES=1
#' @param maxcv starting value for maximum number of causal variants; default 1
#' @param maxcv_stop maximum value to consider for maximum number of causal variants; maxcv_stop >= maxcv.
#' @param NCORES number of cores for parallel computing; recommend NCORES=A, but if on Windows, use NCORES=1
#' @param jam.nM.iter in millions, number of iterations to use in JAM; defailt 5 (5 million)
#' @param r2 r.squared threshold for to find tag SNPs to reduce model search space before expanding to all SNPs; 
#' @param extra.java.arguments A character string to be passed through to the java command line. E.g. to specify a
#' different temporary directory by passing "-Djava.io.tmpdir=/Temp".
#' @return list with 2 components: mpp.pp, a list with 4 components giving the SNP-level results (mpp.pp$PP,mpp.pp$MPP) and SNP group level results (mpp.pp$MPPg, mpp.pp$PPg); and snpGroups, 
#' a list with 2 components giving the SNP groups construced under single-trait (snpGroups\[\[1\]\]) and multi-trait fine-mapping (snpGroups\[\[2\]\])
#' @author Jenn Asimit
#' @import R2BGLiMS
#' @export
FLASHFMwithJAMd <- function (gwas.list, corX, ybar=NULL, N, save.path, TOdds = 1, covY, 
    cpp = 0.99, NCORES=1, maxcv=1, maxcv_stop = 20,jam.nM.iter=5, r2=0.8, extra.java.arguments=NULL) 
{
    maxcv_autocheck = TRUE
    M <- length(gwas.list)
    if(M>6 | M<2) stop("Need at least 2 and at most 6 traits.")
     if(is.null(ybar)) ybar <- rep(0,M)
    for(i in 1:M) gwas.list[[i]] <- as.data.frame(gwas.list[[i]])
    Vy <- diag(covY)
    corX <- as.matrix(corX)
    if(!dir.exists(save.path)) {
     message(c("Directory ",save.path," does not exist. Creating directory ",save.path))
     dir.create(save.path)
     }
    tmpdir <- paste0(save.path, "/tmp",sample(1:1000,1))   
    dir.create(tmpdir) 
    main.input <- JAMmulti2(gwas.list=gwas.list, corX=corX, ybar=ybar, Vy=Vy, N=N, 
            r2 = r2, save.path=save.path,maxcv=maxcv,maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)
    gc(verbose = FALSE)
    ss.stats <- summaryStats(Xmat = FALSE, ybar.all = ybar, main.input = main.input)
    if(M<6) {
     fm.multi <- flashfmU(main.input=main.input, TOdds = TOdds, covY=covY, ss.stats=ss.stats, 
        cpp = cpp, maxmod = NULL, fastapprox = FALSE, NCORES = NCORES)
    } 
    if(M==6) {
     fm.multi <- flashfmU(main.input=main.input, TOdds = TOdds, covY=covY, ss.stats=ss.stats, 
        cpp = cpp, maxmod = NULL, fastapprox = TRUE, NCORES = NCORES)
    }
    snpGroups <- makeSNPgroups2U(main.input=main.input, fm.multi=fm.multi, is.snpmat = FALSE, 
        min.mppi = 0.01, minsnpmppi = 0.01, r2.minmerge = 0.6)
    mpp.pp <- PPsummarise(fm.multi=fm.multi, snpGroups=snpGroups, minPP = 0.01)
    unlink(paste0(tmpdir,"/*"))
    return(list(mpp.pp = mpp.pp, snpGroups = snpGroups))
}



prune_models_by_corr <- function(snpmods, corX, beta, thr = 0.8) {
  pruned <- character(length(snpmods))

  # ensure beta is named; if some SNPs missing in beta, treat their beta as 0
  if (is.null(names(beta))) stop("beta must be a named vector")
  
  for (i in seq_along(snpmods)) {
    snps <- unlist(strsplit(snpmods[i], "%", fixed = TRUE))
    # keep only SNPs that exist in corX
    snps <- snps[snps %in% rownames(corX)]
    if (length(snps) <= 1L) {
      pruned[i] <- paste(snps, collapse = "%")
      next
    }

    # iterative pruning until all pairwise |r| < thr
    repeat {
      subcor <- abs(corX[snps, snps, drop = FALSE])
      diag(subcor) <- 0
      max_r <- max(subcor, na.rm = TRUE)
      if (is.na(max_r) || max_r < thr) break

      # take the first pair that attains the max correlation
      pair_idx <- which(subcor == max_r, arr.ind = TRUE)[1, , drop = TRUE]
      s1 <- snps[pair_idx[1]]
      s2 <- snps[pair_idx[2]]

      # get |beta| (missing names -> treated as 0)
      b1 <- if (!is.na(beta[s1])) abs(beta[s1]) else 0
      b2 <- if (!is.na(beta[s2])) abs(beta[s2]) else 0

      # drop the SNP with the smaller |beta|; if equal, drop s2
      drop_snp <- if (b1 < b2) s1 else s2
      snps <- setdiff(snps, drop_snp)

      if (length(snps) <= 1L) break
    }

    pruned[i] <- paste(snps, collapse = "%")
  }

  # return unique pruned models
  pruned_models <- unique(pruned)
  pruned_models <- pruned_models[nzchar(pruned_models)]
  
  return(pruned_models)
}


#' @title Calculate approximate Bayes' factor (ABF) 
#' @param mod joint SNP model with snps separated by \code{"\%"} e.g. \code{"snp1\%snp2"}
#' @param mbeta joint effect estimates for SNPs in model given by mod; output from multibeta
#' @param SSy sum(y.squared) for trait y
#' @param Sxy vector of sum(xy) for snp x, trait y
#' @param Vy trait variance
#' @param N sample size
#' @return ABF for model mod
#' @author Jenn Asimit
calcABFo <- function(mod,mbeta,SSy,Sxy,Vy,N) {
 if(mod=="1") { msnps <- "one"
 } else { msnps <- c("one", unlist(strsplit(mod, "%", fixed = TRUE))) }
 beta <- mbeta[[mod]]
 num <- SSy-sum(Sxy[msnps]*beta)
 den <- (N-1)*Vy
 k <- length(msnps)-1
 out <- -N*.5*log(num/den)-0.5*k*log(N)
 return(out)
}


calcABF <- compiler::cmpfun(calcABFo)

#' @title Using summary statistics, calculates joint effect estimates
#' @param mod joint SNP model with snps separated by \code{"\%"} e.g. \code{"snp1\%snp2"}
#' @param beta1 named vector single-SNP effect estimates; each effect estimate has name given by SNP id and ids must appear in Gmat columns
#' @param Gmat genotype matrix with SNP columns and individuals rows; could be reference matrix or from sample
#' @param N sample size for trait 
#' @param ybar trait mean
#' @param is.snpmat logical taking value TRUE when Gmat is a genotype matrix and FALSE when Gmat is a SNP covariance matrix
#' @param raf named vector of SNP reference allele frequencies (must match SNP coding used for effect estimates); only needed if Gmat is a covariance matrix
#' @param high_beta_thresh if absolute value of effect estimate < nigh_beta_thresh, the model is re-fit after checking for and removing variants to avoid high correlation
#' @return joint effect estimates for SNPs in model given by mod
#' @author Jenn Asimit
multibetao <- function(mod,beta1,Gmat,N,ybar,is.snpmat,raf=NULL,high_beta_thresh=5) {

 if(mod == "1") {
  mbeta <- ybar; names(mbeta) <- "one"
 } else {
  msnps <- unlist(strsplit(mod, "%", fixed = TRUE))
  if(all(msnps %in% names(beta1))) { 
   B <- beta1[msnps]
   if(is.snpmat) {  
     if(all(msnps %in% colnames(Gmat))){ mbeta <- JAM_PointEstimates_updated(B,X.ref=as.matrix(Gmat[,msnps]),n=N,ybar=ybar)
       } else { mbeta <- NA}
     } else { 
 	       colnames(Gmat) <- names(raf)
 	       rownames(Gmat) <- names(raf)
 	       if(all(msnps %in% colnames(Gmat))) { 
 	       
 	   #  extract covariance submatrix
          Xcov <- as.matrix(Gmat[msnps, msnps, drop = FALSE])
          # ensure exact symmetry
          Xcov <- (Xcov + t(Xcov)) / 2

          # check eigenvalues for positive-definiteness 
          eigvals <- eigen(Xcov, symmetric = TRUE, only.values = TRUE)$values
          tol <- 1e-10  # tolerance for tiny negative eigenvalues

          if (min(eigvals, na.rm = TRUE) < -tol) {
            # matrix is numerically non-PSD — attempt to fix, preferring lqmm if available
            fixed <- FALSE
 	       
 	       try({
                Xcov <- lqmm::make.positive.definite(Xcov)
                Xcov <- (Xcov + t(Xcov)) / 2
                fixed <- TRUE
              }, silent = TRUE)
              
           # as a final guard, add a very small ridge to ensure numerical PD
            Xcov <- Xcov + diag(1e-12, nrow(Xcov))
          } else {
            # if smallest eigenvalue is only slightly negative within tolerance,
            # add a tiny ridge to be safe
            Xcov <- Xcov + diag(1e-12, nrow(Xcov))
          }   
 	       
 	       mbeta <- JAM_PointEstimates_Xcov(B,Xcov=Xcov,raf=raf[msnps],n=N,ybar=ybar)  
 	       
 	       # check for unusually high effects
 	        corX <- cov2cor(Xcov)
          if(!is.na(mbeta)[1] && !is.null(corX) && max(abs(mbeta)) > high_beta_thresh && length(msnps) > 1) {
                      
            pruned_mod <- prune_models_by_corr(snpmods = mod, corX = corX, beta = B, thr = 0.7)
            msnps <- unlist(strsplit(pruned_mod, "%", fixed = TRUE))
            B <- beta1[msnps]
            Xcov <- as.matrix(Gmat[msnps, msnps, drop = FALSE])
            Xcov <- (Xcov + t(Xcov)) / 2
            mbeta <- JAM_PointEstimates_Xcov(B, Xcov = Xcov, raf = raf[msnps], n = N, ybar = ybar)        	       
 	          } 
 	       }else { mbeta <- NA}
 	       }
   } else {mbeta <- NA }
 }
 
 return(mbeta)
} 

multibeta <- compiler::cmpfun(multibetao)


########
#' @title Expanded version of JAM (a single-trait fine-mapping approach) that first runs on thinned SNPs and then expands models on tag SNPs
#' This version starts at a low upper bound for max causal variants and decides on max upper bound based on data
#' @param gwas data.frame with 3 columns named: rsID, beta, EAF
#' @param corX genotype correlation matrix (reference or from sample) 
#' @param ybar trait mean; if trait is transformed to be standard Normal, set ybar = 0, which is default
#' @param Vy trait variance; if trait is transformed to be standard Normal, set Vy = 1, which is default
#' @param N sample size for trait; recommended to give effective sample sizes using GWAS summary statistics in Neff function
#' @param cred probability for credible set; default is 0.99
#' @param save.path path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1"). 
#' @param maxcv starting value for maximum number of causal variants
#' @param maxcv_stop maximum value to consider for maximum number of causal variants; maxcv_stop >= maxcv.
#' @param jam.nM.iter in millions, number of iterations to use in JAM; defailt 1 (1 million)
#' @param r2 r.squared threshold for to find tag SNPs to reduce model search space before expanding to all SNPs; 
#' @param extra.java.arguments A character string to be passed through to the java command line. E.g. to specify a
#' different temporary directory by passing "-Djava.io.tmpdir=/Temp".
#' @return List of credible set variants with their MPP, MPP for all variants, PP for all models 
#' @import R2BGLiMS
#' @export
#' @author Feng Zhou
JAMdynamic <- function(gwas, corX, ybar=0, Vy=1, N, cred=.99, save.path, maxcv=1, maxcv_stop = 20,jam.nM.iter=5, r2=0.80, extra.java.arguments=NULL){
 
  maxcv_autocheck = TRUE
  if(!dir.exists(save.path)) {
     message(c("Directory ",save.path," does not exist. Creating directory ",save.path))
     dir.create(save.path)
     }
 gwas <- as.data.frame(gwas)
 beta1 <- gwas[,"beta"] 
 raf <- gwas[,"EAF"]
 names(beta1) <- names(raf) <- gwas[,"rsID"]
 corX <- as.matrix(corX)
 ksnp <- intersect(colnames(corX),names(raf))
 beta1 <- beta1[ksnp]
 raf <- raf[ksnp]
 corX <- corX[ksnp,ksnp]
 
 fm <- JAMcor.tries.maxcv(beta1=beta1, corX=corX, raf=raf, ybar=ybar, Vy=Vy, N=N, save.path=save.path, maxcv=maxcv, maxcv_stop = maxcv_stop, 
 maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter, r2=r2, extra.java.arguments=extra.java.arguments)
 ind <- order(fm$PP, decreasing = TRUE)
 pp <- fm$PP[ind]
 mod <- rownames(fm)[ind]
 cpp <- cumsum(pp)
 wh <- which(cpp <= cred)
 if (!length(wh)) wh <- 1
 if(cpp[max(wh)] < cred) wh <- c(wh, max(wh) + 1)
 mods <- mod[wh]
 cs <- unique(unlist(strsplit(mods, "%")))
 cs <- cs[!is.na(cs)]
 mpp <- MPPcalc(fm)
 CS <- mpp[cs] 
 CS <- CS[order(CS,decreasing=T)]
 return(list(CS=CS,MPP=mpp,PP=fm))
}  



JAMcor.tries.maxcv <- function(beta1, corX, raf, ybar, Vy, N, save.path, maxcv=1, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1,r2=0.99,extra.java.arguments=NULL){
  
   
  JAM_output <- NULL

  while (maxcv_autocheck == TRUE){
    
    tryCatch({
      
      JAM_output <- JAMexpandedCor2(beta1=beta1, corX=corX, raf=raf, ybar=ybar, Vy=Vy, N=N, 
            r2 = r2, save.path=save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)
        
      maxcv_autocheck = FALSE
      print(paste0('Completed the process when ..., maxcv == ', maxcv))
      
    }, error = function(e){
      
      print(paste0('Keep trying maxcv by adding 1 ..., maxcv == ', maxcv+1))
      
    })
    
    maxcv = maxcv+1
    
    if (maxcv > maxcv_stop & is.null(JAM_output)){
      stop("The maxcv_stop reached, no result! Try increasing maxcv_stop.")
      maxcv_autocheck = FALSE
    }
  
  }
  return(JAM_output)
}


JAMcor.tries <- function (beta1, corX, raf, ybar, Vy, N, save.path,maxcv=2,maxcv_stop =20,maxcv_autocheck =TRUE,jam.nM.iter=1,extra.java.arguments=NULL) 
{
    tryCatch({
        JAMexpandedCor(beta1=beta1, corX=corX, raf=raf, ybar=ybar, Vy=Vy, N=N, r2 = 0.99, 
            save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)
    }, error = function(e) {
        JAMexpandedCor2(beta1=beta1, corX=corX, raf=raf, ybar=ybar, Vy=Vy, N=N, 
            r2 = 0.99, save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)
    })
}



JAMexpandedCor <- function (beta1, corX, raf, ybar, Vy, N, r2 = 0.99, save.path,maxcv=2, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1,extra.java.arguments=NULL) 
{
    covX <- cor2cov(corX, sd = sqrt(2 * raf * (1 - raf)))
   
    snps <- names(beta1)
    snps <- intersect(snps, names(raf))
    beta1 <- beta1[snps]
    if (is.null(colnames(corX))) {
        colnames(corX) <- names(raf)
        rownames(corX) <- names(raf)
    }
    corX <- corX[snps, snps]
    covX <- covX[snps, snps]
    maf <- raf * (raf <= 0.5) + (1 - raf) * (raf > 0.5)
    nsnps <- ncol(corX)
    refG <- cor.refdata.fn(corX, r2)
    corX2 <- corX^2
    taglist <- tagSNP(corX2, threshold = r2)

        BETA <- beta1[colnames(refG)]
        mafs.ref <- maf[colnames(refG)]
		 topmods <- JAM.tries.maxcv(BETA=BETA,Vy=Vy,refG=refG,mafs.ref=mafs.ref,N=N,save.path=save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)
#          jam.results <- R2BGLiMS::JAM(marginal.betas = BETA, trait.variance = Vy, 
#            cor.ref = refG, mafs.ref = mafs.ref, model.space.priors = list(a = 1, 
#                b = length(BETA), Variables = names(BETA)), max.model.dim = maxcv, 
#            n = N, xtx.ridge.term = 0.01, save.path = save.path,n.mil.iter=5)    
#        topmods <- R2BGLiMS::TopModels(jam.results, n.top.models = 1000)
#        if (is.null(ncol(topmods))) {
#            stop("A single model was selected with PP=1. This may mean no convergence because a causal variant is missing from the data and it has no tags in your data.")
#        }
        binout <- as.matrix(topmods[, -ncol(topmods)])
        colnames(binout) <- colnames(topmods)[-ncol(topmods)]
        snpmods <- apply(binout, 1, mod.fn)
        nmod <- apply(binout, 1, sum)
        PP <- topmods[, ncol(topmods)]
        snpPP <- data.frame(rank = 1:length(nmod), size = nmod, 
            logPP = log(PP), PP = PP, str = snpmods, snps = snpmods, 
            stringsAsFactors = FALSE)
        snpPP <- snpPP[order(snpPP$PP, decreasing = TRUE), ]
        expmods <- rlist::list.stack(lapply(snpPP$snps, tagexpand.mod, 
            taglist = taglist))
        wh <- which(duplicated(expmods$snps))
        if (length(wh) > 0) {
            expmods <- expmods[-wh, ]
        }
        row.names(expmods) <- expmods$snps
        check <- sapply(strsplit(expmods[, 2], "%"), function(x) length(x) > 
            length(unique(x)))
        if (sum(check) > 0) 
            expmods <- expmods[-which(check), ]
        mbeta <- lapply(expmods[, 2], multibeta, beta1, 
            covX, N = N, ybar = ybar, is.snpmat = FALSE, 
            raf = raf)
        names(mbeta) <- expmods[, 2]
        SSy <- Vy * (N - 1) + N * ybar^2
        Vx <- diag(covX)
        Mx <- 2 * raf
        Sxy <- c(Sxy.hat(beta1 = beta1, Mx = Mx, N = N, 
            Vx = Vx, muY = ybar), `1` = ybar * N)
        names(Sxy)[length(Sxy)] <- "one"
        lABF <- sapply(expmods$snps, calcABF, mbeta, SSy = SSy, 
            Sxy = Sxy, Vy = Vy, N = N)
        names(lABF) <- expmods$snps
        wh <- which(expmods$snps == "1")
        if (!length(wh)) {
            dd <- data.frame(model = c("1", expmods$snps), 
                tag = c(FALSE, expmods$tag), lBF = c(0, lABF), 
                stringsAsFactors = FALSE)
            l1 <- multibeta("1", beta1, covX, N = N, 
                ybar = ybar, is.snpmat = FALSE, raf = raf)
            mbeta <- rlist::list.append(mbeta, `1` = l1)
        }
        else {
            dd <- data.frame(model = expmods$snps, tag = expmods$tag, 
                lBF = lABF, stringsAsFactors = FALSE)
        }
        SM <- makesnpmod(dd, expected = 2, nsnps = nsnps)
        SM <- best.models.cpp(SM, cpp.thr=.99)[[1]]
        SM <- PP2snpmod(SM)
        tmp <- SM@models[,c("str","PP")]
		rownames(tmp) <- tmp$str
		stfm <- data.frame(PP=tmp[,-1],row.names=tmp$str)	
    return(stfm)
}


JAMexpandedCor2 <- function (beta1, corX, raf, ybar, Vy, N, r2 = 0.99, save.path,maxcv=2, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1,extra.java.arguments=NULL)
{
    covX <- cor2cov(corX, sd = sqrt(2 * raf * (1 - raf)))
   
    snps <- names(beta1)
    snps <- intersect(snps, names(raf))
    beta1 <- beta1[snps]
    if (is.null(colnames(corX))) {
        colnames(corX) <- names(raf)
        rownames(corX) <- names(raf)
    }
    corX <- corX[snps, snps]
    covX <- covX[snps, snps]
    maf <- raf * (raf <= 0.5) + (1 - raf) * (raf > 0.5)
    nsnps <- ncol(corX)
    reftags <- cor.refdata2_mod(corX, beta=beta1[snps], MAF=maf[snps], r2 = r2)
    refGt <- reftags$refG
    taglist <- reftags$taglist
	rm(reftags)
	refG <- lqmm::make.positive.definite(refGt)
		
        BETA <- beta1[colnames(refG)]
        mafs.ref <- maf[colnames(refG)]
		 topmods <- JAM.tries.maxcv(BETA=BETA,Vy=Vy,refG=refG,mafs.ref=mafs.ref,N=N,save.path=save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter,extra.java.arguments=extra.java.arguments)
#         jam.results <- R2BGLiMS::JAM(marginal.betas = BETA, trait.variance = Vy, 
#            cor.ref = refG, mafs.ref = mafs.ref, model.space.priors = list(a = 1, 
#                b = length(BETA), Variables = names(BETA)), max.model.dim = maxcv, 
#            n = N, xtx.ridge.term = 0.01, save.path = save.path,n.mil.iter=5)    
#        topmods <- R2BGLiMS::TopModels(jam.results, n.top.models = 1000)
#        if (is.null(ncol(topmods))) {
#            stop("A single model was selected with PP=1. This may mean no convergence because a causal variant is missing from the data and it has no tags in your data.")
#        }
        binout <- as.matrix(topmods[, -ncol(topmods)])
        colnames(binout) <- colnames(topmods)[-ncol(topmods)]
        snpmods <- apply(binout, 1, mod.fn)
        
         # for any models that have a pair of snps with abs(cor) > 0.8, remove one from the pair with largest mag beta until no pairs have high cor       
        snpmods <- prune_models_by_corr(snpmods=snpmods, corX=refG, beta=BETA, thr = 0.8)
  		nmod <- lengths(strsplit(snpmods, "%", fixed = TRUE))
          
 #       nmod <- apply(binout, 1, sum)
 #       PP <- topmods[, ncol(topmods)]
         snpPP <- data.frame(rank = seq_along(nmod), size = nmod, 
             str = snpmods, snps = snpmods, 
            stringsAsFactors = FALSE)
#        snpPP <- snpPP[order(snpPP$PP, decreasing = TRUE), ]
        expmods <- rlist::list.stack(lapply(snpPP$snps, tagexpand.mod, 
            taglist = taglist))
        wh <- which(duplicated(expmods$snps))
        if (length(wh) > 0) {
            expmods <- expmods[-wh, ]
        }
        row.names(expmods) <- expmods$snps
        check <- sapply(strsplit(expmods[, 2], "%", fixed = TRUE), function(x) length(x) > 
            length(unique(x)))
        if (sum(check) > 0) 
            expmods <- expmods[-which(check), ]
            
        Vy_j <- Vy
		N_j  <- N 
		beta1_j <- beta1   
		ybar_j <- ybar
    
            
#        mbeta <- lapply(expmods[, 2], multibeta, beta1, 
#            covX, N = N, ybar = ybar, is.snpmat = FALSE, 
#            raf = raf)
#        names(mbeta) <- expmods[, 2]

		mbeta_for_mod <- function(mod) {
  			out <- multibeta(mod, beta1_j, covX, N = N_j, ybar = ybar_j, is.snpmat = FALSE, raf = raf)
  			if(mod != "1") modname <- paste(rownames(out)[-1], collapse = "%") # not including intercept in name
  			if(mod == "1") modname <- "1"
  			setNames(list(out), modname)
		}
		
		mbeta <- do.call(c, lapply(expmods[, 2], mbeta_for_mod))


        SSy_j <- Vy_j * (N_j - 1) + N_j * ybar_j^2
        
        Vx <- diag(covX)
        Mx <- 2 * raf
        Sxy_j <- c(Sxy.hat(beta1 = beta1_j, Mx = Mx, N = N_j, 
            Vx = Vx, muY = ybar_j), `1` = ybar_j * N_j)
        names(Sxy_j)[length(Sxy_j)] <- "one"
       
 
		abf_for_mod <- function(mod) {
  			calcABF(mod, mbeta, SSy_j, Sxy_j, Vy_j, N_j)
			}
		expmodnames <- names(mbeta)
		lABF <- vapply(expmodnames, abf_for_mod, numeric(1))
		 names(lABF) <- expmodnames
		        
        wh <- which(expmodnames == "1")
        if (!length(wh)) {
            dd <- data.frame(model = c("1", expmodnames), 
                tag = c(FALSE, expmods$tag), lBF = c(0, lABF), 
                stringsAsFactors = FALSE)
            l1 <- multibeta("1", beta1, covX, N = N, 
                ybar = ybar, is.snpmat = FALSE, raf = raf)
            mbeta <- rlist::list.append(mbeta, `1` = l1)
        }
        else {
            dd <- data.frame(model = expmodnames, tag = expmods$tag, 
                lBF = lABF, stringsAsFactors = FALSE)
        }
        SM <- makesnpmod(dd, expected = 2, nsnps = nsnps)
        SM <- best.models.cpp(SM, cpp.thr=.99)[[1]]
        SM <- PP2snpmod(SM)
        tmp <- SM@models[,c("str","PP")]
		rownames(tmp) <- tmp$str
		stfm <- data.frame(PP=tmp[,-1],row.names=tmp$str)	
    return(stfm)
}


#' @title Construct a credible set for each trait and under each of single and multi-trait fine-mapping, and provide SNP PP details
#' @param mpp.pp object created in flashfm output, e.g. fm$mpp.pp if fm is output from FLASHFMwithJAM or FLASHFMwithFINEMAP
#' @param cred probability for credible set; default is 0.99
#' @return list with 3 components: list of single-trait fine-mapping credible sets for each trait with the SNP PP for each variant, list of multi-trait fine-mapping credible sets for each trait with the SNP PP for each variant, credible set probability
#' @export
#' @author Jenn Asimit
allcredsetsPP <- function (mpp.pp, cred = 0.99) 
{
    M <- length(mpp.pp$PP)
    csfm <- csflfm <- vector("list", M)
    for (i in 1:M) {
        cs <- credsetU(mpp.pp$PP[[i]][, 1], cred)
        mpp <- mpp.pp$MPP[[i]][cs, 1]
        csfm[[i]] <- data.frame(SNP=cs, MPP=mpp)
        csfm[[i]] <-  csfm[[i]][order(csfm[[i]]$MPP,decreasing=TRUE),]
        cs <- credsetU(mpp.pp$PP[[i]][, 2], cred)
        mpp <- mpp.pp$MPP[[i]][cs, 2]
        csflfm[[i]] <- data.frame(SNP=cs, MPP=mpp)
        csflfm[[i]] <-  csflfm[[i]][order(csflfm[[i]]$MPP,decreasing=TRUE),]
    }
    return(list(fm = csfm, flashfm = csflfm, cred = cred))
}

JAM_PointEstimates_updated <- function (marginal.betas = NULL, X.ref = NULL, n = NULL, ybar) 
{
    n.ref <- nrow(X.ref)
    if (is.null(n)) {
        n <- n.ref
    }
    z <- rep(NA, length(marginal.betas))
    names(z) <- names(marginal.betas)
    for (v in 1:length(marginal.betas)) {
        y.pred <- X.ref[, v] * marginal.betas[v]
        y.pred.centered <- y.pred - mean(y.pred)
        z[v] <- X.ref[, v] %*% y.pred.centered
    }
    z <- z * n/n.ref
    xbar <- apply(X.ref, 2, mean)
    for (v in 1:ncol(X.ref)) {
        X.ref[, v] <- X.ref[, v] - xbar[v]
    }
    xtx <- (t(X.ref) %*% X.ref) * n/n.ref
    multivariate.beta.hat <- solve(xtx) %*% z
    b0 <- ybar - sum(xbar * multivariate.beta.hat)
    out <- rbind(b0, multivariate.beta.hat)
    rownames(out) <- c("one", names(marginal.betas))
    return(out)
}


JAM_PointEstimates_Xcov <- function (marginal.betas = NULL, Xcov = NULL, raf, n, ybar) 
{
    xbar <- 2 * raf
    z <- rep(NA, length(marginal.betas))
    names(z) <- names(marginal.betas)
    for (v in 1:length(marginal.betas)) {
        z[v] <- n * (marginal.betas[v] * Xcov[v, v])
    }
    xtx <- n * (Xcov)
    multivariate.beta.hat <- as.vector(solve(xtx) %*% z)
    b0 <- ybar - sum(as.vector(xbar) * multivariate.beta.hat)
    out <- matrix(c(b0, multivariate.beta.hat), ncol = 1)
    rownames(out) <- c("one", names(marginal.betas))
    return(out)
}

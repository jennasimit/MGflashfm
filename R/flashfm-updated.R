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




JAMmaxcv <- function(BETA,Vy,refG,mafs.ref,N,save.path,maxcv=2,jam.nM.iter=1) {
 jam.results <- R2BGLiMS::JAM(marginal.betas = BETA, trait.variance = Vy, 
            cor.ref = refG, mafs.ref = mafs.ref, model.space.priors = list(a = 1, 
                b = length(BETA), Variables = names(BETA)), max.model.dim = maxcv, 
            n = N, xtx.ridge.term = 0.01, save.path = save.path, n.mil.iter=jam.nM.iter) 
 topmods <- R2BGLiMS::TopModels(jam.results, n.top.models = 1000)
        if (is.null(ncol(topmods))  ) {
            stop("A single model was selected with PP=1. This may mean no convergence because a causal variant  is missing from the data and it has no tags in your data.")
        }
 return(topmods)
 }       

JAM.tries.maxcv <- function(BETA,Vy,refG,mafs.ref,N,save.path,maxcv=2, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1){
  
  while (maxcv_autocheck == TRUE){
    
    tryCatch({
      
      JAM_output <-  JAMmaxcv(BETA,Vy,refG,mafs.ref,N,save.path,maxcv=maxcv,jam.nM.iter=jam.nM.iter) 
        
      maxcv_autocheck = FALSE
      print(paste0('Completed the process when ..., maxcv == ', maxcv))
      
    }, error = function(e){
      
      print(paste0('Keep trying maxcv by adding 1 ..., maxcv == ', maxcv+1))
      
    })
    
    maxcv = maxcv+1
    
    if (maxcv == maxcv_stop){
      print("The maxcv_stop reached, no result!")
      maxcv_autocheck = FALSE
    }
  
  }
  return(JAM_output)
}

            
JAMmulti <- function (gwas.list, corX, ybar, Vy, N, r2 = 0.99, save.path, maxcv=2, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1) 
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
		 topmods <- JAM.tries.maxcv(BETA,Vy[j],refG,mafs.ref,N[j],save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter)
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


JAMmulti2 <- function (gwas.list, corX, ybar, Vy, N, r2 = 0.99, save.path, maxcv=2, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1) 
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
    reftags <- cor.refdata2(corX, r2)
    refGt <- reftags$refG
    taglist <- reftags$taglist

 	refG <- lqmm::make.positive.definite(refGt)
 	   
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
		 topmods <- JAM.tries.maxcv(BETA,Vy[j],refG,mafs.ref,N[j],save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter)
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



JAMmulti.tries <- function (gwas.list, corX, ybar, Vy, N, save.path,maxcv=2,maxcv_stop =20,maxcv_autocheck =TRUE,jam.nM.iter=1) 
{
    tryCatch({
        JAMmulti(gwas.list, corX, ybar, Vy, N, r2 = 0.99, 
            save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter)
    }, error = function(e) {
        JAMmulti2(gwas.list, corX, ybar, Vy, N, 
            r2 = 0.99, save.path,maxcv=maxcv,maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter)
    })
}

#' @title Wrapper for flashfm Multi-Trait Fine-Mapping with JAM  - this is the dynamic number of max causal variant version
#' @param gwas.list List of M data.frame objects, where M is the number of traits (at most 6); gwas.list\[\[i\]\] is a data.frame for  trait i with 3 columns named: rsID, beta, EAF
#' @param corX SNP correlation matrix  
#' @param ybar trait mean; if trait is transformed to be standard Normal, set ybar = 0, which is default
#' @param N Vector of length M; Nall\[i\] is the (effective) sample size for  trait i
#' @param save.path Path to save JAM output files; tmp files and could delete these later e.g. save.path=paste0(DIRout,"/tmpJAM/region1").
#' @param TOdds target odds of no sharing to sharing; default is 1
#' @param covY trait covariance matrix (for at most 5 traits and all traits should have a signal in the region, e.g. min p < 1E-6)
#' @param cpp cumulative posterior probability threshold for selecting top models; default 0.99
#' @param NCORES number of cores for parallel computing; recommend NCORES=M, but if on Windows, use NCORES=1
#' @param maxcv starting value for maximum number of causal variants; default 1
#' @param maxcv_stop maximum value to consider for maximum number of causal variants
#' @param maxcv_autocheck Logical for whether to incrementally increase maxcv (TRUE); if do not want to increment maxcv set this to FALSE
#' @param NCORES number of cores for parallel computing; recommend NCORES=A, but if on Windows, use NCORES=1
#' @param jam.nM.iter in millions, number of iterations to use in JAM; defailt 1 (1 million)
#' @return list with 2 components: mpp.pp, a list with 4 components giving the SNP-level results (mpp.pp$PP,mpp.pp$MPP) and SNP group level results (mpp.pp$MPPg, mpp.pp$PPg); and snpGroups, 
#' a list with 2 components giving the SNP groups construced under single-trait (snpGroups\[\[1\]\]) and multi-trait fine-mapping (snpGroups\[\[2\]\])
#' @author Jenn Asimit
#' @export
FLASHFMwithJAMd <- function (gwas.list, corX, ybar, N, save.path, TOdds = 1, covY, 
    cpp = 0.99, NCORES, maxcv=1, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1) 
{
    M <- length(ybar)
    if(M>6 | M<2) stop("Need at least 2 and at most 6 traits.")
    Vy <- diag(covY)
    corX <- as.matrix(corX)
    if(!dir.exists(save.path)) {
     message(c("Directory ",save.path," does not exist. Creating directory ",save.path))
     dir.create(save.path)
     }
    tmpdir <- paste0(save.path, "/tmp",sample(1:1000,1))   
    dir.create(tmpdir) 
    main.input <- JAMmulti2(gwas.list, corX, ybar, Vy, N, 
            r2 = 0.99, save.path,maxcv=maxcv,maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter)
    gc(verbose = FALSE)
    ss.stats <- summaryStats(Xmat = FALSE, ybar.all = ybar, main.input = main.input)
    if(M<6) {
     fm.multi <- flashfmU(main.input, TOdds = TOdds, covY, ss.stats, 
        cpp = cpp, maxmod = NULL, fastapprox = FALSE, NCORES = NCORES)
    } 
    if(M==6) {
     fm.multi <- flashfmU(main.input, TOdds = TOdds, covY, ss.stats, 
        cpp = cpp, maxmod = NULL, fastapprox = TRUE, NCORES = NCORES)
    }
    snpGroups <- makeSNPgroups2(main.input, fm.multi, is.snpmat = FALSE, 
        min.mppi = 0.01, minsnpmppi = 0.01, r2.minmerge = 0.6)
    mpp.pp <- PPsummarise(fm.multi, snpGroups, minPP = 0.01)
    unlink(paste0(tmpdir,"/*"))
    return(list(mpp.pp = mpp.pp, snpGroups = snpGroups))
}






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
#' @param maxcv_stop maximum value to consider for maximum number of causal variants
#' @param maxcv_autocheck Logical for whether to incrementally increase maxcv (TRUE); if do not want to increment maxcv set this to FALSE
#' @param jam.nM.iter in millions, number of iterations to use in JAM; defailt 1 (1 million)
#' @return List of credible set variants with their MPP, MPP for all variants, PP for all models 
#' @import R2BGLiMS
#' @export
#' @author Feng Zhou
JAMdynamic <- function(gwas, corX, ybar=0, Vy=1, N, cred=.99, save.path, maxcv=1, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1){
 
  if(!dir.exists(save.path)) {
     message(c("Directory ",save.path," does not exist. Creating directory ",save.path))
     dir.create(save.path)
     }
 
 beta1 <- gwas[,"beta"] 
 raf <- gwas[,"EAF"]
 names(beta1) <- names(raf) <- gwas[,"rsID"]
 corX <- as.matrix(corX)
 ksnp <- intersect(colnames(corX),names(raf))
 beta1 <- beta1[ksnp]
 raf <- raf[ksnp]
 corX <- corX[ksnp,ksnp]
 
 fm <- JAMcor.tries.maxcv(beta1, corX, raf, ybar, Vy, N, save.path, maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter)
 ind <- order(fm$PP, decreasing = TRUE)
 pp <- fm$PP[ind]
 mod <- rownames(fm)[ind]
 cpp <- cumsum(pp)
 wh <- which(cpp <= cred)
 if (!length(wh)) wh <- 1
 wh <- c(wh, max(wh) + 1)
 mods <- mod[wh]
 cs <- unique(unlist(strsplit(mods, "%")))
 cs <- cs[!is.na(cs)]
 mpp <- MPPcalc(fm)
 CS <- mpp[cs] 
 CS <- CS[order(CS,decreasing=T)]
 return(list(CS=CS,MPP=mpp,PP=fm))
}  



JAMcor.tries.maxcv <- function(beta1, corX, raf, ybar, Vy, N, save.path, maxcv=1, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1){
  
  while (maxcv_autocheck == TRUE){
    
    tryCatch({
      
      JAM_output <- JAMexpandedCor2(beta1, corX, raf, ybar, Vy, N, 
            r2 = 0.99, save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter)
        
      maxcv_autocheck = FALSE
      print(paste0('Completed the process when ..., maxcv == ', maxcv))
      
    }, error = function(e){
      
      print(paste0('Keep trying maxcv by adding 1 ..., maxcv == ', maxcv+1))
      
    })
    
    maxcv = maxcv+1
    
    if (maxcv == maxcv_stop){
      print("The maxcv_stop reached, no result!")
      maxcv_autocheck = FALSE
    }
  
  }
  return(JAM_output)
}


JAMcor.tries <- function (beta1, corX, raf, ybar, Vy, N, save.path,maxcv=2,maxcv_stop =20,maxcv_autocheck =TRUE,jam.nM.iter=1) 
{
    tryCatch({
        JAMexpandedCor(beta1, corX, raf, ybar, Vy, N, r2 = 0.99, 
            save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter)
    }, error = function(e) {
        JAMexpandedCor2(beta1, corX, raf, ybar, Vy, N, 
            r2 = 0.99, save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter)
    })
}



JAMexpandedCor <- function (beta1, corX, raf, ybar, Vy, N, r2 = 0.99, save.path,maxcv=2, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1) 
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
        topmods <- JAM.tries.maxcv(BETA,Vy,refG,mafs.ref,N,save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter)
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


JAMexpandedCor2 <- function (beta1, corX, raf, ybar, Vy, N, r2 = 0.99, save.path,maxcv=2, maxcv_stop = 20, maxcv_autocheck = TRUE,jam.nM.iter=1)
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
    reftags <- cor.refdata2(corX, r2)
    refGt <- reftags$refG
    taglist <- reftags$taglist

		refG <- lqmm::make.positive.definite(refGt)
		
        BETA <- beta1[colnames(refG)]
        mafs.ref <- maf[colnames(refG)]
		 topmods <- JAM.tries.maxcv(BETA,Vy,refG,mafs.ref,N,save.path,maxcv=maxcv, maxcv_stop = maxcv_stop, maxcv_autocheck = maxcv_autocheck,jam.nM.iter=jam.nM.iter)
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


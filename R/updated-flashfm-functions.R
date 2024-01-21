cor.refdata2 <- function (corX, r2 = 0.99) {
    gmat2t <- tagSNP(corX, threshold = sqrt(r2))
    gt <- unlist(gmat2t)
    tg = gt[names(gt) == "tagsnp"]
    refG <- corX[tg, tg]
    return(list(refG = refG, taglist = gmat2t))
}

# edited version of tagSNP in hscovar package
# here all snps in the bin have abs(r)<threshold
tagSNP <- function (mat, threshold = 0.8) 
{
    if (max(abs(mat), na.rm = T) > 1 + 1e-06) 
        stop("Correlation or R-squared matrix is expected.")
    mat <- apply(mat, 1, function(x) {
        x[!is.finite(x)] <- 0
        return(x)
    })
    diag(mat) <- 1
    if ((threshold >= 1) | (threshold <= 0)) 
        stop("Threshold out of range")
    p <- nrow(mat)
    bin <- list()
    counter <- 0
    snpset <- 1:p
    coln <- colnames(mat)
    colnames(mat) <- NULL
    snps <- NULL
    repeat {
        counter <- counter + 1
        ls <- lapply(snpset, function(i) {
            id <- abs(mat[i, snpset]) > threshold
            snpset[id]
        })
        m <- which.max(lapply(ls, length))
        if (length(ls[[m]]) == 1) {
            ts <- ls[[m]]
            candidate <- 1
        } else {
            candidate <- apply(mat[ls[[m]], ls[[m]]], 1, function(x) {
                all(abs(x) > threshold)
            })
            ts <- ls[[m]][candidate][ceiling(sum(candidate)/2)]
        }
        snps <- c(snps,ls[[m]][candidate])
        bin[[counter]] <- list(snps = coln[ls[[m]]][candidate], tagsnp = coln[ts])
        snpset <- setdiff(snpset, snps)
        if ((!length(snpset) > 0) || (counter == p)) 
            break
    }
    z <- length(unlist(rlist::list.select(bin, tp = snps)))
    message(paste(z, "SNPs have been grouped into", counter, 
        "bins"))
    return(bin)
}




best.models.cpp <- function (d, cpp.thr = 0.99, maxmod = NULL) {
	old.cpp = cpp.thr
    if (is.null(maxmod)) {
        d2 <- d
        d2@models <- d@models[order(d@models$PP, decreasing = TRUE), 
            ]
        cpp <- cumsum(d2@models$PP)
        wh <- which(cpp <= cpp.thr)
        if (!length(wh)) 
            wh <- 1
        wh <- c(wh, max(wh) + 1)
        d2@models <- d2@models[wh, ]
        d2@models$PP <- d2@models$PP/sum(d2@models$PP)
        d <- d2
        old.cpp <- cpp.thr
    } else {
    	if (maxmod < nrow(d@models)) {
        d@models <- d@models[order(d@models$PP, decreasing = TRUE), 
            ][1:min(maxmod, nrow(d@models)), ]
        message("Before adjustment, CPP of top", maxmod, "models is", 
            sum(d@models$PP))
        old.cpp <- sum(d@models$PP)
        d@models$PP <- d@models$PP/sum(d@models$PP)
        }
    }
    out <- cbind(d@models, snps = unlist(lapply(strsplit(d@models$str, 
        "%"), makestr)))
    return(list(models = out, old.cpp = old.cpp))
}


flashfmU <- function (main.input, TOdds, covY, ss.stats, cpp = 0.99, maxmod = NULL, 
    fastapprox = FALSE, NCORES) 
{

# check if covY is positive semi-definite and if not, then make it
	minev <- min(eigen(covY)$values)
	A <- covY
	if(minev < 0) {
	 diag(A) <- diag(A)+ abs(minev)+10^(-10)
	 covY <- A
	}

    Nlist <- main.input$Nlist
    Nqq <- as.matrix(Nlist$Nqq)
    Nq3 <- as.vector(Nlist$Nq3)
    Nq4 <- as.vector(Nlist$Nq4)
    N <- Nlist$N
    nsnps <- main.input$nsnps
    mbeta <- main.input$mbeta
    SM <- main.input$SM
    ybar <- ss.stats$ybar
    Sxy <- ss.stats$Sxy
    xcovo <- ss.stats$xcovo
    Mx <- ss.stats$Mx
    nd <- M <- length(SM)
    qt <- names(main.input$SM)
    kappas <- c()
    for (j in 1:length(TOdds)) kappas <- c(kappas, calckappa(nsnps = nsnps, 
        p = 2/nsnps, ndis = nd, target.odds = TOdds[j]))
    kappas <- round(kappas)
    traits <- paste(qt, collapse = "-")
    bestmod.thr <- vector("list", M)
    for (i in 1:M) {
        bm <- best.models.cpp(SM[[i]], cpp.thr = cpp, maxmod)
        bestmod.thr[[i]] <- bm$models
        message("Trait ", i, " (", qt[i], ") ", "has cpp before adjustment: ", 
            bm$old.cpp)
    }
    STR <- lapply(bestmod.thr, "[[", "str")
    PP <- lapply(bestmod.thr, "[[", "PP")
    names(STR) <- qt
    names(PP) <- qt
    SSy <- covY * (Nqq - 1) + Nqq * (ybar %o% ybar)
    for (i in 1:M) mbeta[[i]] <- mbeta[[i]][STR[[i]]]
    pp <- vector("list", length = nd)
    for (kappa in kappas) {
        ret <- marginalpp(STR, PP, mbeta, covY, SSy, Sxy, kappa, 
            N, Nqq, nsnps, Mx, xcovo, Nq3, Nq4, fastapprox, NCORES)
        for (i in 1:nd) pp[[i]] <- cbind(pp[[i]], ret[[i]]$shared.pp)
    }
    for (i in 1:nd) {
        pp[[i]] <- cbind(ret[[i]]$single.pp, pp[[i]])
        colnames(pp[[i]]) <- paste("pp", c("null", round(TOdds, 
            2)), sep = ".")
        rownames(pp[[i]]) <- rownames(ret[[i]])
    }
    mpp <- lapply(pp, MPP.fn)
    names(pp) <- qt
    mpp1 <- lapply(mpp, t)
    MPP <- mpp1[[1]]
    for (k in 2:M) MPP <- gtools::smartbind(MPP, mpp1[[k]], fill = 0)
    return(list(PP = pp, MPP = MPP, sharing = c("null", kappas)))
}



credsetU <- function (modPP, cred = 0.99) 
{
    tmp <- modPP[order(modPP, decreasing = TRUE)]
    cpp <- cumsum(tmp)
    wh <- which(cpp <= cred)
    if (!length(wh)) wh <- 1
    if(cpp[max(wh)] < cred) wh <- c(wh, max(wh) + 1)
    keepmodPP <- tmp[wh]
    mods <- names(keepmodPP)
    cs <- unique(unlist(strsplit(mods, "%")))
    cs <- cs[which(!is.na(cs))]
    return(cs)
}


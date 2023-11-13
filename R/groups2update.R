
#' @title Make SNP groups using fine-mapping information from all of the traits
#' @param main.input output from flashfm.input function
#' @param is.snpmat logical taking value TRUE when genotype matrix is provided and FALSE when covariance matrix is given
#' @param min.mppi trim snp groups with total MPPI < min.mppi in all diseases; default 0.01
#' @param r2.minmerge merge groups with minimum between-group r2 > r2.minmerge; default 0.5
#' @return list of  SNP groups
#' @export
makeSNPgroupsU <- function(main.input,is.snpmat,min.mppi = 0.01,r2.minmerge=0.5) {
snp.data <- main.input$Gmat
SMlist <- main.input$SM
#if(is.snpmat) { Xmat <- new("SnpMatrix",round(snp.data+1)) 
#} else {Xmat <- as.matrix(snp.data) } 
Xmat <- as.matrix(snp.data)
sg <- groupmultiU(SMlist,Xmat,is.snpmat,min.mppi,minsnpmppi=.001,r2.minmerge)

if(length(sg) > 1){
snpgroups <- sg$groups@.Data
ng <- length(snpgroups)
names(snpgroups) <- LETTERS[1:ng] # arbitrary names
if(ng>26) names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
if(ng>52) { 
 names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
 names(snpgroups)[53:min(ng,78)] <- paste0(LETTERS[1:min(26,ng-52)],3)
 }
if(ng>78) {
 names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
 names(snpgroups)[53:min(ng,78)] <- paste0(LETTERS[1:min(26,ng-52)],3)
 names(snpgroups)[79:min(ng,104)] <- paste0(LETTERS[1:min(26,ng-78)],4)
}
Ng <- lapply(snpgroups,length)
sgd <- t(data.frame(Ng)); colnames(sgd) <- "Group Size"
message("SNP group sizes are: "); print(sgd)
} else { snpgroups <- sg}
return(snpgroups)
}


#' @title Make two sets of SNP groups using fine-mapping information from all of the traits using two sets of results and maps the names between them 
#' @param main.input output from flashfm.input function
#' @param fm.multi output from flashfm function
#' @param is.snpmat logical taking value TRUE when genotype matrix is provided and FALSE when covariance matrix is given
#' @param min.mppi trim snp groups with total MPPI < min.mppi in all diseases; default 0.01
#' @param minsnpmppi only group snps with total MPPI > minsnpmppi; default 0.001
#' @param r2.minmerge merge groups with minimum between-group r2 > r2.minmerge; default 0.5
#' @return list of  three objects: groups.fm is a list of SNP groups using the single-trait results; groups.flashfm is a list of SNP groups using the flashfm results; group.sizes is  a table of SNP group sizes for the two sets of groups
#' @export
makeSNPgroups2U <- function(main.input,fm.multi,is.snpmat,min.mppi = 0.01,minsnpmppi=0.001,r2.minmerge=0.5) {
snp.data <- main.input$Gmat
M <- length(fm.multi$PP)
#SMlist <- main.input$SM
#if(is.snpmat) { Xmat <- new("SnpMatrix",round(snp.data+1)) 
#} else {Xmat <- as.matrix(snp.data) } 
SMlist <- vector("list",M)
fmpp  <- fm.multi$PP
for(i in 1:M) {
 ppdf <- data.frame(str=as.character(rownames(fmpp[[i]])),PP=fmpp[[i]][,1], stringsAsFactors = FALSE)
 SMlist[[i]] <- PP2snpmod(ppdf)
 }
Xmat <- as.matrix(snp.data)
sg <- groupmultiU(SMlist,Xmat,is.snpmat,min.mppi,minsnpmppi,r2.minmerge)
if(length(sg) > 1){
snpgroups <- sg$groups@.Data
ng <- length(snpgroups)
names(snpgroups) <- LETTERS[1:ng] # arbitrary names
if(ng>26) names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:(ng-26)],2)
if(ng>52) { 
 names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
 names(snpgroups)[53:min(ng,78)] <- paste0(LETTERS[1:min(26,ng-52)],3)
 }
if(ng>78) {
 names(snpgroups)[27:min(ng,52)] <- paste0(LETTERS[1:min(26,ng-26)],2)
 names(snpgroups)[53:min(ng,78)] <- paste0(LETTERS[1:min(26,ng-52)],3)
 names(snpgroups)[79:min(ng,104)] <- paste0(LETTERS[1:min(26,ng-78)],4)
}
Ng <- lapply(snpgroups,length)
#sgd <- t(data.frame(Ng)); colnames(sgd) <- "Group Size"
#message("SNP group sizes based on single-trait results are: "); print(sgd)
} else { 
 snpgroups <- sg
 ng <- 1
 }

fmpp  <- fm.multi$PP
M <- length(fmpp)
SM2list <- vector("list",M)
for(i in 1:M) {
 ppdf <- data.frame(str=as.character(rownames(fmpp[[i]])),PP=fmpp[[i]][,2], stringsAsFactors = FALSE)
 SM2list[[i]] <- PP2snpmod(ppdf)
 }
sg2 <- groupmultiU(SM2list,Xmat,is.snpmat,min.mppi,minsnpmppi,r2.minmerge)
if(length(sg2) > 1){
snpgroups2 <- sg2$groups@.Data
ng2 <- length(snpgroups2)
names(snpgroups2) <- LETTERS[1:ng2] # arbitrary names
if(ng2>26) names(snpgroups2)[27:min(ng2,52)] <- paste0(LETTERS[1:min(26,ng2-26)],2)
if(ng2>52) { 
 names(snpgroups2)[27:min(ng2,52)] <- paste0(LETTERS[1:min(26,ng2-26)],2)
 names(snpgroups2)[53:min(ng2,78)] <- paste0(LETTERS[1:min(26,ng2-52)],3)
 }
if(ng2>78) {
 names(snpgroups2)[27:min(ng2,52)] <- paste0(LETTERS[1:min(26,ng2-26)],2)
 names(snpgroups2)[53:min(ng2,78)] <- paste0(LETTERS[1:min(26,ng2-52)],3)
 names(snpgroups2)[79:min(ng2,104)] <- paste0(LETTERS[1:min(26,ng2-78)],4)
}

Ng2 <- lapply(snpgroups2,length)
#sgd2 <- t(data.frame(Ng2)); colnames(sgd2) <- "Group Size"
#message("SNP group sizes based on flashfm results are: "); print(sgd2)
} else { 
 snpgroups2 <- sg2
 ng2 <- 1
 }

wh <- NULL
for(i in 1:ng) {
 for(j in 1:ng2) {
  if(length(intersect(snpgroups[[i]],snpgroups2[[j]])) > 0 )wh <- rbind(wh,c(i,j))
 }
}


 
 newsg2 <- snpgroups2  
 newnames <- character(ng2)
 
 ind <- which(1:ng %in% wh[,1])
 wrm <- c() 

 if(any(duplicated(wh[,1]))) {
  dup <- which(duplicated(wh[,1]))
  dups <- c()
  
  for(i in dup) { 
  	dd <- which(wh[,1]==wh[i,1])
  	for(j in 1:length(dd)) newnames[wh[dd[j],2]] <- paste(names(snpgroups)[wh[dd[j],1]],j,sep=".")  			
	wrm <- c(wrm,dd)
	}
   wh <- matrix(wh[-wrm,],ncol=2)
   }
   
   
for(i in 1:nrow(wh)) newnames[wh[i,2]] <- names(snpgroups)[wh[i,1]]
 


nc <- nchar(newnames)
snames <- c(LETTERS[1:26],paste0(LETTERS[1:26],2),paste0(LETTERS[1:26],3),paste0(LETTERS[1:26],4))
if(any(nc==0)){
	ind <- which(nc==0)
	newnames[ind] <-  snames[(ng+1):(ng+length(ind))]	
	}     
names(snpgroups2) <- newnames 
snpgroups2 <- snpgroups2[order(names(snpgroups2))]

Ng <- lapply(snpgroups,length)
sgd <- t(data.frame(Ng)); colnames(sgd) <- "Group Size"
message("SNP group sizes based on single-trait results are: "); print(sgd)

Ng2 <- lapply(snpgroups2,length)
sgd2 <- t(data.frame(Ng2)); colnames(sgd2) <- "Group Size"
message("SNP group sizes based on flashfm results are: "); print(sgd2)

group.sizes <- do.call("smartbind",c(list(t(sgd),t(sgd2)),fill=0))
rownames(group.sizes) <- fm.multi$sharing

return(list(groups.fm=snpgroups,groups.flashfm=snpgroups2, group.sizes=group.sizes))
}





##' @title Group SNPs; adapted from group.multi by Chris Wallace
##' @param SM2 snpmod object
##' @param snp.data SnpMatrix or SNP covariance matrix with snp ids as column names 
##' @param is.snpmat logical taking value TRUE when SnpMatrix is provided and FALSE when covariance matrix is given
##' @param min.mppi trim snp groups with total MPPI < min.mppi in all
##'     diseases
##' @param minsnpmppi only group snps with total MPPI > minsnpmppi
##' @param r2.minmerge merge groups with minimum between-group r2 >
##'     r2.minmerge
##' @return list with three components.
##'
##' First is a data.frame with each row giving summary statistics for
##'     each group.
##'
##' Second is a groups object, each elements ordered according to the rows of the summary
##' 
##' Third is the r2 matrix calculated.
#' @export
groupmultiU <- function (SM2, snp.data, is.snpmat, min.mppi = 0.01, minsnpmppi=0.01, r2.minmerge = 0.6) 
{
    stopifnot(is.list(SM2))
    M <- length(SM2)
    dd <- vector("list",M)
	
    if(is(SM2[[1]], "snpmod")) {
      nsnps <- ncol(snp.data)
      s <- lapply(SM2,function(x) x@snps$var)
	es <- lapply(s,function(s,snps) setdiff(snps,s),snps=colnames(snp.data))
	for(j in 1:M) {
		if(length(es[[j]])>0){
		esmods <- data.frame(str=es[[j]],PP=0,stringsAsFactors = FALSE)
		dd[[j]] <- rbind(SM2[[j]]@models[,c("str","PP")],esmods)
		rownames(dd[[j]]) <- dd[[j]]$str
		SM2[[j]] <- PP2snpmod(dd[[j]])
 				}
 				
		}
    
    }
       
    
    bs <- best.snps(SM2, pp.thr = 0)
    bs <- do.call("rbind", bs)

    snps <- setdiff(unique(bs[bs$Marg_Prob_Incl > minsnpmppi, ]$var), 
        "1")
#    if(length(snps) < 2) {
#     minsnpmppi = 0
#     min.mppi = 0
#     snps <- setdiff(unique(bs[bs$Marg_Prob_Incl > minsnpmppi, ]$var), 
#        "1")
#     }
     if(length(snps) == 1) {
      snpGroups <- list(A=snps)
      out <- snpGroups
     }
     if(length(snps) == 0) {
      minsnpmppi = 0
      min.mppi = 0
      snps <- setdiff(unique(bs[bs$Marg_Prob_Incl > minsnpmppi, ]$var), "1")
     }
     if(length(snps) == 1) {
      snpGroups <- list(A=snps)
      out <- snpGroups
     }
     
     
    if(length(snps) > 1) {    
    if(is.snpmat) {
    	snp.data <- snp.data[, snps]
		r2 <- cor(snp.data)^2
		} else { 
			snp.data <- snp.data[snps,snps]
			r2 <- cov2cor(as.matrix(snp.data))^2 
			}

#	s <- lapply(SM2,function(x) x@snps$var)
	s <- lapply(dd, function(x) rownames(dd))
	es <- lapply(s,function(s,snps) setdiff(snps,s),snps=snps) # prioritised snps not among models for a trait
	

    X <- lapply(SM2, makex)
    mppi <- lapply(X, makemppi)
    
    MPPI <- do.call("smartbind",c(lapply(mppi,function(m) {mm <- matrix(m,nrow=1); colnames(mm) <- names(m); return(mm)}),fill=0))
    MPPI <- t(MPPI)
##    MPPI <- do.call("cbind", lapply(X, makemppi))
    R <- lapply(X, function(x) maker(x)[snps, snps])
    rmax <- rmin <- R[[1]]
    if (length(R) > 1) 
        for (i in 2:length(R)) rmin <- pmin(rmin, R[[i]])
    rmax <- R[[1]]
    if (length(R) > 1) 
        for (i in 2:length(R)) rmax <- pmax(rmax, R[[i]])
    r <- ifelse(abs(rmin) > abs(rmax), rmin, rmax)
    rd <- maked(r, r2)
    h <- hclust(rd, method = "complete")
    d <- as.dendrogram(h)
    r.tol = max(c(quantile(r, 0.9),0))
    
    ## utility functions in local environment
    mem.sum <- function(members) {
        (colSums(MPPI[members, , drop = FALSE]))
    }
    mem.marg <- function(members) {
        sapply(X, function(x) {
            sum(pmin(apply(x$x[, members, drop = FALSE], 1, sum), 
                1) * x$w)
        })
    }
    mem.ab <- function(members) {
        list(a = mem.sum(members), b = mem.marg(members))
    }
    obj.ab <- function(object) {
        members <- labels(object)
        mem.ab(members)
    }
    mem.maxr.minr2 <- function(members) {
        r.sub <- r[members, members, drop = FALSE]
        r2.sub <- r2[members, members, drop = FALSE]
        mx <- max(r.sub[lower.tri(r.sub)], na.rm = TRUE)
        mn <- min(r2.sub[upper.tri(r2.sub)], na.rm = TRUE)
        c(mx, mn)
    }
    cutter <- function(object, mppi.max = 1.01, max.size = 50, 
        marg.sum.ratio = 1.1, max.r = 0, min.r2 = 0.5) {
        if (is.leaf(object)) 
            return(labels(object))
        members <- labels(object)
        ab <- mem.ab(members)
        if (max(ab[[1]]) < min.mppi) 
            return(labels(object))
        if (max(ab[[1]]) > mppi.max) 
            return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
        mxmn <- mem.maxr.minr2(members)
        if (mxmn[1] > r.tol || mxmn[2] < min.r2) 
            return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
        if (min(c(mem.sum(labels(object[[1]])), mem.sum(labels(object[[2]]))) < 
            min.mppi)) 
            return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
        if (max(ab[[1]]) <= mppi.max & all(ab[[1]] < ab[[2]] * 
            marg.sum.ratio)) 
            return(labels(object))
        return(list(cutter(object[[1]],min.r2 =min.r2), cutter(object[[2]],min.r2 =min.r2)))
    }
    mem.summ <- function(members) {
        n <- length(members)
        ab <- mem.ab(members)
        mppi.min <- apply(MPPI[members, , drop = FALSE], 2, min)
        mppi.max <- apply(MPPI[members, , drop = FALSE], 2, max)
        r2.sub <- r2[members, members]
        r2.summ <- summary(r2.sub[upper.tri(r2.sub)])
        r.sub <- r[members, members]
        r.summ <- summary(r.sub[upper.tri(r.sub)])
        c(n = n, sum.mppi = ab[[1]], r2 = r2.summ["Min."], r2 = r2.summ["Max."], 
            r = r.summ["Min."], r = r.summ["Max."], mppi.min = mppi.min, 
            mppi.max = mppi.max)
    }
    ret <- cutter(d,min.r2=r2.minmerge)
    if (!is.list(ret)) 
        ret <- list(ret)
    ret <- LinearizeNestedList(ret)
    ret.mppi <- t(sapply(ret, mem.sum))
    use <- apply(ret.mppi, 1, max) > minsnpmppi
    df <- sapply(ret, mem.summ)
    df <- t(df)
    union.summary <- df[use, , drop = FALSE]
    union.content <- ret[use]
    use <- apply(union.summary[, grep("sum.mppi", colnames(union.summary)), 
        drop = FALSE], 1, max) > min.mppi
    G1 <- union.summary[use, , drop = FALSE]
    G2 <- union.content[use]
    rownames(G1) <- NULL
    merger <- function(G1, G2) {
        maxr2 <- calc.maxmin(r2, G2, fun = max)
        minr2 <- calc.maxmin(r2, G2, fun = min)
        maxr <- calc.maxmin(r, G2, fun = max)
        diag(maxr2) <- 0
        tomerge <- maxr2 > r2.minmerge & maxr < r.tol
        if (any(tomerge, na.rm = TRUE)) {
            wh <- which(tomerge, arr.ind = TRUE)
            wh <- wh[wh[, 1] < wh[, 2], , drop = FALSE]
            wh <- cbind(wh, maxr2[wh])
            wh <- wh[order(wh[, 3], decreasing = TRUE), , drop = FALSE]
            for (k in 1:nrow(wh)) {
                a <- wh[k, 1]
                b <- wh[k, 2]
                sumcols <- grep("sum.mppi", colnames(G1))
                if (any(colSums(G1[c(a, b), sumcols, drop = FALSE]) > 
                  1.01)) 
                  next
                G2[[a]] <- c(G2[[a]], G2[[b]])
                G2[[b]] <- NULL
                for (nm in c(1, sumcols)) G1[a, nm] <- sum(G1[c(a, 
                  b), nm])
                for (nm in setdiff(1:ncol(G1), c(1, sumcols))) G1[a, 
                  nm] <- max(G1[c(a, b), nm])
                G1 <- G1[-b, , drop = FALSE]
                return(merger(G1, G2))
            }
        }
        return(list(G1, G2))
    }
    G <- merger(G1, G2)
    tmp <- G[[2]]
    names(tmp) <- sapply(tmp, "[[", 1)
    newgroups <- new("groups", tmp, tags = names(tmp))
    out <- list(summary = G[[1]], groups = newgroups, r2 = r2)
    }
    
    return(out)
}

##




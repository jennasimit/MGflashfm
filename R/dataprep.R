
#' @title Harmonise list of GWAS data to have same variants, same effect allele, and to remove duplicate snps 
#' @param obsgwas list of GWAS data.frames with column names as specified in the other arguments
#' @param minMAF only variants with MAF > minMAF in all GWAS are retained; default is 0.005
#' @param minINFO only variants with INFO > 0.4 in all GWAS are retained; default is 0.4
#' @param beta_colname text column name for the effect estimates in each obsGWAS data.frame; default "BETA"
#' @param se_colname text column name for the effect estimate standard errors in each obsGWAS data.frame; default "SE"
#' @param snpID_colname text column name for the effect estimates in each obsGWAS data.frame; default "SNP"
#' @param EA_colname text column name for the effect allele in each obsGWAS data.frame; default "EA"
#' @param NEA_colname text column name for the non-effect allele in each obsGWAS data.frame; default "NEA"
#' @param EAfreq_colname text column name for the effect allele frequency (EAF) in each obsGWAS data.frame; default "EAF"
#' @param BP_colname text column name for the base-pair position in each obsGWAS data.frame; default "BP"
#' @param pvalue_colname text column name for the p-value in each obsGWAS data.frame; default "p_value"
#' @param INFO_colname text column name for the INFO score in each obsGWAS data.frame; default "INFO"
#' @param N_colname text column name for the number of individuals with the variant measured in each obsGWAS data frame; default "N".  If there is no sample size column, then there will be no filtering using Nprop
#' @param minNprop Use this for meta-analysis summary statistics - for each trait only variants with N/max(N) > minNprop are retained; default 0.80 (e.g. each variant must be measured in at least 80% of the meta-analysis individuals)
#' @return outputs the list of input GWAS data such that they contain the same variants, have the same effect allele (flip where needed), 
#' and duplicates (by base-pair position) are removed, retaining the variant with highest INFO score
#' @author Jenn Asimit
#' @export
harmoniseGWAS <- function(obsgwas,minMAF=0.005,minINFO=0.4,beta_colname="beta",se_colname="SE",
                       snpID_colname="rsID", EA_colname="EA", NEA_colname="NEA", 
                       EAfreq_colname="EAF", BP_colname="BP", pvalue_colname="p_value", INFO_colname="INFO", N_colname="N", minNprop=0.80) {
# Example for BOLT-LMM output:  
# ooo <- harmoniseGWAS(obsgwas,beta_colname="BETA",EAfreq_colname="A1FREQ",pvalue_colname="P_BOLT_LMM_INF",EA_colname="ALLELE1",NEA_colname="ALLELE0")

 P <- length(obsgwas) # number of observed gwas
 
 for(i in 1:P) {
	obsgwas[[i]] <- as.data.frame(obsgwas[[i]])
	obsgwas[[i]]$INFO=as.numeric(obsgwas[[i]][,INFO_colname])
	obsgwas[[i]]$beta=as.numeric(obsgwas[[i]][,beta_colname])
	obsgwas[[i]]$SE=as.numeric(obsgwas[[i]][,se_colname])
	obsgwas[[i]]$MAF=0.5-abs(0.5-as.numeric(obsgwas[[i]][,EAfreq_colname]))
	obsgwas[[i]]$p_value=as.numeric(obsgwas[[i]][,pvalue_colname])
	obsgwas[[i]]$BP=as.numeric(obsgwas[[i]][,BP_colname])
	obsgwas[[i]]$EA=as.character(obsgwas[[i]][,EA_colname])
	obsgwas[[i]]$NEA=as.character(obsgwas[[i]][,NEA_colname])
	obsgwas[[i]]$EAF=as.numeric(obsgwas[[i]][,EAfreq_colname])
	obsgwas[[i]]$rsID=as.character(obsgwas[[i]][,snpID_colname])
	
	if(N_colname %in% colnames(obsgwas[[i]])) { # If the N column is available in the GWAS, filter out variants that are measured in fewer than minNprop of all individuals 
	 obsgwas[[i]]$N <- as.numeric(obsgwas[[i]][,N_colname])
	 obsgwas[[i]]$Nprop <- obsgwas[[i]]$N/max(obsgwas[[i]]$N)
	 indkeep <- which(obsgwas[[i]]$Nprop >= minNprop)
	 obsgwas[[i]] <- obsgwas[[i]][indkeep,]
	}
	
	indkeep <- which(obsgwas[[i]]$INFO>minINFO & obsgwas[[i]]$MAF>minMAF)
	obsgwas[[i]] <- obsgwas[[i]][indkeep,]
}
 
 for(i in 1:P) { # for each gwas check for and remove any duplicate SNPs (by postion), retaining highest variant with highest INFO score
  	check <- any(duplicated(obsgwas[[i]]$BP)) 
	rmind <- c()
	if(check){
 		dup <- unique(obsgwas[[i]][which(duplicated(obsgwas[[i]]$BP)),"BP"]) # unique duplicate snps
 		for(s in dup){
 			ind <- which(obsgwas[[i]]$BP==s)
  			ksnp <- ind[which.max(obsgwas[[i]]$INFO[ind])] # amongst duplicates, keep snp with max INFO
  			rmind <- c(rmind,setdiff(ind,ksnp))  # store indices of all duplicates that are not having max info
			}
 		obsgwas[[i]] <- obsgwas[[i]][-rmind,] # rm duplicates snps with lower INFO 
	}
 }
# find intersecting snps by BP position
 all_snp_BP <- lapply(obsgwas,function(x) x$BP)
 int_snp_BP <- Reduce("intersect",all_snp_BP)

for(i in 1:P) {
	obsgwas[[i]]<-obsgwas[[i]][which(obsgwas[[i]]$BP %in% int_snp_BP),] # subset to snps present in all gwas
	obsgwas[[i]]<-obsgwas[[i]][order(obsgwas[[i]]$BP),] # order by BP so all GWAS have snps in same order
}

# harmonise gwas so variants are aligned to the same alleles
RPinfo <- obsgwas[[1]] # align to first gwas
# need allele1, allele2, rsID as column names
ind <- which(colnames(RPinfo) %in% c("EA","NEA")) # match these to names in gwas files
colnames(RPinfo)[ind] <- c("allele1","allele2") 

tmp <- obsgwas
obsgwas <- lapply(tmp,alignGWAS,RPinfo=RPinfo)
rm(tmp)

# set same snps in all GWAS
all_snp_BP <- lapply(obsgwas,function(x) x$BP)
int_snp_BP <- Reduce("intersect",all_snp_BP)
for(i in 1:P) {
#	print(i)
	obsgwas[[i]]<-obsgwas[[i]][which(obsgwas[[i]]$BP %in% int_snp_BP),] # subset to snps present in all gwas
	obsgwas[[i]]<-obsgwas[[i]][order(obsgwas[[i]]$BP),] # order by BP so all GWAS have snps in same order
}

return(obsgwas)

}


#' @title Flip GWAS variants to align with SNP correlation matrix (LD) from reference panel or in-sample 
#' @param gwas GWAS data.frame with the following columns (any order, could have others but require these): "rsID" (variant ID), "EA" (effect allele), "NEA" (non-effect allele),
#' "beta" (effect size), "EAF" (effect allele frequency)
#' @param RPinfo reference panel (or in-sample) details data.frame with the following columns (any order): "rsID" (variant ID), "allele1", "allele2" 
#' (this function will align the gwas to allele1, but either allele1 or allele2 may be used -  need consistency for correlation signs)
#' @param details default value is FALSE to return only the aligned GWAS data.frame; TRUE will also return names of excluded snps
#' @return if details=FALSE, return only data.frame of gwas aligned to reference panel (RP) or in-sample LD
#' if details=TRUE: a list with three components: 
#' gwasA = gwas aligned to reference panel (RP) or in-sample LD;
#' excluded = gwas variants that are excluded (not in RP, incompatible alleles (even if flip));
#' ind_excl = indices of removed rows from imput gwas
#' @author Jenn Asimit
#' @export
alignGWAS <- function(gwas,RPinfo,details=FALSE) {

 gwas <- as.data.frame(gwas)
 RPinfo <- as.data.frame(RPinfo)
 snpkeep <- intersect(gwas$rsID,RPinfo$rsID)
 if(length(snpkeep)==0) stop("There is no overlap between the GWAS and reference panel SNP names. Check input.")
 rownames(gwas) <- gwas$rsID
 rownames(RPinfo) <- RPinfo$rsID
 gwas <- gwas[snpkeep,]
 RPinfo <- RPinfo[snpkeep,]
  
 
 flip1 <- which(gwas$EA != RPinfo$allele1 | gwas$NEA != RPinfo$allele2) # check for discrepancies
 check <- excluded <- indrm <- c()
 
 if(length(flip1)>0) check <- which( gwas$EA[flip1] == RPinfo$allele2[flip1] & gwas$NEA[flip1] == RPinfo$allele1[flip1] )  # snps to flip

 if(length(check)>0){
   gwas$NEA[flip1[check]] <- RPinfo$allele2[flip1[check]]
   gwas$EA[flip1[check]] <- RPinfo$allele1[flip1[check]]
   gwas$beta[flip1[check]] <- -gwas$beta[flip1[check]]
   gwas$EAF[flip1[check]] <- 1-gwas$EAF[flip1[check]]

   indrm <-  setdiff(flip1,flip1[check]) # rm these snps from both datasets since allele codings still don't agree if flipped 
   if(length(indrm)>0) {
    excluded <- gwas[indrm,] 
    gwas <- gwas[-indrm,]      
   }
   
 }

if(details) {
out <- list(gwasA=gwas,excluded=excluded,ind_exc=indrm)
} else { out <- gwas}

return(out)

}


#' @title Quality control of dosage genotype matrix prior to LD calculation
#' @param g NxM genotype score matrix for N individuals and M variants
#' @param theta tolerance for dosage, e.g. if theta=0.2, then only dosages d that satisfy d<0.2, 0.8<d<1.2, and d>1.8 are retained 
#' @param BestGuess if TRUE, then return "best guess" genotype matrix (default); if FALSE, return dosage genotype matrix, where dosages that do not meet 
#' tolerance theta, are set to "NA"
#' @return Outputs the "best guess" genotype matrix (if BestGuess=TRUE; default) or the "quality controlled" genotype matrix (if BestGuess=FALSE)
#' @author Jenn Asimit
#' @export
LDqc <- function(g,theta=0.2,BestGuess=TRUE) {
  ind0 <- which(g<=theta)
  ind1 <- which(g>=1-theta & g<=1+theta)
  ind2 <- which(g>=2-theta)
  bg <- rep(NA,length(g))
  if(BestGuess) {
    if(length(ind0)>0) bg[ind0] <- 0
    if(length(ind1)>0) bg[ind1] <- 1
    if(length(ind2)>0) bg[ind2] <- 2 
  } else{
    if(length(ind0)>0) bg[ind0] <- g[ind0]
    if(length(ind1)>0) bg[ind1] <- g[ind1]
    if(length(ind2)>0) bg[ind2] <- g[ind2]   
  }  
  return(bg)
}


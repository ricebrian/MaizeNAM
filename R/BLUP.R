BLUP=function(trait="PlantHeight",family="all",env="all",dereg=TRUE,use_bins=FALSE){
    
    # Line added for debugging purpose
    Z=NamData=B=matrix(NA,2,2)
  
    # Load data
    data(pheno,envir=environment(),package="MaizeNAM")
        
    # Genotypic matrix of lines
    if(use_bins){
      DownloadData()
      data(bins ,envir=environment(),package="MaizeNAM")
      geno = B
      rm(B)
    }else{
      data(geno ,envir=environment(),package="MaizeNAM")
      geno = Z
      rm(Z)
    } 
        
    # FAM
    fam=as.numeric(gsub('Z0|E.+','',rownames(geno)))
    names(fam) = rownames(geno)
    
    # CHR
    tmp = gsub('Bin\\.|S','',colnames(geno))
    tmp = gsub('_.+','',tmp)
    chr=c(table(as.numeric(tmp)))
    
    # Subsetting
    
    if(is.numeric(family)){
      dtaFam = as.numeric(gsub('Z','',as.character(NamData$Family))) 
      NamData = NamData[which(dtaFam%in%family),]
      geno = geno[which(fam%in%family),]
      fam = fam[which(fam%in%family)]
    } 
    
    if(is.numeric(env))  if(length(env)==1) stop("At least two environments where the trait was measured are required")
    
    if(is.numeric(env)){
      E1 = as.numeric(NamData$Environment)
      NamData = NamData[E1%in%env,]
    }
    
    # Model terms
    Y = NamData[,trait]
    G = NamData[,"Genotype"]
    E = NamData[,"Environment"]
    
    # BLUP
    blup = lmer(Y~(1|E)+(1|G))
    BV = coef(blup)$G[,1]
    names(BV) = rownames(coef(blup)$G)
    
    # Matching genotypes
    ge = intersect(names(BV),rownames(geno))
    BV = BV[ge]
    geno = geno[ge,]
    fam = fam[ge]
    
    # Check variance components
    n = c(table(NamData[,"Genotype"]))[ge]
    vc = data.frame(VarCorr(blup))
    vg = vc$vcov[1]
    ve = vc$vcov[3]
    
    cat('Broad-sense H2 of',trait,'is',round(vg/(vg+ve/mean(n)),2),'\n')
    
    # Deregression
    if(dereg){
      h2 = vg/(vg+ve/n)
      BV = BV/h2
    }
        
    LIST = list('Phen'=BV,'Gen'=geno,'Chrom'=chr,'Fam'=fam)
    return(LIST)
  }

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

DownloadData = function( dataset = "all", redownload = FALSE ){
  
  # c("all","panel","bins","gbs")
  
  gbs = paste0(system.file("data","",package="MaizeNAM"),'/gbs.RData')
  panel = paste0(system.file("data","",package="MaizeNAM"),'/panel.RData')
  bins = paste0(system.file("data","",package="MaizeNAM"),'/bins.RData')
  
  if( ! "gbs.RData" %in% dir(system.file("data","",package="MaizeNAM")) | redownload  ){
    cat('Downloading GBS data \n')
    if(dataset=="all"|what=="gbs") download.file("https://urldefense.com/v3/__https://github.com/ricebrian/MaizeNAMfullGenoData/raw/main/MaizeNAMxz.RData__;!!DZ3fjg!vsyljF2rk0eeNeYWTOzjsRjieOd7J_WbZpGb5h-eImKJ_Qw4sGFgJQ2ncbmpsgKnpg$ ", gbs)
  }
  
  if( ! "panel.RData" %in% dir(system.file("data","",package="MaizeNAM")) | redownload ){
    cat('Downloading diversity panel data \n')
    if(dataset=="all"|what=="panel") download.file("https://urldefense.com/v3/__https://github.com/ricebrian/MaizeNAMfullGenoData/raw/main/DiversityPanelData.Rdata__;!!DZ3fjg!vsyljF2rk0eeNeYWTOzjsRjieOd7J_WbZpGb5h-eImKJ_Qw4sGFgJQ2ncbkXlyTsKA$ ", panel)
  }
  
  if( ! "bins.RData" %in% dir(system.file("data","",package="MaizeNAM")) | redownload  ){
    cat('Downloading genotypic bins data \n')
    if(dataset=="all"|what=="bins") download.file("https://urldefense.com/v3/__https://github.com/ricebrian/MaizeNAMfullGenoData/raw/main/bins.RData__;!!DZ3fjg!vsyljF2rk0eeNeYWTOzjsRjieOd7J_WbZpGb5h-eImKJ_Qw4sGFgJQ2ncbm2DvTIyg$ ", bins)
  }
  
  
}

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


LoadData = function( dataset = "panel" ){
  
  if(dataset == "panel"){
    if( "panel.RData" %in% dir(system.file("data","",package="MaizeNAM"))  ){
      dta = paste0(system.file("data","",package="MaizeNAM"),'/panel.RData')
      load(file = dta, envir =  .GlobalEnv )
    }else{
      cat('You need to download the data first! Run this:   DownloadData()   ')
    }
  }
  
  ###
  
  if(dataset == "gbs"){
    if( "gbs.RData" %in% dir(system.file("data","",package="MaizeNAM"))  ){
      dta = paste0(system.file("data","",package="MaizeNAM"),'/gbs.RData')
      load(file = dta, envir =  .GlobalEnv )
    }else{
      cat('You need to download the data first! Run this:   DownloadData()   ')
    }
  }
  
  ###
  
  if(dataset == "bins"){
    if( "bins.RData" %in% dir(system.file("data","",package="MaizeNAM"))  ){
      dta = paste0(system.file("data","",package="MaizeNAM"),'/bins.RData')
      load(file = dta, envir =  .GlobalEnv )
    }else{
      cat('You need to download the data first! Run this:   DownloadData()   ')
    }
  }
  
  
}


# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# Extended Flexible Gwas
gwnam = function(pheno,geno,pop){
  tag = 0; numSnps = ncol(geno)
  pb = txtProgressBar(style = 3)
  sma = apply(geno,2,function(x,y,pop){
    # Print marker under evaluation
    tag <<- tag+1
    TmpDta = data.frame(y=y,f=factor(pop),x=x)
    lvl = levels(TmpDta$f)
    TmpDta = droplevels.data.frame(TmpDta[rowMeans(!is.na(TmpDta))==1,])
    Vy = c(var(TmpDta$y))
    # Null model
    fit0 = lmer(y~(1|f),TmpDta)
    ll0 = logLik(fit0)
    ## Model 1 - Within-family effect
    fit1 = suppressMessages(lmer(y~(1|f)+(1|f):x,TmpDta))
    eff1 = ranef(fit1)$f[,2]
    names(eff1) = rownames(ranef(fit1)$f)
    eff1 = eff1[lvl]
    ll1 = logLik(fit1)
    LRT1 = ll1-ll0
    PVAL1 = -log10(1-pchisq(LRT1,1))
    ## Model 2 - Across-family effect
    eff2 = suppressMessages(lmer(y~x+(1|f),TmpDta))@beta[2]
    ## Coeff of determination
    R2 = 1-c(fit0@devcomp$cmp['sigmaREML'],
             fit1@devcomp$cmp['sigmaREML'])/Vy
    names(R2)=paste0('R2.model',0:1)
    ## Output
    NumObs = nrow(TmpDta)
    out = c(NumObs,P1=PVAL1,R2,Fxd=eff2,eff1)
    setTxtProgressBar(pb, tag/numSnps)
    return(out)},y=pheno,pop=pop)
  rownames(sma) = c('NumObs','MinusLogPvalue','NullModel_R2','AltrModel_R2','OverallSnpEffect',
                    paste('SnpEffPop',sort(unique(pop)),sep=''))
  close(pb); sma[is.na(sma)] = 0
  return(data.frame(t(sma)))}

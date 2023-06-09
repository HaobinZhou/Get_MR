library(mr.raps)
library(TwoSampleMR)
library(dplyr)
library(fs) 
library(ieugwasr)
library(plinkbinr) 
library(tidyr)
library(data.table)
library(MRPRESSO)
library(parallel)
library(doParallel)
library(pbapply)
library(stringr)
library(cause)
library(vroom)
library(MungeSumstats)
library(ldscr)
library(MRlap)


select<-dplyr::select
mcl<-makeCluster


## 2.0更新内容

clean<-function(list,clean_0=FALSE){
  
  l<-pblapply(list,class)
  ll<-c()
  for(i in 1:length(l)){
    if(l[[i]][1]%in% "try-error"){}
    else{ll<-c(ll,i)}
  }
  list<-list[ll]
  
  if(clean_0==TRUE){
    list<-Filter(function(x) nrow(x)!=0, list)
  }
  return(list)
}

get_rsid<-function(chr,pos,version='hg38'){
  library(biomaRt)
  dat<-data.frame(CHR=chr,start=pos,end=pos)
  if(version=="hg38"){
    ver<-listEnsemblArchives()%>%subset(name!='Ensembl GRCh37')
    ur<-ver$url[1]
    snp_mart <- useEnsembl(biomart="ENSEMBL_MART_SNP", 
                           host=ur, 
                           dataset="hsapiens_snp")
    position <- apply(dat, 1, paste, collapse = ":")
    res<-getBM(attributes = c('refsnp_id', 'allele', 'chrom_start'), 
               filters = 'chromosomal_region', 
               values = position, 
               mart = snp_mart)
  }
  if(version=="hg19"){
    snp_mart <- useEnsembl(biomart="ENSEMBL_MART_SNP", 
                           host="https://grch37.ensembl.org", 
                           dataset="hsapiens_snp")
    position <- apply(dat, 1, paste, collapse = ":")
    res<-getBM(attributes = c('refsnp_id', 'allele', 'chrom_start'), 
               filters = 'chromosomal_region', 
               values = position, 
               mart = snp_mart)
  }
  
  return(res)
}

get_exposure = function(id,pval=5e-08,r2=0.001,kb=10000) {
  library(TwoSampleMR)
  dat<-try(data.frame(extract_instruments(id,p1=pval,clump = T,p2=1,r2=r2,kb=kb)))
  if ('try-error'%in%class(dat))   {
    dat=NA
  }
  return(dat)
}

get_outcome=function(id,expo){
  library(TwoSampleMR)
  dat<-try(data.frame(extract_outcome_data(expo$SNP,id,proxies = F,maf_threshold=0.4)))
  if ('try-error'%in%class(dat))   {
    dat=NA
  }
  
  try( return(dat))
}

get_ao<-function(a=NULL){
  ao<-available_outcomes()
  ao<-separate(ao,id,c("a","b","c"),sep="-",remove = FALSE)
  if(is.null(a)==F){ao<-subset(ao,a==a)}
  return(ao)
}

get_exposure_wrong_num<-function(list){
  res<-is.na(list)
  num=0
  foreach (i=1:length(res))%do%{
    if(res[i]==TRUE) num=num+1
  }
  return(num)
}

clean_outcome_from_exposure<-function(expo,outcome){
  snp<-lapply(expo,FUN=function(x)data.frame(SNP=x$SNP))
  snp<-bind_rows(snp)
  snp<-data.frame(SNP=unique(snp$SNP))
  outcome_adj<-merge(outcome,snp)
  return(outcome_adj)
}

clean_GWAS<-function(list,clean=c("bbj","eqtl")){
  clean_GWAS_logit<-function(dat,clean){
    id<-dat$id.exposure[1]
    id<-as.data.frame(id)
    id<-separate(id,"id",c('a','b','c'),sep = "-")
    if(id$a==clean){
      return(data.frame(lg=1))
    }
    else(return(data.frame(lg=0)))
  }
  print("正在获取指定数据集信息")
  for(i in 1:length(clean)){
  lg<-pblapply(list,clean[i],FUN=clean_GWAS_logit)
  
  lg<-bind_rows(lg)
  
  list<-list[which(lg$lg==0)]
  }
  return(list)
  
}



cyclemr<-function(dat,cl_num,type="list"){
clean<-function(list,clean_0=FALSE){
  
  l<-pblapply(list,class)
  ll<-c()
  for(i in 1:length(l)){
    if(l[[i]][1]%in% "try-error"){}
    else{ll<-c(ll,i)}
  }
  list<-list[ll]
  
  if(clean_0==TRUE){
    list<-Filter(function(x) nrow(x)!=0, list)
  }
  return(list)
}
  # base mr
  mr_base<-function(dat){
    library(TwoSampleMR)
    library(dplyr)
    try(mr<-mr(dat) )
    try(mr<-dplyr::select(mr,id.exposure,id.outcome,method,nsnp,b,se,pval)) 
    try(mr_OR<-generate_odds_ratios(mr_res = mr(dat)))
    try(mr_OR<-dplyr::select(mr_OR,id.exposure,id.outcome,method,nsnp,b,se,pval,lo_ci,up_ci,or, or_lci95,or_uci95))
    try(mr_OR<-dplyr::rename(mr_OR,b.OR="b",se.OR="se",pval.OR="pval"))
    try( mr_p_OR<-merge(mr,mr_OR))
    try(return(mr_p_OR))
  }
  
  mr_egger<-function(dat){
    library(TwoSampleMR)
    library(dplyr)
    mr_egger<-mr_pleiotropy_test(dat) 
    try( mr_egger<-dplyr::select(mr_egger,id.exposure,id.outcome,egger_intercept,se,pval))
    try( mr_egger<-dplyr::rename(mr_egger,se.egger="se",pval.egger="pval"))
    try( mr_egger[2:5,]<-NA) 
    try( return(mr_egger))
  }
  
  mr_test<-function(dat){
    library(TwoSampleMR)
    library(dplyr)
    mr_heterogeneity<-mr_heterogeneity(dat)
    try (mr_heterogeneity<-dplyr::select(mr_heterogeneity,id.exposure,id.outcome,method,Q, Q_df,Q_pval))
    try( mr_heterogeneity<-dplyr::rename(mr_heterogeneity,method.he="method"))
    try (mr_heterogeneity[3:5,]<-NA )
    try( return(mr_heterogeneity))
  }
  
  # presso
  cycle_presso<-function(dat){
    library(TwoSampleMR)
    library(MRPRESSO)
    library(dplyr) 
    nsnp_filter=6
    try (mr_presso_res<-mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                  OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat,  
                                  SignifThreshold = 0.05))
    try ( mr_presso_main<-mr_presso_res$`Main MR results`)
    try(mr_presso_main)
    try ( mr_presso_main[3:5,]<-NA)
    try(  return(mr_presso_main))
  }
  
  mr_Presso<-function(dat,num=10000){
    library(TwoSampleMR)
    library(MRPRESSO)
    library(dplyr)
    
    nsnp_filter=6
    set.seed(123)
    try (mr_presso_res<-mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                  OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat,  
                                  SignifThreshold = 0.05,NbDistribution = num))
    return(mr_presso_res)
    
  }
  mr_presso_pval<-function(mr_presso_res){ 
    try ( mr_presso_main<-mr_presso_res$`Main MR results`)
    try ( mr_presso_main[3:5,]<-NA) 
    return(mr_presso_main)
  }
  
  
  mr_presso_snp<-function(mr_presso_res,mr_presso_main,dat,type="list"){
    data_re<-list()
    if(type=="list"){
      for(i in 1:length(mr_presso_res)){
        res<-mr_presso_res$`MR-PRESSO results`[[i]]
        main<-mr_presso_main[[i]]
        data<-dat[[i]]
        try(if(is.na(main[2,6])==FALSE){
          outliers<-which(res$`Outlier Test`$Pvalue<0.05)
          data$mr_keep[outliers]<-FALSE
        })
        data_re[[i]]<-data
        names(data_re)[[i]]<-names(dat)[[i]]
      }
      return(data_re)
    }
    
    if(type=="data"){
      res<-mr_presso_res$`MR-PRESSO results`
      main<-mr_presso_main
      data<-dat
      try(if(is.na(main[2,6])==FALSE){
        outliers<-which(res$`Outlier Test`$Pvalue<0.05)
        data$mr_keep[outliers]<-FALSE
      })
      return(data)
    }
  }
  bind_basemr<-function(base.res,egger.res,test.res){
    res_all<-data.frame()
    for(i in 1:nrow(base.res)){
      if(i%%5==1){
        id<-base.res$id.exposure[i]
        num.bg<-which(egger.res$id.exposure==id)
        num.end<-num.bg+4
        res<-cbind(base.res[i:(i+4),],egger.res[num.bg:num.end,3:5])
        num.bg<-which(test.res$id.exposure==id)[1]
        num.end<-num.bg+4
        res<-cbind(res,test.res[num.bg:num.end,3:6])
        res_all<-rbind(res_all,res)
      }
    }
    return(res_all)
  }
  
  if(type=="list"){
    cl <- makeCluster(cl_num) 
    
    base.res<-pblapply(dat, FUN=mr_base,cl=cl)
    
    egger.res<-pblapply(dat, FUN=mr_egger)
    
    test.res<-pblapply(dat,FUN=mr_test)
    
    stopCluster(cl)
    
    base.res<-clean(base.res)
    
    egger.res<-clean(egger.res)
    
    test.res<-clean(test.res)
    
    base.res<-bind_rows (base.res)
    
    base.res<-subset(base.res,nsnp>2) ## 防止出现wald ratio，导致错误
    
    egger.res<-bind_rows (egger.res)
    
    test.res<-bind_rows (test.res)
    
    res_all<-bind_basemr(base.res,egger.res,test.res)
    
  }
  
  if(type=="data"){
    
    presso_res<-mr_Presso(dat,10000)
    
    presso_pval<-mr_presso_pval(presso_res)
    
    dat_aj<-mr_presso_snp(presso_res,presso_pval,dat,"data")
    
    base.res<-mr_base(dat)
    
    egger.res<-mr_egger(dat)
    
    test.res<-mr_test(dat)
    
    base.res<-bind_rows (base.res)
    
    res_all<-cbind(base.res,presso_pval,egger.res,test.res)
  }
  
  return(res_all)
}


bind_basemr<-function(base.res,egger.res,test.res){
  res_all<-data.frame()
  for(i in 1:nrow(base.res)){
    if(i%%5==1){
      id<-base.res$id.exposure[i]
      num.bg<-which(egger.res$id.exposure==id)
      num.end<-num.bg+4
      res<-cbind(base.res[i:(i+4),],egger.res[num.bg:num.end,3:5])
      num.bg<-which(test.res$id.exposure==id)[1]
      num.end<-num.bg+4
      res<-cbind(res,test.res[num.bg:num.end,3:6])
      res_all<-rbind(res_all,res)
    }
  }
  return(res_all)
}


bind_pressomr<-function(base_res,presso_pval,egger_res,test_res){
  
  presso_id<-names(presso_pval)
  presso_pval<-bind_rows(presso_pval)
  base_res<-bind_rows(base_res)
  egger_res<-bind_rows(egger_res)
  test_res<-bind_rows(test_res)
  
  res_all<-data.frame()
  for(i in 1:nrow(base_res)){
    if(i%%5==1){
      id<-base_res$id.exposure[i]
      
      num.bg<-5*(which(presso_id==id)[1]-1)+1
      num.end<-num.bg+4
      res<-cbind(base_res[i:(i+4),],presso_pval[num.bg:num.end,2:6])
      
      num.bg<-which(egger_res$id.exposure==id)
      num.end<-num.bg+4
      res<-cbind(res,egger_res[num.bg:num.end,3:5])
      
      num.bg<-which(test_res$id.exposure==id)[1]
      num.end<-num.bg+4
      res<-cbind(res,test_res[num.bg:num.end,3:6])
      
      res_all<-rbind(res_all,res)
    }
  }
  return(res_all)
}


bind_trait<-function(res_all,ao){
  res_all<-data.frame(trait=NA,res_all)
  for(i in 1:nrow(res_all)){
    try( num<-which(ao$id==res_all$id.exposure[i]))
    try( res_all$trait[i]<-ao$trait[num])
  }
  return(res_all)
}

bind_trait_outcome<-function(res_all,ao){
  res_all<-data.frame(trait=NA,res_all)
  for(i in 1:nrow(res_all)){
    num<-which(ao$id==res_all$id.outcome[i])
    res_all$trait[i]<-ao$trait[num]
  }
  return(res_all)
}

clean_cl<-function(){
  cl<<-makeCluster(1)
  stopCluster(cl)
  gc()
}

clean_GWAS_id<-function(dat,id){
  lid<-lapply(dat,FUN = function(x)x$id.exposure[1])
  lid<-unlist(lid)
  l<-c()
  for (i in 1:length(id)) {
    l1<-which(lid==id[i])
    l<-c(l,l1)
  }
  dat<-dat[l]
  return(dat)
}


## 1.0内容
# format
format_Mun<-function(file,source="finn_r8",save_path=NULL,lift=F,ref_genome = "hg38",
                     convert_ref_genome = "hg19"){
  library(data.table)
  library(MungeSumstats)
  library(dplyr)
  if(class(file)!="data.frame"){
    dat<-fread(file)
    if(source=="finn_r8"){
      dat<-dat%>%dplyr::select(SNP=rsids,CHR=`#chrom`,BP=pos,A1=ref,
                               A2=alt,FRQ=af_alt,BETA=beta,
                               SE=sebeta,P=pval)
    }
  }
  dat<-as.data.frame(dat)
  dat<-format_sumstats(dat,return_data = TRUE)
  if(is.null(save_path)==F){setwd(save_path)}
  dat<-fwrite(dat,paste0(file,'(Mun_format',ref_genome,').gz'))
  te<-fs::dir_info( tempdir())
  te<-subset(te,size>1e+8)
  fs::file_delete(te$path)
  if(lift==T){
    dat <- MungeSumstats::liftover(sumstats_dt = dat, 
                                   ref_genome = ref_genome,
                                   convert_ref_genome = convert_ref_genome)
    dat<-fwrite(dat,paste0(file,'(Mun_format',convert_ref_genome,').gz'))
  }
  gc()
}

format_cyclemr<-function(data,type="exposure",source="finn_r8"){
  library(TwoSampleMR)
  library(MungeSumstats)
  
  if(source=="finn_r8"){
    data<-format_data(data,type=type,
                      id_col = "id",chr_col ="#chrom",pos="pos",
                      snp_col = "rsids",beta_col = "beta",se_col = "sebeta",
                      effect_allele_col = "alt",other_allele_col = "ref",
                      eaf_col = "af_alt",phenotype_col = "phenotype",pval_col = "pval")
  }
  
  
  
  if(source=="ukb_nosnp"){
    data<-separate(data,variant,c("chr","pos","ref","alt"),sep = ":")
    
    data<-dplyr::select(data,chr,pos,ref,alt,minor_AF,beta,se,pval)
    name<-colnames(data)
    colnames(data) <-c('CHR','BP','A1','A2','FRQ','BETA', 'SE','P')
    
    expo_rs_done<-format_sumstats(data,ref_genome = "GRCh37",return_data = TRUE)
    
    expo_rs_done$BP<-as.character(expo_rs_done$BP)
    colnames(data)<-name
    data<-merge(data,expo_rs_done,by.x=c("chr","pos"),by.y=c("CHR","BP"),all.x=TRUE)
    
    data<-format_data(data,type=type,
                      samplesize_col = "n_complete_samples",
                      snp_col = "SNP",effect_allele_col = "alt",
                      other_allele_col = "ref",eaf_col ="minor_AF",
                      beta_col="beta",se_col = "se",pval_col = "pval",
                      chr_col="chr",pos_col ="pos",id='id',phenotype_col = 'phenotype')
  }
  
  if(source=="Mun"){
    data<-format_data(data,type=type,
                      snp_col = 'SNP',
                      chr_col = 'CHR',
                      pos_col = 'BP',
                      effect_allele_col = 'A2',
                      other_allele_col = "A1",
                      se_col = 'SE',
                      beta_col= 'BETA',
                      eaf_col = 'FRQ',
                      id_col = 'id',
                      phenotype_col = "phenotype",
                      pval_col = "P"
    )
  }
  
  if(source=="covid"){
    data<-format_data(data,type=type,
                      id_col = "id",chr_col ="#CHR",pos="POS",
                      snp_col = "rsid",beta_col = "all_inv_var_meta_beta",se_col = "all_inv_var_meta_sebeta",
                      effect_allele_col = "ALT",other_allele_col = "REF",
                      eaf_col = "all_meta_AF",phenotype_col = "phenotype",
                      ncase_col = "all_inv_var_meta_cases",ncontrol_col = "all_inv_var_meta_controls",
                      pval_col = "all_inv_var_meta_p")
  }
  
  if(source=="outcome"){
    data<-format_data(data,type=type,
                      id_col = "id.outcome",chr_col ="chr.outcome",pos="pos.outcome",
                      snp_col = "SNP",beta_col = "beta.outcome",se_col = "se.outcome",
                      effect_allele_col = "effect_allele.outcome",other_allele_col = "other_allele.outcome",
                      eaf_col = "eaf.outcome",phenotype_col = "outcome",
                      samplesize_col = "samplesize.outcome",
                      pval_col = "pval.outcome")
    
  }
  if(source=="exposure"){
    data<-format_data(data,type=type,
                      id_col = "id.exposure",chr_col ="chr.exposure",pos="pos.exposure",
                      snp_col = "SNP",beta_col = "beta.exposure",se_col = "se.exposure",
                      effect_allele_col = "effect_allele.exposure",other_allele_col = "other_allele.exposure",
                      eaf_col = "eaf.exposure",phenotype_col = "exposure",
                      samplesize_col = "samplesize.exposure",
                      pval_col = "pval.exposure")
    
  }
  
  if(source=="fast_ukb"){
    data<-format_data(data,type=type,id_col="id",
                      phenotype_col = "phenotype",snp_col = "SNP",
                      effect_allele_col = "A1",other_allele_col = "A2",
                      eaf_col ="AF1",beta_col="BETA",se_col = "SE",
                      pval_col = "P",chr_col = "CHR",pos_col = "POS")
  }
  
  
  
  if(source=="bac"){
    data<-format_data(data,type=type,id_col="bac",
                      phenotype_col = "phenotype",snp_col = "rsID",
                      effect_allele_col = "eff.allele",other_allele_col = "ref.allele",
                      beta_col="beta",se_col = "SE",
                      pval_col = "P.weightedSumZ",chr_col = "chr",pos_col = "bp",samplesize_col = "N")
  }
  
  if(source=="finn_r7"){
    data<-format_data(data,type=type,snp_col = "rsids",
                      effect_allele_col = "alt",other_allele_col = "ref",
                      beta_col="beta",se_col = "sebeta",eaf_col = "af_alt",
                      pval_col = "pval",chr_col = "chrom",pos_col = "pos")
  }
  return(data)
}


format_trait<-function(list,short=FALSE,short_num="40"){
  for(i in 1:length(list)){
    expo<-data.frame(exposure=list[[i]]
                     $exposure)
    try(expo<-separate(expo,exposure,c('a','b'),sep='\\|\\|'))
    if(short==TRUE){try(expo$a <- substr(expo$a, 1, short_num))}
    try(expo<-expo$a)
    try(expo<-gsub(":","_",expo))
    try(expo<-gsub("-","_",expo))
    try(expo<-gsub(" ","_",expo))
    try(expo<-gsub(",","_",expo))
    try(expo<-gsub("/","_",expo))
    list[[i]]$exposure<-expo[1]
    
  }
  
  
  return(list)
}

# read
read_vcf_getmr<-function(file_name,nThread = 8,type=".gz"){
  name<-file_name
  
  for(i in 1:nrow(name)){
    dat<-read_sumstats(paste0("./",name[i]),nThread = nThread,nrow=Inf,standardise_headers = FALSE,mapping_file = sumstatsColHeaders)
    
    vroom_write(dat,paste0(name[i],type))
    
    gc()
  }
  print(i)
  
}

read_easy<-function(file_name,pval=5e-08){
  library(data.table)
  dat<-fread(file_name)
  dat<-subset(dat,pval.exposure<pval)
  return(dat)
}


# get
get_eaf_from_1000G<-function(dat,path,type="exposure"){
  
  corrected_eaf_expo<-function(data_MAF){
    effect=data_MAF$effect_allele.exposure
    other=data_MAF$other_allele.exposure
    A1=data_MAF$A1
    A2=data_MAF$A2
    MAF_num=data_MAF$MAF
    EAF_num=1-MAF_num
    
    harna<-is.na(data_MAF$A1)
    harna<-data_MAF$SNP[which(harna==T)]
    
    cor1<-which(data_MAF$effect_allele.exposure !=data_MAF$A1)
    
    
    data_MAF$eaf.exposure=data_MAF$MAF
    data_MAF$type="raw"
    data_MAF$eaf.exposure[cor1]=EAF_num[cor1]
    data_MAF$type[cor1]="corrected"
    cor2<-which(data_MAF$other_allele.exposure ==data_MAF$A1)
    cor21<-setdiff(cor2,cor1)
    cor12<-setdiff(cor1,cor2)
    error<-c(cor12,cor21)
    data_MAF$eaf.exposure[error]=NA
    data_MAF$type[error]="error"
    
    data_MAF<-list(data_MAF=data_MAF,cor1=cor1,harna=harna,error=error)
    
    return(data_MAF)
    
  }
  
  corrected_eaf_out<-function(data_MAF){
    effect=data_MAF$effect_allele.outcome
    other=data_MAF$other_allele.outcome
    A1=data_MAF$A1
    A2=data_MAF$A2
    MAF_num=data_MAF$MAF
    EAF_num=1-MAF_num
    
    harna<-is.na(data_MAF$A1)
    harna<-data_MAF$SNP[which(harna==T)]
    
    cor1<-which(data_MAF$effect_allele.outcome !=data_MAF$A1)
    
    
    data_MAF$eaf.outcome=data_MAF$MAF
    data_MAF$type="raw"
    data_MAF$eaf.outcome[cor1]=EAF_num[cor1]
    data_MAF$type[cor1]="corrected"
    cor2<-which(data_MAF$other_allele.outcome ==data_MAF$A1)
    cor21<-setdiff(cor2,cor1)
    cor12<-setdiff(cor1,cor2)
    error<-c(cor12,cor21)
    data_MAF$eaf.outcome[error]=NA
    data_MAF$type[error]="error"
    
    data_MAF<-list(data_MAF=data_MAF,cor1=cor1,harna=harna,error=error)
    
    return(data_MAF)
    
  }
  
  if(type=="exposure" && (is.na(dat$eaf.exposure[1])==T || is.null(dat$eaf.exposure)==T)){
    r<-nrow(dat)
    
    setwd(path)
    MAF<-fread("fileFrequency.frq",header = T)
    
    dat<-merge(dat,MAF,by.x = "SNP",by.y = "SNP",all.x = T)
    
    dat<-corrected_eaf_expo(dat)
    
    cor1<-dat$cor1
    
    harna<-dat$harna
    
    error<-dat$error
    
    dat<-dat$data_MAF
    
    print(paste0("一共有",(r-length(harna)-length(error)),"个SNP成功匹配EAF,占比",(r-length(harna)-length(error))/r*100,"%"))
    
    print(paste0("一共有",length(cor1),"个SNP是major allele，EAF被计算为1-MAF,在成功匹配数目中占比",length(cor1)/(r-length(harna)-length(error))*100,"%"))
    
    print(paste0("一共有",length(harna),"个SNP在1000G中找不到，占比",length(harna)/r*100,"%"))
    
    print(paste0("一共有",length(error),"个SNP在输入数据与1000G中效应列与参照列，将剔除eaf，占比",length(error)/r*100,"%"))
    
    print("输出数据中的type列说明：")
    print("raw：EAF直接等于1000G里的MAF数值，因为效应列是minor allele")
    print('corrected：EAF等于1000G中1-MAF，因为效应列是major allele')
    print("error：输入数据与1000G里面提供的数据完全不一致，比如这个SNP输入的效应列是C，参照列是G，但是1000G提供的是A-T，这种情况下，EAF会被清空（NA），当成匹配失败")
    
    return(dat)
  }
  
  if(type=="outcome" && (is.na(dat$eaf.outcome[1])==T || is.null(dat$eaf.outcome)==T)){
    r<-nrow(dat)
    
    setwd(path)
    MAF<-fread("fileFrequency.frq",header = T)
    
    dat<-merge(dat,MAF,by.x = "SNP",by.y = "SNP",all.x = T)
    
    dat<-corrected_eaf_out(dat)
    
    cor1<-dat$cor1
    
    harna<-dat$harna
    
    error<-dat$error
    
    dat<-dat$data_MAF
    
    print(paste0("一共有",(r-length(harna)-length(error)),"个SNP成功匹配EAF,占比",(r-length(harna)-length(error))/r*100,"%"))
    
    print(paste0("一共有",length(cor1),"个SNP是major allele，EAF被计算为1-MAF,在成功匹配数目中占比",length(cor1)/(r-length(harna)-length(error))*100,"%"))
    
    print(paste0("一共有",length(harna),"个SNP在1000G找不到，占比",length(harna)/r*100,"%"))
    
    print(paste0("一共有",length(error),"个SNP在输入数据与1000G中效应列与参照列，将剔除eaf，占比",length(error)/r*100,"%"))
    
    print("输出数据中的type列说明：")
    print("raw：EAF直接等于1000G里的MAF数值，因为效应列是minor allele")
    print('corrected：EAF等于1000G中1-MAF，因为效应列是major allele')
    print("error：输入数据与1000G里面提供的数据完全不一致，比如这个SNP输入的效应列是C，参照列是G，但是1000G提供的是A-T，这种情况下，EAF会被清空（NA），当成匹配失败")
    
    return(dat)
  }
  else{return(dat)}
}

get_chr_pos<-function(dat,type="exposure"){
  if(type=="exposure"){
    if(is.na(dat$eaf.exposure[1])==T || is.null(dat$eaf.exposure)==T){
      print("需要先运行get_eaf_from_1000G来匹配eaf，再匹配chr和pos")
    }
    else{
      dat_d<-dat%>%select(SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
      colnames(dat_d) <-c('SNP','A1','A2','FRQ','BETA', 'SE','P')
      dat_d<-format_sumstats(dat_d,ref_genome = "GRCh37",return_data = TRUE)
      dat_d<-format_data(dat_d,type=type,
                         snp_col = 'SNP',
                         chr_col = 'CHR',
                         pos_col = 'BP',
                         effect_allele_col = 'A1',
                         other_allele_col = "A2",
                         se_col = 'SE',
                         beta_col= 'BETA',
                         eaf_col = 'FRQ',
                         pval_col = "P"
      )
      dat_d<-dat_d%>%select(SNP,chr.exposure,pos.exposure)
      dat<-merge(dat,dat_d,by="SNP",all.x=T)
      
      return(dat)
    }
  }
  if(type=="outcome"){
    if(is.na(dat$eaf.outcome[1])==T || is.null(dat$eaf.outcome)==T){
      print("需要先运行get_eaf_from_1000G来匹配eaf，再匹配chr和pos")
    }
    else{
      dat_d<-dat%>%select(SNP,effect_allele.outcome,other_allele.outcome,eaf.outcome,beta.outcome,se.outcome,pval.outcome)
      colnames(dat_d) <-c('SNP','A1','A2','FRQ','BETA', 'SE','P')
      dat_d<-format_sumstats(dat_d,ref_genome = "GRCh37",return_data = TRUE)
      dat_d<-format_data(dat_d,type=type,
                         snp_col = 'SNP',
                         chr_col = 'CHR',
                         pos_col = 'BP',
                         effect_allele_col = 'A1',
                         other_allele_col = "A2",
                         se_col = 'SE',
                         beta_col= 'BETA',
                         eaf_col = 'FRQ',
                         pval_col = "P"
      )
      dat_d<-dat_d%>%select(SNP,chr.outcome,pos.outcome)
      dat<-merge(dat,dat_d,by="SNP",all.x=T)
      
      return(dat)
    }
    
  }
}

get_f<-function(dat,F_value=10){
  log<-is.na(dat$eaf.exposure)
  log<-unique(log)
  if(length(log)==1)
  {if(log==TRUE){
    print("数据不包含eaf，无法计算F统计量")
    return(dat)}
  }
  if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("数据不包含beta，无法计算F统计量")
    return(dat)}
  if(is.null(dat$se.exposure[1])==T || is.na(dat$se.exposure[1])==T){print("数据不包含se，无法计算F统计量")
    return(dat)}
  if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("数据不包含samplesize(样本量)，无法计算F统计量")
    return(dat)}
  
  
  if("FALSE"%in%log && is.null(dat$beta.exposure[1])==F && is.na(dat$beta.exposure[1])==F && is.null(dat$se.exposure[1])==F && is.na(dat$se.exposure[1])==F && is.null(dat$samplesize.exposure[1])==F && is.na(dat$samplesize.exposure[1])==F){
    R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
    F<- (dat$samplesize.exposure-2)*R2/(1-R2)
    dat$R2<-R2
    dat$F<-F
    dat<-subset(dat,F>F_value)
    return(dat)
  }
}


# MR
cause_getmr<-function(expo,outcome,LD_file,r2=0.001,
                      kb=10000,pval=1e-05,cl=NULL){
  format_cause_expo<-function(dat){
    library(cause)
    dat<-gwas_format(dat,snp='SNP',beta_hat ='beta.exposure',
                     se='se.exposure',A1="effect_allele.exposure",
                     A2="other_allele.exposure",chrom='chr.exposure',
                     pos="pos.exposure",p_value = 'pval.exposure'
    )
    return(dat)
  }
  format_cause_out<-function(dat){
    library(cause)
    dat<-gwas_format(dat,snp='SNP',beta_hat ='beta.outcome',
                     se='se.outcome',A1="effect_allele.outcome",
                     A2="other_allele.outcome",chrom='chr.outcome',
                     pos="pos.outcome",p_value = 'pval.outcome'
    )
    return(dat)
  }
  sample_cause<-function(dat,num_snp){
    set.seed(123)
    VAR<-with(dat,sample(snp,size=num_snp,replace=FALSE))
    return(VAR)
  }
  datap_cause<-function(dat){
    dat$p1<-pnorm(abs(dat$beta_hat_1/dat$seb1),lower.tail=F)*2
    return(dat)
  }
  dat_clump_cause<-function(dat){
    for_clump<-data.frame(dat$snp,dat$p1)
    colnames(for_clump)<-c('rsid','pval')
    return(for_clump)
  }
  datap_cause<-function(dat){
    dat$p1<-pnorm(abs(dat$beta_hat_1/dat$seb1),lower.tail=F)*2
    return(dat)
  }
  ld_local<-function(dat,r2,kb,p,LD_file){
    library(plinkbinr)
    
    plink_pathway<-get_plink_exe()
    
    dat<-ld_clump(dat,clump_r2 = r2,clump_kb = kb,
                  clump_p = p,
                  plink_bin =plink_pathway , bfile =LD_file)
    
    return(dat)
  }
  if("list"%in%class(outcome)){
    single_expo=T
    single_outcome=F}
  if("list"%in%class(expo)){
    single_outcome=T
    single_expo=F}
  if(single_expo==T){
    outcome<-pblapply(outcome,FUN=format_cause_out,cl=cl)
    expo<-format_cause_expo(expo)
  }
  if(single_outcome==T){
    expo<-pblapply(expo,FUN=format_cause_expo,cl=cl)
    outcome<-format_cause_out(outcome)
  }
  
  
  if(single_outcome==T){dat<-pblapply(expo,outcome,FUN=gwas_merge,cl=cl)}
  if(single_expo==T){dat<-pblapply(outcome,expo,
                                   FUN=function(outcome,exposure)gwas_merge(exposure,outcome),
                                   cl=cl)}
  
  
  VAR<-pblapply(dat,1000000,FUN=sample_cause)
  
  est<-foreach(i=1:length(VAR)) %do%{
    library(cause)
    data<-est_cause_params(dat[[i]],VAR[[i]])
  }
  dat<-pblapply(dat,FUN=datap_cause)
  for_clump<-lapply(dat, dat_clump_cause)
  clumped<-pblapply(for_clump,0.001,10000,1e-05,LD_file,FUN=ld_local)
  top_vars<-pblapply(clumped,FUN=function(dat) return(dat$rsid))
  cause_res<-list()
  for(i in 1:length(dat)){
    library(cause)
    res<-cause(dat[[i]],est[[i]],top_vars[[i]])
    cause_res[[i]]<-res
    print(i)
  }
  cause_table<-data.frame()
  for (i in 1:length(cause_res)){
    res<-cause_res[[i]]
    elpd<-res$elpd
    elpd$p<-pnorm(-elpd$z,lower.tail = F)
    name<-paste0('file',i)
    elpd$file<-name
    cause_table<-rbind(cause_table,elpd)
  }
  return(cause_table)
}

# RAPS
RAPS_getmr<-function(dat,dir_figure){
  setwd(dir_figure)
  res<-try(mr.raps(dat,over.dispersion = TRUE))
  if(class(res)%in%"try-error"){}
  else{
    exposure<-dat$id.exposure
    dir_create(exposure)
    dir_in<-paste0(dir_figure,'/',exposure)
    setwd(dir_in)
    #plot_name<-paste0(exposure,'raps.pdf')
    ggsave(file='raps_plot.pdf',plot=plot(res),width=9,height=5)
    res<-data.frame(beta.raps=res$beta.hat,se.raps=res$beta.se,
                    eov=res$tau2.hat,se.eov=res$tau2.se,OR.raps=NA,
                    or_lci95.raps=NA,or_uci95.raps=NA)
    res$OR.raps<- exp(res$beta.raps)
    res$or_lci95.raps<-exp(res$beta.raps)-(res$beta.raps*1.96)
    res$or_uci95.raps<-exp(res$beta.raps)+(res$beta.raps*1.96)
    res$pval.raps<-2*pnorm(abs(res$beta.raps/res$se.raps),lower.tail=F) 
    res$pval.eov<-2*pnorm(abs(res$eov/res$se.eov),lower.tail=F)
    if(is.na(res$pval.eov)==FALSE){
      if(res$pval.eov >0.05 ){
        res1<-mr.raps(dat,over.dispersion = F)
        setwd(dir_in)
        #plot_name<-paste0(exposure,'raps.pdf')
        ggsave(file='raps_plot.pdf',plot=plot(res1),width=9,height=5)
        res1<-data.frame(beta.raps=res1$beta.hat,se.raps=res1$beta.se,
                         eov=NA,se.eov=NA)
        res1$eov<-res$eov
        res1$se.eov<-res$se.eov
        res1$OR.raps<- exp(res1$beta.raps)
        res1$or_lci95.raps<-exp(res1$beta.raps)-(res1$beta.raps*1.96)
        res1$or_uci95.raps<-exp(res1$beta.raps)+(res1$beta.raps*1.96)
        res1$pval.raps<-2*pnorm(abs(res1$beta.raps/res1$se.raps),lower.tail=F) 
        res1$pval.eov<-2*pnorm(abs(res1$eov/res1$se.eov),lower.tail=F)
        
        res<-res1
      }
    }
    
    return(res)
  }
}

mr_dircreate_base<-function(root_dir,project_name,date=NULL){  
  library(fs)
  dir_name<-root_dir
  setwd(dir_name)
  
  if(is.null(date)==FALSE){data_u<-date}else{data_u<-Sys.Date()}
  
  dir_name2<-project_name
  dir_name3<-paste0(dir_name2,data_u)
  dir_create(dir_name3)
  setwd(paste0(dir_name,"/",dir_name3))
  dir_name4<-"1.figure"
  dir_name5<-"2.table"
  dir_name6<-"3.figure of sig res"
  dir_name7<-"4.snp with Fval"
  
  paste<-paste0(dir_name,"/",dir_name3,"/")
  
  dir1<-paste0(paste,dir_name4)
  dir.create(dir1)
  
  dir2<-paste0(paste,dir_name5)
  dir.create(dir2)
  
  dir3<-paste0(paste,dir_name6)
  dir.create(dir3)
  
  dir4<-paste0(paste,dir_name7)
  dir.create(dir4)
  
  res<-list(paste=paste,dir1=dir1,dir2=dir2,dir3=dir3,dir4=dir4)
  
  return(res)
}



# PRESSO
mr_Presso<-function(dat,num=10000){
  library(TwoSampleMR)
  library(MRPRESSO)
  library(dplyr)
  
  nsnp_filter=6
  set.seed(123)
  try (mr_presso_res<-mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat,  
                                SignifThreshold = 0.05,NbDistribution = num))
  return(mr_presso_res)
  
}
mr_presso_pval<-function(mr_presso_res){ 
  try ( mr_presso_main<-mr_presso_res$`Main MR results`)
  try ( mr_presso_main[3:5,]<-NA) 
  return(mr_presso_main)
}


mr_presso_snp<-function(mr_presso_res,mr_presso_main,dat,type="list"){
  data_re<-list()
  if(type=="list"){
    for(i in 1:length(mr_presso_res)){
      res<-mr_presso_res[[i]]
      main<-mr_presso_main[[i]]
      data<-dat[[i]]
      try(if(is.na(main[2,6])==FALSE){
        outliers<-which(res$`MR-PRESSO results`$`Outlier Test`$Pvalue<0.05)
        data$mr_keep[outliers]<-FALSE
      })
      data_re[[i]]<-data
      names(data_re)[[i]]<-names(dat)[[i]]
    }
    return(data_re)
  }
  
  if(type=="data"){
    res<-mr_presso_res$`MR-PRESSO results`
    main<-mr_presso_main
    data<-dat
    try(if(is.na(main[2,6])==FALSE){
      outliers<-which(res$`MR-PRESSO results`$`Outlier Test`$Pvalue<0.05)
      data$mr_keep[outliers]<-FALSE
    })
    return(data)
  }
}

#MRlap

mr_lap<-function(expo,outcome,ld,hm3,pval,r2,kb,MR_reverse=1e-03,save_logfiles=F){
  expo<-expo%>%select(rsid=SNP,chr=chr.exposure,pos=pos.exposure,alt=effect_allele.exposure,
                      ref=other_allele.exposure,N=samplesize.exposure,beta=beta.exposure,
                      se=se.exposure)
  expo<-as.data.frame(expo)
  outcome<-outcome%>%select(rsid=SNP,chr=chr.outcome,pos=pos.outcome,alt=effect_allele.outcome,
                            ref=other_allele.outcome,N=samplesize.outcome,beta=beta.outcome,
                            se=se.outcome)
  outcome<-as.data.frame(outcome)
  n_expo<-expo$exposure[1]
  n_out<-outcome$outcome[1]
  
  try(res<-MRlap::MRlap(exposure = expo,
                        exposure_name = n_expo,
                        outcome = outcome,
                        outcome_name = n_out,
                        ld = ld,
                        hm3 = hm3,MR_threshold=pval,MR_pruning_dist=kb,
                        MR_pruning_LD=r2,MR_reverse=MR_reverse,save_logfiles=save_logfiles
  ))
  try(snp<-data.frame(SNP_MRlap=res$MRcorrection$IVs))
  try(res_c<-as.data.frame(res$MRcorrection[-4])%>%select(nsnp=m_IVs,
                                                          beta_MRlap=corrected_effect,
                                                          se_MRlap=corrected_effect_se,
                                                          p_MRlap=corrected_effect_p))
  try(res<-c(res,res_c,snp))
  
  try(return(res))
}


# clean
clean_expo<-function(expo,pval,low_af=0.5,high_af=0.5,
                     clump=TRUE,kb=10000,r2=0.001,LD_file=NULL,af_filter=FALSE){
  library(TwoSampleMR)
  dat<-subset(expo,pval.exposure<pval)
  if(af_filter==TRUE){dat<-subset(dat,eaf.exposure<low_af|eaf.exposure>high_af)}
  
  if(clump==TRUE){
    
    if(is.null(LD_file)==TRUE){
      dat<-clump_data(dat,clump_kb = kb,clump_r2 = r2)}
    else{
      library(plinkbinr)
      plink_pathway<-get_plink_exe()
      snp<-data.frame(rsid=dat$SNP,pval=dat$pval.exposure)
      snp<-try(ld_clump(snp,clump_kb = kb,clump_r2 = r2,plink_bin =plink_pathway , bfile =LD_file))
      if ("try-error"%in% class(snp)){return(dat)}
      else{
        snp<-data.frame(SNP=snp$rsid)
        dat<-merge(snp,dat,by.x="SNP",by.y='SNP')
      }
    }
    
  }
  return(dat)
}

clean_list<-function(list,nrow=10){
  l<-lapply(list,nrow)
  n<-data.frame()
  for(i in 1:length(l)){
    if(is.null(l[[i]])==T || l[[i]]==0 )next
    
    n1<-data.frame(i,l[[i]])
    n<-rbind(n,n1)
  }
  colnames(n)<-c('l','row')
  n<-subset(n,row>nrow)
  list<-list[n$l]
  return(list)
}

clean_IV_from_outsig<-function(dat,MR_reverse=1e-03){
  dat<-subset(dat,pval.outcome>MR_reverse)
  return(dat)
}


LDSC_rg<-function(expo,outcome,an,sample_prev=NA,
                  population_prev=NA,ld,wld,chr_filter=c(1:22),n_blocks=200){
  id.o<-outcome$id.outcome[1]
  id.e<-expo$id.exposure[1]
  
  expo<-expo%>%mutate(Z=beta.exposure/se.exposure)
  expo<-expo%>%select(SNP=SNP,N=samplesize.exposure,Z=Z
                      ,A1=effect_allele.exposure
                      ,A2=other_allele.exposure)
  expo<-as_tibble(expo)
  
  outcome<-outcome%>%mutate(Z=beta.outcome/se.outcome)
  outcome<-outcome%>%select(SNP=SNP,N=samplesize.outcome,Z=Z
                            ,A1=effect_allele.outcome
                            ,A2=other_allele.outcome)
  outcome<-as_tibble(outcome)
  
  
  dat<-list(expo,outcome)
  names(dat)<-c(id.e,id.o)
  
  rm(expo,outcome)
  
  
  res<-try(ldscr::ldsc_rg(dat,ancestry = an,sample_prev=sample_prev,
                          population_prev=population_prev,ld=ld,wld=wld,
                          n_blocks=n_blocks,chr_filter=chr_filter))
  
  return(res)
  
}


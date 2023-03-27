library(mr.raps)
library(TwoSampleMR)
#library(plyr)
library(dplyr)
#library(fs) 
#library(ggplot2) 
#library(lubridate) 
library(ieugwasr)
library(plinkbinr) 
library(tidyr)
#library(progress)
library(data.table)
library(MRPRESSO)
library(parallel)
#library(foreach)
#library(doParallel)
library(pbapply)
library(stringr)
library(cause)
library(vroom)
library(MungeSumstats)
#library(GenomicFiles)
#library(meta)
#library(readr)
#library(readxl)
#library(forestploter)
library(ldscr)
library(MRlap)
#library(meta)
#library(forestplot)
#library(biomaRt)
#library(gwasglue)
#library(coloc)

# format
format_Mun<-function(file,source="finn_r8",save_path=NULL,lift=F,ref_genome = "hg38",
                     convert_ref_genome = "hg19"){
  library(data.table)
  library(MungeSumstats)
  library(dplyr)
  dat<-fread(file)
  if(source=="finn_r8"){
    dat<-dat%>%dplyr::select(SNP=rsids,CHR=`#chrom`,BP=pos,A1=ref,
                             A2=alt,FRQ=af_alt,BETA=beta,
                             SE=sebeta,P=pval)
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
  
  
  
  if(source=="bac_eaf"){
    data<-format_data(data,type=type,id_col="bac",
                      phenotype_col = "phenotype",snp_col = "snp",
                      effect_allele_col = "effect_allele",other_allele_col = "other_allele",
                      eaf_col ="MAF",beta_col="beta",se_col = "se",
                      pval_col = "pval",chr_col = "chr",pos_col = "bp",samplesize_col = "samplesize")
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
  
  if(source=="drug"){
    data<-format_data(data,type=type,snp_col = "SNP",
                      effect_allele_col = "A1",other_allele_col = "A2",
                      beta_col="b",se_col = "se",
                      pval_col = "p",samplesize_col = "n",eaf_col = "freq")
    
  }
  
  
  return(data)
}

format_cause<-function(dat,type="exposure"){
  library(cause)
  
  if(type=='exposure'){
    dat<-gwas_format(dat,snp='SNP',beta_hat ='beta.exposure',
                     se='se.exposure',A1="effect_allele.exposure",
                     A2="other_allele.exposure",chrom='chr.exposure',
                     pos="pos.exposure",p_value = 'pval.exposure'
    )
  }
  
  if(type=='outcome'){
    dat<-gwas_format(dat,snp='SNP',beta_hat ='beta.outcome',
                     se='se.outcome',A1="effect_allele.outcome",
                     A2="other_allele.outcome",chrom='chr.outcome',
                     pos="pos.outcome",p_value = 'pval.outcome'
    )
  }
  
  
  return(dat)
}



# read
read_vcf_cyclemr<-function(file_name,nThread = 8,type=".gz"){
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

# MR
# cause
cause_estimate<-function(dat,num_snp=1000000){
  library(cause)
  cause_sample<-function(dat,num_snp){
    set.seed(123)
    VAR<-with(dat,sample(snp,size=num_snp,replace=FALSE))
    return(VAR)
  }
  VAR<-pblapply(dat,num_snp,FUN=cause_sample)
  
  
  est<-foreach(i=1:length(VAR)) %do%{
    library(cause)
    data<-est_cause_params(dat[[i]],VAR[[i]])
  }
  return(est)
}

cause_add_pval<-function(dat){
  
  cause_pval<-function(dat){
    dat$p1<-pnorm(abs(dat$beta_hat_1/dat$seb1),lower.tail=F)*2
    return(dat)
  }
  dat<-pblapply(dat,FUN=cause_pval)
  return(dat)
}


cause_snp<-function(dat){
  
  get_snp<-function(dat) {
    dat<-dat$rsid
    return(dat)
  }
  
  snp<-lapply(dat,fun=get_snp)
  
  return(snp)
}

mr_cause<-function(dat,estimate,top_variants){
  cause_res<-list()
  for(i in 1:length(dat)){
    library(cause)
    res<-cause(dat[[i]],est[[i]],top_vars[[i]])
    cause_res[[i]]<-res
    print(i)
  }
  return(cause_res)
}
# RAPS
RAPS<-function(dat,dir_figure){
  setwd(dir_figure)
  res<-mr.raps(dat,over.dispersion = TRUE)
  exposure<-paste0(dat$exposure,dat$id.exposure)
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
  
  return(res)
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

# bind

bind_base_raps<-function(table,raps_res,name,save_dir,csv_filename){
  id<-name$V1[which(name$V2==raps_res$exposure)]
  table$beta.raps<-NA
  table$se.raps<-NA
  table$eov<-NA
  table$se.eov<-NA
  table$OR.raps<-NA
  table$or_lci95.raps<-NA
  table$or_uci95.raps<-NA
  table$pval.raps<-NA
  table$pval.eov<-NA
  for(i in 1:nrow(table)){
    if(i%%5==1){
      row<-which(id==table$id.exposure[i])[1]
      table$beta.raps[i] <-raps_res$beta.raps[row]
      table$se.raps[i] <-raps_res$se.raps[row]
      table$eov[i] <-raps_res$eov[row]
      table$se.eov[i] <-raps_res$se.eov[row]
      table$OR.raps[i] <-raps_res$OR.raps[row]
      table$or_lci95.raps[i] <-raps_res$or_lci95.raps[row]
      table$or_uci95.raps[i] <-raps_res$or_uci95.raps[row]
      table$pval.raps[i] <-raps_res$pval.raps[row]
      table$pval.eov[i] <-raps_res$pval.eov[row]
    }
  }
  setwd(save_dir)
  write.csv(table,csv_filename)
  
  return(table)
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

clean_outcome_from_exposure<-function(expo,outcome){
  snp<-lapply(expo,FUN=function(x)data.frame(SNP=x$SNP))
  snp<-bind_rows(snp)
  snp<-data.frame(SNP=unique(snp$SNP))
  outcome_adj<-merge(outcome,snp)
  return(outcome_adj)
}

clean_IV_from_outsig<-function(dat,MR_reverse=1e-03){
  dat<-subset(dat,pval.outcome>MR_reverse)
  return(dat)
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


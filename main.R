#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
test<-read.table(args[1],header=FALSE)
colnames(test)<-c("SAMPLE_ID","SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8","SNP9","SNP10")
print("#1. loading libraries for statistic summary and plot library(pastecs)")
library(Hmisc); library(ggplot2)

print("#2. loading the data for IOP / PCA")
table_all<-read.table(args[2],header=TRUE)
table_all<-merge(table_all,test,by.x="SAMPLE_ID",by.y="SAMPLE_ID")

print("#3 calculate the permulated mean of ansari.test p-value")

#3.1 set up list for pvalue

para_hm_pvalue_list<-c()
para_ht_pvalue_list<-c()
nonpara_hm_pvalue_list<-c()
nonpara_ht_pvalue_list<-c()

for (i0 in 3:snp_i+2)
{
i1<-i0+33
print(paste(i0," ",i1))
wt=sum(table_all[,i1]==0)
ht=sum(table_all[,i1]==1)
hm=sum(table_all[,i1]==2)
print(paste("wt:",wt," ht:",ht," hm:",hm))
if (wt>100 & ht>100 & hm>100) {
Parametric_hm<-var.test(table_all[table_all[,i1]==2,23],table_all[table_all[,i1]==0,23])
Parametric_ht<-var.test(table_all[table_all[,i1]==1,23],table_all[table_all[,i1]==0,23])
Non_Parametric_hm<-ansari.test(table_all[table_all[,i1]==2,23],table_all[table_all[,i1]==0,23])
Non_Parametric_ht<-ansari.test(table_all[table_all[,i1]==1,23],table_all[table_all[,i1]==0,23])
print(paste(Parametric_hm$p.value,":",Parametric_ht$p.value,":",Non_Parametric_hm$p.value,":",Non_Parametric_ht$p.value))
print("#3.3 storing ansari.test p.value into the list")
para_hm_pvalue_list[i0]<-Parametric_hm$p.value
para_ht_pvalue_list[i0]<-Parametric_ht$p.value
nonpara_hm_pvalue_list[i0]<-Non_Parametric_hm$p.value
nonpara_ht_pvalue_list[i0]<-Non_Parametric_ht$p.value
print(nonpara_hm_pvalue_list[i0])
print(i0)
} else{
writeLines(paste(i0," error!!! ","wt:",wt," ht:",ht," hm:",hm,sep=""),"/data5/playyard/taozhang/ukbiobank/genotype_data/ansari.test.error")
print(paste(i0," error!!! ","wt:",wt," ht:",ht," hm:",hm,sep=""))
nonpara_hm_pvalue_list[i0]<-"NA"
nonpara_ht_pvalue_list[i0]<-"NA"
}
}
#4 plot the pvalue distribution
print (para_hm_pvalue_list)
print (para_ht_pvalue_list)
print (nonpara_hm_pvalue_list)
print (nonpara_ht_pvalue_list)
para_hm_pvalue_list<-data.frame(para_hm_pvalue_list,stringsAsFactors =FALSE)
para_ht_pvalue_list<-data.frame(para_ht_pvalue_list,stringsAsFactors =FALSE)
write.table(para_hm_pvalue_list,file=paste(input,".Ftest.hm.pvaluelist",sep=""))
write.table(para_ht_pvalue_list,file=paste(input,".Ftest.ht.pvaluelist",sep=""))
nonpara_hm_pvalue_list<-data.frame(nonpara_hm_pvalue_list,stringsAsFactors =FALSE)
nonpara_ht_pvalue_list<-data.frame(nonpara_ht_pvalue_list,stringsAsFactors =FALSE)
write.table(nonpara_hm_pvalue_list,file=paste(input,".ansari.hm.pvaluelist",sep=""))
write.table(nonpara_ht_pvalue_list,file=paste(input,".ansari.ht.pvaluelist",sep=""))
write.table(table_all,file=paste(input,".combined",sep=""))

#p<-ggplot(data= nonpara_pvalue_list,aes(x= nonpara_pvalue_list))+geom_histogram()
#ggsave(file=paste("trek1.",i0, ".png",sep=""), plot = p, device = "png", path = argv[3],scale = 1, width = NA, height = NA,   dpi = 300, limitsize = TRUE)

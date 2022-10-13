setwd("/scratch/weissvl/shengq2/20210616_human_exomeseq/bam_coverage/result")

library(reshape2)
library(ggplot2)

coverage_files=list.files(".", pattern=".coverage.txt")

f=coverage_files[1]
covs=NULL
for(f in coverage_files){
  cov=read.table(f)$V1[1]
  name=gsub(".coverage.txt", "", f)
  covs=rbind(covs, data.frame("Name"=name, "Coverage"=cov))
}

reads<-read.table("/scratch/weissvl/shengq2/20210616_human_exomeseq/fastqc_raw_summary/result/human_exomeseq.FastQC.reads.tsv", sep="\t", header=T)
reads<-unique(reads[,c("Sample", "Reads")])
reads$Depth<-reads$Reads / 39024143.0 * 150
reads$Coverage<-covs$Coverage
reads<-reads[,-2]

write.csv(reads, file="depth_coverage.csv", row.names=F)

mreads<-melt(reads, id.var="Sample", variable.name="Category")

library(ggplot2)
png("depth_coverage.png", width=2000, height=2000, res=300)
g<-ggplot(mreads,aes(x=Category,y=value)) + geom_violin() + geom_jitter(width=0.2) + theme_bw() + xlab("") + ylab("")
print(g)
dev.off()

reads<-read.csv("depth_coverage.csv", row.names=1, header=T)
mean_v<-apply(reads, 2, mean)
mean_v

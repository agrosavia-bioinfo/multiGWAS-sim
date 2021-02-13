#!/usr/bin/Rscript
library (dplyr)
library (ggplot2)

args = commandArgs (trailingOnly=T)
args = c ("h50")
inDir = args [1]

h2Dirs = list.dirs (inDir, recursive=F)

SNPsType = c ("Significant_SNPs", "Top_SNPs")
for (type in SNPsType) {
	h2DataAll =  NULL
	for (h in h2Dirs) {
		sharedFile = sprintf ("%s/out-statistics-SHARED-%s-table.csv", h, type)
		sharedData = read.csv (sharedFile)
		h2Data     = data.frame (H2=basename (h), sharedData)
		h2DataAll  = if (is.null (h2DataAll)) h2Data else rbind (h2DataAll, h2Data)
	}
	write.csv (h2DataAll, sprintf ("out-statistics-heritabilities-%s-table.csv", type), quote=F, row.names=F)
}

dt    = read.csv ("out-statistics-heritabilities-Top_SNPs-table.csv")
dt1   = filter (dt, SHARED=="Shared2", TOOL=="MultiGWAS")[,-2]
sum1  = dt1 %>% group_by (H2,TOOL) %>% summarize (TPR=median (TPR), FPR=median (FPR))
sum1s = sum1 %>% arrange (TOOL)

plot (sum1s$FPR, sum1s$TPR)
ggplot(sum1s,aes(FPR,TPR,color=TOOL))+geom_line(size = 2, alpha = 0.7)+
      labs(title= "ROC curve", 
           x = "False Positive Rate (1-Specificity)", 
           y = "True Positive Rate (Sensitivity)")

sum1s
dim (dt1)

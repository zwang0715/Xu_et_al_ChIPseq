if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DiffBind")
install.packages("tidyverse")
BiocManager::install("profileplyr")
BiocManager::install("Repitools")

library('Repitools')


library(DiffBind)
library(tidyverse)
library(profileplyr)

setwd("~/diffbind/")

samples_AAR <- read.csv('AAR_sample_sheet.csv')

dbObj_AAR <- dba(sampleSheet=samples_AAR)

dbObj_AAR

dbObj_AAR <- dba.count(dbObj_AAR, bUseSummarizeOverlaps=TRUE)

dbObj_AAR

dba.plotPCA(dbObj_AAR,  attributes=DBA_FACTOR, label=DBA_ID)

plot(dbObj_AAR)

dbObj_AAR <- dba.normalize(dbObj_AAR)

dbObj_AAR <- dba.contrast(dbObj_AAR, categories=DBA_CONDITION, minMembers = 2)

dbObj_AAR <- dba.analyze(dbObj_AAR, method=DBA_ALL_METHODS)

dba.show(dbObj_AAR, bContrasts=T)	

dba.plotPCA(dbObj_AAR, contrast=1, method=DBA_EDGER, attributes=DBA_CONDITION, label=DBA_ID)

dba.plotVenn(dbObj_AAR,contrast=1,method=DBA_ALL_METHODS)

dba.plotMA(dbObj_AAR, method=DBA_EDGER)

dba.plotMA(dbObj_AAR, bXY=TRUE, method = DBA_EDGER)

dba.plotVolcano(dbObj_AAR, contrast=1, method=DBA_EDGER)

plot <- dba.plotVolcano(dbObj_AAR, contrast=1, method=DBA_EDGER)

AAR_report <- dba.report(dbObj_AAR, method=DBA_EDGER, contrast = 1, th=1)

#AAR_report <- dba.report(dbObj_AAR, method=DBA_EDGER, contrast = 1, th=0.05)

df <- annoGR2DF(AAR_report)

#df <- annoGR2DF(plot)

write.table(df, file = "./AAR_diffbind_result.txt", sep = "\t", row.names = TRUE)

profiles <- dba.plotProfile(dbObj_AAR, sites = AAR_report)
dba.plotProfile(profiles, maxSites =2000)

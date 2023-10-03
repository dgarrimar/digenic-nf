#!/usr/bin/env Rscript

library(data.table)
library(bigtabulate)

# Load data
snp1 = fread("one.txt", h = F, data.table=F)[,1]
snp2 = fread("two.txt", h = F, data.table=F)[,1]
snp0 = snp1[!is.na(snp1) & !is.na(snp2)]
snp2 = snp2[!is.na(snp1) & !is.na(snp2)]
snp1 = snp0
n = length(snp1)

# Recode genotypes and obtain freq tables
nms1 = nms2 = as.character(c(0,1,2))

tbl1 = tabulate(snp1+1, 3)
names(tbl1) = nms1
keep1 = tbl1>0
tbl1 = tbl1[keep1]
nms1 = nms1[keep1]
minor1 = nms1[which.min(tbl1)]
major1 = nms1[!nms1 %in% c(1, minor1)]
minor_het1 = unique(c(minor1, "1"))

tbl2 = tabulate(snp2+1, 3)
names(tbl2) = nms2
keep2 = tbl2>0
tbl2 = tbl2[keep2]
nms2 = nms2[keep2]
minor2 = nms2[which.min(tbl2)]
major2 = nms2[!nms2 %in% c(1, minor2)]
minor_het2 = unique(c(minor2, "1"))

tbl1_rec = c(sum(tbl1[minor_het1]), tbl1[major1]) 
tbl2_rec = c(sum(tbl2[minor_het2]), tbl2[major2]) 

tblJ = bigtabulate(cbind(snp1, snp2), 1:2)

# Observed freqs
O1 = sum(tblJ[minor_het1, minor_het2])
O2 = sum(tblJ[minor_het1, major2])
O3 = sum(tblJ[major1, minor_het2])
O4 = sum(tblJ[major1, major2])

# Expected freqs
E1 = tbl1_rec[1]/n*tbl2_rec[1]
E2 = tbl1_rec[1]/n*tbl2_rec[2]
E3 = tbl1_rec[2]/n*tbl2_rec[1]
E4 = tbl1_rec[2]/n*tbl2_rec[2]

# Goodness-of-fit Chi2 test
c2t = chisq.test(c(O1, O2, O3, O4), p = c(E1, E2, E3, E4)/n)

cat(paste(O1,O2,O3,O4),paste(E1,E2,E3,E4),c2t$statistic,c2t$p.value)

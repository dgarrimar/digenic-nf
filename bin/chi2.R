#!/usr/bin/env Rscript

library(data.table)

# Load data
snp1 = fread("one.txt", h = F, data.table=F)[,1]
snp2 = fread("two.txt", h = F, data.table=F)[,1]
snp0 = snp1[!is.na(snp1) & !is.na(snp2)]
snp2 = snp2[!is.na(snp1) & !is.na(snp2)]
snp1 = snp0
n = length(snp1)

# Recode genotypes and obtain freq tables
nms = c(0,1,2)
recode = function(snp, nms = c(0, 1, 2)){
  tbl = tabulate(snp+1, 3)
  names(tbl) = nms
  keep = tbl>0
  tbl = tbl[keep]
  nms = nms[keep]
  minor = as.character(nms[which.min(tbl)])
  major = as.character(nms[!nms %in% c(1, minor)])
  snp = as.character(snp)
  snp = .Internal(gsub(minor, "m", snp, F, F, T, F))
  snp = .Internal(gsub("1", "m", snp, F, F, T, F))
  snp = .Internal(gsub(major, "M", snp, F, F, T, F))
  return(as.factor(snp))
}

snp1 = recode(snp1)
snp2 = recode(snp2)

tbl1 = tabulate(snp1)
names(tbl1) = levels(snp1)

tbl2 = tabulate(snp2)
names(tbl2) = levels(snp2)

# Observed freqs
tblJ = table(snp1, snp2)
O1 = tblJ[1,1]
O2 = tblJ[1,2]
O3 = tblJ[2,1]
O4 = tblJ[2,2]

# Expected freqs
E1 = tbl1[1]/n*tbl2[1]
E2 = tbl1[1]/n*tbl2[2]
E3 = tbl1[2]/n*tbl2[1]
E4 = tbl1[2]/n*tbl2[2]

# Goodness-of-fit Chi2 test
c2t = chisq.test(c(O1, O2, O3, O4), p = c(E1, E2, E3, E4)/n)

cat(paste(O1,O2,O3,O4),paste(E1,E2,E3,E4),c2t$statistic,c2t$p.value)

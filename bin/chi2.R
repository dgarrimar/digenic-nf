#!/usr/bin/env Rscript

library(data.table)

# Load data
snp1 = as.numeric(fread("one.txt", h = F))
snp2 = as.numeric(fread("two.txt", h = F))
snp0 = snp1[!is.na(snp1) & !is.na(snp2)]
snp2 = snp2[!is.na(snp1) & !is.na(snp2)]
snp1 = snp0
n = length(snp1)

# Recode genotypes and obtain freq tables
recode <- function(snp){
    tbl = table(snp)
    minor = names(which.min(tbl))
    major = names(tbl)[!names(tbl) %in% c(1, minor)]
    snp = gsub(minor, "m", snp)
    snp = gsub(1, "m", snp)
    snp = gsub(major, "M", snp)
    return(snp)
}

snp1 = recode(snp1)
snp2 = recode(snp2)
tbl1 = table(snp1)
tbl2 = table(snp2)

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

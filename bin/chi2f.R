#!/usr/bin/env Rscript

library(data.table)
library(bigtabulate)

# Load chunk
args = commandArgs(trailingOnly=TRUE)
pairs = fread(args[1], h = F, data.table = F)

# Iterator
gt.f = args[2]
one_bckp = ''
for (i in 1:nrow(pairs)) {
    
    start=Sys.time()
    # Tabix, load data
    pair = pairs[i, ]
    one = pair[, 1]
    two = pair[, 2]
    if (one == one_bckp) {
        system(sprintf("tabix %s %s | cut -f4- | sed 's/\\t/\\n/g' > two.txt", gt.f, two))
        src2 = fread("two.txt", h = F, data.table = F)[,1]
    } else {
        one_bckp = one
        system(sprintf("tabix %s %s | cut -f4- | sed 's/\\t/\\n/g' > one.txt && tabix %s %s | cut -f4- | sed 's/\\t/\\n/g' > two.txt", gt.f, one, gt.f, two))
        src1 = fread("one.txt", h = F, data.table = F)[,1]
        src2 = fread("two.txt", h = F, data.table = F)[,1]
    }
    snps = cbind(src1, src2)
    snps = snps[complete.cases(snps), ]
    snp1 = snps[,1]
    snp2 = snps[,2]
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
    
    tblJ = bigtabulate(snps, 1:2)
    
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
    
    end=Sys.time()
    cat(one, two, paste(O1,O2,O3,O4), paste(E1,E2,E3,E4), c2t$statistic, c2t$p.value, format(end-start), "\n")
}


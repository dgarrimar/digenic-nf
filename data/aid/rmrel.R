library(data.table)
library(dplyr)
library(tidyr)
library(ukbtools)

setwd("/users/rg/dgarrido/projects/digenic/digenic-nf/data/aid")
king <- fread("king.txt", data.table = F, h = T)
king <- king[king$Kinship>0, ] # remove weird entries
king$ID1 <- as.character(king$ID1)
king$ID2 <- as.character(king$ID2)

dataids <- as.character(fread("keep.txt", data.table = F, h = F)[,1])

write.table(ukb_gen_samples_to_remove(king, ukb_with_data = dataids, cutoff = 0.0442), # 0.0442 up to 3rd degree relationships
            "related2rm.txt", quote = F, row.names = F, col.names =F)

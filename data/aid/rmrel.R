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

cff <- 0.0884 # 0.0442 up to 3rd degree relationships

write.table(ukb_gen_samples_to_remove(king, ukb_with_data = dataids, cutoff = cff), 
            sprintf("related2rm_%s",cff), quote = F, row.names = F, col.names =F)

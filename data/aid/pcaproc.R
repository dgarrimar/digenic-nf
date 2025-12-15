library(data.table)
library(ggplot2)
library(RColorBrewer)
library(microbenchmark)

setwd("/users/rg/dgarrido/projects/digenic/digenic-nf/data/aid")

# Read & cleanup
pcs <- fread("pcs.txt", h = F, data.table = F)
#all(is.na(pcs$V1))
#all(is.na(pcs$V1))
pcs <- pcs[, -c(1,43)]
napc1 <- is.na(pcs[,2])
pcs <- pcs[!napc1, ]
rownames(pcs) <- pcs[,1]
pcs[,1] <- NULL

# Get first 2 PCs
pcs12 <- pcs[, 1:2]
colnames(pcs12) <- c("PC1", "PC2")
# dim(pcs12)

# K-means (3 groups)
set.seed(1)
c <- kmeans(pcs12, 3)
pcs12$c <- c$cluster
centers <- c$centers
# table(pcs12$c)

# Plot PCA
# set.seed(1)
# rndmpick <- sample(rownames(pcs12), 10000)
# dummy1 <- pcs12[rndmpick, ]
p <- ggplot(pcs12, aes(x = PC1, y = PC2, col = as.factor(c))) +
  geom_point(alpha = 0.1) +
  geom_rug(alpha = 0.01) +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw()

D1 <- data.frame(d = sqrt(rowSums(sweep(pcs12[pcs12$c==1,1:2], 2, centers[1,])^2)))    # Distances to centroid (Europeans)

df <- merge(pcs12, D1, by = "row.names")
rownames(df) <- df$Row.names
df$Row.names <- NULL
df$c <- NULL
th <- quantile(df$d, 0.99)  # 1% outliers masked
df$f <- df$d < th
# dim(df)

# set.seed(1)
# rndmpick <- sample(rownames(df), 10000)
# dummy2 <- df[rndmpick, ]
q <- ggplot(df, aes(x = PC1, y = PC2, col = as.factor(f))) +
  geom_point(alpha = 0.1) +
  geom_rug(alpha = 0.01) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw(base_size = 20)

fwrite(data.frame(rownames(df)[df$f]), file = "european.txt", quote = F, sep = "\t", col.names = F, row.names = F) 

# Self reported
sra <- read.delim("ancestry_sr.txt", h = F)
colnames(sra) <- c("ID", "SRA")
rownames(sra) <- sra$ID
sra$ID <- NULL
df2 <- merge(df,sra, by = "row.names")
# dim(df2)

tbl <- table(df2[,c("SRA", "f")])
# tbl
# pheatmap::pheatmap(log10(tbl+1), display_numbers = tbl)
# pheatmap::pheatmap(round(tbl/rowSums(tbl),2), display_numbers = T)

# p
# q

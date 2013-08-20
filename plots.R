# Constants
bar_colors <- c("#333333", "#E6E6E6")

# Length distribution
png("lendis.png", family="Times", pointsize=14)
lendis <- read.table("all.join.cull.lenDis", skip=2, header=TRUE)
lendis.m <- as.matrix(lendis[5:6])
rownames(lendis.m) <- lendis[,1]
par(mar=c(5, 7, 4, 2) + 0.1)
barplot(t(lendis.m), beside=TRUE, col=bar_colors, space=c(0, 0.5), axes=FALSE, xlab="Length", las=2)
axis(2, yaxp=c(0, 22500000, 10), las=2)
title(ylab="Bases (bp)", line=5.5)
legend("topright", c("inter", "intra"), fill=bar_colors)
dev.off()

# 2 kbp length distribution
png("lendis_2K.png", family="Times", pointsize=14)
lendis.2K <- read.table("all.join.cull.lenDis_2K", skip=2, header=TRUE)
lendis.2K.m <- as.matrix(lendis.2K[5:6])
rownames(lendis.2K.m) <- lendis.2K[,1]
par(mar=c(7, 7, 4, 2) + 0.1)
barplot(t(lendis.2K.m), beside=TRUE, col=bar_colors, space=c(0, 0.5), axes=FALSE, las=2)
axis(2, yaxp=c(0, 22500000, 10), las=2)
title(xlab="Length (kbp)", line=4)
title(ylab="Bases (bp)", line=5.5)
legend("topright", c("inter", "intra"), fill=bar_colors)
dev.off()

# Aligned bases by % identity
png("simdis.png", family="Times", pointsize=14)
simdis <- read.table("all.join.cull.simDis", skip=2, header=TRUE)
simdis.m <- as.matrix(simdis[4:5])
par(mar=c(7, 7, 4, 2) + 0.1)
barplot(t(simdis.m), beside=TRUE, col=bar_colors, space=c(0, 0.5), axes=FALSE, las=2)
axis(2, yaxp=c(0, 5150000, 10), las=2)
title(ylab="Aligned Bases (bp)", line=5.5)
title(xlab="Identity (%)", line=4)
legend("topright", c("inter", "intra"), fill=bar_colors)
dev.off()

# Aligned bases by % identity (smaller interval)
png("simdis_05.png", family="Times", pointsize=14)
simdis <- read.table("all.join.cull.simDis_05", skip=2, header=TRUE)
simdis.m <- as.matrix(simdis[5:6])
rownames(simdis.m) <- simdis[,1]
par(mar=c(7, 7, 4, 2) + 0.1)
barplot(t(simdis.m), beside=TRUE, col=bar_colors, space=c(0, 0.5), axes=FALSE, las=2)
axis(2, yaxp=c(0, 3000000, 10), las=2)
title(ylab="Aligned Bases (bp)", line=5.5)
title(xlab="Identity (%)", line=4)
legend("topright", c("inter", "intra"), fill=bar_colors)
dev.off()

# Similarity and Kimura's k by length
sim_kimura <- read.table("length_similarity_kimura", header=TRUE)
png("similarity_by_length.png", family="Times", pointsize=14)
plot(sim_kimura$base_S, sim_kimura$per_sim, xlab="Length (bp)", ylab="Similarity (%)")
dev.off()

png("kimura_by_length.png", family="Times", pointsize=14)
plot(sim_kimura$base_S, sim_kimura$k_kimura, xlab="Length (bp)", ylab="Kimura's k")
dev.off()

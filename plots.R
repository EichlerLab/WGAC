# Load libraries.
library(ggplot2)
library(grid)
options( echo = TRUE )

# Constants
bar_colors <- c("#F8766D", "#00BFC4")
axis_labels <- element_text(family="Times New Roman", size=18)
x_axis_ticks <- element_text(family="Times New Roman", size=18, angle=45, hjust=1, vjust = 1)


library(RColorBrewer)
colors <- brewer.pal(n = 3, name = 'Set1')[1:2]
# colors[1] is red
# colors[2] is blue

# Length distribution
png("lendis.png", width=700)
lendis <- read.table("all.join.cull.lenDis", skip=2, header=TRUE)
sizes <- factor(as.character(lendis$size), levels=unique(lendis$size), ordered=TRUE)
lengths <- rbind(
    data.frame(size=sizes, length=lendis$interlen, type="inter"),
    data.frame(size=sizes, length=lendis$intralen, type="intra")
)
# Convert bases to Megabases by dividing lengths by 1,000,000.
ggplot(lengths, aes(y=length/1000000, x=size, fill=type)) +
geom_bar(position="dodge", stat = "identity") + scale_y_continuous("Duplicated Bases (Mbp)") + xlab("Duplication Length") + theme_bw() + theme(axis.title.x=axis_labels, axis.title.y=axis_labels, axis.text.x=axis_labels, axis.text.y=axis_labels, legend.text=axis_labels, legend.title=axis_labels) +
scale_fill_manual("legend", values = c("intra" = colors[2], "inter" = colors[1] ) )

dev.off()

# 2 kbp length distribution
png("lendis_2K.png", width=700)
lendis.2K <- read.table("all.join.cull.lenDis_2K", skip=2, header=TRUE)
sizes <- factor(as.character(lendis.2K$size), levels=unique(lendis.2K$size), ordered=TRUE)
lengths <- rbind(
    data.frame(size=sizes, length=lendis.2K$interlen, type="inter"),
    data.frame(size=sizes, length=lendis.2K$intralen, type="intra")
)
# Convert bases to Megabases by dividing lengths by 1,000,000.
ggplot(lengths, aes(y=length/1000000, x=size, fill=type)) + geom_bar(position="dodge", stat = "identity") + scale_y_continuous("Duplicated Bases (Mbp)") + xlab("Duplication Length") + theme_bw() + theme(axis.title.x=axis_labels, axis.title.y=axis_labels, axis.text.x=x_axis_ticks, axis.text.y=axis_labels, legend.text=axis_labels, legend.title=axis_labels) +
scale_fill_manual("legend", values = c("intra" = colors[2], "inter" = colors[1] ) )

dev.off()

# Aligned bases by % identity
png("simdis.png", width=700)
simdis <- read.table("all.join.cull.simDis", skip=2, header=TRUE)
identities <- factor(as.character(rownames(simdis)), levels=unique(rownames(simdis)), ordered=TRUE)
lengths <- rbind(
    data.frame(size=identities, length=simdis$TotalinterLength, type="inter"),
    data.frame(size=identities, length=simdis$intraLength, type="intra")
)
# Convert bases to Megabases by dividing lengths by 1,000,000.
ggplot(lengths, aes(y=length/1000000, x=size, fill=type)) + geom_bar(position="dodge", stat = "identity") + scale_y_continuous("Aligned Bases (Mbp)") + xlab("Identity") + theme_bw() + theme(axis.title.x=axis_labels, axis.title.y=axis_labels, axis.text.x=axis_labels, axis.text.y=axis_labels, legend.text=axis_labels, legend.title=axis_labels) +
scale_fill_manual("legend", values = c("intra" = colors[2], "inter" = colors[1] ) )

dev.off()

# Aligned bases by % identity (smaller interval)
png("simdis_05.png", width=700)
simdis <- read.table("all.join.cull.simDis_05", skip=2, header=TRUE)
identities <- factor(as.character(simdis$size), levels=unique(simdis$size), ordered=TRUE)
lengths <- rbind(
    data.frame(size=identities, length=simdis$interLength, type="inter"),
    data.frame(size=identities, length=simdis$intraLength, type="intra")
)
# Convert bases to Megabases by dividing lengths by 1,000,000.
ggplot(lengths, aes(y=length/1000000, x=size, fill=type)) + geom_bar(position="dodge", stat = "identity") + scale_y_continuous("Aligned Bases (Mbp)") + xlab("Identity") + theme_bw() + theme(axis.title.x=axis_labels, axis.title.y=axis_labels, axis.text.x=x_axis_ticks, axis.text.y=axis_labels, legend.text=axis_labels, legend.title=axis_labels) +
scale_fill_manual("legend", values = c("intra" = colors[2], "inter" = colors[1] ) )

dev.off()

# Similarity and Kimura's k by length
png("kimura_by_length.png", width=700)
sim_kimura <- read.table("length_similarity_kimura", header=TRUE)
ggplot(sim_kimura, aes(y=k_kimura, x=base_S)) + geom_point(shape=21) + ylab("Kimura's k") + xlab("Length (bp)") + theme_bw() + theme(axis.title.x=axis_labels, axis.title.y=axis_labels, axis.text.x=axis_labels, axis.text.y=axis_labels, legend.text=axis_labels, legend.title=axis_labels, plot.margin=unit(c(1,1,1,1), "cm"))
dev.off()

# Non-redundant duplication by chromosome/contig.
bar_colors <- c("#F8766D", "#00BFC4", "beige")
png("nonredundant_duplication_by_chromosome.png", width=700)
nr <- read.table("nr_stats.tab", header=TRUE)
chromosomes <- factor(as.character(nr$contig), levels=unique(nr$contig), ordered=TRUE)
lengths <- rbind(
    data.frame(chromosome=chromosomes, length=nr$inter, type="inter"),
    data.frame(chromosome=chromosomes, length=nr$intra, type="intra"),
    data.frame(chromosome=chromosomes, length=nr$total, type="both")
)
# Convert bases to Megabases by dividing Kbp lengths by 1,000.
ggplot(lengths, aes(y=length/1000, x=chromosome, fill=type)) + geom_bar(colour="black", position="dodge", stat = "identity") + scale_y_continuous("Bases (Mbp)") + xlab("Chromosome") + theme_bw() + theme(axis.title.x=axis_labels, axis.title.y=axis_labels, axis.text.x=x_axis_ticks, axis.text.y=axis_labels, legend.text=axis_labels, legend.title=axis_labels) + scale_fill_manual(values=bar_colors)
dev.off()

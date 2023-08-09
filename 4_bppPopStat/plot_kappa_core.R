Species_is="Zymoseptoria_RUN2"
Abb_is="Zymoseptoria"
Abb_slurm="ZT"
Num_chr=485

# plot alignment rate
library(ggplot2)
library(ggforce)
library(dplyr)

# load data
run_loose_kappa <- read.table("2_loose_kapa_list.txt", header = F, sep = "|")
#PAV_genes <- read.table(paste0({Abb_slurm},".GenePAV.txt"), header = TRUE)
core_genes <- read.table("popstat.ZT.4939.kappa3.05.CoreGenes.txt", header = F, sep = "|")

# remove genes in PAV
run_loose_kappa <- run_loose_kappa %>% filter(V1 %in% core_genes$V1)

# get mean and median
kappa_mean=format(round(mean(run_loose_kappa$V3), 2), nsmall = 2)
kappa_median=format(round(median(run_loose_kappa$V3), 2), nsmall = 2)

# number of genes
number_of_genes=(length(run_loose_kappa$V3))

# plot
#ggplot(run_loose_kappa_ZT, aes(x=V3)) + geom_histogram(bins=5000)

# with density line
ggplot(run_loose_kappa) + 
  aes(x=V3) +
  geom_histogram(aes(y = ..density..), alpha = 0.7, bins = 40) +
  geom_vline(xintercept=mean(run_loose_kappa$V3), color = "red") +
  geom_vline(xintercept=median(run_loose_kappa$V3), color = "blue") +
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7), limits=c(0,0.7)) +
  scale_x_continuous("kappa", breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10)) +
  geom_line(aes(y = ..density..), stat = 'density', color = "orange", size = 1) +
  ggtitle(paste(Num_chr,"isolates,",number_of_genes,"genes,","median kappa:",kappa_median,"(blue line)")) +
  theme_classic()

ggsave(file = paste0({Abb_slurm},".",{Num_chr},".",{number_of_genes},"genes.kappa",{kappa_median},".core.jpeg"), width = 9, height = 5, dpi = 300, units = "in", device='jpeg')

write.table(previous_dataset, paste0("popstat.",{Abb_slurm},".",{number_of_genes},".kappa",{kappa_median},".CoreGenes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "|")

library(vegan)
library("dunn.test")
library(ape)
library(cowplot)
library(ggplot2)
library(dplyr)
library(nlme)

# load sputum microbiome and sample metadata
bac_3 <- read.csv("../shell_result/kraken_clean/bac_name.txt",
                               row.names = 1,sep = "\t")
metadata <-  read.csv("../metadata/metadata.csv")

# adjust index match
metadata <- metadata[match(colnames(bac_bck),metadata$sample),]

# shannon index
alpha = "shannon" 
shannon_bac <- diversity(bac3,index=alpha)
# add to metadata 
shan_bac <- data.frame(shannon_bac)
shan_bac <- cbind(shan_bac,metadata)


# -----------------------alpha diversity and plot-------------------------------
# health vs cap 
statis <- wilcox.test(shan_bac1$shannon_bac ~ shan_bac1$group) 
pvalue <- format(statis$p.value, scientific = TRUE, digits = 3)

# plot alpha diversity
e <- ggplot(shan_bac1, aes(x = group, y = shannon_bac,fill=group)) + 
  geom_boxplot(width=0.3) +
  ylim(c(-0.1, 6.5))+
  scale_fill_manual(values = c("#e77d72","#6eab4d"))+
  theme_bw() +theme(panel.grid=element_blank()) +
  ylab("Shannon index") + xlab("")+
  geom_text(x=1,y=6,label=paste0("p=",pvalue))+
  theme(axis.text.x = element_text(size=10), axis.text.y =
          element_text(size=10))
e

# ------------------beta diversity and plot-------------------------------------
bac_bray_dis <- as.matrix(vegdist(bac3, method = "bray"))

#vectors: The principal coordinates with positive eigenvalues.
PCOA <- pcoa(bac_bray_dis)$vectors
var_exp <- pcoa(bac_bray_dis)$values
barplot(var_exp$Eigenvalues) 
beta_dist = as.dist(bac_bray_dis) 


# Run PERMANOVA
ad = adonis(beta_dist ~ metadata$group, permutations=999) 
p_val <- ad$aov.tab[1,6] 
r_sq <- ad$aov.tab[1,5] 


# Run stats on the coordinates
PCOA <- PCOA[metadata$sample,]

# Plot beta diversity PCoA
for(i in 1:ncol(PCOA)){colnames(PCOA)[i] <- paste("PC",i, sep = "")}

PCOA <- cbind(PCOA, rownames(PCOA)) 
colnames(PCOA)[ncol(PCOA)] <- "SampleID"
m <- cbind(metadata, rownames(PCOA))

m <- data.frame(lapply(m, as.character), stringsAsFactors=FALSE) 

colnames(m)[ncol(m)] <- "SampleID"
PCOA <- merge(PCOA, m, by = "SampleID")

PCOA$PC1 <- as.numeric(as.character(PCOA$PC1))
PCOA$PC2 <- as.numeric(as.character(PCOA$PC2))

body_cols = c("#e77d72","#413494","#6eab4d")

pcoa_p <- ggplot(PCOA) +
  geom_point(size = 2, alpha=0.65, aes(x = PC1, y = PC2,color = group))+
  scale_color_manual(values = body_cols) + 
  stat_ellipse(aes(x = PC1, y = PC2, color = group), level = 0.95, 
               alpha = 0.5, show.legend = T) +
  annotate("text", x = -0.45, y = -0.2, label= paste("P = ", p_val), size=2) + 
  annotate("text", x = -0.45, y = -0.25, 
           label= paste("R2 = ", round(r_sq,digits=3)), size=2) + 
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(color=NA), 
        axis.text.y = element_text(color=NA))+
  theme_bw()+theme(panel.grid = element_blank(),legend.position = 'none')


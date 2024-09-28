rm(list=ls())


library(readxl)
library(dplyr)
library(lme4)
library(vegan)
library(ggplot2)
library(gg.gap)
library(patchwork)
library(ape)

# --------------------肠道菌群过滤----------------------------------------------
gut_file <- "gut_species.csv"
gut <- read.csv(gut_file,sep = "\t",row.names = 1) 


match_sample <- intersect(metadata$sample,colnames(gut))

# 提取metadata 和测序数据
gut <- gut[, match_sample]
metadata <- metadata[match(match_sample, metadata$sample),]

# 提取细菌的丰度数据
gut1 <- gut[grepl("k__Bacteria", rownames(gut)),]

# 去除在所有样本中表达量为0的物种
gut1  <- gut1[rowSums(gut1)>0,] 

#只保留在20%人中比例高于0.01%丰度的物种
gut_f <- gut_sp[rowSums(gut_sp>0.01)>0.2*dim(gut_sp)[2],]


# 1. PCoA 图--------------------------------------------------------------------
bac3 <- t(gut_f)
# 计算bray距离，列是物种名字，行是样本名字
bac_bray_dis <- as.matrix(vegdist(bac3, method = "bray"))
PCOA <- pcoa(bac_bray_dis)$vectors
var_exp <- pcoa(bac_bray_dis)$values
# 查看Eigenvalues的分布
barplot(var_exp$Eigenvalues) 

#将 beta_table 转换为距离矩阵。as.dist 函数将相异性矩阵
# 转换为 "dist" 类型的对象，以便进行距离分析。
beta_dist = as.dist(bac_bray_dis) 

# Run PERMANOVA
ad = adonis(beta_dist ~ metadata$group, permutations=999) 
p_val <- ad$aov.tab[1,6] 
r_sq <- ad$aov.tab[1,5] 
p_val
r_sq


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


body_cols=c("#f0999f","#46bbc0")

body_PCOA <- ggplot(PCOA) +
  geom_point(size = 1, alpha=0.85, aes(x = PC1, y = PC2,color = group))+
  scale_color_manual(values = body_cols) + 
  stat_ellipse(aes(x = PC1, y = PC2, color = group), level = 0.95, 
               alpha = 0.5, show.legend = T) +
  annotate("text", x = -0.45, y = -0.2, label= paste("P = ", p_val), size=2) + 
  annotate("text", x = -0.45, y = -0.25, 
           label= paste("R2 = ", round(r_sq,digits=3)), size=2) + 
  labs(x = paste("PC1 (", 
                        round(var_exp$Relative_eig[1], digits=3)*100, "%)", 
                        sep = ""),
      y= paste("PC2 (", 
                        round(var_exp$Relative_eig[2], digits=3)*100, "%)", 
                        sep = ""))+
  theme(axis.text.x = element_text(color=NA), 
        axis.text.y = element_text(color=NA))+
  theme_bw()+theme(panel.grid = element_blank())


body_PCOA

#2. 绘制肠道差异菌群的散点图----------------------------------------------------
# 提取排名前30的物种名字
top30 <- sigdiff_tax_glm[abs(sigdiff_tax_glm$logfc)>2,"taxon"]
gut_f30 <- gut_f[top30$taxon,]

# 按照中位数从大到小排列
gut_f30 <- gut_f30[order(apply(gut_f30,1,median),decreasing = T),]

dat1 <- t(gut_f30)
dat1 <- cbind(metadata[,c("sample", "group")], dat1)
data2<-melt(dat1,id.vals=c("sample", "group")) 
data2$variable=factor(data2$variable,
               levels= rownames(gut_f30))

Fig1b <- ggplot(data2,aes(x=variable,y=value),colour=factor(group))+
  geom_boxplot(aes(fill=group),outlier.colour=NULL,outlier.shape=21,outlier.size=1.5)+
  scale_fill_manual(values=c("#f0999f","#46bbc0"))+
  xlab("") + ylab("Normalized relative abundance") +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   #  panel.border = element_blank(),
                   panel.background = element_blank(),
                   axis.text.x=element_text(angle=45, hjust=1, vjust=1, 
                                            color="black",size=10),
                   axis.text.y=element_text(color="black",size=10),
                   axis.title.x=element_text(),
                   axis.title.y=element_text())

Fig1b


#截断坐标轴
p2 =gg.gap(plot = Fig1b,
           segments = c(40, 65),
           tick_width = 5,
           rel_heights = c(0.8, 0, 0.1),
           ylim = c(0, 100)
)

p2

# 绘制多组学九象限图-----------------------------------------------------------
# 导入呼吸道微生物组差异物种数据
resp_difftax <- read.csv("taxa_different_analysis.csv")
# 获得共有物种
bitax <- intersect(allgut$taxon, resp_difftax$taxon)

# 提取共有物种的logfc 和fdr，然后计算是否显著
gut_bitax <- glm_res %>% filter(taxon %in% bitax) 
# 呼吸道微生物组
res_bitax <- resp_difftax %>% filter(taxon %in% bitax) 
#合并
all_bitax <- merge(gut_bitax,res_bitax,by="taxon") 


group_all <- ifelse(all_bitax$group_gut=="sigdiff" & all_bitax$group_res == "sigdiff","both",
                    ifelse(all_bitax$group_gut=="sigdiff", "gut_only",
                           ifelse(all_bitax$group_res=="sigdiff","resp_only","none")))

all_bitax$all_group <- group_all


# -------------------------------开始绘制九象限图-------------------------------
draw_axis_line <- function(length_x, length_y, 
                           tick_step = NULL, lab_step = NULL){
  axis_x_begin <- -1*length_x
  axis_x_end <- length_x
  
  axis_y_begin  <- -1*length_y
  axis_y_end    <- length_y
  
  if (missing(tick_step))
    tick_step <- 1
  
  if (is.null(lab_step))
    lab_step <- 2
  
  # axis ticks data
  tick_x_frame <- data.frame(ticks = seq(axis_x_begin, axis_x_end, 
                                         by = tick_step))
  
  tick_y_frame <-  data.frame(ticks = seq(axis_y_begin, axis_y_end, 
                                          by = tick_step))
  
  # axis labels data
  lab_x_frame <- subset(data.frame(lab = seq(axis_x_begin, axis_x_end, 
                                             by = lab_step), zero = 0), 
                        lab != 0)
  
  lab_y_frame <- subset(data.frame(lab = seq(axis_y_begin, axis_y_end,
                                             by = lab_step),zero = 0), 
                        lab != 0)
  
  # set tick length
  tick_x_length = 0.05
  tick_y_length = 0.05
  
  # set zero point
  
  data <- data.frame(x = 0, y = 0)
  p <- ggplot(data = data) +
    
    # draw axis line
    geom_segment(y = 0, yend = 0, 
                 x = axis_x_begin, 
                 xend = axis_x_end,
                 linewidth = 0.5) + 
    geom_segment(x = 0, xend = 0, 
                 y = axis_y_begin, 
                 yend = axis_y_end,
                 linewidth = 0.5) +
    # x ticks
    geom_segment(data = tick_x_frame, 
                 aes(x = ticks, xend = ticks, 
                     y = 0, yend = 0 - tick_x_length)) +
    # y ticks
    geom_segment(data = tick_y_frame, 
                 aes(x = 0, xend = 0 - tick_y_length, 
                     y = ticks, yend = ticks)) + 
    
    # labels
    geom_text(data=lab_x_frame, aes(x=lab, y=zero, label=lab), vjust = 1.5) +
    geom_text(data=lab_y_frame, aes(x=zero, y=lab, label=lab), hjust= 1.5) +
    theme_minimal()+
    theme(panel.grid = element_blank(),axis.text = element_blank())
  return(p)
}

p <- draw_axis_line(5,7)
p1 <- p + geom_point(data=all_bitax, aes(gut_logfc, logfc, color = all_group))+
  geom_point(data = all_bitax[all_bitax$all_group == "both",], 
             aes(gut_logfc, logfc),
             color="black",shape=21,size=2)+
  scale_color_manual(values = c("both" = "#dd8653", 
                                "gut_only" = "#59a5d7", 
                                "resp_only" = "#aa65a4", 
                                "#878787"),
                     breaks = c("both","gut_only","resp_only"))+
  geom_text(data = all_bitax[all_bitax$all_group == "both",],
            aes(gut_logfc, logfc,label=taxon,color="#dd8653"),size= 3,vjust=-1,hjust= 0.8)+
  xlab("gut microbiome")+
  ylab("respiratory microbiome")+
  theme(legend.position = c(0.6,0.2))+
  guides(color = guide_legend(title = "", ncol = 1, byrow = TRUE))

p1









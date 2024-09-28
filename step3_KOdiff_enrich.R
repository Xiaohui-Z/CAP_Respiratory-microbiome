rm(list = ls())
options(stringsAsFactors = F)

# 加载包
library(edgeR)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(readxl)
library(statmod)
library(limma)
library(ggrepel)


#-------------------------------------load data---------------------------------

lname <- load(file = "KO_diff.Rdata")

data_cap <- metadata[metadata$group %in% c("health","cap" ),] 
  
data_cap <- data_cap %>% arrange(group,subject,day_group)
data_cap1 <- data_cap[!duplicated(data_cap$subject),] 
table(data_cap1$group) 
huko_cap <- hu_ko[,match(data_cap1$sample, colnames(hu_ko))] 
#过滤低水平ko，去除在10%以下样本中表达的ko
huko_cap <- huko_cap[rowSums(huko_cap>0) > 0.1* ncol(huko_cap),] 
  
# 调整到相同的样本顺序
huko_cap <- huko_cap[,data_cap1$sample]

# limma
group_list <- factor(data_cap1$group, levels = c("cap","health"))
head(group_list)

design <- model.matrix(~0+factor(group_list))
rownames(design) <- colnames(huko_cap)
colnames(design) <- levels(factor(group_list))
head(design)

# 设置对比矩阵
contrast.matrix <- makeContrasts(health ~ cap, levels=colnames(design)) 
print(contrast.matrix)

fit <- lmFit(huko_cap, design)
# 应用对比矩阵
fit2 <- contrasts.fit(fit, contrast.matrix)
# 应用贝叶斯平滑
fit2 <- eBayes(fit2,trend=TRUE)
DEG_lim <- topTable(fit2, adjust="BH", number=Inf)

# 筛选上下调，设定阈值
fc_cutoff <- 2 
pvalue <- 0.1

DEG_lim$regulated <- "normal"

loc_up <- intersect(which(DEG_lim$logFC>log2(fc_cutoff)),
                    which(DEG_lim$FDR<pvalue))
loc_down <- intersect(which(DEG_lim$logFC < (-log2(fc_cutoff))),
                      which(DEG_lim$FDR<pvalue))

DEG_lim$regulated[loc_up] <- "up"
DEG_lim$regulated[loc_down] <- "down"

#-----------------------------绘制火山图----------------------------------------
#提取上下调前10的基因
up10sig <- filter(DEG_wc,regulated=="up") %>%
  distinct(KO,.keep_all = T) %>% top_n(10,abs(logFC))

down10sig <- filter(DEG_wc,regulated=="down") %>%
  distinct(KO,.keep_all = T) %>% top_n(10,abs(logFC))

top10sig <- rbind(up10sig,down10sig)
top10sig$label <- top10sig$gene

# 将top10的label加到总数据中
data <- merge(DEG_wc,top10sig[,c("KO","label"),],
              by = "KO",all.x = TRUE)


#开始画图
valco <- ggplot(data,aes(x=logFC, y = -log10(FDR)))+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "#999999")+
  geom_hline(yintercept = -log10(0.1),linetype = "dashed", color="#999999")+
  geom_point(aes(color=-log10(FDR),size = -log10(FDR)))+
  scale_color_gradientn(values=seq(0,1,0.2),
                        colors = c("#314194","#6cc3de",
                                   "#fae55d","#cd336b","#a33430"))+
  scale_size_continuous(range = c(0.8,1.5))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.05,0.9),
        legend.justification = c(0,1)
  ) +
  guides(col = guide_colorbar(title="-Log10_FDR"),
         size = "none"
  )+
  geom_text_repel(aes(label = label), size = 4) +
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value")

#-----------------------------------基因富集分析--------------------------------   
library(MicrobiomeProfiler)
uplist <- DEG_lim[DEG_lim$regulated=="up","KO"]
downlist <- DEG_lim[DEG_lim$regulated=="down","KO"]

# 开始富集分析
kk.up<- enrichKO(
  uplist,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  #universe = DEG_wc$KO,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2
)
#绘制低表达基因
kk.down <- enrichKO(
  downlist,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  #universe = DEG_wc$KO,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2
)


#整合成图形
kegg_down_dt <- as.data.frame(kk.down)
kegg_up_dt <- as.data.frame(kk.up)

#计算富集的倍数，包装成函数
enrich_fold <- function(dat){
  GeneRatio <- as.numeric(lapply(strsplit(dat$GeneRatio,split="/"),function(x) 
    as.numeric(x[1])/as.numeric(x[2])))
  BgRatio <- as.numeric(lapply(strsplit(dat$BgRatio,split="/"),function(x) 
    as.numeric(x[1])/as.numeric(x[2])  ))
  enrich_factor <- GeneRatio/BgRatio
  return(enrich_factor)
}

kegg_down_dt$enrich_factor <- enrich_fold(kegg_down_dt)
kegg_up_dt$enrich_factor <- enrich_fold(kegg_up_dt)


down_kegg <- kegg_down_dt[(kegg_down_dt$p.adjust<0.1),] 
down_kegg$group=-1
up_kegg <- kegg_up_dt[kegg_up_dt$p.adjust<0.1,]
up_kegg$group=1

kegg_enrich_all <- rbind(up_kegg,down_kegg)

# 导入绘图工具
source('New.functions.R')
kegg_figure=go.kegg_plot(up.data=up_kegg,down.data=down_kegg)
print(kegg_figure)

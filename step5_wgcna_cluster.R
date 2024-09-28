rm(list = ls())

# Load the WGCNA package
library(WGCNA)
library(tidyverse)
library(cowplot)
library(psych)
library(igraph)


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "wgcna_RData.RData");
#The variable lnames contains the names of loaded variables.
lnames

cluster_kraken <- read.csv("respmicro_host_association.csv",
                           row.names = 1)
library(clusterProfiler)
library(org.Hs.eg.db)



tax_cluster5 <- rownames(cluster_kraken[cluster_kraken$cluster == "5",])

#提取salmon的ens id
target_color <- "salmon"
color_id <- colnames(datExpr)[moduleColors==target_color]
allid <- colnames(datExpr)

length(color_id);length(allid)

#### 第一步，从org.Hs.eg.db提取ENSG的ID 和GI号对应关系
keytypes(org.Hs.eg.db)

# bitr in clusterProfiler 
allID <- bitr(allid, fromType = "ENSEMBL", 
              toType = c( "ENTREZID" ), 
              OrgDb = org.Hs.eg.db ) 

degID <- bitr(color_id, fromType = "ENSEMBL", 
              toType = c( "ENTREZID" ), 
              OrgDb = org.Hs.eg.db ) 

enrich <- enrichKEGG(gene = degID[,2],
                     organism='hsa',
                     universe=allID[,2],
                     pvalueCutoff=1,
                     qvalueCutoff=1)

# 计算富集因子
GeneRatio <- as.numeric(lapply(strsplit(enrich$GeneRatio,split="/"),function(x)
  as.numeric(x[1])/as.numeric(x[2])))
head(GeneRatio)

BgRatio <- as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x)
  as.numeric(x[1])/as.numeric(x[2])  ))
head(BgRatio)

enrich_factor <- GeneRatio/BgRatio

out <- data.frame(enrich$ID,
                  enrich$Description,
                  enrich$GeneRatio,
                  enrich$BgRatio,
                  round(enrich_factor,2),
                  enrich$pvalue,
                  enrich$qvalue,
                  enrich$geneID)

colnames(out) <- c("ID","Description","GeneRatio","BgRatio","enrich_factor",
                   "pvalue","qvalue","geneID")

out_sig0.05 <- out[out$qvalue<0.05,]

# barplot
bar <- barplot(enrich,showCategory=10,title="KEGG Pathway",
               colorBy="pvalue")
bar

mutate(enrich, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore",colorBy = "qscore")



# go cc bp bar
DEG <- geneInfo0[geneInfo0$moduleColor==target_color,"geneSymbol"]

ego_CC <- enrichGO(gene=DEG, OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', 
                   ont="CC", pvalueCutoff= 0.05,qvalueCutoff= 1)
ego_MF <- enrichGO(gene=DEG, OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', 
                   ont="MF", pvalueCutoff= 0.05,qvalueCutoff= 1)
ego_BP <- enrichGO(gene=DEG, OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', 
                   ont="BP", pvalueCutoff= 1,qvalueCutoff= 1)


if(F){
  p_BP <- barplot(ego_BP,showCategory = 10,color = "pvalue") + ggtitle("Biological process")+
    scale_fill_gradient(low = "#d31d2c", high = "#626cad",limits = c(0, 0.01))
  p_CC <- barplot(ego_CC,showCategory = 10,color = "pvalue") + ggtitle("Cellular component")+
    scale_fill_gradient(low = "#d31d2c", high = "#626cad",limits = c(0, 0.01))
  p_MF <- barplot(ego_MF,showCategory = 10,color = "pvalue") + ggtitle("Molecular function")+
    scale_fill_gradient(low = "#d31d2c", high = "#626cad",limits = c(0, 0.01))
}

# 使用ggplot2 绘图，根据count从大到小排列
ego_BP_df <- as.data.frame(ego_BP@result) %>%
  arrange(desc(Count))

#绘图
p_BP <- ggplot(ego_BP_df[1:10, ], aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_bar(stat = "identity",width = 0.8) +
  coord_flip() +
  ggtitle("Biological process") +
  scale_fill_gradientn(colors = c("#801f51", "#ffa500", "#4299c8"),
                       values = scales::rescale(c(0, 0.025, 0.05)),
                       limits = c(0, 0.05)) +
  theme_minimal()+
  xlab("")+
  theme(axis.text = element_text(size=8))
p_BP



# 2. cellular component
ego_CC_df <- as.data.frame(ego_CC@result) %>%
  arrange(desc(Count))
# Plot with ggplot2
p_CC <- ggplot(ego_CC_df[1:10, ], aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_bar(stat = "identity",width = 0.8) +
  coord_flip() +
  ggtitle("Biological process") +
  # scale_fill_gradient(low = "#d31d2c", high = "#626cad", limits = c(0, 0.05)) +
  scale_fill_gradientn(colors = c("#801f51", "#ffa500", "#4299c8"),
                       values = scales::rescale(c(0, 0.025, 0.05)),
                       limits = c(0, 0.05)) +
  theme_minimal()+
  xlab("")+
  theme(axis.text = element_text(size=8),legend.position = "none")
p_CC


# 3. molecular function
ego_MF_df <- as.data.frame(ego_MF@result) %>%
  arrange(desc(Count))
# Plot with ggplot2
p_MF <- ggplot(ego_MF_df[1:10, ], aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
  geom_bar(stat = "identity",width = 0.8) +
  coord_flip() +
  ggtitle("Biological process") +
  # scale_fill_gradient(low = "#d31d2c", high = "#626cad", limits = c(0, 0.05)) +
  scale_fill_gradientn(colors = c("#801f51", "#ffa500", "#4299c8"),
                       values = scales::rescale(c(0, 0.025, 0.05)),
                       limits = c(0, 0.05)) +
  theme_minimal()+
  xlab("")+
  theme(axis.text = element_text(size=8),legend.position = "none")
p_MF

# 合并图，后两个图例可以不要
plotc <- p_BP+p_CC+p_MF
plotc

#===============================================================================
# 绘制基因物种互作关系网络，使用cytoscape                                      /
#===============================================================================
# 导入基因注释文件
id2sym <- paste0("all_id2name.txt")

annot = read.csv(file = id2sym,sep="\t",
                 col.names = c("ensemb","protein","symbol"));
head(annot)


target_color <- "salmon" 
target_cluster <- 5

color_id <- colnames(datExpr)[moduleColors==target_color]

# 提取表达量信息
expres_color <- datExpr[,color_id]

color_gene = annot[match(colnames(expres_color), annot$ensemb),"symbol"]
sum(is.na(color_gene)) # 

# 修改列名ens id 为基因名字
colnames(expres_color) <- color_gene

abun_cluster5 <- cluster_kraken %>% filter(cluster == target_cluster) 
abun_cluster5 <- t(abun_cluster5[,1:64])

sum(rownames(abun_cluster5) != rownames(expres_color)) 
cluster_tax_gene <- cbind(abun_cluster5,expres_color) 


# ----------------------------开始进行网络图绘制-------------------------------
# 转换otu-sample矩阵，行为sample，列为物种名字和基因
otu <- cluster_tax_gene

otu <- scale(otu) # 默认按照列
corr_otu <- corr.test(otu,method = "spearman", use = "pairwise")
corr_otu.r = corr_otu$r # 取相关性矩阵R值
corr_otu.p = corr_otu$p # 取相关性矩阵p值

# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
corr_otu.r[corr_otu.p>0.05|abs(corr_otu.r)<0.3] = 0 
s <- corr_otu.r
m <- graph.adjacency(s, weighted = TRUE, mode = "undirected", diag = FALSE)

#去除自相关和孤立的点
m <- simplify(m)
m <- delete.vertices(m, names(degree(m)[degree(m) == 0]))

E(m)$correlation <- E(m)$weight
E(m)$weight <- abs(E(m)$weight)

names(V(m)) 
V(m)$item <- c(rep("microbe",ncol(abun_cluster5)),rep("gene",ncol(expres_color)))


adj <- as.matrix(get.adjacency(m, attr = 'correlation')) 

#保存边的输出
edge <- data.frame(as_edgelist(m, names = T))
df <- as.data.frame(E(m)$correlation)
df[df>0] <- 1
df[df<0] <- -1
colnames(df) <- 'cor'
edge <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(m)$weight,
  correlation = E(m)$correlation,
  cor = df)

# 输出节点的表格
node <- data.frame(
  id = names(V(m))
)

# 增加node表达量和类别信息
node$type = "gene"
node[grepl("s__",node$id),"type"] = "microbe"

# 保存为其他格式,gephi和cytoscape
write.graph(m,
            paste0(folder,dataset,"_network.graphml"),
            format = "graphml")
# gml格式，用于cytoscape
write.graph(m, 
            paste0(folder,dataset,"_network.gml"),
            format = 'gml')

#vasulized using cytoscape 




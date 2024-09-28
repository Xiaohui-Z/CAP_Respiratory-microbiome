
# 用于查看不同物种的丰度比

# 加载包
# library(edgeR)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(ggsci)
library(paletteer)
library(patchwork)



rm(list = ls())
options(stringsAsFactors = F)

# 导入数据
load("KO_diff/kodeg_pathway.Rdata")

#根据day_group分组来绘图
metadata <- as.data.frame(metadata)
col_sort <- c(metadata[grepl("T1",metadata$day_group),"sample"],
              metadata[grepl("T2",metadata$day_group),"sample"],
              metadata[grepl("T3",metadata$day_group),"sample"],
              metadata[metadata$day_group=="health","sample"]
)

metadata_sort <- metadata[match(col_sort,metadata$sample),]
# 初始化一个列表来存储每个图表
plots_list <- list()

# 先取出一个画图
for(path_id in sig_path10$ID){
  #path_id <-sig_path10$ID[3]
  path_name <- sig_path10[sig_path10$ID==path_id,"Description"]
  path_cpm <- map_cpm113gen %>% filter(map==path_id) %>% select(-map)
  
  # 去除不含这个pathway的物种
  col_0 <- colnames(path_cpm[,-1])[colSums(path_cpm[,-1]) ==0]
  path_cpm <- path_cpm %>% select(-all_of(col_0))
  
  # 转换数据格式以适应 ggplot
  df_long <- pivot_longer(path_cpm, cols = -genus, 
                          names_to = "sample", values_to = "abundance")
  df_long$sample <- factor(df_long$sample,levels = col_sort)
  
  # 计算每个 genus 的总丰度
  total_abundance <- df_long %>%
    group_by(genus) %>%
    summarize(total_abundance = sum(abundance, na.rm = TRUE))
  
  # 找出总丰度最高的前10个 genus
  top_genuses <- total_abundance %>%
    top_n(10, total_abundance) %>%
    pull(genus)
  
  # 标记前10个 genus 和 others
  df_long$genus_group <- ifelse(df_long$genus %in% top_genuses, df_long$genus, "Others")

  
  # 绘制堆叠柱状图
  
  # 使用igv调色板
  used_color <- "igv"
  stackbar <- FALSE
  
  d_palettes <- palettes_d_names #获得所有调色板
  fills<- as.character(paletteer_d("ggsci::default_igv", n=20)) #随机选一个查看
    

  p_stack <- ggplot(df_long, aes(x = sample, y = abundance, fill = genus_group)) +
      geom_bar(position="fill",stat = "identity") +
      scale_y_continuous(labels = function(x) sprintf("%.0f", x * 100)) + 
      scale_fill_igv (name = "Genus", labels = c(top_genuses, "Others"), 
                      breaks = c(top_genuses, "Others")) +
      theme_minimal() +
      labs(x = "", y = "Relative abundance(%)", fill = "Genus") + 
      ggtitle(path_name)+
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 4),
            title = element_text(size=8),
            axis.ticks.x = element_blank())  
  }
  
  p_stack

  # 缩小图例
  p_stack <- p_stack + theme(legend.position = "right",
                             legend.title = element_blank(),  # 缩小图例标题字体大小
                             legend.text = element_text(size = 6),  # 缩小图例文本字体大小
                             legend.key.size = unit(0.4, "cm"),  # 缩小图例键的大小
                             legend.spacing = unit(0.1, "cm"),  # 调整图例内部组件的间距
                             legend.margin = margin(3, 3, 3, 3, "pt"))  # 缩小图例边距
  
  plots_list[[path_id]] <- p_stack
  
  ggsave(paste0("KO_diff/single_pathway_stackplot/pathway_plot_",path_id,"_",path_name,".pdf"),
         p_stack,
         height = 3,width = 6)

# 使用 patchwork 将所有图表合并为一个图
combined_plot <- wrap_plots(plots_list, ncol = 1)    
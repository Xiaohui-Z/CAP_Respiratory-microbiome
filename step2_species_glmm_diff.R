# 加载r包
library(lme4)
library(dplyr)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(paletteer)
library(tidyr)
library(tidyverse)

# load microbiome data and metadata--------------------
lnames <- load(file = "taxa_different_data.Rdata")
lnames

dataset_list <- c("kraken_allf_percent")

glm_res_merge <- NULL
for(tax_abun in dataset_list){
  
  abun_file <- eval(parse(text = tax_abun))
  abun_cap <- abun_file %>% select(all_of(data_cap$sample))
  row.names(abun_cap) <- abun_file$Species
  
  # match index of metadata and abundance data
  sum(data_cap$sample==colnames(abun_cap)) 
  
  glm_res <- NULL
  for(i in row.names(abun_cap)){

    data_m = cbind(data_cap,scale(t(abun_cap[i,])))
    data_m$age <- scale(data_m$age)[,1]
    
    # calculate fold chage
    abun_tax <- t(abun_file %>% filter(Species==i) %>% select(-c(1:7)))
    colnames(abun_tax) <- "tax"
    abun_tax <- cbind(metadata[,c("sample","group")],abun_tax)
    
    result <- abun_tax %>% group_by(group) %>%
      summarize(mean_abun = mean(tax))
    
    # formula
    fml <- paste0("group ~ age + antibiotic",i,"+ (1|subject)") 
    mod <- glmer(as.formula(fml), family = binomial, 
                       control = glmerControl(optimizer = "bobyqa",
                                              optCtrl = list(maxfun= 10000)), 
                       data = data_m)
    p <- summary(mod)$coefficients[i,4]

    res <- c("taxon" = i,
             "pvalue" = summary(mod)$coefficients[i,4],
             "fc" = result$mean_abun[1]/result$mean_abun[3])
    glm_res <- bind_rows(glm_res, res)
  }
  
  glm_res$fdr <- p.adjust(glm_res$pvalue, method = "fdr")
}

# filter out significatn difference taxon
diff_tax <- glm_res_merge %>% filter(fdr <0.05 & abun_cap >0.5)

# ----------------------------dfferent taxa plot--------------------------------
# heatmap
abundance <- kraken %>% filter(Species %in% tax) 
abundance <- abundance[,order(metadata$group)] #按照分组顺序画
  
row_color_list <- c("g__Actinomyces"= "#8DD3C7",
                        "g__Bergeyella"=  "#C2E7BD",
                        "g__Campylobacter"="#F7FCB4",
                        "g__Capnocytophaga"= "#E4E3C2",
                        "g__Eikenella" ="#FCCDE5" ,
                        "g__Gemella"="#D9D4AB" ,
                        "g__Granulicatella"="#D2A6B7",
                        "g__Lancefieldella"= "#EE8B86",
                        "g__Lautropia"="#DA8D8B" ,
                        "g__Leptotrichia" ="#A0A3B9" ,
                        "g__Neisseria"="#D3B387",
                        "g__Ottowia"="#99B1BC" ,
                        "g__Prevotella"="#C6C3D4"  ,
                        "g__Rothia" = "#B7DC71",
                        "g__Streptococcus"= "#D0CD66",
                        "g__Veillonella"= "#F3B962")


color_band <- colorRampPalette(c("#62bcc2", "white", "#e99c9d"))(50)

  
# annotation columan
annol_col <- data.frame(group=metadata_t$group)
rownames(annol_col) <- metadata_t$sample
annol_col$group <- factor(annol_col$group)
  
# annotation row
annol_row <- data.frame(genus=diff_info$Genus)
rownames(annol_row) <- diff_info$taxon
annol_row$genus <- factor(annol_row$genus)
  
# 形成行注释颜色注释颜色对应表
d_palete <- palettes_d_names
color_row <- as.character(paletteer_d("ggsci::default_igv",n=16))

  
p1 <- pheatmap(abundance,
               scale = "row",
               show_colnames = T,show_rownames = T,
               cluster_cols = F,
               fontsize_col = 6,
               clustering_method = "ward.D", 
               labels_row = tax,
               annotation_col = annol_col,
               annotation_row = annol_row,
               annotation_colors = list(group = c("cap" = "#f0978c",
                                                  "followup" = "#8cb7f9",
                                                  "health" = "#61d26b"),
                                        genus = row_color_list),
               color = color_band,
               angle_col = 45,
               main = paste0(abun_name,"'s diff taxon"))
p1

# plot fold change hist plot
tax_order <- tax[p1$tree_row$order]
p2_data <- diff_info[match(tax_order,diff_info$Species),]
p2_data$taxon <- factor(p2_data$taxon,levels = rev(p2_data$taxon))
  
p2 <- ggplot(p2_data)+geom_col(aes(x=taxon,y=logfc,fill=fdr))+ coord_flip() +
  scale_fill_continuous(low = "#babbfd",high="#f9c7c5",
                        breaks = seq(0, 1, length.out = 51))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p2
  
#P lot a stacked bar chart for multiple time points 

data_long1 <- as.data.frame(t(abundance)) %>% 
  mutate(sampleid = row.names(.)) 
data_long1$timep <- metadata_t[match(row.names(data_long1),metadata_t$sample),"day_group"]
  
  
data_long2 <- data_long1 %>% pivot_longer(cols = -c(sampleid,timep),
                                            names_to = "species",values_to = "value") %>%
    group_by(sampleid) %>%
    mutate(sum_tax = sum(value)) %>%
    arrange(desc(sum_tax)) %>% 
    filter(sum_tax >20) %>% # 
    arrange(timep,desc(sum_tax))  
  
data_long2$timep <- factor(data_long2$timep,levels = unique(data_long2$timep))
data_long2$species <- factor(data_long2$species,levels = unique(data_long2$species))
data_long2$sampleid <- factor(data_long2$sampleid,levels = unique(data_long2$sampleid))
  
  
fills <- colorRampPalette(brewer.pal(8, "Set3"))(40) 
d_palettes <- palettes_d_names 
fills<- as.character(paletteer_d("ggsci::default_igv", n=40)) 
  
# 绘制堆积图
pt <- ggplot(data_long2)+
  geom_bar(aes(sampleid, value, fill = species),
           stat = "identity", position = "stack")+
  scale_fill_manual(values = rev(fills),
                    name = "",
                    guide = guide_legend(ncol = 2,
                                         keywidth = 0.5,
                                         keyheight = 0.5,
                                         reverse = T))+
  #   xlab("")+
  ylab("realtive abundance")+
  theme_classic()+
  theme(legend.position = c(0.9, 0.9),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank()
  )
pt








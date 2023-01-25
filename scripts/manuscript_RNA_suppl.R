library(ggrepel)
library(ggvenn)
library(ggdendro)
library(grid)
library(gridExtra)
library(cowplot)
library(WriteXLS)

# Load data
#############
load("data/ALK_fusion_data.RData")
source("scripts/functions/helpers.R")

# Volcano plots
#################

res_diff_all <- list("EML4-ALK-V1" = res_diff_expr$RNA$EML_V1_induced_vs_Ctrl_induced,
                     "EML4-ALK-V3" = res_diff_expr$RNA$EML_V3_induced_vs_Ctrl_induced,
                     "KIF5B-ALK" = res_diff_expr$RNA$KIF5B_induced_vs_Ctrl_induced,
                     "TFG-ALK" = res_diff_expr$RNA$TFG_induced_vs_Ctrl_induced,
                     "CUTO8+lorlatinib" = NSCLCL_res$CUTO8.24h_Lorlatinib,
                     "CUTO9+lorlatinib" = NSCLCL_res$CUTO9.24h_Lorlatinib,
                     "CUTO29+lorlatinib" = NSCLCL_res$CUTO29.24h_Lorlatinib,
                     "YU1077+lorlatinib" = NSCLCL_res$YU1077.24h_Lorlatinib,
                     "H2228+lorlatinib" = H2228_res,
                     "H3122+lorlatinib" = H3122_res)

genes_to_label<- c(
  "ALK",
  "EML4",
  "KIF5B",
  "TFG",
  rownames(res_diff_all$`EML4-ALK-V1`)[grep("SERPIN",rownames(res_diff_all$`EML4-ALK-V1`))]
)

p_volcano <- list()
for (cl in names(res_diff_all)) {
  res_proc<- as.data.frame(res_diff_all[[cl]])
  res_proc<- res_proc[!is.na(res_proc$log2FoldChange),]
  p <- plot_volcano(DE_results = res_proc, gene = genes_to_label, labelAll = T, labelCol = "black", useHyperbolicTH = T,logFC_cu = 1.5,curve = 0.1, plotTH=F, p_cu = 0.01, plot_nominal_p = F)
  p <- p +
    # Correcting sample names
    ggtitle(cl) +
    theme(
      plot.title = element_text(hjust = 0.5, size=8), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black", size=0.2),
      axis.ticks = element_line(colour = "black", size = 0.2),
      axis.text = element_text(size=6),  
      axis.title = element_text(size=7)
    ) +
    scale_x_continuous(name = "log2(Fold Change)", limits = c(min(res_proc$log2FoldChange),max(res_proc$log2FoldChange))) +
    scale_y_continuous(name = "-log10(Padj)", limits = c(0,max(-log10(res_proc$padj))))
  p_volcano[[cl]] <- p
}

# SERPINB correlation plot
res_diff_logFC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "HGNC", all.x = TRUE),res_diff_all)
res_diff_logFC <- res_diff_logFC[,grep("log2FoldChange",colnames(res_diff_logFC))]
rownames(res_diff_logFC) <- res_diff_all$`EML4-ALK-V1`$HGNC
colnames(res_diff_logFC) <- names(res_diff_all)

# Extracting SERPINB values
sel_samples <- res_diff_logFC[grep("SERPINB",rownames(res_diff_logFC)),]

# Correlating samples
cor_mat <- cor(sel_samples)

# Running hierarchical clustering on the samples to get clustering
sample_order <- order.dendrogram(as.dendrogram(hclust(d = dist(x = t(sel_samples)))))
cor_mat <- cor_mat[sample_order,sample_order]

# Melting data frame
melted_cormat <- reshape2::melt(cor_mat)
colnames(melted_cormat)[3] <- "Pearson's r"
melted_cormat$Cor <- round(melted_cormat$`Pearson's r`,2)

# Plot
p_cor <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=`Pearson's r`)) + 
  geom_tile()+ 
  theme(
    plot.title = element_text(hjust = 0.5, size=8), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(size=6),
    legend.key.size = unit(3, 'mm'),
    legend.text = element_text(size=5),
    # legend.position = "top"
  ) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red") +
  geom_text(aes(label=Cor),size = 2)


# Alternative to 9 messy volcano plots:
#######################################
library(ggdendro)
res_diff_logFC_padj <- do.call(rbind,res_diff_all)
res_diff_logFC_padj <- res_diff_logFC_padj[res_diff_logFC_padj$HGNC%in%genes_to_label,c("log2FoldChange","padj","HGNC")]
res_diff_logFC_padj <- res_diff_logFC_padj[!is.na(res_diff_logFC_padj$padj),]
res_diff_logFC_padj$sample <- gsub("\\..*","",rownames(res_diff_logFC_padj))
res_diff_logFC_padj <- res_diff_logFC_padj[res_diff_logFC_padj$HGNC!="SERPINA11",] # SERPINA11 and SERPINhas multiple NAs and the rest are 0

res_proc <- res_diff_logFC[rownames(res_diff_logFC) %in% genes_to_label,]
res_proc <- res_proc[!is.na(rowSums(res_proc)),]
df1 <- res_proc
df1$Gene <- rownames(res_proc)
df2 <- reshape2::melt(df1, id = "Gene")

colnames(df2) <- c("Gene","sample","logFC")
dendro1 <- as.dendrogram(hclust(d = dist(x = res_proc)))
dendro2 <- as.dendrogram(hclust(d = dist(x = t(res_proc))))
order1 <- order.dendrogram(dendro1)
order2 <- order.dendrogram(dendro2)

dendro.plot1 <- ggdendrogram(data = dendro1, rotate = TRUE,labels = FALSE,theme_dendro = T)
dendro.plot2 <- ggdendrogram(data = dendro2, rotate = FALSE,labels = FALSE,theme_dendro = T)

# Order rows by clustering from dendro
res_diff_logFC_padj$HGNC <- factor(x = df2$Gene,
                                   levels = df1$Gene[order1],
                                   ordered = TRUE)

p_logFC_heatmap <- ggplot(data = res_diff_logFC_padj, aes(x = factor(sample, levels = colnames(res_proc)[order2]), y = HGNC)) +
  geom_point(aes(fill=log2FoldChange, size=-log10(padj)), colour="black", shape=21, stroke = 0.01)+
  scale_fill_gradient2(low = "blue",mid = "white",high = "red") +
  theme(
    plot.title = element_text(hjust = 0.5, size=8), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(size=6),
    legend.key.size = unit(3, 'mm'),
    legend.text = element_text(size=5),
    legend.key = element_rect(fill = "transparent")
  ) 



# Merge plots
#############

# Volcano plots
p_1<- plot_grid(
  p_volcano$`EML4-ALK-V1`,p_volcano$`KIF5B-ALK`,p_volcano$`TFG-ALK`,
  p_volcano$`CUTO8+lorlatinib`,p_volcano$`CUTO9+lorlatinib`,p_volcano$`CUTO29+lorlatinib`,
  p_volcano$`YU1077+lorlatinib`,p_volcano$`H2228+lorlatinib`,p_volcano$`H3122+lorlatinib`,
  ncol = 3,
  nrow = 3
)

# Correlation plot
p_2<- plot_grid(
  p_cor,NULL,
  ncol = 2,
  nrow = 1,
  rel_widths = c(1.3,1.1)
)

# Put together
p <- plot_grid(
  p_1, p_2,
  nrow = 2,
  ncol = 1,
  rel_heights = c(2,0.5),
  labels = c("a","b"),
  label_size = 7
)

# Save as pdf
ggsave("temp/manuscript_RNA_suppl.pdf",  p, width = 200, height = 270, units = "mm")


# Alternative plot
p2 <- plot_grid(
  p_cor,p_logFC_heatmap,
  nrow = 1,
  rel_widths = c(1.5,1),
  labels = c("a","b"),
  label_size = 7
)
ggsave("temp/manuscript_RNA_suppl_2.pdf",  p2, width = 210, height = 120, units = "mm")


fusions <- c("EML_V1","EML_V3","KIF5B","TFG")
# Getting only DE genes
DE_res <- list()
for (f in fusions) {
  res_proc<- as.data.frame(res_diff_expr$RNA[[paste0(f,"_induced_vs_Ctrl_induced")]])
  res_proc <- res_proc[!is.na(res_proc$log2FoldChange),]
  # temp_DE <- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = res_proc$HGNC, th_logFC = 1.5, th_logP = 2, curve = 0.1)
  # res_proc <- res_proc[res_proc$HGNC %in% temp_DE$DE,]
  DE_res[[paste0(f,"_induced_vs_Ctrl_induced")]] <- res_proc
  print(f)
}

# Supplementary RNA DESeq result excel file
WriteXLS(x = DE_res,
         ExcelFileName = "results/tables/manuscript_RNA_res_suppl.xlsx",
         SheetNames = c("EML4-ALK-V1","EML4-ALK-V3","KIF5B-ALK","TFG-ALK"),
         row.names = TRUE)


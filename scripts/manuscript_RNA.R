library(reshape2)
library(ggrepel)
library(ggvenn)
library(grid)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(gplots)

# Load data
#############
load("data/ALK_fusion_data.RData")
source("scripts/functions/helpers.R")

# Volcano EML_V3
#################

# Load data
res_proc<- as.data.frame(res_diff_expr$RNA$EML_V3_induced_vs_Ctrl_induced)

# Decide which genes to label: ALK, SERPIN, EML4
genes_to_label<- c(
  "ALK",
  "EML4",
  "KIF5B",
  "TFG",
  rownames(res_proc)[grep("SERPINB",rownames(res_proc))]
)

# Plot
p_volc<- plot_volcano(DE_results = res_proc, gene = genes_to_label, labelAll = T, labelCol = "black", useHyperbolicTH = T,logFC_cu = 1.5,curve = 0.1, plotTH=F, p_cu = 0.01, plot_nominal_p = F)
p_volc<- p_volc +
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
  scale_x_continuous(name = "log2(Fold Change)", limits = c(-11,20)) +
  scale_y_continuous(name = "-log10(Padj)", limits = c(0,45))
# Few point wich extremely high logFC not shown!  

# Venn all
###########

# Get DE
genes_DE_ls<- list()
for(f in c("EML_V1", "EML_V3", "KIF5B", "TFG")){
  res_proc<- as.data.frame(res_diff_expr$RNA[[paste0(f,"_induced_vs_Ctrl_induced")]])
  res_proc<- res_proc[!is.na(res_proc$log2FoldChange),]
  genes_DE_ls[[f]]<- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = res_proc$HGNC, th_logFC = 1.5, th_logP = 2, curve = 0.1)
}
genes_DE<- sapply(genes_DE_ls, function(x) x$DE)

# Plot
p_venn<- ggvenn(
  genes_DE, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.2, set_name_size = 7*0.36,
  text_size = 6*0.36,
  show_percentage = F
)

# How many?
sapply(genes_DE_ls, function(x) length(x$DE))
# EML_V1 EML_V3  KIF5B    TFG 
# 1128    887   1707    553 
sapply(genes_DE_ls, function(x) length(x$up))
# EML_V1 EML_V3  KIF5B    TFG 
# 686    597    710    377 
sapply(genes_DE_ls, function(x) length(x$down))
# EML_V1 EML_V3  KIF5B    TFG 
# 442    290    997    176 

# GSEA
######

# Intersect Venn?
genes_is<- sort(intersect(
  intersect(genes_DE$EML_V1, genes_DE$EML_V3),
  intersect(genes_DE$KIF5B, genes_DE$TFG)
))

# GSEA: TNF alpha 1
GSEA<- do_GSEA2(genes_retrieved = genes_is, genes_all = rownames(res_proc), GSEA_db = geneset_ls$Ha_ls, min_genes = 2, isList = T)
GSEA$pathway<- factor(rownames(GSEA),levels=rev(rownames(GSEA)))
GSEA$q<- as.numeric(GSEA$q)
GSEA_sel<- GSEA[GSEA$q<0.01,]

# Barplot
p_GSEA<- ggplot(GSEA_sel, aes(x = pathway, y = -log10(q))) +
  geom_bar(stat="identity", fill="darkblue") +
  coord_flip(expand = F) +
  ylab("-log10(Padj)") +
  xlab("") +
  theme(
    plot.title = element_text(hjust = 0.5, size=8, face = "italic"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
  )

GSEA_sel[,"q",drop=F]
# HALLMARK_TNFA_SIGNALING_VIA_NFKB           4.549819e-14
# HALLMARK_IL2_STAT5_SIGNALING               2.030362e-07
# HALLMARK_KRAS_SIGNALING_UP                 3.237943e-05
# HALLMARK_HYPOXIA                           3.359362e-05
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION 1.605846e-04

# SERPINB4 boxplot
##################

# Dataframe
dataType<- "RNA"
sample_info_tmp<- sample_info[[dataType]]
counts_all<- normalized_counts[[dataType]]
idx_SERPIN<- grep("SERPINB",rownames(counts_all))
for(i in 1:length(idx_SERPIN)){
  idx<- idx_SERPIN[i]
  gene_tmp<- rownames(counts_all)[idx]
  counts_tmp<- counts_all[idx,]
  gene_df_tmp<- data.frame(counts=as.numeric(t(counts_tmp)), gene=gene_tmp, fusion=NA,condition=NA)
  gene_df_tmp$fusion<- sample_info_tmp[names(counts_tmp),"fusion"]
  gene_df_tmp$condition[!sample_info_tmp[names(counts_tmp),"isInduced"]&!sample_info_tmp[names(counts_tmp),"hasLorla"]]<- "Non-induced"
  gene_df_tmp$condition[sample_info_tmp[names(counts_tmp),"isInduced"]&!sample_info_tmp[names(counts_tmp),"hasLorla"]]<- "Induced"
  gene_df_tmp$condition[sample_info_tmp[names(counts_tmp),"isInduced"]&sample_info_tmp[names(counts_tmp),"hasLorla"]]<- "Induced + Lorla"
  gene_df_tmp$condition<- factor(gene_df_tmp$condition, levels=c("Non-induced","Induced","Induced + Lorla"))
  if(i==1) gene_df<- gene_df_tmp
  else gene_df<- rbind(gene_df, gene_df_tmp)
}

# Adding non-induced control for completeness
gene_df$fusion_fullName<- paste0("Ctrl(",gene_df$condition,")")
gene_df$fusion_fullName[gene_df$fusion=="EML_V1"]<- "EML4-ALK-V1"
gene_df$fusion_fullName[gene_df$fusion=="EML_V3"]<- "EML4-ALK-V3"
gene_df$fusion_fullName[gene_df$fusion=="KIF5B"]<- "KIF5B-ALK"
gene_df$fusion_fullName[gene_df$fusion=="TFG"]<- "TFG-ALK"
gene_df$fusion_fullName <- paste0(gene_df$fusion,"_",gene_df$condition)
gene_df <- gene_df[gene_df$gene=="SERPINB4",]

# Separating induced and non-induced
gene_df_non_ind <- gene_df[gene_df$condition=="Non-induced",]
gene_df_ind <- gene_df[gene_df$condition=="Induced",]

# T-test
########

# Induced
t_test_res_ind <- compare_means(counts ~ fusion, data = gene_df_ind, method="t.test",ref.group = "Ctrl")
t_test_res_ind <- as.data.frame(t_test_res_ind)
# .y.    group1 group2        p  p.adj p.format p.signif method
# <chr>  <chr>  <chr>     <dbl>  <dbl> <chr>    <chr>    <chr> 
#   1 counts Ctrl   EML_V1 0.000461 0.0018 0.00046  ***      T-test
# 2 counts Ctrl   KIF5B  0.158    0.16   0.15756  ns       T-test
# 3 counts Ctrl   TFG    0.0311   0.062  0.03109  *        T-test
# 4 counts Ctrl   EML_V3 0.00381  0.011  0.00381  **       T-test

# Non-induced
compare_means(counts ~ fusion, data = gene_df_non_ind, method="t.test",ref.group = "Ctrl")
# .y.    group1 group2     p p.adj p.format p.signif method
# <chr>  <chr>  <chr>  <dbl> <dbl> <chr>    <chr>    <chr> 
#   1 counts Ctrl   EML_V1 0.543     1 0.54     ns       T-test
# 2 counts Ctrl   KIF5B  0.363     1 0.36     ns       T-test
# 3 counts Ctrl   TFG    0.363     1 0.36     ns       T-test
# 4 counts Ctrl   EML_V3 0.892     1 0.89     ns       T-test

# Create annotation for plot
annot_text_t_test_p <- data.frame(
  fusion_fullName = c(2:5),
  counts = rep(3000,4),
  condition = factor("Induced",levels = c("Induced","Non-induced"))
)

# Anova
#######

# Induced
anova_res_ind <- compare_means(counts ~ fusion, data = gene_df_ind, method="anova")
# .y.             p     p.adj p.format p.signif method
# <chr>       <dbl>     <dbl> <chr>    <chr>    <chr> 
#   1 counts 0.00000139 0.0000014 1.4e-06  ****     Anova 

# Non-induced
anova_res_nonInd <- compare_means(counts ~ fusion, data = gene_df_non_ind, method="anova")
# .y.        p p.adj p.format p.signif method
# <chr>  <dbl> <dbl> <chr>    <chr>    <chr> 
#   1 counts 0.477  0.48 0.48     ns       Anova 

# Create annotation for plot
annot_text_anov_p <- data.frame(
  fusion_fullName = c(2,2),
  counts = rep(3250,2),
  condition = factor(c("Induced","Non-induced"),levels = c("Induced","Non-induced"))
)

# Plot
p_SERPINB4 <- ggplot(gene_df, aes(x=fusion_fullName, y=counts)) +
  geom_boxplot(outlier.shape = NA, fill="grey", lwd=0.5) +
  geom_jitter(width = 0.25, size=1) +
  ylim(-2,3250) +
  xlab("") +
  ylab("SERPINB4 expr. (DESeq2 counts)") +
  scale_x_discrete(labels= rep(c("Control","EML4-ALK-V1","EML4-ALK-V3","KIF5B-ALK","TFG-ALK"),2)) + 
  facet_grid(~condition,scales = "free") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
    legend.position = "none",
    axis.text.y = element_text(size=4),
    axis.title.y = element_text(size = 5),
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    strip.background = element_rect(
      color="transparent", 
      fill="transparent"
    ),
    strip.text = element_text(
      size = 5
      )
  ) +
  geom_text(data = annot_text_t_test_p,
            label = t_test_res_ind$p.signif[c(1,4,2,3)],
            size = 2) +
  geom_text(data = annot_text_anov_p,
            label = paste0("Anova,p=",c(anova_res_ind$p.adj,anova_res_nonInd$p.adj)),
            size = 1)

# NSCLCL data
NSCLCL_logFC<- sapply(NSCLCL_res, function(x) x$log2FoldChange)
rownames(NSCLCL_logFC)<- NSCLCL_res[[1]]$HGNC
NSCLCL_SERPINB<- NSCLCL_logFC[grep("SERPINB",rownames(NSCLCL_logFC)),]
NSCLCL_SERPINB<- NSCLCL_SERPINB[,grep("6h",colnames(NSCLCL_SERPINB), invert = T)]
colnames(NSCLCL_SERPINB)<- gsub("\\.6h","",colnames(NSCLCL_SERPINB))

# Merge with ALK fusion data 
res_proc_all<- res_diff_expr$RNA
res_proc_all<- res_proc_all[grep("_induced_vs_Ctrl_induced",names(res_proc_all))]
fusion_logFC<- sapply(res_proc_all, function(x) x$log2FoldChange)
rownames(fusion_logFC)<- res_diff_expr$RNA[[1]]$HGNC
colnames(fusion_logFC)<- gsub("_induced_vs.*","",colnames(fusion_logFC))

# Adding H228 cell lines
H2228_res <- H2228_res[grep("SERPINB",H2228_res$HGNC),c("HGNC","log2FoldChange")]
H2228_res <- H2228_res[order(H2228_res$HGNC),]

# And H3122 cell lines
H3122_res <- H3122_res[grep("SERPINB",H3122_res$HGNC),c("HGNC","log2FoldChange")]
H3122_res <- H3122_res[order(H3122_res$HGNC),]

# Merge
fusion_SERPINB<- cbind(NSCLCL_SERPINB, 
                       fusion_logFC[rownames(NSCLCL_SERPINB),],
                       NSCLCL_res$log2FoldChange,
                       H2228_res$log2FoldChange,
                       H3122_res$log2FoldChange)
colnames(fusion_SERPINB)[c(16,17)] <- c("H2228","H3122")

# heatmap
pdf(file = "results/figs/manuscript_RNA_heatmap.pdf")
heatmap.2(fusion_SERPINB, scale = "none", col=bluered(75), trace="none", dendrogram = "both", xlab="", srtCol=45, margins = c(15,25), symkey=TRUE, key=T, keysize=1, key.title = "", key.xlab = "logFC", density.info = "none", cexCol = 0.5)
dev.off()

fusion_SERPINB[,c("EML_V1", "EML_V3", "KIF5B", "TFG")]
#           EML_V1      EML_V3      KIF5B        TFG
# SERPINB1  -0.9840396  0.61260540 -0.58313506  0.54586062
# SERPINB10  3.6756835  2.67759769  1.21135655  0.04631158
# SERPINB11  5.7887196  8.94096985  1.79727706  7.59644745
# SERPINB12  0.0000000  0.00000000  0.00000000  0.00000000
# SERPINB13  1.4223490  1.62689227  0.86433353  1.95416047
# SERPINB2   6.6779972  6.47459699  2.97494582  5.35456837
# SERPINB3   8.2589220 11.31374463  4.05158587  9.98829564
# SERPINB4  10.0248368 13.29989569  5.79858148 11.95326885
# SERPINB5  -3.7089736 -1.38533682 -5.72198582 -1.77254482
# SERPINB6  -0.5332282 -0.08032204 -0.03233352 -0.16621470
# SERPINB7   2.8723075  4.47639389  3.18109925  3.42769009
# SERPINB8   2.5412205  2.82390021  2.25245612  2.67745139
# SERPINB9   0.0000000  1.71281297  0.00000000  0.00000000

# Merge plots
#############

# Keep space for blot & general figure!

p_venn_GSEA<- plot_grid(
  NULL, p_venn, p_GSEA,   
  labels=c("a","","c"),
  label_size = 12,
  ncol = 3,
  rel_widths =  c(1,1,2)
)

p_SERPIN<- plot_grid(
  p_SERPINB4, NULL, NULL, NULL,
  labels=c("d","f","e"),
  label_size = 12,
  ncol = 2, 
  rel_widths =  c(1,1)
)

p_volc_SERPIN<- plot_grid(
  p_volc, p_SERPIN,
  labels="b",
  ncol = 2,
  rel_widths = c(1,1)
)

p<- plot_grid(
  p_venn_GSEA, p_volc_SERPIN,NULL,
  ncol = 1,
  rel_heights = c(1,3,4)
)
ggsave("results/figs/manuscript_RNA.pdf",  p, width = 178, height = 265, units = "mm")

# Postprocrocessing
# - Import PDF in inkscape
# - Add panel A: exp. outline
# - Venn: change label colors
# - Add panel F: blots
# - Add heatmap & redo labels


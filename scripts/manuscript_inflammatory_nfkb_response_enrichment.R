library(ggplot2)
library(ggrepel)
library(cowplot)

load("data/ALK_fusion_data.RData")
source("scripts/functions/helpers.R")

# Inflammtory response genes, as given in 
# "Inflammatory regulatory network mediated by the joint
# action of NF-kB, STAT3, and AP-1 factors is involved in
# many human cancers"
infl_responos_genes <- c("IL1A","IL1B","IL1R1","IL1R2","IL1RAP","IL1RL1","MYD88",
                         "IRAK2","NFKB1","NFKB2","IL6","LIF","OSMR","JAK2","STAT3",
                         "TNFSF10","TNFRSF10D","TNFRSF11B","TNFRSF21","ATF3","FOS",
                         "FOSL1","FOSL2","JUN","JUNB","MAP3K8","MAP4K4")

inflammatory_response_nfkb <- list(inflammatory_response_nfkb1 = infl_responos_genes,
                                   inflammatory_response_nfkb2 = infl_responos_genes)


# Get all results from lorlatinib treated cell lines
res_diff_all <- list("EML4-ALK-V1" = res_diff_expr$RNA$EML_V1_induced_vs_Ctrl_induced,
                     "EML4-ALK-V3" = res_diff_expr$RNA$EML_V3_induced_vs_Ctrl_induced,
                     "KIF5B-ALK" = res_diff_expr$RNA$KIF5B_induced_vs_Ctrl_induced,
                     "TFG-ALK" = res_diff_expr$RNA$TFG_induced_vs_Ctrl_induced)

# Do enrichment for all samples
GSEA_all <- list()
for (cl in names(res_diff_all)) {
  res_proc<- as.data.frame(res_diff_all[[cl]])
  res_proc<- res_proc[!is.na(res_proc$log2FoldChange),]
  genes_DE_ls<- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = res_proc$HGNC, th_logFC = 1.5, th_logP = 2, curve = 0.1)
  GSEA<- do_GSEA2(genes_retrieved = genes_DE_ls$DE, genes_all = rownames(res_proc), GSEA_db = inflammatory_response_nfkb, min_genes = 2, isList = T)
  GSEA_all[[cl]] <- as.numeric(GSEA$p[1])
}

GSEA <- as.data.frame(unlist(GSEA_all))
GSEA$q <- p.adjust(GSEA[,1],method = "fdr")
GSEA$sample <- gsub(".*\\.","",rownames(GSEA))
GSEA$sample <- factor(GSEA$sample,levels=GSEA[order(GSEA$q,decreasing = T),3])
p <- ggplot(GSEA, aes(x = sample, y = -log10(q))) +
  geom_bar(stat="identity", fill="darkblue") +
  coord_flip(expand = F) +
  ylab("-log10(Padj)") +
  xlab("") +
  theme(
    plot.title = element_text(hjust = 0.5, size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
  )

# Save
ggsave("results/figs/manuscript_inflammatory_response_enrichment_induced_fusions.pdf",p,height = 100,width = 120,units = "mm")


# Get all results from lorlatinib treated cell lines
res_diff_all <- list("CUTO8" = NSCLCL_res$CUTO8.24h_Lorlatinib,
                     "CUTO9" = NSCLCL_res$CUTO9.24h_Lorlatinib,
                     "CUTO29" = NSCLCL_res$CUTO29.24h_Lorlatinib,
                     "YU1077" = NSCLCL_res$YU1077.24h_Lorlatinib,
                     "H2228" = H2228_res,
                     "H3122" = H3122_res)

# Do enrichment for all samples
GSEA_all <- list()
for (cl in names(res_diff_all)) {
  res_proc<- as.data.frame(res_diff_all[[cl]])
  res_proc<- res_proc[!is.na(res_proc$log2FoldChange),]
  genes_DE_ls<- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = res_proc$HGNC, th_logFC = 1.5, th_logP = 2, curve = 0.1)
  GSEA<- do_GSEA2(genes_retrieved = genes_DE_ls$DE, genes_all = rownames(res_proc), GSEA_db = inflammatory_response_nfkb, min_genes = 2, isList = T)
  GSEA_all[[cl]] <- as.numeric(GSEA$p[1])
}

GSEA <- as.data.frame(unlist(GSEA_all))
GSEA$q <- p.adjust(GSEA[,1],method = "fdr")
GSEA$sample <- gsub(".*\\.","",rownames(GSEA))
GSEA$sample <- factor(GSEA$sample,levels=GSEA[order(GSEA$q,decreasing = T),3])
p <- ggplot(GSEA, aes(x = sample, y = -log10(q))) +
  geom_bar(stat="identity", fill="darkblue") +
  coord_flip(expand = F) +
  ylab("-log10(Padj)") +
  xlab("") +
  theme(
    plot.title = element_text(hjust = 0.5, size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
  )

# Save
ggsave("results/figs/manuscript_inflammatory_response_enrichment_CLs.pdf",p,height = 100,width = 120,units = "mm")


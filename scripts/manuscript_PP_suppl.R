library(WriteXLS)
library(ggrepel)
library(cowplot)
library(scales)
library(ggpubr)

# Load data
#############
load("data/ALK_fusion_data.RData")
source("scripts/functions/helpers.R")
SH2andPTBfromSuppExcel2 <- read_excel("data/SH2andPTBfromSuppExcel2.xlsx")

fusions <- c("EML_V1_induced_vs_Ctrl_induced",
             "TFG_induced_vs_Ctrl_induced",
             "KIF5B_induced_vs_Ctrl_induced",
             "EML_V1_induced_lorla_vs_Ctrl_induced",
             "TFG_induced_lorla_vs_Ctrl_induced",
             "KIF5B_induced_lorla_vs_Ctrl_induced")
real_names <- c(EML_V1_induced_vs_Ctrl_induced = "EML4-ALK-V1",
                TFG_induced_vs_Ctrl_induced = "TFG-ALK",
                KIF5B_induced_vs_Ctrl_induced = "KIF5B-ALK",
                EML_V1_induced_lorla_vs_Ctrl_induced = "EML4-ALK-V1+Lorlatinib",
                TFG_induced_lorla_vs_Ctrl_induced = "TFG-ALK+Lorlatinib",
                KIF5B_induced_lorla_vs_Ctrl_induced = "KIF5B-ALK+Lorlatinib")


# Volcano plots
###########
genes_to_label<- rownames(res_diff_expr$PP$EML_V1_induced_lorla_vs_EML_V1_induced)[c(grep("ALK",rownames(res_diff_expr$PP$EML_V1_induced_lorla_vs_EML_V1_induced)),
                                                                                     which(rownames(res_diff_expr$PP$EML_V1_induced_lorla_vs_EML_V1_induced) %in% c("NFKB2_Y77","NFKB1_Y81","STAT3_Y705")))]
p_volcs <- list()
for (f in fusions) {
  res_proc <- res_diff_expr$PP[[f]]
  p_volc <- plot_volcano(DE_results = res_proc, gene = genes_to_label, labelAll = T, labelCol = "black", useHyperbolicTH = T,logFC_cu = 1.5,curve = 0.1, plotTH=F, p_cu = 0.05, plot_nominal_p = F)
  p_volc <- p_volc +
    ggtitle(real_names[[f]]) +
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
    scale_x_continuous(name = "log2(Fold Change)", breaks = pretty_breaks(), limits = c(-max(abs(res_proc$log2FoldChange)),max(abs(res_proc$log2FoldChange)))) +
    scale_y_continuous(name = "-log10(Padj)", breaks = pretty_breaks(), ,limits = c(0,max(-log10(res_proc$padj)))) 
  p_volcs[[f]] <- p_volc
}


# ALK PP Boxplot
################
pp_all <- as.data.frame(normalized_counts$PP)
pp_ALK <- pp_all[grep("ALK",rownames(pp_all)),grep("dox",colnames(pp_all))]
pp_ALK$Site <- rownames(pp_ALK)
pp_ALK <- reshape2::melt(pp_ALK)
colnames(pp_ALK) <- c("Site","Fusion","pp")

# Logging values for better visualization
pp_ALK$pp <- log(pp_ALK$pp)

# Correcting sample names
pp_ALK$Fusion <- gsub("\\.lorla","+Lorlatinib",gsub("\\.dox","_Induced",gsub("TFG","TFG-ALK",gsub("KIF5B","KIF5B-ALK",gsub("V3","EML4-ALK-V3",gsub("V1","EML4-ALK-V1",gsub("\\.1|\\.2","",as.character(pp_ALK$Fusion))))))))

# Adding treatment
pp_ALK$Lorla <- ifelse(grepl("Lorla",pp_ALK$Fusion),"+","-")

# pp_ALK$Fusion <- gsub("_.*","",pp_ALK$Lorla)

# Adding comparison for t-test
my_comparisons <- list( c("Control_Induced", "Control_Induced+Lorlatinib"),
                        c("EML4-ALK-V1_Induced", "EML4-ALK-V1_Induced+Lorlatinib"),
                        c("EML4-ALK-V3_Induced", "EML4-ALK-V3_Induced+Lorlatinib"),
                        c("KIF5B-ALK_Induced", "KIF5B-ALK_Induced+Lorlatinib"),
                        c("TFG-ALK_Induced", "TFG-ALK_Induced+Lorlatinib") )

# Plot
p_ALK_boxplot <- ggplot(pp_ALK, aes(x = Fusion, y=pp,fill=Lorla)) +
  geom_boxplot(lwd=0.5,outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons,method = "t.test",size = 3,label = "p.signif",label.y = c(16,24,26,25,23)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    # legend.position = "none",
    axis.text.y = element_text(size=6),
    axis.title.y = element_text(size = 7),
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size=8),
    legend.key.size = unit(3, 'mm')
  ) +
  labs(x = "", y = "log(normalized counts)") +
  scale_fill_manual(values=c("dodgerblue3", "brown1")) +
  scale_x_discrete(labels= c("Control","","EML4-ALK-V1","","EML4-ALK-V3","","KIF5B-ALK","","TFG-ALK",""))
# Unable to group labels while keeping the p-value indicators, so will move axis text in Inkscape


p1 <- plot_grid(
  p_volcs$EML_V1_induced_vs_Ctrl_induced,
  p_volcs$EML_V1_induced_lorla_vs_Ctrl_induced,
  p_volcs$KIF5B_induced_vs_Ctrl_induced,
  p_volcs$KIF5B_induced_lorla_vs_Ctrl_induced,
  p_volcs$TFG_induced_vs_Ctrl_induced,
  p_volcs$TFG_induced_lorla_vs_Ctrl_induced,
  nrow = 3,
  ncol = 2
)

p2 <- plot_grid(
  p_ALK_boxplot,NULL,
  nrow = 1
)

p <- plot_grid(
  p1,p2,NA,
  nrow = 3,
  rel_heights = c(1,0.3,0.1)
)

# Save plot
ggsave("results/figs/manuscript_PP_suppl.pdf",p, width = 178, height = 270, units = "mm")


PTB <- SH2andPTBfromSuppExcel2$`IRS1-like PTB`[!is.na(SH2andPTBfromSuppExcel2$`IRS1-like PTB`)]
SH2 <- SH2andPTBfromSuppExcel2$`SH2 domains`[!is.na(SH2andPTBfromSuppExcel2$`SH2 domains`)]

res_diff_all <- list("EML4-ALK-V1" = res_diff_expr$PP$EML_V1_induced_lorla_vs_EML_V1_induced,
                     "EML4-ALK-V3" = res_diff_expr$PP$EML_V3_induced_lorla_vs_EML_V3_induced,
                     "KIF5B-ALK" = res_diff_expr$PP$KIF5B_induced_lorla_vs_KIF5B_induced,
                     "TFG-ALK" = res_diff_expr$PP$TFG_induced_lorla_vs_TFG_induced)

res_diff_all_df <- do.call(rbind,res_diff_all)
res_diff_all_df$Protein <- gsub("_.*","",gsub(".*\\.","",rownames(res_diff_all_df)))
res_diff_all_df$Fusion <- gsub("\\..*","",rownames(res_diff_all_df))
res_diff_all_df$Domain <- NA
res_diff_all_df$Domain[res_diff_all_df$Protein%in%SH2] <- "SH2"
res_diff_all_df$Domain[res_diff_all_df$Protein%in%PTB] <- "PTB"
res_diff_all_df$Domain[res_diff_all_df$Protein%in%"ALK"] <- "ALK"
res_diff_all_df$Domain[!(res_diff_all_df$Protein%in%c(SH2,PTB,"ALK"))] <- "Other"
res_diff_all_df$Site <- gsub(".*\\.","",rownames(res_diff_all_df))

# Check distribution of logFCs
hist(res_diff_all_df[,"log2FoldChange"])

fusions <- unique(res_diff_all_df$Fusion)
domains <- unique(res_diff_all_df$Domain)

SEs <- list()
means <- list()
for (fusion in fusions) {
  for (domain in domains) {
    SEs[[fusion]][[domain]] <- sd(res_diff_all_df[res_diff_all_df$Fusion==fusion & res_diff_all_df$Domain==domain,"log2FoldChange"])/sqrt(length(res_diff_all_df[res_diff_all_df$Fusion==fusion & res_diff_all_df$Domain==domain,"log2FoldChange"]))
    means[[fusion]][[domain]] <- mean(res_diff_all_df[res_diff_all_df$Fusion==fusion & res_diff_all_df$Domain==domain,"log2FoldChange"])
  }
}

df <- data.frame(SE = unlist(SEs),mean = unlist(means))
df$Fusion <- gsub("\\..*","",rownames(df))
df$Domain <- gsub(".*\\.","",rownames(df))
df$Fusion <- factor(x = df$Fusion,
                    levels = rev(unique(df$Fusion)),
                    ordered = T)

df$Domain <- factor(x = df$Domain,
                    levels = unique(df$Domain)[c(1,3,2,4)],
                    ordered = T)

df$ID <- paste0(df$Fusion,"_",df$Domain)

# do t-test
stat.test <- res_diff_all_df %>%
  group_by(Fusion) %>%
  t_test(log2FoldChange ~ Domain) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

stat.test
# Fusion            .y. group1 group2   n1   n2   statistic        df        p      p.adj p.adj.signif
# 1  EML4-ALK-V1 log2FoldChange    ALK  Other   12 8967  -7.1112529  11.00736 1.96e-05 4.7040e-04          ***
#   2  EML4-ALK-V1 log2FoldChange    ALK    PTB   12  100  -6.0855060  12.32448 4.87e-05 1.1688e-03           **
#   3  EML4-ALK-V1 log2FoldChange    ALK    SH2   12  185  -6.2072657  11.70441 5.06e-05 1.2144e-03           **
#   4  EML4-ALK-V1 log2FoldChange  Other    PTB 8967  100   3.5021573 100.13207 6.91e-04 1.6584e-02            *
#   5  EML4-ALK-V1 log2FoldChange  Other    SH2 8967  185   4.5248812 187.92330 1.07e-05 2.5680e-04          ***
#   6  EML4-ALK-V1 log2FoldChange    PTB    SH2  100  185  -0.1427306 202.55351 8.87e-01 1.0000e+00           ns
# 7  EML4-ALK-V3 log2FoldChange    ALK  Other   12 8967 -10.3485480  11.00849 5.21e-07 1.2504e-05         ****
#   8  EML4-ALK-V3 log2FoldChange    ALK    PTB   12  100  -9.3863139  12.14020 6.45e-07 1.5480e-05         ****
#   9  EML4-ALK-V3 log2FoldChange    ALK    SH2   12  185  -9.3554940  11.80049 8.34e-07 2.0016e-05         ****
#   10 EML4-ALK-V3 log2FoldChange  Other    PTB 8967  100   3.2263870 100.51266 2.00e-03 4.8000e-02            *
#   11 EML4-ALK-V3 log2FoldChange  Other    SH2 8967  185   4.3594089 187.98896 2.14e-05 5.1360e-04          ***
#   12 EML4-ALK-V3 log2FoldChange    PTB    SH2  100  185   0.3396609 227.18365 7.34e-01 1.0000e+00           ns
# 13   KIF5B-ALK log2FoldChange    ALK  Other   12 8967  -5.5946502  11.00251 1.61e-04 3.8640e-03           **
#   14   KIF5B-ALK log2FoldChange    ALK    PTB   12  100  -4.7703738  11.46098 5.18e-04 1.2432e-02            *
#   15   KIF5B-ALK log2FoldChange    ALK    SH2   12  185  -5.0737935  11.15733 3.43e-04 8.2320e-03           **
#   16   KIF5B-ALK log2FoldChange  Other    PTB 8967  100   5.3659961 100.09195 5.21e-07 1.2504e-05         ****
#   17   KIF5B-ALK log2FoldChange  Other    SH2 8967  185   5.9123347 189.94228 1.54e-08 3.6960e-07         ****
#   18   KIF5B-ALK log2FoldChange    PTB    SH2  100  185  -1.6299043 167.98456 1.05e-01 1.0000e+00           ns
# 19     TFG-ALK log2FoldChange    ALK  Other   12 8967  -4.2760051  11.00226 1.00e-03 2.4000e-02            *
#   20     TFG-ALK log2FoldChange    ALK    PTB   12  100  -3.6803724  11.34087 3.00e-03 7.2000e-02           ns
# 21     TFG-ALK log2FoldChange    ALK    SH2   12  185  -4.2438776  11.18443 1.00e-03 2.4000e-02            *
#   22     TFG-ALK log2FoldChange  Other    PTB 8967  100   4.5605620 100.32584 1.45e-05 3.4800e-04          ***
#   23     TFG-ALK log2FoldChange  Other    SH2 8967  185   0.1595076 188.55377 8.73e-01 1.0000e+00           ns
# 24     TFG-ALK log2FoldChange    PTB    SH2  100  185  -3.5889622 203.36658 4.16e-04 9.9840e-03           ** 

stat.test_subset <- as.data.frame(stat.test)
stat.test_subset <- stat.test_subset[stat.test_subset$group1=="Other"|stat.test_subset$group2=="Other",]
stat.test_subset$comparison <- paste0(stat.test_subset$group1,"_",stat.test_subset$group2)
stat.test_subset$comparison <- gsub("_Other|Other_","",stat.test_subset$comparison)
stat.test_subset$comparison <- paste0(stat.test_subset$Fusion,"_",stat.test_subset$comparison)
df <- merge(df,stat.test_subset,by.x = 5,by.y = 12,all = T)
df$Fusion <- df$Fusion.x
df$p.adj.signif[is.na(df$p.adj.signif)] <- ""

p_bar <- ggplot(df, aes(x=Fusion, y=mean, fill=Domain)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge(0.9),size = 0.3) +
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.3,size = 0.1,
                position=position_dodge(width = 1)) +
  coord_flip() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size=6.5),
    legend.text = element_text(size=5),
    legend.key.size = unit(5, 'mm'),
    legend.position="top",
    plot.margin = unit(c(0, 0, 0, -8), "pt"),
    legend.box.margin=margin(0,-10,-10,-10)
  ) +
  labs(x="", y="mean(log2FoldChange)") + 
  scale_fill_manual(values=c("#999999","#E69F00","#009E73","#56B4E9")) +
  geom_text(aes(label = p.adj.signif, x = Fusion, y = mean), hjust = 1, vjust=1, position = position_dodge(width=1),size = 2)

p_bar


#######

# Construct DB
db_PP <- list(SH2 = res_diff_all_df[res_diff_all_df$Domain=="SH2","Site"],
              PTB = res_diff_all_df[res_diff_all_df$Domain=="PTB","Site"])

fGSEA_RS_p <- list()
for(fusion in fusions) {
  # Get stat and do fGSEA
  stat<- get_GSEA_stat(res_diff_all[[fusion]], isMice = F, MGI_to_HGNC = MGI_to_HGNC)
  res_proc <- do_fGSEA(db = db_PP, stat = stat)
  for(domain in c("SH2","PTB")) {
    # Get padj
    p_pw <- signif(as.numeric(res_proc[res_proc$pathway==domain,"padj"]),2)
    # Generate plot
    fGSEA_RS_p[[fusion]][[domain]] <- plotEnrichment(db_PP[[domain]],stat,ticksSize=.1) +
      # ggtitle(paste0(fusion,": ",domain," (Padj=",p_pw,")")) +
      ggtitle("") +
      geom_line(size=.5, col="green") +
      theme(
        plot.title = element_text(hjust = 0.5, size=7, face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size=0.2),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.text = element_text(size=6),
        axis.title = element_text(size=7),
        plot.margin = unit(c(-10, 5, 0, 5), "pt")
      ) +
      scale_x_continuous(name = "Rank") +
      scale_y_continuous(name = "Enrichment Score") +
      annotate("text", x = 2000, y = -0.1, label = paste0(fusion,": ",domain,"\nPadj = ",p_pw),size = 2)
  }
}

p_fGSEA_RS <- plot_grid(
  fGSEA_RS_p$`EML4-ALK-V1`$SH2,
  fGSEA_RS_p$`EML4-ALK-V1`$PTB,
  fGSEA_RS_p$`EML4-ALK-V3`$SH2,
  fGSEA_RS_p$`EML4-ALK-V3`$PTB,
  fGSEA_RS_p$`KIF5B-ALK`$SH2,
  fGSEA_RS_p$`KIF5B-ALK`$PTB,
  fGSEA_RS_p$`TFG-ALK`$SH2,
  fGSEA_RS_p$`TFG-ALK`$PTB,
  nrow = 4
)

barplot_fGSEA <- plot_grid(
  p_bar,p_fGSEA_RS,
  nrow = 1,
  rel_widths = c(1,1)
)

sites_of_interest <- c("CBL_Y700","CRK_Y136","CRKL_Y207","DOK1_Y449","FRS2_Y306","PTPN11_Y542","NCK1_Y112")

df <- res_diff_all_df[res_diff_all_df$Site%in%sites_of_interest,]
p <- ggplot(df,aes(x = Site, y = log2FoldChange,fill = Site)) +
  geom_bar(stat="identity", color="black",size = 0.3,
           position=position_dodge(),width=0.9) +
  coord_flip() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size=6.5),
    legend.text = element_text(size=5),
    legend.key.size = unit(5, 'mm'),
    legend.position="none",
    plot.margin = unit(c(0, 0, 0, -8), "pt"),
    legend.box.margin=margin(0,-10,-10,-10),
    strip.background =element_rect(fill="white"),
    aspect.ratio = 1/1
  ) +
  scale_fill_manual(values=c("cyan4" ,"deepskyblue3", "dodgerblue4","#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")) 

p_bar_2<- p + facet_grid(~Fusion,scales = "free")

p_out <- plot_grid(
  barplot_fGSEA,p_bar_2,
  nrow = 2
)

# Save
ggsave("results/barplot_fGSEA_barplot.pdf",p_out, width = 180, height = 220, units = "mm")


samples <- c("EML_V1_induced_vs_Ctrl_induced",
             "EML_V3_induced_vs_Ctrl_induced",
             "KIF5B_induced_vs_Ctrl_induced",
             "TFG_induced_vs_Ctrl_induced",
             "EML_V1_induced_lorla_vs_Ctrl_induced",
             "EML_V3_induced_lorla_vs_Ctrl_induced",
             "KIF5B_induced_lorla_vs_Ctrl_induced",
             "TFG_induced_lorla_vs_Ctrl_induced")

DP_res <- list()
for (s in samples) {
  res_proc<- as.data.frame(res_diff_expr$PP[[s]])
  res_proc <- res_proc[!is.na(res_proc$log2FoldChange),]
  # temp_DP <- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = rownames(res_proc), th_logFC = 1.5, th_logP = -log10(0.05), curve = 0)
  # res_proc <- res_proc[rownames(res_proc) %in% temp_DP$DE,]
  # print(temp_DP$DE[grep("ALK_",temp_DP$DE)])
  DP_res[[s]] <- res_proc
  # print(s)
}

# Supplementary PP DEP result excel file
WriteXLS(x = DP_res,
         ExcelFileName = "results/tables/manuscript_PP_res_suppl_all.xlsx",
         SheetNames = c("EML4_ALK_V1","EML4_ALK_V3","KIF5B_ALK","TFG_ALK",
                        "EML4_ALK_V1_lorla","EML4_ALK_V3_lorla","KIF5B_ALK_lorla","TFG_ALK_lorla"),
         row.names = TRUE)

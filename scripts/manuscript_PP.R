library(ggrepel)
library(ggpubr)
library(ggvenn)
library(scales)
library(WriteXLS)
library(cowplot)
library(gridExtra)

load("data/ALK_fusion_data.RData")
source("scripts/functions/helpers.R")

# Volcano plots
###############
fusions <- c("EML_V3_induced_vs_Ctrl_induced",
             "EML_V3_induced_lorla_vs_Ctrl_induced")

real_names <- c("EML4-ALK-V1","EML4-ALK-V3","KIF5B-ALK","TFG-ALK")

genes_to_label<- rownames(res_diff_expr$PP$EML_V1_induced_lorla_vs_EML_V1_induced)[c(grep("ALK",rownames(res_diff_expr$PP$EML_V1_induced_lorla_vs_EML_V1_induced)),
                                                                                     which(rownames(res_diff_expr$PP$EML_V1_induced_lorla_vs_EML_V1_induced) %in% c("NFKB2_Y77","NFKB1_Y81","STAT3_Y705")))]
p_volcs <- list()
for (f in fusions) {
  res_proc <- res_diff_expr$PP[[f]]
  # res_proc <- res_diff_expr$PP$EML_V1_induced_lorla_vs_Ctrl_induced_lorla
  p_volc <- plot_volcano(DE_results = res_proc, gene = genes_to_label, labelAll = T, labelCol = "black", useHyperbolicTH = T,logFC_cu = 1.5,curve = 0.1, plotTH=F, p_cu = 0.05, plot_nominal_p = F)
  p_volc <- p_volc +
    ggtitle("") +
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


# Motif Pie charts
##################
conditions <- c("EML_V1_induced_vs_Ctrl_induced",
                "EML_V1_induced_lorla_vs_Ctrl_induced",
                "EML_V3_induced_vs_Ctrl_induced",
                "EML_V3_induced_lorla_vs_Ctrl_induced",
                "KIF5B_induced_vs_Ctrl_induced",
                "KIF5B_induced_lorla_vs_Ctrl_induced",
                "TFG_induced_vs_Ctrl_induced",
                "TFG_induced_lorla_vs_Ctrl_induced")                

real_names2 <- c(rep("EML4-ALK-V1",2),rep("EML4-ALK-V3",2),rep("KIF5B-ALK",2),rep("TFG-ALK",2))
pie_plots <- list()
data_other_ls <- list()
for (j in 1:length(conditions)) {
  res_proc <- res_diff_expr$PP[[conditions[j]]]
  sites_DP <- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = rownames(res_proc), th_logFC = 1.5, th_logP = -log10(0.05), curve = 0.1)
  res_proc$siteAA <- substr(gsub(".*_","",rownames(res_diff_expr$PP$EML_V3_induced_lorla_vs_Ctrl_induced)),1,1)
  DQ_df <- data.frame(isDP = rownames(res_proc) %in% sites_DP$DE,
                      isDown = rownames(res_proc) %in% sites_DP$down,
                      isUp = rownames(res_proc) %in% sites_DP$up,
                      isY = ifelse(res_proc$siteAA=="Y",TRUE,FALSE))
  
  # Contingency table:
  is_Y_isUp <- nrow(DQ_df[DQ_df$isUp==T & DQ_df$isY==T,])
  not_Y_isUp <- nrow(DQ_df[DQ_df$isUp==T & DQ_df$isY==F,])
  is_Y_notUp <- nrow(DQ_df[DQ_df$isUp==F & DQ_df$isY==T,])
  not_Y_notUp <- nrow(DQ_df[DQ_df$isUp==F & DQ_df$isY==F,])
  
  ft_tab <- matrix(c(is_Y_isUp,is_Y_notUp,
                     not_Y_isUp,not_Y_notUp),nrow = 2)
  
  p_val <- signif(fisher.test(ft_tab)$p.value,digits = 2)
  
  DQ_t_hyper<- table(DQ_df[,"isUp"],DQ_df$isY)
  data_hyper<- data.frame(n=DQ_t_hyper["TRUE",c("TRUE","FALSE")])
  data_hyper$prop<- data_hyper$n/sum(data_hyper$n)
  data_hyper$group<- c("Y", "Other")
  data_hyper$group_label<- paste0(data_hyper$group, " (n=", data_hyper$n, "; ",signif(100*data_hyper$prop,3), "%)")
  data_hyper$ypos<- cumsum(data_hyper$prop)- 0.5*data_hyper$prop
  
  data_other <- data.frame(n = c(is_Y_notUp,not_Y_notUp),
                           prop = c(is_Y_notUp/not_Y_notUp,1-is_Y_notUp/not_Y_notUp),
                           group_label = c("Y","notY"))
  data_other_ls[[conditions[j]]] <- data_other

  n_DP_hyper<- data_hyper["TRUE","n"]
  n_tot_hyper<- sum(data_hyper$n)
  prop_hyper<- 100*round(n_DP_hyper/n_tot_hyper,3)
  n_hyper_not_Y <- data_hyper["FALSE","n"]
  
  # Main pie
  p_pie_fusion <- ggplot(data_hyper, aes(x="", y=prop, fill=group_label)) +
    geom_bar(stat="identity", width=0.5, color="white") +
    coord_polar("y", start=0) +
    labs(
      caption = paste0(prop_hyper,"% (",n_DP_hyper,"/",n_tot_hyper,")",
                       "\n",
                       "P = ",p_val),
      title = real_names2[j]
    ) +
    theme_void() + 
    theme(
      plot.title = element_text(hjust = 0.5, size=8),
      plot.caption = element_text(hjust= 0.5, size=7),
      legend.title = element_blank(),
      legend.position = "none"
    ) + 
    scale_fill_manual(values=c("#999999", "#D55E00")) 
  
  # Background pie 
  p_pie_bg <- ggplot(data_other, aes(x="", y=prop, fill=group_label)) +
    geom_bar(stat="identity", width=0.5, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    labs(
      caption = paste0(data_other[1,2]*100,"%"),
      title = real_names2[j]
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size=8),
      plot.caption = element_text(hjust= 0.5, size=7),
      legend.title = element_blank(),
      legend.position = "none"
    ) + 
    scale_fill_manual(values=c("white", "black")) 
  
  pie_plots[[conditions[j]]] <- list(main = p_pie_fusion,bg = p_pie_bg)
}

p_pie <- grid.arrange(
  pie_plots$EML_V1_induced_vs_Ctrl_induced$main,
  pie_plots$EML_V1_induced_vs_Ctrl_induced$bg,
  pie_plots$EML_V1_induced_lorla_vs_Ctrl_induced$main,
  pie_plots$EML_V1_induced_lorla_vs_Ctrl_induced$bg,
  pie_plots$EML_V3_induced_vs_Ctrl_induced$main,
  pie_plots$EML_V3_induced_vs_Ctrl_induced$bg,
  pie_plots$EML_V3_induced_lorla_vs_Ctrl_induced$main,
  pie_plots$EML_V3_induced_lorla_vs_Ctrl_induced$bg,
  pie_plots$KIF5B_induced_vs_Ctrl_induced$main,
  pie_plots$KIF5B_induced_vs_Ctrl_induced$bg,
  pie_plots$KIF5B_induced_lorla_vs_Ctrl_induced$main,
  pie_plots$KIF5B_induced_lorla_vs_Ctrl_induced$bg,
  pie_plots$TFG_induced_vs_Ctrl_induced$main,
  pie_plots$TFG_induced_vs_Ctrl_induced$bg,
  pie_plots$TFG_induced_lorla_vs_Ctrl_induced$main,
  pie_plots$TFG_induced_lorla_vs_Ctrl_induced$bg,
  nrow = 4,
  ncol = 4
)

# Proportions of non-Y
data_other_ls

# Save pie plot individually
ggsave("results/figs/manuscript_PP_pie.pdf",  p_pie, width = 178, height = 100, units = "mm")


# # Logo plot
# ###########
logo_plot_list <- list()
for (i in 1:length(conditions)) {
  res_proc <- res_diff_expr$PP[[conditions[i]]]
  sites <- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = rownames(res_proc), th_logFC = 1.5, th_logP = -log10(0.05), curve = 0.1)
  logo_plot_df <- PP_motifs[,c(4,5)]
  logo_plot_df$isDE <- ifelse(logo_plot_df$ID%in%sites$DE,"Yes","No")
  logo_plot_df$isHyper <- ifelse(logo_plot_df$ID %in% sites$up,"Yes","No")
  logo_plot_df$isHypo <- ifelse(logo_plot_df$ID %in% sites$down,"Yes","No")
  logo_plot_df <- logo_plot_df[rev(order(logo_plot_df$isHyper)),]
  logo_plot_list[[conditions[i]]] <- logo_plot_df
}

WriteXLS(x = logo_plot_list,
         ExcelFileName = "results/tables/manuscript_logo_plot_data.xlsx",
         SheetNames = gsub("induced","ind",conditions),
         row.names = TRUE)

# Creating logo plot at http://kplogo.wi.mit.edu/submit.html? then adding them in figure in Inkscape

# ALK Tyr logFC Heatmap
#######################
PP_logFC <- res_diff_expr$PP[names(res_diff_expr$PP) %in% c("EML_V1_induced_lorla_vs_EML_V1_induced",
                                                            "EML_V3_induced_lorla_vs_EML_V3_induced",
                                                            "KIF5B_induced_lorla_vs_KIF5B_induced",
                                                            "TFG_induced_lorla_vs_TFG_induced")]

PP_logFC <- Reduce(function(dtf1, dtf2) cbind(dtf1, dtf2),PP_logFC)
PP_logFC <- PP_logFC[grep("ALK",rownames(PP_logFC)),grep("log2FoldChange",colnames(PP_logFC))]
colnames(PP_logFC) <- c("EML4-ALK-V1","EML4-ALK-V3","KIF5B-ALK","TFG-ALK")
PP_logFC$Tyr <- gsub("ALK_","",rownames(PP_logFC))
PP_logFC <- reshape2::melt(PP_logFC)
colnames(PP_logFC) <- c("Tyr","Fusion","log2FC")
PP_logFC$Tyr <- factor(PP_logFC$Tyr, levels=rev(levels(factor(PP_logFC$Tyr))))
level_order <- c("EML4-ALK-V1","EML4-ALK-V3","KIF5B-ALK","TFG-ALK")
PP_logFC$text_col <- ifelse(PP_logFC$log2FC<(-7.5),"white","black")

p_alk_tyr_heatmap <- ggplot(PP_logFC, aes(x = factor(Fusion, levels = level_order), Tyr, fill= log2FC),color = text_col) + 
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",high ="red", 
                       midpoint = 0) +
  geom_text(show.legend = FALSE,aes(label = round(log2FC, 2),color = text_col), size=2) + 
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
    legend.text = element_text(size=7),
    legend.key.size = unit(2.5, 'mm'),
    legend.title = element_text(size=8),
    legend.position="top",
    axis.ticks.x = element_blank(),
  ) + 
  labs(x = "", y = "") +
  scale_colour_manual(values=c("black", "white"))


# HypoPP/HyperPP boxplot
########################
conditions <- c("EML_V1_induced_vs_Ctrl_induced",
                "EML_V1_induced_lorla_vs_Ctrl_induced",
                "EML_V3_induced_vs_Ctrl_induced",
                "EML_V3_induced_lorla_vs_Ctrl_induced",
                "KIF5B_induced_vs_Ctrl_induced",
                "KIF5B_induced_lorla_vs_Ctrl_induced",
                "TFG_induced_vs_Ctrl_induced",
                "TFG_induced_lorla_vs_Ctrl_induced")

DP_sites <- c()
for (i in 1:length(conditions)) {
  sites_DP <- list()
  res_proc <- res_diff_expr$PP[[conditions[i]]]
  sites <- get_DE(logFC = res_proc$log2FoldChange, P = res_proc$padj, genes = rownames(res_proc), th_logFC = 1.5, th_logP = -log10(0.05), curve = 0.1)
  
  for (j in c("up","down")) {
    DP_sites[paste0(conditions[i],"_",j)] <- length(sites[[j]])
  }
}

data <- data.frame(HyperPP = DP_sites[grep("induced_vs_Ctrl_induced_up",names(DP_sites))],
                   HypoPP = DP_sites[grep("induced_vs_Ctrl_induced_down",names(DP_sites))],
                   HyperPP_plus_lorla = DP_sites[grep("lorla_vs_Ctrl_induced_up",names(DP_sites))],
                   HypoPP_plus_lorla = DP_sites[grep("lorla_vs_Ctrl_induced_down",names(DP_sites))])

data$Fusion <- gsub("KIF5B.*","KIF5B-ALK",gsub("TFG.*","TFG-ALK",gsub("EML_V3.*","EML4-ALK-V3",gsub("EML_V1.*","EML4-ALK-V1",rownames(data)))))
data <- reshape2::melt(data,id.vars='Fusion')
data$variable <- paste0(data$Fusion,"_",data$variable)

data$DE <- ifelse(grepl("Hyper",data$variable),"HyperPP","HypoPP")


p_bar <- ggplot(data, aes(variable, value)) +   
  geom_bar(aes(fill = DE),width=0.8, position=position_dodge(width=0.9), stat="identity") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=6),
    axis.title.y = element_text(size = 7),
    axis.line = element_line(colour = "black", size=0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.text = element_text(size=7),
    legend.key.size = unit(2.5, 'mm'),
    legend.title = element_blank(),
    legend.position="top",
    axis.ticks.x = element_blank(),
    strip.background = element_rect(colour="white", fill="white"))  +
  labs(x = "", y = "# significant sites") +
  scale_fill_manual(values=c("brown1","dodgerblue1")) +
  facet_grid(~ Fusion,scales = "free_x") 


# Merge plots
p1 <- plot_grid(
  NULL,p_bar,
  nrow = 1,
  rel_widths = c(0.6,1)
)

p2 <- plot_grid(
  p_volcs$EML_V3_induced_vs_Ctrl_induced,
  p_volcs$EML_V3_induced_lorla_vs_Ctrl_induced,
  nrow = 1
)

p3 <- plot_grid(
  NULL,p_alk_tyr_heatmap,NULL,
  rel_widths = c(0.2,1,2),
  nrow = 1
)

p <- plot_grid(
  p1,NA,p2,p3,
  nrow = 4,
  rel_heights = c(0.4,0.05,0.6,0.5)
)

# Save plot
ggsave("results/figs/manuscript_PP.pdf",  p, width = 178, height = 210, units = "mm")

# Postprocrocessing
# - Import PDF in inkscape
# - Add workflow cartoon 
# - Create logoplots from http://kplogo.wi.mit.edu/submit.html? and merge with volcano plots
# - Edit circle diagram to proper overlap


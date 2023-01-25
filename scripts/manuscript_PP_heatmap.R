library(ggplot2)
load("data/ALK_fusion_data.RData")

sites <- c("STAT3_Y705","FOSL1_S265","NFKB1_Y81","NFKB2_Y77","ALK_Y1604")
fusions <- c("EML_V1_induced_vs_Ctrl_induced",
             "EML_V3_induced_vs_Ctrl_induced",
             "KIF5B_induced_vs_Ctrl_induced",
             "TFG_induced_vs_Ctrl_induced",
             "EML_V1_induced_lorla_vs_Ctrl_induced",
             "EML_V3_induced_lorla_vs_Ctrl_induced",
             "KIF5B_induced_lorla_vs_Ctrl_induced",
             "TFG_induced_lorla_vs_Ctrl_induced")

PP_logFC <- res_diff_expr$PP[names(res_diff_expr$PP) %in% fusions]

PP_logFC <- do.call(rbind,PP_logFC)
PP_logFC$Site <- gsub(".*\\.","",rownames(PP_logFC))
PP_logFC$Fusion <- gsub("_vs.*","",rownames(PP_logFC))
colnames(PP_logFC)[1]<-"log2FC"
PP_logFC <- PP_logFC[PP_logFC$Site%in%sites,c("Site","Fusion","log2FC")]
PP_logFC$Site <- factor(PP_logFC$Site, levels=rev(levels(factor(PP_logFC$Site))))
level_order <- c("FOSL1_S265","STAT3_Y705","NFKB1_Y81","NFKB2_Y77","ALK_Y1604")
PP_logFC$text_col <- ifelse(PP_logFC$log2FC<(-7.5),"white","black")
PP_logFC$Fusion <- factor(PP_logFC$Fusion, levels=levels(factor(PP_logFC$Fusion)))
level_order2 <- c("EML_V1_induced","EML_V3_induced","KIF5B_induced","TFG_induced",
                  "EML_V1_induced_lorla","EML_V3_induced_lorla","KIF5B_induced_lorla","TFG_induced_lorla")


pp_heatmap <- ggplot(PP_logFC, aes(x = factor(Fusion,levels = level_order2), y = factor(Site,levels = level_order), fill= log2FC),color = text_col) + 
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
    legend.key.size = unit(4, 'mm'),
    legend.title = element_text(size=8),
    legend.position="top",
    axis.ticks.x = element_blank(),
    aspect.ratio = 0.9/1.3
  ) + 
  labs(x = "", y = "") +
  scale_colour_manual(values=c("black", "white")) + 
  scale_x_discrete(labels= c("EML4-ALK-V1","EML4-ALK-V3","KIF5B-ALK","TFG-ALK","EML4-ALK-V1","EML4-ALK-V3","KIF5B-ALK","TFG-ALK"))

ggsave("results/figs/manuscript_PP_sites_heatmap.pdf",pp_heatmap)

# Postprocessing
# - Import into Inkscape
# - Add labels for induced and non-induced


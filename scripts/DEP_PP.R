# Load sample_info
sample_info <- readRDS("data/sample_info_PP.rds")
PP_counts <- readRDS("data/ALK_fusion_PP_normalized_counts.rds")

# exp dessign
experimental_design<- data.frame(label=colnames(PP_counts), condition=sample_info[colnames(PP_counts),"cond"],replicate=sample_info[colnames(PP_counts),"rep"])
data<- as.data.frame(PP_counts)
data$id<- rownames(data)
data<- make_unique(data, names = "id", ids = "id") # Data unique already, but required step for DEP
data_se <- make_se(data, 1:ncol(PP_counts), experimental_design) # Create summarized experiment

# Missing values
data_filt <- filter_missval(data_se, thr = 0) # Filter for proteins that are identified in all replicates of at least one condition
data_imp<- data_filt 
data_diff <- test_diff(data_imp, type = "all") # Test all contrasts!
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))
data_results <- get_results(dep)

# Default format
rownames(data_results)<- data_results$name
data_results<- data_results[,grep("centered|significant|name|ID",colnames(data_results),invert = T)]
contrasts<- unique(gsub("_p\\.val|_p\\.adj|_ratio","",colnames(data_results)))
results_DE <- list()
for(cond in contrasts){
  # Contrast
  data_tmp<- data_results[,which(gsub("_ratio|_p.*","",colnames(data_results))==cond)]
  if(length(data_tmp)==0) next # Inverse to be added
  data_tmp<- data_tmp[,c(grep("ratio",colnames(data_tmp)),grep("p.val",colnames(data_tmp)),grep("p.adj",colnames(data_tmp)))]
  colnames(data_tmp)<- c("log2FoldChange", "pvalue", "padj") # Same format as DESeq2
  results_DE[[cond]]<- data_tmp
  # Inverse contrast
  contrast_tmp<- unlist(strsplit(cond, "_vs_"))
  cond<- paste0(contrast_tmp[2],"_vs_",contrast_tmp[1])
  data_tmp$log2FoldChange<- -1*data_tmp$log2FoldChange  
  results_DE[[cond]]<- data_tmp
  # }
}

# Save
saveRDS(results_DE,file="data/ALK_fusion_PP_diff_expression_results.rds")

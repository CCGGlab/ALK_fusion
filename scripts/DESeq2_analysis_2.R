library(DESeq2)

# Get gene annotation map, using gencode.v29.annotation from https://www.gencodegenes.org/
gencode29.annotation <- read.delim("downloads/gencode/gencode.v29.annotation.gtf", header=FALSE, comment.char="#", stringsAsFactors=FALSE)[,9]
gene_info_ls<- strsplit(gencode29.annotation,";")
gene_types<- gsub(" gene_type ","",sapply(gene_info_ls,function(x) x[grep("gene_type",x)]))
gene_names<- gsub(" gene_name ","",sapply(gene_info_ls,function(x) x[grep("gene_name",x)]))
ENSG_names<- gsub("gene_id ","",sapply(gene_info_ls,function(x) x[grep("gene_id",x)]))
gene_info_df<- cbind(HGNC_symbol=gene_names,gene_biotype=gene_types)
rownames(gene_info_df)<- ENSG_names
ENSG_HGNC_map_table<- gene_info_df
  
# Load sample information
sample_info<- readRDS(file="data/sample_info.rds")

# Get quantified files
sampleNames<- rownames(sample_info)
sampleFiles<- paste0("raw/gene_counts/",sampleNames,".gene_counts")

# Add condition from sample info and create sample table
sampleCondition <- as.factor(sample_info[sampleNames,"condition"])
sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles, condition=sampleCondition)

# DESeq2
#########
ddsHTSeq<- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, design=~condition)
colData(ddsHTSeq)$condition<- factor(colData(ddsHTSeq)$condition, levels=levels(sampleTable$condition))
dds<- DESeq(ddsHTSeq)

# Define contrasts
###################

fusions<- c("Ctrl", "EML_V1", "EML_V3", "KIF5B", "TFG")

# All comparisons
conds<- c(fusions, paste0(fusions,"_d"))
contrasts<- data.frame(cond=rep(conds,length(conds)), ref=rep(conds,each=length(conds)))

# Trivial comparisons excluded
contrasts<- contrasts[contrasts$cond!=contrasts$ref,]

# Diff expression analysis for each condition before/after induction
#####################################################################

DE_ls<- list()

for(i in 1:nrow(contrasts)){
  
  cat(i, " ")
  
  # DE
  res_tmp<- results(dds, contrast = c("condition", contrasts$cond[i], contrasts$ref[i]))
  
  # Only protein coding genes
  res_proc<- res_tmp[ENSG_HGNC_map_table[rownames(res_tmp),"gene_biotype"]=="protein_coding"&!is.na(ENSG_HGNC_map_table[rownames(res_tmp),"gene_biotype"]),]
  
  # Convert to HGNC names; No one on one mapping!!!! Multiple gene IDs --> Take the one that has (highest) expression data
  res_proc_HGNC_names<- ENSG_HGNC_map_table[rownames(res_proc),"HGNC_symbol"]
  res_proc$HGNC<- res_proc_HGNC_names
  res_proc<- res_proc[order(res_proc$HGNC,res_proc[,"baseMean"],decreasing=T),]
  res_proc<- res_proc[!duplicated(res_proc$HGNC),]
  res_proc$ENSG<- rownames(res_proc)
  rownames(res_proc)<- res_proc$HGNC
  
  # Sort on genenames
  res_proc<- res_proc[order(res_proc[,"HGNC"]),]
  
  # Redo padj for selected (coding) genes
  res_proc$padj<- p.adjust(res_proc$pvalue,"fdr")
  
  # Data frame
  res_proc<- as.data.frame(res_proc)
  
  # Rename contrasts to default format
  cond_tmp<- contrasts$cond[i]
  if(grepl("_d",cond_tmp)) cond_tmp<- gsub("_d","_induced",cond_tmp)
  else cond_tmp<- paste0(cond_tmp,"_nonInduced")
  ref_tmp<- contrasts$ref[i]
  if(grepl("_d",ref_tmp)) ref_tmp<- gsub("_d","_induced",ref_tmp)
  else ref_tmp<- paste0(ref_tmp,"_nonInduced")
  
  # Put in list
  DE_ls[[paste0(cond_tmp, "_vs_", ref_tmp)]]<- res_proc
} 

# Get normalized counts for each sample & save
###############################################
normalized_counts <- counts(dds, normalized=T)
normalized_counts<- normalized_counts[res_proc$ENSG,]
rownames(normalized_counts)<- res_proc$HGNC
normalized_counts<- as.data.frame(normalized_counts)

# Save all data
###############
saveRDS(DE_ls,file = "data/ALK_fusion_diff_expression_resuls.rds")
saveRDS(normalized_counts,file = "data/ALK_fusion_normalized counts.rds")
saveRDS(FPKMs,file = "data/ALK_fusion_FPKM.rds")

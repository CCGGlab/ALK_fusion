library(readxl)
library(WriteXLS)
library(DEP)

# Y
####

ALK_Y<- read_excel(path = "raw/PP/Gothenburg_Q231531_ST_Y_112420/Table1_Gothenburg_Q231531_pY_FINAL.xlsx",sheet = "Summary",col_names = TRUE,skip = 10)

# Get average normalized abundances
ALK_Y_counts<- data.frame(ALK_Y[,(max(grep("CS",colnames(ALK_Y)))+1):(max(grep("CS",colnames(ALK_Y)))+24)])
colnames(ALK_Y_counts)<- gsub("\\.\\.\\.[[:digit:]]+","",colnames(ALK_Y_counts))

# Get logFC
ALK_Y_logFC<- data.frame(ALK_Y[,121:135])
colnames(ALK_Y_logFC)<- gsub("\\.\\.\\.[[:digit:]]+","",colnames(ALK_Y_logFC))

# Get sites
aaPept<- NULL
for(i in 1:nrow(ALK_Y)){
  cat(i," ")
  aaPept_temp<- strsplit(as.character(ALK_Y[i,"Peptide"]),"\\*")
  aaPept_temp<- aaPept_temp[[1]][-length(aaPept_temp[[1]])]
  aaPept_temp<- substr(aaPept_temp,nchar(aaPept_temp),nchar(aaPept_temp))
  aaPept_temp<- paste(aaPept_temp,collapse = ", ")
  aaPept<- append(aaPept,aaPept_temp)
}
# table(aaPept) # Majority Y, seems ok

# gene + aa position (only take first one if multiple isoforms)
ALK_Y_id<- data.frame(gene=apply(ALK_Y,1,function(x) unlist(strsplit(x["Gene Name"],";"))[1]),aa=aaPept,site=apply(ALK_Y,1,function(x) unlist(strsplit(x["Site"],";"))[1]),ALK_Y_counts,ALK_Y_logFC)
ALK_Y_id<- ALK_Y_id[!is.na(ALK_Y_id[,"gene"]),]

# Split if multiple sites per peptide
ALK_Y_id2<- NULL
for(i in 1:nrow(ALK_Y_id)){
  cat(i," ")
  l_tmp<- length(unlist(strsplit(ALK_Y_id[i,"aa"],",")))
  if(l_tmp==0) next
  else ALK_Y_id_tmp<- as.data.frame(cbind(gene=ALK_Y_id[i,"gene"],aa=unlist(strsplit(ALK_Y_id[i,"aa"],",")),site=unlist(strsplit(ALK_Y_id[i,"site"],",")),ALK_Y_id[i,4:ncol(ALK_Y_id)]))
  ALK_Y_id2<- rbind(ALK_Y_id2,ALK_Y_id_tmp)
}
ALK_Y_id2$aa<- gsub(" *","",ALK_Y_id2$aa) # Remove blank spaces
ALK_Y_id2$site<- gsub(" *","",ALK_Y_id2$site) # Remove blank spaces

# Add numeric value for site
ALK_Y_id2$site_num<- as.numeric(gsub("§","",ALK_Y_id2$site))

# Remove aspecific sites 
ALK_Y_id2<- ALK_Y_id2[ALK_Y_id2$aa=="Y",]

# Median if multiple sites
ALK_Y_id2$id<- paste0(ALK_Y_id2$gene,"_",ALK_Y_id2$aa,ALK_Y_id2$site_num)
ALK_Y_q<- apply(ALK_Y_id2[,grep("\\.",colnames(ALK_Y_id2))],2,function(x) tapply(x,ALK_Y_id2$id,median))

# Finalize
ALK_Y_counts<- ALK_Y_q[,!grepl("\\.\\.\\.",colnames(ALK_Y_q))]
ALK_Y_FC<- ALK_Y_q[,grepl("\\.\\.\\.",colnames(ALK_Y_q))]

# ST
####

ALK_ST<- read_excel(path = "raw/PP/Gothenburg_Q231531_ST_Y_112420/Table3_Gothenburg_Q231531_STMix_FINAL.xlsx",sheet = "Summary",col_names = TRUE,skip = 10)

# Get average normalized abundances
ALK_ST_counts<- data.frame(ALK_ST[,(max(grep("CS",colnames(ALK_ST)))+1):(max(grep("CS",colnames(ALK_ST)))+24)])
colnames(ALK_ST_counts)<- gsub("\\.\\.\\.[[:digit:]]+","",colnames(ALK_ST_counts))

# Get logFC
ALK_ST_logFC<- data.frame(ALK_ST[,121:135])
colnames(ALK_ST_logFC)<- gsub("\\.\\.\\.[[:digit:]]+","",colnames(ALK_ST_logFC))

# Get sites
aaPept<- NULL
for(i in 1:nrow(ALK_ST)){
  cat(i," ")
  aaPept_temp<- strsplit(as.character(ALK_ST[i,"Peptide"]),"\\*")
  aaPept_temp<- aaPept_temp[[1]][-length(aaPept_temp[[1]])]
  aaPept_temp<- substr(aaPept_temp,nchar(aaPept_temp),nchar(aaPept_temp))
  aaPept_temp<- paste(aaPept_temp,collapse = ", ")
  aaPept<- append(aaPept,aaPept_temp)
}
# table(aaPept) # Majority S or T, seems ok

# gene + aa position (only take first one if multiple isoforms)
ALK_ST_id<- data.frame(gene=apply(ALK_ST,1,function(x) unlist(strsplit(x["Gene Name"],";"))[1]),aa=aaPept,site=apply(ALK_ST,1,function(x) unlist(strsplit(x["Site"],";"))[1]),ALK_ST_counts,ALK_ST_logFC)
ALK_ST_id<- ALK_ST_id[!is.na(ALK_ST_id[,"gene"]),]

# Split if multiple sites per peptide
ALK_ST_id2<- NULL
for(i in 1:nrow(ALK_ST_id)){
  cat(i," ")
  l_tmp<- length(unlist(strsplit(ALK_ST_id[i,"aa"],",")))
  if(l_tmp==0) next
  else ALK_ST_id_tmp<- as.data.frame(cbind(gene=ALK_ST_id[i,"gene"],aa=unlist(strsplit(ALK_ST_id[i,"aa"],",")),site=unlist(strsplit(ALK_ST_id[i,"site"],",")),ALK_ST_id[i,4:ncol(ALK_ST_id)]))
  ALK_ST_id2<- rbind(ALK_ST_id2,ALK_ST_id_tmp)
}
ALK_ST_id2$aa<- gsub(" *","",ALK_ST_id2$aa) # Remove blank spaces
ALK_ST_id2$site<- gsub(" *","",ALK_ST_id2$site) # Remove blank spaces

# Add numeric value for site
ALK_ST_id2$site_num<- as.numeric(gsub("§","",ALK_ST_id2$site))

# Remove aspecific sites 
ALK_ST_id2<- ALK_ST_id2[ALK_ST_id2$aa%in%c("S","T"),]

# Median if multiple sites
ALK_ST_id2$id<- paste0(ALK_ST_id2$gene,"_",ALK_ST_id2$aa,ALK_ST_id2$site_num)
ALK_ST_q<- apply(ALK_ST_id2[,grep("\\.",colnames(ALK_ST_id2))],2,function(x) tapply(x,ALK_ST_id2$id,median))

# Finalize
ALK_ST_counts<- ALK_ST_q[,!grepl("\\.\\.\\.",colnames(ALK_ST_q))]
ALK_ST_FC<- ALK_ST_q[,grepl("\\.\\.\\.",colnames(ALK_ST_q))]

# Merge
PP_counts<- as.data.frame(rbind(ALK_Y_counts,ALK_ST_counts))
colnames(PP_counts) <- gsub("3a","3",colnames(PP_counts))
PP_logFC<- rbind(ALK_Y_FC,ALK_ST_FC)


# Generate sample info file
###########################
sample_name<- colnames(PP_counts)

# fusion
fusion<- rep(NA,length(sample_name))
for(f in c("V1","KIF5B","V3","TFG","Control")){
  if(f %in% c("V1","V3")) fusion[grep(f,sample_name)]<- paste0("EML_",f)
  else if(f %in% c("Control")) fusion[grep(f,sample_name)]<- "Ctrl"
  else fusion[grep(f,sample_name)]<- f
} 

# isInduced
isInduced<- rep(F,length(sample_name))
isInduced[grep("\\.dox",sample_name)]<- T

# hasLorla
hasLorla<- rep(F,length(sample_name))
hasLorla[grep("\\.lorla",sample_name)]<- T

# Sample matrix
sample_matrix<- data.frame(fusion,isInduced,hasLorla)
rownames(sample_matrix)<- gsub("3a","3",sample_name)

# DE varname
DE_varname<- colnames(PP_logFC)

for (i in 1:length(DE_varname)) {
  if (grepl("\\.1\\.",DE_varname[i])) {
    DE_varname[i] <- gsub("\\.\\.\\.",".1_",paste0(strsplit(DE_varname[i],"\\.1")[[1]],collapse = ""))
  } else if (grepl("\\.2\\.",DE_varname[i])) {
    DE_varname[i] <- gsub("\\.\\.\\.",".2_",paste0(strsplit(DE_varname[i],"\\.2")[[1]],collapse = ""))
  } else {
    DE_varname[i] <- gsub("\\.\\.\\.","_",DE_varname[i])
  }
}

# DE reference
DE_ref<- sapply(strsplit(DE_varname,"_"),function(x) x[[2]])
DE_samplename<- sapply(strsplit(DE_varname,"_"),function(x) x[[1]])
DE_sample_matrix<- data.frame(DE_ref,DE_varname)
rownames(DE_sample_matrix)<- DE_samplename


# Merge
sample_matrix$DE_varname<- NA 
sample_matrix[rownames(DE_sample_matrix),"DE_varname"]<- DE_varname
sample_matrix$DE_ref<- NA 
sample_matrix[rownames(DE_sample_matrix),"DE_ref"]<- DE_ref

# Arrange sample info
sample_matrix$cond<- "nonInduced"
sample_matrix$cond[sample_matrix$isInduced]<- "induced"
colnames(sample_matrix)[colnames(sample_matrix)=="isLorla"]<- "hasLorla"
sample_matrix$cond[sample_matrix$hasLorla]<- paste(sample_matrix$cond[sample_matrix$hasLorla],"lorla",sep="_")
sample_matrix$cond<- paste(sample_matrix$fusion,sample_matrix$cond,sep="_")
sample_matrix$rep<- NA
for(c in sample_matrix$cond){
  s_tmp<- rownames(sample_matrix[sample_matrix$cond==c,])
  sample_matrix[s_tmp,"rep"]<- 1:length(s_tmp)
}

# Save
saveRDS(sample_matrix,file="data/sample_info_PP.rds")
saveRDS(PP_counts,file="data/ALK_fusion_PP_normalized_counts.rds")


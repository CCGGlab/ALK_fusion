library(protr)
library(readxl)
# Getting Accession names for PP data
ALK_Y<- readxl::read_excel(path = "raw/PP/Gothenburg_Q231531_ST_Y_112420/Table1_Gothenburg_Q231531_pY_FINAL.xlsx",sheet = "Summary",col_names = TRUE,skip = 10)
ALK_Y <- ALK_Y[,grep("Accession|Gene Name",colnames(ALK_Y))]
ALK_ST<- readxl::read_excel(path = "raw/PP/Gothenburg_Q231531_ST_Y_112420/Table3_Gothenburg_Q231531_STMix_FINAL.xlsx",sheet = "Summary",col_names = TRUE,skip = 10)
ALK_ST <- ALK_ST[,grep("Accession|Gene Name",colnames(ALK_ST))]

PP_IDs1 <- rbind(ALK_Y,ALK_ST)
PP_IDs1 <- PP_IDs1[!is.na(PP_IDs1$`Gene Name`) & !is.na(PP_IDs1$Accession),]
PP_IDs1$`Gene Name` <- gsub(";.*","",PP_IDs1$`Gene Name`)
PP_IDs1$Accession <- gsub(";.*","",PP_IDs1$Accession)

PP_counts <- readRDS("data/ALK_fusion_PP_normalized_counts.rds")

PP_IDs2 <- data.frame(Gene_Name = gsub("_.*","",rownames(PP_counts)),
                      Site = gsub(".*_[A-Z]","",rownames(PP_counts)),
                      ID = rownames(PP_counts))

data <- merge(PP_IDs1,PP_IDs2,by.x = 1,by.y = 1)
data <- data[!duplicated(data$ID),]
data$seq5 <- rep(NA,nrow(data))


for (i in 1028:nrow(data)) {
  if (i==1027) next
  uni_seq_temp <- unlist(getUniProt(data$Accession[i]))
  site_temp <- as.numeric(data$Site[i])
  
  seq_dn <- substr(uni_seq_temp,site_temp-5,site_temp-1)
  seq_up <- substr(uni_seq_temp,site_temp+1,site_temp+5)
  
  # Add Xs at beginning or end of seq if cut by start or end of AA seq
  if (nchar(seq_dn)<5) seq_dn <- paste0(paste0(rep("X",5-nchar(seq_dn)),collapse = ""),seq_dn)
  if (nchar(seq_up)<5) seq_up <- paste0(seq_up,paste0(rep("X",5-nchar(seq_up)),collapse = ""))
  
  data$seq5[i] <- paste0(seq_dn,substr(uni_seq_temp,site_temp,site_temp),seq_up,collapse = "")
  print(i)
}

# Save 
saveRDS(data,"data/ALK_fusion_PP_motifs.rds")


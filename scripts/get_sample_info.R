# Create sample info file for RNA-Seq
sample_info_tmp <- data.frame(
  fusion = c("EML_V1","EML_V1","EML_V1","EML_V1","EML_V1","EML_V1",
             "KIF5B","KIF5B","KIF5B","KIF5B","KIF5B","KIF5B","EML_V1",
             "EML_V1","EML_V1","EML_V1","EML_V1","EML_V1","KIF5B",
             "KIF5B","KIF5B","KIF5B","KIF5B","KIF5B","TFG","TFG",
             "TFG","TFG","TFG","TFG","TFG","TFG","TFG","TFG","TFG",
             "TFG","EML_V3","EML_V3","EML_V3","EML_V3","EML_V3",
             "EML_V3","EML_V3","EML_V3","EML_V3","EML_V3","EML_V3",
             "EML_V3","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl",
             "Ctrl","Ctrl","Ctrl","Ctrl","Ctrl"),
  colony = c("C1","C1","D1","D1","F12","F12","E2","E2","E5","E5",
             "G1","G1","C1","C1","D1","D1","F12","F12","E2","E2",
             "E5","E5","G1","G1","F11","F11","C7","C7","D4","D4",
             "F11","F11","C7","C7","D4","D4","A3","A3","B8","B8",
             "C2","C2","A3","A3","B8","B8","C2","C2","A9","A9",
             "B11","B11","C3","C3","A9","A9","B11","B11","C3","C3"),
  isInduced = FALSE,
  condition = c("EML_V1","EML_V1_d","EML_V1","EML_V1_d","EML_V1",
                "EML_V1_d","KIF5B","KIF5B_d","KIF5B","KIF5B_d",
                "KIF5B","KIF5B_d","EML_V1","EML_V1_d","EML_V1",
                "EML_V1_d","EML_V1","EML_V1_d","KIF5B","KIF5B_d",
                "KIF5B","KIF5B_d","KIF5B","KIF5B_d","TFG","TFG_d",
                "TFG","TFG_d","TFG","TFG_d","TFG","TFG_d","TFG",
                "TFG_d","TFG","TFG_d","EML_V3","EML_V3_d","EML_V3",
                "EML_V3_d","EML_V3","EML_V3_d","EML_V3","EML_V3_d",
                "EML_V3","EML_V3_d","EML_V3","EML_V3_d","Ctrl",
                "Ctrl_d","Ctrl","Ctrl_d","Ctrl","Ctrl_d","Ctrl",
                "Ctrl_d","Ctrl","Ctrl_d","Ctrl","Ctrl_d"),
  hasLorla = FALSE,
  row.names = c("V1C1_1","V1C1_1_d","V1D1_1","V1D1_1_d","V1F12_1",
                "V1F12_1d","KIF5BE2_1","KIF5BE21d","KIF5BE5_1",
                "KIF5BE51d","KIF5BG1_1","KIF5BG11d","V1C1_2",
                "V1C1_2d","V1D1_2","V1D1_2d","V1F12_2","V1F12_2d",
                "KIF5BE2_2","KIF5BE22d","KIF5BE5_2","KIF5BE52d",
                "KIF5BG1_2","KIF5BG12d","TFGF11_1","TFGF111d",
                "TFGC7_1","TFGC7_1d","TFGD4_1","TFGD41d","TFGF11_2",
                "TFGF112d","TFGC7_2","TFGC7_2d","TFGD4_2","TFGD4_2d",
                "V3aA3_1","V3aA3_1d","V3aB8_1","V3aB8_1d","V3aC2_1",
                "V3aC2_1d","V3aA3_2","V3aA3_2d","V3aB8_2d","V3aB8_2", # V3aB8_2d and V3aB8_2 are swapped in the wet lab step, correcting for this
                "V3aC2_2","V3aC2_2d","contA9_1","contA9_1d","contB11_1",
                "contB111d","contC3_1","contC3_1d","contA9_2",
                "contA92d","contB11_2","contB112d","contC3_2","contC32d")
)

# d for doxycycline
sample_info_tmp$isInduced[grepl("_d",sample_info_tmp$condition)] <- TRUE

# Save 
saveRDS(sample_matrix,file="data/sample_info.rds")



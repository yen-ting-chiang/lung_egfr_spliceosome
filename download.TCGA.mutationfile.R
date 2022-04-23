library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
update.packages("TCGAbiolinks")

download.TCGA.mutationfile <- function(project_name, 
                     hugo_gene_name){
  #TCGAbiolinks(TCGA mutation data download)-------------------------
  
  query <- GDCquery(sprintf("TCGA-%s", project_name),
                           data.category = "Simple Nucleotide Variation",
                           workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
                           data.type = "Masked Somatic Mutation",
                           legacy = FALSE)
  GDCdownload(query,
              method = "api",
              files.per.chunk = 10)
  
  varscan2.maf <- GDCprepare(query,
                            save = TRUE,
                            save.filename = sprintf("%s_mutation.rdata",
                                                    project_name))
  
  varscan2.maf_target_gene <- varscan2.maf %>% 
    filter(Hugo_Symbol==sprintf("%s", hugo_gene_name)) %>%
    dplyr::select(Tumor_Sample_Barcode,
                  Protein_position,
                  HGVSp_Short,
                  Exon_Number,
                  Variant_Classification,
                  Consequence,
                  IMPACT,
                  SIFT,
                  PolyPhen
    )
  
  #create list of all case IDs with mutation data-------------------------------
  varscan2.maf_case_list=varscan2.maf %>%
    dplyr::select(Tumor_Sample_Barcode)
  
  varscan2.maf_case_list_unique <- varscan2.maf_case_list
  varscan2.maf_case_list_unique=unique(varscan2.maf_case_list_unique)
  # write.csv(PAAD.varscan2.maf_case_list_unique, 
  # file = "PAAD.varscan2.maf_case_list_unique.csv")
  
  
  #create the WT sample Coldata -----------
  
  list_of_samples_containing_MT_target_gene <- varscan2.maf_target_gene[,1, drop = TRUE]
  list_of_samples_containing_MT_target_gene_unique <- 
    unique(list_of_samples_containing_MT_target_gene)
  
  mutation_list <- varscan2.maf_case_list_unique[,1, drop = TRUE]
  DF1 <- data.frame(mutation_list)
  DF1.1 <- DF1 %>% filter(mutation_list %in% 
                            list_of_samples_containing_MT_target_gene_unique == FALSE) %>% 
    mutate(Protein_position = "WT", HGVSp_Short = "WT", IMPACT = "WT")
  colnames(DF1.1)[1] <- "Tumor_Sample_Barcode"
  
  #create the multiple-mutation sample Coldata -----------
  
  DF2 <- data.frame(mutation_list)
  DF2.1 <- varscan2.maf_target_gene %>% 
    filter(duplicated(Tumor_Sample_Barcode) == TRUE)%>% 
    dplyr::select(Tumor_Sample_Barcode, Protein_position, HGVSp_Short, IMPACT)%>% 
    mutate(Protein_position = "multiple", 
           HGVSp_Short = "multiple", 
           IMPACT = "multiple")
  DF2.1 <- unique(DF2.1)
  
  #create the single-mutation sample Coldata -----------
  DF3 <-varscan2.maf_target_gene%>% 
    dplyr::select(Tumor_Sample_Barcode, Protein_position, HGVSp_Short, IMPACT) %>% 
    filter(Tumor_Sample_Barcode %in% DF2.1[,1, drop = TRUE] == FALSE)
  
  
  #combine ColDatas----------------------------------------------------
  DF4 <- full_join(DF1.1, DF2.1, 
                   by = c("Tumor_Sample_Barcode", 
                          "Protein_position", 
                          "HGVSp_Short", 
                          "IMPACT"))
  
  DF5 <- full_join(DF4, DF3, 
                   by = c("Tumor_Sample_Barcode", 
                          "Protein_position", 
                          "HGVSp_Short", 
                          "IMPACT"))
  
  save(DF5, file = sprintf("%s.%s.ColData.rdata", 
                               project_name, 
                               hugo_gene_name))
  write.csv(DF5, file = sprintf("%s.%s.ColData.csv", 
                                    project_name,
                                    hugo_gene_name))
  
  
  
}

download.TCGA.mutationfile (project_name = "LUAD", 
         hugo_gene_name = "EGFR")

ColData <- read.csv("LUAD.EGFR.ColData.csv")
Batch <- read.delim("BatchData.tsv")
ColData <- ColData[,c(2:5)]

Batch <- Batch %>% 
  mutate(joinID = gsub("\\S-\\S\\S\\S\\S-\\S\\S$", 
                               "",
                               Sample))
ColData <- ColData %>% 
  mutate(joinID = gsub("\\S-\\S\\S\\S\\S-\\S\\S$", 
                       "",
                       Tumor_Sample_Barcode))

join_data <- left_join(Batch,
                       ColData, 
                        by= c("joinID" = "joinID"))

join_data_selected <- join_data[,c(1,7,8,9)]

join_data_selected <- join_data_selected %>% distinct()

write.csv(join_data_selected,"join_data_selected.csv")

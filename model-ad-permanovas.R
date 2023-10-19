library(tidyverse)
library(readxl)
library(vegan)

setwd("/media/julio/Storage/mAD_metabolomes/new_lipidomics/")

#SampleIDs
SampleIDs = read_excel("u54-ad-sample-list.xlsx", 
                       col_types = c("text", "text", "skip", 
                                     "skip", "skip", "skip"))

#Data
b1_b2_lipid_df = as.data.frame(t(read.csv("/media/julio/Storage/mAD_metabolomes/new_lipidomics/u54_ad_batch_corrected_data_TL.csv", row.names = 1, check.names = F))) %>%
  rownames_to_column(var = "name") %>% 
  left_join(SampleIDs, by = join_by("name" == "tm_id"))

#Replace NA for zeros
b1_b2_lipid_df2 = b1_b2_lipid_df %>% 
  replace(is.na(.), 0) %>%
  select(!name) %>%
  column_to_rownames(var = "sample_id")

#Import and format metadata.
b1_metadata = read_excel("MODEL-AD MASTER v. 2.xlsx", 
                         sheet = "U54 Mouse Samples", skip = 2) %>%
  select(`Sample ID`, Group, Sex, Genotype) %>%
  distinct()
b1_metadata$`Sample ID` =  as.character(b1_metadata$`Sample ID`)

b2_metadata = read_excel("MODEL-AD MASTER v. 2.xlsx", 
                         sheet = "NEW Samples") %>%
  select(`Sample ID`, Group, Sex, Genotype) %>%
  distinct()

metadata = bind_rows(b1_metadata, b2_metadata) %>% 
  mutate(Age = case_when(str_detect(Group, "4") ~ "4",
                         str_detect(Group, "12") ~ "12",
                         str_detect(Group, "18") ~ "18"),
         Sex = case_when(Sex %in% c("M", "F") ~ Sex,
                         str_detect(Sex, "Male") ~ "M",
                         str_detect(Sex, "Female") ~ "F"),
         Batch = case_when(`Sample ID` %in% b1_metadata$`Sample ID` ~ "One",
                           `Sample ID` %in% b2_metadata$`Sample ID` ~ "Two")
  )

#Cage data
Cage_IDs = read_excel("Julio's 5xFAD GWAS housing data 061423.xlsx") %>%
  select(Name,`Housing ID`, `Prior Housing ID`, Genotype) %>%
  mutate(CageID = case_when(
    !is.na(`Housing ID`) ~ `Housing ID`,
    is.na(`Housing ID`) ~ `Prior Housing ID`
  )) %>%
  select(Name, CageID) %>%
  right_join(metadata, by = join_by("Name" == `Sample ID`))

metadata = Cage_IDs %>%
  filter(!str_detect(Genotype, "PICALM"), !str_detect(Group, "CC 13")) %>%
  mutate(AD_path = case_when(
    str_detect(Genotype, "(?i)hemi") ~ "Yes",
    str_detect(Genotype, "(?i)het") ~ "Yes",
    str_detect(Genotype, "WT") ~ "No",
    str_detect(Genotype, "noncarrier") ~ "No"),
    Background = case_when(
      str_detect(Group, "(?i)TREM") ~ "TREM2",
      !str_detect(Group, "(?i)TREM") ~ "B6J"),
    Name = str_replace_all(Name, "TMF", "")
  )

metadata = metadata %>% 
  unite("Group2", AD_path, Background, remove = F) %>%
  mutate(Group2 = str_replace(Group2, "Yes", "HEMI 5xFAD"), Group2 = str_replace(Group2, "_B6J", ""), Group2 = str_replace(Group2, "No", "WT 5xFAD"))

#Merge data and metadata
combined_merged_data = b1_b2_lipid_df2 %>%
  rownames_to_column(var = "Name") %>%
  inner_join(metadata) %>%
  column_to_rownames(var = "Name") %>%
  filter(complete.cases(.))

combined_merged_data$CageID = as.factor(combined_merged_data$CageID)

#5xFAD vs WT lipidomics
perma_df1 = combined_merged_data %>% filter(Age == 4, Group2 %in% c("WT 5xFAD", "HEMI 5xFAD"))
perma_res1 = adonis2(select(perma_df1, !CageID:Background) ~ Group2 * Sex, 
        data = perma_df1, method = "euclidian")

perma_df2 = combined_merged_data %>% filter(Age == 12, Group2 %in% c("WT 5xFAD", "HEMI 5xFAD"))
perma_res2 = adonis2(select(perma_df2, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df2, method = "euclidian")

perma_df3 = combined_merged_data %>% filter(Age == 18, Group2 %in% c("WT 5xFAD", "HEMI 5xFAD"))
perma_res3 = adonis2(select(perma_df3, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df3, method = "euclidian")

#TREM2 vs WT lipidomics
perma_df4 = combined_merged_data %>% filter(Age == 4, Group2 %in% c("WT 5xFAD", "WT 5xFAD_TREM2"))
perma_res4 = adonis2(select(perma_df4, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df4, method = "euclidian")

perma_df5 = combined_merged_data %>% filter(Age == 12, Group2 %in% c("WT 5xFAD", "WT 5xFAD_TREM2"))
perma_res5 = adonis2(select(perma_df5, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df5, method = "euclidian")

#5xFAD*TREM2 vs 5xFAD
perma_df6 = combined_merged_data %>% filter(Age == 4, Group2 %in% c("HEMI 5xFAD_TREM2", "HEMI 5xFAD"))
perma_res6 = adonis2(select(perma_df6, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df6, method = "euclidian")

perma_df7 = combined_merged_data %>% filter(Age == 12, Group2 %in% c("HEMI 5xFAD_TREM2", "HEMI 5xFAD"))
perma_res7 = adonis2(select(perma_df7, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df7, method = "euclidian")

#Targeted metabolomics

#SampleIDs
tar_SampleIDs = read_excel("u54-ad-sample-list_tar.xlsx", 
                       col_types = c("text", "text", "skip", 
                                     "skip", "skip", "skip"))

# Batch 1
b1_tar_df = read_excel("u54-ad-normalized-data_tar.xlsx", sheet = "u54") %>%
  pivot_longer(!Sample)

b2_tar_df = read_excel("u54-ad-normalized-data_tar.xlsx", sheet = "ad") %>%
  pivot_longer(!Sample)

b1_b2_tar_df = bind_rows(b1_tar_df, b2_tar_df)

#Replace sampleIDs
b1_b2_tar_df = b1_b2_tar_df %>% 
  left_join(SampleIDs, by = join_by("name" == "tm_id")) %>%
  filter(complete.cases(.))

#Combine batch 1 & batch 2
b1_b2_tar_df2 = b1_b2_tar_df %>% 
  group_by(Sample, sample_id) %>% 
  summarise(value = sum(value)) %>%
  pivot_wider(id_cols = Sample, names_from = sample_id, values_from = value) %>%
  replace(is.na(.), 0)

combined_merged_data2 = as.data.frame(t(column_to_rownames(b1_b2_tar_df2, var = "Sample"))) %>%
  rownames_to_column(var = "Name") %>%
  inner_join(metadata)

#write.table(x = select(combined_merged_data2, !CageID:Background), file = "Combined_targeted_mtx_metaboanalyst.txt", quote = F, sep = "\t", row.names = F, col.names = T)
#write.table(x = select(combined_merged_data2, Name, CageID:Background), file = "Tar_mtx_metadata_metaboanalyst.txt", quote = F, sep = "\t", row.names = F, col.names = T)

normalized_data2 = combined_merged_data2 %>%
  filter(!Name %in% c("9379", "8392", "7772", "8393", "8123")) %>%
  column_to_rownames(var = "Name") %>%
  select(!CageID:Background) 

merged_norm_data2 = inner_join(metadata, rownames_to_column(normalized_data2, var = "Name")) %>%
  column_to_rownames(var = "Name") %>%
  filter(complete.cases(.))

merged_norm_data2$CageID = as.factor(merged_norm_data2$CageID)

#5xFAD vs WT tar metabolomics
perma_df8 = merged_norm_data2 %>% filter(Age == 4, Group2 %in% c("WT 5xFAD", "HEMI 5xFAD"))
perma_res8 = adonis2(select(perma_df8, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df8, method = "euclidian")

perma_df9 = merged_norm_data2 %>% filter(Age == 12, Group2 %in% c("WT 5xFAD", "HEMI 5xFAD"))
perma_res9 = adonis2(select(perma_df9, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df9, method = "euclidian") #6.4% batch r^2

perma_df10 = merged_norm_data2 %>% filter(Age == 18, Group2 %in% c("WT 5xFAD", "HEMI 5xFAD"))
perma_res10 = adonis2(select(perma_df10, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df10, method = "euclidian")

#TREM2 vs WT targeted metabolomics
perma_df11 = merged_norm_data2 %>% filter(Age == 4, Group2 %in% c("WT 5xFAD", "WT 5xFAD_TREM2"))
perma_res11 = adonis2(select(perma_df11, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df11, method = "euclidian") #6.9% batch variance

perma_df12 = merged_norm_data2 %>% filter(Age == 12, Group2 %in% c("WT 5xFAD", "WT 5xFAD_TREM2"))
perma_res12 = adonis2(select(perma_df12, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df12, method = "euclidian") #4.6% batch effect

#5xFAD*TREM2 vs 5xFAD
perma_df13 = merged_norm_data2 %>% filter(Age == 4, Group2 %in% c("HEMI 5xFAD_TREM2", "HEMI 5xFAD"), Sex == "F")
perma_res13 = adonis2(select(perma_df13, !CageID:Background) ~ Group2, 
                     data = perma_df13, method = "euclidian") #5.9 % batch effects

perma_df14 = merged_norm_data2 %>% filter(Age == 12, Group2 %in% c("HEMI 5xFAD_TREM2", "HEMI 5xFAD"), Sex == "F")
perma_res14 = adonis2(select(perma_df14, !CageID:Background) ~ Group2, 
                     data = perma_df14, method = "euclidian") #5.3 % variance

#5xFAD vs WT lipidomics with all time points
perma_df15 = merged_norm_data %>% filter(Group2 %in% c("WT 5xFAD", "HEMI 5xFAD"))
perma_res15 = adonis2(select(perma_df15, !CageID:Background) ~ Group2 * Sex, 
                     data = perma_df15, method = "euclidian", strata = perma_df15$Age)

#WT 5xFAD vs. WT 5xFAD_TREM2 lipidomics with all time points
perma_df16 = merged_norm_data %>% filter(Group2 %in% c("WT 5xFAD", "WT 5xFAD_TREM2"))
perma_res16 = adonis2(select(perma_df16, !CageID:Background) ~ Group2 * Sex, 
                      data = perma_df16, method = "euclidian", strata = perma_df16$Age)

#HEMI 5xFAD_TREM2", "HEMI 5xFAD lipidomics with all time points
perma_df17 = merged_norm_data %>% filter(Group2 %in% c("HEMI 5xFAD_TREM2", "HEMI 5xFAD"))
perma_res17 = adonis2(select(perma_df17, !CageID:Background) ~ Group2 * Sex, 
                      data = perma_df17, method = "euclidian", strata = perma_df17$Age)

#5xFAD vs WT lipidomics with all time points
perma_df18 = merged_norm_data2 %>% filter(Group2 %in% c("WT 5xFAD", "HEMI 5xFAD"))
perma_res18 = adonis2(select(perma_df18, !CageID:Background) ~ Group2 * Sex, 
                      data = perma_df18, method = "euclidian", strata = perma_df18$Age)

#WT 5xFAD vs. WT 5xFAD_TREM2 lipidomics with all time points
perma_df19 = merged_norm_data2 %>% filter(Group2 %in% c("WT 5xFAD", "WT 5xFAD_TREM2"))
perma_res19 = adonis2(select(perma_df19, !CageID:Background) ~ Group2 * Sex, 
                      data = perma_df19, method = "euclidian", strata = perma_df19$Age)

#5xFAD vs WT lipidomics with all time points
perma_df20 = merged_norm_data2 %>% filter(Group2 %in% c("HEMI 5xFAD_TREM2", "HEMI 5xFAD"))
perma_res20 = adonis2(select(perma_df20, !CageID:Background) ~ Group2 * Sex, 
                      data = perma_df20, method = "euclidian", strata = perma_df20$Age)

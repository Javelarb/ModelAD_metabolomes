library(tidyverse)
library(MetaboAnalystR)
library(readxl)
library(vegan)
library(pairwiseAdonis)
library(pheatmap)

setwd("/media/julio/Storage/mAD_metabolomes/new_lipidomics/")

z_score <- function(x){
  (x - mean(x)) / sd(x)
}

#Metadata
metadata = read.delim("metadata.txt")
metadata$Name = as.character(metadata$Name)
metadata$Age = as.character(metadata$Age)

#SampleIDs
SampleIDs = read_excel("u54-ad-sample-list.xlsx", 
                       col_types = c("text", "text", "skip", 
                                     "skip", "skip", "skip"))

#Importing and wrangling data
plasma_lip_df = as.data.frame(t(read.csv("u54_ad_batch_corrected_data_TL_raw.csv", row.names = 1, check.names = F))) %>%
  replace(is.na(.), 0) %>%
  rownames_to_column(var = "tm_id") %>%
  inner_join(SampleIDs) %>%
  select(!tm_id) %>%
  pivot_longer(!sample_id) %>%
  mutate(sample_id = paste0(sample_id, "_plasma"))

plasma_tar_df = as.data.frame(t(read.csv("u54_ad_batch_corrected_data_TM_raw.csv", row.names = 1, check.names = F))) %>%
  replace(is.na(.), 0) %>%
  rownames_to_column(var = "tm_id") %>%
  inner_join(SampleIDs) %>%
  select(!tm_id) %>%
  pivot_longer(!sample_id) %>%
  mutate(sample_id = paste0(sample_id, "_plasma"))

liver_df = read_excel("20220519_mapstone_ad_model_liver_107_5500QTRAP_data_combined.xlsx", col_names = F)
names(liver_df) = liver_df[2,]
names(liver_df)[2:3] = c("Common name", "Mode")
liver_df = liver_df[-(1:2),]

liver_lip_df = liver_df %>%
  filter(Mode == "Targeted lipidomics") %>%
  select(!c(`Common name`, Mode)) %>%
  column_to_rownames(var = "Sample IDs")

liver_lip_df = as.data.frame(t(liver_lip_df)) %>%
  replace(is.na(.), 0) %>%
  rownames_to_column(var = "sample_id") %>%
  pivot_longer(!sample_id) %>%
  mutate(sample_id = paste0(sample_id, "_liver"))
liver_lip_df$value = as.numeric(liver_lip_df$value) #Abundance was character

liver_tar_df = liver_df %>%
  filter(Mode == "Tarrgeted Metabolomics") %>%
  select(!c(`Common name`, Mode)) %>%
  column_to_rownames(var = "Sample IDs")

liver_tar_df = as.data.frame(t(liver_tar_df)) %>%
  replace(is.na(.), 0) %>%
  rownames_to_column(var = "sample_id") %>%
  pivot_longer(!sample_id) %>%
  mutate(sample_id = paste0(sample_id, "_liver"))
liver_tar_df$value = as.numeric(liver_tar_df$value) #Abundance was character

brain_df = read_excel("20220519_mapstone_ad_model_cortex_159_5500QTRAP_combined_data.xlsx", col_names = F)
names(brain_df) = brain_df[2,]
brain_df = brain_df[-(1:2),]

brain_lip_df = brain_df %>%
  select(!`9732`) %>% #This sample ID is duplicated for whatever reason
  filter(Mode == "Targeted lipidomics") %>%
  select(!c(`Sample IDs`, Mode)) %>%
  column_to_rownames(var = "Metabolites")

brain_lip_df = as.data.frame(t(brain_lip_df)) %>%
  replace(is.na(.), 0) %>%
  rownames_to_column(var = "sample_id") %>%
  pivot_longer(!sample_id) %>%
  mutate(sample_id = paste0(sample_id, "_brain"))

brain_tar_df = brain_df %>%
  select(!`9732`) %>% #This sample ID is duplicated for whatever reason
  filter(Mode == "Targeted Metabolomics") %>%
  select(!c(`Sample IDs`, Mode)) %>%
  column_to_rownames(var = "Metabolites")

brain_tar_df = as.data.frame(t(brain_tar_df)) %>%
  replace(is.na(.), 0) %>%
  rownames_to_column(var = "sample_id") %>%
  pivot_longer(!sample_id) %>%
  mutate(sample_id = paste0(sample_id, "_brain"))

#Combining dataframes
#Quantile normalization occurs across samples so you should have all the samples needed in your dataframe before normalizing.
#Keeping targeted metabolomics and lipidomics seperate.
lip_df_combined = bind_rows(plasma_lip_df, brain_lip_df, liver_lip_df) %>%
  pivot_wider(id_cols = "sample_id", names_from = "name", values_from = "value") %>%
  #select_if(~ !any(is.na(.))) %>% #Select only metabolites that do not have NA values (Present in all 3)
  separate_wider_delim(cols = sample_id, names = c("Name", "Mode"), delim = "_", cols_remove = F) %>%
  inner_join(metadata) %>%
  filter(!sample_id %in% c("9379_plasma", "7772_plasma", "20740_brain"), !Age == "18") %>%
  relocate(CageID:Background, .after = sample_id)

lip_df_combined[,-(1:12)] = log2(lip_df_combined[,-(1:12)]) %>% #Normalization
  replace(is.na(.), 0)
#write.table(x = lip_df_combined, file = "BLP_combined_lipidomics.txt", quote = F, sep = "\t", row.names = F, col.names = T)

tar_df_combined = bind_rows(plasma_tar_df, brain_tar_df, liver_tar_df) %>%
  pivot_wider(id_cols = "sample_id", names_from = "name", values_from = "value") %>%
  #select_if(~ !any(is.na(.))) %>%
  separate_wider_delim(cols = sample_id, names = c("Name", "Mode"), delim = "_", cols_remove = F) %>%
  inner_join(metadata) %>%
  filter(!Age == "18") %>%
  relocate(CageID:Background, .after = sample_id)

tar_df_combined[,-(1:12)] = log2(tar_df_combined[,-(1:12)]) %>% #Normalization
  replace(is.na(.), 0)
#write.table(x = tar_df_combined, file = "BLP_combined_targmet.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#PCA Plots
x = data.frame(Mode_var = c("brain", "brain", "liver", "liver", "plasma", "plasma"),
           Age_var = c("4", "12", "4", "12", "4", "12"))

for (i in 1:6) {
lipid_dm = lip_df_combined %>% 
  filter(Mode == x$Mode_var[i], Age == x$Age_var[i]) %>%
  column_to_rownames(var = "sample_id") %>%
  select(!Name:Background) %>%
  vegdist(., method = "euclidean")
lip_pcoa = cmdscale(lipid_dm, eig = T, k = nrow(filter(lip_df_combined, Mode == x$Mode_var[i], Age == x$Age_var[i]))-1, add = T)
lip_pcoa_eig = eigenvals(lip_pcoa)
lip_pcoa_var = lip_pcoa_eig/sum(lip_pcoa_eig)
lip_pcoa_var[1:3] #First 3 axes

lip_pcoa_plot = as.data.frame(lip_pcoa$points[,1:3]) %>%
  rownames_to_column(var = "Name") %>%
  separate_wider_delim(cols = Name, names = c("sample_id", "Mode"), delim = "_", cols_remove = F) %>%
  inner_join(metadata, by = join_by("sample_id" == "Name")) %>%
  unite("Combined", Group2, Mode, remove = F)

ggplot(data = lip_pcoa_plot) +
  aes(x = V1, y = V2) +
  stat_ellipse(geom = "polygon", aes(fill = Group2), alpha = 0.1) +
  geom_point(aes(pch = Sex, fill = Group2), size = 4, alpha = .5) +
  geom_text(aes(label = Sex), size = 2, alpha = .25) +
  scale_shape_manual(values = c(21, 22, 23)) +
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw() +
  scale_fill_manual(values = c("firebrick3", "steelblue1", "gold", "forestgreen")) +
  labs(x = bquote("PC1:"~.(round(lip_pcoa_var[1]*100, digits = 1))~"%"), 
       y = bquote("PC2:"~.(round(lip_pcoa_var[2]*100, digits = 1))~"%"),
       subtitle = bquote(.(x$Mode_var[i])~"-"~.(x$Age_var[i])~"mo."), fill = "Genotype", pch = "Sex", 
       title = "Lipidomics")

#Targeted metabolomics
tar_dm = tar_df_combined %>%
  filter(Mode == x$Mode_var[i], Age == x$Age_var[i]) %>%
  column_to_rownames(var = "sample_id") %>%
  select(!Name:Background) %>%
  vegdist(., method = "euclidean")
tar_pcoa = cmdscale(tar_dm, eig = T, k = nrow(filter(tar_df_combined, Mode == x$Mode_var[i], Age == x$Age_var[i]))-1, add = T)
tar_pcoa_eig = eigenvals(tar_pcoa)
tar_pcoa_var = tar_pcoa_eig/sum(tar_pcoa_eig)
tar_pcoa_var[1:3] #First 3 axes

tar_pcoa_plot = as.data.frame(tar_pcoa$points[,1:3]) %>%
  rownames_to_column(var = "Name") %>%
  separate_wider_delim(cols = Name, names = c("sample_id", "Mode"), delim = "_", cols_remove = F) %>%
  inner_join(metadata, by = join_by("sample_id" == "Name")) %>%
  unite("Combined", Group2,Mode, remove = F)

ggplot(data = tar_pcoa_plot) +
  aes(x = V1, y = V2) +
  stat_ellipse(geom = "polygon", alpha = .1, aes(fill = Group2)) +
  geom_point(aes(pch = Sex, fill = Group2), size = 4, alpha = .5) +
  geom_text(aes(label = Sex), size = 2, alpha = .25) +
  scale_shape_manual(values = c(21, 22, 23)) +
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw() +
  scale_fill_manual(values = c("firebrick3", "steelblue1", "gold", "forestgreen")) +
  labs(x = bquote("PC1:"~.(round(tar_pcoa_var[1]*100, digits = 1))~"%"), 
       y = bquote("PC2:"~.(round(tar_pcoa_var[2]*100, digits = 1))~"%"),
       subtitle = bquote(.(x$Mode_var[i])~"-"~.(x$Age_var[i])~"mo."), fill = "Genotype", pch = "Sex", 
       title = "Targeted metabolomics")

#PERMANOVA
lip_perma_df = lip_df_combined %>% 
  filter(Mode == x$Mode_var[i], Age == x$Age_var[i]) %>%
  column_to_rownames(var = "sample_id")

pairwise.adonis2(select(lip_perma_df, !Name:Background) ~ Group2*Sex, data = lip_perma_df, nperm = 999, method = "euclidean")

tar_perma_df = tar_df_combined %>% 
  filter(Mode == x$Mode_var[i], Age == x$Age_var[i]) %>%
  column_to_rownames(var = "sample_id")

pairwise.adonis2(select(tar_perma_df, !Name:Background) ~ Group2*Sex, data = tar_perma_df, nperm = 999, method = "euclidean")

#Heatmap
#Lipidomics
lip_zscore_df = apply(X = select(lip_perma_df, !Name:Background), MARGIN = 2, FUN = z_score)
lip_zscore_df = as.data.frame(lip_zscore_df)
lip_zscore_df = replace(lip_zscore_df, is.na(lip_zscore_df), 0)

lip_merged_zscores = lip_zscore_df %>% 
  select(which(!colSums(.) == 0)) %>%
  rownames_to_column(var = "sample_id") %>%
  separate_wider_delim(cols = sample_id, names = c("Name", "Mode"), delim = "_", cols_remove = F) %>%
  inner_join(metadata) %>% 
  relocate(CageID:Background, .after = sample_id) %>%
  column_to_rownames(var = "sample_id") %>%
  arrange(Group2, Sex)

pheatmap(t(select(lip_merged_zscores, !Name:Background)), 
         annotation_col = data.frame(rownames_to_column(lip_merged_zscores, var = "sample_id") %>% 
                                       select(sample_id, Sex) %>% column_to_rownames(var = "sample_id"),
                                     rownames_to_column(lip_merged_zscores, var = "sample_id") %>% 
                                       select(sample_id, Group2) %>% column_to_rownames(var = "sample_id")),
         show_colnames = F, fontsize_row = 1, cluster_cols = T, main = bquote("Lipidomics -"~.(x$Mode_var[i])~"-"~.(x$Age_var[i])~"mo."))

#Targeted metabolomics
tar_zscore_df = apply(X = select(tar_perma_df, !Name:Background), MARGIN = 2, FUN = z_score)
tar_zscore_df = as.data.frame(tar_zscore_df)
tar_zscore_df = replace(tar_zscore_df, is.na(tar_zscore_df), 0)

tar_merged_zscores = tar_zscore_df %>% 
  select(which(!colSums(.) == 0)) %>%
  rownames_to_column(var = "sample_id") %>%
  separate_wider_delim(cols = sample_id, names = c("Name", "Mode"), delim = "_", cols_remove = F) %>%
  inner_join(metadata) %>% 
  relocate(CageID:Background, .after = sample_id) %>%
  column_to_rownames(var = "sample_id") %>%
  arrange(Group2, Sex)

pheatmap(t(select(tar_merged_zscores, !Name:Background)), 
         annotation_col = data.frame(rownames_to_column(tar_merged_zscores, var = "sample_id") %>% 
                                       select(sample_id, Sex) %>% column_to_rownames(var = "sample_id"),
                                     rownames_to_column(tar_merged_zscores, var = "sample_id") %>% 
                                       select(sample_id, Group2) %>% column_to_rownames(var = "sample_id")),
         show_colnames = F, fontsize_row = 1, cluster_cols = T, main = bquote("Targ. metabolomics -"~.(x$Mode_var[i])~"-"~.(x$Age_var[i])~"mo."))
}

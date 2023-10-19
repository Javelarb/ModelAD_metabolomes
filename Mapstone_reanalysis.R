library(tidyverse)
library(MetaboAnalystR)
library(readxl)
library(vegan)
library(pheatmap)

setwd("/media/julio/Storage/mAD_metabolomes/new_lipidomics/")

#SampleIDs which translates from what MS core used
SampleIDs = read_excel("u54-ad-sample-list.xlsx", 
                       col_types = c("text", "text", "skip", 
                                     "skip", "skip", "skip"))

#Data
b1_b2_lipid_df = as.data.frame(t(read.csv("/media/julio/Storage/mAD_metabolomes/new_lipidomics/u54_ad_batch_corrected_data_TL.csv", row.names = 1, check.names = F))) %>%
  rownames_to_column(var = "name")

#Merge with sampleIDs
b1_b2_lipid_df = b1_b2_lipid_df %>% 
  left_join(SampleIDs, by = join_by("name" == "tm_id")) #%>%
  #filter(complete.cases(.))

#Replace NAs with zero
b1_b2_lipid_df2 = b1_b2_lipid_df %>% 
  replace(is.na(.), 0) %>%
  select(!name) %>%
  column_to_rownames(var = "sample_id")

#Importing and formatting metadata for batches 1 & 2
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

#Merge data with metadata
combined_merged_data = rownames_to_column(b1_b2_lipid_df2, var = "Name") %>%
  #rownames_to_column(var = "Name") %>%
  inner_join(metadata)

#Write and export for metaboanalyst
#write.table(x = metadata, file = "metaboanalyst_metadata.txt", quote = F, sep = "\t", row.names = F, col.names = T)
#write.table(x = test2, file = "metaboanalyst_data.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#Metaboanalyst section (This is old code for when the batch effects where not fixed yet)
#mSet = InitDataObjects(data.type = "conc", anal.type = "stat", paired = F)
#mSet = Read.TextData(mSet, "metaboanalyst_data.txt", format = "rowu", lbl.type = "disc")
#mSet = ReplaceMin(mSet)
#mSet = PreparePrenormData(mSet)
#mSet = Normalization(mSet, "QuantileNorm", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)

#normalized_data = mSet[["dataSet"]][["norm"]] #%>% 
  #rownames_to_column(var = "Name") %>%
  #inner_join(metadata)

#Continuation
normalized_data = combined_merged_data %>%
  #filter(Group2 %in% c("HEMI 5xFAD_TREM2", "HEMI 5xFAD"), Age %in% c("4")) %>%
  column_to_rownames(var = "Name") %>%
  select(!CageID:Background) 

#PCoA
dist_matrix = vegdist(normalized_data, method = "euclidean")
pcoa = cmdscale(dist_matrix, eig = T, k = nrow(normalized_data)-1, add = T)
pcoa_eig = eigenvals(pcoa)
pcoa_var = pcoa_eig/sum(pcoa_eig)
pcoa_var[1:3] #First 3 axes

pcoa_plot_df = as.data.frame(pcoa$points[,1:2]) %>%
  rownames_to_column(var = "Name") %>%
  inner_join(metadata)

#Batch effect fixed
ggplot(data = pcoa_plot_df) +
  aes(x = V1, y = V2, color = Batch) +
  geom_text(aes(label = Name)) +
  theme_bw() +
  labs(x = bquote("PC1:"~.(round(pcoa_var[1]*100, digits = 1))~"%"), 
       y = bquote("PC2:"~.(round(pcoa_var[2]*100, digits = 1))~"%"),
       title = "New normalized lipidomics")

#Permanovas
merged_norm_data = inner_join(metadata, rownames_to_column(normalized_data, var = "Name")) %>%
  column_to_rownames(var = "Name") #%>%
  #filter(complete.cases(.))

merged_norm_data$CageID = as.factor(merged_norm_data$CageID)


#### Old code ####
#Single factor PERMA, Batch effects ~3% of the variance. 
#adonis2(select(merged_norm_data, !CageID:Background) ~ Batch, data = merged_norm_data, method = "euclidian")          

#What about the other factors?
#adonis2(select(merged_norm_data, !CageID:Background) ~ Age + Sex + AD_path + Background + Batch + CageID, 
#        data = merged_norm_data, method = "euclidian") #Maybe AD path will become sig with specific time points

#Differential abundance
#t_test_df = merged_norm_data %>%
#  select(!CageID:Background) %>%
#  map_df(~ broom::tidy(t.test(. ~ merged_norm_data$Group2)), .id = 'var') %>%
#  select(var, p.value)

#t_test_df$padj = p.adjust(t_test_df$p.value)

#Need to calc log2 fold change
#log2fc = merged_norm_data %>%
#  group_by(Group2) %>%
#  summarise(across(!CageID:Background, mean)) %>%
#  column_to_rownames(var = "Group2") %>%
#  t(.)

#log2fc = as.data.frame(log2fc)
#log2fc$Ratio = log2fc$`HEMI 5xFAD_TREM2`/log2fc$`HEMI 5xFAD`
#log2fc$Log2FC = log(log2fc$Ratio, base = 2)

#log2fc = log2fc %>%
#  rownames_to_column(var = "var") %>%
#  left_join(t_test_df)

#log2fc = log2fc %>% 
#  mutate(color = case_when(
#  padj < 0.05 & Log2FC > 1 ~ "color1",
#  padj < 0.05 & Log2FC < -1 ~ "color2",
#  padj > 0.05 ~ "color3",
#  Log2FC > -1 ~ "color3",
#  Log2FC < 1 ~ "color3"
#))

#how many are infinity?
#sum(log2fc$Log2FC == Inf) #360
#sum(log2fc$Log2FC == -Inf) #126

#ggplot(data = log2fc) +
#  aes(x = Log2FC, y = -log10(padj), color = color) +
#  geom_hline(yintercept = -log10(0.05), lty = 2) +
#  geom_vline(xintercept = 1, lty = 3) +
#  geom_vline(xintercept = -1, lty = 3) +
#  geom_point(alpha = 0.5, size = 2) + 
#  scale_color_manual(values = c("steelblue", "#BF0D3E", "gray"), labels = c("Increased", "Decreased", "N.S.")) +
#  labs(title = "5xFAD*TREM2 vs. 5xFAD", subtitle = "4 mo. lipidomics",
#       y = expression("-log"[10]*"(p-adj)"), color = "") +
#  theme_bw() +
#  ggrepel::geom_text_repel(aes(x = Log2FC, y = -log10(padj), label = var), color = "black", size = 2, max.overlaps = 8) +
#  annotate("text", x = -Inf, y = Inf, size = 2.5, hjust = 0, vjust = 1, label = bquote("Negative:"~.(sum(log2fc$padj < 0.05 & log2fc$Log2FC < -1)))) +
#  annotate("text", x = Inf, y = Inf, size = 2.5, hjust = 1, vjust = 1, label = bquote("Positive:"~.(sum(log2fc$padj < 0.05 & log2fc$Log2FC > 1))))

#Heatmap
#z_score <- function(x){
#  (x - mean(x)) / sd(x)
#}

#zscore_df = apply(X = select(merged_norm_data, !CageID:Background), MARGIN = 2, FUN = z_score)
#zscore_df = as.data.frame(zscore_df)

#Top N
#top_100 = slice_min(log2fc, order_by = padj, n = 100)
#zscore_df = select(zscore_df, top_100$var)

#merged_zscore_data = inner_join(metadata, rownames_to_column(zscore_df, var = "Name")) %>%
#  column_to_rownames(var = "Name") %>%
#  filter(complete.cases(.)) %>%
#  arrange(desc(Age), Group2, Sex)

#pheatmap(t(select(merged_zscore_data, !(CageID:Background))), 
#         annotation_col = data.frame(merged_zscore_data %>% rownames_to_column(var = "Name") %>% select(Name, Sex) %>% column_to_rownames(var = "Name"),
#                                     merged_zscore_data %>% rownames_to_column(var = "Name") %>% select(Name, Group2) %>% column_to_rownames(var = "Name"),
#                                     merged_zscore_data %>% rownames_to_column(var = "Name") %>% select(Name, Age) %>% column_to_rownames(var = "Name")),
#         show_colnames = F, main = "Z-scores - Top 100 features", cluster_cols = T, fontsize_row = 5)
#### Old code ####


#Specific AD path groups
perma_df1 = merged_norm_data %>% filter(Age == 18, Group2 %in% c("WT 5xFAD", "HEMI 5xFAD"))
#perma_df1 = perma_df1[!rownames(perma_df1) == "8123",]

adonis2(select(perma_df1, !CageID:Background) ~ Group2 + Sex + CageID, 
        data = perma_df1, method = "euclidian")

dist_matrix = perma_df1 %>%
  select(!CageID:Background) %>%
  vegdist(., method = "euclidean")

pcoa = cmdscale(dist_matrix, eig = T, k = nrow(perma_df1)-1, add = T)
pcoa_eig = eigenvals(pcoa)
pcoa_var = pcoa_eig/sum(pcoa_eig)
pcoa_var[1:3] #First 3 axes

pcoa_plot_df = as.data.frame(pcoa$points[,1:2]) %>%
  rownames_to_column(var = "Name") %>%
  inner_join(metadata)

loadings = perma_df1 %>%
  select(!CageID:Background) %>%
  envfit(ord = pcoa, env = ., na.rm = TRUE)

#make the loadings
loadings_df = data.frame(scores(x = loadings, display = "vectors"), R = loadings$vectors$r, pval = loadings$vectors$pvals) %>%
  filter(pval < 0.05) %>%
  rownames_to_column(var = "Metabolite")

ggplot(data = pcoa_plot_df) +
  aes(x = V1, y = V2, fill = Group2, pch = Batch, label = Sex) +
  geom_segment(data = loadings_df, inherit.aes = F, aes(x = 0, xend = Dim1*max(abs(pcoa_plot_df$V1)), y = 0, yend = Dim2*max(abs(pcoa_plot_df$V2))), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey90") +
  #stat_ellipse(aes(group = Sex, lty = Sex), alpha = 1/2, show.legend = F) +
  geom_point(size = 5) +
  #ggrepel::geom_text_repel(data = loadings_df, inherit.aes = FALSE, aes(x = Dim1*max(abs(pcoa_plot_df$V1)), 
  #                                                                      y = Dim2*max(abs(pcoa_plot_df$V2)), label = Metabolite), size = 4, max.overlaps = 20) + 
  geom_text(size = 3) +
  theme_bw() +
  guides(color = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  scale_shape_manual(values = c(21, 22)) +
  labs(x = bquote("PC1:"~.(round(pcoa_var[1]*100, digits = 1))~"%"), 
       y = bquote("PC2:"~.(round(pcoa_var[2]*100, digits = 1))~"%"),
       title = "WT vs. WT TREM2 - 12 mo.")

#Heatmap
zscore_df = apply(X = select(perma_df1, `CE(12:0)`:`Tridecanoylcarnitine_AC(13:0)`), MARGIN = 2, FUN = z_score)
zscore_df = as.data.frame(zscore_df)

merged_zscore_data = inner_join(metadata, rownames_to_column(zscore_df, var = "Name")) %>%
  column_to_rownames(var = "Name") %>%
  select_if(~ !any(is.na(.)))

merged_zscore_data = arrange(merged_zscore_data, Group2, Sex)

pheatmap(t(select(merged_zscore_data, `CE(12:0)`:`Tridecanoylcarnitine_AC(13:0)`)), 
         annotation_col = data.frame(merged_zscore_data %>% rownames_to_column(var = "Name") %>% select(Name, Sex) %>% column_to_rownames(var = "Name"),
                                     merged_zscore_data %>% rownames_to_column(var = "Name") %>% select(Name, Group2) %>% column_to_rownames(var = "Name")),
         show_colnames = F, fontsize_row = 1, main = "WT vs. 5xFAD Z-scores - 18 mo.", cluster_cols = T)

#Note: its named t test but its using a wilcox test
t_test_df = perma_df1 %>%
  select(`CE(12:0)`:`Tridecanoylcarnitine_AC(13:0)`) %>%
  map_df(~ broom::tidy(wilcox.test(. ~ perma_df1$Group2)), .id = 'var') %>%
  select(var, p.value)

t_test_df$padj = p.adjust(t_test_df$p.value)

#Need to calc log2 fold change
log2fc = perma_df1 %>%
  group_by(Group2) %>%
  summarise(across(`CE(12:0)`:`Tridecanoylcarnitine_AC(13:0)`, mean)) %>%
  column_to_rownames(var = "Group2") %>%
  t(.)

log2fc = as.data.frame(log2fc)
log2fc$Ratio = log2fc$`HEMI 5xFAD`/log2fc$`WT 5xFAD`
log2fc$Log2FC = log(log2fc$Ratio, base = 2)

log2fc = log2fc %>%
  rownames_to_column(var = "var") %>%
  left_join(t_test_df)

log2fc = log2fc %>% 
  mutate(color = case_when(
    padj < 0.05 & Log2FC > 1 ~ "color1",
    padj < 0.05 & Log2FC < -1 ~ "color2",
    padj > 0.05 ~ "color3",
    Log2FC > -1 ~ "color3",
    Log2FC < 1 ~ "color3"
  ))

#how many are infinity?
sum(log2fc$Log2FC == Inf) #360
sum(log2fc$Log2FC == -Inf) #126

ggplot(data = log2fc) +
  aes(x = Log2FC, y = -log10(padj), color = color) +
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  geom_vline(xintercept = 1, lty = 3) +
  geom_vline(xintercept = -1, lty = 3) +
  geom_point(alpha = 0.5, size = 2) + 
  scale_color_manual(values = c("steelblue", "#BF0D3E", "gray"), labels = c("Increased", "Decreased", "N.S.")) +
  labs(title = "Significantly different metabolites", subtitle = "Unnormalized data",
       y = expression("-log"[10]*"(p-adj)"), color = "Abundance in\nBatch 1") +
  theme_bw() +
  #ggrepel::geom_text_repel(aes(x = Log2FC, y = -log10(padj), label = var), color = "black", size = 2, max.overlaps = 8) +
  annotate("text", x = -Inf, y = Inf, size = 2.5, hjust = 0, vjust = 1, label = bquote("Negative:"~.(sum(log2fc$padj < 0.05 & log2fc$Log2FC < -1)))) +
  annotate("text", x = Inf, y = Inf, size = 2.5, hjust = 1, vjust = 1, label = bquote("Positive:"~.(sum(log2fc$padj < 0.05 & log2fc$Log2FC > 1))))

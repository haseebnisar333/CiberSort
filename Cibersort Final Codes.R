sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
results <- cibersort(sig_matrix, mixture_file)

#my dataset mixture file
mixture_file2 <- "/Users/Admin/Desktop/data/df1.txt"
results2 <- cibersort(sig_matrix, mixture_file2)

rownames(LM22) <- LM22$Gene.symbol
LM22 <- LM22[,-1]
# Assuming results2 is your data with immune cell types as columns
# Remove non-numeric columns (like 'P-value', 'Correlation', and 'RMSE') for correlation calculation
immune_cells_data <- results2[, 1:22]


# Calculate the Spearman correlation between immune cells and signature genes (you need to have signature genes data)
cor_matrix <- cor(immune_cells_data, LM22, method = "spearman")

#read gene expression matrix 
input <- df1
cibersort_perm = 100
#Quantile normalization of input mixture, default = FALSE for RNA-Seq data
cibersort_qn = TRUE
#whether to apply absolute mode in cibersort
cibersort_abs = TRUE

#sig.score = for each mixture sample, define S as the median expression,level of all genes in the signature matrix divided by the median expression level of all genes in the mixture. Multiple cell subset fractions by S.
cibersort_abs_method = "sig.score"
res_ciber <- cibersort(sig_matrix, mixture_file2, perm = cibersort_perm, QN = cibersort_qn)
head(res_ciber)
res_ciber <- res_ciber[,1:22]
res_ciber


M = cor(immune_cells_data)
corrplot(M, order = 'hclust', addrect = 2)

results3 <- cbind(sampletype = rownames(immune_cells_data), immune_cells_data)

res_ciber_df <- data.frame(results3[,1:23]) %>% 
  pivot_longer(cols = 2:23, names_to = "Celltype", values_to = "Proportion") %>% 
  mutate(Proportion = as.numeric(Proportion)) 

## plot proportion estimates in a heat map, rows indexed by samples and
## columns indexed by cell types
## color filled with estimated proportion of tissue of composition per cell

ggplot(res_ciber_df, 
       aes(x = Celltype, y = sampletype, fill = Proportion)) + geom_tile() + 
  scale_fill_viridis_c(breaks = seq(0, 1, 0.1), 
                       limits = c(0, 1), 
                       label = seq(0, 100, 10)) + 
  labs(title = "Reference Based Deconvolution with CIBERSORT",
       x = "Cell-Type",
       y = "Sample",
       fill = "% of Tissue") + 
  theme(text = element_text(family = "Fira Sans"),
        plot.background = element_rect(fill = "white"), 
        axis.text.x = element_text(angle = 45, vjust = 0.65), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key.size = unit(3,"line"))


# Ensure proportions sum to 1 per sample
res_ciber_df2 <- res_ciber_df %>%
  group_by(sampletype) %>%
  mutate(Proportion = Proportion / sum(Proportion)) %>%
  ungroup()
# Stacked barplot: Immune cell distribution per sample
ggplot(res_ciber_df, aes(x = sampletype, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity") +
  labs(title = "Distribution of Immune Cells in Each Sample",
       x = "Sample", y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(family = "Fira Sans")) +
  scale_fill_viridis_d()


# Grouped barplot: Comparison between groups
group_means <- res_ciber_df %>%
  group_by(sampletype, Celltype) %>%
  summarise(Mean_Proportion = mean(Proportion), .groups = "drop")

ggplot(group_means, aes(x = Celltype, y = Mean_Proportion, fill = sampletype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of Immune Cell Expression Between Groups",
       x = "Cell Type", y = "Mean Proportion", fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(family = "Fira Sans")) +
  scale_fill_manual(values = c("Control" = "#1f78b4", "Patient" = "#e31a1c"))


immune_cells_data <- as.data.frame(immune_cells_data)
immune_cells_data2 <- cbind(Sample = rownames(immune_cells_data), immune_cells_data)
# Replace the labels in the "Sample" column
immune_cells_data2 <- immune_cells_data2 %>%
  mutate(Sample = gsub("^Control.*", "Control", Sample),  # Replace anything starting with "Control"
         Sample = gsub("^Patient.*", "Patient", Sample))  # Replace anything starting with "Patient"
control_data <- immune_cells_data2 %>% filter(Sample == "Control")
patient_data <- immune_cells_data2 %>% filter(Sample == "Patient")

# Step 3: Compute correlations for each cell type
# Remove the "Sample" column for correlation computation
celltype_cols <- colnames(immune_cells_data2)[2:23] # Columns for cell types

# Aggregate data by group
group_means <- immune_cells_data2 %>%
  group_by(Sample) %>% 
  summarise(across(everything(), mean, na.rm = TRUE)) 

# Transpose the data for easier correlation computation
transposed_data <- as.data.frame(t(group_means[-1]))  # Remove Sample column and transpose
# Set proper column names (Control, Patient) for clarity
colnames(transposed_data) <- c("Control", "Patient")
transposed_data2 <- cbind(ENSEMBLID = rownames(transposed_data), transposed_data)
colnames(transposed_data2) <- gsub("ENSEMBLID", "Celltype", colnames(transposed_data2))
write_xlsx(immune_cells_data2, "immune_cells_data2.xlsx")
write_xlsx(transposed_data2, "transposed_data2.xlsx")



# Step 1: Reshape the data into long format
long_data <- transposed_data2 %>%
  pivot_longer(cols = c("Control", "Patient"), 
               names_to = "Group", 
               values_to = "Proportion")


# Step 2: Create the bar plot
ggplot(long_data, aes(x = Celltype, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 'dodge' places bars side by side
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "Proportion of Each Cell Type in Control vs Patient",
    x = "Cell Type",
    y = "Proportion",
    fill = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for readability
    plot.title = element_text(hjust = 0.5)  # Center plot title
  )
data <- transposed_data2








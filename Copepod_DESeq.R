library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(scales)
library(pheatmap)

options(bitmapType = "cairo")

setwd("/gpfs1/cl/pbio3990/GroupProjects/Acartia/")

################################################
#
#normalize hudsonica data 
#
###############################################
ahuddata <- read.delim("ahud_samples_R.txt")
ahuddata_subset <- subset(ahuddata,treatment=="AM" & generation=="F0")
ahud_subset <- subset(ahuddata, subset = (treatment %in% c("AM","OW")& generation == "F0"))


countsTable <- read.table("/gpfs1/cl/pbio3990/GroupProjects/Acartia/salmon.isoform.counts.matrix.filteredAssembly",
                          header = TRUE, row.names = 1)

##########################
#
#subset the tonsa data files
#
#########################

tonsadata <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                        header = TRUE,
                        stringsAsFactors = TRUE,
                        row.names = 1)
head(tonsadata)
treatment <- tonsadata[which(tonsadata$FinalTemp =="BASE"), ]
print(treatment)

#keeping these colums because they are the correct conditions for this study 
columns_to_keep <- grep("^AA_F0|HA_F0", colnames(countsTable), value=TRUE)
subset.counts_hudsonica <- countsTable[, columns_to_keep]

#round the subsetted data becase DeSeq does not like decimals
subset.counts_hudsonica <- round(subset.counts_hudsonica)

head(subset.counts_hudsonica)


################################################
#
#subsetting the tonsa counts matrix
#
###############################################
countsTable2 <- read.table("/gpfs1/cl/pbio3990/GroupProjects/Acartia/tonsa_counts.txt",
                           header = TRUE, row.names = 1)
subset.counts_tonsa <- countsTable2[, c("N1C3","N2C2","N1C4","N1C2","N2C3","N2C1")]


#round the subsetted data becase DeSeq does not like decimals
subset.counts_tonsaround<- round(subset.counts_tonsa)

head(subset.counts_tonsaround)


##################################################
#
# filtering BLAST query genes for best evalue 
#
##################################################

library(dplyr)

blast_data <- read.table("acartia_blast_max1", header = FALSE, sep = "\t", 
                         stringsAsFactors = FALSE)

colnames(blast_data) <- c("query", "subject", "identity", "length", "mismatch", 
                          "gapopen", "q.start", "q.end", "s.start", "s.end", 
                          "evalue", "bitscore")
best_blast_matches <- blast_data %>%
  group_by(query)%>%
  filter(evalue == min(evalue))%>%
  ungroup()

best_blast_matches <- best_blast_matches %>%
  group_by(query) %>%
  slice(1) %>%
  ungroup()

print(best_blast_matches)


############################################################
#
#call genes from counts matrix and BLAST 
#merge the filtered BLAST file with the counts matricies 
#
###########################################################

merged_ahud1 <- merge(best_blast_matches, subset.counts_hudsonica, by.x = "query", by.y = "row.names", all = TRUE)

merged_tonsa <- merge(best_blast_matches, subset.counts_tonsaround, by.x = "subject", by.y = "row.names", all = TRUE)

merged_all <- merge (merged_ahud1, merged_tonsa, by.x = "query", by.y = "query")

############################################################
#
# subset merged table
#
#############################################################

subset.merged_all <- merged_all[, c("subject.x","query","AA_F0_Rep1_","AA_F0_Rep2_","AA_F0_Rep3_","HA_F0_Rep1_","HA_F0_Rep2_","HA_F0_Rep3_","N2C2","N1C4","N1C2","N2C3","N2C1","N1C3")]

subset.merged_all$combined <- paste(subset.merged_all[,1], subset.merged_all[,2], sep = "_")

write.csv(subset.merged_all, "/gpfs1/cl/pbio3990/GroupProjects/Acartia/subset.merged_all", row.names = FALSE)

#make column 1 the row names 
rownames(subset.merged_all) <- subset.merged_all$combined
#merge the names of the genes so they are all unique and then make them the row names 

subset.merged_all <- subset.merged_all[, -c(1,2, ncol(subset.merged_all))]

################################
#
#Normalize subsetted merged file
#
##############################
head(subset.merged_all)
dim(subset.merged_all)
# 5669     12
# subset.merged_all[acartia_counts < 0] <- 0
# any(acartia_counts < 0)
# which(acartia_counts < 0, arr.ind = TRUE)

acartia_conds <- read.delim("~/Projects/eco_genomics_2/final_project/scripts/acartia_conditions.txt")
head(acartia_conds)

#create dds object 
dds <- DESeqDataSetFromMatrix(countData = subset.merged_all, colData = acartia_conds,
                              design = ~ Treatment + Species + Treatment:Species)
dds <- DESeq(dds)

#final normalized counts matrix 
normalized_counts <- counts(dds, normalized = TRUE)
file.path <- "/gpfs1/cl/pbio3990/GroupProjects/Acartia/normalized_counts"
write.csv(normalized_counts, file = file.path, row.names = TRUE)

head(normalized_counts)
dim(normalized_counts)
#[1] 5669   12

#DESeq results, use these names to run contrasts 
resultsNames(dds)
#[1] "Intercept"                 "Treatment_OW_vs_AM"        "Species_Atonsa_vs_Ahud"    "TreatmentOW.SpeciesAtonsa"

###################################################
#
#Constrasts and Results!!
#
###################################################
#
#Contrast 1: species
#
##################################################

res_species <- results(dds, name = "Species_Atonsa_vs_Ahud", alpha = 0.05)
res_species <- res_species[order(res_species$padj),] #order by significance 
head(res_species)
summary(res_species)
# out of 5655 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 2601, 46%
#LFC < 0 (down)     : 1901, 34%
#outliers [1]       : 70, 1.2%
#low counts [2]     : 0, 0%
#(mean count < 0

#look at counts of a specific top gene that we're interested in to validate that the model is working 
d <- plotCounts(dds, gene = "TRINITY_DN134966_c2_g1_i1_TRINITY_DN369_c1_g1::TRINITY_DN369_c1_g1_i7::g.2127::m.2127", 
                int=(c("Treatment", "Species")), returnData = TRUE)
d
#count number in each of our samples in a specific gene, plot these data 

p <- ggplot(d, aes(x=Treatment, y= count, color=Treatment, shape=Species))+
  theme_minimal()+theme(text=element_text(size=20), panel.grid.major=element_line(color="grey"))
p <- p+geom_point(position=position_jitter(w=0.2,h=0), size=3)
p

#
#PCA top 5 differentially expressed genes
#

# Step 1: Extract the top 5 differentially regulated genes from your DESeq results
res_species <- results(dds, name = "Species_Atonsa_vs_Ahud", alpha = 0.05)

# Order by padj (adjusted p-value) and select top 5 most significant genes
res_species <- res_species[order(res_species$padj), ]

# Extract the top 5 genes with the smallest adjusted p-value (padj) and a log2FoldChange > 1 (optional)
top_genes <- head(res_species, 5)

# Step 2: Extract expression data for the top 5 genes
# Extract normalized counts for the top genes from your DESeqDataSet
top_gene_expression <- counts(dds, normalized = TRUE)[rownames(top_genes), ]

# Step 3: Run PCA on the expression data for these top genes (transpose to get samples as rows)
pca_result <- prcomp(t(top_gene_expression), scale. = TRUE)

# Step 4: Create a data frame with PCA results
pca_df <- data.frame(pca_result$x)

# Add sample information (assuming Treatment and Species are in the colData of `dds`)
pca_df$Sample <- colnames(top_gene_expression)
pca_df$Treatment <- colData(dds)$Treatment
pca_df$Species <- colData(dds)$Species

# Step 5: Plot the PCA results using ggplot2
library(ggplot2)

# Choose colors and shapes for the plot (modify according to your conditions)
Treatment_colors <- c("OW" = "tomato", "AM" = "cornflowerblue")
shapes_choose <- c("Ahud" = 16, "Atonsa" = 18)

# Create the PCA plot
PCA_species <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Species)) + 
  geom_point(size = 5) + 
  scale_shape_manual(values = shapes_choose) + 
  scale_color_manual(values = Treatment_colors) + 
  labs(title = "PCA of Top 5 Differentially Regulated Genes (Species: Tonsa vs Hudsonica)",
       x = paste0('PC1: ', round(100 * summary(pca_result)$importance[2, 1], 1), '%'),
       y = paste0('PC2: ', round(100 * summary(pca_result)$importance[2, 2], 1), '%')) + 
  theme_bw(base_size = 16)

# Step 6: Display the plot
print(PCA_species)

#save the plot
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/species_pca_topgenes", device = "pdf")

###########################################
#
# heatmap of top 5 genes 
#
###########################################

vsd <- vst(dds, blind=FALSE)

topgenes <- head(rownames(res_species), 100)
mat <- assay(vsd)[topgenes,]
mat <- mat-rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("Treatment", "Species")])
heatmap_species <- pheatmap(mat,annotation_col = df, show_rownames = FALSE, cluster_cols = F, cluster_rows = T)

ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/species_heatmap_topgenes", device = "pdf")

########################################
#
#second contrast, treatment
#
########################################
res_treatment <- results(dds, name = "Treatment_OW_vs_AM", alpha = 0.05)
res_treatment <- res_treatment[order(res_treatment$padj),]
head(res_treatment)
summary(res_treatment)
# out of 5655 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 454, 8%
# LFC < 0 (down)     : 362, 6.4%
# outliers [1]       : 70, 1.2%
# low counts [2]     : 0, 0%
# (mean count < 0)

#look at counts of a specific top gene that we're interested in to validate that the model is working 
d <- plotCounts(dds, gene = "TRINITY_DN137257_c0_g2_i1_TRINITY_DN502_c0_g1::TRINITY_DN502_c0_g1_i1::g.3527::m.3527", 
                int=(c("Treatment", "Species")), returnData = TRUE)
d
#count number in each of our samples in a specific gene, plot these data 

p <- ggplot(d, aes(x=Treatment, y=count, color=Treatment, shape=Species))+
  theme_minimal()+theme(text=element_text(size=20), panel.grid.major=element_line(color="grey"))
p <- p+geom_point(position=position_jitter(w=0.2,h=0), size=3)
p

#
#PCA of top differentially regulated genes 
#

# Step 1: Extract the top 5 differentially regulated genes from your DESeq results
res_treatment <- results(dds, name = "Treatment_OW_vs_AM", alpha = 0.05)

# Order by padj (adjusted p-value) and select top 5 most significant genes
res_treatment <- res_treatment[order(res_treatment$padj), ]

# Extract the top 5 genes with the smallest adjusted p-value (padj) and a log2FoldChange > 1 (optional)
top_genes <- head(res_treatment, 5)

# Step 2: Extract expression data for the top 5 genes
# Extract normalized counts for the top genes from your DESeqDataSet
top_gene_expression <- counts(dds, normalized = TRUE)[rownames(top_genes), ]

# Step 3: Run PCA on the expression data for these top genes (transpose to get samples as rows)
pca_result <- prcomp(t(top_gene_expression), scale. = TRUE)

# Step 4: Create a data frame with PCA results
pca_df <- data.frame(pca_result$x)

# Add sample information (assuming Treatment and Species are in the colData of `dds`)
pca_df$Sample <- colnames(top_gene_expression)
pca_df$Treatment <- colData(dds)$Treatment
pca_df$Species <- colData(dds)$Species

# Step 5: Plot the PCA results using ggplot2
library(ggplot2)

# Choose colors and shapes for the plot (modify according to your conditions)
Treatment_colors <- c("OW" = "tomato", "AM" = "cornflowerblue")
shapes_choose <- c("Ahud" = 16, "Atonsa" = 18)

# Create the PCA plot
PCA_treatment <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Species)) + 
  geom_point(size = 5) + 
  scale_shape_manual(values = shapes_choose) + 
  scale_color_manual(values = Treatment_colors) + 
  labs(title = "PCA of Top 5 Differentially Regulated Genes (Treatment_OW_vs_AM)",
       x = paste0('PC1: ', round(100 * summary(pca_result)$importance[2, 1], 1), '%'),
       y = paste0('PC2: ', round(100 * summary(pca_result)$importance[2, 2], 1), '%')) + 
  theme_bw(base_size = 16)

# Step 6: Display the plot
print(PCA_treatment)
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/treatment_pca_topgenes", device = "pdf")

##################################################
#
# heatmap of top 5 differentially expressed genes
#
##################################################

topgenes <- head(rownames(res_treatment), 100)
mat <- assay(vsd)[topgenes,]
mat <- mat-rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("Treatment", "Species")])
heatmap_treatment <- pheatmap(mat,annotation_col = df, show_rownames = FALSE, cluster_cols = F, cluster_rows = T)
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/treatment_heatmap_topgenes", device = "pdf")

#############################################################
#
# third contrast, interaction of species and treatment
#
#############################################################

res_interaction <- results(dds, name = "TreatmentOW.SpeciesAtonsa", alpha = 0.05)
res_interaction <- res_interaction[order(res_interaction$padj),]
head(res_interaction)
summary(res_interaction)
#out of 5655 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 287, 5.1%
# LFC < 0 (down)     : 411, 7.3%
# outliers [1]       : 70, 1.2%
# low counts [2]     : 0, 0%
# (mean count < 0)

#look at counts of a specific top gene that we're interested in to validate that the model is working 
d <- plotCounts(dds, gene = "TRINITY_DN137257_c0_g2_i1_TRINITY_DN502_c0_g1::TRINITY_DN502_c0_g1_i1::g.3527::m.3527", 
                int=(c("Treatment", "Species")), returnData = TRUE)
d
#count number in each of our samples in a specific gene, plot these data 

p <- ggplot(d, aes(x=Treatment, y=count, color=Treatment, shape=Species))+
  theme_minimal()+theme(text=element_text(size=20), panel.grid.major=element_line(color="grey"))
p <- p+geom_point(position=position_jitter(w=0.2,h=0), size=3)
p

################################################
#
#heatmap of top 5 differentially expressed genes
#
###############################################

library(pheatmap)
vsd <- vst(dds, blind=FALSE)

topgenes <- head(rownames(res_interaction), 100)
mat <- assay(vsd)[topgenes,]
mat <- mat-rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("Treatment", "Species")])
pheatmap(mat,annotation_col = df, show_rownames = FALSE, cluster_cols = F, cluster_rows = T)
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/interaction_heatmap_topgenes", device = "pdf")

#PCA of top differentially regulated genes 

# Step 1: Extract the top 5 differentially regulated genes from your DESeq results
res_interaction <- results(dds, name = "TreatmentOW.SpeciesAtonsa", alpha = 0.05)

# Order by padj (adjusted p-value) and select top 5 most significant genes
res_interaction <- res_interaction[order(res_interaction$padj), ]

# Extract the top 5 genes with the smallest adjusted p-value (padj) and a log2FoldChange > 1 (optional)
top_genes <- head(res_interaction, 5)

# Step 2: Extract expression data for the top 5 genes
# Extract normalized counts for the top genes from your DESeqDataSet
top_gene_expression <- counts(dds, normalized = TRUE)[rownames(top_genes), ]

# Step 3: Run PCA on the expression data for these top genes (transpose to get samples as rows)
pca_result <- prcomp(t(top_gene_expression), scale. = TRUE)

# Step 4: Create a data frame with PCA results
pca_df <- data.frame(pca_result$x)

# Add sample information (assuming Treatment and Species are in the colData of `dds`)
pca_df$Sample <- colnames(top_gene_expression)
pca_df$Treatment <- colData(dds)$Treatment
pca_df$Species <- colData(dds)$Species

# Step 5: Plot the PCA results using ggplot2
library(ggplot2)

# Choose colors and shapes for the plot (modify according to your conditions)
Treatment_colors <- c("OW" = "tomato", "AM" = "cornflowerblue")
shapes_choose <- c("Ahud" = 16, "Atonsa" = 18)

# Create the PCA plot
PCA_interaction <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Species)) + 
  geom_point(size = 5) + 
  scale_shape_manual(values = shapes_choose) + 
  scale_color_manual(values = Treatment_colors) + 
  labs(title = "PCA of Top 5 Differentially Regulated Genes (TreatmentOW.SpeciesAtonsa)",
       x = paste0('PC1: ', round(100 * summary(pca_result)$importance[2, 1], 1), '%'),
       y = paste0('PC2: ', round(100 * summary(pca_result)$importance[2, 2], 1), '%')) + 
  theme_bw(base_size = 16)

# Step 6: Display the plot
print(PCA_interaction)
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/interaction_pca_topgenes", device = "pdf")

#######################################
#
#PCA of all genes
#
#######################################

#PCA!!
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Species"), returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))

Treatment_colors <- c("OW"="tomato", "AM" = "cornflowerblue")
shapes_choose <- c("Ahud" = 16, "Atonsa"= 18)

p <- ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Species))+
  geom_point(size=5)+
  scale_shape_manual(values = shapes_choose)+
  scale_color_manual(values = Treatment_colors)+
  labs(x = paste0('PC1:', percentVar[1], '%'),
       y = paste0('PC2:', percentVar[2], '%'))+
  theme_bw(base_size = 16)
p
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/pca_allgenes", device = "pdf")

library(gridExtra)

combined_plot_PCA <- grid.arrange(PCA_species, PCA_treatment, PCA_interaction,ncol=3 )
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(scales)
library(pheatmap)

options(bitmapType = "cairo")

setwd("/gpfs1/cl/pbio3990/GroupProjects/Acartia/")

################################################
#
#normalize hudsonica data 
#
###############################################
ahuddata <- read.delim("ahud_samples_R.txt")
ahuddata_subset <- subset(ahuddata,treatment=="AM" & generation=="F0")
ahud_subset <- subset(ahuddata, subset = (treatment %in% c("AM","OW")& generation == "F0"))


countsTable <- read.table("/gpfs1/cl/pbio3990/GroupProjects/Acartia/salmon.isoform.counts.matrix.filteredAssembly",
                          header = TRUE, row.names = 1)

##########################
#
#subset the tonsa data files
#
#########################

tonsadata <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                        header = TRUE,
                        stringsAsFactors = TRUE,
                        row.names = 1)
head(tonsadata)
treatment <- tonsadata[which(tonsadata$FinalTemp =="BASE"), ]
print(treatment)

#keeping these colums because they are the correct conditions for this study 
columns_to_keep <- grep("^AA_F0|HA_F0", colnames(countsTable), value=TRUE)
subset.counts_hudsonica <- countsTable[, columns_to_keep]

#round the subsetted data becase DeSeq does not like decimals
subset.counts_hudsonica <- round(subset.counts_hudsonica)

head(subset.counts_hudsonica)


################################################
#
#subsetting the tonsa counts matrix
#
###############################################
countsTable2 <- read.table("/gpfs1/cl/pbio3990/GroupProjects/Acartia/tonsa_counts.txt",
                           header = TRUE, row.names = 1)
subset.counts_tonsa <- countsTable2[, c("N1C3","N2C2","N1C4","N1C2","N2C3","N2C1")]


#round the subsetted data becase DeSeq does not like decimals
subset.counts_tonsaround<- round(subset.counts_tonsa)

head(subset.counts_tonsaround)


##################################################
#
# filtering BLAST query genes for best evalue 
#
##################################################

library(dplyr)

blast_data <- read.table("acartia_blast_max1", header = FALSE, sep = "\t", 
                         stringsAsFactors = FALSE)

colnames(blast_data) <- c("query", "subject", "identity", "length", "mismatch", 
                          "gapopen", "q.start", "q.end", "s.start", "s.end", 
                          "evalue", "bitscore")
best_blast_matches <- blast_data %>%
  group_by(query)%>%
  filter(evalue == min(evalue))%>%
  ungroup()

best_blast_matches <- best_blast_matches %>%
  group_by(query) %>%
  slice(1) %>%
  ungroup()

print(best_blast_matches)


############################################################
#
#call genes from counts matrix and BLAST 
#merge the filtered BLAST file with the counts matricies 
#
###########################################################

merged_ahud1 <- merge(best_blast_matches, subset.counts_hudsonica, by.x = "query", by.y = "row.names", all = TRUE)

merged_tonsa <- merge(best_blast_matches, subset.counts_tonsaround, by.x = "subject", by.y = "row.names", all = TRUE)

merged_all <- merge (merged_ahud1, merged_tonsa, by.x = "query", by.y = "query")

############################################################
#
# subset merged table
#
#############################################################

subset.merged_all <- merged_all[, c("subject.x","query","AA_F0_Rep1_","AA_F0_Rep2_","AA_F0_Rep3_","HA_F0_Rep1_","HA_F0_Rep2_","HA_F0_Rep3_","N2C2","N1C4","N1C2","N2C3","N2C1","N1C3")]

subset.merged_all$combined <- paste(subset.merged_all[,1], subset.merged_all[,2], sep = "_")

write.csv(subset.merged_all, "/gpfs1/cl/pbio3990/GroupProjects/Acartia/subset.merged_all", row.names = FALSE)

#make column 1 the row names 
rownames(subset.merged_all) <- subset.merged_all$combined
#merge the names of the genes so they are all unique and then make them the row names 

subset.merged_all <- subset.merged_all[, -c(1,2, ncol(subset.merged_all))]

################################
#
#Normalize subsetted merged file
#
##############################
head(subset.merged_all)
dim(subset.merged_all)
# 5669     12
# subset.merged_all[acartia_counts < 0] <- 0
# any(acartia_counts < 0)
# which(acartia_counts < 0, arr.ind = TRUE)

acartia_conds <- read.delim("~/Projects/eco_genomics_2/final_project/scripts/acartia_conditions.txt")
head(acartia_conds)

#create dds object 
dds <- DESeqDataSetFromMatrix(countData = subset.merged_all, colData = acartia_conds,
                              design = ~ Treatment + Species + Treatment:Species)
dds <- DESeq(dds)

#final normalized counts matrix 
normalized_counts <- counts(dds, normalized = TRUE)
file.path <- "/gpfs1/cl/pbio3990/GroupProjects/Acartia/normalized_counts"
write.csv(normalized_counts, file = file.path, row.names = TRUE)

head(normalized_counts)
dim(normalized_counts)
#[1] 5669   12

#DESeq results, use these names to run contrasts 
resultsNames(dds)
#[1] "Intercept"                 "Treatment_OW_vs_AM"        "Species_Atonsa_vs_Ahud"    "TreatmentOW.SpeciesAtonsa"

###################################################
#
#Constrasts and Results!!
#
###################################################
#
#Contrast 1: species
#
##################################################

res_species <- results(dds, name = "Species_Atonsa_vs_Ahud", alpha = 0.05)
res_species <- res_species[order(res_species$padj),] #order by significance 
head(res_species)
summary(res_species)
# out of 5655 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 2601, 46%
#LFC < 0 (down)     : 1901, 34%
#outliers [1]       : 70, 1.2%
#low counts [2]     : 0, 0%
#(mean count < 0

#look at counts of a specific top gene that we're interested in to validate that the model is working 
d <- plotCounts(dds, gene = "TRINITY_DN134966_c2_g1_i1_TRINITY_DN369_c1_g1::TRINITY_DN369_c1_g1_i7::g.2127::m.2127", 
                int=(c("Treatment", "Species")), returnData = TRUE)
d
#count number in each of our samples in a specific gene, plot these data 

p <- ggplot(d, aes(x=Treatment, y= count, color=Treatment, shape=Species))+
  theme_minimal()+theme(text=element_text(size=20), panel.grid.major=element_line(color="grey"))
p <- p+geom_point(position=position_jitter(w=0.2,h=0), size=3)
p

#
#PCA top 5 differentially expressed genes
#

# Step 1: Extract the top 5 differentially regulated genes from your DESeq results
res_species <- results(dds, name = "Species_Atonsa_vs_Ahud", alpha = 0.05)

# Order by padj (adjusted p-value) and select top 5 most significant genes
res_species <- res_species[order(res_species$padj), ]

# Extract the top 5 genes with the smallest adjusted p-value (padj) and a log2FoldChange > 1 (optional)
top_genes <- head(res_species, 5)

# Step 2: Extract expression data for the top 5 genes
# Extract normalized counts for the top genes from your DESeqDataSet
top_gene_expression <- counts(dds, normalized = TRUE)[rownames(top_genes), ]

# Step 3: Run PCA on the expression data for these top genes (transpose to get samples as rows)
pca_result <- prcomp(t(top_gene_expression), scale. = TRUE)

# Step 4: Create a data frame with PCA results
pca_df <- data.frame(pca_result$x)

# Add sample information (assuming Treatment and Species are in the colData of `dds`)
pca_df$Sample <- colnames(top_gene_expression)
pca_df$Treatment <- colData(dds)$Treatment
pca_df$Species <- colData(dds)$Species

# Step 5: Plot the PCA results using ggplot2
library(ggplot2)

# Choose colors and shapes for the plot (modify according to your conditions)
Treatment_colors <- c("OW" = "tomato", "AM" = "cornflowerblue")
shapes_choose <- c("Ahud" = 16, "Atonsa" = 18)

# Create the PCA plot
PCA_species <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Species)) + 
  geom_point(size = 5) + 
  scale_shape_manual(values = shapes_choose) + 
  scale_color_manual(values = Treatment_colors) + 
  labs(title = "PCA of Top 5 Differentially Regulated Genes (Species: Tonsa vs Hudsonica)",
       x = paste0('PC1: ', round(100 * summary(pca_result)$importance[2, 1], 1), '%'),
       y = paste0('PC2: ', round(100 * summary(pca_result)$importance[2, 2], 1), '%')) + 
  theme_bw(base_size = 16)

# Step 6: Display the plot
print(PCA_species)

#save the plot
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/species_pca_topgenes", device = "pdf")

###########################################
#
# heatmap of top 5 genes 
#
###########################################

vsd <- vst(dds, blind=FALSE)

topgenes <- head(rownames(res_species), 100)
mat <- assay(vsd)[topgenes,]
mat <- mat-rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("Treatment", "Species")])
heatmap_species <- pheatmap(mat,annotation_col = df, show_rownames = FALSE, cluster_cols = F, cluster_rows = T)

ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/species_heatmap_topgenes", device = "pdf")

########################################
#
#second contrast, treatment
#
########################################
res_treatment <- results(dds, name = "Treatment_OW_vs_AM", alpha = 0.05)
res_treatment <- res_treatment[order(res_treatment$padj),]
head(res_treatment)
summary(res_treatment)
# out of 5655 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 454, 8%
# LFC < 0 (down)     : 362, 6.4%
# outliers [1]       : 70, 1.2%
# low counts [2]     : 0, 0%
# (mean count < 0)

#look at counts of a specific top gene that we're interested in to validate that the model is working 
d <- plotCounts(dds, gene = "TRINITY_DN137257_c0_g2_i1_TRINITY_DN502_c0_g1::TRINITY_DN502_c0_g1_i1::g.3527::m.3527", 
                int=(c("Treatment", "Species")), returnData = TRUE)
d
#count number in each of our samples in a specific gene, plot these data 

p <- ggplot(d, aes(x=Treatment, y=count, color=Treatment, shape=Species))+
  theme_minimal()+theme(text=element_text(size=20), panel.grid.major=element_line(color="grey"))
p <- p+geom_point(position=position_jitter(w=0.2,h=0), size=3)
p

#
#PCA of top differentially regulated genes 
#

# Step 1: Extract the top 5 differentially regulated genes from your DESeq results
res_treatment <- results(dds, name = "Treatment_OW_vs_AM", alpha = 0.05)

# Order by padj (adjusted p-value) and select top 5 most significant genes
res_treatment <- res_treatment[order(res_treatment$padj), ]

# Extract the top 5 genes with the smallest adjusted p-value (padj) and a log2FoldChange > 1 (optional)
top_genes <- head(res_treatment, 5)

# Step 2: Extract expression data for the top 5 genes
# Extract normalized counts for the top genes from your DESeqDataSet
top_gene_expression <- counts(dds, normalized = TRUE)[rownames(top_genes), ]

# Step 3: Run PCA on the expression data for these top genes (transpose to get samples as rows)
pca_result <- prcomp(t(top_gene_expression), scale. = TRUE)

# Step 4: Create a data frame with PCA results
pca_df <- data.frame(pca_result$x)

# Add sample information (assuming Treatment and Species are in the colData of `dds`)
pca_df$Sample <- colnames(top_gene_expression)
pca_df$Treatment <- colData(dds)$Treatment
pca_df$Species <- colData(dds)$Species

# Step 5: Plot the PCA results using ggplot2
library(ggplot2)

# Choose colors and shapes for the plot (modify according to your conditions)
Treatment_colors <- c("OW" = "tomato", "AM" = "cornflowerblue")
shapes_choose <- c("Ahud" = 16, "Atonsa" = 18)

# Create the PCA plot
PCA_treatment <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Species)) + 
  geom_point(size = 5) + 
  scale_shape_manual(values = shapes_choose) + 
  scale_color_manual(values = Treatment_colors) + 
  labs(title = "PCA of Top 5 Differentially Regulated Genes (Treatment_OW_vs_AM)",
       x = paste0('PC1: ', round(100 * summary(pca_result)$importance[2, 1], 1), '%'),
       y = paste0('PC2: ', round(100 * summary(pca_result)$importance[2, 2], 1), '%')) + 
  theme_bw(base_size = 16)

# Step 6: Display the plot
print(PCA_treatment)
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/treatment_pca_topgenes", device = "pdf")

##################################################
#
# heatmap of top 5 differentially expressed genes
#
##################################################

topgenes <- head(rownames(res_treatment), 100)
mat <- assay(vsd)[topgenes,]
mat <- mat-rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("Treatment", "Species")])
heatmap_treatment <- pheatmap(mat,annotation_col = df, show_rownames = FALSE, cluster_cols = F, cluster_rows = T)
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/treatment_heatmap_topgenes", device = "pdf")

#############################################################
#
# third contrast, interaction of species and treatment
#
#############################################################

res_interaction <- results(dds, name = "TreatmentOW.SpeciesAtonsa", alpha = 0.05)
res_interaction <- res_interaction[order(res_interaction$padj),]
head(res_interaction)
summary(res_interaction)
#out of 5655 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 287, 5.1%
# LFC < 0 (down)     : 411, 7.3%
# outliers [1]       : 70, 1.2%
# low counts [2]     : 0, 0%
# (mean count < 0)

#look at counts of a specific top gene that we're interested in to validate that the model is working 
d <- plotCounts(dds, gene = "TRINITY_DN137257_c0_g2_i1_TRINITY_DN502_c0_g1::TRINITY_DN502_c0_g1_i1::g.3527::m.3527", 
                int=(c("Treatment", "Species")), returnData = TRUE)
d
#count number in each of our samples in a specific gene, plot these data 

p <- ggplot(d, aes(x=Treatment, y=count, color=Treatment, shape=Species))+
  theme_minimal()+theme(text=element_text(size=20), panel.grid.major=element_line(color="grey"))
p <- p+geom_point(position=position_jitter(w=0.2,h=0), size=3)
p

################################################
#
#heatmap of top 5 differentially expressed genes
#
###############################################

library(pheatmap)
vsd <- vst(dds, blind=FALSE)

topgenes <- head(rownames(res_interaction), 100)
mat <- assay(vsd)[topgenes,]
mat <- mat-rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("Treatment", "Species")])
heatmap_interaction <- pheatmap(mat,annotation_col = df, show_rownames = FALSE, cluster_cols = F, cluster_rows = T)
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/interaction_heatmap_topgenes", device = "pdf")

#PCA of top differentially regulated genes 

# Step 1: Extract the top 5 differentially regulated genes from your DESeq results
res_interaction <- results(dds, name = "TreatmentOW.SpeciesAtonsa", alpha = 0.05)

# Order by padj (adjusted p-value) and select top 5 most significant genes
res_interaction <- res_interaction[order(res_interaction$padj), ]

# Extract the top 5 genes with the smallest adjusted p-value (padj) and a log2FoldChange > 1 (optional)
top_genes <- head(res_interaction, 5)

# Step 2: Extract expression data for the top 5 genes
# Extract normalized counts for the top genes from your DESeqDataSet
top_gene_expression <- counts(dds, normalized = TRUE)[rownames(top_genes), ]

# Step 3: Run PCA on the expression data for these top genes (transpose to get samples as rows)
pca_result <- prcomp(t(top_gene_expression), scale. = TRUE)

# Step 4: Create a data frame with PCA results
pca_df <- data.frame(pca_result$x)

# Add sample information (assuming Treatment and Species are in the colData of `dds`)
pca_df$Sample <- colnames(top_gene_expression)
pca_df$Treatment <- colData(dds)$Treatment
pca_df$Species <- colData(dds)$Species

# Step 5: Plot the PCA results using ggplot2
library(ggplot2)

# Choose colors and shapes for the plot (modify according to your conditions)
Treatment_colors <- c("OW" = "tomato", "AM" = "cornflowerblue")
shapes_choose <- c("Ahud" = 16, "Atonsa" = 18)

# Create the PCA plot
PCA_interaction <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Species)) + 
  geom_point(size = 5) + 
  scale_shape_manual(values = shapes_choose) + 
  scale_color_manual(values = Treatment_colors) + 
  labs(title = "PCA of Top 5 Differentially Regulated Genes (TreatmentOW.SpeciesAtonsa)",
       x = paste0('PC1: ', round(100 * summary(pca_result)$importance[2, 1], 1), '%'),
       y = paste0('PC2: ', round(100 * summary(pca_result)$importance[2, 2], 1), '%')) + 
  theme_bw(base_size = 16)

# Step 6: Display the plot
print(PCA_interaction)
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/interaction_pca_topgenes", device = "pdf")

#######################################
#
#PCA of all genes
#
#######################################

#PCA!!
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Species"), returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))

Treatment_colors <- c("OW"="tomato", "AM" = "cornflowerblue")
shapes_choose <- c("Ahud" = 16, "Atonsa"= 18)

p <- ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Species))+
  geom_point(size=5)+
  scale_shape_manual(values = shapes_choose)+
  scale_color_manual(values = Treatment_colors)+
  labs(x = paste0('PC1:', percentVar[1], '%'),
       y = paste0('PC2:', percentVar[2], '%'))+
  theme_bw(base_size = 16)
p
ggsave("/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/pca_allgenes", device = "pdf")

library(gridExtra)

combined_plot_PCA <- grid.arrange(PCA_species, PCA_treatment, PCA_interaction,ncol=3 )
pdf(file = "/gpfs1/cl/pbio3990/GroupProjects/Acartia/Acartia_figures/combined_plot_PCA.pdf",
    width = 600, height = 350)

combined_plot_heatmap <- grid.arrange(heatmap_species, heatmap_treatment, heatmap_interaction, ncol=3)



























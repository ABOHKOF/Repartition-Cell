############################################################
# Repartition of Immune Cells - Seurat
############################################################

# ===============================
# 1. Libraries
# ===============================

library(Seurat)
library(patchwork)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(here)

# ===============================
# 2. Data import
# ===============================

data <- readRDS(
  here("data/raw/atlas_for_deconv_CD10neg_annot_figs2.rds")
)

# ===============================
# 3. Subset immune cells
# ===============================

data_immuno <- subset(
  data,
  subset = annot_atlas_low == "Immune"
)

# ===============================
# 4. Normalization & feature selection
# ===============================

data_immuno <- NormalizeData(
  data_immuno,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

data_immuno <- FindVariableFeatures(
  data_immuno,
  selection.method = "vst",
  nfeatures = 2000
)

all.genes <- rownames(data_immuno)
top20 <- head(VariableFeatures(data_immuno), 20)

VariableFeaturePlot(data_immuno) +
  LabelPoints(points = top20, repel = TRUE)

# ===============================
# 5. Scaling, PCA, clustering
# ===============================

data_immuno <- ScaleData(data_immuno, features = all.genes)

data_immuno <- RunPCA(data_immuno, features = VariableFeatures(data_immuno))

data_immuno <- FindNeighbors(data_immuno, dims = 1:20)
data_immuno <- FindClusters(data_immuno, resolution = 0.5)

data_immuno <- RunUMAP(
  data_immuno,
  dims = 1:20,
  n.neighbors = 30
)

# ===============================
# 6. UMAP visualisation
# ===============================

DimPlot(data_immuno, reduction = "umap", label = TRUE)
DimPlot(data_immuno, group.by = "DiseaseID.x")
DimPlot(data_immuno, group.by = "annot_atlas", label = TRUE)
DimPlot(data_immuno, group.by = "patientID", label = TRUE)

# ===============================
# 7. Cell proportion per condition
# ===============================

meta_condition <- data_immuno@meta.data %>%
  select(DiseaseID.x, annot_atlas)

cell_counts <- meta_condition %>%
  group_by(DiseaseID.x, annot_atlas) %>%
  summarise(n_cells = n(), .groups = "drop")

cell_props <- cell_counts %>%
  group_by(DiseaseID.x) %>%
  mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells * 100
  )

ggplot(
  cell_props,
  aes(x = annot_atlas, y = proportion, fill = DiseaseID.x)
) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_classic() +
  labs(
    x = "",
    y = "Proportion of cells (%)",
    fill = "Condition"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ===============================
# 8. Cell proportion per patient
# ===============================

meta_patient <- data_immuno@meta.data %>%
  select(patientID, annot_atlas)

cell_counts_patient <- meta_patient %>%
  group_by(patientID, annot_atlas) %>%
  summarise(n_cells = n(), .groups = "drop")

cell_props_patient <- cell_counts_patient %>%
  group_by(patientID) %>%
  mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells * 100
  )

ggplot(
  cell_props_patient,
  aes(x = annot_atlas, y = proportion, fill = patientID)
) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_classic() +
  labs(
    x = "",
    y = "Proportion of cells (%)",
    fill = "Patient"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

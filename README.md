# Repartition-Cell
# importation des library

library(Seurat)
library(patchwork)
library(dplyr)
library(tidyverse)

# chemin des données
data <- readRDS("C:/Users/RICHMOND/Documents/Richmond/singlecell_rein/atlas_for_deconv_CD10neg_annot_figs2.rds")
View(data@meta.data)
summary(data@meta.data)

data_immuno <- subset(data, subset = annot_atlas_low == "Immune") #considerons que les cellules immunitaires

#normalisation des données
data_immuno <- NormalizeData(data_immuno, normalization.method = "LogNormalize",scale.factor = 10000)
data_immuno <- FindVariableFeatures(data_immuno, selection.method = "vst",  nfeatures = 2000)

all.genes_immuno <- rownames(data_immuno)
top20_immuno <- head(VariableFeatures(data_immuno), 20)
plot1_immuno <- VariableFeaturePlot(data_immuno)
plot1_immuno
plot2_immuno <- LabelPoints(plot = plot1_immuno, points = top20_immuno, repel = TRUE)
plot1_immuno + plot2_immuno

data_immuno <- ScaleData(data_immuno, features = all.genes_immuno)
print(data_immuno[["pca"]], dims = 1:5, nfeatures = 5)

data_immuno <- FindNeighbors(data_immuno, dims = 1:20)
data_immuno <- FindClusters(data_immuno, resolution = 0.5)

head(Idents(data_immuno), 5)

data_immuno <- RunUMAP(data_immuno, dims = 1:20, n.neighbors = 30) #regroupement non superviser en fonction des 20 dimensions


DimPlot(data_immuno, reduction = "umap", label = TRUE)
DimPlot(data_immuno, group.by = "DiseaseID.x")
DimPlot(data_immuno, group.by = "annot_atlas",
        label = TRUE)

DimPlot(data_immuno, group.by = "patientID",
        label = TRUE)
#________________________________________________________________
#           regarder la proportion de distribution des cellules
#----------------------------------------------------------------

# Extraire les métadonnées utiles
meta <- data_immuno@meta.data %>%
  select(DiseaseID.x, annot_atlas)

# Comptage des cellules
cell_counts <- meta %>%
  group_by(DiseaseID.x, annot_atlas) %>%
  summarise(n_cells = n(), .groups = "drop")

# Calcul des proportions par condition
cell_props <- cell_counts %>%
  group_by(DiseaseID.x) %>%
  mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells * 100
  )

cell_props

library(ggplot2)

ggplot(cell_props,
       aes(x = annot_atlas,
           y = proportion,
           fill = DiseaseID.x)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8)) +
  theme_classic() +
  labs(
    x = "",
    y = "Proportion of cells (%)",
    fill = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
#-------------------------------------------
#       pour chaque patient 
#-------------------------------------------
library(ggpubr)
meta_patient <- data_immuno@meta.data %>%
  select(patientID, annot_atlas)

cell_counts_patient <- meta_patient %>%
  group_by(patientID, annot_atlas) %>%
  summarise(n_cells = n(), .groups = "drop")

# Calcul des proportions par condition
cell_props_patient <- cell_counts_patient %>%
  group_by(patientID) %>%
  mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells * 100
  )
cell_props_patient
ggplot(cell_props_patient,
       aes(x = annot_atlas,
           y = proportion,
           fill = patientID)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8))  +
  labs(
    x = "",
    y = "Proportion of cells (%)",
    fill = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


##################################################
# Script for motivation figure - top HVGs and SVGs
# Lukas Weber, Nov 2021
##################################################

# to run in interactive session on cluster:
# module load conda_R/4.1.x
# Rscript filename.R


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggsci)


# ------------
# load results
# ------------

# DLPFC dataset

list_DLPFC <- list(
  spe_DLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG.rds"))), 
  spe_DLPFC_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_DLPFC_HVGs.rds")))
)

for (i in seq_along(list_DLPFC)){
  stopifnot(all(list_DLPFC[[i]]$gene_name == list_DLPFC[[1]]$gene_name))
}

names_DLPFC <- gsub("^spe_", "", names(list_DLPFC))
names_DLPFC

names(list_DLPFC) <- names_DLPFC


# SlideSeqHippo dataset

list_SlideSeqHippo <- list(
  spe_SlideSeqHippo_nnSVG_clusters = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_clusters.rds"))), 
  spe_SlideSeqHippo_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_SlideSeqHippo_HVGs.rds"))), 
  spe_SlideSeqHippo_deviance = rowData(readRDS(here("outputs", "results", "deviance", "spe_SlideSeqHippo_deviance.rds"))), 
  spe_SlideSeqHippo_deviance_clusters = rowData(readRDS(here("outputs", "results", "deviance", "spe_SlideSeqHippo_deviance_clusters.rds"))),
  spe_SlideSeqHippo_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_SlideSeqHippo_SPARKX.rds"))), 
  spe_SlideSeqHippo_SPARKX_clusters = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_SlideSeqHippo_SPARKX_clusters.rds")))
)

for (i in seq_along(list_SlideSeqHippo)){
  stopifnot(all(list_SlideSeqHippo[[i]]$gene_name == list_SlideSeqHippo[[1]]$gene_name))
}

names_SlideSeqHippo <- gsub("^spe_", "", names(list_SlideSeqHippo))
names_SlideSeqHippo

names(list_SlideSeqHippo) <- names_SlideSeqHippo


# ----------------------------------
# identify top genes - DLPFC dataset
# ----------------------------------

# top SVG from nnSVG
ix_top_SVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$rank == 1)
nm_top_SVG <- list_DLPFC[["DLPFC_nnSVG"]]$gene_name[ix_top_SVG]
nm_top_SVG  ## MBP

# top HVG
ix_top_HVG <- which(list_DLPFC[["DLPFC_HVGs"]]$rank == 1)
nm_top_HVG <- list_DLPFC[["DLPFC_HVGs"]]$gene_name[ix_top_HVG]
nm_top_HVG  ## PLP1


# rank of PCP4 in nnSVG results
ix_PCP4_nnSVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$gene_name == "PCP4")
rank_PCP4_nnSVG <- list_DLPFC[["DLPFC_nnSVG"]]$rank[ix_PCP4_nnSVG]
rank_PCP4_nnSVG  ## 81

# rank of PCP4 in HVGs
ix_PCP4_HVGs <- which(list_DLPFC[["DLPFC_HVGs"]]$gene_name == "PCP4")
rank_PCP4_HVGs <- list_DLPFC[["DLPFC_HVGs"]]$rank[ix_PCP4_HVGs]
rank_PCP4_HVGs  ## 88


# rank of HBB in nnSVG results
ix_HBB_nnSVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$gene_name == "HBB")
rank_HBB_nnSVG <- list_DLPFC[["DLPFC_nnSVG"]]$rank[ix_HBB_nnSVG]
rank_HBB_nnSVG  ## 121

# rank of HBB in HVGs
ix_HBB_HVGs <- which(list_DLPFC[["DLPFC_HVGs"]]$gene_name == "HBB")
rank_HBB_HVGs <- list_DLPFC[["DLPFC_HVGs"]]$rank[ix_HBB_HVGs]
rank_HBB_HVGs  ## 25


# rank of IGKC in nnSVG results
ix_IGKC_nnSVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$gene_name == "IGKC")
rank_IGKC_nnSVG <- list_DLPFC[["DLPFC_nnSVG"]]$rank[ix_IGKC_nnSVG]
rank_IGKC_nnSVG  ## 49

# rank of IGKC in HVGs
ix_IGKC_HVGs <- which(list_DLPFC[["DLPFC_HVGs"]]$gene_name == "IGKC")
rank_IGKC_HVGs <- list_DLPFC[["DLPFC_HVGs"]]$rank[ix_IGKC_HVGs]
rank_IGKC_HVGs  ## 13


# rank of NPY in nnSVG results
ix_NPY_nnSVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$gene_name == "NPY")
rank_NPY_nnSVG <- list_DLPFC[["DLPFC_nnSVG"]]$rank[ix_NPY_nnSVG]
rank_NPY_nnSVG  ## 870

# rank of NPY in HVGs
ix_NPY_HVGs <- which(list_DLPFC[["DLPFC_HVGs"]]$gene_name == "NPY")
rank_NPY_HVGs <- list_DLPFC[["DLPFC_HVGs"]]$rank[ix_NPY_HVGs]
rank_NPY_HVGs  ## 51


# --------------------------------------------------
# identify top genes - Slide-seq mouse hippo dataset
# --------------------------------------------------

# Cpne9

# rank of Cpne9 in nnSVG results
ix_Cpne9_nnSVG_clusters <- which(list_SlideSeqHippo[["SlideSeqHippo_nnSVG_clusters"]]$gene_name == "Cpne9")
rank_Cpne9_nnSVG_clusters <- list_SlideSeqHippo[["SlideSeqHippo_nnSVG_clusters"]]$rank[ix_Cpne9_nnSVG_clusters]
rank_Cpne9_nnSVG_clusters  ## 141

# rank of Cpne9 in HVGs
ix_Cpne9_HVGs <- which(list_SlideSeqHippo[["SlideSeqHippo_HVGs"]]$gene_name == "Cpne9")
rank_Cpne9_HVGs <- list_SlideSeqHippo[["SlideSeqHippo_HVGs"]]$rank[ix_Cpne9_HVGs]
rank_Cpne9_HVGs  ## 7975

# rank of Cpne9 in deviance
ix_Cpne9_deviance <- which(list_SlideSeqHippo[["SlideSeqHippo_deviance"]]$gene_name == "Cpne9")
rank_Cpne9_deviance <- list_SlideSeqHippo[["SlideSeqHippo_deviance"]]$rank[ix_Cpne9_deviance]
rank_Cpne9_deviance  ## 2151

# rank of Cpne9 in deviance_clusters
ix_Cpne9_deviance_clusters <- which(list_SlideSeqHippo[["SlideSeqHippo_deviance_clusters"]]$gene_name == "Cpne9")
rank_Cpne9_deviance_clusters <- list_SlideSeqHippo[["SlideSeqHippo_deviance_clusters"]]$rank[ix_Cpne9_deviance_clusters]
rank_Cpne9_deviance_clusters  ## 2969

# rank of Cpne9 in SPARKX
ix_Cpne9_SPARKX <- which(list_SlideSeqHippo[["SlideSeqHippo_SPARKX"]]$gene_name == "Cpne9")
rank_Cpne9_SPARKX <- list_SlideSeqHippo[["SlideSeqHippo_SPARKX"]]$rank[ix_Cpne9_SPARKX]
rank_Cpne9_SPARKX  ## 102

# rank of Cpne9 in SPARKX_clusters
ix_Cpne9_SPARKX_clusters <- which(list_SlideSeqHippo[["SlideSeqHippo_SPARKX_clusters"]]$gene_name == "Cpne9")
rank_Cpne9_SPARKX_clusters <- list_SlideSeqHippo[["SlideSeqHippo_SPARKX_clusters"]]$rank[ix_Cpne9_SPARKX_clusters]
rank_Cpne9_SPARKX_clusters  ## 174


# Rgs14

# rank of Rgs14 in nnSVG results
ix_Rgs14_nnSVG_clusters <- which(list_SlideSeqHippo[["SlideSeqHippo_nnSVG_clusters"]]$gene_name == "Rgs14")
rank_Rgs14_nnSVG_clusters <- list_SlideSeqHippo[["SlideSeqHippo_nnSVG_clusters"]]$rank[ix_Rgs14_nnSVG_clusters]
rank_Rgs14_nnSVG_clusters  ## 19

# rank of Rgs14 in HVGs
ix_Rgs14_HVGs <- which(list_SlideSeqHippo[["SlideSeqHippo_HVGs"]]$gene_name == "Rgs14")
rank_Rgs14_HVGs <- list_SlideSeqHippo[["SlideSeqHippo_HVGs"]]$rank[ix_Rgs14_HVGs]
rank_Rgs14_HVGs  ## 2636

# rank of Rgs14 in deviance
ix_Rgs14_deviance <- which(list_SlideSeqHippo[["SlideSeqHippo_deviance"]]$gene_name == "Rgs14")
rank_Rgs14_deviance <- list_SlideSeqHippo[["SlideSeqHippo_deviance"]]$rank[ix_Rgs14_deviance]
rank_Rgs14_deviance  ## 5025

# rank of Rgs14 in deviance_clusters
ix_Rgs14_deviance_clusters <- which(list_SlideSeqHippo[["SlideSeqHippo_deviance_clusters"]]$gene_name == "Rgs14")
rank_Rgs14_deviance_clusters <- list_SlideSeqHippo[["SlideSeqHippo_deviance_clusters"]]$rank[ix_Rgs14_deviance_clusters]
rank_Rgs14_deviance_clusters  ## 6516

# rank of Rgs14 in SPARKX
ix_Rgs14_SPARKX <- which(list_SlideSeqHippo[["SlideSeqHippo_SPARKX"]]$gene_name == "Rgs14")
rank_Rgs14_SPARKX <- list_SlideSeqHippo[["SlideSeqHippo_SPARKX"]]$rank[ix_Rgs14_SPARKX]
rank_Rgs14_SPARKX  ## 663

# rank of Rgs14 in SPARKX_clusters
ix_Rgs14_SPARKX_clusters <- which(list_SlideSeqHippo[["SlideSeqHippo_SPARKX_clusters"]]$gene_name == "Rgs14")
rank_Rgs14_SPARKX_clusters <- list_SlideSeqHippo[["SlideSeqHippo_SPARKX_clusters"]]$rank[ix_Rgs14_SPARKX_clusters]
rank_Rgs14_SPARKX_clusters  ## 37


# ------------------------------------
# identify SVGs that are not also HVGs
# ------------------------------------

df_all <- full_join(as.data.frame(list_DLPFC[["DLPFC_nnSVG"]]), 
                    as.data.frame(list_DLPFC[["DLPFC_HVGs"]]), 
                    by = "gene_name", 
                    suffix = c(".nnSVG", ".HVGs"))

head(df_all)

df_SVGs_not_HVGs <- df_all[df_all$rank.nnSVG <= 100 & df_all$rank.HVGs > 100, ]
head(df_SVGs_not_HVGs)
View(df_SVGs_not_HVGs)  ## see FTH1

df_SVGs_not_HVGs[df_SVGs_not_HVGs$gene_name == "FTH1", ]  ## rank 14 for nnSVG, rank 1666 for HVGs


View(df_all[df_all$rank.nnSVG <= 50 & df_all$rank.HVGs > 50, ])


# plot FTH1

# load full object
spe_DLPFC_nnSVG <- readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG.rds"))

expr_FTH1 <- counts(spe_DLPFC_nnSVG)[rowData(spe_DLPFC_nnSVG)$gene_name == "FTH1", ]

df_plot <- as.data.frame(cbind(
  spatialCoords(spe_DLPFC_nnSVG), 
  expr_FTH1 = expr_FTH1))

ggplot(df_plot, aes(x = x, y = y, color = expr_FTH1)) + 
  geom_point(size = 0.8) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray95", high = "blue") + 
  ggtitle("FTH1") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())


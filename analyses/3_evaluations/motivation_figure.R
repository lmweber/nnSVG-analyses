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
  DLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG.rds"))), 
  DLPFC_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_DLPFC_HVGs.rds"))), 
  DLPFC_deviance = rowData(readRDS(here("outputs", "results", "deviance", "spe_DLPFC_deviance.rds"))), 
  DLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_DLPFC_SPARKX.rds")))
)

# check row order of gene names is consistent
for (i in seq_along(list_DLPFC)){
  stopifnot(all(list_DLPFC[[i]]$gene_name == list_DLPFC[[1]]$gene_name))
}


# SlideSeqHippo dataset

list_SlideSeqHippo <- list(
  SlideSeqHippo_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG.rds"))), 
  SlideSeqHippo_nnSVG_clusters = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_clusters.rds"))), 
  SlideSeqHippo_nnSVG_onevsall = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_onevsall.rds"))), 
  SlideSeqHippo_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_SlideSeqHippo_HVGs.rds"))), 
  SlideSeqHippo_deviance = rowData(readRDS(here("outputs", "results", "deviance", "spe_SlideSeqHippo_deviance.rds"))), 
  SlideSeqHippo_deviance_clusters = rowData(readRDS(here("outputs", "results", "deviance", "spe_SlideSeqHippo_deviance_clusters.rds"))),
  SlideSeqHippo_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_SlideSeqHippo_SPARKX.rds"))), 
  SlideSeqHippo_SPARKX_clusters = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_SlideSeqHippo_SPARKX_clusters.rds")))
)

# check row order of gene names is consistent
for (i in seq_along(list_SlideSeqHippo)){
  stopifnot(all(list_SlideSeqHippo[[i]]$gene_name == list_SlideSeqHippo[[1]]$gene_name))
}


# -----------------------------------------------
# ranks of top and/or known genes - DLPFC dataset
# -----------------------------------------------

# top genes

# top SVG from nnSVG
ix_top_nnSVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$rank == 1)
nm_top_nnSVG <- list_DLPFC[["DLPFC_nnSVG"]]$gene_name[ix_top_nnSVG]
nm_top_nnSVG  ## MBP
# top HVG
ix_top_HVGs <- which(list_DLPFC[["DLPFC_HVGs"]]$rank == 1)
nm_top_HVGs <- list_DLPFC[["DLPFC_HVGs"]]$gene_name[ix_top_HVGs]
nm_top_HVGs  ## PLP1
# top deviance
ix_top_deviance <- which(list_DLPFC[["DLPFC_deviance"]]$rank == 1)
nm_top_deviance <- list_DLPFC[["DLPFC_deviance"]]$gene_name[ix_top_deviance]
nm_top_deviance  ## MBP
# top SPARK-X
ix_top_SPARKX <- which(list_DLPFC[["DLPFC_SPARKX"]]$rank == 1)
nm_top_SPARKX <- list_DLPFC[["DLPFC_SPARKX"]]$gene_name[ix_top_SPARKX]
nm_top_SPARKX  ## MBP


# known laminar genes

# rank of MOBP in nnSVG
ix_MOBP_nnSVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$gene_name == "MOBP")
rank_MOBP_nnSVG <- list_DLPFC[["DLPFC_nnSVG"]]$rank[ix_MOBP_nnSVG]
rank_MOBP_nnSVG  ## 6
# rank of MOBP in HVGs
ix_MOBP_HVGs <- which(list_DLPFC[["DLPFC_HVGs"]]$gene_name == "MOBP")
rank_MOBP_HVGs <- list_DLPFC[["DLPFC_HVGs"]]$rank[ix_MOBP_HVGs]
rank_MOBP_HVGs  ## 8
# rank of MOBP in deviance
ix_MOBP_deviance <- which(list_DLPFC[["DLPFC_deviance"]]$gene_name == "MOBP")
rank_MOBP_deviance <- list_DLPFC[["DLPFC_deviance"]]$rank[ix_MOBP_deviance]
rank_MOBP_deviance  ## 12
# rank of MOBP in SPARK-X
ix_MOBP_SPARKX <- which(list_DLPFC[["DLPFC_SPARKX"]]$gene_name == "MOBP")
rank_MOBP_SPARKX <- list_DLPFC[["DLPFC_SPARKX"]]$rank[ix_MOBP_SPARKX]
rank_MOBP_SPARKX  ## 12

# rank of PCP4 in nnSVG
ix_PCP4_nnSVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$gene_name == "PCP4")
rank_PCP4_nnSVG <- list_DLPFC[["DLPFC_nnSVG"]]$rank[ix_PCP4_nnSVG]
rank_PCP4_nnSVG  ## 81
# rank of PCP4 in HVGs
ix_PCP4_HVGs <- which(list_DLPFC[["DLPFC_HVGs"]]$gene_name == "PCP4")
rank_PCP4_HVGs <- list_DLPFC[["DLPFC_HVGs"]]$rank[ix_PCP4_HVGs]
rank_PCP4_HVGs  ## 88
# rank of PCP4 in deviance
ix_PCP4_deviance <- which(list_DLPFC[["DLPFC_deviance"]]$gene_name == "PCP4")
rank_PCP4_deviance <- list_DLPFC[["DLPFC_deviance"]]$rank[ix_PCP4_deviance]
rank_PCP4_deviance  ## 76
# rank of PCP4 in SPARK-X
ix_PCP4_SPARKX <- which(list_DLPFC[["DLPFC_SPARKX"]]$gene_name == "PCP4")
rank_PCP4_SPARKX <- list_DLPFC[["DLPFC_SPARKX"]]$rank[ix_PCP4_SPARKX]
rank_PCP4_SPARKX  ## 1186

# rank of SNAP25 in nnSVG
ix_SNAP25_nnSVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$gene_name == "SNAP25")
rank_SNAP25_nnSVG <- list_DLPFC[["DLPFC_nnSVG"]]$rank[ix_SNAP25_nnSVG]
rank_SNAP25_nnSVG  ## 25
# rank of SNAP25 in HVGs
ix_SNAP25_HVGs <- which(list_DLPFC[["DLPFC_HVGs"]]$gene_name == "SNAP25")
rank_SNAP25_HVGs <- list_DLPFC[["DLPFC_HVGs"]]$rank[ix_SNAP25_HVGs]
rank_SNAP25_HVGs  ## 23
# rank of SNAP25 in deviance
ix_SNAP25_deviance <- which(list_DLPFC[["DLPFC_deviance"]]$gene_name == "SNAP25")
rank_SNAP25_deviance <- list_DLPFC[["DLPFC_deviance"]]$rank[ix_SNAP25_deviance]
rank_SNAP25_deviance  ## 16
# rank of SNAP25 in SPARK-X
ix_SNAP25_SPARKX <- which(list_DLPFC[["DLPFC_SPARKX"]]$gene_name == "SNAP25")
rank_SNAP25_SPARKX <- list_DLPFC[["DLPFC_SPARKX"]]$rank[ix_SNAP25_SPARKX]
rank_SNAP25_SPARKX  ## 38


# known non-laminar genes

# rank of HBB in nnSVG
ix_HBB_nnSVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$gene_name == "HBB")
rank_HBB_nnSVG <- list_DLPFC[["DLPFC_nnSVG"]]$rank[ix_HBB_nnSVG]
rank_HBB_nnSVG  ## 121
# rank of HBB in HVGs
ix_HBB_HVGs <- which(list_DLPFC[["DLPFC_HVGs"]]$gene_name == "HBB")
rank_HBB_HVGs <- list_DLPFC[["DLPFC_HVGs"]]$rank[ix_HBB_HVGs]
rank_HBB_HVGs  ## 25
# rank of HBB in deviance
ix_HBB_deviance <- which(list_DLPFC[["DLPFC_deviance"]]$gene_name == "HBB")
rank_HBB_deviance <- list_DLPFC[["DLPFC_deviance"]]$rank[ix_HBB_deviance]
rank_HBB_deviance  ## 28
# rank of HBB in SPARK-X
ix_HBB_SPARKX <- which(list_DLPFC[["DLPFC_SPARKX"]]$gene_name == "HBB")
rank_HBB_SPARKX <- list_DLPFC[["DLPFC_SPARKX"]]$rank[ix_HBB_SPARKX]
rank_HBB_SPARKX  ## 6538

# rank of IGKC in nnSVG
ix_IGKC_nnSVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$gene_name == "IGKC")
rank_IGKC_nnSVG <- list_DLPFC[["DLPFC_nnSVG"]]$rank[ix_IGKC_nnSVG]
rank_IGKC_nnSVG  ## 49
# rank of IGKC in HVGs
ix_IGKC_HVGs <- which(list_DLPFC[["DLPFC_HVGs"]]$gene_name == "IGKC")
rank_IGKC_HVGs <- list_DLPFC[["DLPFC_HVGs"]]$rank[ix_IGKC_HVGs]
rank_IGKC_HVGs  ## 13
# rank of IGKC in deviance
ix_IGKC_deviance <- which(list_DLPFC[["DLPFC_deviance"]]$gene_name == "IGKC")
rank_IGKC_deviance <- list_DLPFC[["DLPFC_deviance"]]$rank[ix_IGKC_deviance]
rank_IGKC_deviance  ## 11
# rank of IGKC in SPARK-X
ix_IGKC_SPARKX <- which(list_DLPFC[["DLPFC_SPARKX"]]$gene_name == "IGKC")
rank_IGKC_SPARKX <- list_DLPFC[["DLPFC_SPARKX"]]$rank[ix_IGKC_SPARKX]
rank_IGKC_SPARKX  ## 3584

# rank of NPY in nnSVG
ix_NPY_nnSVG <- which(list_DLPFC[["DLPFC_nnSVG"]]$gene_name == "NPY")
rank_NPY_nnSVG <- list_DLPFC[["DLPFC_nnSVG"]]$rank[ix_NPY_nnSVG]
rank_NPY_nnSVG  ## 870
# rank of NPY in HVGs
ix_NPY_HVGs <- which(list_DLPFC[["DLPFC_HVGs"]]$gene_name == "NPY")
rank_NPY_HVGs <- list_DLPFC[["DLPFC_HVGs"]]$rank[ix_NPY_HVGs]
rank_NPY_HVGs  ## 51
# rank of NPY in deviance
ix_NPY_deviance <- which(list_DLPFC[["DLPFC_deviance"]]$gene_name == "NPY")
rank_NPY_deviance <- list_DLPFC[["DLPFC_deviance"]]$rank[ix_NPY_deviance]
rank_NPY_deviance  ## 7
# rank of NPY in SPARK-X
ix_NPY_SPARKX <- which(list_DLPFC[["DLPFC_SPARKX"]]$gene_name == "NPY")
rank_NPY_SPARKX <- list_DLPFC[["DLPFC_SPARKX"]]$rank[ix_NPY_SPARKX]
rank_NPY_SPARKX  ## 13583


# bandwidth parameter estimates for genes of interest

list_DLPFC[["DLPFC_nnSVG"]][c(
  ix_top_nnSVG, 
  ix_MOBP_nnSVG, ix_PCP4_nnSVG, ix_SNAP25_nnSVG, 
  ix_HBB_nnSVG, ix_IGKC_nnSVG, ix_NPY_nnSVG
), ]


# plot histogram of bandwidth parameters

df <- as.data.frame(list_DLPFC[["DLPFC_nnSVG"]])
df$is_svg <- df$rank <= 1000
genes_hvgs <- list_DLPFC[["DLPFC_HVGs"]][list_DLPFC[["DLPFC_HVGs"]]$rank <= 1000, ]$gene_id
df$is_hvg <- df$gene_id %in% genes_hvgs
genes_deviance <- list_DLPFC[["DLPFC_deviance"]][list_DLPFC[["DLPFC_deviance"]]$rank <= 1000, ]$gene_id
df$is_deviance <- df$gene_id %in% genes_deviance
genes_sparkx <- list_DLPFC[["DLPFC_SPARKX"]][list_DLPFC[["DLPFC_SPARKX"]]$rank <= 1000, ]$gene_id
df$is_sparkx <- df$gene_id %in% genes_sparkx

df_sub <- df[df$is_svg, ]

ggplot(df, aes(x = phi, group = is_hvg, fill = is_hvg)) + 
  geom_histogram(bins = 40, position = "identity", color = "black", size = 0.25) + 
  scale_fill_manual(values = c("#69b3a2", "#404080")) + 
  labs(x = "phi") + 
  ggtitle("Inverse bandwidth phi, nnSVG, top 1000 SVGs, DLPFC dataset") + 
  theme_bw()

fn <- here("plots", "bandwidth", paste0("bandwidth_nnSVG_top1000SVGs_top1000HVGs"))
ggsave(paste0(fn, ".pdf"), width = 6, height = 4.5)
ggsave(paste0(fn, ".png"), width = 6, height = 4.5)


ggplot(df, aes(x = phi, group = is_deviance, fill = is_deviance)) + 
  geom_histogram(bins = 40, position = "identity", color = "black", size = 0.25) + 
  scale_fill_manual(values = c("#69b3a2", "#404080")) + 
  labs(x = "phi") + 
  ggtitle("Inverse bandwidth phi, nnSVG, top 1000 SVGs, DLPFC dataset") + 
  theme_bw()

fn <- here("plots", "bandwidth", paste0("bandwidth_nnSVG_top1000SVGs_top1000deviance"))
ggsave(paste0(fn, ".pdf"), width = 6, height = 4.5)
ggsave(paste0(fn, ".png"), width = 6, height = 4.5)


ggplot(df, aes(x = phi, group = is_sparkx, fill = is_sparkx)) + 
  geom_histogram(bins = 40, position = "identity", color = "black", size = 0.25) + 
  scale_fill_manual(values = c("#69b3a2", "darkred")) + 
  labs(x = "phi") + 
  ggtitle("Inverse bandwidth phi, nnSVG, top 1000 SVGs, DLPFC dataset") + 
  theme_bw()

fn <- here("plots", "bandwidth", paste0("bandwidth_nnSVG_top1000SVGs_top1000SPARKX"))
ggsave(paste0(fn, ".pdf"), width = 6, height = 4.5)
ggsave(paste0(fn, ".png"), width = 6, height = 4.5)



# -----------------------------------------------------------------
# ranks of top and/or known genes - Slide-seqV2 mouse hippo dataset
# -----------------------------------------------------------------

# Cpne9

# rank of Cpne9 in nnSVG
ix_Cpne9_nnSVG <- which(list_SlideSeqHippo[["SlideSeqHippo_nnSVG"]]$gene_name == "Cpne9")
rank_Cpne9_nnSVG <- list_SlideSeqHippo[["SlideSeqHippo_nnSVG"]]$rank[ix_Cpne9_nnSVG]
rank_Cpne9_nnSVG  ## 144
# rank of Cpne9 in nnSVG_clusters
ix_Cpne9_nnSVG_clusters <- which(list_SlideSeqHippo[["SlideSeqHippo_nnSVG_clusters"]]$gene_name == "Cpne9")
rank_Cpne9_nnSVG_clusters <- list_SlideSeqHippo[["SlideSeqHippo_nnSVG_clusters"]]$rank[ix_Cpne9_nnSVG_clusters]
rank_Cpne9_nnSVG_clusters  ## 141
# rank of Cpne9 in nnSVG_onevsall
ix_Cpne9_nnSVG_onevsall <- which(list_SlideSeqHippo[["SlideSeqHippo_nnSVG_onevsall"]]$gene_name == "Cpne9")
rank_Cpne9_nnSVG_onevsall <- list_SlideSeqHippo[["SlideSeqHippo_nnSVG_onevsall"]]$rank[ix_Cpne9_nnSVG_onevsall]
rank_Cpne9_nnSVG_onevsall  ## 180
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

# rank of Rgs14 in nnSVG
ix_Rgs14_nnSVG <- which(list_SlideSeqHippo[["SlideSeqHippo_nnSVG"]]$gene_name == "Rgs14")
rank_Rgs14_nnSVG <- list_SlideSeqHippo[["SlideSeqHippo_nnSVG"]]$rank[ix_Rgs14_nnSVG]
rank_Rgs14_nnSVG  ## 171
# rank of Rgs14 in nnSVG_clusters
ix_Rgs14_nnSVG_clusters <- which(list_SlideSeqHippo[["SlideSeqHippo_nnSVG_clusters"]]$gene_name == "Rgs14")
rank_Rgs14_nnSVG_clusters <- list_SlideSeqHippo[["SlideSeqHippo_nnSVG_clusters"]]$rank[ix_Rgs14_nnSVG_clusters]
rank_Rgs14_nnSVG_clusters  ## 19
# rank of Rgs14 in nnSVG_onevsall
ix_Rgs14_nnSVG_onevsall <- which(list_SlideSeqHippo[["SlideSeqHippo_nnSVG_onevsall"]]$gene_name == "Rgs14")
rank_Rgs14_nnSVG_onevsall <- list_SlideSeqHippo[["SlideSeqHippo_nnSVG_onevsall"]]$rank[ix_Rgs14_nnSVG_onevsall]
rank_Rgs14_nnSVG_onevsall  ## 127
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


# bandwidth parameter estimates for genes of interest

list_SlideSeqHippo[["SlideSeqHippo_nnSVG_clusters"]][c(
  ix_Cpne9_nnSVG, ix_Rgs14_nnSVG
), ]


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


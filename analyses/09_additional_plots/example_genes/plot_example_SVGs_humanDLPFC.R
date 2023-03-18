#############################
# Script to plot example SVGs
# Lukas Weber, May 2022
#############################

library(SpatialExperiment)
library(STexampleData)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)


# directory to save plots
dir_plots <- here(file.path("plots", "example_genes"))


# ------------------
# Visium human DLPFC
# ------------------

# Maynard and Collado-Torres et al. (2021)

spe <- Visium_humanDLPFC()

# genes of interest
known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")

ix_known <- which(rowData(spe)$gene_name %in% known_genes)
names(ix_known) <- rowData(spe)$gene_name[ix_known]
stopifnot(length(ix_known) == length(known_genes))

ix_known <- ix_known[known_genes]
ix_known

# additional genes
additional_genes <- c("CALM1", "CST3")

ix_additional <- which(rowData(spe)$gene_name %in% additional_genes)
names(ix_additional) <- rowData(spe)$gene_name[ix_additional]
stopifnot(length(ix_additional) == length(additional_genes))

ix_additional <- ix_additional[additional_genes]
ix_additional


# UMI counts
colData(spe)$MOBP <- counts(spe)[ix_known["MOBP"], ]
colData(spe)$PCP4 <- counts(spe)[ix_known["PCP4"], ]
colData(spe)$SNAP25 <- counts(spe)[ix_known["SNAP25"], ]
colData(spe)$HBB <- counts(spe)[ix_known["HBB"], ]
colData(spe)$IGKC <- counts(spe)[ix_known["IGKC"], ]
colData(spe)$NPY <- counts(spe)[ix_known["NPY"], ]

colData(spe)$CALM1 <- counts(spe)[ix_additional["CALM1"], ]
colData(spe)$CST3 <- counts(spe)[ix_additional["CST3"], ]


# expression plots

# known genes

df <- as.data.frame(cbind(colData(spe), spatialCoords(spe))) %>% 
  select(c("barcode_id", "in_tissue", known_genes, 
           "pxl_col_in_fullres", "pxl_row_in_fullres")) %>% 
  filter(in_tissue == 1) %>% 
  pivot_longer(., cols = known_genes, 
               names_to = "gene", values_to = "counts") %>% 
  mutate(gene = factor(gene, levels = known_genes)) %>% 
  mutate(bandwidth = factor(
    ifelse(gene %in% c("MOBP", "PCP4", "SNAP25"), "large length scale", "small length scale")))


ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = counts)) + 
  facet_wrap(~ gene, nrow = 2) + 
  geom_point(size = 0.05) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradientn(trans = "log1p", 
                        colors = c("gray90", mid = "blue", high = "black"), 
                        breaks = c(0, 500, 600), labels = c("0", "", "600")) + 
  ggtitle("Selected SVGs: human DLPFC") + 
  theme_bw() + 
  guides(color = guide_colorbar(ticks = FALSE)) + 
  theme(strip.text = element_text(face = "italic"), 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "example_genes_humanDLPFC")
ggsave(paste0(fn, ".png"), width = 5.75, height = 4.5)


# additional genes

df <- as.data.frame(cbind(colData(spe), spatialCoords(spe))) %>% 
  select(c("barcode_id", "in_tissue", additional_genes, 
           "pxl_col_in_fullres", "pxl_row_in_fullres")) %>% 
  filter(in_tissue == 1) %>% 
  pivot_longer(., cols = additional_genes, 
               names_to = "gene", values_to = "counts") %>% 
  mutate(gene = factor(gene, levels = additional_genes))


ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = counts)) + 
  facet_wrap(~ gene, nrow = 1) + 
  geom_point(size = 0.05) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradientn(trans = "log1p", 
                        colors = c("gray90", high = "blue"), 
                        breaks = c(0, 50), labels = c("0", "50")) + 
  ggtitle("Additional SVGs: human DLPFC") + 
  theme_bw() + 
  guides(color = guide_colorbar(ticks = FALSE)) + 
  theme(strip.text = element_text(face = "italic"), 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "additional_genes_humanDLPFC")
ggsave(paste0(fn, ".png"), width = 4.25, height = 2.5)


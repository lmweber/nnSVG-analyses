#############################
# Script to plot example SVGs
# Lukas Weber, Feb 2022
#############################


library(SpatialExperiment)
library(STexampleData)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)


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


# UMI counts
colData(spe)$MOBP <- counts(spe)[ix_known["MOBP"], ]
colData(spe)$PCP4 <- counts(spe)[ix_known["PCP4"], ]
colData(spe)$SNAP25 <- counts(spe)[ix_known["SNAP25"], ]
colData(spe)$HBB <- counts(spe)[ix_known["HBB"], ]
colData(spe)$IGKC <- counts(spe)[ix_known["IGKC"], ]
colData(spe)$NPY <- counts(spe)[ix_known["NPY"], ]


# expression plots

df <- as.data.frame(cbind(colData(spe), spatialCoords(spe))) %>% 
  select(c("barcode_id", "in_tissue", known_genes, "x", "y")) %>% 
  filter(in_tissue == 1) %>% 
  pivot_longer(., cols = known_genes, 
               names_to = "gene", values_to = "counts") %>% 
  mutate(gene = factor(gene, levels = known_genes))

ggplot(df, aes(x = x, y = y, color = counts)) + 
  facet_wrap(~ gene, nrow = 2) + 
  geom_point(size = 0.05) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradientn(trans = "log1p", 
                        colors = c("gray90", mid = "blue", high = "black"), 
                        breaks = c(0, 500, 600), labels = c("0", "", "600")) + 
  ggtitle("Example SVGs: DLPFC") + 
  theme_bw() + 
  guides(color = guide_colorbar(ticks = FALSE)) + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave(here("plots", "example_svgs", "humanDLPFC_6known.png"), width = 6.25, height = 4.75)


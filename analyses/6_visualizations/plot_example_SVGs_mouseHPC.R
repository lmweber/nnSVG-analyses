#############################
# Script to plot example SVGs
# Lukas Weber, Feb 2022
#############################


library(SpatialExperiment)
library(STexampleData)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)


# -----------------------------
# Slide-seqV2 mouse hippocampus
# -----------------------------

# Stickels et al. (2020) and Cable et al. (2021)

spe <- SlideSeqV2_mouseHPC()

# cell type labels are stored in colData
colData(spe)

# genes of interest
ix_Cpne9 <- which(rowData(spe)$gene_name == "Cpne9")
ix_Rgs14 <- which(rowData(spe)$gene_name == "Rgs14")

colData(spe)$Cpne9 <- counts(spe)[ix_Cpne9, ]
colData(spe)$Rgs14 <- counts(spe)[ix_Rgs14, ]

df <- as.data.frame(cbind(colData(spe), spatialCoords(spe)))
levs <- c("none", sort(unique(df$celltype)))
df$celltype[is.na(df$celltype)] <- "none"
df$celltype <- factor(df$celltype, levels = levs)

pal <- unname(palette.colors(36, "Polychrome 36"))
pal <- pal[c(1, 3:36)]


# plot all cells

# add points in two steps to avoid overplotting NA points
ggplot() + 
  geom_point(data = df[df$celltype == "none", ], 
             aes(x = xcoord, y = ycoord), 
             size = 0.01, color = "gray90") + 
  geom_point(data = df[df$celltype != "none", ], 
             aes(x = xcoord, y = ycoord, color = celltype), 
             size = 0.01) + 
  scale_color_manual(values = pal) + 
  coord_fixed() + 
  ggtitle("Mouse HPC") + 
  guides(color = guide_legend(override.aes = list(size = 2.5))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave(here("plots", "example_svgs", "mouseHPC_celltypes.png"), width = 6.75, height = 5)


# plot expression of genes of interest in CA3 region

# subset CA3 region
df_sub <- 
  filter(df, celltype == "CA3") %>% 
  pivot_longer(., cols = c("Cpne9", "Rgs14"), 
               names_to = "gene", values_to = "counts") %>% 
  mutate(gene = factor(gene, levels = c("Cpne9", "Rgs14")))


ggplot(df_sub, aes(x = xcoord, y = ycoord, color = counts)) + 
  facet_wrap(~ gene, nrow = 2) + 
  geom_point(size = 0.1, alpha = 0.5) + 
  coord_fixed() + 
  scale_color_gradientn(trans = "log1p", 
                        colors = c("gray90", mid = "red", high = "black"), 
                        breaks = c(0, 7, 8), labels = c("0", "", "8")) + 
  ggtitle("Example SVGs: mouse HPC") + 
  theme_bw() + 
  guides(color = guide_colorbar(ticks = FALSE)) + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave(here("plots", "example_svgs", "mouseHPC_known.png"), width = 4, height = 4.5)


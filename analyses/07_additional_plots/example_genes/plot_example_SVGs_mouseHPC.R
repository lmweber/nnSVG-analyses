#############################
# Script to plot example SVGs
# Lukas Weber, May 2022
#############################

library(SpatialExperiment)
library(STexampleData)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)


# directory to save plots
dir_plots <- here(file.path("plots", "example_genes"))


# ---------------------
# Slide-seqV2 mouse HPC
# ---------------------

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


# plot all spots

# add points in two steps to avoid overplotting NA points
ggplot() + 
  geom_point(data = df[df$celltype == "none", ], 
             aes(x = xcoord, y = ycoord), 
             size = 0.01, color = "gray90") + 
  geom_point(data = df[df$celltype != "none", ], 
             aes(x = xcoord, y = ycoord, color = celltype), 
             size = 0.01) + 
  scale_color_manual(values = pal, name = "cell type") + 
  coord_fixed() + 
  ggtitle("Mouse HPC") + 
  guides(color = guide_legend(override.aes = list(size = 2.5))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "mouseHPC_celltypes")
ggsave(paste0(fn, ".png"), width = 6.5, height = 4.75)


# plot expression of example genes in CA3 cell type

# subset CA3 cell type
df_sub <- 
  filter(df, celltype == "CA3") %>% 
  pivot_longer(., cols = c("Cpne9", "Rgs14"), 
               names_to = "gene", values_to = "counts") %>% 
  mutate(gene = factor(gene, levels = c("Cpne9", "Rgs14")))

# add points in two steps to avoid overplotting NA points
ggplot(df_sub, aes(x = xcoord, y = ycoord, color = counts)) + 
  facet_wrap(~ gene, nrow = 2) + 
  geom_point(size = 0.01) + 
  geom_point(data = df_sub[df_sub$counts > 0, ], size = 0.01) + 
  coord_fixed() + 
  scale_color_gradientn(trans = "log1p", 
                        colors = c("gray90", mid = "red", high = "black"), 
                        breaks = c(0, 7, 8), labels = c("0", "", "8")) + 
  ggtitle("Selected SVGs: mouse HPC") + 
  theme_bw() + 
  guides(color = guide_colorbar(ticks = FALSE)) + 
  theme(strip.text = element_text(face = "italic"), 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "example_genes_mouseHPC")
ggsave(paste0(fn, ".png"), width = 3.75, height = 4.25)


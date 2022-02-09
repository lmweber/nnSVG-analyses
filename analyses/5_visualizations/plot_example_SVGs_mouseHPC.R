#############################
# Script to plot example SVGs
# Lukas Weber, Feb 2022
#############################


library(SpatialExperiment)
library(STexampleData)
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
  ggtitle("Slide-seqV2 mouse hippocampus") + 
  guides(color = guide_legend(override.aes = list(size = 2.5))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave(here("plots", "example_svgs", "mousehippo_all.png"), width = 6.75, height = 5)


# plot expression of genes of interest in CA3 region

# subset CA3 region
df_sub <- df[df$celltype == "CA3", ]


ggplot(df_sub, aes(x = xcoord, y = ycoord, color = Cpne9)) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_color_gradient(trans = "sqrt", low = "gray85", high = "red") + 
  ggtitle("Slide-seqV2 mouse hippocampus") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave(here("plots", "example_svgs", "mousehippo_Cpne9.png"), width = 6, height = 3.5)


ggplot(df_sub, aes(x = xcoord, y = ycoord, color = Rgs14)) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_color_gradient(trans = "sqrt", low = "gray85", high = "red") + 
  ggtitle("Slide-seqV2 mouse hippocampus") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave(here("plots", "example_svgs", "mousehippo_Rgs14.png"), width = 6, height = 3.5)



ggplot(df_sub, aes(x = xcoord, y = ycoord, color = Cpne9)) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_color_viridis_c(trans = "sqrt") + 
  ggtitle("Slide-seqV2 mouse hippocampus") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave(here("plots", "example_svgs", "mousehippo_Cpne9.png"), width = 6, height = 3.5)



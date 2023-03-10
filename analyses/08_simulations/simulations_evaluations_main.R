#########################################################
# Simulations: performance evaluations - main simulations
# Lukas Weber, Mar 2023
#########################################################


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


dir_sims <- here("outputs", "simulations")
dir_plots <- here("plots", "simulations")


# ------------
# load results
# ------------

res_list <- readRDS(here(dir_sims, "main", "res_simulations_main.rds"))

length(res_list)

sim_names <- gsub("^sim_", "", names(res_list))
sim_names


# -----------------------------
# calculate performance metrics
# -----------------------------

# compare sensitivity and specificity in each simulation

# p-value thresholds
thresholds <- c(0.01, 0.05, 0.1)

metrics <- c("TPR", "FPR")

res_perf <- data.frame(
  simulation = rep(sim_names, each = length(metrics) * length(thresholds)), 
  metric = rep(rep(metrics, each = length(thresholds)), length(sim_names)), 
  threshold = rep(thresholds, length(sim_names) * length(metrics)), 
  value = rep(NA, length(sim_names) * length(metrics) * length(thresholds))
)

for (s in seq_along(sim_names)) {
  for (p in seq_along(thresholds)) {
    
    is_true <- res_list[[s]]$expressed
    is_pos <- res_list[[s]]$pval <= thresholds[p]
    
    true_pos <- sum(is_true & is_pos)
    false_neg <- sum(is_true & !is_pos)
    false_pos <- sum(!is_true & is_pos)
    true_neg <- sum(!is_true & !is_pos)
    
    # calculate performance metrics
    tpr <- true_pos / (true_pos + false_neg)
    fpr <- false_pos / (false_pos + true_neg)
    
    # store in correct rows in results table
    ix_tpr <- which(
      res_perf$simulation == sim_names[s] & 
      res_perf$metric == "TPR" & 
      res_perf$threshold == thresholds[p]
    )
    ix_fpr <- which(
      res_perf$simulation == sim_names[s] & 
      res_perf$metric == "FPR" & 
      res_perf$threshold == thresholds[p]
    )
    
    res_perf[ix_tpr, "value"] <- tpr
    res_perf[ix_fpr, "value"] <- fpr
  }
}

res_perf


# ----------------
# plot performance
# ----------------

# reshape data frame
df_plot <- res_perf
df_plot$simulation <- factor(df_plot$simulation, levels = sim_names)
df_plot$metric <- factor(df_plot$metric, levels = metrics)
df_plot$threshold <- factor(df_plot$threshold, levels = thresholds)

df_plot <- pivot_wider(df_plot, names_from = metric, values_from = value)


# create plot
ggplot(df_plot, aes(x = FPR, y = TPR, color = threshold)) + 
  facet_wrap(~ simulation, nrow = 3) + 
  geom_point(size = 2.25) + 
  scale_color_manual(values = c("orange1", "firebrick1", "purple3"), 
                     name = "p-value\nthreshold") + 
  #geom_vline(xintercept = thresholds[1], color = "orange1", lty = 2) + 
  #geom_vline(xintercept = thresholds[2], color = "firebrick1", lty = 2) + 
  #geom_vline(xintercept = thresholds[3], color = "purple3", lty = 2) + 
  xlim(c(0, 0.1)) + 
  ylim(c(0, 1)) + 
  ggtitle("nnSVG performance") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

fn <- file.path(dir_plots, "simulations_performance")
ggsave(paste0(fn, ".pdf"), width = 6, height = 5.5)
ggsave(paste0(fn, ".png"), width = 6, height = 5.5)


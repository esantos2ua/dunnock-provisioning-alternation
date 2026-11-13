# Figures for the manuscript (vector PDF, read by Typst)
#
# Colours: categorical slots 1-2 of the validated reference palette (dataviz
# skill), monogamous = blue, polyandrous = orange; nulls in neutral grey.
# Behavioral Ecology: all labels in black; Helvetica is used as the metric
# equivalent of Arial (the base pdf device has no Arial).
# Run from the repo root: Rscript analysis/05_figures.R

suppressMessages({ library(ggplot2); library(patchwork) })
source("R/models.R")

m <- readRDS("output/results/models.rds")
dir.create("manuscript/figures", recursive = TRUE, showWarnings = FALSE)

col_system <- c(monogamous = "#2a78d6", polyandrous = "#eb6834")
ink <- "black"; ink2 <- "#52514e"; grid_col <- "#e6e5e0"; null_col <- "#8f8e89"
scen_system <- c(mono = "monogamous", poly_fm = "polyandrous", poly_fab = "polyandrous",
                 poly_ab = "polyandrous")
scen_short <- c(mono = "Pairs: female vs male", poly_fm = "Polyandrous: female vs males",
                poly_fab = "Trios: among individuals", poly_ab = "Polyandrous: alpha vs beta")

theme_paper <- theme_minimal(base_size = 8.5, base_family = "Helvetica") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.3),
        axis.text = element_text(colour = ink), axis.title = element_text(colour = ink),
        strip.text = element_text(colour = ink, face = "bold", hjust = 0),
        legend.position = "none", plot.margin = margin(4, 6, 4, 4))

# ---- Figure 1: observed vs null alternation rates -----------------------------
v <- m$vs_null
obs <- m$observed_pooled
f1 <- rbind(
  data.frame(scenario = obs$scenario, dataset = "Observed (A)", est = obs$rate,
             lo = obs$lo, hi = obs$hi, kind = "observed"),
  data.frame(scenario = v$scenario, dataset = null_labels[v$null], est = v$null_mean,
             lo = ifelse(v$null == "sorted", NA, v$null_lo),
             hi = ifelse(v$null == "sorted", NA, v$null_hi), kind = "null")
)
f1$dataset <- factor(f1$dataset, levels = rev(c("Observed (A)", null_labels)))
f1$panel <- factor(scen_short[f1$scenario], levels = scen_short)
f1$system <- scen_system[f1$scenario]
ref <- data.frame(panel = factor(scen_short[obs$scenario], levels = scen_short), x = obs$rate)

p1 <- ggplot(f1, aes(est, dataset)) +
  geom_vline(data = ref, aes(xintercept = x), linetype = "22", colour = ink2, linewidth = 0.35) +
  geom_linerange(aes(xmin = lo, xmax = hi, colour = ifelse(kind == "observed", system, "null")),
                 linewidth = 0.8, na.rm = TRUE) +
  geom_point(aes(fill = ifelse(kind == "observed", system, "null")), shape = 21,
             colour = "white", stroke = 0.6, size = 2.4) +
  scale_colour_manual(values = c(col_system, null = null_col)) +
  scale_fill_manual(values = c(col_system, null = null_col)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~panel, ncol = 2) +
  labs(x = "Alternation rate (alternations / maximum possible)", y = NULL) +
  theme_paper
ggsave("manuscript/figures/fig1_observed_vs_null.pdf", p1, width = 174, height = 95,
       units = "mm", device = grDevices::pdf, useDingbats = FALSE)
ggsave("manuscript/figures/fig1_observed_vs_null.png", p1, width = 174, height = 95,
       units = "mm", dpi = 300)

# ---- Figure 2: odds ratios, observed vs each null ------------------------------
v$panel <- factor(scen_short[v$scenario], levels = scen_short)
v$null_lab <- factor(null_labels[v$null], levels = rev(null_labels))
v$system <- scen_system[v$scenario]
p2 <- ggplot(v, aes(exp(log_or), null_lab)) +
  geom_vline(xintercept = 1, colour = ink2, linewidth = 0.35) +
  geom_linerange(aes(xmin = exp(lo), xmax = exp(hi), colour = system), linewidth = 0.8) +
  geom_point(aes(fill = system), shape = 21, colour = "white", stroke = 0.6, size = 2.4) +
  scale_colour_manual(values = col_system) + scale_fill_manual(values = col_system) +
  scale_x_log10(breaks = c(0.8, 1, 1.25, 1.5, 2, 2.5, 3)) +
  facet_wrap(~panel, ncol = 2) +
  labs(x = "Odds ratio, observed vs null alternation (log scale)", y = NULL) +
  theme_paper
ggsave("manuscript/figures/fig2_odds_ratios.pdf", p2, width = 174, height = 85,
       units = "mm", device = grDevices::pdf, useDingbats = FALSE)
ggsave("manuscript/figures/fig2_odds_ratios.png", p2, width = 174, height = 85,
       units = "mm", dpi = 300)

# ---- Figure 3: provisioning rates by role --------------------------------------
d <- readRDS("data/derived/analysis_data.rds")
birds <- bird_visit_table(d$watches)
birds$rate <- birds$visits / birds$duration_h
rt <- m$rate_tab
rt$group <- factor(rt$group, levels = levels(birds$group))
rt$system <- ifelse(grepl("pair", rt$group), "monogamous", "polyandrous")
birds$system_c <- as.character(birds$system)
set.seed(1)
p3 <- ggplot() +
  geom_jitter(data = birds, aes(group, rate, colour = system_c), width = 0.18, height = 0,
              alpha = 0.25, size = 0.7, stroke = 0) +
  # Capped bars and a small dot, so the narrow beta-male interval is not hidden
  geom_errorbar(data = rt, aes(group, ymin = lo, ymax = hi), colour = ink, linewidth = 0.6,
                width = 0.25) +
  geom_point(data = rt, aes(group, rate, fill = system), shape = 21, colour = ink,
             stroke = 0.4, size = 1.8) +
  scale_colour_manual(values = col_system) + scale_fill_manual(values = col_system) +
  coord_cartesian(ylim = c(0, quantile(birds$rate, 0.99))) +
  labs(x = NULL, y = "Provisioning visits per hour") +
  theme_paper + theme(panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3))
ggsave("manuscript/figures/fig3_provisioning_rates.pdf", p3, width = 114, height = 70,
       units = "mm", device = grDevices::pdf, useDingbats = FALSE)
ggsave("manuscript/figures/fig3_provisioning_rates.png", p3, width = 114, height = 70,
       units = "mm", dpi = 300)
message("Figures written to manuscript/figures/")

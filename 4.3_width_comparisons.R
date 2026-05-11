library(tidyverse)
library(lubridate)
library(dplyr)
library(ggtext)

# ---------------------------------------------------------------------------------------------------------------------------
# SWOT WSE validation
# ---------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------
# Node level
# ---------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------
# READ IN DATA
# ---------------------------------------------------------------------------------------------------------------------------

node_SWOT_ortho_vC <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16/node_width_SWOT_Ortho.csv') %>%
  mutate(source = "PIC0") %>%
  rename(old_node_id = node_id) %>%
  filter(abs(residuals) < 1500) %>%
  filter(dark_frac < 0.5)

node_SWOT_ortho_vD <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v17b/node_width_SWOT_Ortho.csv') %>%
  mutate(source = "PGD0") %>%
  filter(abs(residuals) < 1500) %>%
  filter(dark_frac < 0.5)


# ---------------------------------------------------------------------------------------------------------------------------
# Get all data to the same SWORD version (v17b)
# ---------------------------------------------------------------------------------------------------------------------------

# SWORD translator to change Version C data to SWORD v17b naming convention
SWORD_translator <- read_csv("/Users/camryn/Desktop/SWORD_translation/NA_NodeIDs_v17b_vs_v16.csv")

# translate the vC SWORD v16 data to SWORD v17b
node_SWOT_ortho_vC <- node_SWOT_ortho_vC %>%
  left_join(SWORD_translator %>% 
              select(v16_node_id, v17_node_id),
            by = c("old_node_id" = "v16_node_id")) %>%
  rename(node_id = v17_node_id)

# merge all dataframes together
node_SWOT_ortho <- bind_rows(node_SWOT_ortho_vC, node_SWOT_ortho_vD) %>%
  filter(dark_frac < 0.5)

# compute version inclusion
# -1 = only in vC, 0 = both, 1 = only in vD
all_nodes <- node_SWOT_ortho %>%
  distinct(node_id, source) %>%         
  group_by(node_id) %>%
  summarise(
    has_RiverSP   = any(source == "PIC0"),
    has_RiverTile = any(source == "PGD0"),
    .groups = "drop") %>%
  mutate(
    version_inclusion = case_when(has_RiverSP & has_RiverTile ~ 0L, has_RiverSP & !has_RiverTile ~ -1L, !has_RiverSP & has_RiverTile ~ 1L, TRUE ~ NA_integer_)) %>%
  select(node_id, version_inclusion)

node_SWOT_ortho <- node_SWOT_ortho %>%
  left_join(all_nodes, by = "node_id")


# ---------------------------------------------------------------------------------------------------------------------------
# TABLES -- SUMMARY STATS
# ---------------------------------------------------------------------------------------------------------------------------

# ABSOLUTE NODE WIDTH TABLE BY VERSION
# -----------------------------------------------------
table_absolute_node_width <- node_SWOT_ortho %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_abs_68ile = round(quantile(abs(residuals), 0.68, na.rm = TRUE), 1),
    error_abs_50ile = round(quantile(abs(residuals), 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias), na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE), 1),
    error_perdiff_68ile = round(quantile(percent_diff, 0.68, na.rm = TRUE), 2),
    error_perdiff_50ile = round(quantile(percent_diff, 0.50, na.rm = TRUE), 2),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_ortho %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(width, ortho_width_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(width, ortho_width_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_node_width <- table_absolute_node_width %>%
  left_join(cor_table, by = "source")


# ABSOLUTE NODE WIDTH TABLE BY RIVER vD
# -----------------------------------------------------
table_absolute_node_width <- node_SWOT_ortho_vD %>%
  group_by(river) %>%
  summarise(
    med_swot_wd = round(median(width), 1),
    med_ortho_wd = round(median(ortho_width_m), 1),
    # error metrics
    error_abs_68ile = round(quantile(abs(residuals), 0.68, na.rm = TRUE), 1),
    error_abs_50ile = round(quantile(abs(residuals), 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias), na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE), 1),
    error_perdiff_68ile = round(quantile(abs(percent_diff), 0.68, na.rm = TRUE), 2),
    error_perdiff_50ile = round(quantile(percent_diff, 0.50, na.rm = TRUE), 2),
    # min_error_percentdiff_68ile = round(min(abs(percent_diff)), 1),
    # max_error_percentdiff_50ile = round(max(abs(percent_diff)), 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_ortho_vD %>%
  group_by(river) %>%
  summarise(
    r_value = round(cor(width, ortho_width_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(width, ortho_width_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_node_width <- table_absolute_node_width %>%
  left_join(cor_table, by = "river")




# RELATIVE NODE WSE TABLE BY VERSION INCLUSION
# -----------------------------------------------------

# NODES UNIQUE TO vC & vD
table_absolute_node_width <- node_SWOT_ortho %>%
  group_by(version_inclusion) %>%
  summarise(
    # error metrics
    error_abs_68ile = round(quantile(abs(residuals), 0.68, na.rm = TRUE), 1),
    error_abs_50ile = round(quantile(abs(residuals), 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias), na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE), 1),
    error_perdiff_68ile = round(quantile(abs(percent_diff), 0.68, na.rm = TRUE), 2),
    error_perdiff_50ile = round(quantile(percent_diff, 0.50, na.rm = TRUE), 2),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_ortho %>%
  group_by(version_inclusion) %>%
  summarise(
    r_value = round(cor(width, ortho_width_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(width, ortho_width_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_node_width <- table_absolute_node_width %>%
  left_join(cor_table, by = "version_inclusion") %>%
  filter(version_inclusion != 0) %>% # drop 0, which are obs in both C&D
  mutate(version_inclusion = factor(version_inclusion, labels = c("vC0", "vD0")))


# !!!!!!!!!! need to fix this chunk once I have vD from riverSP so I can sort based on pass/cycle id too

# # SAME SUBSET OF NODES FOR vC & vD
same_version_subset_node_SWOT_insitu <- node_SWOT_ortho %>%
  filter(version_inclusion == 0) %>%                            # keep only reaches present in both versions
  group_by(node_id, cycle_id, pass_id) %>%
  filter(all(c("PIC0", "PGD0") %in% source)) %>%        # require both sources initially
  mutate(
    RiverSP_resid_na   = any(source == "PIC0"   & is.na(residuals_nobias)),
    RiverTile_resid_na = any(source == "PGD0" & is.na(residuals_nobias))) %>%
  # drop the partner row when the counterpart has NA residuals_nobias
  filter(
    !(source == "PGD0" & RiverSP_resid_na),
    !(source == "PIC0"   & RiverTile_resid_na)) %>%
  # after removals, keep only triples that still contain both sources
  filter(all(c("PIC0", "PGD0") %in% source)) %>%
  ungroup() %>%
  select(-RiverSP_resid_na, -RiverTile_resid_na)

table_absolute_node_width <- same_version_subset_node_SWOT_insitu %>%
  filter(version_inclusion == 0) %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_abs_68ile = round(quantile(abs(residuals), 0.68, na.rm = TRUE), 1),
    error_abs_50ile = round(quantile(abs(residuals), 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias), na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE), 1),
    error_perdiff_68ile = round(quantile(abs(percent_diff), 0.68, na.rm = TRUE), 2),
    error_perdiff_50ile = round(quantile(percent_diff, 0.50, na.rm = TRUE), 2),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_ortho %>%
  filter(version_inclusion == 0) %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(width, ortho_width_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(width, ortho_width_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_node_width <- table_absolute_node_width %>%
  left_join(cor_table, by = "source")













# ---------------------------------------------------------------------------------------------------------------------------
# PLOTS -- SUMMARY STATS
# ---------------------------------------------------------------------------------------------------------------------------

# RELATIVE C vs D
# --------------------------------------------------

# Compute n
n_relative_df <- node_SWOT_ortho %>%
  group_by(source) %>%
  summarise(n_unique_nodes = n_distinct(node_id), # count of non-NA residuals
            n = sum(!is.na(residuals)), .groups = "drop") # count of unique nodes

# CDF plot By Absolute Difference
ggplot(node_SWOT_ortho, aes(x = abs(residuals), color = source, linetype = source)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = expression("| SWOT -" ~ italic("in situ") ~ "Width | (m)"), y = "Cumulative Probability",
       title = "By absolute difference") +
  annotate("text", x = 240, y = 0.71, hjust = 0,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PIC0", ]$residuals), 0.68, na.rm = TRUE), 1),
                         "m, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PGD0", ]$residuals), 0.68, na.rm = TRUE), 1),
                         "m"),
           color = "#222222", size = 5) +
  annotate("text", x = 240, y = 0.53, hjust = 0,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PIC0", ]$residuals), 0.5, na.rm = TRUE), 1),
                         "m, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PGD0", ]$residuals), 0.5, na.rm = TRUE), 1),
                         "m"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1, vjust = 0,
           label = paste0("Version C: ", 
                          n_relative_df[n_relative_df$source == "PIC0", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "PIC0", ]$n, " total"),
           color = "#E69F00", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1, vjust = 0, 
           label = paste0("Version D: ", 
                          n_relative_df[n_relative_df$source == "PGD0", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "PGD0", ]$n, " total"), 
           color = "#0072B2", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("PIC0" = "#E69F00", "PGD0" = "#0072B2")) +
  scale_linetype_manual(values = c("PIC0" = "solid", "PGD0" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 500))
# width 7.17 height 6.35



# CDF plot By Percent Difference
ggplot(node_SWOT_ortho, aes(x = abs(percent_diff), color = source, linetype = source)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = expression("| SWOT -" ~ italic("in situ") ~ "Width | (% difference)"), y = "Cumulative Probability", 
       title = "By percent difference") +
  annotate("text", x = 70, y = 0.71, hjust = 0,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PIC0", ]$percent_diff), 0.68, na.rm = TRUE), 1),
                         "%, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PGD0", ]$percent_diff), 0.68, na.rm = TRUE), 1),
                         "%"),
           color = "#222222", size = 5) +
  annotate("text", x = 70, y = 0.53, hjust = 0,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PIC0", ]$percent_diff), 0.5, na.rm = TRUE), 1),
                         "%, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "PGD0", ]$percent_diff), 0.5, na.rm = TRUE), 1),
                         "%"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1, vjust = 0,
           label = paste0("Version C: ", 
                          n_relative_df[n_relative_df$source == "PIC0", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "PIC0", ]$n, " total"),
           color = "#E69F00", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1, vjust = 0, 
           label = paste0("Version D: ", 
                          n_relative_df[n_relative_df$source == "PGD0", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "PGD0", ]$n, " total"), 
           color = "#0072B2", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("PIC0" = "#E69F00", "PGD0" = "#0072B2")) +
  scale_linetype_manual(values = c("PIC0" = "solid", "PGD0" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# width 7.17 height 6.35




# ---------------------------------------------------------------------------------------------------------------------------
# INTER-RIVER COMPARISON PLOTS
# ---------------------------------------------------------------------------------------------------------------------------

# Reorder the factor levels by river to line up with color palette
node_SWOT_ortho_vD$river <- factor(
  node_SWOT_ortho_vD$river,
  levels = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"))

# compute n for each river
counts <- node_SWOT_ortho_vD %>%
  group_by(river) %>%
  summarise(n = n()) %>%
  ungroup()

# now make sure facor level of counts match with color palette
counts$river <- factor(counts$river, levels = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"))

color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")

# violin! percent difference
ggplot(node_SWOT_ortho_vD, aes(x = river, y = abs(percent_diff), fill = river)) + 
  geom_violin(alpha = 0.8, color = NA) +
  xlab("River") +
  ylab(expression(atop("SWOT -" ~ italic("in situ") ~ "Width", "(% difference)"))) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  # add counts below each violin
  geom_text(data = counts,
            aes(x = river, y = -2, label = paste0("n=", n)),
            inherit.aes = FALSE,
            vjust = 1, size = 6) +
  theme_minimal(base_size = 25) +
  scale_fill_manual(
    values = color_palette,
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  scale_x_discrete(
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 0.9),
        # give a little extra bottom margin so the -1 labels aren't cut off
        plot.margin = margin(t = 5, r = 5, b = 20, l = 5)) +
  coord_cartesian(ylim = c(-5, 150))
# 9.44, 6.01


# violin absolute
ggplot(node_SWOT_ortho_vD, aes(x = river, y = abs(residuals), fill = river)) + 
  geom_violin(alpha = 0.8, color = NA) +
  xlab("River") +
  ylab(expression(atop("| SWOT -" ~ italic("in situ") ~ "Width |(m)"))) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  # add counts below each violin
  geom_text(data = counts,
            aes(x = river, y = -2, label = paste0("n=", n)),
            inherit.aes = FALSE,
            vjust = 1, size = 6) +
  theme_minimal(base_size = 25) +
  scale_fill_manual(
    values = color_palette,
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  scale_x_discrete(
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 0.9),
        # give a little extra bottom margin so the -1 labels aren't cut off
        plot.margin = margin(t = 5, r = 5, b = 20, l = 5)) +
  coord_cartesian(ylim = c(-5, 750))
# 9.44, 6.01


# SCATTER PLOT

# Reorder the factor levels by river to line up with color palette
node_SWOT_ortho_vD$river <- factor(
  node_SWOT_ortho_vD$river,
  levels = c("CL", "CD", "PR", "upperYR", "lowerYR", "SJ"))

color_palette <- c("#F2C14E","#3B6064", "#F4845F", "#DA627D", "#9A348E", "#8EAD7A")
# correlation test
cor_test <- cor.test(node_SWOT_ortho_vD$width, node_SWOT_ortho_vD$ortho_width_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# scatterplot SWOT vs ortho width
ggplot() +
  geom_point(data = subset(node_SWOT_ortho_vD, river != "SJ"), aes(x = ortho_width_m, y = width, color = river), size = 2.5) +
  # plot SJ last so we can see the points
  geom_point(data = subset(node_SWOT_ortho_vD, river == "SJ"), aes(x = ortho_width_m, y = width, color = river), size = 2.5) +
  scale_color_manual(values = color_palette) +
  xlab(expression(atop(~ italic("In situ") ~ "Width (m)"))) +
  ylab("SWOT Width (m)") +
  theme_minimal(base_size = 25) +
  geom_abline(linetype = "dashed", color = "gray") +
  ylim(0,2700) +
  xlim(0,2700) +
  annotate("text",
           x = min(node_SWOT_ortho_vD$ortho_width_m, na.rm = TRUE),
           y = 2700,
           label = paste0("r = ", round(r_value, 4),
                          "\np value = ", round(signif(p_value, 3), 4),
                          "\nn = ", nrow(node_SWOT_ortho_vD)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none")
# 6.66, 6.01






#-------------------------------------------------------------------------------------
# Cross track bias explore

# Filter residuals to 2.5–97.5% per river
node_SWOT_ortho_vD_filt <- node_SWOT_ortho_vD %>%
  filter(!is.na(residuals), !is.na(xtrk_dist)) %>%
  group_by(river) %>%
  filter(
    residuals >= quantile(residuals, 0.025, na.rm = TRUE),
    residuals <= quantile(residuals, 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

# Correlation stats by river after filtering
cor_stats <- node_SWOT_ortho_vD_filt %>%
  group_by(river) %>%
  summarise(
    n = n(),
    r_value = ifelse(n >= 3, cor(residuals, abs(xtrk_dist), method = "pearson"), NA_real_),
    p_value = ifelse(n >= 3, cor.test(residuals, abs(xtrk_dist))$p.value, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "r = ", round(r_value, 4),
      "\np value = ", signif(p_value, 3),
      "\nn = ", n
    )
  )

# Plot
ggplot(
  node_SWOT_ortho_vD_filt,
  aes(x = abs(xtrk_dist / 1000), y = residuals, color = river)
) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  facet_wrap(~ river, scales = "free") +
  scale_color_manual(values = color_palette) +
  xlab(expression(atop("|Cross-track distance| (km)"))) +
  ylab("SWOT Width residuals (m)") +
  geom_text(
    data = cor_stats,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 6
  ) +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none")






# Overall correlation stats across all rivers
cor_stats_all <- node_SWOT_ortho_vD_filt %>%
  summarise(
    n = n(),
    r_value = ifelse(n >= 3, cor(residuals, abs(xtrk_dist), method = "pearson"), NA_real_),
    p_value = ifelse(n >= 3, cor.test(residuals, abs(xtrk_dist))$p.value, NA_real_)
  ) %>%
  mutate(
    label = paste0(
      "All rivers\n",
      "r = ", round(r_value, 4),
      "\np = ", signif(p_value, 3),
      "\nn = ", n
    )
  )

ggplot(
  node_SWOT_ortho_vD_filt,
  aes(x = abs(xtrk_dist / 1000), y = residuals, color = river)
) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = FALSE,
    linewidth = 1,
    color = "black"
  ) +
  scale_color_manual(values = color_palette) +
  xlab("|Cross-track distance| (km)") +
  ylab("SWOT Width residuals (m)") +
  geom_text(
    data = cor_stats_all,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 6
  ) +
  theme_minimal(base_size = 25)














# Expand to long df for width type comparisons
long_df <- node_SWOT_ortho_vD %>%
  select(river, ortho_width_m, width) %>%
  pivot_longer(cols = c(ortho_width_m, width),
               names_to = "width_type",
               values_to = "width_value") %>%
  mutate(width_type = recode(width_type,
                             ortho_width_m = "Ortho",
                             width = "SWOT"))

# Set river order
long_df$river <- factor(long_df$river, 
                        levels = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"))

# split violin comparisons
devtools::install_github("psyteachr/introdataviz")
ggplot(long_df, aes(x = river, y = width_value, fill = width_type)) +
  introdataviz::geom_split_violin(alpha = .5, trim = FALSE, color= NA) +
  geom_boxplot(width = .2, alpha = .8, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  scale_x_discrete(labels = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")) +
  scale_fill_manual(values=c("lightblue","darkblue")) +
  theme_minimal(base_size = 30) +
  ylab("width (m)") +
  coord_cartesian(ylim = c(0, 3500))
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 0.9))
  
  


  
  
  library(dplyr)
  library(ggplot2)
  library(introdataviz)
  
  rivers <- c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
  
  plots <- lapply(rivers, function(riv) {
    
    df_riv <- long_df %>% filter(river == riv)
    
    ggplot(df_riv, aes(x = river, y = width_value, fill = width_type)) +
      introdataviz::geom_split_violin(
        alpha = 0.5, trim = FALSE, color = NA
      ) +
      geom_boxplot(
        width = 0.2, alpha = 0.8, fatten = NULL,
        show.legend = FALSE
      ) +
      stat_summary(
        fun.data = "mean_se",
        geom = "pointrange",
        position = position_dodge(0.175),
        show.legend = FALSE
      ) +
      scale_fill_manual(values = c("lightblue", "darkblue")) +
      theme_minimal(base_size = 30) +
      labs(
        title = riv,
        y = "Width (m)",
        x = NULL
      ) +
      coord_cartesian(ylim = c(0, 300)) +
      theme(
        legend.position = "none"
      )
  })
  
  names(plots) <- rivers
  plots$CL
  


# save node means for data viz

node_means <- node_SWOT_ortho_vD %>%
  group_by(node_id) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))


write.csv(node_means, "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverTile_v17b/node_avg_width_SWOT_Ortho.csv", row.names = FALSE)























stats <- node_SWOT_ortho_vD %>%
  group_by(river) %>%
  summarise(
    n = n(),
    q68 = quantile(abs(percent_diff), 0.68, na.rm = TRUE)
  ) %>%
  ungroup()

stats$river <- factor(stats$river,
                      levels = c("CL","SJ","CD","PR","upperYR","lowerYR"))

ggplot(node_SWOT_ortho_vD, aes(x = river, y = abs(percent_diff), fill = river)) +
  geom_violin(alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  # combined multi-line label
  geom_text(data = stats,
            aes(x = river, y = -3,
                label = paste0("q=", round(q68,1), "%\n", "n=", n)),
            inherit.aes = FALSE,
            vjust = 1,
            size = 5,
            lineheight = 0.9) +
  coord_cartesian(ylim = c(-8, 150)) +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none", panel.clip = "off")

















# ---------------------------------------------------------------------------------------------------------------------------
# OLD SCRATCH
# ---------------------------------------------------------------------------------------------------------------------------




# AGU MESS
# --------------------------------------------------

# Compute n
n_relative_df <- node_SWOT_ortho %>%
  group_by(source) %>%
  summarise(n_unique_nodes = n_distinct(node_id), # count of non-NA residuals
            n = sum(!is.na(residuals)), .groups = "drop") # count of unique nodes

# CDF plot By Absolute Difference
ggplot(node_SWOT_ortho, aes(x = abs(residuals), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT -" ~ italic("in situ") ~  "Width (m)", y = "Cumulative Probability", 
       title = "Node Absolute Difference") +
  annotate("text", x = 345, y = 0.71,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverSP", ]$residuals), 0.68, na.rm = TRUE), 1),
                         "m, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverTile", ]$residuals), 0.68, na.rm = TRUE), 1),
                         "m"),
           color = "#222222", size = 6) +
  annotate("text", x = 345, y = 0.53,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverSP", ]$residuals), 0.5, na.rm = TRUE), 1),
                         "m, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverTile", ]$residuals), 0.5, na.rm = TRUE), 1),
                         "m"),
           color = "#222222", size = 6) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1, vjust = 0,
           label = paste0("Version C: ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n, " total"),
           color = "#E97132", size = 6) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1, vjust = 0, 
           label = paste0("Version D: ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n, " total"), 
           color = "darkblue", size = 6) +
  theme_minimal(base_size = 20) +
  scale_color_manual(values = c("RiverSP" = "#E97132", "RiverTile" = "darkblue")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 500))
# width 610 height 550



# CDF plot By Percent Difference
ggplot(node_SWOT_ortho, aes(x = abs(percent_diff), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT -" ~ italic("in situ") ~ "Width (% difference)", y = "Cumulative Probability", 
       title = "Node Absolute Width Percent Difference") +
  annotate("text", x = 100, y = 0.71,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverSP", ]$percent_diff), 0.68, na.rm = TRUE), 1),
                         "%, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverTile", ]$percent_diff), 0.68, na.rm = TRUE), 1),
                         "%"),
           color = "#222222", size = 6) +
  annotate("text", x = 100, y = 0.53,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverSP", ]$percent_diff), 0.5, na.rm = TRUE), 1),
                         "%, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverTile", ]$percent_diff), 0.5, na.rm = TRUE), 1),
                         "%"),
           color = "#222222", size = 6) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1.1, vjust = 0,
           label = paste0("Version C: ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n, " total"),
           color = "#E97132", size = 6) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1.1, vjust = 0, 
           label = paste0("Version D: ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n, " total"), 
           color = "darkblue", size = 6) +
  theme_minimal(base_size = 20) +
  scale_color_manual(values = c("RiverSP" = "#E97132", "RiverTile" = "darkblue")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# width 610 height 550










AOI <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/SWOT/YF24_AOI_SWOT_tiles.csv')


length(unique(AOI$pass))


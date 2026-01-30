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
  mutate(source = "RiverSP") %>%
  filter(abs(residuals) < 1500) %>%
  filter(dark_frac < 0.5)

node_SWOT_ortho_vD <- read_csv('/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverTile_v17b/node_width_SWOT_Ortho.csv') %>%
  mutate(source = "RiverTile") %>%
  filter(abs(residuals) < 1500) %>%
  filter(dark_frac < 0.5)

# merge all dataframes together
node_SWOT_ortho <- bind_rows(node_SWOT_ortho_vC, node_SWOT_ortho_vD) %>%
  filter(dark_frac < 0.5)


# ---------------------------------------------------------------------------------------------------------------------------
# TABLES -- SUMMARY STATS
# ---------------------------------------------------------------------------------------------------------------------------

# ABSOLUTE NODE WIDTH TABLE BY VERSION
# -----------------------------------------------------
table_absolute_node_width <- node_SWOT_ortho %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_absolute_68ile = quantile(abs(residuals), 0.68, na.rm = TRUE),
    error_absolute_50ile = quantile(abs(residuals), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals_nobias), na.rm = TRUE),
    error_percentdiff_68ile = quantile(percent_diff, 0.68, na.rm = TRUE),
    error_percentdiff_50ile = quantile(percent_diff, 0.50, na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_ortho %>%
  group_by(source) %>%
  summarise(
    r_value = cor(width, ortho_width_m, use = "complete.obs", method = "pearson"),
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
    median_swot_width = median(width),
    median_ortho_width = median(ortho_width_m),
    # error metrics
    error_absolute_68ile = quantile(abs(residuals), 0.68, na.rm = TRUE),
    error_absolute_50ile = quantile(abs(residuals), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals_nobias), na.rm = TRUE),
    error_percentdiff_68ile = quantile(abs(percent_diff), 0.68, na.rm = TRUE),
    error_percentdiff_50ile = quantile(percent_diff, 0.50, na.rm = TRUE),
    min_error_percentdiff_68ile = min(abs(percent_diff)),
    max_error_percentdiff_50ile = max(abs(percent_diff)),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    
    # !!!!!!!!!!!!!! how to quantify count when this will vary from vC to vD
    # right now I'm taking RiverTile which has more unique obs
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id[source == "RiverTile"]))

# Add correlations
cor_table <- node_SWOT_ortho_vD %>%
  group_by(river) %>%
  summarise(
    r_value = cor(width, ortho_width_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(width, ortho_width_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_node_width <- table_absolute_node_width %>%
  left_join(cor_table, by = "river")


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
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = expression("| SWOT -" ~ italic("in situ") ~ "Width | (m)"), y = "Cumulative Probability",
       title = "By absolute difference") +
  annotate("text", x = 240, y = 0.71, hjust = 0,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverSP", ]$residuals), 0.68, na.rm = TRUE), 1),
                         "m, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverTile", ]$residuals), 0.68, na.rm = TRUE), 1),
                         "m"),
           color = "#222222", size = 5) +
  annotate("text", x = 240, y = 0.53, hjust = 0,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverSP", ]$residuals), 0.5, na.rm = TRUE), 1),
                         "m, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverTile", ]$residuals), 0.5, na.rm = TRUE), 1),
                         "m"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1, vjust = 0,
           label = paste0("Version C: ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n, " total"),
           color = "#E97132", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1, vjust = 0, 
           label = paste0("Version D: ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n, " total"), 
           color = "darkblue", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "#E97132", "RiverTile" = "darkblue")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 500))
# width 7.17 height 6.35



# CDF plot By Percent Difference
ggplot(node_SWOT_ortho, aes(x = abs(percent_diff), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = expression("| SWOT -" ~ italic("in situ") ~ "Width | (% difference)"), y = "Cumulative Probability", 
       title = "By percent difference") +
  annotate("text", x = 70, y = 0.71, hjust = 0,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverSP", ]$percent_diff), 0.68, na.rm = TRUE), 1),
                         "%, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverTile", ]$percent_diff), 0.68, na.rm = TRUE), 1),
                         "%"),
           color = "#222222", size = 5) +
  annotate("text", x = 70, y = 0.53, hjust = 0,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverSP", ]$percent_diff), 0.5, na.rm = TRUE), 1),
                         "%, vD:", 
                         round(quantile(abs(node_SWOT_ortho[node_SWOT_ortho$source == "RiverTile", ]$percent_diff), 0.5, na.rm = TRUE), 1),
                         "%"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1, vjust = 0,
           label = paste0("Version C: ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n, " total"),
           color = "#E97132", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1, vjust = 0, 
           label = paste0("Version D: ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n, " total"), 
           color = "darkblue", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "#E97132", "RiverTile" = "darkblue")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# width 610 height 550




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
  ylab("|SWOT - Ortho Width| (m)") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  theme_minimal(base_size = 25) +
  scale_fill_manual(
    values = color_palette,
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  scale_x_discrete(
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 0.9)) + 
  coord_cartesian(ylim = c(0, 750))



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








# all clusters comparison
# ---------------------------------------------------------------------------------------------------------------------------

# Set working directory to the orthomosaics directory with the combined_ortho_SWOT_ csvs
# wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16"
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverTile_v17b"


setwd(wd)

# Get list of all CSV files in working directory that start with combined_ortho_SWOT
csv_files <- list.files(wd, pattern = "^combined_ortho_SWOT_.*\\.csv$", full.names = TRUE)

# Merge all PT files into a combined dataframe with filename column
data_list <- lapply(csv_files, function(file) {
  df <- read.csv(file)
  return(df)
})

# Combine all dataframes into one
combined_df <- bind_rows(data_list)


# add river names to df
combined_df <- combined_df %>%
  mutate(river_code = substr(reach_id, 1, 6),
         river = case_when(
           # putting the reach id first ensures case_when won't overwrite SJ/BL labels
           # SWORD v16: "81260300061", "81260300231", "81260300241", "81260300251"
           # SWORD v17b: "81260300181", "81260300191", "81260300201", "81260300211"
           # only SJ reaches need to be adjusted here
           reach_id %in% c("81260300181", "81260300191", "81260300201", "81260300211") ~ "SJ", 
           reach_id %in% c("81270100111", "81270100121", "81270100131", "81270100141", "81270100151", "81270100161", "81270200011", "81270200021") ~ "BL",
           river_code == "812701" ~ "lowerYR", # until the Circle bifurcation
           river_code == "812509" ~ "lowerYR", # past the PR confluence
           river_code == "812705" ~ "upperYR", # Circle up
           river_code == "812508" ~ "CD",
           river_code == "812603" ~ "PR",
           river_code == "812605" ~ "PR",
           river_code == "812604" ~ "CL",
           TRUE ~ NA_character_))


# ---------------------------------------------------------------------------------------------------------------------------
# data viz for RiverSP / RiverTile
# ---------------------------------------------------------------------------------------------------------------------------

# "CD", "CL", "lowerPR", "lowerYR", "SJ", "upperPR", "upperYR"
color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#F4845F", "#DA627D")


# Combo CDF plot
ggplot(combined_ortho_SWORD_df, aes(x = percent_diff, color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability", 
       title = "CDF of SWOT - Ortho Width %diff") +
  annotate("text", x = 125, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error, 2), 
                         ", Version D:", round(percentile_68_error_RiverTile, 2)), 
           color = "#222222", size = 5) +
  annotate("text", x = 125, y = 0.53, 
           label = paste("|50%ile| Version C:", round(percentile_50_error, 2), 
                         ", Version D:", round(percentile_50_error_RiverTile, 2)), 
           color = "#222222", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "darkblue", "RiverTile" = "#E97132")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "longdash")) +
  xlim(0, 200)

# Combo CDF plot
ggplot(combined_ortho_SWORD_df, aes(x = abs(residuals), color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - Ortho Width (m)", y = "Cumulative Probability", 
       title = "CDF of SWOT - Ortho Width (m)") +
  annotate("text", x = 525, y = 0.71, 
           label = paste("|68%ile| Version C:", round(percentile_68_error, 2), 
                         ", Version D:", round(percentile_68_error_RiverTile, 2)), 
           color = "#222222", size = 5) +
  annotate("text", x = 525, y = 0.53, 
           label = paste("|50%ile| Version C:", round(percentile_50_error, 2), 
                         ", Version D:", round(percentile_50_error_RiverTile, 2)), 
           color = "#222222", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "darkblue", "RiverTile" = "#E97132")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "longdash")) +
  xlim(0, 1000)

# ---------------------------------------------------------------------------------------------------------------------------
# data viz for one SWOT data type

# c("CL" = "#F2C14E", "upper_PR" = "#F4845F", "upper_YR" ="#DA627D", lower_YR = "#9A348E"))
color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#F4845F", "#DA627D")

# color_palette <- c("#F2C14E", "#F4845F", "#DA627D", "#9A348E", )


# plot SWOT vs GNSS width
ggplot(combined_ortho_SWORD_df, aes(x = ortho_width_m, y = width, color = factor(river))) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_palette) +
  xlab("Ortho width (m)") +
  ylab("SWOT width (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", x = min(combined_ortho_SWORD_df$ortho_width_m, na.rm = TRUE), 
           y = max(combined_ortho_SWORD_df$width, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", round(signif(p_value, 3), 4),
                          "\nn = ", nrow(combined_ortho_SWORD_df)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none") + 
  labs(color = "Reach ID") 


# CDF plot of absolute wse difference
ggplot(combined_ortho_SWORD_df, aes(x = abs(width - ortho_width_m))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - Ortho Width (m)", y = "Cumulative Probability", title = "CDF of SWOT Width - Ortho Width") +
  annotate("text", x = 350, y = 0.71, label = paste("68% abs diff:", round(percentile_68_error, 4)), color = "#222222", size = 6) +
  annotate("text", x = 350, y = 0.53, label = paste("50% abs diff:", round(percentile_50_error, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20) +
  xlim(0, 1500)

# CDF plot of percent difference
ggplot(combined_ortho_SWORD_df, aes(x = abs(percent_diff))) +
  stat_ecdf(geom = "step", color = "darkblue", size = 1) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability", title = "CDF of SWOT-Ortho Width %diff") +
  annotate("text", x = 75, y = 0.71, label = paste("68% abs diff:", round(percentile_68_percent, 4)), color = "#222222", size = 6) +
  annotate("text", x = 75, y = 0.53, label = paste("50% abs diff:", round(percentile_50_percent, 4)), color = "#222222", size = 6) +
  theme_minimal(base_size = 20) +
  xlim(0, 100)



# ---------------------------------------------------------------------------------------------------------------------------
# Comparison with SWOT parameters plots

# correlation test
cor_test <- cor.test(abs(combined_ortho_SWORD_df$dark_frac), combined_ortho_SWORD_df$percent_diff)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001

# plot percent diff vs swot df variables
ggplot(combined_ortho_SWORD_df, aes(x = dark_frac, y = percent_diff, color = river)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_palette) +
  xlab("dark_frac") +
  ylab("percent_diff (%)") +
  theme_minimal(base_size = 30) +
  annotate("text", x = max(combined_ortho_SWORD_df$dark_frac-.28, na.rm = TRUE), 
           y = max(combined_ortho_SWORD_df$percent_diff, na.rm = TRUE), 
           label = paste0("r = ", round(r_value, 4), "\np value = ", round(signif(p_value, 3),5)),
           hjust = 0, vjust = 1, size = 8) +
  theme(legend.position = "none")

color_palette <- c("#D86A1A", "#6D398B", "#00429D", "#F8A31B", "#2E7D32")
color_palette <- c("#6D398B", "#00429D","#2E7D32")

# violin comparison
# factor(node_q_b), factor(pass_id)
ggplot(combined_ortho_SWORD_df, aes(x = factor(substr(time_utc, 0, 10)), y = percent_diff_nobias, fill=factor(substr(time_utc, 0, 10))), drop = FALSE) + 
  geom_violin(alpha = 0.8, color= NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  xlab("time_utc") +
  ylab("Percent Difference (%)") +
  scale_fill_manual(values=c("lightblue","darkblue")) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none")
# violin comparision by river
ggplot(combined_ortho_SWORD_df, aes(x = river, y = percent_diff, fill=river), drop = FALSE) + 
  geom_violin(alpha = 0.8, color= NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  xlab("River") +
  ylab("Percent Difference (%)") +
  scale_fill_manual(values=color_palette) +
  theme_minimal(base_size = 30)+
  scale_x_discrete(
    breaks = c("CL", "upper_PR", "upper_YR"),
    labels = c("Coleen", "Upper Porcupine", "Upper Yukon")
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 0.9)) +
  ylim(0, 270)

# "CD", "CL", "lowerPR", "lowerYR", "SJ", "upperPR", "upperYR"
color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#F4845F", "#DA627D")


# Combo CDF plot by river
ggplot(node_SWOT_ortho_vC, aes(x = percent_diff, color = river, linetype = river)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability", 
       title = "CDF of SWOT - Ortho width %diff") +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("CD" = "#3B6064", "CL" = "#F2C14E", "upper_PR" = "#F4845F", "upper_YR" ="#DA627D", "lower_YR" = "#9A348E", "SJ" = "#8EAD7A")) +
  #theme(legend.position = "none") +
  xlim(0, 100)

# Combo CDF plot by river
ggplot(RiverTile_df, aes(x = abs(residuals), color = river, linetype = river)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - Ortho Width (m)", y = "Cumulative Probability", 
       title = "CDF of SWOT - Ortho width %diff") +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("CL" = "#F2C14E", "upper_PR" = "#F4845F", "upper_YR" ="#DA627D", lower_YR = "#9A348E")) +
  theme(legend.position = "none") 

# "CD", "CL", "lowerPR", "lowerYR", "SJ", "upperPR", "upperYR"
color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#F4845F", "#DA627D")

# Expand to long df for width type comparisons
long_df <- RiverTile_df %>%
  select(river, ortho_width_m, width) %>%
  pivot_longer(cols = c(ortho_width_m, width),
               names_to = "width_type",
               values_to = "width_value") %>%
  mutate(width_type = recode(width_type,
                             ortho_width_m = "Ortho",
                             width = "SWOT"))

# Set river order
long_df$river <- factor(long_df$river, 
                        levels = c("CL", "upper_PR", "upper_YR", "lower_YR"))

# split violin comparisons
devtools::install_github("psyteachr/introdataviz")
ggplot(long_df, aes(x = river, y = width_value, fill = width_type)) +
  introdataviz::geom_split_violin(alpha = .5, trim = FALSE, color= NA) +
  geom_boxplot(width = .2, alpha = .8, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  scale_x_discrete(labels = c("Coleen", 'Porcupine', "Single-channel Yukon", "Braided Yukon")) +
  scale_fill_manual(values=c("lightblue","darkblue")) +
  theme_minimal(base_size = 30) +
  ylab("width (m)") +
  ylim(0, 3500) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 0.9))

# Expand to long df for width type comparisons by day
long_df <- combined_ortho_SWORD_df %>%
  select(time_utc, ortho_width_m, width) %>%
  pivot_longer(cols = c(ortho_width_m, width),
               names_to = "width_type",
               values_to = "width_value") %>%
  mutate(width_type = recode(width_type,
                             ortho_width_m = "Ortho",
                             width = "SWOT"))
# split violin comparisons
devtools::install_github("psyteachr/introdataviz")
ggplot(long_df, aes(x = factor(substr(time_utc, 0, 10)), y = width_value, fill = width_type)) +
  introdataviz::geom_split_violin(alpha = .5, trim = FALSE, color= NA) +
  geom_boxplot(width = .2, alpha = .8, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  scale_fill_manual(values=c("lightblue","darkblue")) +
  theme_minimal(base_size = 30) +
  ylab("width (m)") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 0.9))

# CDF plot by day (used to look at upper PR differences)
ggplot(combined_ortho_SWORD_df, aes(x = percent_diff, color = factor(substr(time_utc, 0, 10)), linetype = factor(substr(time_utc, 0, 10)))) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "Width Percent Difference (%)", y = "Cumulative Probability", 
       title = "CDF of SWOT - Ortho width upper_PR") +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("2024-07-10" = "lightblue", "2024-07-16" = "#00429D")) +
  theme(legend.position = "none")





















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


library(tidyverse)
library(lubridate)
library(dplyr)
library(ggtext)

# ---------------------------------------------------------------------------------------------------------------------------
# SWOT Slope validation
# ---------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------
# read in data
# ---------------------------------------------------------------------------------------------------------------------------

# PT
reach_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_PT.csv") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_slope_SWOT_PT.csv") %>%
  filter(dark_frac < 0.5)

# GNSS
reach_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_GNSS.csv") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_slope_SWOT_GNSS.csv") %>%
  filter(dark_frac < 0.5)

# ---------------------------------------------------------------------------------------------------------------------------
# Get all data to the same SWORD version (v17b)
# ---------------------------------------------------------------------------------------------------------------------------

# SWORD translator to change Version C data to SWORD v17b naming convention
SWORD_translator <- read_csv("/Users/camryn/Desktop/SWORD_translation/NA_ReachIDs_v17b_vs_v16.csv")

# translate the vC SWORD v16 data to SWORD v17b
reach_SWOT_PT_vC <- reach_SWOT_PT_vC %>%
  left_join(SWORD_translator %>% 
              select(v16_reach_id, v17_reach_id),
            by = c("old_reach_id" = "v16_reach_id")) %>%
  rename(reach_id = v17_reach_id)

reach_SWOT_GNSS_vC <- reach_SWOT_GNSS_vC %>%
  left_join(SWORD_translator %>% 
              select(v16_reach_id, v17_reach_id),
            by = c("old_reach_id" = "v16_reach_id")) %>%
  rename(reach_id = v17_reach_id)

# merge all dataframes together
reach_SWOT_full_insitu <- bind_rows(reach_SWOT_PT_vC, reach_SWOT_PT_vD, reach_SWOT_GNSS_vC, reach_SWOT_GNSS_vD) %>%
  mutate(insitu_slope_m_m = coalesce(slope_m_m_abs, reach_drift_slope_m_m_abs)) %>%
  mutate(insitu_slope_nobias_m_m = coalesce(mean_reach_PT_slope_no_bias_m_m, reach_drift_slope_m_m_abs_nobias)) %>%
  mutate(insitu_time_utc = coalesce(pt_time_UTC, wse_drift_midpoint_UTC))

# compute version inclusion
# -1 = only in vC, 0 = both, 1 = only in vD
all_reaches <- reach_SWOT_full_insitu %>%
  distinct(reach_id, source, insitu_time_utc) %>%         
  group_by(reach_id) %>%
  summarise(
    has_RiverSP   = any(source == "RiverSP"),
    has_RiverTile = any(source == "RiverTile"),
    .groups = "drop") %>%
  mutate(
    version_inclusion = case_when(has_RiverSP & has_RiverTile ~ 0L, has_RiverSP & !has_RiverTile ~ -1L, !has_RiverSP & has_RiverTile ~ 1L, TRUE ~ NA_integer_)) %>%
  select(reach_id, version_inclusion)

reach_SWOT_full_insitu <- reach_SWOT_full_insitu %>%
  left_join(all_reaches, by = "reach_id")

# ---------------------------------------------------------------------------------------------------------------------------
# TABLES -- SUMMARY STATS
# ---------------------------------------------------------------------------------------------------------------------------

# RELATIVE REACH WSE TABLE BY VERSION
# -----------------------------------------------------
table_relative_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_68ile = quantile(abs(slope_residuals_nobias)*100000, 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(slope_residuals_nobias)*100000, 0.50, na.rm = TRUE),
    MAE = mean(abs(slope_residuals_nobias)*100000, na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals_nobias)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(slope_abs, insitu_slope_nobias_m_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = "source")

# relabel and reorder
table_relative_reach_slope <- table_relative_reach_slope %>%
  mutate(source = factor(source,
                         levels = c("RiverTile", "RiverSP"),   # swapped order
                         labels = c("vD0", "vC0")))            # relabels

# bar chart of count of residuals_nobias by version
ggplot(table_relative_reach_slope, aes(x = source, y = n, fill = source)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n),
            vjust = -0.5,
            size = 8) +
  ylab("Count") +
  coord_cartesian(ylim = c(11, 150)) +
  scale_fill_manual(values = c("vC0" = "#E97132",
                               "vD0" = "darkblue")) +
  theme_classic(base_size = 34) +
  theme(axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")
# 3.16, 6.54

# RELATIVE REACH SLOPE TABLE BY GNSS/PT
# -----------------------------------------------------
table_relative_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    # error metrics (use all data)
    error_68ile = quantile(abs(slope_residuals_nobias)*100000, 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(slope_residuals_nobias)*100000, 0.50, na.rm = TRUE),
    MAE = mean(abs(slope_residuals_nobias), na.rm = TRUE)*100000,
    
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals_nobias)),
    
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id)
  )

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    r_value = cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(slope_abs, insitu_slope_nobias_m_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = "insitu_type")


# RELATIVE REACH WSE TABLE BY RIVER
# -----------------------------------------------------
table_relative_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(river) %>%
  summarise(
    # error metrics
    error_68ile = quantile(abs(slope_residuals_nobias)*100000, 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(slope_residuals_nobias)*100000, 0.50, na.rm = TRUE),
    MAE = mean(abs(slope_residuals_nobias)*100000, na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals_nobias)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(river) %>%
  summarise(
    r_value = cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(slope_abs, insitu_slope_nobias_m_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = "river")


# ABSOLUTE WSE TABLE BY VERSION
# -----------------------------------------------------
table_absolute_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_68ile = quantile(abs(slope_residuals)*100000, 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(slope_residuals)*100000, 0.50, na.rm = TRUE),
    MAE = mean(abs(slope_residuals)*100000, na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = cor(slope_abs, insitu_slope_m_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(slope_abs, insitu_slope_m_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_reach_slope <- table_absolute_reach_slope %>%
  left_join(cor_table, by = "source")

# ABSOLUTE WSE TABLE BY GNSS/PT
# -----------------------------------------------------
table_absolute_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    # error metrics
    error_68ile = quantile(abs(slope_residuals)*100000, 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(slope_residuals)*100000, 0.50, na.rm = TRUE),
    MAE = mean(abs(slope_residuals)*100000, na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    r_value = cor(slope_abs, insitu_slope_m_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(slope_abs, insitu_slope_m_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_reach_slope <- table_absolute_reach_slope %>%
  left_join(cor_table, by = "insitu_type")


# ---------------------------------------------------------------------------------------------------------------------------
# PLOTS -- SUMMARY STATS
# ---------------------------------------------------------------------------------------------------------------------------

# RELATIVE C vs D
# --------------------------------------------------

# Compute n
n_relative_df <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(n_unique_reaches = n_distinct(reach_id), # count of non-NA residuals
            n = sum(!is.na(slope_residuals_nobias)), .groups = "drop") # count of unique reaches

# CDF plot (fixed cm placement)
ggplot(reach_SWOT_full_insitu, aes(x = abs(slope_residuals_nobias)*100000, color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = expression("|SWOT -" ~ italic("in situ") ~ "Slope| (cm/km)"), y = "Cumulative Probability", 
       title = "By SWOT version") +
  annotate("text", x = 4, y = 0.71, hjust = 0,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverSP", ]$slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
                         "cm/km, vD:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverTile", ]$slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
                         "cm/km"),
           color = "#222222", size = 5) +
  annotate("text", x = 4, y = 0.53, hjust = 0,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverSP", ]$slope_residuals_nobias)*100000, 0.5, na.rm = TRUE), 2),
                         "cm/km, vD:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverTile", ]$slope_residuals_nobias)*100000, 0.5, na.rm = TRUE), 2),
                         "cm/km"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1, vjust = 0,
           label = paste0("Version C: ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n, " total"),
           color = "#E97132", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1, vjust = 0, 
           label = paste0("Version D: ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n, " total"), 
           color = "darkblue", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "#E97132", "RiverTile" = "darkblue")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 13))
# width 7.17 height 6.35



# RELATIVE GNSS / PT with lines for C/D split
# --------------------------------------------------

# Compute n
n_relative_df <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(n_unique_reaches = n_distinct(reach_id[source == "RiverTile"]), # count of non-NA residuals
            n = sum(!is.na(slope_residuals_nobias)), .groups = "drop") # count of unique reaches

# CDF plot (fixed cm placement)
ggplot(reach_SWOT_full_insitu, aes(x = abs(slope_residuals_nobias)*100000, color = insitu_type, linetype = insitu_type)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = expression("|SWOT -" ~ italic("in situ") ~ "Slope| (cm/km)"), y = "Cumulative Probability", 
       title = expression("By" ~ italic("in situ") ~ "measurement type")) +
  annotate("text", x = 4, y = 0.71, hjust = 0,
           label = paste("|68%ile| PT:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "PT", ]$slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
                         "cm/km, GNSS:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "GNSS", ]$slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
                         "cm/km"),
           color = "#222222", size = 5) +
  annotate("text", x = 4, y = 0.53, hjust = 0,
           label = paste("|50%ile| PT:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "PT", ]$slope_residuals_nobias)*100000, 0.5, na.rm = TRUE), 2),
                         "cm/km, GNSS:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "GNSS", ]$slope_residuals_nobias)*100000, 0.5, na.rm = TRUE), 2),
                         "cm/km"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1, vjust = 0,
           label = paste0("PT: ", 
                          n_relative_df[n_relative_df$insitu_type == "PT", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$insitu_type == "PT", ]$n, " total"),
           color = "#C03F61", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1, vjust = 0, 
           label = paste0("GNSS: ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n, " total"), 
           color = "#404A22", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("PT" = "#C03F61", "GNSS" = "#404A22")) +
  scale_linetype_manual(values = c("PT" = "solid", "GNSS" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 13))
# width 7.17 height 6.35



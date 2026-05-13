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
reach_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v17b/reach_slope_SWOT_PT.csv") %>%
  filter(dark_frac < 0.5)

# GNSS
reach_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_GNSS.csv") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v17b/reach_slope_SWOT_GNSS.csv") %>%
  filter(dark_frac < 0.5)

# ---------------------------------------------------------------------------------------------------------------------------
# Get all data to the same SWORD version (v17b)
# ---------------------------------------------------------------------------------------------------------------------------

# SWORD translator to change Version C data to SWORD v17b naming convention
SWORD_translator <- read_csv("/Users/camryn/Desktop/SWORD_translation/NA_ReachIDs_v17b_vs_v16.csv")

# translate the vC SWORD v16 data to SWORD v17b
reach_SWOT_PT_vC <- reach_SWOT_PT_vC %>%
  left_join(SWORD_translator %>% 
              dplyr::select(v16_reach_id, v17_reach_id),
            by = c("old_reach_id" = "v16_reach_id")) %>%
  rename(reach_id = v17_reach_id)

reach_SWOT_GNSS_vC <- reach_SWOT_GNSS_vC %>%
  left_join(SWORD_translator %>% 
              dplyr::select(v16_reach_id, v17_reach_id),
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
    has_RiverSP   = any(source == "PIC0"),
    has_RiverTile = any(source == "PGD0"),
    .groups = "drop") %>%
  mutate(
    version_inclusion = case_when(has_RiverSP & has_RiverTile ~ 0L, has_RiverSP & !has_RiverTile ~ -1L, !has_RiverSP & has_RiverTile ~ 1L, TRUE ~ NA_integer_)) %>%
  dplyr::select(reach_id, version_inclusion)

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
    error_68ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
    error_50ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.50, na.rm = TRUE), 2),
    MAE = round(mean(abs(slope_residuals_nobias)*100000, na.rm = TRUE), 2),
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals_nobias)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(slope_abs, insitu_slope_nobias_m_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = "source")

# relabel and reorder
table_relative_reach_slope <- table_relative_reach_slope %>%
  mutate(source = factor(source,
                         levels = c("PGD0", "PIC0"),   # swapped order
                         labels = c("vD0", "vC0")))            # relabels

# bar chart of count of residuals_nobias by version
ggplot(table_relative_reach_slope, aes(x = source, y = n, fill = source)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n),
            vjust = -0.5,
            size = 8) +
  ylab("Count") +
  coord_cartesian(ylim = c(11, 157)) +
  scale_fill_manual(values = c("vC0" = "#E69F00",
                               "vD0" = "#0072B2")) +
  theme_classic(base_size = 34) +
  theme(axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")
# 3.16, 6.54

# RELATIVE REACH SLOPE TABLE BY GNSS/PT
# -----------------------------------------------------
table_relative_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
    error_50ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.50, na.rm = TRUE), 2),
    MAE = round(mean(abs(slope_residuals_nobias)*100000, na.rm = TRUE), 2),
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals_nobias)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id)
  )

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(slope_abs, insitu_slope_nobias_m_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = c("source", "insitu_type"))


# RELATIVE REACH SLOPE TABLE BY RIVER
# -----------------------------------------------------
table_relative_reach_slope <- reach_SWOT_full_insitu %>%
  filter(source == "PGD0") %>%
  filter(insitu_type == "PT") %>%
  mutate(river = case_when(river %in% c("lowerPR", "upperPR") ~ "PR",TRUE ~ river)) %>%
  group_by(river) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
    error_50ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.50, na.rm = TRUE), 2),
    MAE = round(mean(abs(slope_residuals_nobias)*100000, na.rm = TRUE), 2),
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals_nobias)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(river) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(slope_abs, insitu_slope_nobias_m_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = "river")


# ABSOLUTE SLOPE TABLE BY VERSION
# -----------------------------------------------------
table_absolute_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(slope_residuals)*100000, 0.68, na.rm = TRUE), 2),
    error_50ile = round(quantile(abs(slope_residuals)*100000, 0.50, na.rm = TRUE), 2),
    MAE = round(mean(abs(slope_residuals)*100000, na.rm = TRUE), 2),
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(slope_abs, insitu_slope_m_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_reach_slope <- table_absolute_reach_slope %>%
  left_join(cor_table, by = "source")

# ABSOLUTE SLOPE TABLE BY GNSS/PT
# -----------------------------------------------------
table_absolute_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(slope_residuals)*100000, 0.68, na.rm = TRUE), 2),
    error_50ile = round(quantile(abs(slope_residuals)*100000, 0.50, na.rm = TRUE), 2),
    MAE = round(mean(abs(slope_residuals)*100000, na.rm = TRUE), 2),
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(slope_abs, insitu_slope_m_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_reach_slope <- table_absolute_reach_slope %>%
  left_join(cor_table, by = "insitu_type")





# RELATIVE SLOPE TABLE BY VERSION INCLUSION
# -----------------------------------------------------

# REACHES UNIQUE TO vC & vD
table_relative_reach_slope <- reach_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
    error_50ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.50, na.rm = TRUE), 2),
    MAE = round(mean(abs(slope_residuals_nobias)*100000, na.rm = TRUE), 2),
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals_nobias)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

cor_table <- reach_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(n = sum(complete.cases(slope_abs, insitu_slope_nobias_m_m)),
            r_value = if (n > 1) {round(cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"), 4)} else {
              NA_real_}, 
            p_value = if (n > 1) {cor.test(slope_abs, insitu_slope_nobias_m_m,
                                           method = "pearson")$p.value} else {NA_real_},
            .groups = "drop") %>%
  dplyr::select(-n)


# Join everything to one table
table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = "version_inclusion") %>%
  filter(version_inclusion != 0) %>% # drop 0, which are obs in both C&D
  mutate(version_inclusion = factor(version_inclusion, labels = c("vC0", "vD0")))

# SAME SUBSET
same_version_subset_reach_SWOT_insitu <- reach_SWOT_full_insitu %>%
  filter(version_inclusion == 0) %>%                            # keep only reaches present in both versions
  group_by(reach_id, insitu_time_utc, insitu_type) %>%
  filter(all(c("PIC0", "PGD0") %in% source)) %>%        # require both sources initially
  mutate(
    RiverSP_resid_na   = any(source == "PIC0"   & is.na(slope_residuals_nobias)),
    RiverTile_resid_na = any(source == "PGD0" & is.na(slope_residuals_nobias))
  ) %>%
  # drop the partner row when the counterpart has NA slope_residuals_nobias
  filter(
    !(source == "PGD0" & RiverSP_resid_na),
    !(source == "PIC0"   & RiverTile_resid_na)
  ) %>%
  # after removals, keep only triples that still contain both sources
  filter(all(c("PIC0", "PGD0") %in% source)) %>%
  ungroup() %>%
  dplyr::select(-RiverSP_resid_na, -RiverTile_resid_na)

# Add correlations
cor_table <- same_version_subset_reach_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(slope_abs, insitu_slope_nobias_m_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(slope_abs, insitu_slope_nobias_m_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

table_relative_reach_slope <- same_version_subset_reach_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    error_68ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
    error_50ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.50, na.rm = TRUE), 2),
    MAE = round(mean(abs(slope_residuals_nobias)*100000, na.rm = TRUE), 2),
    n = sum(!is.na(slope_residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id),
    .groups = "drop")

# Join everything to one table
table_relative_reach_slope <- table_relative_reach_slope %>%
  left_join(cor_table, by = c("source"))























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

# CDF plot 
ggplot(reach_SWOT_full_insitu, aes(x = abs(slope_residuals_nobias)*100000, color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = expression("|SWOT -" ~ italic("in situ") ~ "Slope| (cm/km)"), y = "Cumulative Probability", 
       title = "By SWOT version") +
  annotate("text", x = 4, y = 0.71, hjust = 0,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PIC0", ]$slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
                         "cm/km, vD:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PGD0", ]$slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
                         "cm/km"),
           color = "#222222", size = 5) +
  annotate("text", x = 4, y = 0.53, hjust = 0,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PIC0", ]$slope_residuals_nobias)*100000, 0.5, na.rm = TRUE), 2),
                         "cm/km, vD:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PGD0", ]$slope_residuals_nobias)*100000, 0.5, na.rm = TRUE), 2),
                         "cm/km"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1, vjust = 0,
           label = paste0("Version C: ", 
                          n_relative_df[n_relative_df$source == "PIC0", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$source == "PIC0", ]$n, " total"),
           color = "#E69F00", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1, vjust = 0, 
           label = paste0("Version D: ", 
                          n_relative_df[n_relative_df$source == "PGD0", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$source == "PGD0", ]$n, " total"), 
           color = "#0072B2", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("PIC0" = "#E69F00", "PGD0" = "#0072B2")) +
  scale_linetype_manual(values = c("PIC0" = "solid", "PGD0" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 13))
# width 7.17 height 6.35



# RELATIVE GNSS / PT with lines for C/D split
# --------------------------------------------------

# Compute n
n_relative_df <- reach_SWOT_full_insitu %>%
  filter(source == "PGD0") %>%
  group_by(insitu_type) %>%
  summarise(n_unique_reaches = n_distinct(reach_id[source == "PGD0"]), # count of non-NA residuals
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
           color = "#009E73", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1, vjust = 0, 
           label = paste0("GNSS: ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n, " total"), 
           color = "#CC79A7", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("PT" = "#009E73", "GNSS" = "#CC79A7")) +
  scale_linetype_manual(values = c("PT" = "solid", "GNSS" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 13))
# width 7.17 height 6.35




problems <- reach_SWOT_PT_vD %>%
  filter(river == "PR") 


ggplot(problems, aes(x = factor(reach_id), y = slope_residuals_nobias*100000)) + 
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  theme_minimal(base_size = 25)
  
  
problems_table <- problems %>%
  group_by(reach_id) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(slope_residuals_nobias)*100000, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(slope_residuals_nobias)*100000, na.rm = TRUE), 1),
    # count of non-NA residuals
    n = sum(!is.na(slope_residuals_nobias)),
    # count of unique reaches
    n_unique_reaches = n_distinct(reach_id))


# %>%
#   filter(slope_residuals_nobias*100000 > 3)
#   
#   filter(reach_id != '81270100061')


# Calculate the 68th percentile error
percentile_68_error <- quantile(abs(problems$slope_residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(problems$slope_residuals_nobias), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error:", percentile_68_error*100000))
print(paste("50th Percentile Error:", percentile_50_error*100000))


color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")


# correlation plot explorations:

reach_SWOT_PT_vD <- reach_SWOT_PT_vD %>%
  mutate(river = case_when(river %in% c("lowerPR", "upperPR") ~ "PR",TRUE ~ river))
reach_SWOT_GNSS_vD <- reach_SWOT_GNSS_vD %>%
  mutate(river = case_when(river %in% c("lowerPR", "upperPR") ~ "PR",TRUE ~ river))

# plot layovr_val
ggplot(reach_SWOT_PT_vD, aes(x = layovr_val, y = abs(slope_residuals_nobias)*100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("layovr_val") +
  ylab("slope error") +
  theme_minimal(base_size = 30) +
  labs(color = "River") 
ggplot(reach_SWOT_GNSS_vD, aes(x = layovr_val, y = abs(slope_residuals_nobias)*100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("layovr_val") +
  ylab("slope error") +
  ylim(0,7) +
  theme_minimal(base_size = 30) +
  labs(color = "River") 

# plot abs(xtrk_dist)
ggplot(reach_SWOT_PT_vD, aes(x = abs(xtrk_dist), y = abs(slope_residuals_nobias)*100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("abs(xtrk_dist)") +
  ylab("slope error") +
  theme_minimal(base_size = 30) +
  labs(color = "River") 
ggplot(reach_SWOT_GNSS_vD, aes(x = abs(xtrk_dist), y = abs(slope_residuals_nobias)*100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("abs(xtrk_dist)") +
  ylab("slope error") +
  theme_minimal(base_size = 30) +
  ylim(0,7) +
  labs(color = "River") 

# plot dark_frac
ggplot(reach_SWOT_PT_vD, aes(x = dark_frac, y = abs(slope_residuals_nobias)*100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("dark_frac") +
  ylab("slope error") +
  theme_minimal(base_size = 30) +
  labs(color = "River") 
ggplot(reach_SWOT_GNSS_vD, aes(x = dark_frac, y = abs(slope_residuals_nobias)*100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("dark_frac") +
  ylab("slope error") +
  ylim(0,7) +
  theme_minimal(base_size = 30) +
  labs(color = "River")

# plot slope_m_m_abs
ggplot(reach_SWOT_PT_vD, aes(x = slope_m_m_abs, y = abs(slope_residuals_nobias)*100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("slope_m_m") +
  ylab("slope error") +
  theme_minimal(base_size = 30) +
  labs(color = "River") 
ggplot(reach_SWOT_GNSS_vD, aes(x = reach_drift_slope_m_m_abs_nobias, y = abs(slope_residuals_nobias)*100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("slope_m_m") +
  ylab("slope error") +
  ylim(0,7) +
  theme_minimal(base_size = 30) +
  labs(color = "River") 

# plot width
ggplot(reach_SWOT_PT_vD, aes(x = width, y = abs(slope_residuals_nobias)*100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("width") +
  ylab("slope error") +
  theme_minimal(base_size = 30) +
  labs(color = "River") 
ggplot(reach_SWOT_GNSS_vD, aes(x = width, y = abs(slope_residuals_nobias)*100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("width") +
  ylab("slope error") +
  ylim(0,7) +
  theme_minimal(base_size = 30) +
  labs(color = "River") 

# plot n_good_nod
ggplot(reach_SWOT_GNSS_vD, aes(x = n_good_nod, y = abs(slope_residuals_nobias)*100000, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  xlab("n_good_nod") +
  ylab("slope error") +
  ylim(0,7) +
  theme_minimal(base_size = 30) +
  labs(color = "River") 




library(tidyverse)
library(lubridate)
library(dplyr)
library(ggtext)

# ---------------------------------------------------------------------------------------------------------------------------
# SWOT WSE & slope validation
# ---------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------
# Node level
# ---------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------
# READ IN DATA
# ---------------------------------------------------------------------------------------------------------------------------

# PT
node_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT")  %>%
  mutate(source = "PIC0") %>%
  rename(old_node_id = node_id) %>%
  filter(dark_frac < 0.5)
node_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "PGD0") %>%
  filter(dark_frac < 0.5)

# GNSS
node_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS") %>%
  mutate(source = "PIC0") %>%
  rename(old_node_id = node_id) %>%
  filter(dark_frac < 0.5)
node_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v17b/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS") %>%
  mutate(source = "PGD0") %>%
  filter(dark_frac < 0.5)

# ---------------------------------------------------------------------------------------------------------------------------
# Get all data to the same SWORD version (v17b)
# ---------------------------------------------------------------------------------------------------------------------------

# SWORD translator to change Version C data to SWORD v17b naming convention
SWORD_translator <- read_csv("/Users/camryn/Desktop/SWORD_translation/NA_NodeIDs_v17b_vs_v16.csv")

# translate the vC SWORD v16 data to SWORD v17b
node_SWOT_PT_vC <- node_SWOT_PT_vC %>%
  left_join(SWORD_translator %>% 
      select(v16_node_id, v17_node_id),
    by = c("old_node_id" = "v16_node_id")) %>%
  rename(node_id = v17_node_id)

node_SWOT_GNSS_vC <- node_SWOT_GNSS_vC %>%
  left_join(SWORD_translator %>% 
      select(v16_node_id, v17_node_id),
    by = c("old_node_id" = "v16_node_id")) %>%
  rename(node_id = v17_node_id)

# merge all dataframes together
node_SWOT_full_insitu <- bind_rows(node_SWOT_PT_vC, node_SWOT_PT_vD, node_SWOT_GNSS_vC, node_SWOT_GNSS_vD) %>%
  mutate(insitu_wse_m = coalesce(pt_wse_m, mean_node_drift_wse_m)) %>%
  mutate(insitu_wse_nobias_m = coalesce(pt_wse_nobias_m, mean_node_drift_wse_no_bias_m)) %>%
  mutate(insitu_time_utc = coalesce(pt_time_UTC, time_UTC))

# compute version inclusion
# -1 = only in vC, 0 = both, 1 = only in vD
all_nodes <- node_SWOT_full_insitu %>%
  distinct(node_id, source) %>%         
  group_by(node_id) %>%
  summarise(
    has_RiverSP   = any(source == "PIC0"),
    has_RiverTile = any(source == "PGD0"),
    .groups = "drop") %>%
  mutate(
    version_inclusion = case_when(has_RiverSP & has_RiverTile ~ 0L, has_RiverSP & !has_RiverTile ~ -1L, !has_RiverSP & has_RiverTile ~ 1L, TRUE ~ NA_integer_)) %>%
  select(node_id, version_inclusion)

node_SWOT_full_insitu <- node_SWOT_full_insitu %>%
  left_join(all_nodes, by = "node_id")

# ---------------------------------------------------------------------------------------------------------------------------
# TABLES -- SUMMARY STATS
# ---------------------------------------------------------------------------------------------------------------------------

# RELATIVE NODE WSE TABLE BY VERSION
# -----------------------------------------------------

table_relative_node_WSE <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals_nobias)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias)*100, na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE)*100, 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = "source")

# library(flextable)
# ft <- flextable(table_relative_node_WSE)
# save_as_docx(ft, path = "my_table.docx")

# BAR CHART of OBS COUNT

# relabel and reorder
table_relative_node_WSE <- table_relative_node_WSE %>%
  mutate(source = factor(source,
                         levels = c("PGD0", "PIC0"),  
                         labels = c("vD0", "vC0")))

# bar chart of count of residuals_nobias by version
ggplot(table_relative_node_WSE, aes(x = source, y = n, fill = source)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n), vjust = -0.5, size = 8) +
  ylab("Count") +
  coord_cartesian(ylim = c(2000, 7150)) +
  scale_fill_manual(values = c("vC0" = "#E69F00",
                               "vD0" = "#0072B2")) +
  theme_classic(base_size = 34) +
  theme(axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")
# 3.16, 6.54

# RELATIVE NODE WSE TABLE BY VERSION INCLUSION
# -----------------------------------------------------

# NODES UNIQUE TO vC & vD
table_relative_node_WSE <- node_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals_nobias)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias)*100, na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE)*100, 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = "version_inclusion") %>%
  filter(version_inclusion != 0) %>% # drop 0, which are obs in both C&D
  mutate(version_inclusion = factor(version_inclusion, labels = c("PIC0", "PGD0")))


# SAME SUBSET OF NODES FOR vC & vD
same_version_subset_node_SWOT_insitu <- node_SWOT_full_insitu %>%
  filter(version_inclusion == 0) %>%                            # keep only reaches present in both versions
  group_by(node_id, insitu_time_utc, insitu_type) %>%
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

table_relative_node_WSE <- same_version_subset_node_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals_nobias)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias)*100, na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE)*100, 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- same_version_subset_node_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = "source")




# RELATIVE NODE WSE TABLE BY GNSS/PT
# -----------------------------------------------------
table_relative_node_WSE <- node_SWOT_full_insitu %>%
  group_by(insitu_type, source) %>% 
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals_nobias)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias)*100, na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE)*100, 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_full_insitu %>%
  group_by(insitu_type, source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = c("insitu_type", "source"))


# RELATIVE NODE WSE TABLE BY RIVER
# -----------------------------------------------------
# for vD PT data

table_relative_node_WSE <- node_SWOT_full_insitu %>%
  filter(source == "PGD0") %>%
  filter(insitu_type == "PT") %>%
  mutate(river = case_when(river %in% c("lowerPR", "upperPR") ~ "PR",TRUE ~ river)) %>%
  group_by(river, insitu_type) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals_nobias)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias)*100, na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE)*100, 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_full_insitu %>%
  group_by(river) %>%
  mutate(river = case_when(river %in% c("lowerPR", "upperPR") ~ "PR", TRUE ~ river)) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = c("river"))

# ABSOLUTE WSE TABLE BY VERSION
# -----------------------------------------------------
table_absolute_node_WSE <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals)*100, na.rm = TRUE), 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_node_WSE <- table_absolute_node_WSE %>%
  left_join(cor_table, by = "source")

# ABSOLUTE WSE TABLE BY GNSS/PT
# -----------------------------------------------------
table_absolute_node_WSE <- node_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals)*100, na.rm = TRUE), 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_node_WSE <- table_absolute_node_WSE %>%
  left_join(cor_table, by = c("source", "insitu_type"))



# ---------------------------------------------------------------------------------------------------------------------------
# PLOTS -- SUMMARY STATS
# ---------------------------------------------------------------------------------------------------------------------------

# RELATIVE C vs D
# --------------------------------------------------

# Compute n
n_relative_df <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(n_unique_nodes = n_distinct(node_id), # count of non-NA residuals
            n = sum(!is.na(residuals_nobias)), .groups = "drop") # count of unique nodes

# CDF plot: Node Relative WSE differences
ggplot(node_SWOT_full_insitu, aes(x = abs(residuals_nobias)*100, color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = expression("|SWOT -" ~ italic("in situ") ~ "WSE| (cm)"), y = "Cumulative Probability",
    title = "By SWOT version") +
  annotate("text", x = 40, y = 0.71, hjust = 0,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "PIC0", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm, vD:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "PGD0", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  annotate("text", x = 40, y = 0.53, hjust = 0,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "PIC0", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm, vD:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "PGD0", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm"),
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



# RELATIVE GNSS / PT for vD
# --------------------------------------------------

node_SWOT_PGD0_insitu <- node_SWOT_full_insitu %>%
  filter(source == 'PGD0')

# Compute n
n_relative_df <- node_SWOT_PGD0_insitu %>%
  group_by(insitu_type) %>%
  summarise(n_unique_nodes = n_distinct(node_id[source == "PGD0"]), # count of non-NA residuals
            n = sum(!is.na(residuals_nobias)), .groups = "drop") # count of unique nodes

# CDF plot
ggplot(node_SWOT_PGD0_insitu, aes(x = abs(residuals_nobias)*100, color = insitu_type, linetype = insitu_type)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = expression("|SWOT -" ~ italic("in situ") ~ "WSE| (cm)"), y = "Cumulative Probability", 
       title = expression("By" ~ italic("in situ") ~ "measurement type")) +
  annotate("text", x = 40, y = 0.71, hjust = 0,
           label = paste("|68%ile| PT:", 
                         round(quantile(abs(node_SWOT_PGD0_insitu[node_SWOT_PGD0_insitu$insitu_type == "PT", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm, GNSS:", 
                         round(quantile(abs(node_SWOT_PGD0_insitu[node_SWOT_PGD0_insitu$insitu_type == "GNSS", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  annotate("text", x = 40, y = 0.53, hjust = 0,
           label = paste("|50%ile| PT:", 
                         round(quantile(abs(node_SWOT_PGD0_insitu[node_SWOT_PGD0_insitu$insitu_type == "PT", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm, GNSS:", 
                         round(quantile(abs(node_SWOT_PGD0_insitu[node_SWOT_PGD0_insitu$insitu_type == "GNSS", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1, vjust = 0,
           label = paste0("PT: ", 
                          n_relative_df[n_relative_df$insitu_type == "PT", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$insitu_type == "PT", ]$n, " total"),
           color = "#009E73", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1, vjust = 0, 
           label = paste0("GNSS: ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n, " total"), 
           color = "#CC79A7", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("PT" = "#009E73", "GNSS" = "#CC79A7")) +
  scale_linetype_manual(values = c("PT" = "solid", "GNSS" = "solid")) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 150))
# width 7.17 height 6.35






# ---------------------------------------------------------------------------------------------------------------------------
# ***************************************************************************************************************************
# Reach level
# ***************************************************************************************************************************
# ---------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------
# read in data
# ---------------------------------------------------------------------------------------------------------------------------

# PT
reach_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "PIC0") %>%
  rename(old_reach_id = reach_id) %>%
  filter(dark_frac < 0.5)
reach_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "PGD0") %>%
  filter(dark_frac < 0.5)

# GNSS
reach_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_GNSS.csv") %>%
  rename(old_reach_id = reach_id) %>%
  mutate(source = "PIC0") %>%
  filter(dark_frac < 0.5)
reach_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_wse_SWOT_GNSS.csv") %>%
  mutate(source = "PGD0") %>%
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
  mutate(insitu_wse_m = coalesce(mean_reach_pt_wse_m, mean_reach_drift_wse_m)) %>%
  mutate(insitu_wse_nobias_m = coalesce(pt_wse_nobias_m, mean_reach_drift_wse_no_bias_m)) %>%
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
  select(reach_id, version_inclusion)

reach_SWOT_full_insitu <- reach_SWOT_full_insitu %>%
  left_join(all_reaches, by = "reach_id")


# ---------------------------------------------------------------------------------------------------------------------------
# TABLES -- SUMMARY STATS
# ---------------------------------------------------------------------------------------------------------------------------

# RELATIVE REACH WSE TABLE BY VERSION
# -----------------------------------------------------
table_relative_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals_nobias)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias)*100, na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE)*100, 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = "source")

# relabel and reorder
table_relative_reach_WSE <- table_relative_reach_WSE %>%
  mutate(source = factor(source,
                         levels = c("PGD0", "PIC0"),   # swapped order
                         labels = c("vD0", "vC0")))            # relabels

# bar chart of count of residuals_nobias by version
ggplot(table_relative_reach_WSE, aes(x = source, y = n, fill = source)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n),
            vjust = -0.5,
            size = 8) +
  ylab("Count") +
  coord_cartesian(ylim = c(11, 185)) +
  scale_fill_manual(values = c("vC0" = "#E69F00",
                               "vD0" = "#0072B2")) +
  theme_classic(base_size = 34) +
  theme(axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")
# 3.16, 6.54


# RELATIVE REACH WSE TABLE BY VERSION INCLUSION
# -----------------------------------------------------

# REACHES UNIQUE TO vC & vD
table_relative_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals_nobias)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias)*100, na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE)*100, 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

cor_table <- reach_SWOT_full_insitu %>%
  group_by(version_inclusion) %>%
  summarise(n = sum(complete.cases(wse, insitu_wse_nobias_m)),
    r_value = if (n > 1) {round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4)} else {
      NA_real_}, 
    p_value = if (n > 1) {cor.test(wse, insitu_wse_nobias_m,
               method = "pearson")$p.value} else {NA_real_},
    .groups = "drop")

# Join everything to one table
table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = "version_inclusion") %>%
  filter(version_inclusion != 0) %>% # drop 0, which are obs in both C&D
  mutate(version_inclusion = factor(version_inclusion, labels = c("vC0", "vD0")))

# SAME SUBSET
same_version_subset_reach_SWOT_insitu <- reach_SWOT_full_insitu %>%
  filter(version_inclusion == 0) %>%                            # keep only reaches present in both versions
  group_by(reach_id, insitu_time_utc, insitu_type) %>%
  filter(all(c("PIC0", "PGD0") %in% source)) %>%        # require both sources initially
  mutate(
    RiverSP_resid_na   = any(source == "PIC0"   & is.na(residuals_nobias)),
    RiverTile_resid_na = any(source == "PGD0" & is.na(residuals_nobias))
  ) %>%
  # drop the partner row when the counterpart has NA residuals_nobias
  filter(
    !(source == "PGD0" & RiverSP_resid_na),
    !(source == "PIC0"   & RiverTile_resid_na)
  ) %>%
  # after removals, keep only triples that still contain both sources
  filter(all(c("PIC0", "PGD0") %in% source)) %>%
  ungroup() %>%
  select(-RiverSP_resid_na, -RiverTile_resid_na)

# Add correlations
cor_table <- same_version_subset_reach_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

table_relative_reach_WSE <- same_version_subset_reach_SWOT_insitu %>%
  group_by(source) %>%
  summarise(
    error_68ile = round(quantile(abs(residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals_nobias)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias)*100, na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE)*100, 2),
    n = sum(!is.na(residuals_nobias)),
    n_unique_reaches = n_distinct(reach_id),
    .groups = "drop")

# Join everything to one table
table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = c("source"))



# RELATIVE REACH WSE TABLE BY GNSS/PT
# -----------------------------------------------------
table_relative_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(insitu_type, source) %>%
  summarise(
    # error metrics (use all data)
    error_68ile = round(quantile(abs(residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals_nobias)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias)*100, na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE)*100, 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(insitu_type, source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = c("insitu_type", "source"))





# RELATIVE REACH WSE TABLE BY RIVER
# -----------------------------------------------------
table_relative_reach_WSE <- reach_SWOT_full_insitu %>%
  filter(source == "PGD0") %>%
  filter(insitu_type == "PT") %>%
  mutate(river = case_when(river %in% c("lowerPR", "upperPR") ~ "PR",TRUE ~ river)) %>%
  group_by(river) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals_nobias)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals_nobias)*100, na.rm = TRUE), 1),
    bias = round(median(bias, na.rm = TRUE)*100, 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique reaches
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(river) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = "river")

# ABSOLUTE WSE TABLE BY VERSION
# -----------------------------------------------------
table_absolute_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals)*100, na.rm = TRUE), 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_reach_WSE <- table_absolute_reach_WSE %>%
  left_join(cor_table, by = "source")

# ABSOLUTE WSE TABLE BY GNSS/PT
# -----------------------------------------------------
table_absolute_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    # error metrics
    error_68ile = round(quantile(abs(residuals)*100, 0.68, na.rm = TRUE), 1),
    error_50ile = round(quantile(abs(residuals)*100, 0.50, na.rm = TRUE), 1),
    MAE = round(mean(abs(residuals)*100, na.rm = TRUE), 1),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source, insitu_type) %>%
  summarise(
    r_value = round(cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"), 4),
    p_value = tryCatch(cor.test(wse, insitu_wse_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_reach_WSE <- table_absolute_reach_WSE %>%
  left_join(cor_table, by = c("source", "insitu_type"))


# ---------------------------------------------------------------------------------------------------------------------------
# PLOTS -- SUMMARY STATS
# ---------------------------------------------------------------------------------------------------------------------------

# RELATIVE C vs D
# --------------------------------------------------

# Compute n
n_relative_df <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(n_unique_reaches = n_distinct(reach_id), # count of non-NA residuals
            n = sum(!is.na(residuals_nobias)), .groups = "drop") # count of unique nodes

# Reach Relative WSE Differences
ggplot(reach_SWOT_full_insitu, aes(x = abs(residuals_nobias)*100, color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = expression("|SWOT -" ~ italic("in situ") ~ "WSE| (cm)"), y = "Cumulative Probability", 
       title = "By SWOT version") +
  annotate("text", x = 40, y = 0.71, hjust = 0,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PIC0", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm, vD:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PGD0", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  annotate("text", x = 40, y = 0.53, hjust = 0,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PIC0", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm, vD:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "PGD0", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm"),
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
  coord_cartesian(xlim = c(0, 150))
# width 7.17 height 6.35



# RELATIVE GNSS / PT with lines for C/D split
# --------------------------------------------------

reach_SWOT_PGD0_insitu <- reach_SWOT_full_insitu %>%
  filter(source == 'PGD0')

# Compute n
n_relative_df <- reach_SWOT_PGD0_insitu %>%
  group_by(insitu_type) %>%
  summarise(n_unique_reaches = n_distinct(reach_id[source == "PGD0"]), # count of non-NA residuals
            n = sum(!is.na(residuals_nobias)), .groups = "drop") # count of unique nodes

# CDF plot (fixed cm placement)
ggplot(reach_SWOT_PGD0_insitu, aes(x = abs(residuals_nobias)*100, color = insitu_type, linetype = insitu_type)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = expression("|SWOT -" ~ italic("in situ") ~ "WSE| (cm)"), y = "Cumulative Probability", 
       title = expression("By" ~ italic("in situ") ~ "measurement type")) +
  annotate("text", x = 40, y = 0.71, hjust = 0,
           label = paste("|68%ile| PT:", 
                         round(quantile(abs(reach_SWOT_PGD0_insitu[reach_SWOT_PGD0_insitu$insitu_type == "PT", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm, GNSS:", 
                         round(quantile(abs(reach_SWOT_PGD0_insitu[reach_SWOT_PGD0_insitu$insitu_type == "GNSS", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  annotate("text", x = 40, y = 0.53, hjust = 0,
           label = paste("|50%ile| PT:", 
                         round(quantile(abs(reach_SWOT_PGD0_insitu[reach_SWOT_PGD0_insitu$insitu_type == "PT", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm, GNSS:", 
                         round(quantile(abs(reach_SWOT_PGD0_insitu[reach_SWOT_PGD0_insitu$insitu_type == "GNSS", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm"),
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
  theme(legend.position = "none")  +
  coord_cartesian(xlim = c(0, 150))
# width 7.17 height 6.35






# ---------------------------------------------------------------------------------------------------------------------------
# INTER-RIVER COMPARISON PLOTS
# ---------------------------------------------------------------------------------------------------------------------------


# WSE
# -------------------------------------------

node_SWOT_PT_vD <- node_SWOT_PT_vD %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",   # merge into one
    TRUE ~ river                              # keep all others unchanged
  ))

# Reorder the factor levels for river
node_SWOT_PT_vD$river <- factor(
  node_SWOT_PT_vD$river,
  levels = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
)

color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")

# compute n for each river
counts <- node_SWOT_PT_vD %>%
  group_by(river) %>%
  summarise(n = n()) %>%
  ungroup()

# Replot
ggplot(node_SWOT_PT_vD, aes(x = river, y = abs(residuals_nobias)*100, fill = river)) + 
  geom_violin(alpha = 0.8, color = NA) +
  xlab("River") +
  ylab("| SWOT - PT WSE | (cm)") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  # add counts below each violin
  geom_text(data = counts,
            aes(x = river, y = -1, label = paste0("n=", n)),
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
  coord_cartesian(ylim = c(-5, 60))
# 9.44, 6.01
# **********************************

# Replot
ggplot(node_SWOT_PT_vD, aes(x = river, y = bias*100, fill = river)) + 
  geom_violin(alpha = 0.8, color = NA) +
  xlab("River") +
  ylab("| SWOT - PT WSE | (cm)") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  # add counts below each violin
  geom_text(data = counts,
            aes(x = river, y = -1, label = paste0("n=", n)),
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
  coord_cartesian(ylim = c(-25, 50))




# Slope
# -------------------------------------------
reach_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v17b/reach_slope_SWOT_GNSS.csv") %>%
  filter(dark_frac < 0.5)
reach_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v17b/reach_slope_SWOT_PT.csv") %>%
  filter(dark_frac < 0.5)

# merge all dataframes together
reach_SWOT_full_insitu <- bind_rows(reach_SWOT_PT_vD, reach_SWOT_GNSS_vD) %>%
  mutate(insitu_slope_m_m = coalesce(slope_m_m_abs, reach_drift_slope_m_m_abs)) %>%
  mutate(insitu_slope_nobias_m_m = coalesce(mean_reach_PT_slope_no_bias_m_m, reach_drift_slope_m_m_abs_nobias)) %>%
  mutate(insitu_time_utc = coalesce(pt_time_UTC, wse_drift_midpoint_UTC))

reach_SWOT_full_insitu <- reach_SWOT_full_insitu %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",   # merge into one
    TRUE ~ river                              # keep all others unchanged
  ))

# Reorder the factor levels for river
reach_SWOT_full_insitu$river <- factor(
  reach_SWOT_full_insitu$river,
  levels = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
)

color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")

# compute n for each river
counts <- reach_SWOT_full_insitu %>%
  group_by(river) %>%
  summarise(n = n()) %>%
  ungroup()

# Replot
ggplot(reach_SWOT_full_insitu, aes(x = river, y = abs(slope_residuals_nobias)*100000, fill = river)) + 
  geom_violin(alpha = 0.8, color = NA) +
  xlab("River") +
  ylab("| SWOT -" ~ italic("in situ") ~ "Slope | (cm/km)") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  # add counts below each violin
  geom_text(data = counts,
            aes(x = river, y = -0.2, label = paste0("n=", n)),
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
  coord_cartesian(ylim = c(-0.5, 12))
# 9.44, 6.01







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

# PT
node_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT")
node_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/node_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT")

# GNSS
node_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverSP_v16/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS")
node_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/node/RiverTile_v17b/node_SWOT_GNSS_3mdiff.csv") %>%
  mutate(insitu_type = "GNSS")

# merge all dataframes together
node_SWOT_full_insitu <- bind_rows(node_SWOT_PT_vC, node_SWOT_PT_vD, node_SWOT_GNSS_vC, node_SWOT_GNSS_vD) %>%
  mutate(insitu_wse_m = coalesce(pt_wse_m, mean_node_drift_wse_m)) %>%
  mutate(insitu_wse_nobias_m = coalesce(pt_wse_nobias_m, mean_node_drift_wse_no_bias_m))

# ---------------------------------------------------------------------------------------------------------------------------
# TABLES -- SUMMARY STATS
# ---------------------------------------------------------------------------------------------------------------------------

# RELATIVE NODE WSE TABLE BY VERSION
# -----------------------------------------------------
table_relative_node_WSE <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_68ile = quantile(abs(residuals_nobias), 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(residuals_nobias), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals_nobias), na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = "source")

# RELATIVE NODE WSE TABLE BY GNSS/PT
# -----------------------------------------------------
table_relative_node_WSE <- node_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    # error metrics (use all data)
    error_68ile = quantile(abs(residuals_nobias), 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(residuals_nobias), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals_nobias), na.rm = TRUE),

    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    
    # !!!!!!!!!!!!!! how to quantify count when this will vary from vC to vD
    # right now I'm taking RiverTile which has more unique obs
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id[source == "RiverTile"])
  )

# Add correlations
cor_table <- node_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    r_value = cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = "insitu_type")

# RELATIVE NODE WSE TABLE BY RIVER
# -----------------------------------------------------
table_relative_node_WSE <- node_SWOT_full_insitu %>%
  group_by(river) %>%
  summarise(
    # error metrics
    error_68ile = quantile(abs(residuals_nobias), 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(residuals_nobias), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals_nobias), na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    
    # !!!!!!!!!!!!!! how to quantify count when this will vary from vC to vD
    # right now I'm taking RiverTile which has more unique obs
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id[source == "RiverTile"]))

# Add correlations
cor_table <- node_SWOT_full_insitu %>%
  group_by(river) %>%
  summarise(
    r_value = cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_node_WSE <- table_relative_node_WSE %>%
  left_join(cor_table, by = "river")

# ABSOLUTE WSE TABLE BY VERSION
# -----------------------------------------------------
table_absolute_node_WSE <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_68ile = quantile(abs(residuals), 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(residuals), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals), na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(wse, insitu_wse_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_node_WSE <- table_absolute_node_WSE %>%
  left_join(cor_table, by = "source")

# ABSOLUTE WSE TABLE BY GNSS/PT
# -----------------------------------------------------
table_absolute_node_WSE <- node_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    # error metrics
    error_68ile = quantile(abs(residuals), 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(residuals), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals), na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_nodes = n_distinct(node_id))

# Add correlations
cor_table <- node_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    r_value = cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(wse, insitu_wse_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_node_WSE <- table_absolute_node_WSE %>%
  left_join(cor_table, by = "insitu_type")



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

# CDF plot (fixed cm placement)
ggplot(node_SWOT_full_insitu, aes(x = abs(residuals_nobias)*100, color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - in situ WSE (cm)", y = "Cumulative Probability", 
       title = "Node Relative WSE Differences") +
  annotate("text", x = 95, y = 0.71,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "RiverSP", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm, vD:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "RiverTile", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  annotate("text", x = 95, y = 0.53,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "RiverSP", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm, vD:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$source == "RiverTile", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1.1, vjust = 0,
           label = paste0("Version C: ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n, " total"),
           color = "#E97132", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1.1, vjust = 0, 
           label = paste0("Version D: ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n, " total"), 
           color = "darkblue", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "#E97132", "RiverTile" = "darkblue")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "solid")) +
  theme(legend.position = "none") +
  xlim(0, 150)
# width 610 height 550



# RELATIVE GNSS / PT with lines for C/D split
# --------------------------------------------------

# Compute n
n_relative_df <- node_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(n_unique_nodes = n_distinct(node_id[source == "RiverTile"]), # count of non-NA residuals
            n = sum(!is.na(residuals_nobias)), .groups = "drop") # count of unique nodes

# CDF plot (fixed cm placement)
ggplot(node_SWOT_full_insitu, aes(x = abs(residuals_nobias)*100, color = insitu_type, linetype = insitu_type)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - in situ WSE (cm)", y = "Cumulative Probability", 
       title = "Node Relative WSE Differences") +
  annotate("text", x = 95, y = 0.71,
           label = paste("|68%ile| PT:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$insitu_type == "PT", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm, GNSS:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$insitu_type == "GNSS", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  annotate("text", x = 95, y = 0.53,
           label = paste("|50%ile| PT:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$insitu_type == "PT", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm, GNSS:", 
                         round(quantile(abs(node_SWOT_full_insitu[node_SWOT_full_insitu$insitu_type == "GNSS", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1.1, vjust = 0,
           label = paste0("PT: ", 
                          n_relative_df[n_relative_df$insitu_type == "PT", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$insitu_type == "PT", ]$n, " total"),
           color = "#C03F61", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1.1, vjust = 0, 
           label = paste0("GNSS: ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n_unique_nodes, 
                          " unique nodes, ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n, " total"), 
           color = "#404A22", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("PT" = "#C03F61", "GNSS" = "#404A22")) +
  scale_linetype_manual(values = c("PT" = "solid", "GNSS" = "solid")) +
  theme(legend.position = "none") +
  xlim(0, 150)
# width 610 height 550





# ---------------------------------------------------------------------------------------------------------------------------
# Reach level
# ---------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------
# read in data
# ---------------------------------------------------------------------------------------------------------------------------

# PT
reach_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverSP")
reach_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_wse_SWOT_PT.csv") %>%
  mutate(insitu_type = "PT") %>%
  mutate(source = "RiverTile")

# GNSS
reach_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_SWOT_GNSS.csv") 
reach_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_SWOT_GNSS.csv")

# merge all dataframes together
reach_SWOT_full_insitu <- bind_rows(reach_SWOT_PT_vC, reach_SWOT_PT_vD, reach_SWOT_GNSS_vC, reach_SWOT_GNSS_vD) %>%
  mutate(insitu_wse_m = coalesce(mean_reach_pt_wse_m, mean_reach_drift_wse_m)) %>%
  mutate(insitu_wse_nobias_m = coalesce(pt_wse_nobias_m, mean_reach_drift_wse_no_bias_m))


# ---------------------------------------------------------------------------------------------------------------------------
# TABLES -- SUMMARY STATS
# ---------------------------------------------------------------------------------------------------------------------------

# RELATIVE REACH WSE TABLE BY VERSION
# -----------------------------------------------------
table_relative_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    # error metrics
    error_68ile = quantile(abs(residuals_nobias), 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(residuals_nobias), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals_nobias), na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = "source")

# RELATIVE REACH WSE TABLE BY GNSS/PT
# -----------------------------------------------------
table_relative_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    # error metrics (use all data)
    error_68ile = quantile(abs(residuals_nobias), 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(residuals_nobias), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals_nobias), na.rm = TRUE),
    
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    
    # !!!!!!!!!!!!!! how to quantify count when this will vary from vC to vD
    # right now I'm taking RiverTile which has more unique obs
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id[source == "RiverTile"])
  )

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    r_value = cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(wse, insitu_wse_nobias_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_relative_reach_WSE <- table_relative_reach_WSE %>%
  left_join(cor_table, by = "insitu_type")

# RELATIVE REACH WSE TABLE BY RIVER
# -----------------------------------------------------
table_relative_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(river) %>%
  summarise(
    # error metrics
    error_68ile = quantile(abs(residuals_nobias), 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(residuals_nobias), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals_nobias), na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(residuals_nobias)),
    
    # !!!!!!!!!!!!!! how to quantify count when this will vary from vC to vD
    # right now I'm taking RiverTile which has more unique obs
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id[source == "RiverTile"]))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(river) %>%
  summarise(
    r_value = cor(wse, insitu_wse_nobias_m, use = "complete.obs", method = "pearson"),
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
    error_68ile = quantile(abs(residuals), 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(residuals), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals), na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(source) %>%
  summarise(
    r_value = cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(wse, insitu_wse_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_reach_WSE <- table_absolute_reach_WSE %>%
  left_join(cor_table, by = "source")

# ABSOLUTE WSE TABLE BY GNSS/PT
# -----------------------------------------------------
table_absolute_reach_WSE <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    # error metrics
    error_68ile = quantile(abs(residuals), 0.68, na.rm = TRUE),
    error_50ile = quantile(abs(residuals), 0.50, na.rm = TRUE),
    MAE = mean(abs(residuals), na.rm = TRUE),
    # count of non-NA residuals
    n = sum(!is.na(residuals)),
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id))

# Add correlations
cor_table <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(
    r_value = cor(wse, insitu_wse_m, use = "complete.obs", method = "pearson"),
    p_value = tryCatch(cor.test(wse, insitu_wse_m)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Join everything to one table
table_absolute_reach_WSE <- table_absolute_reach_WSE %>%
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
            n = sum(!is.na(residuals_nobias)), .groups = "drop") # count of unique nodes

# CDF plot (fixed cm placement)
ggplot(reach_SWOT_full_insitu, aes(x = abs(residuals_nobias)*100, color = source, linetype = source)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - in situ WSE (cm)", y = "Cumulative Probability", 
       title = "Reach Relative WSE Differences") +
  annotate("text", x = 95, y = 0.71,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverSP", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm, vD:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverTile", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  annotate("text", x = 95, y = 0.53,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverSP", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm, vD:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverTile", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1.1, vjust = 0,
           label = paste0("Version C: ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n, " total"),
           color = "#E97132", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1.1, vjust = 0, 
           label = paste0("Version D: ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n, " total"), 
           color = "darkblue", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "#E97132", "RiverTile" = "darkblue")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "solid")) +
  theme(legend.position = "none") +
  xlim(0, 150)
# width 610 height 550



# RELATIVE GNSS / PT with lines for C/D split
# --------------------------------------------------

# Compute n
n_relative_df <- reach_SWOT_full_insitu %>%
  group_by(insitu_type) %>%
  summarise(n_unique_reaches = n_distinct(reach_id[source == "RiverTile"]), # count of non-NA residuals
            n = sum(!is.na(residuals_nobias)), .groups = "drop") # count of unique nodes

# CDF plot (fixed cm placement)
ggplot(reach_SWOT_full_insitu, aes(x = abs(residuals_nobias)*100, color = insitu_type, linetype = insitu_type)) +
  stat_ecdf(geom = "step", size = 1.2) +
  geom_hline(yintercept = 0.68, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey") +
  labs(x = "SWOT - in situ WSE (cm)", y = "Cumulative Probability", 
       title = "Reach Relative WSE Differences") +
  annotate("text", x = 95, y = 0.71,
           label = paste("|68%ile| PT:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "PT", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm, GNSS:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "GNSS", ]$residuals_nobias)*100, 0.68, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  annotate("text", x = 95, y = 0.53,
           label = paste("|50%ile| PT:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "PT", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm, GNSS:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "GNSS", ]$residuals_nobias)*100, 0.5, na.rm = TRUE), 1),
                         "cm"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1.1, vjust = 0,
           label = paste0("PT: ", 
                          n_relative_df[n_relative_df$insitu_type == "PT", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$insitu_type == "PT", ]$n, " total"),
           color = "#C03F61", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1.1, vjust = 0, 
           label = paste0("GNSS: ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n, " total"), 
           color = "#404A22", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("PT" = "#C03F61", "GNSS" = "#404A22")) +
  scale_linetype_manual(values = c("PT" = "solid", "GNSS" = "solid")) +
  theme(legend.position = "none") +
  xlim(0, 150)
# width 610 height 550





# ---------------------------------------------------------------------------------------------------------------------------
# SWOT Slope validation
# ---------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------
# read in data
# ---------------------------------------------------------------------------------------------------------------------------

# PT
reach_SWOT_PT_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_slope_SWOT_PT.csv") 
reach_SWOT_PT_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_slope_SWOT_PT.csv") 

# GNSS
reach_SWOT_GNSS_vC <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverSP_v16/reach_SWOT_GNSS.csv") 
reach_SWOT_GNSS_vD <- read_csv("/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_SWOT_GNSS.csv")

# merge all dataframes together
reach_SWOT_full_insitu <- bind_rows(reach_SWOT_PT_vC, reach_SWOT_PT_vD, reach_SWOT_GNSS_vC, reach_SWOT_GNSS_vD) %>%
  mutate(insitu_slope_m_m = coalesce(slope_m_m_abs, reach_drift_slope_m_m_abs)) %>%
  mutate(insitu_slope_nobias_m_m = coalesce(mean_reach_PT_slope_no_bias_m_m, reach_drift_slope_m_m_abs_nobias))


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
    
    # !!!!!!!!!!!!!! how to quantify count when this will vary from vC to vD
    # right now I'm taking RiverTile which has more unique obs
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id[source == "RiverTile"])
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
    
    # !!!!!!!!!!!!!! how to quantify count when this will vary from vC to vD
    # right now I'm taking RiverTile which has more unique obs
    # count of unique nodes
    n_unique_reaches = n_distinct(reach_id[source == "RiverTile"]))

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
  labs(x = "SWOT - in situ Slope (cm/km)", y = "Cumulative Probability", 
       title = "By SWOT Version") +
  annotate("text", x = 8.5, y = 0.71,
           label = paste("|68%ile| vC:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverSP", ]$slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
                         "cm/km, vD:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverTile", ]$slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
                         "cm/km"),
           color = "#222222", size = 5) +
  annotate("text", x = 8.5, y = 0.53,
           label = paste("|50%ile| vC:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverSP", ]$slope_residuals_nobias)*100000, 0.5, na.rm = TRUE), 2),
                         "cm/km, vD:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$source == "RiverTile", ]$slope_residuals_nobias)*100000, 0.5, na.rm = TRUE), 2),
                         "cm/km"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1.1, vjust = 0,
           label = paste0("Version C: ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$source == "RiverSP", ]$n, " total"),
           color = "#E97132", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1.1, vjust = 0, 
           label = paste0("Version D: ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$source == "RiverTile", ]$n, " total"), 
           color = "darkblue", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("RiverSP" = "#E97132", "RiverTile" = "darkblue")) +
  scale_linetype_manual(values = c("RiverSP" = "solid", "RiverTile" = "solid")) +
  theme(legend.position = "none") +
  xlim(0, 15)
# width 610 height 550



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
  labs(x = "SWOT - in situ Slope (cm/km)", y = "Cumulative Probability", 
       title = "By in situ measurement type") +
  annotate("text", x = 8.5, y = 0.71,
           label = paste("|68%ile| PT:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "PT", ]$slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
                         "cm/km, GNSS:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "GNSS", ]$slope_residuals_nobias)*100000, 0.68, na.rm = TRUE), 2),
                         "cm/km"),
           color = "#222222", size = 5) +
  annotate("text", x = 8.5, y = 0.53,
           label = paste("|50%ile| PT:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "PT", ]$slope_residuals_nobias)*100000, 0.5, na.rm = TRUE), 2),
                         "cm/km, GNSS:", 
                         round(quantile(abs(reach_SWOT_full_insitu[reach_SWOT_full_insitu$insitu_type == "GNSS", ]$slope_residuals_nobias)*100000, 0.5, na.rm = TRUE), 2),
                         "cm/km"),
           color = "#222222", size = 5) +
  # Add counts in lower right
  annotate("text", x = Inf, y = 0.08,
           hjust = 1.1, vjust = 0,
           label = paste0("PT: ", 
                          n_relative_df[n_relative_df$insitu_type == "PT", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$insitu_type == "PT", ]$n, " total"),
           color = "#C03F61", size = 5) +
  annotate("text", x = Inf, y = 0.02, 
           hjust = 1.1, vjust = 0, 
           label = paste0("GNSS: ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n_unique_reaches, 
                          " unique reaches, ", 
                          n_relative_df[n_relative_df$insitu_type == "GNSS", ]$n, " total"), 
           color = "#404A22", size = 5) +
  theme_minimal(base_size = 18) +
  scale_color_manual(values = c("PT" = "#C03F61", "GNSS" = "#404A22")) +
  scale_linetype_manual(values = c("PT" = "solid", "GNSS" = "solid")) +
  theme(legend.position = "none") +
  xlim(0, 15)
# width 610 height 550





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

# "CD", "CL", "lowerPR", "lowerYR", "SJ", "upperPR", "upperYR"
color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#DA627D")

color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")

# Replot
ggplot(node_SWOT_PT_vD, aes(x = river, y = abs(residuals_nobias)*100, fill = river)) + 
  geom_violin(alpha = 0.8, color = NA) +
  xlab("River") +
  ylab("|SWOT - PT WSE| (cm)") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  theme_minimal(base_size = 25) +
  scale_fill_manual(
    values = color_palette,
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  scale_x_discrete(
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 0.9)
  ) + 
  coord_cartesian(ylim = c(0, 60))

# **********************************




# Slope
# -------------------------------------------

reach_SWOT_PT_vD <- reach_SWOT_PT_vD %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",   # merge into one
    TRUE ~ river                              # keep all others unchanged
  ))

# Reorder the factor levels for river
reach_SWOT_PT_vD$river <- factor(
  reach_SWOT_PT_vD$river,
  levels = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR")
)

# "CD", "CL", "lowerPR", "lowerYR", "SJ", "upperPR", "upperYR"
color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#DA627D")

color_palette <- c("#F2C14E", "#8EAD7A", "#3B6064", "#F4845F", "#DA627D", "#9A348E")

# Replot
ggplot(reach_SWOT_PT_vD, aes(x = river, y = abs(slope_residuals_nobias)*100000, fill = river)) + 
  geom_violin(alpha = 0.8, color = NA) +
  xlab("River") +
  ylab("|SWOT - PT Slope| (cm/km)") +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 3, lwd = 1) +
  theme_minimal(base_size = 25) +
  scale_fill_manual(
    values = color_palette,
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  scale_x_discrete(
    breaks = c("CL", "SJ", "CD", "PR", "upperYR", "lowerYR"),
    labels = c("Coleen", "Sheenjek", "Chandalar", "Porcupine", "Single-channel Yukon", "Braided Yukon")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 0.9)
  ) + 
  coord_cartesian(ylim = c(0, 6))
# 1000, 515

# **********************************







# correlation test
cor_test_nobias <- cor.test(node_SWOT_full_insitu$wse, node_SWOT_full_insitu$pt_wse_nobias_m)

# Extract r and p-value
r_value_nobias <- cor_test_nobias$estimate # Pearson correlation coefficient
p_value_nobias <- cor_test_nobias$p.value # 





node_SWOT_full_insitu <- node_SWOT_full_insitu %>%
  mutate(river = case_when(
    river %in% c("lowerPR", "upperPR") ~ "PR",   # merge into one
    TRUE ~ river                              # keep all others unchanged
  ))


# "CD", "CL", "lowerPR", "lowerYR", "SJ", "upperPR", "upperYR"
color_palette <- c("#3B6064", "#F2C14E", "#F4845F", "#9A348E", "#8EAD7A", "#DA627D", "red")

# plot SWOT vs PT wse
ggplot(node_SWOT_full_insitu, 
       aes(x = pt_wse_nobias_m, y = wse, color = factor(river))) +
  geom_point(size = 4) +
  scale_color_manual(
    values = color_palette,
    breaks = c("CD", "CL", "lowerPR", "lowerYR", "SJ", "upperYR"),
    labels = c("Chandalar", "Coleen", "Porcupine", "Braided Yukon", "Sheenjek", "Single-channel Yukon")) +
  xlab("PT WSE (m)") +
  ylab("SWOT WSE (m)") +
  theme_minimal(base_size = 30) +
  geom_abline(linetype = "dashed", color = "gray") +  # 1:1 line
  annotate("text", 
           x = min(node_SWOT_full_insitu$pt_wse_nobias_m, na.rm = TRUE), 
           y = max(node_SWOT_full_insitu$wse, na.rm = TRUE), 
           label = paste0("r = ", round(r_value_nobias, 4), 
                          "\np value = ", signif(p_value_nobias, 3),
                          "\nn = ", nrow(node_SWOT_full_insitu)),
           hjust = 0, vjust = 1, size = 8) +
  labs(color = "River")


# need to label each version & add river id


# add river names to df
reach_SWOT_GNSS_vD <- reach_SWOT_GNSS_vD %>%
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
# Check to make sure rivers are labeled as expected
# reach_SWOT_GNSS_vD %>%
#   filter(reach_id %in% c("81270100111", "81270100121", "81270100131", "81270100141", "81270100151", "81270100161", "81270200011", "81270200021")) %>%
#   select(reach_id, river_code, river)

# save to CSV
# write.csv(reach_SWOT_GNSS_vD, file = '/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/wse/reach/RiverTile_v17b/reach_SWOT_GNSS.csv', row.names = FALSE)


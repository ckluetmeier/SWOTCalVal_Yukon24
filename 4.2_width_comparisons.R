# all clusters comparison
# ---------------------------------------------------------------------------------------------------------------------------

# Set working directory to the orthomosaics directory with the combined_ortho_SWOT_ csvs
wd <- "/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats/CalVal_dataframes/width/node/RiverSP_v16"
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
combined_ortho_SWORD_df <- combined_df %>%
  mutate(
    river_code = substr(reach_id, 1, 6),
    river = case_when(
      river_code == "812701" ~ "lower_YR",
      river_code == "812705" ~ "upper_YR",
      river_code == "812508" ~ "CD",
      river_code == "812603" ~ "upper_PR",
      river_code == "812605" ~ "upper_PR",
      river_code == "812604" ~ "CL",
      TRUE ~ NA_character_
    )
  )

# optional filter to look at isolated groups (e.g. by river)
# combined_ortho_SWORD_df <- combined_ortho_SWORD_df %>%
#   filter(river == "upper_YR")

# Create a source column to identify which dataset each row comes from
RiverSP_df <- combined_ortho_SWORD_df %>%
  mutate(source = "RiverSP")

RiverTile_df <- combined_ortho_SWORD_df %>%
  mutate(source = "RiverTile")

# optional filter to look at isolated groups (e.g. by river)
# RiverTile_df_subset <- RiverTile_df %>%
#   filter(river == "lower_YR")

percentile_68_error_RiverTile <- quantile(abs(RiverTile_df_subset$percent_diff), 0.68, na.rm=TRUE)
percentile_68_error_RiverTile <- quantile(abs(RiverTile_df_subset$residuals), 0.68, na.rm=TRUE)
print(percentile_68_error_RiverTile)

# Combine both dataframes
combined_ortho_SWORD_df <- bind_rows(RiverSP_df, RiverTile_df)

# ---------------------------------------------------------------------------------------------------------------------------
# Stats & plots
# ---------------------------------------------------------------------------------------------------------------------------

# Calculate the 68th percentile error
percentile_68_error <- quantile(abs(combined_ortho_SWORD_df$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(combined_ortho_SWORD_df$residuals), 0.50, na.rm=TRUE)
# for % diff
percentile_68_percent <- quantile(abs(combined_ortho_SWORD_df$percent_diff), 0.68, na.rm=TRUE)
percentile_50_percent <- quantile(abs(combined_ortho_SWORD_df$percent_diff), 0.50, na.rm=TRUE)

#print the result
print(paste("68th Percentile Error:", percentile_68_percent))
print(paste("50th Percentile Error:", percentile_50_percent))

# correlation test
cor_test <- cor.test(combined_ortho_SWORD_df$width, combined_ortho_SWORD_df$ortho_width_m)

# Extract r and p-value
r_value <- cor_test$estimate # Pearson correlation coefficient
p_value <- cor_test$p.value # highly statistically significant is P < 0.001


# Stats for RiverSP / RiverTile
# Calculate the 68th & 50th percentile error
percentile_68_error <- quantile(abs(RiverSP_df$residuals), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(RiverSP_df$residuals), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile <- quantile(abs(RiverTile_df$residuals), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile <- quantile(abs(RiverTile_df$residuals), 0.50, na.rm=TRUE)

# Calculate the 68th &50th percentile error, no bias
percentile_68_error_nobias <- quantile(abs(RiverSP_df$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_nobias <- quantile(abs(RiverSP_df$residuals_nobias), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile_nobias <- quantile(abs(RiverTile_df$residuals_nobias), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile_nobias <- quantile(abs(RiverTile_df$residuals_nobias), 0.50, na.rm=TRUE)

# Calculate the 68th & 50th percentile error
percentile_68_error <- quantile(abs(RiverSP_df$percent_diff), 0.68, na.rm=TRUE)
percentile_50_error <- quantile(abs(RiverSP_df$percent_diff), 0.50, na.rm=TRUE)

percentile_68_error_RiverTile <- quantile(abs(RiverTile_df$percent_diff), 0.68, na.rm=TRUE)
percentile_50_error_RiverTile <- quantile(abs(RiverTile_df$percent_diff), 0.50, na.rm=TRUE)

summary <- group_by(RiverTile_df, river) %>% summarise(
  count = n(),
  mean = mean(abs(residuals), na.rm = TRUE),
  sd = sd(abs(residuals), na.rm = TRUE),
  median = median(abs(residuals), na.rm = TRUE),
  IQR = IQR(abs(residuals), na.rm= TRUE),
  min =min(abs(residuals), na.rm = TRUE),
  max =max(abs(residuals), na.rm= TRUE),
  quant68 = quantile(abs(residuals), 0.68, na.rm=TRUE),
  quant68_nobais = quantile(abs(residuals_nobias), 0.68, na.rm=TRUE)
)



summary <- group_by(RiverTile_df, river) %>% summarise(
  count = n(),
  mean = mean(abs(ortho_width_m), na.rm = TRUE),
  sd = sd(abs(ortho_width_m), na.rm = TRUE),
  median = median(abs(ortho_width_m), na.rm = TRUE),
  IQR = IQR(abs(ortho_width_m), na.rm= TRUE),
  min =min(abs(ortho_width_m), na.rm = TRUE),
  max =max(abs(ortho_width_m), na.rm= TRUE),
  quant68 = quantile(abs(ortho_width_m), 0.68, na.rm=TRUE),
  quant68_nobais = quantile(abs(ortho_width_m), 0.68, na.rm=TRUE)
)


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
ggplot(RiverTile_df, aes(x = percent_diff, color = river, linetype = river)) +
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

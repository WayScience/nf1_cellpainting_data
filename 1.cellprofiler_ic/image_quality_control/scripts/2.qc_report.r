suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(platetools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(RColorBrewer))

# Paths to CSV files to generate QC report
path_to_qc_results <- file.path("../Corrected_Images/Corrected_Plate_5/IC_QC_RunImage.csv")

# Output path for bar chart
barchart_output_file <- file.path(paste0("./qc_figures/Plate_5_flagged_fov_per_well_chart.png"))

# Output path for platemap
platemap_output_file <- file.path(paste0("./qc_figures/Plate_5_platemap_flagged_fov_per_well.png"))

# Output path for site bar plot
site_FOV_output_file <- file.path(paste0("./qc_figures/Plate_5_per_site_flagged_fov.png"))

# Read in CSV files
qc_df <- read.csv(path_to_qc_results)

dim(qc_df)
head(qc_df)

# Count the number of flagged image sets
flagged_image_sets <- sum(qc_df$Image_Quality_Control_QC_Flag == 1)

well_flag_counts <- qc_df %>%
  group_by(Metadata_Well) %>%
  summarise(FlaggedFOVs = sum(Image_Quality_Control_QC_Flag))

# Print the number of image sets flagged and not processed
print(flagged_image_sets)
# Print the table
head(well_flag_counts)


# Extract the first letter from 'Metadata_Well' to create a color palette
unique_starting_letters <- unique(substr(well_flag_counts$Metadata_Well, 1, 1))

# Create a color palette based on unique starting letters
color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique_starting_letters))

# Map each unique starting letter to a color in a named vector
color_dict <- setNames(color_palette, unique_starting_letters)

# Create a new column 'Color' in the data frame based on the starting letter
well_flag_counts$Color <- color_dict[substr(well_flag_counts$Metadata_Well, 1, 1)]

# Increase the figure size to extend the chart horizontally
options(repr.plot.width=14, repr.plot.height=6)

# Create a bar chart using ggplot2 with the 'Color' column for colors
fov_chart <- ggplot(well_flag_counts, aes(x = Metadata_Well, y = FlaggedFOVs, fill = Color)) +
  geom_bar(stat = 'identity') +
  labs(x = 'Metadata_Well', y = 'Number of Flagged FOVs', title = 'Count of Flagged FOVs per Well for Plate 5') +
  scale_fill_identity() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot to qc_figures
ggsave(
        barchart_output_file,
        fov_chart,
        dpi = 500,
        height = 6,
        width = 10
    )


# Add genotype to data frame to label on plot
well_flag_counts <- well_flag_counts %>%
  mutate(
    Metadata_genotype = case_when(
      as.numeric(gsub("\\D", "", Metadata_Well)) %in% 1:4 ~ "WT",
      as.numeric(gsub("\\D", "", Metadata_Well)) %in% 5:8 ~ "HET",
      as.numeric(gsub("\\D", "", Metadata_Well)) %in% 9:12 ~ "Null",
    )
  )

fov_platemap <- platetools::raw_map(
    data = well_flag_counts$FlaggedFOVs,
    well = well_flag_counts$Metadata_Well,
    plate = 96,
    size = 8
    ) +

    ggtitle(paste("Platemap of Flagged FOV Count Per Well in Plate 5")) +
    geom_point(aes(shape = well_flag_counts$Metadata_genotype)) +
    scale_shape_discrete(name = "Genotype") +
    scale_fill_gradient(
      name = "Flagged FOV Count",
      low = "white",
      high = "red",
    ) +
    theme(
    plot.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),  # Adjust the text size
    legend.title = element_text(size = 8),  # Adjust the title size
    legend.key.size = unit(0.5, "cm")      # Adjust the key size
  )

    ggsave(
    platemap_output_file,
    fov_platemap,
    dpi = 500,
    height = 3.5,
    width = 6
    )

# Replace 'path_to_qc_results' with your actual data frame name
flagged_qc_df <- qc_df[qc_df$Image_Quality_Control_QC_Flag == 1, ]

dim(flagged_qc_df)

# Example using viridis palette
site_FOV_plot <- ggplot(flagged_qc_df, aes(x = factor(Metadata_Site))) +
  geom_bar(fill = "steelblue") +
  labs(title = "Count of Flagged FOVs per Site Integer in Plate 5",
       x = "Site",
       y = "Count")

# Save plot to qc_figures
ggsave(
        site_FOV_output_file,
        site_FOV_plot,
        dpi = 500,
        height = 6,
        width = 10
    )

# Load the parquet file into misclassed_df
misclassed_df <- arrow::read_parquet("./misclassified_cells.parquet")

dim(misclassed_df)

# Convert Metadata_Site to character in flagged_qc_df
flagged_qc_df <- flagged_qc_df %>%
  mutate(Metadata_Site = as.character(Metadata_Site))

# Convert Metadata_Site to character in misclassed_df
misclassed_df <- misclassed_df %>%
  mutate(Metadata_Site = as.character(Metadata_Site))

# Perform anti-join
filtered_misclassed_df <- misclassed_df %>%
  anti_join(flagged_qc_df, by = c("Metadata_Well", "Metadata_Site"))

dim(filtered_misclassed_df)

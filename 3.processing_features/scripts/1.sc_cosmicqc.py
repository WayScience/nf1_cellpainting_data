#!/usr/bin/env python
# coding: utf-8

# # Perform single-cell quality control
# 
# In this notebook, we perform single-cell quality control using [coSMicQC](https://github.com/WayScience/coSMicQC). 
# We filter the single cells by identifying outliers with z-scores, and use either combinations of features or one feature for each condition. For more information on the functions, please refer [**here**](https://github.com/WayScience/coSMicQC/blob/main/src/cosmicqc/analyze.py).
# 
# NOTE: As this method uses z-scoring, there is an automatic assumption that the distribution of the features are normal. If they are not, this method will be less likely to work. We have confirmed the features we are using are of normal distribution.
# 
# We use the feature(s) below to assess the technical quality of the segmented single-cells:
# 
# ### Assessing poor nuclei segmentation
# 
# To identify nuclei segmentations that include multiple nuclei, we use the following feature(s) as one condition:
# 
# - **Nuclei Solidity:** This metric quantifies how irregularly shaped a nuclei segmentation is. A value of 1 indicates that the segmentation is perfectly round and lower values indicate a very irregularly shaped nuclei (e.g., lot of indentation or protrusions).
# - **Nuclei Integrated Intensity:** This metric quantifies the total intensity of all pixels in a nuclei segmentation. Higher values mean that there are likely more pixels which can indicate multiple nuclei within the segmentation.
# 
# To identify nuclei segmentations where the respective nuclei is over-saturated or undergoing mitosis, we use the following feature(s) as two more conditions.
# 
# - **Nuclei Upper Quartile Intensity:** This metric quantifies the distribution of the pixel intensities in which 75% of pixels fall below the value (capturing the top 25% of pixel intensities). Values near 1 from this metric indicate that most pixel intensities are at a very high value.
# - **Nuclei MAD (Median Absolute Deviation) Intensity:** This metric quantifies the variability of pixel intensities within the nuclei. Higher values of this metric indicate more texture while very low values indicate very little differences between pixel intensities (e.g., low variance).

# ## Import libraries

# In[1]:


import pathlib
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
import seaborn as sns
import yaml

from cosmicqc import find_outliers
from cytodataframe import CytoDataFrame
from pycytominer import annotate


# ## Define helper function

# In[2]:


# helper function for CytoDataFrame to not change the brightness of the image
def do_not_adjust_image_brightness(image: Image.Image) -> Image.Image:
    """
    Function for CytoDataFrame to do nothing to the image and return it unchanged.

    Args:
        image (Image.Image): The input PIL Image.

    Returns:
        Image.Image: The unchanged PIL Image.
    """
    # Return the input image as is without any modifications
    return image


# ## Set paths and variables

# In[3]:


# Directory with data
data_dir = pathlib.Path("./data/converted_data/").resolve(strict=True)

# Directory to save cleaned data
cleaned_dir = pathlib.Path("./data/cleaned_profiles/")
cleaned_dir.mkdir(exist_ok=True)

# Directory to save qc figures
qc_fig_dir = pathlib.Path("./qc_figures")
qc_fig_dir.mkdir(exist_ok=True)

# metadata columns to include in output data frame
metadata_columns = [
    "Image_Metadata_Plate",
    "plate_alias",
    "Image_Metadata_Well",
    "Image_Metadata_Site",
    "Metadata_Nuclei_Location_Center_X",
    "Metadata_Nuclei_Location_Center_Y",
    "Image_FileName_DAPI",
    "Image_PathName_DAPI",
    "Nuclei_AreaShape_BoundingBoxMaximum_X",
    "Nuclei_AreaShape_BoundingBoxMaximum_Y",
    "Nuclei_AreaShape_BoundingBoxMinimum_X",
    "Nuclei_AreaShape_BoundingBoxMinimum_Y",
]

# Path to dictionary
dictionary_path = pathlib.Path("./plate_info_dictionary.yaml")


# ## Load in dictionary of plates to process

# In[4]:


# Load plate information from YAML file
with open(dictionary_path, "r") as file:
    plate_info = yaml.safe_load(file)

# Load in specific plates relevant to manuscript
plates = ["Plate_3_prime", "Plate_3", "Plate_5", "Plate_6"]

# Load in specified plates from plate_info_dictionary
plates = plate_info.keys()
dfs = {plate: pd.read_parquet(plate_info[plate]["dest_path"]) for plate in plates}

# Concatenate all dataframes into a single dataframe
combined_df = pd.concat(dfs.values(), ignore_index=True)

# Create a mapping for plate aliases
plate_alias_mapping = {
    "Plate_3": "Plate A",
    "Plate_3_prime": "Plate B",
    "Plate_5": "Plate C",
    "Plate_6": "Plate D",
}

# Add the plate_alias column
combined_df["plate_alias"] = combined_df["Image_Metadata_Plate"].map(
    plate_alias_mapping
)

# Print output
print(combined_df.shape)
combined_df.head()


# ## Get count of all plates without HET cells

# In[5]:


# Map plate names to their corresponding platemap CSV filenames
platemap_files = {
    "Plate_3": "platemap_NF1_plate3.csv",
    "Plate_3_prime": "platemap_NF1_plate3.csv",  # Plate_3_prime uses the same platemap as Plate_3
    "Plate_5": "platemap_NF1_plate5.csv",
    "Plate_6": "platemap_NF1_plate6.csv",
}

# List of plate names
plates = list(platemap_files.keys())

# Path to genotype annotation directory
genotype_annot_dir = pathlib.Path("../0.download_data/metadata")

# Annotate and filter each plate, then collect all into a list
all_annotated = []
for plate in plates:
    plate_df = dfs[plate]
    genotype_annot_path = genotype_annot_dir / platemap_files[plate]

    # Annotate with genotype info
    annotated_plate_df = annotate(
        profiles=plate_df,
        platemap=genotype_annot_path,
        join_on=["Metadata_well_position", "Image_Metadata_Well"],
    )

    # Remove HET rows
    annotated_plate_df = annotated_plate_df[
        annotated_plate_df["Metadata_genotype"] != "HET"
    ]

    all_annotated.append(annotated_plate_df)

# Combine all annotated data
total_annotated_df = pd.concat(all_annotated, ignore_index=True)

# Add Image_ prefix to match outlier DataFrames
total_annotated_df = total_annotated_df.rename(
    columns={
        "Metadata_Plate": "Image_Metadata_Plate",
        "Metadata_Well": "Image_Metadata_Well",
    }
)

# Print total number of cells (excluding HET)
print(
    f"Total number of cells across all plates (excluding HET): {total_annotated_df.shape[0]}"
)


# In[6]:


# Define the QC features
qc_features = [
    "Nuclei_Intensity_UpperQuartileIntensity_DAPI",
    "Nuclei_Intensity_IntegratedIntensity_DAPI",
    "Nuclei_Intensity_MADIntensity_DAPI",
    "Nuclei_AreaShape_Solidity",
]

# Filter combined_df to only include metadata columns and QC features
filtered_combined_df = combined_df[metadata_columns + qc_features]

filtered_combined_df_cdf = CytoDataFrame(
    data=filtered_combined_df, image_adjustment=do_not_adjust_image_brightness
)[
    [
        "Image_FileName_DAPI",
        "Nuclei_Intensity_MADIntensity_DAPI",
        "Nuclei_AreaShape_Solidity",
        "Nuclei_Intensity_UpperQuartileIntensity_DAPI",
    ]
]

# Sample 2 rows from the dataframe
filtered_combined_df_cdf.sample(n=2, random_state=0)


# ## Over-saturated nuclei (mitosis/debris)
# 
# NOTE: Threshold was determined with trial and error to find where the cutoff for good to bad quality or mitosis-ing single-cells are.

# In[7]:


# Set outlier threshold that maximizes removing most technical outliers and minimizes good cells
outlier_threshold = 2

# find nuclei with overly high intensity (over-saturated)
feature_thresholds = {
    "Nuclei_Intensity_UpperQuartileIntensity_DAPI": outlier_threshold,
}

nuclei_high_int_outliers = find_outliers(
    df=filtered_combined_df,
    metadata_columns=metadata_columns,
    feature_thresholds=feature_thresholds,
)

nuclei_high_int_outliers_cdf = CytoDataFrame(
    data=nuclei_high_int_outliers, image_adjustment=do_not_adjust_image_brightness
)[
    [
        "Nuclei_Intensity_UpperQuartileIntensity_DAPI",
        "Image_FileName_DAPI",
    ]
]

# Sort the outliers
nuclei_high_int_outliers_cdf = nuclei_high_int_outliers_cdf.sort_values(
    by="Nuclei_Intensity_UpperQuartileIntensity_DAPI", ascending=True
)

print(nuclei_high_int_outliers_cdf.shape)
nuclei_high_int_outliers_cdf.head(2)


# In[8]:


nuclei_high_int_outliers_cdf.sample(n=2, random_state=0)


# In[9]:


# Create a new column 'qc_status' in filtered_combined_df
copy_filtered_combined_df = filtered_combined_df.copy()
copy_filtered_combined_df["qc_status"] = copy_filtered_combined_df.index.map(
    lambda idx: "Failed" if idx in nuclei_high_int_outliers.index else "Passed"
)

# Set figure size
plt.figure(figsize=(6, 6))

# Histogram plot for nuclei high intensity outliers, grouped by qc_status
sns.histplot(
    data=copy_filtered_combined_df,
    x="Nuclei_Intensity_UpperQuartileIntensity_DAPI",
    hue="qc_status",  # Color bins by qc_status
    multiple="layer",  # Stack the bars for better comparison
    bins=30,  # Adjust the number of bins as needed
    kde=False,  # Disable KDE
    palette={
        "Passed": "#006400",
        "Failed": "#CC0099",
    },  # Custom colors for Passed and Failed
    legend=False,
)

# Customize the plot
plt.xlabel("Nuclei upper quartile intensity (DAPI)", fontsize=14)
plt.ylabel("Count", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(
    title="QC Status",
    loc="upper right",
    labels=["Failed", "Passed"],
    fontsize=12,
    title_fontsize=14,
)

plt.savefig(
    qc_fig_dir / "high_intensity_nuclei_outliers_histogram.png",
    dpi=500,
    bbox_inches="tight",
)

# Show the plot
plt.show()


# In[10]:


# Print out the number of outliers across plates
outlier_counts = nuclei_high_int_outliers["Image_Metadata_Plate"].value_counts()

# Calculate the percentage of outliers
total_counts = combined_df["Image_Metadata_Plate"].value_counts()
outlier_percentages = (outlier_counts / total_counts) * 100

# Print the counts and percentages
for plate, count in outlier_counts.items():
    print(f"{plate}: {count} outliers ({outlier_percentages[plate]:.2f}%)")


# In[11]:


# Set outlier threshold that maximizes removing most technical outliers and minimizes good cells
outlier_threshold = 2

# find nuclei with overly high intensity (over-saturated)
feature_thresholds = {
    "Nuclei_Intensity_MADIntensity_DAPI": outlier_threshold,
}

blurry_nuclei_outliers = find_outliers(
    df=filtered_combined_df,
    metadata_columns=metadata_columns,
    feature_thresholds=feature_thresholds,
)

blurry_nuclei_outliers_cdf = CytoDataFrame(
    data=blurry_nuclei_outliers, image_adjustment=do_not_adjust_image_brightness
)[
    [
        "Nuclei_Intensity_MADIntensity_DAPI",
        "Image_FileName_DAPI",
    ]
]

# Sort the outliers
blurry_nuclei_outliers_cdf = blurry_nuclei_outliers_cdf.sort_values(
    by="Nuclei_Intensity_MADIntensity_DAPI", ascending=True
)

print(blurry_nuclei_outliers_cdf.shape)
blurry_nuclei_outliers_cdf.head(2)


# In[12]:


blurry_nuclei_outliers_cdf.sample(n=2, random_state=0)


# In[13]:


# Create a new column 'qc_status' in filtered_combined_df
copy_filtered_combined_df = filtered_combined_df.copy()
copy_filtered_combined_df["qc_status"] = copy_filtered_combined_df.index.map(
    lambda idx: "Failed" if idx in blurry_nuclei_outliers.index else "Passed"
)

# Set figure size
plt.figure(figsize=(6, 6))

# Histogram plot for blurry nuclei outliers, grouped by qc_status
sns.histplot(
    data=copy_filtered_combined_df,
    x="Nuclei_Intensity_MADIntensity_DAPI",
    hue="qc_status",  # Color bins by qc_status
    multiple="layer",  # Stack the bars for better comparison
    bins=30,  # Adjust the number of bins as needed
    kde=False,  # Disable KDE
    palette={
        "Passed": "#006400",
        "Failed": "#CC0099",
    },  # Custom colors for Passed and Failed
    legend=False,
)

# Customize the plot
plt.xlabel("Nuclei MAD intensity (DAPI)", fontsize=14)
plt.ylabel("Count", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(
    title="QC Status",
    loc="upper right",
    labels=["Passed", "Failed"],
    fontsize=12,
    title_fontsize=14,
)

plt.savefig(
    qc_fig_dir / "blurry_nuclei_outliers_histogram.png",
    dpi=500,
    bbox_inches="tight",
)

# Show the plot
plt.show()


# In[14]:


# Print out the number of outliers across plates
outlier_counts = blurry_nuclei_outliers["Image_Metadata_Plate"].value_counts()

# Calculate the percentage of outliers
total_counts = combined_df["Image_Metadata_Plate"].value_counts()
outlier_percentages = (outlier_counts / total_counts) * 100

# Print the counts and percentages
for plate, count in outlier_counts.items():
    print(f"{plate}: {count} outliers ({outlier_percentages[plate]:.2f}%)")


# ## Over-segmented nuclei (reflected by irregular, non-circular shape)
# 
# NOTE: Threshold was determined with trial and error to find where the cutoff for good to bad quality single-cell are.

# In[15]:


# find irregular shaped nuclei
feature_thresholds = {
    # outlier threshold for only solidity (values less than the mean = poor quality)
    "Nuclei_AreaShape_Solidity": -1.25,
    # outlier threshold for only intensity (values above than the mean = highly intense)
    "Nuclei_Intensity_IntegratedIntensity_DAPI": 2,
}

irregular_nuclei_outliers = find_outliers(
    df=filtered_combined_df,
    metadata_columns=metadata_columns,
    feature_thresholds=feature_thresholds,
)

irregular_nuclei_outliers_cdf = CytoDataFrame(
    data=irregular_nuclei_outliers, image_adjustment=do_not_adjust_image_brightness
)[
    [
        "Nuclei_AreaShape_Solidity",
        "Nuclei_Intensity_IntegratedIntensity_DAPI",
        "Image_Metadata_Plate",
        "Image_Metadata_Well",
        "Image_Metadata_Site",
        "Image_FileName_DAPI",
    ]
]

print(irregular_nuclei_outliers_cdf.shape)
irregular_nuclei_outliers_cdf.sort_values(
    by="Nuclei_AreaShape_Solidity", ascending=False
).head(2)


# In[16]:


irregular_nuclei_outliers_cdf.sample(n=2, random_state=0)


# In[17]:


# Create a new column 'qc_status' in filtered_combined_df
copy_filtered_combined_df = filtered_combined_df.copy()
copy_filtered_combined_df["qc_status"] = copy_filtered_combined_df.index.map(
    lambda idx: "Failed" if idx in irregular_nuclei_outliers.index else "Passed"
)

# Set figure size
plt.figure(figsize=(12, 6))

# Scatterplot for intensity and solidity in the DataFrame
ax = sns.scatterplot(
    data=copy_filtered_combined_df,
    x="Nuclei_AreaShape_Solidity",
    y="Nuclei_Intensity_IntegratedIntensity_DAPI",
    hue="qc_status",  # Color points by qc_status
    alpha=0.5,  # Set transparency for better visibility
    palette={
        "Passed": "#006400",
        "Failed": "#CC0099",
    },  # Custom colors for Passed and Failed
)

# Customize the plot
plt.xlabel("Nuclei solidity", fontsize=14)
plt.ylabel("Nuclei sum intensity (DAPI)", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Retrieve handles and labels correctly
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles=handles,
    labels=labels,
    title="QC Status",
    loc="upper left",
    fontsize=12,
    title_fontsize=14,
)

plt.savefig(
    qc_fig_dir / "irregular_nuclei_outliers_scatterplot.png",
    dpi=500,
    bbox_inches="tight",
)

# Show the plot
plt.show()


# In[18]:


# Print out the number of outliers across plates
outlier_counts = irregular_nuclei_outliers["Image_Metadata_Plate"].value_counts()

# Calculate the percentage of outliers
total_counts = combined_df["Image_Metadata_Plate"].value_counts()
outlier_percentages = (outlier_counts / total_counts) * 100

# Print the counts and percentages
for plate, count in outlier_counts.items():
    print(f"{plate}: {count} outliers ({outlier_percentages[plate]:.2f}%)")


# ### Get the count of outliers across all plates for each condition

# In[19]:


# Filtering out the HET cells for calculation
nuclei_high_int_outliers_filtered = nuclei_high_int_outliers[
    ~(
        (
            (nuclei_high_int_outliers["Image_Metadata_Plate"] == "Plate_6")
            & (
                nuclei_high_int_outliers["Image_Metadata_Well"].str.startswith(
                    ("C", "F")
                )
            )
        )
        | (
            (
                nuclei_high_int_outliers["Image_Metadata_Plate"].isin(
                    ["Plate_3", "Plate_3_prime", "Plate_5"]
                )
            )
            & (
                nuclei_high_int_outliers["Image_Metadata_Well"].str.endswith(
                    ("05", "06", "07", "08")
                )
            )
        )
    )
]

blurry_nuclei_outliers_filtered = blurry_nuclei_outliers[
    ~(
        (
            (blurry_nuclei_outliers["Image_Metadata_Plate"] == "Plate_6")
            & (blurry_nuclei_outliers["Image_Metadata_Well"].str.startswith(("C", "F")))
        )
        | (
            (
                blurry_nuclei_outliers["Image_Metadata_Plate"].isin(
                    ["Plate_3", "Plate_3_prime", "Plate_5"]
                )
            )
            & (
                blurry_nuclei_outliers["Image_Metadata_Well"].str.endswith(
                    ("05", "06", "07", "08")
                )
            )
        )
    )
]

# Define the unique cell ID columns
merge_cols = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Site",
    "Metadata_Nuclei_Location_Center_X",
    "Metadata_Nuclei_Location_Center_Y",
]

# Combine both filtered outlier DataFrames and drop duplicates
combined_outliers_filtered = pd.concat(
    [nuclei_high_int_outliers_filtered, blurry_nuclei_outliers_filtered]
).drop_duplicates(subset=merge_cols)

# Merge with total_annotated_df to flag outliers
is_combined_outlier = (
    total_annotated_df.merge(
        combined_outliers_filtered[merge_cols].assign(is_combined_outlier=True),
        on=merge_cols,
        how="left",
    )["is_combined_outlier"]
    .astype("boolean")
    .fillna(False)
)

# Count and print
num_failed_combined = is_combined_outlier.sum()
percent_failed_combined = (num_failed_combined / total_annotated_df.shape[0]) * 100
print(
    f"Combined high intensity and blurry nuclei outliers: {num_failed_combined} cells failed ({percent_failed_combined:.2f}% failed)"
)


# In[20]:


# Filter irregular nuclei outliers
irregular_nuclei_outliers_filtered = irregular_nuclei_outliers[
    ~(
        (
            (irregular_nuclei_outliers["Image_Metadata_Plate"] == "Plate_6")
            & (
                irregular_nuclei_outliers["Image_Metadata_Well"].str.startswith(
                    ("C", "F")
                )
            )
        )
        | (
            (
                irregular_nuclei_outliers["Image_Metadata_Plate"].isin(
                    ["Plate_3", "Plate_3_prime", "Plate_5"]
                )
            )
            & (
                irregular_nuclei_outliers["Image_Metadata_Well"].str.endswith(
                    ("05", "06", "07", "08")
                )
            )
        )
    )
]

# Define ID columns
merge_cols = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Site",
    "Metadata_Nuclei_Location_Center_X",
    "Metadata_Nuclei_Location_Center_Y",
]

# Combine blurry and high intensity, drop duplicates
previous_outliers = pd.concat(
    [nuclei_high_int_outliers_filtered, blurry_nuclei_outliers_filtered]
).drop_duplicates(subset=merge_cols)

# Remove any overlapping cells from the irregular set
irregular_cleaned = (
    irregular_nuclei_outliers_filtered.merge(
        previous_outliers[merge_cols], on=merge_cols, how="left", indicator=True
    )
    .query("_merge == 'left_only'")
    .drop(columns="_merge")
)

# Then count as before
unique_irregular_cells = (
    total_annotated_df.merge(
        irregular_cleaned[merge_cols].assign(is_irregular=True),
        on=merge_cols,
        how="left",
    )["is_irregular"]
    .astype("boolean")
    .fillna(False)
)

# Print counts
num_irregular = unique_irregular_cells.sum()
percent_irregular = (num_irregular / total_annotated_df.shape[0]) * 100
print(
    f"Unique irregular nuclei outliers (excluding blurry/high intensity): {num_irregular} ({percent_irregular:.2f}%)"
)


# ## Save cleaned profiles (filtered to remove outliers) per plate

# In[21]:


# Collect the indices of the outliers
outlier_indices = pd.Index(
    pd.concat(
        [nuclei_high_int_outliers, irregular_nuclei_outliers, blurry_nuclei_outliers]
    ).index.unique()
)

# Remove rows with outlier indices from combined_df
combined_df_cleaned = combined_df.drop(outlier_indices)

# Save cleaned data for each plate and update the dictionary with cleaned paths
for plate in plates:
    plate_df = combined_df[
        combined_df["Image_Metadata_Plate"] == plate
    ]  # Original plate data
    plate_df_cleaned = combined_df_cleaned[
        combined_df_cleaned["Image_Metadata_Plate"] == plate
    ]  # Cleaned plate data

    # Calculate number of failed cells (rows removed)
    failed_cells = plate_df.shape[0] - plate_df_cleaned.shape[0]

    # Calculate percentage of failed cells
    failed_percentage = (failed_cells / plate_df.shape[0]) * 100

    # Print the number of failed cells, percentage, and remaining single cells
    print(f"{plate}: {failed_cells} cells failed ({failed_percentage:.2f}% failed)")
    print(f"{plate}: {plate_df_cleaned.shape[0]} single cells in cleaned dataset")

    # Clean the plate data
    plate_df_cleaned = plate_df_cleaned.drop(
        columns=["plate_alias"]
    )  # Remove plate_alias column

    # Save cleaned data for each plate
    cleaned_path = f"{cleaned_dir}/{plate}_cleaned.parquet"
    plate_df_cleaned.to_parquet(cleaned_path)
    plate_info[plate]["cleaned_path"] = cleaned_path
    print(plate, ":", plate_df_cleaned.shape)


# ## Determine count of the number of outliers detected across all plates
# 
# This count excludes HET cells.

# In[22]:


# Define columns that uniquely identify each cell
merge_cols = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Site",
    "Metadata_Nuclei_Location_Center_X",
    "Metadata_Nuclei_Location_Center_Y",
]

# Combine outlier DataFrames (already have Image_ prefixes)
combined_outliers = pd.concat(
    [
        nuclei_high_int_outliers_filtered,
        irregular_nuclei_outliers_filtered,
        blurry_nuclei_outliers_filtered,
    ]
)

# Drop duplicates just in case (some cells may appear in more than one category)
combined_outliers = combined_outliers.drop_duplicates(subset=merge_cols)

# Merge to flag outliers in the annotated data
is_outlier_col = total_annotated_df.merge(
    combined_outliers[merge_cols].assign(is_outlier=True), on=merge_cols, how="left"
)["is_outlier"].astype("boolean")
total_annotated_df["is_outlier"] = is_outlier_col.fillna(False)

# Count and percent
num_outliers = total_annotated_df["is_outlier"].sum()
percent_outliers = (num_outliers / total_annotated_df.shape[0]) * 100
print(
    f"Total outliers detected across all plates (total = {total_annotated_df.shape[0]}): {num_outliers} ({percent_outliers:.2f}% of all cells)"
)


# In[23]:


# Filter for Plates 3, 3_prime, and 5
plates_to_check = ["Plate_3", "Plate_3_prime", "Plate_5"]

filtered_df = total_annotated_df[
    total_annotated_df["Image_Metadata_Plate"].isin(plates_to_check)
]

# Count and percent for these plates
num_outliers_plates_3_5 = filtered_df["is_outlier"].sum()
percent_outliers_plates_3_5 = (num_outliers_plates_3_5 / filtered_df.shape[0]) * 100

print(
    f"Outliers detected in Plates 3, 3_prime, and 5 (total = {filtered_df.shape[0]}): {num_outliers_plates_3_5} ({percent_outliers_plates_3_5:.2f}% of these plates' cells)"
)


# In[24]:


# Filter for Plate_6
plate_6_df = total_annotated_df[total_annotated_df["Image_Metadata_Plate"] == "Plate_6"]

# Count and percent for Plate_6
num_outliers_plate_6 = plate_6_df["is_outlier"].sum()
percent_outliers_plate_6 = (num_outliers_plate_6 / plate_6_df.shape[0]) * 100

print(
    f"Outliers detected in Plate 6 (total = {plate_6_df.shape[0]}): {num_outliers_plate_6} ({percent_outliers_plate_6:.2f}% of Plate 6 cells)"
)


# ## Dump the new cleaned path to the dictionary for downstream processing

# In[25]:


with open(dictionary_path, "w") as file:
    yaml.dump(plate_info, file)


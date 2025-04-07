##################################################################################
#                                                                                #
#             Parkinson's Disease Clinical Study II synthetic data               #  
#                                                                                #
##################################################################################

# Set the working directory to the location of your project files
setwd("~/Documents/Parkinson")

# --------------------- Run Required Previous Scripts --------------------------
# Source previous scripts that load data and perform preliminary processing.
# These scripts load clinical metadata, create a shifted matrix, and process data.
source("clinical_metadata_v4_2023/DataReading_v4.R")
source("data_processing.R") 
# ------------------------- 1. Collect and Filter Data ---------------------------
# Merge clinical_data and metadata by participant_id using a left join. 
# This keeps only the clinical information for patients who have corresponding records.
master_matrix <- merge(clinical_data, metadata, by = "participant_id", all.x = TRUE)
dim(master_matrix)  # Check dimensions to inspect the number of records

# Reorder and select the variables used in the analysis, including metadata.
master_matrix <- master_matrix %>% select(
  participant_id, visit_month, study, diagnosis_at_baseline, age_at_baseline, sex, race, 
  case_control_other_at_baseline, global_famhistory,
  on_levodopa, on_dopamine_agonist, on_other_pd_medications, 
  has_known_GBA_mutation_in_WGS, has_known_LRRK2_mutation_in_WGS,
  has_known_SNCA_mutation_in_WGS, has_known_APOE_E4_mutation_in_WGS, 
  has_known_PD_Mutation_in_WGS,
  # test_abeta, test_ptau, test_tau, guid, study_participant_id  # Uncomment if necessary
  mds_updrs_part_i_summary_score, mds_updrs_part_ii_summary_score, 
  mds_updrs_part_iii_summary_score, mds_updrs_part_iv_summary_score, 
  code_upd2hy_hoehn_and_yahr_stage, moca_total_score,
  pdq39_mobility_score, pdq39_adl_score, pdq39_emotional_score, pdq39_stigma_score, 
  pdq39_social_score, pdq39_cognition_score, pdq39_communication_score,
  pdq39_discomfort_score, Schwad_ADL_score, REM_score, ess_sleepiness_score, 
  upsit_total_score
)

# ------------------------- Change Node Names ----------------------------------
# Define new node names for clarity
node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", 
                "age_at_baseline", "sex", "race", "case_control_other_at_baseline", 
                "global_famhistory", "GBA_mut", "LRRK2_mut", "SNCA_mut", "APOE_E4_mut", 
                "PD_Mut", "levodopa", "dopamine", "other_med", "UPDRS1", "UPDRS2", 
                "UPDRS3", "UPDRS4", "HOEHN_Stage", "MOCA", "Mobility39", "ADL39",
                "Emotional39", "Stigma39", "Social39", "Cognition39", "Communication39",
                "Discomfort39", "Schwad_ADL", "REM", "ESS", "UPSIT")
# Rename the columns in master_matrix using the new node names
colnames(master_matrix) <- node_names

# ----------------- Induce Community Structure for the Network ------------------
# Define a community structure for grouping variables:
# - "Demographics": first 9 variables
# - "Mutations": next 5 variables
# - "Medication": following 3 variables
# - "General PD Severity": next 5 variables
# - "Cognitive": 1 variable
# - "Disability": next 9 variables
# - "Sleep": 2 variables
# - "Smell": 1 variable
community_structure <- c(
  rep("Demographics", 9),
  rep("Mutations", 5),
  rep("Medication", 3),
  rep("General PD Severity", 5),
  rep("Cognitive", 1),
  rep("Disability", 9),
  rep("Sleep", 2),
  rep("Smell", 1)
)
structure <- data.frame(node_names, community_structure)

# --------------------- Filter by NA Values ------------------------------------
# Convert master_matrix to a data frame and filter out rows with more than 7 missing values.
multiple_visit <- as.data.frame(master_matrix)
multiple_visit <- multiple_visit %>% filter(rowSums(is.na(.)) <= 7)
dim(multiple_visit)  # Check the dimensions (e.g., 4096 patients)

# Calculate the percentage of missing values for each variable.
porcentaje_nas_variables <- as.data.frame(colMeans(is.na(multiple_visit)))

# Define a threshold (umbral) for missing values (70% in this example).
umbral <- 0.7
# Select columns with missing percentage below the threshold.
columnas_filtradas <- colnames(multiple_visit)[porcentaje_nas_variables < umbral]
print(columnas_filtradas)  # Display the selected variables

# Retain only the filtered columns in the final dataset.
multiple_visit <- multiple_visit[, columnas_filtradas]

# Apply an additional mutation: for male patients with UPDRS3 > 10, set UPSIT to 40.
multiple_visit <- multiple_visit %>%
  mutate(UPSIT = ifelse(sex == "Male" & UPDRS3 > 10, 40, UPSIT))
# Update the structure data frame to keep only filtered node names.
structure <- subset(structure, structure$node_names %in% columnas_filtradas)

# ------------------- List Variables to Scale ----------------------------------
# Define the variables that will be scaled for network analysis.
variables_to_scale <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39",
                        "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39",
                        "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT")

# ----------------- Network Analysis using MariNET and EBICglasso ---------------
# Load the MariNET package.
library("MariNET")

# Run a linear mixed model analysis (MariNET) including sex as a covariate.
model_1 <- lmm_analysis(multiple_visit, variables_to_scale, random_effects = "(1|sex/participant_id)")
# Run a MariNET analysis without including sex as a covariate.
model <- lmm_analysis(multiple_visit, variables_to_scale)

# Estimate a network using the EBICglasso method with Spearman correlation.
network <- estimateNetwork(multiple_visit[variables_to_scale],
                           default = "EBICglasso",
                           tuning = 0.25,
                           threshold = FALSE,
                           corMethod = "spearman")

# Compare the models:
# - difference_models_1: Difference between model_1 (with covariate) and model (without covariate).
# - difference_models_2: Difference between model_1 and the EBICglasso network.
difference_models_1 <- differentiation(model_1, model)
difference_models_2 <- differentiation(model_1, network$graph)

# ----------------- Extract Community Structure for Plotting ---------------------
# Obtain the community structure for the variables to scale from the structure dataset.
structure_plot <- structure %>%
  filter(node_names %in% variables_to_scale) %>%
  pull(community_structure)

# ----------------------- Define Custom Labels and Layout -------------------------
# Define custom labels for the network graph.
custom_labels <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mob39",
                   "ADL39", "Emot39", "Stig39", "Social39", "Cogn39",
                   "Commu39", "Discomf39", "Schw_ADL", "ESS", "UPSIT")

# Define a custom layout matrix for plotting multiple graphs.
# Layout matrix: 3 rows x 2 columns.
# Graphs will be arranged as follows:
# 1: Top-left, 2: Top-right, 3: Middle-left, 4: Middle-right, 5: Bottom-left, 6: Bottom-right (Legend or additional plot)
layout_matrix <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2, byrow = TRUE)

# Uncomment the next line to save the plot as a PNG file.
# png("third_use_case_full.png", width = 800, height = 1000, res = 150)  # Adjust dimensions and resolution as needed

# Set up the custom layout for plotting.
layout(layout_matrix)

# Adjust graphical parameters: set smaller margins and larger title font.
par(mar = c(0, 0, 0, 0))  # Set margins to zero
par(cex.main = 1)         # Increase title font size

# ------------------------- Plot 1: MariNET with Sex as Covariate -------------------
qgraph(model_1,
       title = "(A) MariNET using sex as covariate",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen", "lightblue", "orange", "pink", "grey"),
       labels = custom_labels,
       legend = FALSE,
       vsize = 11,
       label.font = 1.3,
       label.cex = 1.4)

# ----------------------------- Plot 2: Legend --------------------------------------
# Create an empty plot area for the legend.
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# Add a legend using the unique community structure labels and defined colors.
legend("center",
       legend = unique(structure_plot),
       fill = c("lightgreen", "lightblue", "orange", "pink", "grey"),
       title = "",
       cex = 1.1,   # Increase legend font size
       bty = "n")   # Remove legend border

# ------------------ Plot 3: MariNET without Sex as Covariate -------------------------
qgraph(model,
       title = "(B) MariNET not including sex as covariate",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen", "lightblue", "orange", "pink", "grey"),
       labels = custom_labels,
       legend = FALSE,
       vsize = 11,
       label.font = 1.3,
       label.cex = 1.4)

# ----------------- Plot 4: Difference between (A) and (B) ----------------------------
qgraph(difference_models_1,
       title = "(C) Difference between (A) and (B)",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen", "lightblue", "orange", "pink", "grey"),
       labels = custom_labels,
       legend = FALSE,
       vsize = 11,
       label.font = 1.3,
       label.cex = 1.4)

# ----------------- Plot 5: EBICGlasso Network -----------------------------------------
qgraph(network$graph,
       title = "(D) EBICGlasso",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen", "lightblue", "orange", "pink", "grey"),
       labels = custom_labels,
       legend = FALSE,
       vsize = 11,
       label.font = 1.3,
       label.cex = 1.4)

# -------------- Plot 6: Difference between (A) and (D) --------------------------------
qgraph(difference_models_2,
       title = "(E) Difference between (A) and (D)",
       layout = "spring",
       groups = structure_plot,
       color = c("lightgreen", "lightblue", "orange", "pink", "grey"),
       labels = custom_labels,
       legend = FALSE,
       vsize = 11,
       label.font = 1.3,
       label.cex = 1.4)

# Reset the layout to the default single plot layout.
layout(1)

# If saving the plot as PNG, close the device.
dev.off()

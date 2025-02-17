################################################################################

#Automatized script for longitudinal analysis by linear mixed models

#Marina Vargas-Fernández


library(dplyr)      # filtering datasets
library(rstatix)    # summary statistics
library(ggpubr)     # convenient summary statistics and plots
library(tidyverse)  # not sure if usefull
library(lme4)       # linear mixed effects model
library(qgraph)     # network package

setwd("~/Documents/Parkinson")

#Run required previous scripts

source("clinical_metadata_v4_2023/DataReading_v4.R")
source("shifted_matrix.R")
source("data_processing.R") #previous NetworkComparison_3.R
#source("Long_analysis_v4.R")

#1. Collect and filter data 

# View(metadata)

master_matrix <- merge(clinical_data, metadata,by = "participant_id", all.x = T) 
# left join, only need clinical information on patients who have registers 
dim(master_matrix)

# All variables used in the analysis, including metadata information, reordered
master_matrix <- master_matrix %>% select(participant_id, visit_month, study, diagnosis_at_baseline, age_at_baseline, sex, race, case_control_other_at_baseline, global_famhistory
                                          ,on_levodopa, on_dopamine_agonist, on_other_pd_medications, has_known_GBA_mutation_in_WGS, has_known_LRRK2_mutation_in_WGS
                                          ,has_known_SNCA_mutation_in_WGS, has_known_APOE_E4_mutation_in_WGS, has_known_PD_Mutation_in_WGS
                                          #, test_abeta, test_ptau, test_tau, guid, study_participant_id
                                          ,mds_updrs_part_i_summary_score, mds_updrs_part_ii_summary_score, mds_updrs_part_iii_summary_score, mds_updrs_part_iv_summary_score, code_upd2hy_hoehn_and_yahr_stage, moca_total_score
                                          ,pdq39_mobility_score, pdq39_adl_score, pdq39_emotional_score, pdq39_stigma_score, pdq39_social_score, pdq39_cognition_score, pdq39_communication_score
                                          ,pdq39_discomfort_score, Schwad_ADL_score, REM_score, ess_sleepiness_score, upsit_total_score)
# Change node names
node_names <- c("participant_id", "visit_month", "study", "diagnosis_at_baseline", "age_at_baseline", "sex", "race", "case_control_other_at_baseline", "global_famhistory", "GBA_mut", "LRRK2_mut"
                ,"SNCA_mut", "APOE_E4_mut", "PD_Mut", "levodopa", "dopamine", "other_med", "UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "HOEHN_Stage", "MOCA", "Mobility39", "ADL39","Emotional39", "Stigma39", "Social39", "Cognition39","Communication39","Discomfort39", "Schwad_ADL", "REM", "ESS", "UPSIT")
colnames(master_matrix) <- node_names

# Induce community structure, and create dataset which contains this information
community_structure <- c(rep("Demographics", 9),rep("Mutations", 5),rep("Medication", 3), rep("General PD Severity", 5), rep("Cognitive", 1) , rep("Disability", 9), rep("Sleep", 2), rep("Smell",1) )
structure <- data.frame(node_names, community_structure)


# Filter by NA (7 max)
multiple_visit <- as.data.frame(master_matrix)
multiple_visit <- multiple_visit %>%
  filter(rowSums(is.na(.)) <= 7)
#test for unfiltered: multiple_visit<-master_matrix
dim(multiple_visit) #4096 patients

# Percentage of NA is calculated by variable
porcentaje_nas_variables <- as.data.frame(colMeans(is.na(multiple_visit)))

# And later filtered under selected value
umbral <- 0.7
columnas_filtradas <- colnames(multiple_visit)[porcentaje_nas_variables < umbral]
# Show selected variables
print(columnas_filtradas)

# Only filtered variables are included on final dataset
multiple_visit <- multiple_visit[, columnas_filtradas]

structure <- subset(structure, structure$node_names %in% columnas_filtradas)


#2. Scale selected variables and order by ID

# Ordered by participant id
long_visit <- multiple_visit %>%
  arrange(participant_id) %>%
  group_by(participant_id)

# Scale to z-score continuous variables
# List of variables to scale
variables_to_scale <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA", "Mobility39", 
                        "ADL39", "Emotional39", "Stigma39", "Social39", "Cognition39", 
                        "Communication39", "Discomfort39", "Schwad_ADL", "ESS", "UPSIT")

# Loop through each variable and scale it
for (variable in variables_to_scale) {
  long_visit[[variable]] <- scale(long_visit[[variable]])
}

#Check, patients are not filtered by idiopathic PD, all of them are included
table(long_visit$diagnosis_at_baseline) 
long_visit <- as.data.frame(long_visit)
 
#write.table(long_visit, file = "clinical_database.csv", col.names = T, row.names = F, sep = ",")

#3. Build multiple Linear Mixed Model 

#Select variables to include on network
dependent_variables <- variables_to_scale 
models <- list()

variables_string <- paste(variables_to_scale, collapse = " + ")
# Loop through each dependent variable
for (dependent_variable in dependent_variables) {
  # Create formula for the model
  formula <-as.formula(paste(dependent_variable, "~ ", variables_string, "+ (1|participant_id)"))
  
  # Create and fit the model
  model <- lmer(formula, data = long_visit)
  
  # Print the model summary
  #print(summary(model))
  
  # Optionally, you can store the model in a list or perform further analysis
  models[[dependent_variable]] <- summary(model)$coefficients[,3]
}

# Convert list of lists into matrix
score<- do.call(rbind, models)
colnames(score)[1] <- variables_to_scale[1]

# Change format of score matrix using external functions 
source("shifted_matrix.R")

score_matrix <- score_matrix(score, shift_matrix(score))

# Get node degree
grado <- rowSums(score_matrix) + colSums(score_matrix)

# Normalize degree
grado_normalizado <- grado / max(grado)

node_size <- 5 + 3.5 * grado_normalizado  # Adjust normalization to size

# Get community structure
structure_plot <- structure %>%
  filter(node_names %in% variables_to_scale) %>%
  pull(community_structure)

# Plot adjusted node size
qgraph(score_matrix, vsize = node_size, groups = structure_plot,  layout="spring", color=c("lightgreen", "lightblue","orange", "pink","grey") ) 

# Addition ( M + t(M) )
symm_score <- (score_matrix + t(score_matrix) )/2

# Plot new network (symmetrical)
qgraph(symm_score, vsize = node_size, groups = structure_plot,  layout="spring", color=c("lightgreen", "lightblue","orange", "pink","grey") ) 

#To save image!!
#setwd("~/Documents/Parkinson/images_HD/")

# Open a PNG graphics device
#png("qgraph_plot_2.png", width = 1800, height = 1000, res = 500)

# Generate the qgraph plot
qgraph(symm_score, vsize = node_size, groups = structure_plot,  layout="spring", color=c("lightgreen", "lightblue","orange", "pink","grey"), 
       legend = T)#, legend.cex = 0.2) 

# Close the graphics device
#dev.off()

###############################################################################################
################Test to see if network estimation produces same result#########################
library(bootnet)

network <- estimateNetwork(long_visit[variables_to_scale], default = "EBICglasso", tuning = 0.25, threshold = FALSE, corMethod ="spearman")
filtered_structure <- structure[structure$node_names %in% variables_to_scale, ]

qgraph(network$graph, groups = filtered_structure$community_structure,  layout="spring", color=c("lightgreen", "lightblue","orange", "pink","grey") ) 

###############################################################################################
################Test to see if mlVAR produces same result######################################
library(mlVAR)

fit1 <- mlVAR(long_visit, variables_to_scale, "participant_id", estimator = "lmer")
#Summary of all parameter estimates:
  summary(fit1)
# Compare temporal relationships:
#layout(t(1:2))
plot(fit1, "temporal", title = "Estimated temporal relationships", layout = "spring")
# Compare contemporaneous partial correlations:
plot(fit1, "contemporaneous", title = "Estimated contemporaneous relationships",
     layout = "spring")
# Compare between-subjects partial correlations:
plot(fit1, "between", title = "Estimated between-subjects relationships",
     layout = "spring")



###############################################################################################
#check with mutations
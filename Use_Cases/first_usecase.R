##################################################################################
#                                                                                #
#                      COVID-19 student mental health study                      #  
#                                                                                #
##################################################################################

# ----------------- 1. Loading Packages & Data ------------------------------

# Load required libraries:
library(lme4)       # For fitting linear mixed-effects models
library(dplyr)      # For data manipulation (filtering, selecting, etc.)
library(tidyr)      # For reshaping/tidying data
library(qgraph)     # For network visualization
library(bootnet)    # For bootstrapping network estimates
library(mlVAR)      # For multi-level vector autoregression (mlVAR) models
library(MariNET)    # For MariNET analysis (specific network analysis method)

# Load pre-cleaned dataset (assumes that the file clean_network.RData contains a data object named Data2)
# Data for this use case can be obtained from https://osf.io/erp7v.
#load("~/Documents/Data/clean_network.RData")

# Define the original variable names to investigate in the dataset
vars <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q11","Q12","Q13","Q14","Q15","Q16","Q17","Q18")

# Copy the dataset for further manipulation
Data5b <- Data2

# Create new, more descriptive labels for the variables
varLabs <- c("Relax","Irritable","Worry","Nervous","Future","Anhedonia",
             "Tired","Hungry","Alone","Angry","Social_offline","Social_online","Music",
             "Procrastinate","Outdoors","C19_occupied","C19_worry","Home")

# Rename the columns in Data5b that match the names in 'vars' to the new labels in 'varLabs'
names(Data5b)[names(Data5b) %in% vars] <- varLabs

# Remove items that are not needed for the analysis:
# Here, "Hungry", "Angry", "Music", and "Procrastinate" are excluded.
Data5b <- Data5b %>% select(-Hungry, -Angry, -Music, -Procrastinate)
# Update the labels to reflect the removed variables.
varLabs <- varLabs[!varLabs %in% c("Hungry", "Angry", "Music", "Procrastinate")]

# Assign the cleaned dataset to 'master_matrix' for use in later analyses
master_matrix <- Data5b

# --------------- 3. Fit mlVAR Network Models: Orthogonal Estimation -------------

# Run the multi-level VAR model using the mlVAR package.
# - 'vars': the list of variables to include in the network (using updated labels)
# - 'idvar': identifier for subjects
# - 'dayvar': variable representing the day (repeated measures)
# - 'beepvar': variable representing measurement occasions within a day
# - 'lags': number of time lags (1 in this case)
# - 'temporal' & 'contemporaneous': estimation method ("orthogonal")
# - 'nCores': number of cores for parallel processing to speed up computation
res <- mlVAR(master_matrix,   
             vars = varLabs, 
             idvar = "id",
             dayvar = "day", 
             beepvar = "beep", 
             lags = 1,
             temporal = "orthogonal", 
             contemporaneous = "orthogonal",
             nCores = 8)

# Uncomment the line below if you want to save the mlVAR result object to file for future use:
# save(res, file = "network_orthogonal.RData")

# Reset master_matrix to Data5b (if any changes were made) 
master_matrix <- Data5b

# Create a community structure for the network.
# In this example, the first 7 variables are labeled "Stress", the next 3 "Social", and the final 4 "COVID-19".
community_structure <- c(rep("Stress", 7), rep("Social", 3), rep("COVID-19", 4))
structure <- data.frame(varLabs, community_structure)

# Perform a MariNET analysis using lmm_analysis (from the MariNET package)
# This function likely computes a symmetrical score based on a linear mixed model specified by random effects "(1|id/day)"
symm_score <- lmm_analysis(master_matrix, varLabs, random_effects = "(1|id/day)")

# Estimate a network using the EBICglasso method via the estimateNetwork function.
# - This uses a subset of columns (2:15) from master_matrix
# - 'tuning' parameter is set to 0.25
# - 'corMethod' is set to "spearman" for Spearman correlation
network <- estimateNetwork(master_matrix[, c(2:15)], 
                           default = "EBICglasso", 
                           tuning = 0.25, 
                           threshold = FALSE, 
                           corMethod = "spearman")

# Define node labels for plotting the network graphs
names <- c("Relax","Irritable","Worry","Nervous","Future",
           "Anhedonia","Tired","Alone",
           "Social-offline", "Social-online", "Outdoors",
           "C19-occupied", "C19-worry", "Home")

# Define grouping of nodes (communities) for visualization purposes
gr <- list('Stress' = c(1:7), 'Social' = c(8:10), 'COVID-19' = c(11:14))

# Extract the contemporaneous network from the mlVAR results.
# 'layout' sets the node arrangement, 'nonsig' determines how non-significant edges are handled,
# and 'rule' "and" indicates that an edge is retained only if significant in both directions.
cont <- getNet(res, "contemporaneous", layout = "spring", nonsig = "show", rule = "and")

# Define colors for nodes corresponding to their community membership.
node_colors <- c(rep("#A6CEE3", 7), rep("#FDBF6F", 3), rep("#B2DF8A", 4))

# (Optional) Extract community structure for plotting purposes (if needed for later customizations)
structure_plot <- structure %>%
  filter(varLabs %in% varLabs) %>%  # Here, adjust filter condition if necessary
  pull(community_structure)

# (Optional) To save high resolution png
png("first_use_case.png", 
    width = 2400,       
    height = 1800,       
    res = 300,
    type = "cairo",
    bg = "transparent")  
# Adjust graphical layout: define a 2x2 grid to display three graphs and a legend
layout(matrix(c(1, 2,
                3, 4), byrow = TRUE, ncol = 2), 
       widths = c(1, 1), heights = c(1, 1))

# -------------------- 1. First Graph (Top-Left): MariNET -----------------------
edge_labels <- ifelse(abs(symm_score) > 16, round(symm_score, 2), NA)
edge_label_colors <- ifelse(edge_labels > 0, "blue",
                            ifelse(edge_labels < 0, "red", "black"))


graph1 <- qgraph(symm_score,
                 title = "(A) MariNET",
                 layout = "spring",
                 theme = 'colorblind',    # Use a colorblind-friendly palette
                 negDashed = FALSE,       # Do not dash negative edges
                 edge.labels = edge_labels,      # To add label weight
                 edge.label.cex = 1,    # Label weight size
                 edge.label.font = 2,         # negrita
                 edge.label.color = edge_label_colors,
                 groups = gr,             # Group nodes based on pre-defined community structure
                 nodeNames = names,       # Node labels for the graph
                 legend = FALSE,          # Legend is plotted separately
                 color = node_colors,     # Use pre-defined node colors
                 vsize = 7,              # Node size
                 label.cex = 1.2)           

# ---------------------- 2. Legend (Top-Right) ----------------------------------
# Create an empty plot space for the legend
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

# Add legend at the center of the empty plot
legend("center", legend = unique(structure$community_structure), fill = unique(node_colors), 
       title = "",        # No title for the legend
       cex = 1.3,         # Increase font size of legend text
       bty = "n")         # Remove border around legend

# ------------------ 3. Second Graph (Bottom-Left): mlVAR ------------------------
edge_labels <- ifelse(abs(cont) > 0.3, round(cont, 2), NA)
edge_label_colors <- ifelse(edge_labels > 0, "blue",
                            ifelse(edge_labels < 0, "red", "black"))
graph2 <- qgraph(cont, 
                 layout = "spring",
                 title = "(B) MlVAR", 
                 theme = 'colorblind', 
                 negDashed = FALSE, 
                 edge.labels = edge_labels,     
                 edge.label.cex = 1,   
                 edge.label.font = 2,         # negrita
                 edge.label.color = edge_label_colors,
                 groups = gr, 
                 legend = FALSE, 
                 nodeNames = names,
                 color = node_colors, 
                 vsize = 7, 
                 label.cex = 1.2)           

# ----------------- 4. Third Graph (Bottom-Right): EBICGlasso ---------------------
edge_labels <- ifelse(abs(network$graph) > 0.3, round(network$graph, 2), NA)
edge_label_colors <- ifelse(edge_labels > 0, "blue",
                            ifelse(edge_labels < 0, "red", "black"))
graph3 <- qgraph(network$graph,
                 title = "(C) EBICGlasso",
                 layout = "spring",
                 theme = 'colorblind', 
                 negDashed = FALSE,
                 edge.labels = edge_labels,     
                 edge.label.cex = 1,        # tamaño del texto
                 edge.label.font = 2,         # negrita
                 edge.label.color = edge_label_colors,
                 groups = gr, 
                 nodeNames = names,
                 legend = FALSE,
                 color = node_colors, 
                 vsize = 7, 
                 label.cex = 1.2)           

# Reset the plotting layout back to the default (single plot layout)
dev.off()
layout(1)


# Comparison metrics 

W_mari   <- normalization(symm_score)
W_mlvar  <- cont
W_glasso <- network$graph


metrics <- compare_three_networks(W_mari, W_mlvar, W_glasso)
    
print(metrics$PairwiseMetrics)
  
print(metrics$TopologicalMetrics)
    

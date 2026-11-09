
# ============================================================
# MICE PROTEIN EXPRESSION PROJECT
# 01 - DATA EXPLORATION
# ============================================================


# ------------------------------------------------------------
# 1. INSTALL AND LOAD PACKAGES
# ------------------------------------------------------------

# Run install.packages() only once if you don't have them.
# install.packages(c("readxl", "dplyr", "ggplot2"))

library(readxl)
library(dplyr)
library(ggplot2)


# ------------------------------------------------------------
# 2. LOAD THE DATASET
# ------------------------------------------------------------

mice_data <- read_excel("data/Data_Cortex_Nuclear.xls")

# Display the first 6 rows
head(mice_data)

# Display number of rows and columns
dim(mice_data)

# Display column names
names(mice_data)

# Display the structure and data types
str(mice_data)


# ------------------------------------------------------------
# 3. UNDERSTAND THE DATASET
# ------------------------------------------------------------

# Number of observations
nrow(mice_data)

# Number of variables
ncol(mice_data)

# Number of quantitative variables
sum(sapply(mice_data, is.numeric))

# Number of qualitative variables
sum(sapply(mice_data, function(x) is.character(x) ||
             is.factor(x)))


# ------------------------------------------------------------
# 4. INVESTIGATE MISSING VALUES
# ------------------------------------------------------------

# Count missing values in each column
missing_values <- colSums(is.na(mice_data))

# Display only columns containing missing values
missing_values[missing_values > 0]

# Total number of missing values in the dataset
sum(is.na(mice_data))

# Percentage of missing values per column
missing_percentage <- colMeans(is.na(mice_data)) * 100

# Display columns with missing values, sorted from highest
# to lowest percentage
sort(missing_percentage[missing_percentage > 0],
     decreasing = TRUE)


# ------------------------------------------------------------
# 5. UNDERSTAND THE TARGET VARIABLE: CLASS
# ------------------------------------------------------------

# Number of observations in each class
class_distribution <- table(mice_data$class)

class_distribution

# Percentage of observations in each class
round(prop.table(class_distribution) * 100, 2)

# Visualize class distribution
barplot(
  class_distribution,
  main = "Distribution of Mouse Classes",
  xlab = "Class",
  ylab = "Number of observations",
  col = "steelblue",
  las = 2
)


# ------------------------------------------------------------
# 6. UNDERSTAND MOUSE IDENTIFIERS
# ------------------------------------------------------------

# Extract the mouse ID before the underscore
mouse_base_id <- sub("_.*", "", mice_data$MouseID)

# Add the extracted ID to the dataset
mice_data$mouse_base_id <- mouse_base_id

# Number of distinct mice
length(unique(mouse_base_id))

# Number of observations belonging to each mouse
observations_per_mouse <- table(mouse_base_id)

# How many mice have each number of observations?
table(observations_per_mouse)


# ------------------------------------------------------------
# 7. CHECK WHETHER EACH MOUSE HAS ONE CLASS
# ------------------------------------------------------------

mouse_class_check <- mice_data %>%
  group_by(mouse_base_id) %>%
  summarise(
    number_of_classes = n_distinct(class),
    .groups = "drop"
  )

# Display the number of mice with each number of classes
table(mouse_class_check$number_of_classes)


# ------------------------------------------------------------
# 8. UNDERSTAND THE EXPERIMENTAL VARIABLES
# ------------------------------------------------------------

# Genotype distribution
table(mice_data$Genotype)

# Treatment distribution
table(mice_data$Treatment)

# Behavior distribution
table(mice_data$Behavior)

# Target class distribution by experimental variables
table(mice_data$Genotype, mice_data$class)
table(mice_data$Treatment, mice_data$class)
table(mice_data$Behavior, mice_data$class)


# ------------------------------------------------------------
# 9. SUMMARY OF PROTEIN VARIABLES
# ------------------------------------------------------------

# Select numeric columns (protein measurements)
protein_data <- mice_data %>%
  select(where(is.numeric))

# Summary statistics for each protein
summary(protein_data)

# Mean of each protein
sapply(protein_data, mean, na.rm = TRUE)

# Standard deviation of each protein
sapply(protein_data, sd, na.rm = TRUE)


# ------------------------------------------------------------
# 10. FIRST LOOK AT PROTEIN DISTRIBUTIONS
# ------------------------------------------------------------

# Distribution of one protein: DYRK1A_N
hist(
  mice_data$DYRK1A_N,
  main = "Distribution of DYRK1A_N",
  xlab = "Protein expression",
  col = "lightblue",
  breaks = 30
)


# ------------------------------------------------------------
# END OF INITIAL DATA EXPLORATION
# ------------------------------------------------------------


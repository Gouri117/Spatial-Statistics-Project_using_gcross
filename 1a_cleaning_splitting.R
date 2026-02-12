################################
# Code to rename Tumor 1-7 into a general category called "Tumor"
# The dataset is from /nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc//nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc
library(dplyr)
library(stringr)
library(readr)

data <- read.csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/subset_data_samples_with_outcome.csv", header = TRUE)

data <- data %>%
  mutate(CELL_TYPE = ifelse(str_detect(CELL_TYPE, "Tumor"), "Tumor", CELL_TYPE))

write.csv(data, "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/subset_data_samples_with_outcome_tumor_category.csv")

# Step 1: Get unique sample IDs in order
unique_samples <- unique(data$sample)

# Step 2: Split into two sets of 98
split_1_samples <- unique_samples[1:98]
split_2_samples <- unique_samples[99:196]

# Step 3: Filter the data
data_split1 <- data %>% filter(sample %in% split_1_samples)
data_split2 <- data %>% filter(sample %in% split_2_samples)

# Step 4: Save to CSV
write.csv(data_split1, "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/subset_data_samples_with_outcome_split1.csv")
write.csv(data_split2, "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/subset_data_samples_with_outcome_split2.csv")

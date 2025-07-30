# Load necessary libraries
library(dplyr)
library(tidyr)

data = get(load("/home/mohammad/Desktop/New_folder_10_08_2024/fin/data/count_vivo.RData"))

filtered_data <- subset(data, grepl("II", treat) & !grepl("III", treat))
filtered_data  <- filtered_data[!(grepl("Pre", filtered_data$Sample)), ]
filtered_data  <-filtered_data[!(grepl("t0", filtered_data$Sample)), ]
data1=filtered_data[,c(1,3,7)]
# Create a unique numbering within each treatment group using group_indices()
data1 <- data1 %>%
  mutate(group_id = group_indices(., treat, Sample)) %>%  # Ensures unique group numbering
  mutate(treat_group = paste0(treat, "_", group_id))  # Create unique group label

# Reshape data: Create a wide format while handling missing values
matrix_data <- data1 %>%
  group_by(treat_group) %>%
  mutate(row_num = row_number()) %>%  # Assigns row numbers within each group
  pivot_wider(names_from = row_num, values_from = Count, names_prefix = "values_") %>%
  ungroup()

# Ensure unique row names
rownames(matrix_data) <- make.names(matrix_data$treat_group, unique = TRUE)  # Unique row names
matrix_data <- matrix_data %>% select(-treat_group)  # Remove treat_group column

# Convert dataframe to matrix
matrix_result <- as.matrix(matrix_data)

# Print the final matrix
print(matrix_result, quote = FALSE)
countt = matrix_result[,c(4:24)]

countt <- apply(countt, 2, as.numeric)

# Ensure it remains a matrix (apply can return a dataframe sometimes)
countt <- as.matrix(numeric_matrix)

colnames(countt)= c("Hom.25", "Hom.26", "Hom.27", 
                                "PTEN.15", "PTEN.16", "PTEN.17", 
                                "TSC2.18", "TSC2.19", "TSC2.20", 
                                "A829P.21", "A829P.22", "A829P.28", 
                                "D816A.10", "D816A.11", "D816A.13", 
                                "D816E.23", "D816E.24", "D816E.36", 
                                "V654A.29", "V654A.30", "V654A.39")

my_vector <- c("DMSO_1", "DMSO_2", "DMSO_3", "DMSO_4", "DMSO_5", 
               "Imatinib_1", "Imatinib_2", "Imatinib_3", "Imatinib_4", "Imatinib_5", 
               "Sunitinib_1", "Sunitinib_2", "Sunitinib_3", "Sunitinib_4")

rownames(countt)= my_vector

vivo_II_count=countt
vivo_II_count_ad <- sapply(seq(1, ncol(vivo_II_count), by=3), function(i) rowSums(vivo_II_count[, i:(i+2)]))
new_vec <- unique(sub("\\.\\d+$", "", colnames(vivo_II_count)))
colnames(vivo_II_count_ad)=new_vec
save(vivo_II_count,file = "/home/mohammad/Desktop/data_gist_paper/in_vivo/vivo_II_count.RData")
save(vivo_II_count_ad,file = "/home/mohammad/Desktop/data_gist_paper/in_vivo/vivo_II_count_ad.RData")




#--- volume II
data = get(load("/home/mohammad/Desktop/New_folder_10_08_2024/fin/data/volume_vivo.RData"))

filtered_data <- subset(data, grepl("II", treat) & !grepl("III", treat))
filtered_data  <-filtered_data[!(grepl("t0", filtered_data$treat)), ]
filtered_data$treat = rep(c("DMSO","Imatinib","Sunitinib"),times=c(5,5,4))
vivo_II_volume = filtered_data
save(vivo_II_volume,file = "/home/mohammad/Desktop/data_gist_paper/in_vivo/vivo_II_volume.RData")

#---------------------------     
#---------------------------          VIVO III
#---------------------------     



# Load necessary libraries
library(dplyr)
library(tidyr)

data = get(load("/home/mohammad/Desktop/New_folder_10_08_2024/fin/data/count_vivo.RData"))
filtered_data <- subset(data, grepl("III", treat))  # Ensures it ends with "III"

filtered_data  <- filtered_data[!(grepl("Pre", filtered_data$Sample)), ]
filtered_data  <-filtered_data[!(grepl("t0", filtered_data$Sample)), ]

data1=filtered_data[,c(1,3,7)]
# Create a unique numbering within each treatment group using group_indices()
data1 <- data1 %>%
  mutate(group_id = group_indices(., treat, Sample)) %>%  # Ensures unique group numbering
  mutate(treat_group = paste0(treat, "_", group_id))  # Create unique group label

# Reshape data: Create a wide format while handling missing values
matrix_data <- data1 %>%
  group_by(treat_group) %>%
  mutate(row_num = row_number()) %>%  # Assigns row numbers within each group
  pivot_wider(names_from = row_num, values_from = Count, names_prefix = "values_") %>%
  ungroup()

# Ensure unique row names
rownames(matrix_data) <- make.names(matrix_data$treat_group, unique = TRUE)  # Unique row names
matrix_data <- matrix_data %>% select(-treat_group)  # Remove treat_group column
matrix_data$Sample
# Convert dataframe to matrix
matrix_result <- as.matrix(matrix_data)

# Print the final matrix
print(matrix_result, quote = FALSE)
countt = matrix_result[,c(4:18)]

countt <- apply(countt, 2, as.numeric)

# Ensure it remains a matrix (apply can return a dataframe sometimes)
#countt <- as.matrix(numeric_matrix)

colnames(countt)= c("T1parental.4","T1parental.6","T1parental.9",
                    "A829P.21", "A829P.22", "A829P.28", 
                    "D816A.10", "D816A.11", "D816A.13", 
                    "D816E.23", "D816E.24", "D816E.36", 
                    "V654A.29", "V654A.30", "V654A.39")

my_vector <- c("DMSO_1", "DMSO_2", "DMSO_3", "DMSO_4", 
               "Imatinib_1", "Imatinib_2", "Imatinib_3", "Imatinib_4", "Imatinib_5",  "Imatinib_5",
               "Sunitinib_1", "Sunitinib_2", "Sunitinib_3", "Sunitinib_4", "Sunitinib_5")

rownames(countt)= my_vector

vivo_III_count=countt
vivo_III_count_ad <- sapply(seq(1, ncol(vivo_III_count), by=3), function(i) rowSums(vivo_III_count[, i:(i+2)]))
new_vec <- unique(sub("\\.\\d+$", "", colnames(vivo_III_count)))
colnames(vivo_III_count_ad)=new_vec
save(vivo_III_count,file = "/home/mohammad/Desktop/data_gist_paper/in_vivo/vivo_III_count.RData")
save(vivo_III_count_ad,file = "/home/mohammad/Desktop/data_gist_paper/in_vivo/vivo_III_count_ad.RData")




#--- volume II
data = get(load("/home/mohammad/Desktop/New_folder_10_08_2024/fin/data/volume_vivo.RData"))

filtered_data <- subset(data, grepl("III", treat))  # Ensures it ends with "III"

filtered_data  <-filtered_data[!(grepl("t0", filtered_data$treat)), ]
filtered_data$treat = rep(c("DMSO","Imatinib","Sunitinib"),times=c(4,6,5))
vivo_III_volume = filtered_data
save(vivo_III_volume,file = "/home/mohammad/Desktop/data_gist_paper/in_vivo/vivo_III_volume.RData")

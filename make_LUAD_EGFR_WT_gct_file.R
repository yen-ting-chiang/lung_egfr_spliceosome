RowData <- read.csv("LUAD_batch_corrected_version_removed.csv",
                    check.names = FALSE)
RowData <- RowData %>% 
  select(!c(1:2))

ColData <- read.csv("join_data_selected.csv")
ColData <- ColData %>% 
  select(!c(1))

library(dplyr)
RowData_t <- as.data.frame(t(RowData))
library(tibble)
RowData_t <- tibble::rownames_to_column(RowData_t, "VALUE")
joined_data <- full_join(RowData_t, 
                         ColData, 
                         by = c("VALUE"="Sample"))

joined_data_filtered <- joined_data %>% 
  filter(IMPACT == "WT")
joined_data_filtered_t <- as.data.frame(t(joined_data_filtered))
joined_data_filtered_t <-na.omit(joined_data_filtered_t)
write.csv(joined_data_filtered_t, 
         file = "LUAD_batch_corrected_EGFR_WT.csv")

joined_data_filtered <- joined_data %>% 
  filter(IMPACT == "MODERATE")
joined_data_filtered_t <- as.data.frame(t(joined_data_filtered))
joined_data_filtered_t <-na.omit(joined_data_filtered_t)
write.csv(joined_data_filtered_t, 
          file = "LUAD_batch_corrected_EGFR_MODERATE.csv")


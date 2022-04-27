RowData <- read.csv("LUAD_batch_corrected_version_removed.csv",
                    check.names = FALSE)
RowData <- RowData %>% 
  select(!c(1:2))

ColData <- read.csv("join_data_selected_exon19del_L858R.csv")
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
  filter(exon19del.729.761._L858R == "p.L858R")
joined_data_filtered_t <- as.data.frame(t(joined_data_filtered))
joined_data_filtered_t <-na.omit(joined_data_filtered_t)
write.csv(joined_data_filtered_t, 
          file = "LUAD_batch_corrected_EGFR_L858R.csv")

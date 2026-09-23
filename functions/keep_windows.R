#libraries
library(tidyverse)
library(here)


### ID large gaps (> 5 days) in the data to remove after regularization ####
keep_windows <- function(sp_df){

  #ID gaps greater than 5 days
  all_keep_df <- data.frame()

  for(i in 1:length(unique(sp_df$id))){

  curr_id = unique(sp_df$id)[i]
  windows_df <- sp_df %>%
    filter(id == curr_id) %>%
    arrange(date) %>%
    mutate(
      # Calculate time diff in days
      gap = as.numeric(difftime(date, lag(date, default = first(date)), units = "days")),
      # Increment burst ID whenever a gap > threshold occurs
      burst_id = cumsum(gap >= 5)) %>%
    ungroup()
    
  keep_windows <- windows_df %>%
    group_by(id, burst_id) %>%
    summarise(start_time = min(date), 
              end_time = max(date), 
              .groups = "drop")
   
  all_keep_df <- rbind(all_keep_df, keep_windows)
    
  } #end tag id loop
  
  return(all_keep_df)

} #end function

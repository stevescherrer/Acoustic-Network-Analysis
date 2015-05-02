#TODO: stripchart for detections?
# 57451 detected before tagged?
# which tag has > 150000 detections?
# are all BRFA in outs PB tags?
# station swapping for greg

# merge_all
# Rebuild dataframe from which boxplots are made

# Comparing Bottomfish tag returns to Acoustic Web App Model
# Written 14 January 2015 by Stephen Scherrer

#### Clearning Workspace and setting directory ---------
rm(list=ls()) # Clear workspace
setwd('/Users/stephenscherrer/Documents/Work/UH/Projects/dissertation work/Spacial Ecology/Acoustic-Network-Analysis/Output Files/')

# Import any principle dependencies----------
# install.packages('wesanderson') # color palett for plotting
# install.packages('matlab')
# install.packages('maps')
# install.packages('mapdata')
# install.packages('maptools')
# install.packages('scales')
# install.packages('ggmap')
library("matlab", lib.loc)
library('maps')
library('mapdata')
library('maptools')  # for shapefiles
library('scales')  # for transparency
library('ggmap')
source('/Users/stephenscherrer/Documents/Programming/R/utility_functions.R')
library('reshape') # merge_all



## Functions-----------------------

#### Functions-------------------------------------------------------------------

### experiment_dates
experiment_dates = function(vue_data, bottomfish_tag_ids = FALSE){
  ## Function takes a vue_data dataframe and returns the start date, end date and
  ## elapsed time of study. If provided a list of tag IDs, study starts on 
  ## first detection from a tag, and ends on last detection
  ## Arguments: 
  ## vue_data: A dataframe from a vue export called with function load_vemco
  ## bottomfish_tag_ids: an optional list of specific tags for analysis
  ## Returns:
  ## Dataframe of 3 clolumns
  temp_data = clean_vue(vue_data, bottomfish_tag_ids)
  date_range = as.POSIXct(range(strftime(temp_data$datetime, format = "%Y-%m-%d")))
  elapsed_dates = date_range[2] - date_range[1]
  dates = as.data.frame(cbind(as.character(date_range[1]), as.character(date_range[2]), elapsed_dates))
  colnames(dates)[1] = 'Start Date'
  colnames(dates)[2] = 'End Date'
  colnames(dates)[3] = 'Length of Study'
  return (dates)
}

### Adding/Adjusting study dates

adjust_vue_study_dates = function(vue_data, tagging_data, bottomfish_tag_ids = FALSE){
  ## Function turns POSIXct date into a single number indicating the number of
  ## days of the study that have elapsed at the time of a detection
  ## Arguments:
  ## vue_data: a data frame from a vue export, called with function load_vemco
  ## tagging_data: a dataframe of tagging metadata, called with function 
  ##load_tagging_data
  ## bottomfish_tag_ids = an optional list of specific tags to analyze
  ## Returns:
  ## Dataframe similar to vue_data but with additional study date column
  if(bottomfish_tag_ids[1] != FALSE){
    vue_data = vue_data[vue_data$tag_id %in% bottomfish_tag_ids, ]
    tagging_data = tagging_data[tagging_data$vem_tag_id %in% bottomfish_tag_ids, ]
  }
  earliest_date = min(tagging_data$datetime)
  vue_data$study_date = as.numeric(difftime(vue_data$datetime, earliest_date), units = 'days')
  #tagging_data$study_date = as.numeric(difftime(tagging_data$datetime, earliest_date, na.rm = TRUE), units = 'days')
  return(vue_data)
}

adjust_tagging_study_dates = function(vue_data, tagging_data, bottomfish_tag_ids = FALSE){
  ## Function turns POSIXct date into a single number indicating the number of
  ## days of the study that have elapsed at the time of a detection
  ## Arguments:
  ## vue_data: a data frame from a vue export, called with function load_vemco
  ## tagging_data: a dataframe of tagging metadata, called with function 
  ##load_tagging_data
  ## bottomfish_tag_ids = an optional list of specific tags to analyze
  ## Returns:
  ## Dataframe similar to tagging_data but with additional study date column
  if(bottomfish_tag_ids[1] != FALSE){
    vue_data = vue_data[vue_data$tag_id %in% bottomfish_tag_ids, ]
    tagging_data = tagging_data[tagging_data$vem_tag_id %in% bottomfish_tag_ids, ]
  }
  earliest_date = min(tagging_data$datetime)
  #vue_data$study_date = as.numeric(difftime(vue_data$datetime, earliest_date, na.rm = TRUE), units = 'days')
  tagging_data$study_date = as.numeric(difftime(tagging_data$datetime, earliest_date), units = 'days')
  return(tagging_data)
}


adjust_receiver_study_dates = function(receiver_data, tagging_data, bottomfish_tag_ids = FALSE){
  ## Function turns POSIXct date into a single number indicating the number of
  ## days of the study that have elapsed at the time of a detection
  ## Arguments:
  ## receiver_data: a data frame from a receiver export, called with function load_receiver
  ## tagging_data: a dataframe of tagging metadata, called with function 
  ##load_tagging_data
  ## bottomfish_tag_ids = an optional list of specific tags to analyze
  ## Returns:
  ## Dataframe similar to tagging_data but with additional study date column
  if(bottomfish_tag_ids[1] != FALSE){
    #vue_data = vue_data[vue_data$tag_id %in% bottomfish_tag_ids, ]
    tagging_data = tagging_data[tagging_data$vem_tag_id %in% bottomfish_tag_ids, ]
  }
  earliest_date = min(tagging_data$datetime)
  #vue_data$study_date = as.numeric(difftime(vue_data$datetime, earliest_date, na.rm = TRUE), units = 'days')
  receiver_data$deployment_study_date = as.numeric(difftime(receiver_data$deployment_date, earliest_date), units = 'days')
  receiver_data$recovery_study_date = as.numeric(difftime(receiver_data$recovery_date, earliest_date), units = 'days')
  return(receiver_data)
}

### tagging_date
tagging_date = function(tagging_data, bottomfish_tag_ids = FALSE){
  ## Function to retrieve date transmitters were deployed
  ## Arguments:
  ## tagging_data: a data dataframe of tag ids and tagging POSIXct datetimes
  ## called by function load_tagging_data
  ## bottomfish_tag_ids: an optonal list of specific tags to analyze
  ## Returns:
  ## dataframe with three columns. first column is tag id, second is date 
  ## tag was deployed, third is study date of tag deployment
  tagging_data = tagging_data[which(is.na(tagging_data$vem_tag_id) == 0), ]
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(tagging_data$vem_tag_id))[tagging_data$vem_tag_id])}
  tagging_date = matrix(0,length(bottomfish_tag_ids), 3)
  tagging_date[ ,1] = as.numeric(as.character(bottomfish_tag_ids))
  for (i in 1:length(bottomfish_tag_ids)){
    tagging_date[i,2] = as.character(tagging_data$datetime[tagging_data$vem_tag_id == bottomfish_tag_ids[i]])
    tagging_date[i,3] = as.numeric(tagging_data$study_date[tagging_data$vem_tag_id == bottomfish_tag_ids[i]])
  }
  tag_date_out = as.data.frame(tagging_date)
  colnames(tag_date_out)[1] =  'tag_id'
  colnames(tag_date_out)[2] =  'tagging_date'
  colnames(tag_date_out)[3] =  'study_date'
  return (tag_date_out)
}


first_transmission = function(vue_data, bottomfish_tag_ids = FALSE){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  first_transmissions = matrix(0, length(bottomfish_tag_ids), 3)
  first_transmissions[ ,1] = bottomfish_tag_ids
  for (i in 1:length(bottomfish_tag_ids)){
    indv = vue_data[vue_data$tag_id == bottomfish_tag_ids[i], ]
    first_transmissions[i,2] = as.character(min(indv$datetime[as.character(indv$station) != 'Tagging Location']))
    first_transmissions[i,3] = as.numeric(min(indv$study_date[as.character(indv$station) != 'Tagging Location']))
  }
  first_transmissions = as.data.frame(first_transmissions)
  colnames(first_transmissions) = c('tag_id', 'first_transmission_datetime', 'first_transmission_study_date')
  return (first_transmissions)
}


last_transmission = function(vue_data, bottomfish_tag_ids = FALSE){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  last_transmissions = matrix(0, length(bottomfish_tag_ids), 3)
  last_transmissions[ ,1] = bottomfish_tag_ids
  for (i in 1:length(bottomfish_tag_ids)){
    indv = vue_data[vue_data$tag_id == bottomfish_tag_ids[i], ]
    last_transmissions[i,2] = as.character(max(indv$datetime))
    last_transmissions[i,3] = as.numeric(max(indv$study_date))
  }
  last_transmissions = as.data.frame(last_transmissions)
  colnames(last_transmissions) = c('tag_id', 'last_transmission_datetime', 'last_transmission_study_date')
  return (last_transmissions)
}

time_between_first_transmission_and_tagging_date = function(vue_data, tagging_data, bottomfish_tag_ids){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  time_diff = as.numeric(as.character(first_transmission(vue_data, bottomfish_tag_ids)$first_transmission_study_date)) - as.numeric(as.character(tagging_date(tagging_data, bottomfish_tag_ids)$study_date))
  time_between_first_and_tag_date = as.data.frame(cbind(bottomfish_tag_ids, time_diff))
  colnames(time_between_first_and_tag_date) = c('tag_id', 'time_diff_between_tagging_and_first_detection')
  return(time_between_first_and_tag_date)
}

time_at_liberty = function(vue_data, tagging_data, bottomfish_tag_ids){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  time_diff = as.numeric(as.character(last_transmission(vue_data, bottomfish_tag_ids)$last_transmission_study_date)) - as.numeric(as.character(tagging_date(tagging_data, bottomfish_tag_ids)$study_date))
  time_at_lib = as.data.frame(cbind(bottomfish_tag_ids, time_diff))
  colnames(time_at_lib) = c('tag_id', 'time_at_liberty')
  return(time_at_lib)
}

time_beteween_first_last_transmission = function(vue_data, bottomfish_tag_ids = FALSE){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  time_between_first_last_transmission = matrix(0,length(bottomfish_tag_ids), 2)
  time_between_first_last_transmission[ ,1] = bottomfish_tag_ids
  time_between_first_last_transmission[ ,2] = (as.numeric(as.character(last_transmission(vue_data, bottomfish_tag_ids)$last_transmission_study_date)) - 
                                                 as.numeric(as.character(first_transmission(vue_data, bottomfish_tag_ids)$first_transmission_study_date)))
  #   for (i in 1:length(bottomfish_tag_ids)){
  #     indv_data = vue_data[vue_data$tag_id == bottomfish_tag_ids[i], ]
  #     time_between_first_last_transmission[i,2] = as.numeric(as.character(max(indv_data$study_date))) - )
  #   }
  time_between_first_last = as.data.frame(time_between_first_last_transmission)
  colnames(time_between_first_last) = c('tag_id', 'time_between_first_last_transmission')
  return(time_between_first_last)
}




transmission_stats = function(vue_data, tagging_data, bottomfish_tag_ids = FALSE){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  t_stats_mat = matrix(0, length(bottomfish_tag_ids), 6)
  # tag_id
  t_stats_mat[ ,1] = as.numeric(as.character(bottomfish_tag_ids))
  # study date of tagging
  t_stats_mat[ ,2] = ceiling(as.numeric(as.character(tagging_date(tagging_data, bottomfish_tag_ids)$study_date)))
  # study date of first transmission
  t_stats_mat[ ,3] = ceiling(as.numeric(as.character(first_transmission(vue_data, bottomfish_tag_ids)$first_transmission_study_date)))
  # study date of last transmission
  t_stats_mat[ ,4] = ceiling(as.numeric(as.character(last_transmission(vue_data, bottomfish_tag_ids)$last_transmission_study_date)))
  # time at liberty
  t_stats_mat[ ,5] = ceiling(as.numeric(as.character(time_at_liberty(vue_data, tagging_data, bottomfish_tag_ids)$time_at_liberty)))
  # time between first and last detection
  t_stats_mat[ ,6] = ceiling(as.numeric(as.character(time_beteween_first_last_transmission(vue_data, bottomfish_tag_ids)$time_between_first_last_transmission)))
  indv_date_stats = as.data.frame(t_stats_mat)
  colnames(indv_date_stats)[1] =  'tag_id'
  colnames(indv_date_stats)[2] =  'tagging_study_date'
  colnames(indv_date_stats)[3] =  'first_transmission'  
  colnames(indv_date_stats)[4] =  'last_transmission'
  colnames(indv_date_stats)[5] =  'time_at_liberty'
  colnames(indv_date_stats)[6] =  'time_between_first_and_last_transmissions'
  return(indv_date_stats)
}

unique_days_detected = function(vue_data, bottomfish_tag_ids = FALSE){
  ## Returns the number of unique days each fish was detected on the array
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  unique_days_detected = matrix(0,length(bottomfish_tag_ids),1)
  for (i in 1:length(bottomfish_tag_ids)){
    indv_data = vue_data[vue_data$tag_id == bottomfish_tag_ids[i],]
    unique_days_detected[i] = length(unique(floor(indv_data$study_date)))
  }
  print(fivenum(unique_days_detected))
  
  pdf('Unique_days_detected_boxplot.pdf')
  boxplot(unique_days_detected[1:length(bottomfish_tag_ids)], main = 'Unique Days Detected', ylab = 'Days' )
  dev.off()
  
  unique_detection_days = as.data.frame(cbind(bottomfish_tag_ids, unique_days_detected))
  colnames(unique_detection_days)[1] =  'tag_id'
  colnames(unique_detection_days)[2] =  'unique_days_detected'
  return(unique_detection_days)
}

## Number of detections
number_of_detections = function(vue_data, bottomfish_tag_ids = FALSE){
  ## Determines the number of detections for each tag on the array. Also determines 
  # the percentage of all detections those individuals account for. 
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  number_of_detections = matrix(0,length(bottomfish_tag_ids),1)
  for (a in 1:length(bottomfish_tag_ids)){
    number_of_detections[a] = length(vue_data$tag_id[vue_data$tag_id == bottomfish_tag_ids[a]])
  }
  percentage_of_detections = (number_of_detections / length(vue_data$datetime))*100
  #print(length(vue_data))
  fivenum(number_of_detections)
  #pdf('Total Number of Detections Boxplot.pdf')
  boxplot(number_of_detections[1:length(bottomfish_tag_ids)], main = 'Total Number of Detections', ylab = 'Number of Detections' )
  #dev.off()
  number_and_percent_detections = as.data.frame(cbind(bottomfish_tag_ids, number_of_detections, percentage_of_detections))
  colnames(number_and_percent_detections)[1] = 'tag_id'
  colnames(number_and_percent_detections)[2] = '#_of_detections'
  colnames(number_and_percent_detections)[3] = '%_of_all_detections'
  return(number_and_percent_detections)
}

detections_per_day = function(vue_data, tagging_data, bottomfish_tag_ids = FALSE){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  tol = floor(transmission_stats(vue_data, tagging_data, bottomfish_tag_ids)$time_at_liberty)+1
  detections_per_day = number_of_detections(vue_data, bottomfish_tag_ids)[,2]/
    tol
  detections_per_day_out = as.data.frame(cbind(bottomfish_tag_ids, detections_per_day))
  colnames(detections_per_day_out)[1] = 'tag_id'
  colnames(detections_per_day_out)[2] = 'detections / day'
  return(detections_per_day_out)
}



count_detections_per_day = function(vue_data, bottomfish_tag_ids = FALSE){
  ## Function returns a vector that is the same length as the number of total days
  ## in the study. Each bin represents a day and the number is the total number
  ## of transmissions detected for that day.
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  counting_detections = rep(0, max(ceiling(vue_data$study_date)))
  for (i in 1:length(vue_data$study_date)){
    counting_detections[floor(vue_data$study_date[i])] = counting_detections[floor(vue_data$study_date[i])]+1
  }
  return(counting_detections)
}

plot_transmission_frequency = function(vue_data, bottomfish_tag_ids = FALSE){
  
  ### Function to plot transmission recovered per day for all tags and individual
  ## tags
  
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  # counting_detections = rep(0, max(ceiling(vue_data$study_date))) 
  
  pdf('Detections Per Day For All Tags Histogram.pdf')
  par(mfrow = c(1,1), oma=c(0,0,2,0))
  hist(floor(vue_data[,1]), xlim = c(0, max(vue_data[,1])), breaks = c(0, max(vue_data[,1])+1), xlab = 'Study Date', ylab = 'Detection Frequency', main = 'Detections per Day for All Tags')
  dev.off
  
  for (tag_id in bottomfish_tag_ids){
    indv_data = vue_data[vue_data$tag_id == tag_id, ]
    pdf(paste(tag_id,' Detection Per Day Histogram.pdf'))
    par(mfrow = c(1,1), oma=c(0,0,2,0))
    hist(floor(indv_data$study_date),  
         breaks = c(max(indv_data$study_date)+1), 
         xlim = c(0, max(vue_data$study_date)), xlab = 'Study Date', 
         ylab = 'Detection Frequency', 
         main = paste(tag_id, 'Detections per Day'))
    # Plotting abline for tagging date
    abline(v = min(indv_data$study_date), col = zissou.red) # since tagging date
    # represented as first detection, use min study date
    dev.off
  }
  
  for (i in 1:length(bottomfish_tag_ids)){
    indv_data = vue_data[vue_data$tag_id == tag_id, ]
    
    ## Determining if a detection represents a tag at a new station
    date_of_station_change = c()
    for (a in 2:length(indv_data$station)){
      if (indv_data$station[a] != indv_data$station[a-1]){
        date_of_station_change = c(date_of_station_change, indv_data$study_date[a])
      }
    }
    pdf(paste(tag_id, 'Detections Per Day Histogram.pdf'))
    par(mfrow = c(1,1), oma=c(0,0,2,0))
    hist(floor(indv_data$study_date), 
         breaks = max(vue_data[,1])+1, 
         xlim = c(0, max(vue_data[,1])), xlab = 'Study Date',
         ylab = 'Daily Detections', 
         main = paste(tag_id, 'Detections per Day')) 
    ## plotting abline for tagging date
    abline(v = min(indv_data$study_date), col = zissou.red) # since tagging date
    # represented as first detection, use min study date
    ## Plotting abline for dates when transmitter is detected at a new station
    abline(v = date_of_station_change, col = zissou.blue)
    dev.off()
  }  
}


### distance traveled
distance_traveled = function(vue_data, bottomfish_tag_ids = FALSE){
  ## Function to measure the total distance in km between subsequent receivers
  ## that a fish was detected
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  individual_distance = matrix(0, length(bottomfish_tag_ids), 1)
  for (i in 1:length(bottomfish_tag_ids)){
    subset = vue_data[vue_data$tag_id == bottomfish_tag_ids[i],]
    for (a in 2:length(subset$lon)){
      individual_distance[i] = individual_distance[i] + lldist(point1 = c(subset$lon[a-1], subset$lat[a-1]), point2 = c(subset$lon[a], subset$lat[a]))
    }
  }
  total_distance_tracked = sum(individual_distance)
  print (total_distance_tracked)
  #pdf('distance_tracked_boxplot.pdf')
  #frame()
  #boxplot(individual_distance, main = 'Distance Tracked (km)', ylab = 'km' )
  #dev.off()
  #individual_distance_5num = fivenum(individual_distance[1:4])
  #print (individual_distance_5num)
  distance_traveled_out = as.data.frame(cbind(bottomfish_tag_ids,individual_distance))
  colnames(distance_traveled_out) = c('tag_id', 'distance_tracked')
  return(distance_traveled_out)
}

distance_per_day = function(vue_data, tagging_data, bottomfish_tag_ids =FALSE){
  ## Function to measure the total distance / day at liberty in km between 
  ## subsequent receivers that a fish was detected at
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  distance_per_day = distance_traveled(vue_data, bottomfish_tag_ids)[,2] / 
    transmission_stats(vue_data, tagging_data, bottomfish_tag_ids)$time_at_liberty
  print(fivenum(distance_per_day))
  boxplot(distance_per_day, main = 'Distance Standardized by time at liberty', ylab = 'km/day')
  distance_per_day_out = as.data.frame(cbind(bottomfish_tag_ids, distance_per_day))
  colnames(distance_per_day_out) = c('tag_id', 'distance/day')
  return(distance_per_day_out)
}

number_of_movements = function(vue_data, bottomfish_tag_ids = FALSE){
  ## Function to return the number of movements between stations a transmitter was
  ## detected making
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  movements = matrix(-1, length(bottomfish_tag_ids), 2) ## This begins at -1 because tagging location shouldn't count as a movement between stations
  movements[ ,1] = bottomfish_tag_ids
  for (i in 1:length(bottomfish_tag_ids)){
    indv_data =  vue_data[vue_data$tag_id == bottomfish_tag_ids[i],]
    for (a in 2:length(indv_data$station)){
      if (isTRUE(indv_data$station[a] != indv_data$station[a-1])){
        movements[i,2] = movements[i,2]+1
      }
    }
  }
  movements = as.data.frame(movements)
  colnames(movements) = c('tag_id', 'detected_movements')
  return (movements)
}

number_of_movements_per_day = function(vue_data, tagging_data, bottomfish_tag_ids =FALSE){
  ## Function to return the number of movements between stations a transmitter was
  ## detected making standardized by dividing by transmitter's days at liberty
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  move_per_day = number_of_movements(vue_data, bottomfish_tag_ids)$detected_movements / transmission_stats(vue_data, tagging_data, bottomfish_tag_ids)$time_at_liberty
  move_per_day_out = as.data.frame(cbind(bottomfish_tag_ids, move_per_day))
  colnames(move_per_day_out)[1] = 'tag_id'
  colnames(move_per_day_out)[2] = 'detected movements per day at liberty'
  return(move_per_day_out)
}

brfa_crossings = function(vue_data, tagging_data, bottomfish_tag_ids =FALSE){
  ## Function to return the number of movements between stations a transmitter was
  ## detected making across BRFA boundaries
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  movement_matrix = as.data.frame(matrix(0, length(bottomfish_tag_ids), 3))
  movement_matrix[ ,1] = bottomfish_tag_ids
  for (i in 1:length(bottomfish_tag_ids)){
    indv_data = vue_data[vue_data$tag_id == bottomfish_tag_ids[i], ]
    if (isTRUE(length(indv_data) > 1)){
      for (a in 2:length(indv_data$station)){
        if (in_brfa(indv_data$lon[a], indv_data$lat[a]) != in_brfa(indv_data$lon[a-1], indv_data$lat[a-1])){
          movement_matrix[i,2] = movement_matrix[i,2]+1
        }
      }
    }
    movement_matrix[i,3] = movement_matrix[i,2] / 
      transmission_stats(vue_data, tagging_data, bottomfish_tag_ids[i])$time_at_liberty
  }
  colnames(movement_matrix)[1] = 'tag_id'
  colnames(movement_matrix)[2] = 'BRFA_crossings'
  colnames(movement_matrix)[3] = 'BRFA_crossings_per_day'
  return (movement_matrix)
}

in_brfa_e = function(lon, lat){
  if((lon >= -157.68333333 && lon <= -157.53333333) && 
       (lat >= 21.28333333 && lat <= 21.4166666)){
    return(in_bfra = TRUE)
  }else{ 
    return(FALSE)
  }
}

in_brfa_f = function(lon, lat){
  if((lon >= -157.5666667 && lon <= -157.3666667) && 
       (lat >= 20.9666667 && lat <= 21.333333333)){
    return(in_bfra = TRUE)
  }else{ 
    return(FALSE)
  } 
}

in_brfa = function(lon, lat){
  if(isTRUE(in_brfa_e(lon, lat)) || isTRUE(in_brfa_f(lon,lat))){
    return (TRUE)
  }else{
    return(FALSE)
  }
}

brfa_crossings_per_day = function(vue_data, tagging_data, bottomfish_tag_ids = FALSE){
  ## Function to return the number of movements between stations a transmitter was
  ## detected making across BRFA boundaries standardized by dividing by 
  ## transmitter's days at liberty
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  cross_per_day = brfa_crossings(vue_data, bottomfish_tag_ids)[,2] /
    transmission_stats(vue_data, tagging_data, bottomfish_tag_ids)$time_at_liberty
  cross_per_day_out = as.data.frame(cbind(bottomfish_tag_ids, cross_per_day))
  colnames(cross_per_day_out)[1] = 'tag_id'
  colnames(cross_per_day_out)[2] = 'BRFA_crossings_per_day'
  return(cross_per_day_out)
}

dates_of_location_switching = function(vue_data, bottomfish_tag_ids =FALSE){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  date_changes_all = c()
  for (i in 1:length(bottomfish_tag_ids)){
    indv_data = vue_data[vue_data$tag_id == bottomfish_tag_ids[i],]
    for (a in 2:length(indv_data$station)){
      if (isTRUE(indv_data$station[a] != indv_data$station[a-1])){
        date_changes_all = rbind(date_changes_all, indv_data[a-1, ], indv_data[a, ])
      }
    }
  }
  date_changes_all = as.data.frame(date_changes_all)
  #colnames(date_changes_all)[1] = 'tag_id'
  #colnames(date_changes_all)[2] = 'receiver'
  #colnames(date_changes_all)[3] = 'location'
  #colnames(date_changes_all)[4] = 'date_of_detection'
  #colnames(date_changes_all)[5] = 'arrive/depart' # 0 = depature at a receiver, 1 = arrival
  return (date_changes_all)
}

arrive_date = function(vue_data, bottomfish_tag_ids =FALSE){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  arrive_dates = c()
  for (i in 1:length(bottomfish_tag_ids)){
    indv_data = vue_data[vue_data$tag_id == bottomfish_tag_ids[i], ]
    for (a in 2:length(indv_data$station)){
      if(isTRUE(indv_data$station[a] != indv_data$station[a-1])){
        arrive_dates = rbind(c(bottomfish_tag_ids[i], indv_data$station[a], indv_data$datetime[a]))}
    }
    arrive_dates = as.data.frame(arrive_dates)
    colnames(arrive_dates) = c('tag_id', 'station_arrived', 'arrival_study_date')
  }
  return(arrive_dates)
}

depart_date = function(vue_data, bottomfish_tag_ids =FALSE){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  depart_dates = c()
  for (i in 1:length(bottomfish_tag_ids)){
    indv_data = vue_data[vue_data$tag_id == bottomfish_tag_ids[i], ]
    for (a in 2:length(indv_data$station)){
      if(isTRUE(indv_data$station[a] != indv_data$station[a-1])){
        depart_dates = rbind(c(bottomfish_tag_ids[i], indv_data$station[a-1], indv_data$study_date[a-1]))}
    }
    depart_dates = as.data.frame(depart_dates)
    colnames(depart_dates) = c('tag_id', 'station_departed', 'departure_study_date')
  }
  return(depart_dates)
}

all_stations_detected = function(vue_data, bottomfish_tag_ids = FALSE){
  ## Function to return a list of each subsequent station a transmitter was detected
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])
    to_return = FALSE
  }else{
    to_return = TRUE}
  for (i in 1:length(bottomfish_tag_ids)){
    indv_data = vue_data[vue_data$tag_id == bottomfish_tag_ids[i], ]
    all_stations = c(paste('Detection History for Tag no. ', bottomfish_tag_ids[i], ': '), indv_data$station[1])
    if (length(indv_data)>1){
      for (a in 2:length(indv_data$station)){
        if (isTRUE(indv_data$station[a] != indv_data$station[a-1])){
          all_stations = as.character(c(all_stations, ',', indv_data$station[a]))
        }
      }
    }
    if(to_return == FALSE){
      print(c(all_stations))
    }
  }
  if(to_return == TRUE){
    return(cat(all_stations))
  }
}

unique_stations_detected = function(vue_data, bottomfish_tag_ids = FALSE){
  ## Function to return a list of each unique station a transmitter was detected
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  for (i in 1:length(bottomfish_tag_ids)){
    indv_data = vue_data[vue_data$tag_id == bottomfish_tag_ids[i], ]
    unique_station_detections = unique(indv_data$station)
    unique_stations = c(paste('Unique stations visited for tag no.', bottomfish_tag_ids[i], ': '), unique_station_detections)
    print(unique_stations)
  }
}

### Plotting_movements
plotting_movements = function(vue_data, receiver_data, bottomfish_tag_ids =FALSE){
  ## Function to return the number of movements between stations a transmitter was
  ## detected making across BRFA boundaries
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  
  ## Defining map shapefile
  pngMAP_df<- get_map(location = c(lon = -157.75, lat = 21.251), 
                      source = "google", zoom = 9,color='bw')
  
  ## Cycling through tags and plotting
  for (tag_id in bottomfish_tag_ids){
    lonlats = as.matrix(cbind(vue_data$lon[vue_data$tag_id == tag_id], vue_data$lat[vue_data$tag_id == tag_id]))
    ll_plot_index = 1
    if (length(lonlats[ ,1]) > 1){
      for (i in 2:length(lonlats[ ,1])){
        if (lonlats[i,1] != lonlats[i-1,1] || lonlats[i,2] != lonlats[i-1,2]){
          ll_plot_index = c(ll_plot_index, i)
        }
      }
      plot_lonlats = as.data.frame(lonlats[ll_plot_index, ])
      colnames(plot_lonlats) = c('lon', 'lat')
    }
    
    plot.new()
    
    png(sprintf('%s Movement Map.png', tag_id))
    print(ggmap(pngMAP_df) + 
            geom_point(color = zissou.gold, size = 2, data = receiver_data,
                       aes(x = lon, y = lat)) +
            geom_path(color = zissou.blue, mapping = aes(x = c(-157.566, -157.566), y = c(21.0333, 20.9166))) +
            geom_path(color = zissou.blue, mapping = aes(x = c(-157.566, -157.366), y = c(20.9166, 20.9166))) +
            geom_path(color = zissou.blue, mapping = aes(x = c(-157.366, -157.366), y = c(20.9166, 21.0333))) +
            geom_path(color = zissou.blue, mapping = aes(x = c(-157.366, -157.566), y = c(21.0333, 21.0333))) +
            geom_path(color = zissou.blue, mapping = aes(x = c(-157.683, -157.533), y = c(21.4166, 21.4166))) +
            geom_path(color = zissou.blue, mapping = aes(x = c(-157.533, -157.533), y = c(21.4166, 21.2833))) +
            geom_path(color = zissou.blue, mapping = aes(x = c(-157.683, -157.533), y = c(21.2833, 21.2833))) +
            geom_point(color = zissou.red, size = 2, data = plot_lonlats, aes(x = lon, y = lat)) +
            geom_path(color = zissou.red, size = 2, data = plot_lonlats, aes(x = lon, y = lat)) + 
            labs(title = sprintf('Detections of Tag %s', tag_id)))
    dev.off()
  }
}

get_tagging_metadata = function(tagging_data, bottomfish_tag_ids = FALSE){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  tagging_meta = as.data.frame(matrix(0, length(bottomfish_tag_ids), 4))
  for (i in 1:length(bottomfish_tag_ids)){
    tagging_meta[i,1] = as.numeric(as.character(tagging_data$vem_tag_id[which(tagging_data$vem_tag_id == bottomfish_tag_ids[i])]))
    tagging_meta[i,2] = as.character(tagging_data$"species"[which(tagging_data$vem_tag_id == bottomfish_tag_ids[i])])
    tagging_meta[i,3] = as.numeric(as.character(tagging_data$"fork_length(cm)"[which(tagging_data$vem_tag_id == bottomfish_tag_ids[i])]))
    tagging_meta[i,4] = as.character(tagging_data$"area_of_capture"[which(tagging_data$vem_tag_id == bottomfish_tag_ids[i])])
  }
  colnames(tagging_meta) = c('tag_id', 'species', 'fork_length', 'tagging_location')
  return(tagging_meta)
}

unique_receivers = function(vue_data, tagging_data, bottomfish_tag_ids = FALSE){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])}
  receivers_visited = as.data.frame(matrix(0, length(bottomfish_tag_ids), 3))
  for (i in 1:length(bottomfish_tag_ids)){
    indv_data = vue_data[which(vue_data$tag_id == bottomfish_tag_ids[i]), ]
    receivers_visited[i,1] = as.numeric(as.character(tagging_data$vem_tag_id[which(tagging_data$vem_tag_id == bottomfish_tag_ids[i])]))
    receivers_visited[i,2] = length(unique(indv_data$station)) - 1 # Subtract to account for tagging location
    receivers_visited[i,3] = receivers_visited[i,2] / transmission_stats(vue_data, tagging_data, bottomfish_tag_ids = bottomfish_tag_ids[i])$time_at_liberty
  }
  colnames(receivers_visited) = c('tag_id', 'receivers_detected', 'receivers_detected_per_day')
  return(receivers_visited)
}

## Making strip chart of detections
vemco_stripchart = function(vue_data, bottomfish_tag_ids = FALSE){
  plot_data = clean_vue(vue_data, bottomfish_tag_ids)
  pdf('Stripchart of Detections.pdf')
  par(las = 1)
  stripchart(study_date~tag_id, data = plot_data,
             #ylab = 'Tag ID',
             xlab = 'Study Date',
             main = 'Detections During Study Period')
  dev.off()
}

spatial_evenness = function(vue_data, receiver_data){
  ### function to calculate spacitail eveness based on Pielou 1966 from TinHan 2014
  ## outputs a dataframe first column is tag id, second column is spacial evenness index
  spacial_evenness_df = as.data.frame(matrix(data = 0, nrow = length(unique(vue_data$tag_id)), ncol = 2))
  colnames(spacial_evenness_df) = c('tag_id', 'spacial_evenness_metric')
  spacial_evenness_df$tag_id = unique(vue_data$tag_id)
  R = length(unique(receiver_data$station_name)) #could replaces "station_name" with "zone"
  for(i in 1:length(unique(vue_data$tag_id))){
    indv_data = vue_data[which(vue_data$tag_id == unique(vue_data$tag_id)[i]), ]
    spatial_sum = c()
    for(a in 1:length(unique(receiver_data$station_name))){
      rho_i = length(indv_data$datetime[which(as.character(indv_data$station) == as.character(receiver_data$station_name[a]))]) / length(indv_data$station)
      spatial_sum = c(spatial_sum, (rho_i * log(rho_i)))
    }
    print(spatial_sum)
    spacial_evenness_df[i, 2] = (-1 * sum((spatial_sum[spatial_sum != 'NaN']))) / log(R)
  }
  return(spacial_evenness_df)
}

calculate_time_between_detections = function(vue_data, bottomfish_tag_ids == FALSE){
  if (bottomfish_tag_ids[1] == FALSE){
    bottomfish_tag_ids = unique(as.numeric(levels(vue_data$tag_id))[vue_data$tag_id])
  }
  time_between_df = as.data.frame(matrix(0, length(bottomfish_tag_ids), 7))
  colnames(time_between_df) = c('tag_id', 'min', '1st_quartile', 'median', 'third_quartile', 'max', 'mean')
  for (i in 1:length(bottomfish_tag_ids)){
    indv_data = vue_data[vue_data$tag_id == bottomfish_tag_ids[i], ]
    diff_time = c()
    for (a in 3:length(indv_data$datetime)){ # Begin indexing at 3 because first detection is generated tagging date
      diff_time = c(diff_time, difftime(time1 = indv_data$datetime[a], time2 = indv_data$datetime[a-1], units = 'days'))
    }
    time_between_df[i, 2:6] = fivenum(diff_time)
    time_between_df[i, 7]   = mean(diff_time)
  }
  return(time_between_df)
}

#### USAGE -------------------------------------


### Importing Data -----------------------------

## Importing receiver data

receiver_data = load_receiver(filename = 'DEPLOYMENT_RECOVERY_LOG.csv', filepath = '/Users/stephenscherrer/Dropbox/Lab Folder/Oahu Receiver Data Files/')

## Import tagging log
tagging_data = load_tagging_data('/Users/stephenscherrer/Dropbox/Lab Folder/Oahu Receiver Data Files/Bottomfish_Tag_Master.csv')

## Importing VUE Data
vue_data = load_vemco('/Users/stephenscherrer/Dropbox/Lab Folder/Oahu Receiver Data Files/VUE_Export_2015-Jan-20.csv')

### Cleaning data

## Fixing missing Lat and Lons
vue_data = clean_vue_lat_lon(vue_data, receiver_data)

## Adding tagging location for all tags present
vue_data = generate_tagging_detection(tagging_data, vue_data)

## Removing detections from botcam crew
vue_data = remove_location(vue_data, 'With BotCam Crew')

## Pulling out bottomfish tag ids
ws_tag_ids = c(633, 68295, 8893, 8832, 1425, 29700, 2139, 8892, 9095, 
            3430, 3580, 3407, 18616, 3265, 52933, 68291, 630, 62142, 
            88758, 83056, 83054, 88767, 88765, 46863, 46862, 46861, 
            58432, 100521, 29711, 3339, 3281, 3402, 18596, 29710, 
            41761, 41753, 41760, 41739, 41748, 41773, 52928, 49019, 
            52930, 52934, 62015, 61910, 62016, 61907, 61902, 62012, 
            62019, 62011, 62010, 61899, 61904, 62020, 61898, 62014, 
            62013, 52873, 62021, 60982, 60983, 62017, 68281, 68308, 
            642, 68321, 68286, 68320, 68319, 68290, 624, 637, 622, 
            640, 635, 623, 62121, 46109, 618, 52446, 52456, 52458, 
            111194, 87568, 68317, 10565, 61911, 61914, 62018, 61909, 
            6841, 627, 6847, 6840, 645, 78173, 655, 654, 6843, 52426, 
            52422, 52927, 61900, 111190, 58446, 617, 52428, 46114, 
            52455, 659, 46111, 46118, 58460, 58453, 58445, 78159, 646, 
            62145, 52436, 87568, 653, 616, 52452, 46116, 52454, 52425, 
            68272, 12827, 46108, 52435, 46105, 52434, 52433, 52442, 
            87597, 52419, 62127, 62146, 46110, 62140, 46115, 46106, 
            46097, 59108, 58438, 58431, 58448, 58456, 58458, 6858, 
            62122, 62125, 62126, 62131, 12828, 68251, 68315, 46107, 
            52449, 58447, 52039, 61906, 59112, 52929, 46117, 46094, 
            52440, 631, 78172, 62137, 68253, 629, 68282, 68314, 68305, 
            61908, 68252, 68306, 68309, 641, 648, 651, 657, 6844, 62138, 
            6857, 62133, 62124, 46104, 46100, 61905, 52450, 78168, 6845, 
            58455, 52429, 78145, 6851, 52439, 46103, 52432, 46099, 58452, 
            58434, 46102, 2994, 111182, 58437, 59114, 108314, 59106, 33560, 
            33543, 102061, 656, 625, 68307, 59105, 62147, 660, 52041, 87562, 
            58435, 6852, 58459, 111193, 58457, 58442, 33549, 59110, 33539, 
            33535, 33533, 33531, 68310, 103300, 33530, 33540, 33541, 33538, 
            59109, 33532, 33537, 33542, 33548, 68307, 59113, 87602, 658, 52457, 
            46095, 33546, 58433, 6859, 78153, 68249, 62149, 58441, 78176, 52421, 
            46098, 638, 68316, 6848, 59107, 6850, 58443, 78177, 111195, 
            628, 78175, 52420, 639, 62148, 62022, 46101, 52438, 52423, 
            58444, 62134, 78171, 68312, 6855, 58440, 62150, 52431, 636, 
            68287, 78152, 643, 62130, 87569, 6853, 78141, 62132, 68311, 
            62144, 52424, 62139, 6846, 78169, 78136, 103307, 52437, 46112, 
            58449, 49018, 6856, 78167, 650, 68250, 52427, 58439, 632, 
            58451, 52451, 49075, 87563, 85532, 59111, 33545, 33544, 33547, 
            62123, 58450, 87570, 6854, 46096, 52430, 52447, 62128, 62136, 
            62141, 46113, 62143, 62129, 52444, 52441, 6849, 58454, 58436, 
            6842 )
##NOTE: Removed tag 57451 because it was detected at cross and moved a lot and fucked up analysis

all_bottomfish_tags = ws_tag_ids
# Removing data from unwanted tags
vue_data = clean_vue(vue_data, ws_tag_ids)
dead_fish = clean_vue(vue_data, c(37969, 57459))
vue_data = clean_vue(vue_data, c(57451, 37954, 37969, 57459), exclude = TRUE) # first 2 tags are at cross seamount. cannot find tags in hard copy logs. probably monchong tagged with BF tags. Third and fourth tags are very dead fish right next to a recevier
opakapaka_tag_ids = ws_tag_ids[!(ws_tag_ids %in% c(57451, 37954, 37969, 57459))]

# Removing botcam detections for which I don't have time to look up lon/lats
vue_data = remove_location(vue_data, 'With BotCam Crew')

# Add column to vue_data data frame for what day of study a detection occurred
vue_data = adjust_vue_study_dates(vue_data, tagging_data, opakapaka_tag_ids)
dead_fish = adjust_vue_study_dates(dead_fish, tagging_data, opakapaka_tag_ids)

# Add column to tagging_data frame for what day of study tagging occurred
tagging_data = adjust_tagging_study_dates(vue_data, tagging_data, opakapaka_tag_ids)

# add columns to receiver_data frame for what day of study deployments and recoveries occured
receiver_data = adjust_receiver_study_dates(receiver_data, tagging_data)
write.csv(receiver_data, file = "receiver_data_with_study_dates.csv")
# Determining first, last and elapsed number of days in study
exp_dates = experiment_dates(vue_data)
print(exp_dates)

# pulling first detection, last detection, and time at liberty for each tag
all_transmission_stats = transmission_stats(vue_data, tagging_data)
print(all_transmission_stats)

# pulling out tags with > 7 days detection history
tags_of_interest = all_transmission_stats$tag_id[all_transmission_stats$time_at_liberty>7]
toi_transmission_stats = transmission_stats(vue_data, tagging_data, tags_of_interest)
print(toi_transmission_stats)

# getting meta data associated with transmitters
all_meta_data = get_tagging_metadata(tagging_data)
toi_meta_data = get_tagging_metadata(tagging_data, tags_of_interest)

# getting tagging dates
all_tagging_date = tagging_date(tagging_data, unique(as.numeric(vue_data$tag_id)))
toi_tagging_date = tagging_date(tagging_data, tags_of_interest)

# number of detections and percentage of total detections for each individual
all_number_of_detections = number_of_detections(vue_data)
print(all_number_of_detections)

toi_number_of_detections = number_of_detections(vue_data, tags_of_interest)
print(toi_number_of_detections)

# number of detections/day for each individual
#### IS THIS ALL DAYS OF STUDY OR ALL DAYS AT LIBERTY?
all_detections_per_day = detections_per_day(vue_data, tagging_data)
print(all_detections_per_day)

toi_detections_per_day = detections_per_day(vue_data, tagging_data, tags_of_interest)
print(toi_detections_per_day)

# unique days detected for each individual
all_unique_days = unique_days_detected(vue_data)
print(all_unique_days)

toi_unique_days = unique_days_detected(vue_data, tags_of_interest)
print(toi_unique_days)

# distance traveled by individual
all_distance_traveled = distance_traveled(vue_data)
print(all_distance_traveled)

toi_distance_traveled = distance_traveled(vue_data, tags_of_interest)
print(toi_distance_traveled)

# distance/day at liberty
all_distance_per_day = distance_per_day(vue_data, tagging_data)
print(all_distance_per_day)

toi_distance_per_day = distance_per_day(vue_data, tagging_data, tags_of_interest)
print(toi_distance_per_day)

# All/Unique stations detected per individual
all_stations_detected(vue_data)
unique_stations_detected(vue_data, opakapaka_tag_ids)

all_stations_detected(vue_data, tags_of_interest)
unique_stations_detected(vue_data, tags_of_interest)

# number of movements
all_movements = number_of_movements(vue_data)
print(all_movements)

toi_movements = number_of_movements(vue_data, tags_of_interest)
print(toi_movements)

#number of movements / day at liberty
all_movements_per_day = number_of_movements_per_day(vue_data, tagging_data)
print(all_movements_per_day)

toi_movements_per_day = number_of_movements_per_day(vue_data, tagging_data, tags_of_interest)
print(toi_movements_per_day)

# number of BRFA crossings and crossings/day
all_brfa_crossings = brfa_crossings(vue_data, tagging_data)
print (all_brfa_crossings)

toi_brfa_crossings = brfa_crossings(vue_data, tagging_data, tags_of_interest)
print(toi_brfa_crossings)


# number of receivers visited
all_receivers_detected = unique_receivers(vue_data, tagging_data)
print(all_receivers_detected)

toi_receivers_detected = unique_receivers(vue_data, tagging_data, tags_of_interest)
print(toi_receivers_detected)


### Combining all data:
all_data_out = cbind(all_meta_data,
                     all_transmission_stats, 
                     all_tagging_date[ ,-1],
                     all_number_of_detections[ ,-1],
                     all_detections_per_day[ ,-1],
                     all_unique_days[ ,-1],
                     all_movements[ ,-1],
                     all_movements_per_day[ ,-1],
                     all_brfa_crossings[ ,-1],
                     all_distance_traveled[ ,-1],
                     all_distance_per_day[ ,-1],
                     all_receivers_detected[ ,-1]
)

toi_data_output = cbind(toi_meta_data,
                        toi_transmission_stats, 
                        toi_tagging_date[ ,-1],
                        toi_number_of_detections[ ,-1],
                        toi_detections_per_day[ ,-1],
                        toi_unique_days[ ,-1],
                        toi_movements[ ,-1],
                        toi_movements_per_day[ ,-1],
                        toi_brfa_crossings[ ,-1],
                        toi_distance_traveled[ ,-1],
                        toi_distance_per_day[ ,-1],
                        toi_receivers_detected[ ,-1]
)


sink("bottomfish_analysis_output.txt")
for(i in 1:length(all_data_out$tag_id)){
  cat(paste("Transmitter ID: ", all_data_out$tag_id[i], "\n"))
  cat(paste("Species: ", all_data_out$species[i], "\n"))
  cat(paste("Fork Length: ", all_data_out$fork_length[i], " cm", "\n"))
  cat(paste("Tagging Date: ", all_data_out$tagging_date[i], "\n"))
  cat(paste("Tagging Location: ", all_data_out$tagging_location[i], "\n"))  
  cat(paste("Day of First Detected: ", all_data_out$first_transmission[i], "\n"))
  cat(paste("Day of Last Detection: ", all_data_out$last_transmission[i], "\n"))
  cat(paste("Days at Liberty: ", all_data_out$time_at_liberty[i], "\n"))
  cat(paste("Unique Days Detected: ", all_data_out$"all_unique_days[, -1]"[i], "\n"))
  cat(paste("Number of Detections (Total): ", all_data_out$"#_of_detections"[i], "\n"))
  cat(paste("Number of Detections / Day: ", all_data_out$"all_detections_per_day[, -1]"[i], "\n"))
  cat(paste("Number of Receivers Detected (Total): ", all_data_out$"receivers_detected"[i], "\n"))
  cat(paste("Number of receivers Detected / Day: ", all_data_out$receivers_detected_per_day[i], "\n"))
  cat(paste("Movements Into/Out of BRFA (Total): ", all_data_out$BRFA_crossings[i], "\n"))
  cat(paste("Movements Into/Out of BRFA / Day: ", all_data_out$BRFA_crossings_per_day[i], "\n"))
  cat(paste("Approximate Distance Tracked (Total): ", all_data_out$all_distance_traveled[i], "km","\n"))
  cat(paste("Approximate Distance Tracked / Day: ", all_data_out$all_distance_per_day[i],"km", "\n")) 
  cat(paste(all_stations_detected(vue_data, all_data_out$tag_id[i]), "\n"))
  cat("\n")
}
sink()

write.csv(all_data_out, file = 'bottomfish_summary.csv')

plot.new()
pdf('Number of Detections by Tag Boxplot.pdf')
par(mfrow = c(1,1), oma=c(0,0,2,0))  
median = signif(fivenum(toi_data_output$"#_of_detections")[3], 3)
first_q = signif(fivenum(toi_data_output$"#_of_detections")[2], 3)
third_q = signif(fivenum(toi_data_output$"#_of_detections")[4], 3)
boxplot(toi_data_output$"#_of_detections", main = 'Number of Detections by Tag', ylab = 'Detections')
legend('topright', cex = .75, legend = c(sprintf('1st Quartile: %s', first_q), sprintf('Median: %s', median), sprintf('3rd Quartile: %s', third_q)),
       col = 'black')
dev.off()

plot.new()
pdf('Percentage of All Detections By Tag Boxplot.pdf')
median = signif(fivenum(toi_data_output$"%_of_all_detections" )[3], 3)
first_q = signif(fivenum(toi_data_output$"%_of_all_detections" )[2], 3)
third_q = signif(fivenum(toi_data_output$"%_of_all_detections" )[4], 3)
boxplot(toi_data_output$"%_of_all_detections" , main = 'Percentage of All Detections by Tag', ylab = '%')
legend('topright', cex = .75, legend = c(sprintf('1st Quartile: %s', first_q), sprintf('Median: %s', median), sprintf('3rd Quartile: %s', third_q)),
       col = 'black')
dev.off()

plot.new()
pdf('Mean Detections Per Day Boxplot.pdf')
median = signif(fivenum(toi_data_output$"toi_detections_per_day[, -1]")[3], 3)
first_q = signif(fivenum(toi_data_output$"toi_detections_per_day[, -1]")[2], 3)
third_q = signif(fivenum(toi_data_output$"toi_detections_per_day[, -1]")[4], 3)
boxplot(toi_data_output$"toi_detections_per_day[, -1]", main = 'Mean Detections By Day', ylab = '# Detections/Day')
legend('topright', cex = .75, legend = c(sprintf('1st Quartile: %s', first_q), sprintf('Median: %s', median), sprintf('3rd Quartile: %s', third_q)),
       col = 'black')
dev.off()

plot.new()
pdf('Unique Days Detected Boxplot.pdf')
median = signif(fivenum(toi_data_output$"toi_unique_days[, -1]")[3], 3)
first_q = signif(fivenum(toi_data_output$"toi_unique_days[, -1]")[2], 3)
third_q = signif(fivenum(toi_data_output$"toi_unique_days[, -1]")[4], 3)
boxplot(toi_data_output$"toi_unique_days[, -1]" , main = 'Unique Days Detected', ylab = 'Days')
legend('topright', cex = .75, legend = c(sprintf('1st Quartile: %s', first_q), sprintf('Median: %s', median), sprintf('3rd Quartile: %s', third_q)),
       col = 'black')
dev.off()

plot.new()
pdf('Total Distance Tracked Boxplot.pdf')
median = signif(fivenum(toi_data_output$"toi_distance_traveled[, -1]")[3], 3)
first_q = signif(fivenum(toi_data_output$"toi_distance_traveled[, -1]")[2], 3)
third_q = signif(fivenum(toi_data_output$"toi_distance_traveled[, -1]")[4], 3)
boxplot(toi_data_output$"toi_distance_traveled[, -1]", main = 'Total Distance Tracked', ylab = 'km')
legend('topright', cex = .75, legend = c(sprintf('1st Quartile: %s', first_q), sprintf('Median: %s', median), sprintf('3rd Quartile: %s', third_q)),
       col = 'black')
dev.off()

plot.new()
pdf('Distance Tracked Per Day Boxplot.pdf')
median = signif(fivenum(toi_data_output$"toi_distance_per_day[, -1]")[3], 3)
first_q = signif(fivenum(toi_data_output$"toi_distance_per_day[, -1]")[2], 3)
third_q = signif(fivenum(toi_data_output$"toi_distance_per_day[, -1]")[4], 3)
boxplot(toi_data_output$"toi_distance_per_day[, -1]", main = 'Distance Tracked/Day', ylab = 'km/Day')
legend('topright', cex = .75, legend = c(sprintf('1st Quartile: %s', first_q), sprintf('Median: %s', median), sprintf('3rd Quartile: %s', third_q)),
       col = 'black')
dev.off()

plot.new()
pdf('Number of Unique Receivers Visited Boxplot.pdf')
median = signif(fivenum(toi_data_output$receivers_detected)[3], 3)
first_q = signif(fivenum(toi_data_output$receivers_detected)[2], 3)
third_q = signif(fivenum(toi_data_output$receivers_detected)[4], 3)
boxplot(toi_data_output$receivers_detected, main = 'Receivers Visited', ylab = '# Receivers')
legend('topright', cex = .75, legend = c(sprintf('1st Quartile: %s', first_q), sprintf('Median: %s', median), sprintf('3rd Quartile: %s', third_q)),
       col = 'black')
dev.off()

plot.new()
pdf('Movements Between Receivers By Day Boxplot.pdf')
median = signif(fivenum(toi_data_output$"toi_movements_per_day[, -1]")[3], 3)
first_q = signif(fivenum(toi_data_output$"toi_movements_per_day[, -1]")[2], 3)
third_q = signif(fivenum(toi_data_output$"toi_movements_per_day[, -1]")[4], 3)
boxplot(toi_data_output$"toi_movements_per_day[, -1]", main = 'Movements Between Receivers/Day', ylab = 'Movements/Day')
legend('topright', cex = .75, legend = c(sprintf('1st Quartile: %s', first_q), sprintf('Median: %s', median), sprintf('3rd Quartile: %s', third_q)),
       col = 'black')
dev.off()

plot.new()
pdf('Total BRFA Crossings Boxplot.pdf')
median = signif(fivenum(toi_data_output$BRFA_crossings)[3], 3)
first_q = signif(fivenum(toi_data_output$BRFA_crossings)[2], 3)
third_q = signif(fivenum(toi_data_output$BRFA_crossings)[4], 3)
boxplot(toi_data_output$BRFA_crossings, main = 'Total BRFA Crossings', ylab = 'BRFA Crossings')
legend('topright', cex = .75, legend = c(sprintf('1st Quartile: %s', first_q), sprintf('Median: %s', median), sprintf('3rd Quartile: %s', third_q)),
       col = 'black')
dev.off()

plot.new()
pdf('BRFA Crossings Per Day Boxplot.pdf')
median = signif(fivenum(toi_data_output$BRFA_crossings_per_day)[3], 3)
first_q = signif(fivenum(toi_data_output$BRFA_crossings_per_day)[2], 3)
third_q = signif(fivenum(toi_data_output$BRFA_crossings_per_day)[4], 3)
boxplot(toi_data_output$BRFA_crossings_per_day, main = 'BRFA Crossings/Day', ylab = 'BRFA Crossings/Day')
legend('topright', cex = .75, legend = c(sprintf('1st Quartile: %s', first_q), sprintf('Median: %s', median), sprintf('3rd Quartile: %s', third_q)),
       col = 'black')
dev.off()


## Plotting out detections by study date
for (i in 1:length(tags_of_interest)){
  tag_id = tags_of_interest[i]
  title = sprintf('%s Study Date of Detections.pdf', tag_id)
  pdf(title)
  par(mfrow = c(1,1))
  plot_title = sprintf('Detections of Tag %s', tag_id)
  indv_data = vue_data[which(vue_data$tag_id == tag_id), ]
  hist(indv_data$study_date,  
       breaks = ceiling(max(vue_data$study_date)), 
       main = plot_title, xlab = 'Days at Liberty (Binned Daily)', 
       ylab = 'Number of Transmissions', ylim = c(0,30), 
       xlim = c(0, max(vue_data$study_date))+5)
  abline(v = tagging_data$study_date[which(tagging_data$vem_tag_id == tag_id)], col = 'red')
  dates = dates_of_location_switching(vue_data, tag_id)
  abline(v = dates[dates[,5] == 1, 4], col = 258)
  dev.off()
}


#### For Dead Tags
for (i in 1:length(unique(dead_fish$tag_id))){
  tag_id = dead_fish$tag_id[i]
  title = sprintf('%s Study Date of Detections.pdf', tag_id)
  pdf(title)
  par(mfrow = c(1,1))
  plot_title = sprintf('Detections of Tag %s', tag_id)
  indv_data = dead_fish[which(dead_fish$tag_id == tag_id), ]
  hist(indv_data$study_date,  
       breaks = ceiling(max(vue_data$study_date)), 
       main = plot_title, xlab = 'Days at Liberty (Binned Daily)', 
       ylab = 'Number of Transmissions', ylim = c(0,30), 
       xlim = c(0, max(vue_data$study_date))+5)
  abline(v = tagging_data$study_date[which(as.character(tagging_data$vem_tag_id) == as.character(tag_id))], col = 'red')
  dates = dates_of_location_switching(dead_fish, tag_id)
  abline(v = dates[dates[,5] == 1, 4], col = 258)
  dev.off()
}


## All Detections
pdf('All Detections During Study.pdf')
par(mfrow = c(1,1))
hist(vue_data$study_date, breaks = ceiling(max(vue_data$study_date)/30), 
     main = 'Total Transmissions Detected During Study', 
     xlab = 'Days (Binned Every 30 Days)', xlim = c(0, max(vue_data$study_date)+5),
     ylab = 'Transmissions Detected', ylim = c(0,7000))
dev.off()

## Detection Stripchart
vemco_stripchart(vue_data)
vemco_stripchart(dead_fish)

## Making Maps
## loading in receivers to plot
plot_receivers = load_receiver("/Users/stephenscherrer/Documents/Work/UH/Projects/dissertation work/Spacial Ecology/Bottomfish Analysis Feb 2015/plotting_files/phase_1_deployment_recovery_log.csv")
plotting_movements(vue_data, plot_receivers)




### Plotting Receiver Positions

## Phase 1 acoustic array
plot_receivers_phase_1 = load_receiver("/Users/stephenscherrer/Documents/Work/UH/Projects/dissertation work/Spacial Ecology/Bottomfish Analysis Feb 2015/plotting_files/phase_1_deployment_recovery_log.csv")
## Defining map shapefile
pngMAP_df<- get_map(location = c(lon = -157.75, lat = 21.251), 
                    source = "google", zoom = 9,color='color')

plot.new()

png("Locations of Phase 1 Receivers.png")
print(ggmap(pngMAP_df) + 
        geom_point(color = zissou.red, size = 2, data = plot_receivers_phase_1,
                   aes(x = lon, y = lat)) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.566, -157.566), y = c(21.0333, 20.9166))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.566, -157.366), y = c(20.9166, 20.9166))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.366, -157.366), y = c(20.9166, 21.0333))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.366, -157.566), y = c(21.0333, 21.0333))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.683, -157.533), y = c(21.4166, 21.4166))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.533, -157.533), y = c(21.4166, 21.2833))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.683, -157.533), y = c(21.2833, 21.2833))) +
        labs(title = "Locations of Phase 1 Receivers"))
dev.off()


## Phase 2 acoustic array
plot_receivers_phase_2 = load_receiver("/Users/stephenscherrer/Documents/Work/UH/Projects/dissertation work/Spacial Ecology/Bottomfish Analysis Feb 2015/plotting_files/phase_2_deployment_recovery_log.csv")
## Defining map shapefile
pngMAP_df<- get_map(location = c(lon = -157.75, lat = 21.251), 
                    source = "google", zoom = 9,color='color')

plot.new()

png("Locations of Phase 2 Receivers.png")
print(ggmap(pngMAP_df) + 
        geom_point(color = zissou.red, size = 2, data = plot_receivers_phase_2,
                   aes(x = lon, y = lat)) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.566, -157.566), y = c(21.0333, 20.9166))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.566, -157.366), y = c(20.9166, 20.9166))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.366, -157.366), y = c(20.9166, 21.0333))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.366, -157.566), y = c(21.0333, 21.0333))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.683, -157.533), y = c(21.4166, 21.4166))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.533, -157.533), y = c(21.4166, 21.2833))) +
        geom_path(color = zissou.gold, mapping = aes(x = c(-157.683, -157.533), y = c(21.2833, 21.2833))) +
        labs(title = "Locations of Phase 2 Receivers"))
dev.off()

paka_tagging_data = tagging_data[as.character(tagging_data$species) == "Opakapaka", ]
fivenum(as.numeric(paka_tagging_data$"fork_length(cm)"))

## For each/all fish:
# How long tag active/time at liberty
# Number of detections
# Recovery rate
# Date of last tag
# Histogram of detections per time bin
# Weekly
# By time bin
# Number of receivers detected at.
# Determine receiver deployment configuration during those periods
# Receiver name for each recorded movement along with dates.
# Plot of receiver movements by date
# Determine distance swam
# Number of BRFA crossings
# Number of channel crossings
## For all fish
# Total study durration.
# histogram of all detections by time
# Remove last 2 weeks of tag records to account for predation
# Determine which tags meet criteria for consideration
# Export receiver deployment configurations for processing through web app
# Import web app recovery rates
# t-test tagged fish rates with Ho: recovery rate = model output
# plot recovery rate vs model projections

distances_between_receivers = matrix(0, length(plot_receivers_phase_1$lat), length(plot_receivers_phase_1$lat))
for(i in 1:length(plot_receivers_phase_1$lat)){
  for (a in 1:length(plot_receivers_phase_1$lat)){
    distances_between_receivers[i, a] = lldist(point1 = c(plot_receivers_phase_1$lon[i], plot_receivers_phase_1$lat[i]), point2 = c(plot_receivers_phase_1$lon[a], plot_receivers_phase_1$lat[a]))
  }
}



minimum_distance_between_receivers = min(distances_between_receivers[distances_between_receivers != 0])
maximum_distance_between_receivers = max(distances_between_receivers)
mean_distance_between_receivers = mean(distances_between_receivers)


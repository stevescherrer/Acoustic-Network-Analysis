#### Utility Functions
#### A Collection of useful functions

### installing principle dependencies ----------------
#install.packages('wesanderson')
library('wesanderson')



#### Functions ---------------------------------------

### load_vemco -> Importing VUE Data files -------------------------
get_vemco = function(filename, filepath = FALSE, format = '%Y-%m-%d %H:%M:%S'){
  proj_dir = getwd()
  if (isTRUE(filepath != FALSE)) {setwd(filepath)}
  vue_data_raw = read.csv(filename)
  vue_data_cleaned = vue_col_names(vue_data_raw)
  vue_data_cleaned$datetime = strptime(vue_data_cleaned$datetime, 
                                             format = '%Y-%m-%d %H:%M:%S',
                                             tz = "GMT")
  vue_data_cleaned$datetime = convert_tz(vue_data_cleaned$datetime, new.tz = 'HST')
  vue_data_cleaned$tag_id = clean_tag_id(vue_data_cleaned$tag_id)
  vue_data_cleaned$receiver = clean_receiver(vue_data_cleaned$receiver)
  setwd(proj_dir)
  return (vue_data_cleaned)
}

convert_tz = function(datetime, new.tz = 'HST'){
  #Function to convert GMT/UTC times to HST time
  datetime.new.tz = strptime(datetime, format = '%Y-%m-%d %H:%M:%S', tz = new.tz)
  dateoffset = datetime-datetime.new.tz
  datetime.new.tz = datetime.new.tz + dateoffset
  return(datetime.new.tz)
}

vue_col_names = function(vue_data_raw){
  colnames(vue_data_raw)[1]  <- 'datetime'
  colnames(vue_data_raw)[2]  <- 'receiver'
  colnames(vue_data_raw)[3]  <- 'tag_id'
  colnames(vue_data_raw)[4]  <- 'name'
  colnames(vue_data_raw)[5]  <- 'tag_serial'
  colnames(vue_data_raw)[6]  <- 'sensor_value'
  colnames(vue_data_raw)[7]  <- 'sensor.unit'
  colnames(vue_data_raw)[8]  <- 'station'
  colnames(vue_data_raw)[9]  <- 'lat'
  colnames(vue_data_raw)[10] <- 'lon'
  return (vue_data_raw)
}

clean_tag_id = function(tag_id){
  ## Returns tag ID number as a factor, removing the 'A69-####-' prefix
  cleaned_id = as.factor(substring(tag_id, 10, ))
  return (cleaned_id)
}

clean_receiver = function(receiver){
  ## Returns receiver number as a factor, remvoing the 'VR2W-' Prefix
  cleaned_receiver = as.factor(substring(receiver, 6))
  return (cleaned_receiver)
}

load_vemco = function(filename, filepath = FALSE){
  get_vemco(filename, filepath)
}

### load_receiver_data -> Importing Receiver Data files ---------------------
get_receiver_data = function(filename, filepath = FALSE){
  proj_dir = getwd()
  if (filepath != FALSE){
    setwd(filepath)
  }
  receiver_dates = receiver_col_names(read.csv(filename))
  receiver_dates$deployment_date = strptime(receiver_dates$deployment_date, format = '%m/%d/%y %H:%M', tz = 'HST')
  receiver_dates$recovery_date = strptime(receiver_dates$recovery_date, format = '%m/%d/%y %H:%M', tz = 'HST')
  receiver_dates$lat = convert_lat_lon(receiver_dates$lat_deg, receiver_dates$lat_min)
  receiver_dates$lon = convert_lat_lon(receiver_dates$lon_deg, receiver_dates$lon_min)
  setwd(proj_dir)
  return (receiver_dates)
}

receiver_col_names = function(receiver_file_raw){
  receiver_file = receiver_file_raw
  colnames(receiver_file)[1] = 'serviced'
  colnames(receiver_file)[2] = 'station_name'
  colnames(receiver_file)[3] = 'consecutive_deployment_number'
  colnames(receiver_file)[4] = 'deployment_date'
  colnames(receiver_file)[5] = 'recovery_date'
  colnames(receiver_file)[6] = 'downloaded'
  colnames(receiver_file)[7] = 'in_data_set'
  colnames(receiver_file)[8] = 'lat_deg'
  colnames(receiver_file)[9] = 'lat_min'
  colnames(receiver_file)[10] = 'lon_deg'
  colnames(receiver_file)[11] = 'lon_min'
  colnames(receiver_file)[12] = 'depth'
  colnames(receiver_file)[13] = 'vr2w_serial'
  colnames(receiver_file)[14] = 'acoustic_release_serial'
  colnames(receiver_file)[15] = 'acoustic_release_battery_life'
  colnames(receiver_file)[16] = 'acoustic_release_voltage_at_deployment'
  colnames(receiver_file)[17] = 'acoustic_release_serial_code'
  colnames(receiver_file)[18] = 'temperature_logger_serial'
  colnames(receiver_file)[19] = 'location_code'
  colnames(receiver_file)[20] = 'deployed_by'
  colnames(receiver_file)[21] = 'recovered_by'
  colnames(receiver_file)[22] = 'comments_deployment'
  colnames(receiver_file)[23] = 'comments_recovery'
  return (receiver_file)
}

load_receiver = function(filename, filepath = FALSE){
  return(get_receiver_data(filename, filepath))
}

### get_tagging_data -> Importing Tagging Data ---------------------
get_tagging_data = function(filename, filepath = FALSE){
  proj_dir = getwd()
  if (isTRUE(filepath != FALSE)){
    setwd(filepath)
  }
  tagging_data_raw = read.csv(filename)
  tagging_meta_data = meta_data_col_names(tagging_data_raw)
  tagging_meta_data$datetime = strptime(tagging_meta_data$datetime, format = '%m/%d/%y %H:%M', tz = 'HST')
  tagging_meta_data$vem_tag_id = as.factor(tagging_meta_data$vem_tag_id)
  tagging_meta_data$lat = convert_lat_lon(tagging_meta_data$lat_deg, tagging_meta_data$lat_min)
  tagging_meta_data$lon = convert_lat_lon(tagging_meta_data$lon_deg, tagging_meta_data$lon_min)
  setwd(proj_dir)
  return (tagging_meta_data)
}

meta_data_col_names = function(data_frame){
  colnames(data_frame)[1]  <- 'unique_id'
  colnames(data_frame)[2]  <- 'datetime'
  colnames(data_frame)[3]  <- 'species'
  colnames(data_frame)[4]  <- 'conventional_tag_id'
  colnames(data_frame)[5]  <- 'vem_tag_type'
  colnames(data_frame)[6]  <- 'vem_tag_serial'
  colnames(data_frame)[7]  <- 'vem_tag_id'
  colnames(data_frame)[8]  <- 'fork_length(cm)'
  colnames(data_frame)[9] <- 'precaudal_length(cm)'
  colnames(data_frame)[10] <- 'cohort'
  colnames(data_frame)[11] <- 'area_of_capture'
  colnames(data_frame)[12] <- 'lat_deg'
  colnames(data_frame)[13] <- 'lat_min'
  colnames(data_frame)[14] <- 'lon_deg'
  colnames(data_frame)[15] <- 'lon_min'
  colnames(data_frame)[16] <- 'stomach_everted'
  colnames(data_frame)[17] <- 'eyes_popped'
  colnames(data_frame)[18] <- 'bladder_vented'
  colnames(data_frame)[19] <- 'point_of_incision'
  colnames(data_frame)[20] <- 'dna_clip'
  colnames(data_frame)[21] <- 'cannulation'
  colnames(data_frame)[22] <- 'sex'
  colnames(data_frame)[23] <- 'video'
  colnames(data_frame)[24] <- 'photo'
  colnames(data_frame)[25] <- 'photo_name'
  colnames(data_frame)[26] <- 'audio_log_file'
  colnames(data_frame)[27] <- 'dropshot'
  colnames(data_frame)[28] <- 'tissue_sample'
  colnames(data_frame)[29] <- 'gut_sample'
  colnames(data_frame)[30] <- 'tagger'
  colnames(data_frame)[31] <- 'notes'
  colnames(data_frame)[32] <- 'recaptured'
  colnames(data_frame)[33] <- 'detections'
  colnames(data_frame)[34] <- 'comments'
  return(data_frame)
}

## Usage
load_tagging_data = function(filename, filepath = FALSE){
  get_tagging_data(filename, filepath)
}

### clean_vue_lat_lon -> Fixing VUE data file receiver metadata-----
clean_vue_lat_lon = function(vue_data_df, receiver_data_df){
  station = rep(NA, times = length(vue_data_df$datetime))
  for (i in 1:length(receiver_data_df$station_name)){
    receiver_subset_index = which(vue_data_df$receiver == receiver_data_df$vr2w_serial[i])
    deploy_subset_index = which(vue_data_df$datetime >= receiver_data_df$deployment_date[i])
    recover_subset_index = which(vue_data_df$datetime < na.omit(receiver_data_df$recovery_date[i]))
    ind = Reduce(intersect, list(receiver_subset_index, deploy_subset_index, recover_subset_index))
    vue_data_df$lat[ind] = receiver_data_df$lat[i]
    vue_data_df$lon[ind] = receiver_data_df$lon[i]
    station[ind] = as.character(receiver_data_df$station_name[i])
  }
  vue_data_df$station = as.factor(station)
  return(vue_data_df)
}
    
### remove_location -> Removing station location from database --------------
remove_location = function(vue_data, location_to_remove = FALSE){
  if (location_to_remove == FALSE){
    return (vue_data)
  }else{
    keep_data = vue_data[vue_data$station != location_to_remove, ]
    return (keep_data)}
}

### clean_vue -> Removing database enteries not associated with tag ids ------------------
clean_vue = function(vue_data, bottomfish_tag_ids = FALSE, exclude = FALSE){
  ##Function for removing all tags from a vue DB not explicitly kept by the tag_ids input.
  if (bottomfish_tag_ids[1] == FALSE){
    return (vue_data)}
  if (exclude == FALSE){
  keep_data = vue_data[vue_data$tag_id %in% bottomfish_tag_ids, ]
  return(keep_data)
  } else if (exclude == TRUE){
    keep_data = vue_data[!(vue_data$tag_id %in% bottomfish_tag_ids), ]
    return(keep_data)
  }
}

### clean_tags -> Removes tag ids that are not in the database-----
clean_tags = function(bottomfish_tag_ids, vue_data){
  ## function to remove tags from tag id list that are not in vue database
  keep_tags = bottomfish_tag_ids[bottomfish_tag_ids %in% vue_data$tag_id]
  print (paste('The Following Tags Were Not In the VUE Data:', as.character(bottomfish_tag_ids[!(bottomfish_tag_ids %in% vue_data$tag_id)]),sep = ' '))
  return(keep_tags)
}

### generate_tagging_detection -> for each tag in dataset, -----------
###   adds a detection corrosponding to the tagging location data
generate_tagging_detection = function(tagging_data_df, vue_df){
  tag_ids = na.exclude(tagging_data_df$vem_tag_id[tagging_data_df$vem_tag_id %in% vue_df$tag_id])
  for (i in 1:length(tag_ids)){
    false_record     = vue_df[1, ]
    false_record[1]  = as.POSIXct(as.character(tagging_data_df$datetime[which(tagging_data_df$vem_tag_id == tag_ids[i])]))
    false_record[2]  = 0
    false_record[3]  = tag_ids[i]
    false_record[4]  = NA
    false_record[5]  = NA
    false_record[6]  = NA
    false_record[7]  = 1
    false_record[8]  = 'Tagging Location'
    false_record[9]  = convert_lat_lon(tagging_data_df$lat_deg[which(tagging_data_df$vem_tag_id == tag_ids[i])], tagging_data_df$lat_min[which(tagging_data_df$vem_tag_id == tag_ids[i])])
    false_record[10] = convert_lat_lon(tagging_data_df$lon_deg[which(tagging_data_df$vem_tag_id == tag_ids[i])], tagging_data_df$lon_min[which(tagging_data_df$vem_tag_id == tag_ids[i])])
    vue_data = rbind(false_record, vue_data)
    }
  # Ordering data by date, then tag, then receiver
  vue_data = vue_data[order(vue_data$datetime, vue_data$tag_id, vue_data$receiver), ]
  return (vue_data)
  }

### convert_lat_lon -> Converting Latitude and Logitude ----------------
convert_lat_lon = function(ll_deg, ll_min = FALSE){
  ## Converts latitude and longitude between ll minutes and ll decimal degrees
  # 2 usages:
  # Convert decimal degrees to degree minutes
  # 1 argument
  # ll_pref is a single argument of latitude or longitude in decimal degrees
  # Returns a prefix and decimal for that argument
  # Convert degree minutes to decimal degrees
  # 2 arguments
  # ll_pref is the latitude or longitude's degree
  # ll_min is the degree minutes
  # returns a single float of ll in decimal degrees
  if (ll_min[1] == FALSE){ #then we are going from one number to two
    ll_bin = matrix(0, length(ll_deg), 2)
    for (r in 1:length(ll_deg)){
      if (isTRUE(ll_deg[r] >= 0)){
        ll_dec = ll_deg[r] - floor(ll_deg[r])
        ll_bin[r, ] = c(floor(ll_deg[r]), (ll_dec)*60)
      } else {
        ll_dec = (ll_deg[r] - ceiling(ll_deg[r]))*-1
        ll_bin[r, ] = c(ceiling(ll_deg[r]), (ll_dec)*60)
      }
    }
  }else{ #if we are converting from two numbers to one
    ll_bin = matrix(0, length(ll_deg), 1)
    for (r in 1:length(ll_deg)){
    ll_dec_deg = abs(ll_deg[r]) + (abs(ll_min[r])/60)
    if (isTRUE(ll_deg[r] < 0)){
      ll_dec_deg = ll_dec_deg*(-1)
    }
    ll_bin[r] = ll_dec_deg
  }
  }
return (ll_bin)
}

### lldist -> Getting distance (km) between -------------------------
lldist <- function(point1, point2){
  #The following program computes the distance on the surface of the earth 
  # between two points point1 and point2 in KM
  # Both the points are of the form (Longitude, Latitude)
  #From: http://www.biostat.umn.edu/~sudiptob/Software/distonearth.R
  R <- 6371
  p1rad <- point1 * pi/180
  p2rad <- point2 * pi/180
  d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))  
  d <- acos(d)
  R*d
}

### Plotting Colors -------------------------------
zissou.teal   = wes_palette("Zissou", 5)[1]
zissou.blue   = wes_palette("Zissou", 5)[2]
zissou.gold   = wes_palette("Zissou", 5)[3]
zissou.yellow = wes_palette("Zissou", 5)[4]
zissou.red    = wes_palette("Zissou", 5)[5]


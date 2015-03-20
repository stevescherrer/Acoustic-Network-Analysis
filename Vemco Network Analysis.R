#### Vemco Network Analysis ------------------------------------
#### Written 19 Mar 2015 by Steve Scherrer and Greg Burgess
#### Last Modified: 19 March 2015
#### Email: scherrer@hawaii.edu, gburgess@hawaii.edu

#### R script for analysis of Vemco acoustic data using network
  #### Analysis for identification of homerange and migratory
  #### pathways

#### Clearing Workspace and Assigning Working Directory ------------
rm(list=ls())
setwd('/Users/stephenscherrer/Documents/Work/UH/Projects/dissertation work/Spacial Ecology/Acoustic-Network-Analysis/')

#### Installing Principle Dependancies -----------------------------
#install.packages('dplyr')
library('dplyr')
#install.packages('wes_anderson')
library('wesanderson') # Functions: wes_palette
source('/Users/stephenscherrer/Documents/Programming/R/utility_functions.R')
#install.packages('igraph')
library('igraph')

#### Loading and Cleaning Datafiles --------------------------------
## Importing receiver data
receiver_data = load_receiver(filename = 'DEPLOYMENT_RECOVERY_LOG.csv', 
                              filepath = '/Users/stephenscherrer/Dropbox/Lab Folder/Oahu Receiver Data Files/')

## Import tagging log
tagging_data = load_tagging_data(filename = 'Bottomfish_Tag_Master.csv',
                                 filepath = '/Users/stephenscherrer/Dropbox/Lab Folder/Oahu Receiver Data Files/')

## Importing VUE Data
vue_data = load_vemco(filename = 'VUE_Export_2015-Jan-20.csv',
                      filepath = '/Users/stephenscherrer/Dropbox/Lab Folder/Oahu Receiver Data Files/')

#### Cleaning data -------------------------------------------------

## Fixing missing Lat and Lons
vue_data = clean_vue_lat_lon(vue_data, receiver_data)

## Removing detections from botcam crew
vue_data = remove_location(vue_data, 'With BotCam Crew')

## Removing tags in bottomfish tagging database that are not, 
  ## or likely not actual tagged bottomfish (maybe monchong?)


# Pulling all tag ids associated with opakapaka species
tag_ids = na.exclude(tagging_data$vem_tag_id[
  tagging_data$species == 'Opakapaka']) 

# Converting tag ids from factor to numerics
tag_ids = unique(as.numeric(levels(tag_ids))[tag_ids])

## Removing Vue data for tags not associated with tag_ids of interest
vue_data = clean_vue(vue_data, tag_ids, exclude = FALSE)

## Removing vue data from unwanted tags
# unwanted_tags = c()
# vue_data = clean_vue(vue_data, unwanted_tags, exclude = TRUE)
  
## Removing botcam detections from botcam drops because receivers 
  ## were only deployed for 45 min.
vue_data = remove_location(vue_data, 'With BotCam Crew')


#### Creating Networks -------------------------------------------

get_adjacency = function(vue_df, 
                         tag_id = FALSE, 
                         time_period = FALSE){
  ### A function to create an adjacency matrix from a vue dataset
  ### Arguments: 
    ### vue_df = a dataframe containing a vue export with 
      ### columns for datetime, station, and tag_id. Imported
      ### with load_vemco function from 'R_Utility_Functions.R' 
      ### script
    ### tag_id = a tag id (or multiple ids) corrosponding to the
      ### unique transmitter number for a vemco tag (or tags)
    ### time_period = a list containing two POSIXct objects 
      ### corrosponding to the begining of the time for analysis
      ### and the up to but not including the end time for analysis
  
  ### TO DO: Add functionality to determine which stations had 
  ### receivers during a given time_period. Replace vue_df$station
  ### with these stations.
  
  if (tag_id == FALSE){ # if no tag ids specified, all tag ids used
    tag_id = unique(vue_df$tag_id)
  }
  if (time_period == FALSE){ # if no dates specified, all dates used
    time_period = c(min(vue_df$datetime), max(vue_df$datetime))
  }
  # Pulling out detections for a specific tag
  vue_df = clean_vue(vue_data, tag_id, exclude = FALSE)
  
  # To Do: pull out receivers present during a specific time period
  
  # Build adjacency matrix, with receivers as nodes and fish 
    # movements as edges
  adj_matrix = matrix(0, length(unique(vue_df$station)), 
                      length(unique(vue_df$station)))
  colnames(adj_matrix) = to
  rownames(adj_matrix) = from
  # If station changes, increase adjacency matrix value by one
  for (i in 2:length(vue_df$station)){
    if (vue_df$station[i] != vue_df$station[i-1]){
      adj_matrix[as.numeric(vue_df$station[i-1]), 
                 as.numeric(vue_df$station[i])] = 
        adj_matrix[as.numeric(vue_df$station[i-1]), 
                   as.numeric(vue_df$station[i])] + 1
    }
  }
  # Convert adjacency matrix into a graph object
  vemco_graph = graph.adjacency(adj_matrix, mode = 'directed', 
                                weighted = TRUE)
  # Return graph object
  return(vemco_graph)   
}

#### Make graphs for each fish over appropriate time interval
  #### Monthly? Total study durration for such limited data?

#### Plot graphs

#### Test for nonrandom associations of receivers 
  #### (see Mourier et al 2012, Bejder et al 2005)

#### Calculate community membership (group number of the community/
  #### cluster in the network) from observed matrix for group size
  #### and number of communities in network 
    #### CAN THIS BE DONE USING GRAPH PARTITIONING INSTEAD??

#### Determine if receiver association are significantly different 
  #### from random using coefficient of variation and likelihood 
  #### ratio testing
    #### Note: if normality assumption violated, use non-parametric

#### Monthly (or weekly?) graph for each fish. Compare centrality 
  #### metrics for each receiver. Centrality metrics indicate Core
  #### Use Areas (CUR)
    ## Node Strength: total number of incoming/outgoing movements
      ## (Barrat et al 2004)
    ## Closeness: How central a receiver's position in network space
      ## (Urban et al 2009)
    ## Eigenvector Centrality: how strategically placed a receiver
      ## is. high value, has high node strength (Bodin et all 2011)
  
#### Use structural equivilance graph to identify receivers that are
  #### similar across centrality metrics. These are CURs

#### Principle Components Analysis (PCA) to determine which metric 
  #### to identify CUR.

#### Spearman correlation analysis to remove collinearity. Centrality
  #### metrics where rho > 0.75 removed according to:
    #### eigenvector centrality > node strength > closeness
    #### with metric of lowest strength removed.

#### PCA remaining centrality metrics to determine most influential
  #### in determining network shape. Only PCA with values > 1 used
  #### and only keep those that account for 80% of variance.

#### Centrality metrics with highest absolute loading values retained.

#### Structural equivilance to select receivers with similar 
  #### characteristics to core receivers get added to core group
  #### Core use receivers only identified for time periods where enough
  #### data to produce structural equivilance graph.

#### Percentage approach: Identify CUR where 50% of movements occur
  #### Select CUR one at a time based on number of detections at and
  #### movements to a receiver starting at receiver with highest number 
  #### movements. Select until 50% total movements reached.
#### Select until 95% of total movements reached. This gives General
  #### Use receivers. 

#### Receiver removal analysis of CUR
  #### Method 1: independent removal of each CUR from network
  #### Method 2: successive removal
  #### Then construct new network after each removal. 
  #### Visually examine network and metrics
  #### Ave path length: measure of ease of movement
    #### between pairs of receivers, aka how many receivers an
    #### individual passes when moving between any 2 receivers. Low ave
    #### path length indicates individual passed through few receivers.
    #### (Rayfield et al 2011)
  #### Density: Measures route selection. (ranges between 0 and 1). 
    #### When all receivers connected to all others, density = 1. 
    #### Individual has more routes to choose with larger density
    #### (Rayfield et al 2011)
  #### Components: number of subnetworks or isolates. Represents
    #### network fragmentation. (Bodin et al 2011)
  #### Compare network metrics before and after removal with 
    #### Mann - Whitney U test. if removal does not result in signif
    #### changes to network, receiver is not "Core Use"

#### Number and Frequency of routes between two receivers taken
  #### by an individual. classify as one-way or two-way. 
  #### Classify by frequency of use (high use, low use). With such
  #### poor data, likely to be mostly low use. 


      

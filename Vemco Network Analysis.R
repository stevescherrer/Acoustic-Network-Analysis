#### Vemco Network Analysis ------------------------------------
#### Written 19 Mar 2015 by Steve Scherrer and Greg Burgess
#### Last Modified: 19 March 2015
#### Email: scherrer@hawaii.edu, gburgess@hawaii.edu

#### R script for analysis of Vemco acoustic data using network
  #### Analysis for identification of homerange and migratory
  #### pathways

#### Clearing Workspace and Assigning Working Directory ------------
rm(list=ls())
numStations = 22
user = function(username = FALSE){
  if (username == FALSE){
    print('Please enter a user name (Greg/Steve)')
  }
    if (username == 'Steve' || username == 'steve'){
      setwd('/Users/stephenscherrer/Documents/Work/UH/Projects/dissertation work/Spacial Ecology/Acoustic-Network-Analysis/')
      }
    if (username == 'Greg' || username == 'greg'|| username == 'g'){
      setwd('D:/Workspace/Acoustic-Network-Analysis/')
    }
  }

#### Installing Principle Dependancies -----------------------------

library('dplyr')
library('wesanderson') # Functions: wes_palette
source('utility_functions.R')
library('igraph')
library('intergraph')
#### Loading and Cleaning Datafiles --------------------------------

### Loading data from saved dataframes
load('vue_data.Rda')
load('receiver_data.Rda')

load_all_data = function(){
  #### YOU MUST NEVER RUN THIS GREG. OR YOU WILL BE HAUNTED. BY A SEXY
  #### GHOST. YOU THINK ITS COOL, BUT SHE'S A GHOST. THATS MESSED UP.
  #### ALSO, SHE HAS A GHOST PENIS...
  
  ### Constructing dataframes 
  ## Importing receiver data
  receiver_data = load_receiver(filename = 'DEPLOYMENT_RECOVERY_LOG.csv', 
                                filepath = '/Users/stephenscherrer/Dropbox/Lab Folder/Oahu Receiver Data Files/')
  
  ## Import tagging log
  tagging_data = load_tagging_data(filename = 'Bottomfish_Tag_Master.csv',
                                   filepath = '/Users/stephenscherrer/Dropbox/Lab Folder/Oahu Receiver Data Files/')
  
  ## Importing VUE Data
  vue_data = load_vemco(filename = 'VUE_Export_2015-Jan-20.csv',
                        filepath = '/Users/stephenscherrer/Dropbox/Lab Folder/Oahu Receiver Data Files/')
  
}
  #### Cleaning data -------------------------------------------------
  processData = function() {
  ## Fixing missing Lat and Lons
  vue_data = clean_vue_lat_lon(vue_data, receiver_data)
  
  ## Removing detections from botcam crew
  vue_data = remove_location(vue_data, 'With BotCam Crew')
  
  ## Removing tags in bottomfish tagging database that are not, 
    ## or likely not actual tagged bottomfish (at cross 
    ## maybe monchong?)
   vue_data = clean_vue(vue_data, c(57451, 37954), exclude = TRUE)
  
  # Pulling all tag ids associated with opakapaka species
  #tag_ids = na.exclude(tagging_data$vem_tag_id[
    #tagging_data$species == 'Opakapaka']) 
  
  ## Converting tag ids from factor to numerics
  #tag_ids = unique(as.numeric(levels(tag_ids))[tag_ids])
  
  ## Pulling all tag ids associated with white shark species
  # Tag IDs for White Sharks
  tag_ids = c(633, 68295, 8893, 8832, 1425, 29700, 2139, 8892, 9095, 
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
  
  ## Removing Vue data for tags not associated with tag_ids of interest
  vue_data = clean_vue(vue_data, tag_ids, exclude = FALSE)
  
  ## Removing vue data from unwanted tags
  # unwanted_tags = c()
  # vue_data = clean_vue(vue_data, unwanted_tags, exclude = TRUE)
    
  ## Removing botcam detections from botcam drops because receivers 
    ## were only deployed for 45 min.
  vue_data = remove_location(vue_data, 'With BotCam Crew')
  
  ## Saving Data Frames
    save(vue_data, file = 'vue_data.Rda')
    save(receiver_data, file = 'receiver_data.Rda')
 }

 ## Pulling only those tag IDS that appear in the data
 tag_ids = unique(vue_data$tag_id)

#### Creating Networks -------------------------------------------

get_graph = function(vue_df, 
                     tag_id = FALSE, 
                     time_period = FALSE,
                     igraph=TRUE,
                     binary=FALSE,
                     removeLoops=FALSE){
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
  
  # Order vue_df by tag id
    vue_df = vue_df[order(vue_df$tag_id, vue_df$datetime), ]
  
  # To Do: pull out receivers present during a specific time period
  
  # Build adjacency matrix, with receivers as nodes and fish 
    # movements as edges
    # Rows indicate movement from a receiver
    # Columns indicate movement to a receiver
  adj_matrix = matrix(0, numStations,numStations)
  # If station changes, increase adjacency matrix value by one
  for (i in 2:length(vue_df$station)){
    if(vue_df$tag_id[i] == vue_df$tag_id[i-1]){
      prevLoc = as.numeric(vue_df$station[i-1])
      newLoc = as.numeric(vue_df$station[i])
      if(binary) {
        adj_matrix[prevLoc, newLoc] = 1
      }
      else {
        adj_matrix[prevLoc, newLoc] =  adj_matrix[prevLoc, newLoc] + 1
      }
    }
  }
  if (removeLoops) {
    for (i in 1:numStations) {
      adj_matrix[i,i]=0
    }
  }
  if(igraph) {
    # Convert adjacency matrix into a graph object
    vemco_graph = graph.adjacency(adj_matrix, mode = 'directed', 
                                  weighted = TRUE)
    return (vemco_graph)
  }
  # Return graph object
  return(adj_matrix)   
} 

subset_time = function(vue_data, start = min(vue_data$datetime), end = max(vue_data$datetime)){
  new_vue_df = vue_data[which(vue_data$datetime >= as.POSIXct(start) & vue_data$datetime < as.POSIXct(end)), ]
  return (new_vue_df)
}


subset_tag = function(vue_data, tag_id){
  new_vue_df = vue_data[as.character(vue_data$tag_id) %in% as.character(tag_id), ]
  return (new_vue_df)
}

station_ids = function(vue_data){
  vue_data = vue_data[order(vue_data$station), ]
  station_number = rep(0, length(vue_data$station))
  for (i in 2:length(vue_data$station)){
    if (vue_data$station[i] != vue_data$station[i-1]){
      station_number[i] = station_number[i-1] + 1
    } else {
      station_number[i] = station_number[i-1]
    }
  }
  vue_data$station_number = station_number+1
  vue_data = vue_data[order(vue_data$tag_id, vue_data$datetime), ]
  return(vue_data)
}

station_ids_map = function(vue_data){
  vue_data = vue_data[order(vue_data$station), ]
  station_number = as.data.frame(matrix(0, length(vue_data$station), 2))
  colnames(station_number) = c('Station Number', 'Station Name')
  for (i in 2:length(vue_data$station)){
    if (vue_data$station[i] != vue_data$station[i-1]){
      station_number[i, 1] = station_number[i-1, 1] + 1
    } else {
      station_number[i, 1] = station_number[i-1, 1]
    }
    station_number[i,2] = as.character(vue_data$station[i])
  }
  station_number[ ,1] = station_number[ ,1] + 1
  station_number = station_number[2:length(station_number$'Station Number'), ]
  station_code_map = unique(station_number)
  return(station_code_map)
}



anlysis1 = function() {
  #for each fish, for each time delta, graph the fish movement,
  #convert a weighted matrix to a binary matirx
  #add up all the binary matricies (across the time deltas), and see which ones have the highest values
  
  # Create 3d matrix [tagID, timedelta, reciever, receiver]
  # for each fish
  #  for each time delta
  #    create adj matrix
  # total += matrix we just made
}


anlysis2 = function() {
  #How important is each node?
  
  #results = dict
  # for each node in the graph:
    # temp = delete.vertices(graph, v) //remove a node v
    # results[node] = no.clusters(temp)  //counts the number of isolates
  #print results
}


analysis3 = function() {
  #Kernel based approach????
}


analysis4 = function() {
  #How important is each edge?
  #4 and 18 are the top left most nad bottom right most stations.
  
  #F = max flow from a to b
  #k = # of ways to get from a to b
  #F/k = importance of an edge
  # F = graph.maxflow(graph, source, target, capacity=NULL) //gives max flow from source to target
  # k = vertex.connectivity(graph, source=NULL, target=NULL, checks=TRUE)
  # F/k = importance of edge from a to b
  # edges with high F/k are important
  
  graph = get_graph(vue_data,igraph=T,binary=F,removeLoops=F)
  result = matrix(,numStations,numStations);
  for(source in 1:numStations) {
    for(target in 1:numStations) {
      if(source!=target) {
        flow = graph.maxflow(graph, source, target, capacity=NULL)
        numPaths = vertex.connectivity(graph, source=NULL, target=NULL, checks=TRUE)
        if(numPaths> 0) {
          result[source,target] = flow$value/numPaths
        }
        else {
          result[source,target] = 0
        }
      }
    }
  }
  return (result)
}

analysis5 = function() {
  # Which paths are imprtant to a species?
  # for each fish, make a binary matrix for whole experiment
  # sum all fish matricies
  # divide resulting matrix by # of fish
  #//low numbers mean a path is only used by some individuals
  #//high numbers mean a path is important to a species.
  result = matrix(,numStations,numStations);
  tagIds = unique(vue_data$tag_id)
  first = TRUE
  for (tagId in tagIds) {
    #mat = make a binary matrix for the whole experiment
    mat = get_graph(vue_data, tag_id=tagId,igraph=FALSE,binary=TRUE)
    if (first) {
      result = matrix
      result = mat
      first=FALSE
    }
    else {
      result = result + mat
    }
  }
  return (result)
}

analysis6 = function() {
#community detection
#steve has code for this
  node_groupings = function(graph_obj){
    bf = fastgreedy.community(as.undirected(graph_obj)) # blondel et al ocmmnity.
    #summary(bf)
    plot(bf, graph_obj)  
    library(ape)
    dendPlot(bf, mode = 'phylo')
    bflap = graph.laplacian(graph_obj)
    eig.anal = eigen(bflap)
    plot(eig.anal$values, col = zissou.blue, ylab = 'eigenvalue of graph laplacian')
    faction = get.vertex.attribute(graph_obj, 'Community')
    f.colors = as.character(length(faction))
    #plot(faction, xlab = 'actor number', ylab = 'fiedler vector entry', col = f.colors)
    #abline(0,0,lwd = 2)
  }

}










#### Analysis from other paper ----------------------------
#### Make graphs for each shark over appropriate time interval
  #### Yearly?

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


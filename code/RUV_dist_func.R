RUVdist.prep <- function(spatial_dat,fish_dat,vid_name,species,site,cam_height,height_scaler,nframes,object_id_per_species = TRUE){

# Function to calculate the plane equation coefficients
calculate_plane <- function(p1, p2, p3) {
  v1 <- p2 - p1    
  v2 <- p3 - p1
  normal <- c(
    v1[2] * v2[3] - v1[3] * v2[2],
    v1[3] * v2[1] - v1[1] * v2[3],
    v1[1] * v2[2] - v1[2] * v2[1]
  ) #normal vector to the plane
  a <- normal[1]
  b <- normal[2]
  c <- normal[3]
  d <- - (a * p1[1] + b * p1[2] + c * p1[3])
  return(c(a, b, c, d))
}

# Function to project a point onto a plane
project_point_onto_plane <- function(point, plane_coefficients) {
  a <- plane_coefficients[1]
  b <- plane_coefficients[2]
  c <- plane_coefficients[3]
  d <- plane_coefficients[4]
  
  # Calculate the normal vector from the point to the plane
  t <- -(a * point[1] + b * point[2] + c * point[3] + d) / (a^2 + b^2 + c^2)
  
  # Calculate the projection
  projected_point <- point + t * c(a, b, c)
  
  return(projected_point)
}

# Function to find distance from a point onto a plane
dist_point_onto_plane <- function(point, plane_coefficients) {
  a <- plane_coefficients[1]
  b <- plane_coefficients[2]
  c <- plane_coefficients[3]
  d <- plane_coefficients[4]
  
  # Calculate the normal vector from the point to the plane
  dist <- abs(a * point[1] + b * point[2] + c * point[3] + d) / sqrt((a^2 + b^2 + c^2))
  
  return(dist)
}

# Function to calculate the angle between three points
calculate_angle <- function(p1, p2, p3) {
  # Vectors from p2 to p1 and p2 to p3
  v1 <- p1 - p2
  v2 <- p3 - p2
  
  # Dot product of v1 and v2
  dot_product <- sum(v1 * v2)
  
  # Magnitudes of v1 and v2
  mag_v1 <- sqrt(sum(v1^2))
  mag_v2 <- sqrt(sum(v2^2))
  
  # Cosine of the angle
  cos_theta <- dot_product / (mag_v1 * mag_v2)
  
  # Angle in radians
  angle_radians <- acos(cos_theta)
  
  # Convert angle to degrees
  angle_degrees <- angle_radians * (180 / pi)
  
  return(angle_degrees)
}

# function to move closer to a plane by a distance d
move_point_closer_to_plane <- function(A, B, C, D, x1, y1, z1, d) {
  # Calculate the magnitude of the normal vector
  normal_magnitude <- sqrt(A^2 + B^2 + C^2)
  
  # Calculate the unit normal vector
  unit_normal <- c(A, B, C) / normal_magnitude
  
  # Calculate the new coordinates
  new_x <- x1 - d * unit_normal[1]
  new_y <- y1 - d * unit_normal[2]
  new_z <- z1 - d * unit_normal[3]
  
  return(c(new_x, new_y, new_z))
}

fish_dat_final = data.frame()

vid_names = unique(spatial_dat[[vid_name]]) 

#loop through each video, calculating plane equation, distance on plane, theta, gamma and height.
for(j in (1:length(vid_names))){
  
  dat_fish_temp = fish_dat %>% filter(fish_dat[[vid_name]] == vid_names[j])
  dat_spatial_temp = spatial_dat %>% filter(spatial_dat[[vid_name]] == vid_names[j])
  
  # ground points
  p1 <- unlist(c(0, 0, -cam_height))
  p2 <- unlist(dat_spatial_temp %>% filter(Species=="BOTTOM") %>%
                 filter(row_number() == 1) %>% dplyr::select(X,Z,Y))
  p3 <- unlist(dat_spatial_temp %>% filter(Species=="BOTTOM") %>%
                 filter(row_number() == 2) %>% dplyr::select(X,Z,Y))
  
  # Calculate the plane equation coefficients
  plane_coefficients <- calculate_plane(p1, p2, p3)
  
  
  #find the height and distance of each fish to the camera along the seafloor
  dat_fish_temp$height = NA
  dat_fish_temp$height_alt = NA
  dat_fish_temp$plane_dist = NA
  dat_fish_temp$original_dist = NA
  

  #find vertical angle gamma (proper way)
  
  #reference point 
  top_point <- dat_spatial_temp %>% filter(Species == "UPPER") %>% dplyr::select(X,Z,Y)
  
  
  if(dim(top_point)[1]==1){
    #projected point
    #find distance of point to plane
    height_temp <- dist_point_onto_plane(top_point, plane_coefficients)/1000
    
    #find distance from the point to the camera along the plane
    projected_point_top <- project_point_onto_plane(unlist(top_point), plane_coefficients)
    
    top_dist_temp <- sqrt(sum((projected_point_top[1:2])^2) + (p1[3] - projected_point_top[3])^2)/1000
    
    # Find vertical angle gamma
    gamma_angle =  atan((height_temp-cam_height/1000)/top_dist_temp)*180/pi
    
    
    dat_fish_temp$gamma = unlist(gamma_angle)
    
  }else{
    dat_fish_temp$gamma = NA
    }
  
  
  #find horizontal angle theta
  
  #reference point 
  left_point <- dat_spatial_temp %>% filter(Species == "LEFT") %>% dplyr::select(X,Z,Y)
  right_point <- dat_spatial_temp %>% filter(Species == "RIGHT") %>% dplyr::select(X,Z,Y)
  
  #projected point
  if(dim(left_point)[1]==1 & dim(right_point)[1]==1){
    projected_point_left <- project_point_onto_plane(unlist(left_point), plane_coefficients)
    projected_point_right <- project_point_onto_plane(unlist(right_point), plane_coefficients)
    
    #horizontal angle 
    theta_angle <- calculate_angle(unlist(projected_point_left), p1, unlist(projected_point_right))
    theta_angle
    
    dat_fish_temp$theta = theta_angle
  }else{
    dat_fish_temp$theta = NA
  }
  
  for (i in 1:nrow(dat_fish_temp)){
    new_point = c(dat_fish_temp$X[i],dat_fish_temp$Z[i],dat_fish_temp$Y[i])
    
    # Project the point onto the plane
    projected_point <- project_point_onto_plane(new_point, plane_coefficients)
    
    #find the distance (m) along the plane to this point
    dat_fish_temp$plane_dist[i] = sqrt(sum((projected_point[1:2])^2) + (p1[3] - projected_point[3])^2)/1000
    
    dat_fish_temp$original_dist[i] = sqrt(sum(new_point^2))/1000
    
    #find alternative height metric
    dat_fish_temp$height_alt[i] = sqrt(sum((new_point-projected_point)^2))/1000
    
    #perpendicular distance (m) from the point to the plane
    dat_fish_temp$height[i] <- dist_point_onto_plane(new_point, plane_coefficients)/1000
    
    #calculate theta angle of fish to left and right point
    dat_fish_temp$theta_angle_L[i] <- calculate_angle(unlist(projected_point_left), p1, unlist(projected_point))
    dat_fish_temp$theta_angle_R[i] <- calculate_angle(unlist(projected_point), p1, unlist(projected_point_right))
    
  }
  
  fish_dat_final <-rbind(fish_dat_final,dat_fish_temp)
}

#distribution of thetas
min_theta = min(fish_dat_final$theta,na.rm=T) #51.99647
min_gamma = min(fish_dat_final$gamma,na.rm=T)

#filter out videos with no angles
fish_dat_final = fish_dat_final %>% filter(is.na(gamma)==FALSE,is.na(theta)==FALSE)

#now filter out fish with a horizontal viewing angle greater then min_theta

#first check if any fish were observed outside of the viewing window
#fish_dat_final$theta_angle_L[fish_dat_final$theta_angle_L<0]
#fish_dat_final$theta_angle_R[fish_dat_final$theta_angle_R<0]

#fish_dat_final$theta_angle_L[fish_dat_final$theta_angle_L>90]
#fish_dat_final$theta_angle_R[fish_dat_final$theta_angle_R>90]

fish_dat_final$theta_angle_LR_min = apply(data.frame(fish_dat_final$theta_angle_R,fish_dat_final$theta_angle_L),c(1),min)

fish_dat_final$theta_cut_off = (fish_dat_final$theta-min_theta)/2

#dat_fish_pre_filt = fish_dat_final
fish_dat_final = fish_dat_final %>% filter(theta_angle_LR_min>theta_cut_off)

#lose 341 fish...
#dim(WA_dat_fish_pre_filt)[1]-dim(fish_dat_final)[1]



#to get truncations
height_buffer_dat = fish_dat_final %>% group_by(fish_dat_final[[species]]) %>% summarise(height_buffer = max(height)*height_scaler)

colnames(height_buffer_dat)[1] = species

fish_dat_final <- left_join(fish_dat_final,height_buffer_dat,by=species)

fish_dat_final$left_trunc_gamma_min =  (fish_dat_final$height_buffer-cam_height/1000)/tan(min(fish_dat_final$gamma)*pi/180)

#update to the smallest values possible assuming angle down is the same as angle up
fish_dat_final$left_trunc_gamma_min[fish_dat_final$left_trunc_gamma_min<(cam_height/1000)/tan(min(fish_dat_final$gamma)*pi/180)] = (cam_height/1000)/tan(min(fish_dat_final$gamma)*pi/180)


fish_dat_final <- fish_dat_final %>% filter(plane_dist>left_trunc_gamma_min)


fish_dat_final$distance = fish_dat_final$plane_dist


distance_dat = data.frame(Region.Label=fish_dat_final[[site]],Area=0,
                             Sample.Label=fish_dat_final[[vid_name]],distance=fish_dat_final$distance,
                             Effort=nframes,species=fish_dat_final[[species]],object=as.factor(1:length(fish_dat_final[[vid_name]])),
                          left_trunc=fish_dat_final$left_trunc_gamma_min,theta=fish_dat_final$theta)


video_codes = rbind(spatial_dat,fish_dat) %>%
  dplyr::select(all_of(site),all_of(vid_name)) %>% distinct()

species_list = fish_dat %>%
  dplyr::select(all_of(species)) %>% distinct() %>% c()


distance_data_final = data.frame()

for (j in 1:length(species_list[[1]])){

  distance_data_temp = distance_dat %>% filter(species == species_list[[1]][j])

  if(dim(distance_data_temp)[1]>0){

    if(object_id_per_species){
      distance_data_temp$object = 1:dim(distance_data_temp)[1]
    }

    for (i in 1:dim(video_codes)[1]){
      #condition to only add in NAs for videos where the fish was observed at other videos in that site
      if(distance_data_temp %>% filter(Region.Label==video_codes[[site]][i]) %>% nrow()>0){

        if((video_codes[[vid_name]][i] %in% distance_data_temp$Sample.Label)==FALSE){
          dist_dat_empty = data.frame(Region.Label=video_codes[[site]][i],Area=0,
                                      Sample.Label=video_codes[[vid_name]][i],distance=NA,
                                      Effort=40,species=species_list[[1]][j],
                                      object=NA,left_trunc=NA,theta=NA)
          distance_data_temp = rbind(distance_data_temp,dist_dat_empty)
        }
      }
    }
    distance_data_final = rbind(distance_data_final,distance_data_temp)
  }
}



return(distance_data_final)
}


#update1: Fix lower angle equivalence
#update2: Code to handle the base of the camera stand not being on the plane equation







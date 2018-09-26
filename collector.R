library(ggplot2)
library(dplyr)
library(purrr)

#Function library

############################
#High-Level FFT Functions
############################

#This identifies the three highest values in the first half of an fft, ignoring the first entry.
#Input: a transformed fft
#Output: a dataframe with three entries and coordinates (x,y)
identifyThreeHighestValues <- function(transformed_vector){
  transformed_vector = Re(Mod(transformed_vector))
  transformed_vector = transformed_vector[2:floor((length(transformed_vector)/2))]
  transformed_data = data.frame("foo" =transformed_vector)
  largest = t(apply(transformed_data, 2, order, decreasing=TRUE)[ 1:3, ])[1,]
  data.frame("x"= largest+1, "y"= transformed_data[largest,])
}


# Inputs a data.frame with an "x" and "y" column, with three entries in each, corresponding to three points.
# Outputs a vector with two entries: the "true peak" frequency and the magnitude of said frequency
#If the true peak is not within the interval obtained from the largest three magnitudes, then the raw peak is used instead.
truePeak <- function(data){
  x<-data$x
  y<-data$y
  peakModel<-lm(y~ poly(x,2,raw=TRUE))
  a<-as.data.frame(coef(peakModel))[3,1]
  if (a>=0) {
    return(c(data$x[1],data$y[1])) 
  }
  else{
    b<-as.data.frame(coef(peakModel))[2,1]
    c<-as.data.frame(coef(peakModel))[1,1]
    
    p=-b/(2*a)
    k=c-(b*b)/(4*a)
    #c(a,b,c)
  }
  if(p>min(data$x) && p<max(data$x)){
    return(c(p,k))
  }
  else{
    data_y = data.frame(data[,2])
    index = t(apply(data_y, 2, order, decreasing=TRUE)[1,])[1,]
    as.matrix(data)[index,]
  }
}


#Produces the True Peak data
#Input: a transformed vector
#Output: a numerical vector, [True Frequency, True Magnitude]
getTruePeakVector <- function(transformed_vector){
  truePeak(identifyThreeHighestValues(transformed_vector))
}


#Produces the true frequency
#Input: a transformed vector, an optional window size
#Output: A frequency
TrueFrequency <- function(transformed_vector, windowSize = 4){
  getTruePeakVector(transformed_vector)[1]/windowSize
}

#Produces the true magnitude
#Input: a transformed vector
#Output: A magnitude
TrueMagnitude <- function(transformed_vector){
  getTruePeakVector(transformed_vector)[2]
}


#Produces the energy of a fourier transform with or without DC term
#Input: a transformed vector, the DC term is preset to off
#Output: A number
energy = function(transformed_vector,includeDC=FALSE){
  if(includeDC){
    (Re(sum(transformed_vector*Conj(transformed_vector))))^(1/2)
  }
  else{
    transformed_vector = transformed_vector[-c(1)]
    (Re(sum(transformed_vector*Conj(transformed_vector))))^(1/2)
  }
}

#Computes entropy from a tranformed vector
#Input: Takes in a transformed vector, option for including DC signal
#Output: Returns an entropy value
entropy = function(transformed_vector,includeDC=FALSE){
  if(includeDC){
    transformed_vector=(Re(Mod(transformed_vector)))
  }
  else{
    transformed_vector = transformed_vector[-c(1)]
    transformed_vector=(Re(Mod(transformed_vector)))
  }
  
  #Exception Handling
  if(sum(transformed_vector == 0)){
    return(1)
  }
  probabilities = (transformed_vector/sum(transformed_vector))
  if(any(probabilities == 0)){
    return(1)
  }
  -sum(probabilities*log10(probabilities))
}

###########################
#End of FFT Library
###########################


# Filtering function to desired frequency
filter_time <- function(data,freq){
  # frequency is the desired frequency, taken_freq is the frequency that we used
  times <- data[,1] # the times
  final_time <- times[length(times)]
  initial_time <- times[1]
  taken_freq <- length(times)/(final_time - initial_time)
  n <- 1/freq
  m <- taken_freq*n # approximate distance to next index
  iterations <- floor((final_time - initial_time)/n) # total number of times at the end
  indices <- vector(mode="numeric", length=iterations) # initialize vector
  indices[1] <- 1
  for(i in 1:(iterations-1)){
    actual_indices <- c(max(m*i-10,1):min(m*i+10,length(times)))
    potential_times <- vector(mode="numeric", length=length(actual_indices))
    for(j in 1:length(actual_indices)){
      potential_times[j] <- abs(times[actual_indices[j]]-(initial_time + i*n))
    } # here we made two lists: one with differences to initial_time+i*n, the other with the actual range of indices
    next_time <- min(potential_times) # smallest distance to initial_time+i*n among indices
    is_there <- function(i){potential_times[i] == next_time}
    current_index <- 1:length(potential_times) %>% detect_index(is_there) # finding index so that time is closest to initial_time+i*n
    indices[i+1] <- actual_indices[current_index]
  }   
  new_data <- data[indices,]
  for(i in 2:(iterations)){
    new_data[i,1] <- new_data[1,1] + (i-1)*n
  }
  return(new_data)
}

#Get average of a column
extractavg <- function(dframe, col){
  foo = mean(dframe[,col])
  return(foo)
}

#Get maximum of a column
extractmax <- function(dframe, col){
  foo = max(dframe[,col])
  return(foo)
}

#Get minimum of a column
extractmin <- function(dframe, col){
  foo = min(dframe[,col])
  return(foo)
}


#Get stdev of a column
extractsd <- function(dframe, col){
  foo = sd(dframe[,col])
  return(foo)
}

#Get data (magnitude, significant freq) from the fast fourier transform
extractfft <- function(dframe, col, seglength){
  fft_0 <- fft(dframe[,col])
  tfreq <- TrueFrequency(fft_0)
  tamp <- TrueMagnitude(fft_0)
  en <- energy(fft_0)
  ent <- entropy(fft_0)
  fft_1 <- Mod(fft_0[2:floor(length(fft_0)/2)])
  return(c(max(fft_1), tamp, tfreq, en, ent, which.max(fft_1)/seglength))
  
}

#Set working directory if needed (depends on individual enviornment)
#setwd("C:/Users/SangJoon/Downloads/sensor/Test/")

#Read in the raw data
master_data <- read.csv(file.choose())

#add total linear acc (euclidean norm)
master_data <- master_data %>%
  mutate(aTT = sqrt(ax^2+ay^2+az^2))

#crop the data (seconds)
crop = 5
master_data <- master_data[master_data$time > crop & master_data$time < (max(master_data$time) - crop),]

#filter the data (hz)
filterfreq = 8
filtered_data <- filter_time(master_data,filterfreq)

#Make the output data frame
output_data <- data.frame("Timestamp" = NULL,
                          "gforcexavg" = NULL, "gforcexstdev" = NULL, "gforcexmax" = NULL, "gforcexmin" = NULL,
                          "gforcexmaxamp" = NULL, "gforcexmaxfreq" = NULL, "tgforcexmaxamp" = NULL, "tgforcexmaxfreq" = NULL, "gforcexenergy" = NULL, "gforcexentropy" = NULL,
                          "gforceyavg" = NULL, "gforceystdev" = NULL, "gforceymax" = NULL, "gforceymin" = NULL,
                          "gforceymaxamp" = NULL, "gforceymaxfreq" = NULL, "tgforceymaxamp" = NULL, "tgforceymaxfreq" = NULL, "gforceyenergy" = NULL, "gforceyentropy" = NULL,
                          "gforcezavg" = NULL, "gforcezstdev" = NULL, "gforcezmax" = NULL, "gforcezmin" = NULL,
                          "gforcezmaxamp" = NULL, "gforcezmaxfreq" = NULL, "tgforcezmaxamp" = NULL, "tgforcezmaxfreq" = NULL, "gforcezenergy" = NULL, "gforcezentropy" = NULL,
                          "linaccxavg" = NULL, "linaccxstdev" = NULL, "linaccxmax" = NULL, "linaccxmin" = NULL,
                          "linaccxmaxamp" = NULL, "linaccxmaxfreq" = NULL, "tlinaccxmaxamp" = NULL, "tlinaccxmaxfreq" = NULL, "linaccxenergy" = NULL, "linaccxentropy" = NULL,
                          "linaccyavg" = NULL, "linaccystdev" = NULL, "linaccymax" = NULL, "linaccymin" = NULL,
                          "linaccymaxamp" = NULL, "linaccymaxfreq" = NULL, "tlinaccymaxamp" = NULL, "tlinaccymaxfreq" = NULL, "linaccyenergy" = NULL, "linaccyentropy" = NULL,
                          "linacczavg" = NULL, "linacczstdev" = NULL, "linacczmax" = NULL, "linacczmin" = NULL,
                          "linacczmaxamp" = NULL, "linacczmaxfreq" = NULL, "tlinacczmaxamp" = NULL, "tlinacczmaxfreq" = NULL, "linacczenergy" = NULL, "linacczentropy" = NULL,
                          "linaccavg" = NULL, "linaccstdev" = NULL, "linaccmax" = NULL, "linaccmin" = NULL,
                          "linaccmaxamp" = NULL, "linaccmaxfreq" = NULL, "tlinaccmaxamp" = NULL, "tlinaccmaxfreq" = NULL, "linaccenergy" = NULL, "linaccentropy" = NULL,
                          "speedavg" = NULL, "speedstdev" = NULL, "speedmax" = NULL, "speedmin" = NULL,
                          "gyroxavg" = NULL, "gyroxstdev" = NULL, "gyroxmax" = NULL, "gyroxmin" = NULL, "gyroxmaxamp" = NULL, "gyroxmaxfreq" = NULL, "tgyroxmaxamp" = NULL, "tgyroxmaxfreq" = NULL, "gyroxenergy" = NULL, "gyroxentropy" = NULL,
                          "gyroyavg" = NULL, "gyroystdev" = NULL, "gyroymax" = NULL, "gyroymin" = NULL, "gyroymaxamp" = NULL, "gyroymaxfreq" = NULL, "tgyroymaxamp" = NULL, "tgyroymaxfreq" = NULL, "gyroyenergy" = NULL, "gyroyentropy" = NULL,
                          "gyrozavg" = NULL, "gyrozstdev" = NULL, "gyrozmax" = NULL, "gyrozmin" = NULL, "gyrozmaxamp" = NULL, "gyrozmaxfreq" = NULL, "tgyrozmaxamp" = NULL, "tgyrozmaxfreq" = NULL, "gyrozenergy" = NULL, "gyrozentropy" = NULL,
                          "sampleno" = NULL)


#Determine size of segments (seconds)
segsize = 4

#Determine number of segments
numsegments <- 2*floor((max(filtered_data$time)-crop)/segsize)-2

#For each data segments
for (i in 0:numsegments){
  
  #Extract the segment of raw data
  data_temp <- filtered_data[filtered_data$time>crop+segsize*i/2
                             & filtered_data$time <crop+segsize*(i+2)/2,]
  
  #Build data tuple of statistics of this segment
  ftx <- extractfft(data_temp, "ax", segsize)
  fty <- extractfft(data_temp, "ay", segsize)
  ftz <- extractfft(data_temp, "az", segsize)
  x <- extractfft(data_temp, "aTT", segsize)
  gx <- extractfft(data_temp, "wx", segsize)
  gy <- extractfft(data_temp, "wy", segsize)
  gz <- extractfft(data_temp, "wz", segsize)
  fx <- extractfft(data_temp, "gFx", segsize)
  fy <- extractfft(data_temp, "gFy", segsize)
  fz <- extractfft(data_temp, "gFz", segsize)
  data_tuple <- data.frame("Timestamp" = crop+segsize*(i+1)/2,
                           "gforcexavg" = extractavg(data_temp, "gFx"),
                           "gforcexstdev" = extractsd(data_temp, "gFx"),
                           "gforcexmax" = extractmax(data_temp, "gFx"),
                           "gforcexmin" = extractmin(data_temp, "gFx"),
                           "gforcexmaxamp" = fx[1],
                           "gforcexmaxfreq" = fx[6],
                           "tgforcexmaxamp" = fx[2],
                           "tgforcexmaxfreq" = fx[3],
                           "gforcexenergy" = fx[4],
                           "gforcexentropy" = fx[5],
                           "gforceyavg" = extractavg(data_temp, "gFy"),
                           "gforceystdev" = extractsd(data_temp, "gFy"),
                           "gforceystmax" = extractmax(data_temp, "gFy"),
                           "gforceystmin" = extractmin(data_temp, "gFy"),
                           "gforceymaxamp" = fy[1],
                           "gforceymaxfreq" = fy[6],
                           "tgforceymaxamp" = fy[2],
                           "tgforceymaxfreq" = fy[3],
                           "gforceyenergy" = fy[4],
                           "gforceyentropy" = fy[5],
                           "gforcezavg" = extractavg(data_temp, "gFz"),
                           "gforcezstdev" = extractsd(data_temp, "gFz"),
                           "gforcezstmax" = extractmax(data_temp, "gFz"),
                           "gforcezstmin" = extractmin(data_temp, "gFz"),
                           "gforcezmaxamp" = fz[1],
                           "gforcezmaxfreq" = fz[6],
                           "tgforcezmaxamp" = fz[2],
                           "tgforcezmaxfreq" = fz[3],
                           "gforcezenergy" = fz[4],
                           "gforcezentropy" = fz[5],
                           "linaccxavg" = extractavg(data_temp, "ax"),
                           "linaccxstdev" = extractsd(data_temp, "ax"),
                           "linaccxmax" = extractmax(data_temp, "ax"),
                           "linaccxmin" = extractmin(data_temp, "ax"),
                           "linaccxmaxamp" = ftx[1],
                           "linaccxmaxfreq" = ftx[6],
                           "tlinaccxmaxamp" = ftx[2],
                           "tlinaccxmaxfreq" = ftx[3],
                           "linaccxenergy" = ftx[4],
                           "linaccxentropy" = ftx[5],
                           "linaccyavg" = extractavg(data_temp, "ay"),
                           "linaccystdev" = extractsd(data_temp, "ay"),
                           "linaccymax" = extractmax(data_temp, "ay"),
                           "linaccymin" = extractmin(data_temp, "ay"),
                           "linaccymaxamp" = fty[1],
                           "linaccymaxfreq" = fty[6],
                           "tlinaccymaxamp" = fty[2],
                           "tlinaccymaxfreq" = fty[3],
                           "linaccyenergy" = fty[4],
                           "linaccyentropy" = fty[5],
                           "linacczavg" = extractavg(data_temp, "az"),
                           "linacczstdev" = extractsd(data_temp, "az"),
                           "linacczmax" = extractmax(data_temp, "az"),
                           "linacczmin" = extractmin(data_temp, "az"),
                           "linacczmaxamp" = ftz[1],
                           "linacczmaxfreq" = ftz[6],
                           "tlinacczmaxamp" = ftz[2],
                           "tlinacczmaxfreq" = ftz[3],
                           "linacczenergy" = ftz[4],
                           "linacczentropy" = ftz[5],
                           "linaccavg" = extractavg(data_temp, "aTT"),
                           "linaccstdev" = extractsd(data_temp, "aTT"),
                           "linaccmax" = extractmax(data_temp, "aTT"),
                           "linaccmin" = extractmin(data_temp, "aTT"),
                           "linaccmaxamp" = x[1],
                           "linaccmaxfreq" = x[6],
                           "tlinaccmaxamp" = x[2],
                           "tlinaccmaxfreq" = x[3],
                           "linaccenergy" = x[4],
                           "linaccentropy" = x[5],
                           "speedavg" = extractavg(data_temp, "Speed..m.s."),
                           "speedstdev" = extractsd(data_temp, "Speed..m.s."),
                           "speedmax" = extractmax(data_temp, "Speed..m.s."),
                           "speedmin" = extractmin(data_temp, "Speed..m.s."),
                           "gyroxavg" = extractavg(data_temp, "wx"),
                           "gyroxstdev" = extractsd(data_temp, "wx"),
                           "gyroxmax" = extractmax(data_temp, "wx"),
                           "gyroxmin" = extractmin(data_temp, "wx"),
                           "gyroxmaxamp" = gx[1],
                           "gyroxmaxfreq" = gx[6],
                           "tgyroxmaxamp" = gx[2],
                           "tgyroxmaxfreq" = gx[3],
                           "gyroxenergy" = gx[4],
                           "gyroxentropy" = gx[5],
                           "gyroyavg" = extractavg(data_temp, "wy"),
                           "gyroystdev" = extractsd(data_temp, "wy"),
                           "gyroymax" = extractmax(data_temp, "wy"),
                           "gyroymin" = extractmin(data_temp, "wy"),
                           "gyroymaxamp" = gy[1],
                           "gyroymaxfreq" = gy[6],
                           "tgyroymaxamp" = gy[2],
                           "tgyroymaxfreq" = gy[3],
                           "gyroyenergy" = gy[4],
                           "gyroyentropy" = gy[5],
                           "gyrozavg" = extractavg(data_temp, "wz"),
                           "gyrozstdev" = extractsd(data_temp, "wz"),
                           "gyrozmax" = extractmax(data_temp, "wz"),
                           "gyrozmin" = extractmin(data_temp, "wz"),
                           "gyrozmaxamp" = gz[1],
                           "gyrozmaxfreq" = gz[6],
                           "tgyrozmaxamp" = gz[2],
                           "tgyrozmaxfreq" = gz[3],
                           "gyrozenergy" = gz[4],
                           "gyrozentropy" = gy[5],
                           "sampleno" = nrow(data_temp))
  
  #attach the tuple to the output frame
  output_data <- rbind(output_data, data_tuple)
  
}

output_data <- output_data %>%
  mutate(Action="")

if (readline(prompt = "Vehicle filter? (Y/N) : ") == "Y") {
  output_data <- filter(output_data,speedavg>1)
}

output_data$Action <- readline(prompt = "Enter Label : ")

existing_data <- read.csv("trainingset.csv")

existing_data <- rbind(existing_data, output_data)


write.csv(existing_data, file="trainingset.csv", row.names = FALSE)
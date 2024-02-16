library(dplyr)
## library(tidyverse)
library(tidyr)
library(rhdf5)
library(stringi)
library(lubridate)
library(readr)
library(doParallel)
library(foreach)



## Function Usage
## directory: Any directory stores waveform datasets with an extension ".h5". It doesn't have to be the root directory. 
## full-time: Boolean, whether or not you would like to get a full duration of each waveform/trend for each participant, default is 'TRUE'.
## time_window: Numeric (unit = hrs), the duration of interest, use only when "full_time" is set to "FALSE", default is NULL. 
## time_break: Numeric (unit = mins), the time interval used for calculcating moving averages.
## clinic_file: The full path of the clinical dataset.

get_n_min_avg <- function(directory, full_time = T, time_window = NULL, time_break, clinic_file) {
     comb <- function(d1, d2) {
          if (nrow(d1) != nrow(d2)) {as.data.frame(t(bind_rows(as.data.frame(t(d1)), as.data.frame(t.data.frame(d2)))))} 
          else {cbind(d1, d2)}
          }    ## create a function to inform foreach loop how to combine result from each iteration at level 2
     subject_files <- list.files(directory, recursive = T, pattern = '*.h5')    ## search all files with an extension ".h5" under the directory, including sub-folders
     full_path <- paste0(directory, subject_files)
     clinic <- read_csv(clinic_file)    ## load clinical data
     if (full_time == T) {    ## Below is the code for extracting complete non-eeg data for each patient
          if (!missing(time_window)) {
               stop('argument "full_time" must be set to "FALSE" when "time_window" has a value')
               } else {
                    registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))    ## create do parallel cluster based on # of tasks per node
                    aggr_data <- foreach(a = full_path, .combine = plyr::rbind.fill,     ## level 1 foreach loop to search over the HDF5 file(s) for each patient
                                    .packages = c('stringi', 'plyr', 'dplyr', 'rhdf5', 'tidyr', 'readr', 'lubridate', 'doParallel')) %:% 
                         foreach(b = iter(h5ls(a, datasetinfo = F) %>%    ## level 2 foreach loop to return the field names of non-eeg forms/trends and of their timestamps under each patient, then iterating each pair of names by row.
                                          filter(group %in% c('/numerics', '/waves') & !grepl('Timestamps|Quality', name, ignore.case = T)) %>% 
                                          select(1:2) %>% 
                                          unite(field_name, group, name, sep = '/', remove = F) %>% 
                                          mutate(field_time = paste0(field_name, ' ', 'Timestamps')) %>% 
                                          select(3, 1, 4), by = 'row'), 
                                 i = icount(),    ## get the real-time iteration # 
                                 .combine = comb) %dopar% {
                                      subject_num <- stri_extract_first_regex(a, '(?<=_)[\\d]+(?=_)')    ## get subject # from file name
                                      site_num <- stri_extract_first_regex(a, '(?<=/)[\\d]+(?=_)')    ## get site # from file name
                                      field <- tryCatch(h5read(a, b[, 2], bit64conversion='bit64'), error = function(err) NULL)    ## read each non-eeg waveform/trend
                                      time <- as.POSIXct(tryCatch(h5read(a, b[, 3], bit64conversion='bit64'), error = function(err) NULL) / 1e6, origin = '1970-01-01', tz = 'GMT')    ## read each timestamp
                                      full_field <- aggregate(field, list(cut(time[1:length(field)], breaks = paste0(time_break, ' min'))), mean)[, 'x']    ## calculate moving averages for each xx-min interval
                                      full_window <- aggregate(time[1:length(field)], list(cut(time[1:length(field)], breaks = paste0(time_break, ' min'))), first)[, 'x']    ## get start time for each xx-min interval
                                      means_time <- data.frame(full_field, full_window)
                                      colnames(means_time)[1:2] <- c(b[, 1], paste0(b[, 1], '_start'))
                                      EMS <- clinic[clinic$zSubjectID == subject_num, c('F502Q02', 'F502Q03', 'F502Q06', 'F502Q07')][1,]    ## get EMS info from clinical dataset for each patient
                                      call_time <- as.POSIXct(paste(EMS$F502Q02, EMS$F502Q03), format = "%d%b%Y %H:%M:%S",
                                                              tz = case_when(site_num %in% c('1753') ~ 'America/Chicago', 
                                                                             site_num %in% c('1125', '1435') ~ 'America/Los_Angeles', 
                                                                             T ~ 'America/New_York'))
                                      arrival_time <- as.POSIXct(paste(EMS$F502Q06, EMS$F502Q07), format = "%d%b%Y %H:%M:%S",
                                                                 tz = case_when(site_num %in% c('1753') ~ 'America/Chicago', 
                                                                                site_num %in% c('1125', '1435') ~ 'America/Los_Angeles', 
                                                                                T ~ 'America/New_York'))
                                      if (i == 1) {    ## inlcude EMS info into the first iteration, instead of every iteration
                                           n_hr_subject <- bind_cols(subject.id = subject_num, 
                                                                     site = site_num, 
                                                                     ems_call = call_time, 
                                                                     ems_arrival = arrival_time, 
                                                                     means_time)
                                           print(subject_num)
                                           } else {
                                                n_hr_subject <- means_time
                                                }
                                      n_hr_subject
                                      }
                    getDoParWorkers()
                    gc()
                    h5closeAll()
                    }
          } else {    ## Below is the code for extracting non-eeg data with any specified duration, which is similar to the section above
               if (missing(time_window)) {
                    stop('argument "time_window" must have a value when "full_time" is set to "FALSE"')
                    } else {
                         registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
                         aggr_data <- foreach(a = full_path, .combine = plyr::rbind.fill, 
                                              .packages = c('stringi', 'plyr', 'dplyr', 'rhdf5', 'tidyr', 'readr', 'lubridate', 'doParallel')) %:% 
                              foreach(b = iter(h5ls(a, datasetinfo = F) %>% 
                                                    filter(group %in% c('/numerics', '/waves') & !grepl('Timestamps|Quality', name, ignore.case = T)) %>% 
                                                    select(1:2) %>% 
                                                    unite(field_name, group, name, sep = '/', remove = F) %>% 
                                                    mutate(field_time = paste0(field_name, ' ', 'Timestamps')) %>% 
                                                    select(3, 1, 4), by = 'row'),
                                      i = icount(), 
                                      .combine = cbind) %dopar% {
                                           subject_num <- stri_extract_first_regex(a, '(?<=_)[\\d]+(?=_)')
                                           site_num <- stri_extract_first_regex(a, '(?<=/)[\\d]+(?=_)')
                                           field <- tryCatch(h5read(a, b[, 2], bit64conversion='bit64'), error = function(err) NULL)
                                           time <- as.POSIXct(tryCatch(h5read(a, b[, 3], bit64conversion='bit64'), error = function(err) NULL) / 1e6, origin = '1970-01-01', tz = 'GMT')
                                           n_hr_window <- time[time < min(time) + hours(time_window)]    ## subset of timestamps from the beginning based on the time window specified by user
                                           means_time <- data.frame(aggregate(field[1:length(n_hr_window)],    ## calculate moving averages for each xx-min interval
                                                                              list(cut(n_hr_window, breaks = paste0(time_break, ' min'))),
                                                                              mean)[1:(time_window*60/time_break), 'x'],
                                                                    aggregate(n_hr_window,    ## get start time for each xx-min interval
                                                                              list(cut(n_hr_window, breaks = paste0(time_break, ' min'))),
                                                                              first)[1:(time_window*60/time_break), 'x'])
                                           colnames(means_time)[1:2] <- c(b[, 1], paste0(b[, 1], '_start'))
                                           EMS <- clinic[clinic$zSubjectID == subject_num, c('F502Q02', 'F502Q03', 'F502Q06', 'F502Q07')][1, ]
                                           call_time <- as.POSIXct(paste(EMS$F502Q02, EMS$F502Q03), format = "%d%b%Y %H:%M:%S",
                                                                   tz = case_when(site_num %in% c('1753') ~ 'America/Chicago', 
                                                                                  site_num %in% c('1125', '1435') ~ 'America/Los_Angeles', 
                                                                                  T ~ 'America/New_York'))
                                           arrival_time <- as.POSIXct(paste(EMS$F502Q06, EMS$F502Q07), format = "%d%b%Y %H:%M:%S",
                                                                      tz = case_when(site_num %in% c('1753') ~ 'America/Chicago', 
                                                                                     site_num %in% c('1125', '1435') ~ 'America/Los_Angeles', 
                                                                                     T ~ 'America/New_York'))
                                           if (i == 1) {
                                                n_hr_subject <- bind_cols(subject.id = subject_num, 
                                                                          site = site_num, 
                                                                          ems_call = call_time, 
                                                                          ems_arrival = arrival_time, 
                                                                          means_time)
                                                print(subject_num)
                                           } else {
                                                n_hr_subject <- means_time
                                           }
                                           n_hr_subject
                                      }
                         getDoParWorkers()
                         gc()
                         h5closeAll()
                    }
               }
     aggr_data <- aggr_data %>%
          fill(1:4, .direction = 'down')
     aggr_data
     }



## first_hr_window <- get_n_min_avg('/oak/stanford/groups/zihuai/precicecap_data_main/', 1, 5, '/oak/stanford/groups/zihuai/precicecap_data_main/Clinic Data/Output/baseline_CRF.csv')

## write_csv(first_hr_window, '/oak/stanford/groups/zihuai/precicecap_data_main/Output Data/1hr_5min_avg.csv')

## first_72hr_window <- get_n_min_avg('/oak/stanford/groups/zihuai/precicecap_data_main/', full_time = F, time_window = 72, time_break = 5, '/oak/stanford/groups/zihuai/precicecap_data_main/Clinic Data/Output/baseline_CRF.csv')

## write_csv(first_72hr_window, '/oak/stanford/groups/zihuai/precicecap_data_main/Output Data/72hr_5min_avg.csv')


first_full_window <- get_n_min_avg('/oak/stanford/groups/zihuai/precicecap_data_main/', full_time = T, time_break = 5, clinic_file = '/oak/stanford/groups/zihuai/precicecap_data_main/Clinic Data/Output/baseline_CRF.csv')

write_csv(first_full_window, '/oak/stanford/groups/zihuai/precicecap_data_main/Output Data/full_5min_avg.csv')

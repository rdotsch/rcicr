
# Load reverse correlation toolbox
library(rcicr)
library(stringr)

# R Studio에서 .r 파일이 있는 경로 찾기
rstudioapi::getActiveDocumentContext
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


ci_each_participant <- 0

dir <- "analysis/"
outputdir <- paste0(dir, "output/")
filenames <- list.files(dir, pattern = "\\.csv$")
filenum <- length(filenames)

female_rdata <- "stimuli_female/rcic_seed_1_time_8_31_2024_18_30.Rdata"
male_rdata <- "stimuli_male/rcic_seed_1_time_Aug_31_2024_21_32.Rdata"

data <- NULL;
ci_list <- list(); ci_anti_list <- list()

for (this in 1:filenum){
  
  # load data
  filename <- paste0(dir, filenames[this])
  tem_data <- read.csv(filename)
  
  # organize data 1: Fill in all rows of survey response
  tem_demograph <- tem_data[, c("gender", "age", "nationality", "ses.response")]
  tem_demograph <- subset(tem_demograph, ses.response != "")
  tem_data$gender <- tem_demograph$gender
  tem_data$age <- tem_demograph$age
  tem_data$nationality <- tem_demograph$nationality
  tem_data$ses.response <- tem_demograph$ses.response
  
  # organize data 2: Remove empty cells
  tem_data <- subset(tem_data, selectedstim != "")


  # generate CIs -------------------------------
  if (tem_data$imgpath[1] == "stimuli_male/"){
    tem_rdata <- male_rdata
  } 
  if (tem_data$imgpath[1] == "stimuli_female/"){
    tem_rdata <- female_rdata
  } 
  
  # Extract stimulus number based on trial number
  tem_data$stim <- tem_data$imgset
  
  # Extract ori/inv from selectedstim
  tem_data$oriinv <- str_match(tem_data$selectedstim, "_([invori]+).")[,2]
  
  # Recode left/right selection to weights in CI
  tem_data$response[tem_data$oriinv == "ori"] <- 1
  tem_data$response[tem_data$oriinv == "inv"] <- -1
  
  if (ci_each_participant == 1){
    # CI for each participant ----------------------------
    tem_id <- tem_data$id[1]
    tem_name <- paste0("p", tem_id)
    
    tem_ci <- generateCI2IFC(tem_data$stim, tem_data$response, "base", tem_rdata, scaling="matched",
                             filename = tem_name, targetpath = outputdir)
    tem_ci_list <- setNames(list(tem_ci), tem_name)
    ci_list <- c(ci_list, tem_ci_list)
    
    tem_ci_anti <- generateCI2IFC(tem_data$stim, tem_data$response, "base", tem_rdata, scaling="matched", 
                                  antiCI = TRUE, filename = tem_name, targetpath = outputdir)
    tem_ci_anti_list <- setNames(list(tem_ci_anti), tem_name)
    ci_anti_list <- c(ci_anti_list, tem_ci_anti_list)  
  }

  data <- rbind(data, tem_data)
  tem_data = NULL;
}
if (ci_each_participant == 1){
  scaled_cis <- autoscale(ci_list, save_as_pngs = TRUE, targetpath = outputdir)
  scaled_cis_anti <- autoscale(ci_anti_list, save_as_pngs = TRUE, targetpath = outputdir)
}
# CI for all participants ----------------------------

data_f <- subset(data, imgpath == "stimuli_female/")
data_m <- subset(data, imgpath == "stimuli_male/")

filename_f <- "all_F"; filename_m <- "all_M"

if (nrow(data_f) >= 1){
  all_ci_anti_f <- generateCI2IFC(data_f$stim, data_f$response, "base", female_rdata, scaling="matched", 
                                  antiCI = TRUE, filename = filename_f, targetpath = outputdir)
  all_ci_f <- generateCI2IFC(data_f$stim, data_f$response, "base", female_rdata, scaling="matched",
                             filename = filename_f, targetpath = outputdir)
}

if (nrow(data_m) >= 1){
  all_ci_m <- generateCI2IFC(data_m$stim, data_m$response, "base", male_rdata, scaling="matched",
                             filename = filename_m, targetpath = outputdir)
  all_ci_anti_m <- generateCI2IFC(data_m$stim, data_m$response, "base", male_rdata, scaling="matched", 
                                antiCI = TRUE, filename = filename_m, targetpath = outputdir)
}

write.table(data, paste0(dir, "output/data_all.csv"), sep = ",", row.names = FALSE, col.names = TRUE, fileEncoding = "cp949")
write.table(data_f, paste0(dir, "output/data_f.csv"), sep = ",", row.names = FALSE, col.names = TRUE, fileEncoding = "cp949")
write.table(data_m, paste0(dir, "output/data_m.csv"), sep = ",", row.names = FALSE, col.names = TRUE, fileEncoding = "cp949")
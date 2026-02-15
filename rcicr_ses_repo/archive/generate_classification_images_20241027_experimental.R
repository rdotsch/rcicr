
library(dplyr)  # bind_rows 함수 포함된 패키지

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

## female_rdata <- "stimuli_female/rcic_seed_1_time_8_31_2024_18_30.Rdata"
female_rdata_anti <- "rcic_seed_1_time_8_31_2024_18_30_ANTI.Rdata"
female_rdata_anti_mod <- "rcic_seed_1_time_8_31_2024_18_30_MODIFIED.Rdata"
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
    tem_rdata <- female_rdata_anti_mod
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
  
  data <- bind_rows(data, tem_data)   # 열 개수 맞지 않을 때 나타나는 error 해결할 수 있음
  
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
  all_ci_anti_f <- generateCI2IFC(data_f$stim, data_f$response, "base", female_rdata_anti_mod, scaling="matched", 
                                  antiCI = TRUE, filename = filename_f, targetpath = outputdir)
  ##all_ci_f <- generateCI2IFC(data_f$stim, data_f$response, "base", female_rdata, scaling="matched",
    ##                         filename = filename_f, targetpath = outputdir)
  
  all_ci_f_zmap <- generateCI(data_f$stim, data_f$response, "base", female_rdata_anti_mod,
                              filename = paste0(filename_f, "_zmap"),
                              targetpath = dir, zmap = TRUE, zmapmethod = "t.test",sigma = 1, threshold = 1 )
  
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




## 테ㅐ스트
# 🔹 1. 기존 Rdata 파일 로드

load("stimuli_female/rcic_seed_1_time_8_31_2024_18_30.Rdata")
load("rcic_seed_1_time_8_31_2024_18_30_ANTI.Rdata")

# 🔹 2. 기존 Base 이미지 확인

print("📌 기존 Base 얼굴 데이터 확인:")
print(dim(base_faces))  # base_faces 구조 확인

# 🔹 3. Base 얼굴을 Anti-CI 이미지로 교체

##밑에 내용 확인
base_faces <- all_ci_anti_f  # Anti-CI를 새로운 베이스 이미지로 설정

# 🔹 4. 기존 데이터의 모든 변수 저장

save(list = ls(), file = "rcic_seed_1_time_8_31_2024_18_30_ANTI.Rdata")





# 🔹 2. base_faces에서 특정 요소 삭제
if ("ci" %in% names(base_faces)) base_faces$ci <- NULL
if ("scaled" %in% names(base_faces)) base_faces$scaled <- NULL
if ("combined" %in% names(base_faces)) base_faces$combined <- NULL

# 🔹 3. 변경된 base_faces 확인
print(names(base_faces))  # ✅ "base"만 남아 있는지 확인

# 🔹 4. 변경된 Rdata 파일 저장
save(list = ls(), file = "rcic_seed_1_time_8_31_2024_18_30_MODIFIED.Rdata")



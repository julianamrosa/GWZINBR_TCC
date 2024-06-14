library(readr)
#korea_base_artigo <- read_csv("UnB/2024/TCC2/korea_base_artigo.csv")
korea_base_artigo <- read_csv("C:/Users/Juliana Rosa/OneDrive/Documents/TCC2/GWZINBR_TCC-main/korea_base_artigo.csv")
#korea_base_artigo <- read_csv("C:/Juliana/TCC/GWZINBR_TCC-main/korea_base_artigo.csv")

if(!require(sp)){
  install.packages("sp")
  library(sp)
}

## GOLDEN ##

## adaptive

#zinb, adaptive, aic
startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = c("Healthcare_access", "Crowding"), lat = "y", long = "x", offset = "ln_total",
       model = "zinb", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#36 mins

#zip, adaptive, aic
startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = c("Healthcare_access", "Crowding"), lat = "y", long = "x", offset = "ln_total",
       model = "zip", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#10 mins - diferente do sas

#zinb, adaptive, cv
startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = c("Healthcare_access", "Crowding"), lat = "y", long = "x", offset = "ln_total",
       model = "zinb", method = "adaptive_bsq", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#16 mins

#zip, adaptive, cv
startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = c("Healthcare_access", "Crowding"), lat = "y", long = "x", offset = "ln_total",
       model = "zip", method = "adaptive_bsq", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#1.4 min

## fixed

#zinb, fixed, aic
startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = c("Healthcare_access", "Crowding"), lat = "y", long = "x", offset = "ln_total",
       model = "zinb", method = "fixed_g", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#12 mins

#zip, fixed, aic
startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = c("Healthcare_access", "Crowding"), lat = "y", long = "x", offset = "ln_total",
       model = "zip", method = "fixed_g", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#10 mins - diferente do sas

#zinb, fixed, cv
startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = c("Healthcare_access", "Crowding"), lat = "y", long = "x", offset = "ln_total",
       model = "zinb", method = "fixed_g", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#14 mins

#zip, fixed, cv
startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = c("Healthcare_access", "Crowding"), lat = "y", long = "x", offset = "ln_total",
       model = "zip", method = "fixed_g", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#1.1 min

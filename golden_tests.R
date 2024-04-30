library(readr)
#korea_base_artigo <- read_csv("UnB/2024/TCC2/korea_base_artigo.csv")
#korea_base_artigo <- read_csv("C:/Users/Juliana Rosa/OneDrive/Documents/TCC2/GWZINBR-main/korea_base_artigo.csv")
korea_base_artigo <- read_csv("C:/Juliana/TCC/GWZINBR-main/korea_base_artigo.csv")

if(!require(sp)){
  install.packages("sp")
  library(sp)
}

### TESTES ###

## ZINB ##

# Teste 1: adaptive_bsq - cv #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = "Healthcare_access", weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zinb", method = "adaptive_bsq", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#7 mins

# Teste 2: adaptive_bsq - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = "Healthcare_access", weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zinb", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#14 mins

# Teste 3: fixed_g - cv #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = "Healthcare_access", weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zinb", method = "fixed_g", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#20 mins

# Teste 4: fixed_g - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = "Healthcare_access", weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zinb", method = "fixed_g", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#34 mins

## ZIP ##

# Teste 9: adaptive_bsq - cv #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = "Healthcare_access", weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zip", method = "adaptive_bsq", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#2.2 mins

# Teste 10: adaptive_bsq - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = "Healthcare_access", weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zip", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#9 mins --> cvs um pouco diferentes

# Teste 11: fixed_g - cv #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = "Healthcare_access", weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zip", method = "fixed_g", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#1 min

# Teste 12: fixed_g - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zip", method = "fixed_g", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#43 mins --> errado

## NEGBIN ##

# Teste 17: adaptive_bsq - cv #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "negbin", method = "adaptive_bsq", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime

# Teste 18: adaptive_bsq - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "negbin", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#57 segs

# Teste 19: fixed_g - cv #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "negbin", method = "fixed_g", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#1.2 mins

# Teste 20: fixed_g - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "negbin", method = "fixed_g", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#45 segs

## POISSON ##

# Teste 25: adaptive_bsq - cv #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "poisson", method = "adaptive_bsq", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime

# Teste 26: adaptive_bsq - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "poisson", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime

# Teste 27: fixed_g - cv #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "poisson", method = "fixed_g", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#17 segs

# Teste 28: fixed_g - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "poisson", method = "fixed_g", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime

#lembretes:
#3- rodar todos os testes
#4- return --> invisible

#Obs.: implementar gráfico;

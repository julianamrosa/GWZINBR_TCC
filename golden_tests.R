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
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zinb", method = "adaptive_bsq", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE)
endTime <- Sys.time()
endTime-startTime
#1.5 mins

# Teste 2: adaptive_bsq - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zinb", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE)
endTime <- Sys.time()
endTime-startTime
#2.2 mins

# Teste 3: fixed_g - cv #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zinb", method = "fixed_g", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE)
endTime <- Sys.time()
endTime-startTime
#erro

# Teste 4: fixed_g - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "zinb", method = "fixed_g", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE)
endTime <- Sys.time()
endTime-startTime
#erro

## NEGBIN ##

# Teste 18: adaptive_bsq - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "negbin", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE)
endTime <- Sys.time()
endTime-startTime
#56 segs

# Teste 19: fixed_g - cv #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "negbin", method = "fixed_g", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE)
endTime <- Sys.time()
endTime-startTime
#erro

# Teste 20: fixed_g - aic #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "negbin", method = "fixed_g", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE)
endTime <- Sys.time()
endTime-startTime
#erro

## POISSON ##

# Teste 27: fixed_g - cv #

startTime <- Sys.time()
Golden(data = korea_base_artigo,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = NULL, weight = NULL, lat = "x", long = "y", offset = "ln_total",
       model = "poisson", method = "fixed_g", bandwidth = "cv", globalmin = FALSE, distancekm = TRUE)
endTime <- Sys.time()
endTime-startTime
#erro: depois de 22 mins

#dúvidas: band deve conter os valores iniciais ou não?
#         verificar nomes dos objetos na saída (e elementos num geral)
#         decimais em CV são importantes?
#         anotar o que é cada parâmetro

#lembretes:
#3- rodar todos os testes
#4- return --> invisible

#possibilidades: tirar min_bandwidth da saída, manter só gss_results (trocar o nome para min_bandwidth e trocar o label xmin para bandwidth)
#                tirar h? mas preciso de acesso fácil ao seu valor...
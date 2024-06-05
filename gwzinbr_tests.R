
library(readr)

#paths Jess
korea_base_artigo <- read_csv("C:/Users/jehhv/OneDrive/Documentos/UnB/2024/TCC2/korea_base_artigo.csv")
korea_base_artigo <- read_csv("D:/Users/jessica.abreu/Documents/UnB/tcc/korea_base_artigo.csv")

#paths Ju
korea_base_artigo <- read_csv("C:/Users/Juliana Rosa/OneDrive/Documents/TCC2/GWZINBR_TCC-main/korea_base_artigo.csv")
korea_base_artigo <- read_csv("C:/Juliana/TCC/GWZINBR_TCC-main/korea_base_artigo.csv")

## TESTE ARTIGO ##

startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
                diff_sd+Crowding+Migration+Health_behavior,
        xvarinf = c("Healthcare_access", "Crowding"),
        lat = "y", long = "x", offset = "ln_total", method = "adaptive_bsq",
        model = "zinb", distancekm = TRUE, h=82, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#1.4 min

## FIXED ##

#Teste 1 do relatório
startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        xvarinf = "Healthcare_access",
        lat = "x", long = "y", offset = "ln_total", method = "fixed_g",
        model = "zinb", distancekm = TRUE, h=226.73, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#1.6 mins

#Teste 2 do relatório
startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        xvarinf = "Healthcare_access",
        lat = "x", long = "y", offset = "ln_total", method = "fixed_g",
        model = "zip", distancekm = TRUE, h=733.70, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#6 segs

#Teste 3 do relatório
startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        lat = "x", long = "y", offset = "ln_total", method = "fixed_g",
        model = "negbin", distancekm = TRUE, h=189.74, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#5 segs

#Teste 4 do relatório
startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        lat = "x", long = "y", offset = "ln_total", method = "fixed_g",
        model = "poisson", distancekm = TRUE, h=733.70, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#5 segs

## ADAPTIVE ##

#Teste 5
startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        xvarinf = "Healthcare_access",
        lat = "x", long = "y", offset = "ln_total", method = "adaptive_bsq",
        model = "zinb", distancekm = TRUE, h=230, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#1.1 min

#Teste 6
startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        xvarinf = "Healthcare_access",
        lat = "x", long = "y", offset = "ln_total", method = "adaptive_bsq",
        model = "zip", distancekm = TRUE, h=230, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#7 segs

#Teste 7
startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        lat = "x", long = "y", offset = "ln_total", method = "adaptive_bsq",
        model = "negbin", distancekm = TRUE, h=230, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#2 segs

#Teste 8
startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        lat = "x", long = "y", offset = "ln_total", method = "adaptive_bsq",
        model = "poisson", distancekm = TRUE, h=48, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#1 seg

#decidir teste do relatório
#mandar dúvidas pro Alan
#rodar teste escolhido, comparar e tirar prints
#subir para o CRAN
#relatório!!!

#antes de subir para o cran, (1) tirar views e transformar em outputs e (2) trocar return por invisible
#estamos deixando assim por enquanto para rodar os testes e tirar os prints

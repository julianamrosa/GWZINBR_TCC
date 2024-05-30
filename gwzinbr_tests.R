
library(readr)

#paths Jess
korea_base_artigo <- read_csv("C:/Users/jehhv/OneDrive/Documentos/UnB/2024/TCC2/korea_base_artigo.csv")
korea_base_artigo <- read_csv("D:/Users/jessica.abreu/Documents/UnB/tcc/korea_base_artigo.csv")

#paths Ju
korea_base_artigo <- read_csv("C:/Users/Juliana Rosa/OneDrive/Documents/TCC2/GWZINBR-main/korea_base_artigo.csv")
korea_base_artigo <- read_csv("C:/Juliana/TCC/GWZINBR-main/korea_base_artigo.csv")

#Teste 1 do relatório
startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        xvarinf = "Healthcare_access",
        lat = "x", long = "y", offset = "ln_total", method = "fixed_g",
        model = "zinb", distancekm = TRUE, h=226.73)
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
        model = "zip", distancekm = TRUE, h=733.70)
endTime <- Sys.time()
endTime-startTime
#6 segs

#Teste 3 do relatório
startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        lat = "x", long = "y", offset = "ln_total", method = "fixed_g",
        model = "negbin", distancekm = TRUE, h=189.74)
endTime <- Sys.time()
endTime-startTime
#5 segs

#Teste 4 do relatório
startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        lat = "x", long = "y", offset = "ln_total", method = "fixed_g",
        model = "poisson", distancekm = TRUE, h=733.70)
endTime <- Sys.time()
endTime-startTime
#5 segs

#fazer alterações do pibic
#decidir teste do relatório
#rodar e comparar certinho
#subir para o CRAN
#relatório!!!

#antes de subir para o cran, (1) tirar views e transformar em outputs e (2) trocar return por invisible
#estamos deixando assim por enquanto para rodar os testes e tirar os prints

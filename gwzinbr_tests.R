
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

# alterar para vetor: dbb dlb e todos os outros
# linha 744  varabetalambda: trocar matrizes por rep(), se necessario
#tirar NULL dos defaults

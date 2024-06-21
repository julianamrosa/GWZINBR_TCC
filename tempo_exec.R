library(readr)
library(gwzinbr)
korea_df <- read_csv("C:/Users/Juliana Rosa/OneDrive/Documents/TCC2/GWZINBR_TCC-main/korea_base_artigo.csv")

#POISSON

startTime <- Sys.time()
Golden(data = korea_df,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, lat = "y", long = "x", offset = "ln_total",
       model = "poisson", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime

startTime <- Sys.time()
gwzinbr(data = korea_df, formula = n_covid1~Morbidity+high_sch_p
        +Healthcare_access+diff_sd+Crowding+Migration+Health_behavior,
        lat = "y", long = "x", offset = "ln_total",
        method = "adaptive_bsq", model = "poisson", distancekm = TRUE,
        h=79, force=TRUE)
endTime <- Sys.time()
endTime-startTime

#NEGBIN

startTime <- Sys.time()
Golden(data = korea_df,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, lat = "y", long = "x", offset = "ln_total",
       model = "negbin", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime

startTime <- Sys.time()
gwzinbr(data = korea_df, formula = n_covid1~Morbidity+high_sch_p
        +Healthcare_access+diff_sd+Crowding+Migration+Health_behavior,
        lat = "y", long = "x", offset = "ln_total",
        method = "adaptive_bsq", model = "negbin", distancekm = TRUE,
        h=82, force=TRUE)
endTime <- Sys.time()
endTime-startTime

#PIZ

startTime <- Sys.time()
Golden(data = korea_df,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = c("Healthcare_access", "Crowding"), lat = "y", long = "x", offset = "ln_total",
       model = "zip", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime

startTime <- Sys.time()
gwzinbr(data = korea_df, formula = n_covid1~Morbidity+high_sch_p
        +Healthcare_access+diff_sd+Crowding+Migration+Health_behavior,
        xvarinf = c("Healthcare_access", "Crowding"), lat = "y",
        long = "x", offset = "ln_total", method = "adaptive_bsq",
        model = "zip", distancekm = TRUE, h=56, force=TRUE)
endTime <- Sys.time()
endTime-startTime
#3.8 mins

#BNIZ

startTime <- Sys.time()
Golden(data = korea_df,formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
         diff_sd+Crowding+Migration+Health_behavior, xvarinf = c("Healthcare_access", "Crowding"), lat = "y", long = "x", offset = "ln_total",
       model = "zinb", method = "adaptive_bsq", bandwidth = "aic", globalmin = FALSE, distancekm = TRUE, force=TRUE)
endTime <- Sys.time()
endTime-startTime

startTime <- Sys.time()
gwzinbr(data = korea_df, formula = n_covid1~Morbidity+high_sch_p
        +Healthcare_access+diff_sd+Crowding+Migration+Health_behavior,
        xvarinf = c("Healthcare_access", "Crowding"), lat = "y",
        long = "x", offset = "ln_total", method = "adaptive_bsq",
        model = "zinb", distancekm = TRUE, h=82, force=TRUE)
endTime <- Sys.time()
endTime-startTime
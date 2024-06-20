##### CARREGANDO PACOTES E DADOS #####

library(leaflet)
library(readr)
library(dplyr)
korea_df <- read_csv("C:/Users/Juliana Rosa/OneDrive/Documents/TCC2/GWZINBR_TCC-main/korea_base_artigo.csv")
#korea_df <- read_csv("C:/Users/jehhv/OneDrive/Documentos/UnB/2024/TCC2/GWZINBR_TCC-main/GWZINBR_TCC-main/korea_base_artigo.csv")


##### 1- MAPA VAR RESPOSTA #####

#Criando escala de cores

korea_df <- korea_df %>% mutate(n_covid1_col=case_when(
  n_covid1 == 0 ~ "#F0F0F0",     
  n_covid1 == 1 ~ "#D7EAF3",     
  n_covid1 == 2 ~ "#B0D5E6",     
  n_covid1 %in% 3:4 ~ "#89C0D9", 
  n_covid1 == 5 ~ "#62ABCC",     
  n_covid1 %in% 6:8 ~ "#3A96BF", 
  n_covid1 %in% 9:19 ~ "#1381B2",
  n_covid1 >= 20 ~ "#006C9E"     
))

#Criando popup do mapa

state_popup <- paste0("<strong>Local: </strong>", 
                      korea_df$IDNAME, 
                      "<br><strong>Casos: </strong>", 
                      korea_df$n_covid1)

#Mapa

mapa_resposta <- korea_df %>%
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%  
  addCircleMarkers(
    color = ~n_covid1_col,
    fillColor = ~n_covid1_col,
    fillOpacity = 0.8,
    weight = 1,
    radius = 7,
    stroke = TRUE,
    popup = state_popup,
    lat = ~y,
    lng = ~x
  ) %>%
  addLegend(
    labels = c("0", "1", "2", "3-4", "5", "6-8", "9-19", ">19"),
    colors = c("#F0F0F0", "#D7EAF3", "#B0D5E6", "#89C0D9", "#62ABCC", "#3A96BF",
               "#1381B2", "#006C9E"),
    title = "Casos de COVID-19"
  )

mapa_resposta


##### 2- MAPA VAR EXPLICATIVA - Crowding #####

#Definindo classes

k <- floor(1+3.3*log10(nrow(korea_df))) #número de classes
minimum <- min(korea_df$Crowding)
maximum <- max(korea_df$Crowding)
h <- round((maximum - minimum)/k, 1) #amplitude das classes
min_trunc <- trunc(minimum*10)/10
seq(min_trunc, min_trunc+8*h, h) #classes

#Criando escala de cores

korea_df <- korea_df %>% mutate(crowding_col=case_when(
  Crowding >= 9.9 & Crowding < 11.5 ~ "#F0F0F0",     
  Crowding >= 11.5 & Crowding < 13.1 ~ "#D7EAF3",     
  Crowding >= 13.1 & Crowding < 14.7 ~ "#B0D5E6",     
  Crowding >= 14.7 & Crowding < 16.3 ~ "#89C0D9", 
  Crowding >= 16.3 & Crowding < 17.9 ~ "#62ABCC",     
  Crowding >= 17.9 & Crowding < 19.5 ~ "#3A96BF", 
  Crowding >= 19.5 & Crowding < 21.1 ~ "#1381B2",
  Crowding >= 21.1 & Crowding < 22.7 ~ "#006C9E"     
))

#Criando popup do mapa

state_popup <- paste0("<strong>Local: </strong>", 
                      korea_df$IDNAME, 
                      "<br><strong>Nível de Aglomeração: </strong>", 
                      round(korea_df$Crowding, 3))

#Mapa

#sugestão: usar razão de chances para trazer infos mais interessantes 
mapa_explicativa <- korea_df %>%
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%  
  addCircleMarkers(
    color = ~crowding_col,
    fillColor = ~crowding_col,
    fillOpacity = 0.8,
    weight = 1,
    radius = 7,
    stroke = TRUE,
    popup = state_popup,
    lat = ~y,
    lng = ~x
  ) %>%
  addLegend(
    labels = c("9.9 |- 11.5", "11.5 |- 13.1", "13.1 |- 14.7", "14.7 |- 16.3",
               "16.3 |- 17.9", "17.9 |- 19.5", "19.5 |- 21.1", "21.1 |- 22.7"),
    colors = c("#F0F0F0", "#D7EAF3", "#B0D5E6", "#89C0D9", "#62ABCC", "#3A96BF",
               "#1381B2", "#006C9E"),
    title = "Nível de Aglomeração"
  )

mapa_explicativa

##### 3- MAPA BETA - Crowding #####

#Pegando as estimativas

#obs.: rodar definição da função gwzinbr antes
mod <- gwzinbr(data = korea_df, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        xvarinf = c("Healthcare_access", "Crowding"),
        lat = "y", long = "x", offset = "ln_total", method = "adaptive_bsq",
        model = "zinb", distancekm = TRUE, h=82, force=TRUE)

korea_df$Crowding_est <- round(as.numeric((mod$parameter_estimates)[,"Crowding"]), 3)

#Definindo classes

minimum <- min(korea_df$Crowding_est)
maximum <- max(korea_df$Crowding_est)
h <- round((maximum - minimum)/k, 3) #amplitude das classes
min_trunc <- trunc(minimum*1000)/1000
seq(min_trunc, min_trunc+k*h, h) #classes

#Criando escala de cores

korea_df <- korea_df %>% mutate(crowding_est_col=case_when(
  Crowding_est >= (-0.744) & Crowding_est < (-0.577) ~ "#F0F0F0",     
  Crowding_est >= (-0.577) & Crowding_est < (-0.410) ~ "#D7EAF3",     
  Crowding_est >= (-0.410) & Crowding_est < (-0.243) ~ "#B0D5E6",     
  Crowding_est >= (-0.243) & Crowding_est < (-0.076) ~ "#89C0D9", 
  Crowding_est >= (-0.076) & Crowding_est < 0.091 ~ "#62ABCC",     
  Crowding_est >= 0.091 & Crowding_est < 0.258 ~ "#3A96BF", 
  Crowding_est >= 0.258 & Crowding_est < 0.425 ~ "#1381B2",
  Crowding_est >= 0.425 & Crowding_est < 0.592 ~ "#006C9E"     
))

#Criando popup do mapa

state_popup <- paste0("<strong>Local: </strong>", 
                      korea_df$IDNAME, 
                      "<br><strong>Estimativa para Crowding: </strong>", 
                      korea_df$Crowding_est)

#Mapa

#mudar escala de cores para incluir o aspecto negativo->positivo
mapa_beta <- korea_df %>%
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%  
  addCircleMarkers(
    color = ~crowding_est_col,
    fillColor = ~crowding_est_col,
    fillOpacity = 0.8,
    weight = 1,
    radius = 7,
    stroke = TRUE,
    popup = state_popup,
    lat = ~y,
    lng = ~x
  ) %>%
  addLegend(
    labels = c("-0.744 |- -0.577", "-0.577 |- -0.410", "-0.410 |- -0.243",
               "-0.243 |- -0.076", "-0.076 |- 0.091", "0.091 |- 0.258",
               "0.258 |- 0.425", "0.425 |- 0.592"),
    colors = c("#F0F0F0", "#D7EAF3", "#B0D5E6", "#89C0D9", "#62ABCC",
               "#3A96BF", "#1381B2", "#006C9E"),
    title = "Estimativa para Parâmetro Crowding"
  )

mapa_beta

##### 3.1- MAPA BETA - Crowding #####

#Pegando as estimativas

#obs.: rodar definição da função gwzinbr antes
mod <- gwzinbr(data = korea_df, 
               formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
                 diff_sd+Crowding+Migration+Health_behavior,
               xvarinf = c("Healthcare_access", "Crowding"),
               lat = "y", long = "x", offset = "ln_total", method = "adaptive_bsq",
               model = "zinb", distancekm = TRUE, h=82, force=TRUE)

korea_df$Crowding_est <- round(as.numeric((mod$parameter_estimates)[,"Crowding"]), 3)

#Definindo classes

minimum <- min(korea_df$Crowding_est)
maximum <- max(korea_df$Crowding_est)
minimum
maximum

#Criando escala de cores

korea_df <- korea_df %>% mutate(crowding_est_col=case_when(
  Crowding_est >= (-0.8) & Crowding_est < (-0.4) ~ "#ef746f",     
  Crowding_est >= (-0.4) & Crowding_est < 0 ~ "#fbb4b9",     
  Crowding_est ==0 ~ "#ffffff",     
  Crowding_est > 0 & Crowding_est <= 0.4 ~ "#9ecae1", 
  Crowding_est > 0.4 & Crowding_est <= 0.8 ~ "#4292c6"
))
#scale_color <- c("#ef746f", "#fdd0a2", "#ffffff", "#9ecae1", "#4292c6")

#Criando popup do mapa

state_popup <- paste0("<strong>Local: </strong>", 
                      korea_df$IDNAME, 
                      "<br><strong>Estimativa para Crowding: </strong>", 
                      korea_df$Crowding_est)

#Mapa

#mudar escala de cores para incluir o aspecto negativo->positivo
mapa_beta <- korea_df %>%
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%  
  addCircleMarkers(
    color = ~crowding_est_col,
    fillColor = ~crowding_est_col,
    fillOpacity = 0.8,
    weight = 1,
    radius = 7,
    stroke = TRUE,
    popup = state_popup,
    lat = ~y,
    lng = ~x
  ) %>%
  addLegend(
    labels = c("-0.8 |- 0.4", "-0.4 |- 0", "0", "0 -| 0.4", "0.4 -| 0.8"),
    colors = c("#ef746f", "#fbb4b9", "#ffffff", "#9ecae1", "#4292c6"),
    title = "Estimativa para Parâmetro Crowding"
  )

mapa_beta

##### 4- MAPA ALPHA #####

#Pegando as estimativas

korea_df$alpha <- round(as.numeric((mod$alpha_estimates)[,"alpha"]), 3)

#Definindo classes

minimum <- min(korea_df$alpha)
maximum <- max(korea_df$alpha)
h <- ceiling(1000*(maximum - minimum)/k)/1000 #amplitude das classes
min_trunc <- trunc(minimum*1000)/1000
seq(min_trunc, min_trunc+k*h, h) #classes

#Criando escala de cores

korea_df <- korea_df %>% mutate(alpha_col=case_when(
  alpha >= 0.000 & alpha < 0.509 ~ "#F0F0F0",     
  alpha >= 0.509 & alpha < 1.018 ~ "#D7EAF3",     
  alpha >= 1.018 & alpha < 1.527 ~ "#B0D5E6",     
  alpha >= 1.527 & alpha < 2.036 ~ "#89C0D9", 
  alpha >= 2.036 & alpha < 2.545 ~ "#62ABCC",     
  alpha >= 2.545 & alpha < 3.054 ~ "#3A96BF", 
  alpha >= 3.054 & alpha < 3.563 ~ "#1381B2",
  alpha >= 3.563 & alpha < 4.072 ~ "#006C9E"     
))

#Criando popup do mapa

state_popup <- paste0("<strong>Local: </strong>", 
                      korea_df$IDNAME, 
                      "<br><strong>Estimativa para alpha: </strong>", 
                      korea_df$alpha)

#Mapa

#mudar escala de cores para incluir o aspecto negativo->positivo
mapa_alpha <- korea_df %>%
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%  
  addCircleMarkers(
    color = ~alpha_col,
    fillColor = ~alpha_col,
    fillOpacity = 0.8,
    weight = 1,
    radius = 7,
    stroke = TRUE,
    popup = state_popup,
    lat = ~y,
    lng = ~x
  ) %>%
  addLegend(
    labels = c("0.000 |- 0.509", "0.509 |- 1.018", "1.018 |- 1.527",
               "1.527 |- 2.036", "2.036 |- 2.545", "2.545 |- 3.054",
               "3.054 |- 3.563", "3.563 |- 4.072"),
    colors = c("#F0F0F0", "#D7EAF3", "#B0D5E6", "#89C0D9", "#62ABCC",
               "#3A96BF", "#1381B2", "#006C9E"),
    title = "Estimativa para Parâmetro alpha"
  )

mapa_alpha

##### 5- MAPA LAMBDA #####

#Pegando as estimativas

korea_df$lambda <- round(as.numeric((mod$parameter_estimates)[,"Inf_Crowding"]), 1)

#Definindo classes

minimum <- min(korea_df$lambda)
maximum <- max(korea_df$lambda)
h <- ceiling(10*(maximum - minimum)/k)/10 #amplitude das classes
min_trunc <- trunc(minimum*10)/10
seq(min_trunc, min_trunc+k*h, h) #classes

#Criando escala de cores

korea_df <- korea_df %>% mutate(lambda_col=case_when(
  lambda >= (-296.5) & lambda < (-232.6) ~ "#F0F0F0",     
  lambda >= (-232.6) & lambda < (-168.7) ~ "#D7EAF3",     
  lambda >= (-168.7) & lambda < (-104.8) ~ "#B0D5E6",     
  lambda >= (-104.8) & lambda < (-40.9) ~ "#89C0D9", 
  lambda >= (-40.9) & lambda < 23.0 ~ "#62ABCC",     
  lambda >= 23.0 & lambda < 86.9 ~ "#3A96BF", 
  lambda >= 86.9 & lambda < 150.8 ~ "#1381B2",
  lambda >= 150.8 & lambda < 214.7 ~ "#006C9E"     
))

#Criando popup do mapa

state_popup <- paste0("<strong>Local: </strong>", 
                      korea_df$IDNAME, 
                      "<br><strong>Estimativa para Crowding (IZ): </strong>", 
                      korea_df$lambda)

#Mapa

#mudar escala de cores para incluir o aspecto negativo->positivo
mapa_lambda <- korea_df %>%
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%  
  addCircleMarkers(
    color = ~lambda_col,
    fillColor = ~lambda_col,
    fillOpacity = 0.8,
    weight = 1,
    radius = 7,
    stroke = TRUE,
    popup = state_popup,
    lat = ~y,
    lng = ~x
  ) %>%
  addLegend(
    labels = c("-296.5 |- -232.6",     
                 "-232.6 |- -168.7",     
                 "-168.7 |- -104.8",     
                 "-104.8 |- -40.9", 
                 "-40.9 |- 23.0",     
                 "23.0 |- 86.9", 
                 "86.9 |- 150.8",
                 "150.8 |- 214.7"),
    colors = c("#F0F0F0", "#D7EAF3", "#B0D5E6", "#89C0D9", "#62ABCC",
               "#3A96BF", "#1381B2", "#006C9E"),
    title = "Estimativa para Parâmetro Crowding IZ"
  )

mapa_lambda

##### 5.1- MAPA LAMBDA #####

#Pegando as estimativas

korea_df$lambda <- round(as.numeric((mod$parameter_estimates)[,"Inf_Crowding"]), 1)

#Definindo classes

minimum <- min(korea_df$lambda)
maximum <- max(korea_df$lambda)
minimum
maximum

#Criando escala de cores

korea_df <- korea_df %>% mutate(lambda_col=case_when(
  lambda >= (-300) & lambda < (-150) ~ "#ef746f",     
  lambda >= (-150) & lambda < 0 ~ "#fbb4b9",     
  lambda == 0 ~ "#ffffff",     
  lambda > 0 & lambda <= 150 ~ "#9ecae1", 
  lambda > 150 & lambda <= 300 ~ "#4292c6"
))
#scale_color <- c("#ef746f", "#fdd0a2", "#ffffff", "#9ecae1", "#4292c6")

#Criando popup do mapa

state_popup <- paste0("<strong>Local: </strong>", 
                      korea_df$IDNAME, 
                      "<br><strong>Estimativa para Crowding (IZ): </strong>", 
                      korea_df$lambda)

#Mapa

#mudar escala de cores para incluir o aspecto negativo->positivo
mapa_lambda <- korea_df %>%
  leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%  
  addCircleMarkers(
    color = ~lambda_col,
    fillColor = ~lambda_col,
    fillOpacity = 0.8,
    weight = 1,
    radius = 7,
    stroke = TRUE,
    popup = state_popup,
    lat = ~y,
    lng = ~x
  ) %>%
  addLegend(
    labels = c("-300 |- -150", "-150 |- 0", "0", "0 -| 150", "150 -| 300"),
    colors = c("#ef746f", "#fbb4b9", "#ffffff", "#9ecae1", "#4292c6"),
    title = "Estimativa para Parâmetro Crowding IZ"
  )

mapa_lambda

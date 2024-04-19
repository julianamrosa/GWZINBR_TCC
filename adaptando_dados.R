library(ggplot2)
library(dplyr)

#ZIP

ggplot(korea_base_artigo)+
  geom_histogram(aes(n_covid1), bins=100)

korea_base_zip <- korea_base_artigo%>%
  filter(n_covid1<500)

ggplot(korea_base_zip)+
  geom_histogram(aes(n_covid1), bins=100)

write.csv(korea_base_zip, "C:/Users/Juliana Rosa/OneDrive/Documents/TCC2/GWZINBR-main/korea_base_zip.csv", row.names=FALSE)

#NEGBIN

ggplot(korea_base_artigo)+
  geom_histogram(aes(n_covid1), bins=100)

korea_base_negbin <- korea_base_artigo[-sample(which(korea_base_artigo$n_covid1==0), 51), ]

ggplot(korea_base_negbin)+
  geom_bar(aes(n_covid1), bins=100)

#POISSON

ggplot(korea_base_artigo)+
  geom_histogram(aes(n_covid1), bins=100)

korea_base_poisson <- korea_base_artigo[-sample(which(korea_base_artigo$n_covid1==0), 51), ]%>%filter(n_covid1<10)

ggplot(korea_base_poisson)+
  geom_histogram(aes(n_covid1), bins=100)

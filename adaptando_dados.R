library(ggplot2)
library(dplyr)

korea_base_zip <- korea_base_artigo%>%
  filter(n_covid1<500)

ggplot(korea_base_zip)+
  geom_boxplot(aes(n_covid1))

ggplot(korea_base_zip)+
  geom_histogram(aes(n_covid1), bins=100)

ggplot(korea_base_artigo)+
  geom_histogram(aes(n_covid1), bins=100)

write.csv(korea_base_zip, "C:/Juliana/TCC/GWZINBR-main/korea_base_zip.csv", row.names=FALSE)



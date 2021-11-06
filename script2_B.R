# IBE 875 - Modelagem de Distribuição de Espécies
# PPGE/PPGBio
# Professores: Rodrigo Tardin, Maria Lucia Lorini, Mariana Vasconcellos
# Script 2 - Processamento e exploração visual de preditores ambientais

####################################
# Script 2 - Preditores ambientais #
####################################

#### Definindo diretório de trabalho ####
setwd("D:/MDE") #Mude para o endereço da pasta da disciplina no seu computador
getwd()

#### Carregando pacotes ####
library(biomod2)
library(raster)
library(rgdal) #Carregar shapefiles
library(sdmpredictors) #Carregar as camadas online
library(maptools)
library(usdm) #Teste de Inflação de Variância (VIF) para 'stacks'
library(ecospat)
library(CoordinateCleaner)
library(rgbif)
library(spocc)
library(spThin)
library(dplyr)
library(dismo)
library(gridExtra)
library(rgdal)
library(rgeos)
library(dplyr)

#### Listando e carregando os datasets pre-existentes ####

list_datasets()

#Listando as camadas do WorldClim
list_layers("WorldClim")

#Obtendo camadas do Worldclim a partir da pasta da disciplina (se vc fez o download completo da pasta do GoogleDrive o caminho abaixo deve funcionar)
#Camadas topográficas
alt=raster("./Camadas/Presente/WC_alt_lonlat.tif")
alt
#Altitude

#Camadas climáticas
bio1=raster("./Camadas/Presente/WC_bio1_lonlat.tif")
#Temperatura Anual Média

bio3=raster("./Camadas/Presente/WC_bio3_lonlat.tif")
#Isotermalidade

bio4=raster("./Camadas/Presente/WC_bio4_lonlat.tif")
#Sazonalidade da temperatura

bio7=raster("./Camadas/Presente/WC_bio7_lonlat.tif")
#Variação Anual da Temperatura

bio12=raster("./Camadas/Presente/WC_bio12_lonlat.tif")
#Precipitação Anual

bio15=raster("./Camadas/Presente/WC_bio15_lonlat.tif")
#Sazonalidade da precipitação

#### OPCIONAL ####
#Obtendo camadas do WorldClim online em caso de não possuir/querer baixar

#alt=load_layers("WC_alt", datadir = "./Camadas")
#Altitude

#bio1=load_layers("WC_bio1", datadir = "./Camadas")
#Temperatura Média Anual

#bio3=load_layers("WC_bio3", datadir = "./Camadas")
#Isotermalidade

#bio4=load_layers("WC_bio4", datadir = "./Camadas")
#Sazonalidade da temperatura (SD x)

#bio7=load_layers("WC_bio7", datadir = "./Camadas")
#Variação Anual de Temperatura 

#bio12=load_layers("WC_bio12", datadir = "./Camadas")
#Precipitação Anual

#bio15=load_layers("WC_bio15", datadir = "./Camadas")
#Sazonalidade da precipitação (CV)

#### Explorando os detalhes das camadas climáticas ####

#Obtendo os limites e detalhes das camadas
bbox(alt); ncol(alt); nrow(alt) ; res(alt)
bbox(bio1); ncol(bio1); nrow(bio1) ; res(bio1)
bbox(bio3); ncol(bio3); nrow(bio3) ; res(bio3)
bbox(bio4); ncol(bio4); nrow(bio4) ; res(bio4)
bbox(bio7); ncol(bio7); nrow(bio7) ; res(bio7)
bbox(bio12); ncol(bio12); nrow(bio12) ; res(bio12)
bbox(bio15); ncol(bio15); nrow(bio15) ; res(bio15)
compareRaster(alt, bio1, bio3, bio4, bio7, bio12, bio15)

#### Visualizando as camadas climáticas ####
plot(alt, main="Altitude")
par(mfrow=c(3,2),mar=c(3,3,2,0)) #ajustar parâmetros gráficos para particionar a área do plot (mfrow) em 3 linhas e 2 colunas e ajustar as margens (mar)
plot(bio1, col=topo.colors(255), main="Temperatura Média Anual (bio1)")
plot(bio3, col=topo.colors(255), main="Isotermalidade (bio3)")
plot(bio4, col=topo.colors(255), main="Sazonalidade da temperatura (bio4)")
plot(bio7, col=topo.colors(255), main="Variação Anual de Temperatura (bio7)")
plot(bio12, col=topo.colors(255), main="Precipitação Anual (bio12)")
plot(bio15, col=topo.colors(255), main="Sazonalidade da precipitação (bio15)")
dev.off() #retorna os parâmetros gráficos para o default
###

#### Processando as camadas climáticas ####

#Cortando as camadas para o Eurafrica
#Carregando o mapa do Eurafrica
# AQUI EU PEGUEI O RECORTE DE PORTUGAL PARA O TREINO

Portugal=getData('GADM', country= 'Portugal', level=1, path = "./Camadas")
Portugal <- gUnaryUnion(Portugal)
plot(Portugal)
Eurafrica=readOGR("./Eurafrica/eurafrica.shp")
plot(Eurafrica)
#Cortando cada camada de uma vez para o Eurafrica e 'perfumando' minimamente os mapas

altEurafrica=mask(crop(alt,Eurafrica),Eurafrica)
plot(altEurafrica, main = "Altitude Eurafrica",xlab = "Longitude", ylab = "Latitude")

bio1Eurafrica=mask(crop(bio1,Eurafrica),Eurafrica)
plot(bio1Eurafrica, main = "bio1 Eurafrica",xlab = "Longitude", ylab = "Latitude")

bio3Eurafrica=mask(crop(bio3,Eurafrica),Eurafrica)
plot(bio3Eurafrica, main = "bio3 Eurafrica",xlab = "Longitude", ylab = "Latitude")

bio4Eurafrica=mask(crop(bio4,Eurafrica),Eurafrica)
plot(bio4Eurafrica, main = "bio4 Eurafrica",xlab = "Longitude", ylab = "Latitude")

bio7Eurafrica=mask(crop(bio7,Eurafrica),Eurafrica)
plot(bio7Eurafrica, main = "bio7 Eurafrica",xlab = "Longitude", ylab = "Latitude")

bio12Eurafrica=mask(crop(bio12,Eurafrica),Eurafrica)
plot(bio12Eurafrica, main = "bio12 Eurafrica",xlab = "Longitude", ylab = "Latitude")

bio15Eurafrica=mask(crop(bio15,Eurafrica),Eurafrica)
plot(bio15Eurafrica, main = "bio15 Eurafrica",xlab = "Longitude", ylab = "Latitude")


#Obtendo os limites e detalhes das camadas cortadas para o Eurafrica
bbox(altEurafrica); ncol(altEurafrica); nrow(altEurafrica) ; res(altEurafrica)
bbox(bio1Eurafrica); ncol(bio1Eurafrica); nrow(bio1Eurafrica) ; res(bio1Eurafrica)
bbox(bio3Eurafrica); ncol(bio3Eurafrica); nrow(bio3Eurafrica) ; res(bio3Eurafrica)
bbox(bio4Eurafrica); ncol(bio4Eurafrica); nrow(bio4Eurafrica) ; res(bio4Eurafrica)
bbox(bio7Eurafrica); ncol(bio7Eurafrica); nrow(bio7Eurafrica) ; res(bio7Eurafrica)
bbox(bio12Eurafrica); ncol(bio12Eurafrica); nrow(bio12Eurafrica) ; res(bio12Eurafrica)
bbox(bio15Eurafrica); ncol(bio15Eurafrica); nrow(bio15Eurafrica) ; res(bio15Eurafrica)
compareRaster(altEurafrica, bio1Eurafrica, bio3Eurafrica, bio4Eurafrica, bio7Eurafrica, bio12Eurafrica, bio15Eurafrica)

#### Comparando as dimensões das camadas Eurafrica e Mundo
bbox(bio1); ncol(bio1); nrow(bio1) ; res(bio1)
bbox(bio1Eurafrica); ncol(bio1Eurafrica); nrow(bio1Eurafrica) ; res(bio1Eurafrica)

#### Unindo as variáveis em um unico 'stack' ####

biostack=stack(altEurafrica,bio1Eurafrica,bio3Eurafrica,bio4Eurafrica,bio7Eurafrica,bio12Eurafrica,bio15Eurafrica)
plot(biostack)

#### Cortando as variáveis para o polígono da Araponga ####

lusitanicapolygon=Portugal
plot(lusitanicapolygon)

#Plotando o polígono da distribuição da Araponga no mapa do Eurafrica
plot(altEurafrica,  main="Altitude Eurafrica")
plot(lusitanicapolygon, add=T)

#Cortando cada camada de uma vez para o poligono da Araponga e 'perfumando' minimamente os mapas

altlusitanica=mask(crop(alt,lusitanicapolygon),lusitanicapolygon)
plot(altlusitanica, main = "Altitude lusitanica",xlab = "Longitude", ylab = "Latitude")
plot(lusitanicapolygon, add=T)

bio1lusitanica=mask(crop(bio1,lusitanicapolygon),lusitanicapolygon)
plot(bio1lusitanica, main = "bio1 lusitanica",xlab = "Longitude", ylab = "Latitude")

bio3lusitanica=mask(crop(bio3,lusitanicapolygon),lusitanicapolygon)
plot(bio3lusitanica, main = "bio3 lusitanica",xlab = "Longitude", ylab = "Latitude")

bio4lusitanica=mask(crop(bio4,lusitanicapolygon),lusitanicapolygon)
plot(bio4lusitanica, main = "bio4 lusitanica",xlab = "Longitude", ylab = "Latitude")

bio7lusitanica=mask(crop(bio7,lusitanicapolygon),lusitanicapolygon)
plot(bio7lusitanica, main = "bio7 lusitanica",xlab = "Longitude", ylab = "Latitude")

bio12lusitanica=mask(crop(bio12,lusitanicapolygon),lusitanicapolygon)
plot(bio12lusitanica, main = "bio12 lusitanica",xlab = "Longitude", ylab = "Latitude")

bio15lusitanica=mask(crop(bio15,lusitanicapolygon),lusitanicapolygon)
plot(bio15lusitanica, main = "bio15 lusitanica",xlab = "Longitude", ylab = "Latitude")

#Cortando o 'stack' com as variáveis para o polígono da Araponga
biostack1=mask(crop(biostack,lusitanicapolygon),lusitanicapolygon)
biostack1
summary(biostack1)

#Plotando o 'stack' com as camadas climáticas dentro do poligono da  Araponga
plot(biostack1)

#Checando existência de correlação e multicolinearidade entre as variáveis

#Correlação entre as camadas usando a função 'pairs'
pairs(biostack1)

#Visualizar a correlação entre as camadas usando o pacote corrplot
set.seed(1963) #seed number para sempre gerar os mesmos pontos abaixo
backgr <- randomPoints(biostack1, 10000)
absclim <- data.frame(raster::extract(biostack1, backgr)) # o comando :: forca o uso da funcao especifica daquele pacote, evitando conflito entre pacotes -> raster::extract = utilize a funcao 'extract' do pacote 'raster' 
absclim.std <- data.frame(scale(absclim)) # Normalizar as variáveis

library(corrplot)
M <- cor(absclim.std)
corrplot.mixed(M, upper = "ellipse", lower = "number",number.cex = 0.8,tl.cex = 0.8)

# Checagem de multicolinearidade entre as camadas - remoção stepwise com threshold = 0.7

vifcor(biostack1,th=0.7)
#th = valor de corte para correlação

vifstep(biostack1, th = 10)
#th = valor de corte para VIF


#Unindo as variáveis em um unico 'stack' após as análises de correlação e multicolinearidade

biostack2=stack(altlusitanica,bio3lusitanica, bio7lusitanica,bio12lusitanica,bio15lusitanica)
biostack2

#Plotando o novo 'stack' com as camadas climáticas não correlacionadas e não-colinares dentro do poligono da  Araponga
plot(biostack2, main=c("Altitude", "Isotermalidade (bio3)", "Variação Anual de Temperatura (bio7)", "Precipitação Anual (bio12)","Sazonalidade da precipitação (bio15)"))

#Salvando o espaço de trabalho com todos os objetos num documento RData que pode ser carregado posteriormente.
save.image(file="script2_B.RData")

## Fim do script2
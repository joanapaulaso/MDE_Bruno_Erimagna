# IBE 875 - Modelagem de Distribuição de Espécies
# PPGE/PPGBio
# Professores: Rodrigo Tardin, Maria Lucia Lorini, Mariana Vasconcellos
# Script 6 - Projeção, transferabilidade e incerteza

######################################################
# Script 6 - Projeção, transferabilidade e incerteza #
######################################################

# Geração dos modelos usando os algoritmos de classificação (GBM, RF) e aprendizado de máquina (MaxEnt)

#### Definindo area de trabalho ####
setwd("D:/MDE") #Mude para o endereco da pasta da disciplina no seu computador
getwd()

#### Carregando pacotes ####
library(biomod2)
library(raster)

#### Carregando os objetos do script 2, 3 e 4 ####
load("script2_B.RData")
load("script3_B.RData")
load("script4_B.RData")
load("script5_B.RData")

#### Projetando as predições dos modelos no espaço geográfico - presente ####

#Gerando o objeto com as predições espaciais de todos os modelos e rodadas separadamente para GLM, GAM, SRE
projec_procnias1 <- BIOMOD_Projection(
  modeling.output = procnias1model, #objeto com os modelos criados
  new.env = biostack2, #variaveis ambientais onde os modelos vão ser projetados
  proj.name = "GLM_GAM_SRE", #nome da pasta com as projeções
  selected.models = "all",
  binary.meth = "TSS", #metodo para binarização
    output.format = ".img", compress = ' xz ' ,
  clamping.mask = F, do.stack = T
)

# Sumário do objeto com as projeções no espaço geográfico
projec_procnias1


#Gerando o objeto com as predições espaciais de todos os modelos e rodadas separadamente para GBM, RF, MaxEnt
projec_procnias2 <- BIOMOD_Projection(
  modeling.output = procnias2model, #objeto com os modelos criados
  new.env = biostack2, #variaveis ambientais onde os modelos vão ser projetados
  proj.name = "GBM_RF_MAX", #nome da pasta com as projeções
  selected.models = "all",
  binary.meth = "TSS", #metodo para binarização
  output.format = ".img", compress = ' xz ' ,
  clamping.mask = F, do.stack = T
)


# Sumário do objeto com as projeções no espaço geográfico
projec_procnias2

# Projetando as predições dos modelos individuais por algoritmo
plot(projec_procnias1, str.grep = "GLM")
plot(projec_procnias1, str.grep = "GAM")
plot(projec_procnias1, str.grep = "SRE")
plot(projec_procnias2, str.grep = "GBM")
plot(projec_procnias2, str.grep = "RF")
plot(projec_procnias2, str.grep = "MAXENT.Phillips")

#### Preparando as camadas para transferir as predições do modelo no tempo ####

#Camadas obtidas do WorldClim (CMIP 6) dentro dos Shared Socio-Economic Pathway 585 para os anos de 2060-2080 (https://www.carbonbrief.org/cmip6-the-next-generation-of-climate-models-explained).

#Carregando as camadas futuras correspondentes as que foram usadas para modelar no tempo presente  (bio 3, bio 7, bio 12, bio 15)
#Nomes dos arquivos das camadas na pasta tem que ser os mesmos do presente para 'dar match' e ser possível projetar no futuro

bio3_2060=raster("./Camadas/Futuro/WC_bio3_lonlat.tif") 
#Isotermalidade 2060

bio7_2060=raster("./Camadas/Futuro/WC_bio7_lonlat.tif")
#Variação Anual da Temperatura 2060

bio12_2060=raster("./Camadas/Futuro/WC_bio12_lonlat.tif")
#Precipitação Anual 2060

bio15_2060=raster("./Camadas/Futuro/WC_bio15_lonlat.tif")
#Sazonalidade da precipitação 2060

#### Cortando as camadas futuras para o polígono da Araponga ####
bio3procnias_futuro=mask(crop(bio3_2060,procniaspolygon),procniaspolygon)

bio7procnias_futuro=mask(crop(bio7_2060,procniaspolygon),procniaspolygon)

bio12procnias_futuro=mask(crop(bio12_2060,procniaspolygon),procniaspolygon)

bio15procnias_futuro=mask(crop(bio15_2060,procniaspolygon),procniaspolygon)

#Formatando a bio3_2060 para ficar no mesmo 'range' de valores do que o bio 3 no tempo presente
bio3procnias_futuro=bio3procnias_futuro/100

#### Unindo as variáveis futuras em um unico 'stack' ####
biostack_futuro=stack(altprocnias,bio3procnias_futuro,bio7procnias_futuro,bio12procnias_futuro,bio15procnias_futuro)

#### Plotando o 'stack' das variáveis futuras e comparando com o 'stack' do presente (melhorar com os códigos da Mari) ####
plot(biostack_futuro, col=topo.colors(255))
plot(biostack2, col=topo.colors(255))

#### Projetando no espaço geográfico os modelos do presente nas condições climáticas futuras ####

projec_procnias1_futuro <- BIOMOD_Projection(
  modeling.output = procnias1model, #modelo do presente_envelope e regressao
  new.env = biostack_futuro, #aqui você indica as camadas futuras que vão ser usadas para projetar
  proj.name = "Modelos Futuro_GLM_GAM_SRE", #nome da pasta
  selected.models = "all",
  binary.meth = "TSS",
  output.format = ".img", compress = ' xz ' ,
  clamping.mask = T, do.stack = T
)


projec_procnias2_futuro <- BIOMOD_Projection(
  modeling.output = procnias2model, #modelo do presente_classificação e maxima entropia
  new.env = biostack_futuro, #aqui você indica as camadas futuras que vão ser usadas para projetar
  proj.name = "Modelos Futuro_RF_BRT_MAXENT", #nome da pasta
  selected.models = "all",
  binary.meth = "TSS",
  output.format = ".img", compress = ' xz ' ,
  clamping.mask = T, do.stack = T
)

#### Plotando os modelos individuais projetados no futuro ####
plot(projec_procnias1_futuro, str.grep = "GLM")
plot(projec_procnias1_futuro, str.grep = "GAM")
plot(projec_procnias1_futuro, str.grep = "SRE")
plot(projec_procnias2_futuro, str.grep = "GBM")
plot(projec_procnias2_futuro, str.grep = "RF")
plot(projec_procnias2_futuro, str.grep = "MAXENT.Phillips")

#Salvando o espaço de trabalho com todos os objetos num documento RData que pode ser carregado posteriormente.
save.image(file="script6.RData")

# Fim do Script 6
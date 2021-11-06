# IBE 875 - Modelagem de Distribuição de Especies
# PPGE/PPGBio
# Professores: Rodrigo Tardin, Maria Lucia Lorini, Mariana Vasconcellos
# Script 3 - Registros de ocorrencia, pseudo-ausencias e background.

########################################################
# Script 3 - Ocorrencia, pseudo-ausencias e background #
########################################################

#### Definindo area de trabalho ####
setwd("D:/MDE") #Mude para o endereco da pasta da disciplina no seu computador
getwd()

#### Carregando pacotes ####
library(tidyr)
library(maptools)
library(spThin)
library(raster)
library(biomod2)

#### Carregando os objetos do script 2 ####
load("script2_B.RData")

#### Carregando os registros de ocorrencia ####

#### Obtendo os registros do diretorio ####
lusitanic=read.csv("lusitanic.csv", header=T,  sep = ";")
View(lusitanic)

#### Checando inconsistencias relacionadas aos registros de ocorrencias - problemas, filtro e limpeza ####
# Caso existam linhas nas colunas das coordenadas sem informacoes (NA), a funcao as remove
lusitanic_na <- lusitanic %>% tidyr::drop_na(decimalLongitude, decimalLatitude)
View(lusitanic_na)
#Quantos registros foram removidos por nao conter coordenadas?

# Limpando registros com outros problemas referente as coordenadas 
# Marcando os registros potencialmente problematicos 
flags_spatial <- CoordinateCleaner::clean_coordinates(
  x = lusitanic_na, 
  species = "species",
  lon = "decimalLongitude", 
  lat = "decimalLatitude",
  tests = c("capitals", # raio ao redor de capitais
            "centroids", # raio ao redor de centroides de paises e provincias
            "duplicates", # duplicatas
            "equal", # coordenadas iguais
            "gbif", # raio ao redor da sede da GBIF
            "institutions", # raio ao redor de instituicoes de pesquisa em biodiversidade
            "seas", # pontos no mar
            "urban", # pontos dentro de areas urbanas
            "validity", # ponto de fora do sistema de coordenadas
            "zeros" # zeros e pontos onde lat = lon 
  )
)

# Resultado dos registros que foram problematicos
# TRUE = coordenadas 'limpas'
# FALSE = coordenadas potencialmente problematicas
head(flags_spatial)
summary(flags_spatial)

# Excluindo os pontos considerados como problematicos na analise anterior
lusitanic_f <- lusitanic_na %>% 
  dplyr::filter(flags_spatial$.summary == TRUE)
#Quantos registros ficaram retidos após a remoção dos problematicos?

#### Limpando e selecionando as colunas de interesse ####
lusitanic_f = lusitanic_f %>%
  dplyr::select(species, decimalLongitude, decimalLatitude)

#Checando a extensão do arquivo limpo
nrow(lusitanic_f)
colnames(lusitanic_f)

#### Plotando os registros de ocorrencia filtrados para visualizacao ####
plot(altlusitanica, xlim=c(-20,30),ylim=c(10,60), main="Pontos de ocorrencia por atitude")
data(wrld_simpl)
plot(wrld_simpl, add=TRUE, border='darkgray', lwd=1)
points(lusitanic_f$decimalLongitude, lusitanic_f$decimalLatitude, col='blue', cex=0.3)

#### Filtrando os registros de ocorrencia baseado em uma distancia 'x' para evitar potenciais problemas de autocorrelacao espacial ####
dir.create("./thinned") # criando um diretorio para salvar os dados a serem filtrados

lusitanic_thin <- 
  thin( loc.data = lusitanic_f, 
        lat.col = "decimalLatitude", long.col = "decimalLongitude", 
        spec.col = "species", 
        thin.par = 50 #distancia em km que os registros vao ficar separados#  
        ,reps = 5, 
        locs.thinned.list.return = TRUE, 
        write.files = TRUE, 
        max.files = 2, 
        out.dir = "./thinned", out.base = "lusitanic_thinned", 
        write.log.file = TRUE,
        log.file = "./thinned/lusitanic_thinned_log.txt")

#Explorando visualmente o efeito do processo de filtragem espacial dos dados
plotThin(lusitanic_thin) # Clique 'enter' no console varias vezes para visualizar diferentes graficos

#Valores chegaram no valor máximo com as repeticoes usadas. De uma olhada no arquivo 'log' com todos os detalhes

#Carregando o novo arquivo com o processo de filtragem espacial realizado.
lusitanic_thin=read.csv("./thinned/lusitanic_thinned_thin1.csv", sep=",")
View(lusitanic_thin)

#### Plotando os registros de ocorrencia filtrados em relacao aos antigos ####

plot(alt, xlim=c(-20,30),ylim=c(20,50), main="Pontos de ocorrencia por atitude")
plot(wrld_simpl, add=TRUE, border='darkgray', lwd=1)
points(lusitanic_f$decimalLongitude, lusitanic_f$decimalLatitude, col='blue', cex=0.3)
points(lusitanic_thin$decimalLongitude, lusitanic_thin$decimalLatitude, col='red', cex=0.2)

#Inserindo uma coluna contendo todos os registros de ocorrencia igual a 1
resp.occ=as.numeric("resp.occ")
str(resp.occ)
lusitanic_thin[,"resp.occ"] = 1
resp.occ=lusitanic_thin$resp.occ
is.numeric(resp.occ)

#Checando se houve adicao da coluna
colnames(lusitanic_thin)
View(lusitanic_thin)

# Agora que os registros foram filtrados quanto a inconsistencias e potenciais problemas de autocorrelacao espacial, estamos prontos para gerar as pseudo-ausencias

#### Formatando os registros de ocorrencia para a modelagem no ambiente BIOMOD ####

#Salvando um objeto apenas com o nome da especie
resp.name <- 'lusitanic'

# dados de presença da nossa especie 
is.numeric(lusitanic_thin$resp.occ)
myResp <- as.numeric(lusitanic_thin[,"resp.occ"])

#Criando um objeto apenas com as coordenadas da Araponga e concatenando o valor 1 para todas as coordenadas
myRespXY <- lusitanic_thin[which(myResp==1),c("decimalLongitude", "decimalLatitude")]

# Criando uma camada raster baseada nas camadas ambientais do script 2 para nossa resposta variavel 
myResp <- reclassify(subset(biostack2,1,drop=TRUE), c(-Inf,Inf,0))
myResp[cellFromXY(myResp,myRespXY)] <- 1

#Formatando os dados no modelo BIOMOD para depois iniciar o procedimento de modelagem com diferentes algoritmos
#Dado a falta de ausencia no conjunto de dados, temos que informar a criacao das pseudo-ausencias
#Estrategia aleatoria
#100 pseudo ausencias e 4 sets
bm.lusitanic100 = BIOMOD_FormatingData(resp.var = resp.occ,
                                      expl.var = biostack2,
                                      resp.xy = myRespXY,
                                      resp.name = resp.name,
                                      PA.nb.rep = 4,
                                      PA.nb.absences = 100, 
                                      PA.strategy = "random")

#sumario do novo objeto criado
bm.lusitanic100


#Visualizacao da configuracao espacial das pseudo-ausencias
plot(bm.lusitanic100) 


#Agora vamos mudar a estrategia de criacao das pseudo-ausencias para "disk", onde estabelecemos a distancia minima e maxima entre as pseudo-ausencias e os registros de ocorrencia 


#100 pseudo ausencias com estrategia 'disk'
bm.lusitanic100disk = BIOMOD_FormatingData(resp.var = resp.occ,
                                          expl.var = biostack2,
                                          resp.xy = myRespXY,
                                          resp.name = resp.name,
                                          PA.nb.rep = 4,
                                          PA.nb.absences = 100, 
                                          PA.strategy = "disk", 
                                          PA.dist.min = 100000) #distancia em metros

#sumario do novo objeto criado
bm.lusitanic100disk

#Visualização da configuração espacial das pseudo-ausencias
plot(bm.lusitanic100disk) 

#Agora vamos pre-estabelecer que a nossa area onde serao geradas as pseudo-ausencias sera fora do envelope climatico gerado com o algoritmo GAM. 
# JOANA: Ajustei esse algoritimo para GAM. Para esse conjunto de dados o SRE não pôde gerar pseudo-ausências
# Vejam o que é melhor para vocês quando introduzirem os shp corretos.
bm.lusitanic100sre = BIOMOD_FormatingData(resp.var = resp.occ,
                                         expl.var = biostack2,
                                         resp.xy = myRespXY,
                                         resp.name = resp.name,
                                         PA.nb.rep = 5,
                                         PA.nb.absences = 500, 
                                         PA.strategy = "sre", 
                                         PA.sre.quant = 0.025, 
                                         na.rm = T)

#sumario do novo objeto criado
bm.lusitanic100sre                                          
#Visualização da configuração espacial das Pseudo Ausências
plot(bm.lusitanic100sre)

#Salvando o espaço de trabalho com todos os objetos num documento RData que pode ser carregado posteriormente.
save.image(file="script3_B.RData")

# Quais foram as principais diferencas na geracao das pseudo ausencias com as diferentes estrategias?
# Explore mais opcoes, alterando o numero de repeticoes, numero de pseudo-ausencia, distancia minima e quantil do SRE

## Fim do script3 ##
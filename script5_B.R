# IBE 875 - Modelagem de Distribuição de Espécies
# PPGE/PPGBio
# Professores: Rodrigo Tardin, Maria Lucia Lorini, Mariana Vasconcellos
# Script 5 - Algoritmos de classificação e aprendizado de máquina

###############################################################
# Script 5 - Algoritmos classificacao, aprendizado de maquina #
###############################################################

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

#### Checando as opções 'default' de cada algoritmo ####
BIOMOD_ModelingOptions()

# Modificando as opções do Algoritmo GBM em relaçao a fração aleatória dos dados usados para ajustar cada árvore
modelo_op2=BIOMOD_ModelingOptions(GBM = list( distribution = 'bernoulli',
                                              n.trees = 2500,
                                              interaction.depth = 7,
                                              n.minobsinnode = 5,
                                              shrinkage = 0.001,
                                              bag.fraction = 0.7,
                                              train.fraction = 1,
                                              cv.folds = 3,
                                              keep.data = FALSE,
                                              verbose = FALSE,
                                              perf.method = 'cv',
                                              n.cores = 1),
                                  MAXENT.Phillips = list( path_to_maxent.jar = 'D:/MDE',
                                              memory_allocated = 512,
                                              background_data_dir = 'default',
                                              maximumbackground = 'default',
                                              maximumiterations = 200,
                                              visible = FALSE,
                                              linear = TRUE,
                                              quadratic = TRUE,
                                              product = TRUE,
                                              threshold = TRUE,
                                              hinge = FALSE,
                                              lq2lqptthreshold = 80,
                                              l2lqthreshold = 10,
                                              hingethreshold = 15,
                                              beta_threshold = -1,
                                              beta_categorical = -1,
                                              beta_lqp = -1,
                                              beta_hinge = -1,
                                              betamultiplier = 1,
                                              defaultprevalence = 0.5))

#### Formatando os dados, camadas e fazendo conexão para o formato do Maxent ####
# Para trabalhar com o MAXENT nós temos que criar um caminho para um diretório contendo nossas variaveis explanatórias em ascii

maxent.background.dat.dir <- "maxent_bg"
dir.create(maxent.background.dat.dir, showWarnings = FALSE, recursive = TRUE)

## salvando as variáveis preditoras em .ascii
for(var_ in names(biostack2)){
  cat("\n> saving", paste0(var_, ".asc"))
  writeRaster(subset(biostack2, var_), 
                      filename = file.path(maxent.background.dat.dir,
                      paste0(var_, ".asc")),
                      overwrite = TRUE)
}

## definindo o camimnho para o arquivo maxent.jar
path.to.maxent.jar <- file.path(getwd(), "maxent.jar")

#### Definindo os parâmetros dos modelos, incluindo a escolha dos algoritmos, o numero de rodadas e a divisão do conjunto de dados entre treino e teste. ####

lusitanic2model = BIOMOD_Modeling(
  bm.lusitanic100disk, #objeto das pseudo-ausências
  models = c("GBM", "RF","MAXENT.Phillips"), #algoritmos de escolha
  models.options = modelo_op2, #opções personalizadas ou padrões dos algoritmos
  NbRunEval = 3, #Numero de rodadas
  DataSplit = 70, #Divisão treino/teste
  Prevalence = 0.5,
  VarImport = 3, #permutações para gerar o valor de importância das variáveis
  models.eval.meth = c("TSS", "ROC"), #método de avaliação do desempenho dos modelos
  SaveObj = TRUE, #se os modelos serão salvos ou não
  rescal.all.models = F,
  do.full.models = FALSE,
  modeling.id = "lusitanic_classif_maxent")

#Sumário do objeto com os modelos criados, onde é possível ver quais modelos foram gerados (set de pseudoausencia + rodada + algoritmo).
lusitanic2model

#### Obtendo a importância de cada variável usada nos modelos ####
var_import_lusitanic2=get_variables_importance(lusitanic2model) 

# obtendo os valores de importância de cada variável para cada modelo
var_import_lusitanic2

# obtendo os valores médios de importância de cada variável para cada algoritmo
var_import_lusitanic2=apply(var_import_lusitanic2,c(1,2),mean)

# Salvando os valores de importancia das variáveis em um arquivo .csv
write.csv(var_import_lusitanic2,paste0("./Outputs/", "var_import_lusitanic2.csv"))

#Qual variável foi a mais importante para cada algoritmo? Houveram diferenças entre algoritmos?

#### Criando curvas de resposta para cada algoritmo ####

#Carregando os modelos individuais que foram gerados acuma
lusitanic2_gbm=BIOMOD_LoadModels(lusitanic2model, models = 'GBM')
lusitanic2_rf=BIOMOD_LoadModels(lusitanic2model, models = 'RF')
lusitanic2_MaxEnt=BIOMOD_LoadModels(lusitanic2model, models = 'MAXENT.Phillips')

#Plotando as curvas de resposta para cada modelo para cada variável, incluindo todas as rodadas.

biomod2::response.plot2(models = lusitanic2_gbm, Data = get_formal_data(lusitanic2model, 'expl.var'),show.variables = get_formal_data(lusitanic2model, 'expl.var.names'), do.bivariate = F, save.file = 'jpeg', ImageSize = 720, name= "Outputs/GBM_curva_resposta", fixed.var.metric = 'median', legend = F, display_title = T, data_species = get_formal_data(lusitanic2model, 'resp.var') , plot = T)

biomod2::response.plot2(models = lusitanic2_rf, Data = get_formal_data(lusitanic2model, 'expl.var'),show.variables = get_formal_data(lusitanic2model, 'expl.var.names'), do.bivariate = F, save.file = 'jpeg', ImageSize = 720, name= "Outputs/RF_curva_resposta", fixed.var.metric = 'median', legend = F, display_title = T, data_species = get_formal_data(lusitanic2model, 'resp.var') )

biomod2::response.plot2(models = lusitanic2_MaxEnt, Data = get_formal_data(lusitanic2model, 'expl.var'),show.variables = get_formal_data(lusitanic2model, 'expl.var.names'), do.bivariate = F, save.file = 'jpeg', ImageSize = 720, name= "Outputs/MaxEnt_curva_resposta", fixed.var.metric = 'median', legend = F, display_title = T, data_species = get_formal_data(lusitanic2model, 'resp.var'))

#Salvando o espaço de trabalho com todos os objetos num documento RData que pode ser carregado posteriormente.
save.image(file="script5_B.RData")

# Fim do Script 5
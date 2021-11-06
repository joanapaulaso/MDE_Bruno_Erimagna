# IBE 875 - Modelagem de Distribuicao de Especies
# PPGE/PPGBio
# Professores: Rodrigo Tardin, Maria Lucia Lorini, Mariana Vasconcellos
# Script 4 - Algoritmos de envelopes climaticos e regressao

##########################################################
# Script 4 - Algoritmos envelopes climaticos e regressao #
##########################################################

# Geracao dos modelos usando os algoritmos de envelope (SRE) e regressao  (GLM e GAM)

#### Definindo area de trabalho ####
setwd("D:/MDE") #Mude para o endereco da pasta da disciplina no seu computador
getwd()

#### Carregando pacotes ####
library(biomod2)

#### Carregando os objetos do script 2 e 3 ####
load("script2_B.RData")
load("script3_B.RData")

#### Criando um diretorio 'Outputs' para guardar tabelas com resultados das modelagens ####
dir.create("Outputs")

#### Checando as opcoes 'default' de cada algoritmo ####
BIOMOD_ModelingOptions()

# Modificando as opcoes do Algoritmo GAM quanto ao numero de nós (k = 4) para evitar sobreajuste 
modelo_op <- BIOMOD_ModelingOptions(GAM = list( algo = 'GAM_mgcv', type = 's_smoother', k = 4, interaction.level = 0, myFormula = NULL,  family = binomial(link = 'logit'), method = 'GCV.Cp', optimizer = c('outer','newton'), select = FALSE, knots = NULL,   paraPen = NULL, control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07, maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15, rank.tol = 1.49011611938477e-08, nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0), optim = list(factr=1e+07), newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0), outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE
, efs.lspmax = 15, efs.tol = 0.1, keepData = FALSE , scale.est = "fletcher", edge.correct = FALSE) ))

#### Definindo os parametros dos modelos, incluindo a escolha dos algoritmos, o numero de rodadas e a divisao do conjunto de dados entre treino e teste. ####

lusitanic1model = BIOMOD_Modeling(
  bm.lusitanic100disk, #objeto das pseudo-ausencias
  models = c("GLM", "GAM", "SRE"), #algoritmos de escolha
  models.options = modelo_op, #opções personalizadas ou padroes dos algoritmos
  NbRunEval = 3, #Numero de rodadas
  DataSplit = 70, #Divisao treino/teste
  Prevalence = 0.5,
  VarImport = 3, #permutações para gerar o valor de importancia das variaveis
  models.eval.meth = c("TSS", "ROC"), #metodo de avaliação do desempenho dos modelos
  SaveObj = TRUE, #se os modelos serao salvos ou não
  rescal.all.models = F,
  do.full.models = FALSE,
  modeling.id = "lusitanic_envelope_regres")


#Sumario do objeto com os modelos criados, onde é possível ver quais modelos foram gerados (conjunto de pseudoausencia + rodada + algoritmo).
lusitanic1model

#### Obtendo a importancia de cada variavel usada nos modelos ####
var_import_lusitanic1=get_variables_importance(lusitanic1model) 

# obtendo os valores de importancia de cada variavel para cada modelo
var_import_lusitanic1

# obtendo os valores medios de importancia de cada variavel para cada algoritmo
var_import_lusitanic1=apply(var_import_lusitanic1, c(1,2), mean)

# obtendo os valores medios de importancia de cada variavel para cada algoritmo
var_import_lusitanic1

# Salvando os valores de importancia das variaveis em um arquivo .csv
write.csv(var_import_lusitanic1,paste0("./Outputs/", "var_import_lusitanic1.csv"))

#Qual variavel foi a mais importante para cada algoritmo? Houve diferenças entre algoritmos?

#### Criando curvas de resposta para cada algoritmo ####

#Carregando os modelos individuais que foram gerados acima
lusitanic1_glm=BIOMOD_LoadModels(lusitanic1model, models = 'GLM')
lusitanic1_gam=BIOMOD_LoadModels(lusitanic1model, models = 'GAM')
lusitanic1_sre=BIOMOD_LoadModels(lusitanic1model, models = 'SRE')

#Plotando as curvas de resposta para cada modelo para cada variavel, incluindo todas as rodadas.

biomod2::response.plot2(models = lusitanic1_glm, Data = get_formal_data(lusitanic1model, 'expl.var'),show.variables = get_formal_data(lusitanic1model, 'expl.var.names'), do.bivariate = F, name= "GLM_curva_resposta", fixed.var.metric = 'median', legend = F, display_title = T, data_species = get_formal_data(lusitanic1model, 'resp.var') )

biomod2::response.plot2(models = lusitanic1_gam, Data = get_formal_data(lusitanic1model, 'expl.var'),show.variables = get_formal_data(lusitanic1model, 'expl.var.names'), do.bivariate = F, name= "GAM_curva_resposta", fixed.var.metric = 'median', legend = F, display_title = T, data_species = get_formal_data(lusitanic1model, 'resp.var') )

biomod2::response.plot2(models = lusitanic1_sre, Data = get_formal_data(lusitanic1model, 'expl.var'),show.variables = get_formal_data(lusitanic1model, 'expl.var.names'), do.bivariate = F, name= "SRE_curva_resposta", fixed.var.metric = 'median', legend = F, display_title = T, data_species = get_formal_data(lusitanic1model, 'resp.var') )

#Salvando o espaço de trabalho com todos os objetos num documento RData que pode ser carregado posteriormente.
save.image(file="script4_B.RData")

# Fim do script 4
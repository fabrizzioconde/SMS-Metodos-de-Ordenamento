# Nome dos pacotes
packages <- c("randomForest", "care","dplyr")

# Instalando os pacotes que ainda nao foram instalados
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Carregando os pacotes
invisible(lapply(packages, library, character.only = TRUE))


# Importando a base de dados simulada 1
nome_do_arquivo<-"dados_simul_1.csv"
dados_1 <- read.csv(nome_do_arquivo)
View(dados_1)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df<-dados_1


#Rank do Valor p

valor.p<-function(genotipo_fenotipo)
{
  valor_p_bruto<-vector();
  valor_p_ajustado<-vector();
  y<-genotipo_fenotipo[,ncol(genotipo_fenotipo)];
  m<-ncol(genotipo_fenotipo)-1;
  model_regression<-list();
  saida<-data.frame();
  
  for (i in 1:(ncol(genotipo_fenotipo)-1))
  {
    x<-genotipo_fenotipo[,i];
    model_regression[[i]]<-lm(y~x,data=genotipo_fenotipo);
    valor_p_bruto[i]<-ifelse(all(x==3)|all(x==2)|all(x==1),1,summary(model_regression[[i]])[[4]][2,4]);
    
  }
  valor_p_ajustado<- m*valor_p_bruto;  
  valor_p<-data.frame();
  valor_p<-cbind(valor_p_bruto,valor_p_ajustado);
  colnames(valor_p)<-c("Valor p bruto","Valor p ajustado")
  rownames(valor_p)<-names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo)-1)];
  return(valor_p)
}

valor_p<-as.data.frame(valor.p(df))
valor_p_order<-valor_p[order(valor_p[,1], decreasing = FALSE),]


# Rank ISCORE

# Função para calcular o I score de uma única variável preditiva
calculate_I_score_single <- function(data, response_var, predictor_var) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  # Criar partições baseadas na variável preditora
  partitions <- data %>%
    group_by(!!sym(predictor_var)) %>%
    summarise(n_j = n(),
              y_j_mean = mean(get(response_var)))
  
  # Calcular I score
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Função para calcular o I score de todas as variáveis preditivas
calculate_all_I_scores <- function(data, response_var) {
  predictor_vars <- setdiff(names(data), response_var)
  I_scores <- sapply(predictor_vars, function(var) calculate_I_score_single(data, response_var, var))
  return(I_scores)
}


#Rank da Random Forest
ntree<-4000 #Numero de arvores dentro da floresta aleatoria.

set.seed(1) #Semente aleatoria para floresta aleatoria.

RF <- randomForest(fenotipo~., data=df,
                   ntree=ntree,
                   mtry=ncol(df)-1,
                   importance=TRUE)    

rank_RF<-sort(importance(RF)[,1],decreasing=TRUE)

#Rank do CAR Score

dados <- list() # Criando uma lista de dataframes
dados$x <- df[-length(names(df))] # Criando o primeiro elemento (x) da lista
dados$y <- df$fenotipo # Criando o primeiro elemento (y) da lista
xnames = names(df[-length(names(df))])

# Base de dados Completa com o pacote CAR
car = carscore(dados$x, dados$y, diagonal = FALSE) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_car <- xnames[order(car^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_car # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE) # Correlação Marginal.
rank_mcorr <- xnames[order(mcorr^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_mcorr # Rank dos SNPs ordenados de forma decrescente pelo Correlação Marginal.

# Rank da Correlação Parcial
# Base de dados Completa com o pacote CAR
pcor = pcor.shrink(cbind(dados$x,dados$y), lambda=0, verbose=FALSE)[-1,1]
rank_pcorr <- xnames[order(pcor^2, decreasing=TRUE)]
rank_pcorr

# Rank ISCORE
# Base de dados Completa com o pacote I SCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_I_scores_sorted<- sort(I_scores, decreasing = TRUE) # Ordena de forma decrescente o I SCORE
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores) # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.


# Construindo o Dataframe com todos os Ranks
rank_todos_1<-cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_1)<-c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_1)

# Exportando o dataframe com os primeiros 20 SNPs para todos Ranks
write.csv(rank_todos_1[1:20,], file = "Ranks_dados_1.csv", row.names = FALSE)

# Instalar e carregar o pacote xtable
if (!require("xtable")) {
  install.packages("xtable")
  library(xtable)
}

# Selecionando apenas as 20 primeiras linhas do dataframe
rank_todos_1_top20 <- head(rank_todos_1, 20)

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(rank_todos_1_top20, caption = "Rank dos SNPs por seis métodos de ornamento para os 20 primeiros SNPs na simulação 1.", label = "tab:ranks")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Ranks_dados_1.tex", include.rownames = FALSE)


################# Base de Dados Simulada 2 ###################
# Importando a base de dados simulada 2
nome_do_arquivo<-"dados_simul_2.csv"
dados_2 <- read.csv(nome_do_arquivo)
View(dados_2)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df<-dados_2


#Rank do Valor p

valor.p<-function(genotipo_fenotipo)
{
  valor_p_bruto<-vector();
  valor_p_ajustado<-vector();
  y<-genotipo_fenotipo[,ncol(genotipo_fenotipo)];
  m<-ncol(genotipo_fenotipo)-1;
  model_regression<-list();
  saida<-data.frame();
  
  for (i in 1:(ncol(genotipo_fenotipo)-1))
  {
    x<-genotipo_fenotipo[,i];
    model_regression[[i]]<-lm(y~x,data=genotipo_fenotipo);
    valor_p_bruto[i]<-ifelse(all(x==3)|all(x==2)|all(x==1),1,summary(model_regression[[i]])[[4]][2,4]);
    
  }
  valor_p_ajustado<- m*valor_p_bruto;  
  valor_p<-data.frame();
  valor_p<-cbind(valor_p_bruto,valor_p_ajustado);
  colnames(valor_p)<-c("Valor p bruto","Valor p ajustado")
  rownames(valor_p)<-names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo)-1)];
  return(valor_p)
}

valor_p<-as.data.frame(valor.p(df))
valor_p_order<-valor_p[order(valor_p[,1], decreasing = FALSE),]


# Rank ISCORE

# Função para calcular o I score de uma única variável preditiva
calculate_I_score_single <- function(data, response_var, predictor_var) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  # Criar partições baseadas na variável preditora
  partitions <- data %>%
    group_by(!!sym(predictor_var)) %>%
    summarise(n_j = n(),
              y_j_mean = mean(get(response_var)))
  
  # Calcular I score
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Função para calcular o I score de todas as variáveis preditivas
calculate_all_I_scores <- function(data, response_var) {
  predictor_vars <- setdiff(names(data), response_var)
  I_scores <- sapply(predictor_vars, function(var) calculate_I_score_single(data, response_var, var))
  return(I_scores)
}


#Rank da Random Forest
ntree<-4000 #Numero de arvores dentro da floresta aleatoria.

set.seed(1) #Semente aleatoria para floresta aleatoria.

RF <- randomForest(fenotipo~., data=df,
                   ntree=ntree,
                   mtry=ncol(df)-1,
                   importance=TRUE)    

rank_RF<-sort(importance(RF)[,1],decreasing=TRUE)

#Rank do CAR Score

dados <- list() # Criando uma lista de dataframes
dados$x <- df[-length(names(df))] # Criando o primeiro elemento (x) da lista
dados$y <- df$fenotipo # Criando o primeiro elemento (y) da lista
xnames = names(df[-length(names(df))])

# Base de dados Completa com o pacote CAR
car = carscore(dados$x, dados$y, diagonal = FALSE) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_car <- xnames[order(car^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_car # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE) # Correlação Marginal.
rank_mcorr <- xnames[order(mcorr^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_mcorr # Rank dos SNPs ordenados de forma decrescente pelo Correlação Marginal.

# Rank da Correlação Parcial
# Base de dados Completa com o pacote CAR
pcor = pcor.shrink(cbind(dados$x,dados$y), lambda=0, verbose=FALSE)[-1,1]
rank_pcorr <- xnames[order(pcor^2, decreasing=TRUE)]
rank_pcorr

# Rank ISCORE
# Base de dados Completa com o pacote I SCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_I_scores_sorted<- sort(I_scores, decreasing = TRUE) # Ordena de forma decrescente o I SCORE
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores) # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.


# Construindo o Dataframe com todos os Ranks
rank_todos_2<-cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_2)<-c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_2)

# Exportando o dataframe com os primeiros 20 SNPs para todos Ranks
write.csv(rank_todos_2[1:20,], file = "Ranks_dados_2.csv", row.names = FALSE)

# Instalar e carregar o pacote xtable
if (!require("xtable")) {
  install.packages("xtable")
  library(xtable)
}

# Selecionando apenas as 20 primeiras linhas do dataframe
rank_todos_2_top20 <- head(rank_todos_2, 20)

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(rank_todos_2_top20, caption = "Rank dos SNPs por seis métodos de ornamento para os 20 primeiros SNPs na simulação 2.", label = "tab:ranks")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Ranks_dados_2.tex", include.rownames = FALSE)



################# Base de Dados Simulada 3 ###################
# Importando a base de dados simulada 3
nome_do_arquivo<-"dados_simul_3.csv"
dados_3 <- read.csv(nome_do_arquivo)
View(dados_3)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df<-dados_3


#Rank do Valor p

valor.p<-function(genotipo_fenotipo)
{
  valor_p_bruto<-vector();
  valor_p_ajustado<-vector();
  y<-genotipo_fenotipo[,ncol(genotipo_fenotipo)];
  m<-ncol(genotipo_fenotipo)-1;
  model_regression<-list();
  saida<-data.frame();
  
  for (i in 1:(ncol(genotipo_fenotipo)-1))
  {
    x<-genotipo_fenotipo[,i];
    model_regression[[i]]<-lm(y~x,data=genotipo_fenotipo);
    valor_p_bruto[i]<-ifelse(all(x==3)|all(x==2)|all(x==1),1,summary(model_regression[[i]])[[4]][2,4]);
    
  }
  valor_p_ajustado<- m*valor_p_bruto;  
  valor_p<-data.frame();
  valor_p<-cbind(valor_p_bruto,valor_p_ajustado);
  colnames(valor_p)<-c("Valor p bruto","Valor p ajustado")
  rownames(valor_p)<-names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo)-1)];
  return(valor_p)
}

valor_p<-as.data.frame(valor.p(df))
valor_p_order<-valor_p[order(valor_p[,1], decreasing = FALSE),]


# Rank ISCORE

# Função para calcular o I score de uma única variável preditiva
calculate_I_score_single <- function(data, response_var, predictor_var) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  # Criar partições baseadas na variável preditora
  partitions <- data %>%
    group_by(!!sym(predictor_var)) %>%
    summarise(n_j = n(),
              y_j_mean = mean(get(response_var)))
  
  # Calcular I score
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Função para calcular o I score de todas as variáveis preditivas
calculate_all_I_scores <- function(data, response_var) {
  predictor_vars <- setdiff(names(data), response_var)
  I_scores <- sapply(predictor_vars, function(var) calculate_I_score_single(data, response_var, var))
  return(I_scores)
}


#Rank da Random Forest
ntree<-4000 #Numero de arvores dentro da floresta aleatoria.

set.seed(1) #Semente aleatoria para floresta aleatoria.

RF <- randomForest(fenotipo~., data=df,
                   ntree=ntree,
                   mtry=ncol(df)-1,
                   importance=TRUE)    

rank_RF<-sort(importance(RF)[,1],decreasing=TRUE)

#Rank do CAR Score

dados <- list() # Criando uma lista de dataframes
dados$x <- df[-length(names(df))] # Criando o primeiro elemento (x) da lista
dados$y <- df$fenotipo # Criando o primeiro elemento (y) da lista
xnames = names(df[-length(names(df))])

# Base de dados Completa com o pacote CAR
car = carscore(dados$x, dados$y, diagonal = FALSE) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_car <- xnames[order(car^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_car # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE) # Correlação Marginal.
rank_mcorr <- xnames[order(mcorr^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_mcorr # Rank dos SNPs ordenados de forma decrescente pelo Correlação Marginal.

# Rank da Correlação Parcial
# Base de dados Completa com o pacote CAR
pcor = pcor.shrink(cbind(dados$x,dados$y), lambda=0, verbose=FALSE)[-1,1]
rank_pcorr <- xnames[order(pcor^2, decreasing=TRUE)]
rank_pcorr

# Rank ISCORE
# Base de dados Completa com o pacote I SCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_I_scores_sorted<- sort(I_scores, decreasing = TRUE) # Ordena de forma decrescente o I SCORE
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores) # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.


# Construindo o Dataframe com todos os Ranks
rank_todos_3<-cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_3)<-c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_3)

# Exportando o dataframe com os primeiros 20 SNPs para todos Ranks
write.csv(rank_todos_3[1:20,], file = "Ranks_dados_3.csv", row.names = FALSE)

# Instalar e carregar o pacote xtable
if (!require("xtable")) {
  install.packages("xtable")
  library(xtable)
}

# Selecionando apenas as 20 primeiras linhas do dataframe
rank_todos_3_top20 <- head(rank_todos_3, 20)

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(rank_todos_3_top20, caption = "Rank dos SNPs por seis métodos de ornamento para os 20 primeiros SNPs na simulação 3.", label = "tab:ranks")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Ranks_dados_3.tex", include.rownames = FALSE)


################# Base de Dados Simulada 4 ###################
# Importando a base de dados simulada 4
nome_do_arquivo<-"dados_simul_4.csv"
dados_4 <- read.csv(nome_do_arquivo)
View(dados_4)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df<-dados_4


#Rank do Valor p

valor.p<-function(genotipo_fenotipo)
{
  valor_p_bruto<-vector();
  valor_p_ajustado<-vector();
  y<-genotipo_fenotipo[,ncol(genotipo_fenotipo)];
  m<-ncol(genotipo_fenotipo)-1;
  model_regression<-list();
  saida<-data.frame();
  
  for (i in 1:(ncol(genotipo_fenotipo)-1))
  {
    x<-genotipo_fenotipo[,i];
    model_regression[[i]]<-lm(y~x,data=genotipo_fenotipo);
    valor_p_bruto[i]<-ifelse(all(x==3)|all(x==2)|all(x==1),1,summary(model_regression[[i]])[[4]][2,4]);
    
  }
  valor_p_ajustado<- m*valor_p_bruto;  
  valor_p<-data.frame();
  valor_p<-cbind(valor_p_bruto,valor_p_ajustado);
  colnames(valor_p)<-c("Valor p bruto","Valor p ajustado")
  rownames(valor_p)<-names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo)-1)];
  return(valor_p)
}

valor_p<-as.data.frame(valor.p(df))
valor_p_order<-valor_p[order(valor_p[,1], decreasing = FALSE),]


# Rank ISCORE

# Função para calcular o I score de uma única variável preditiva
calculate_I_score_single <- function(data, response_var, predictor_var) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  # Criar partições baseadas na variável preditora
  partitions <- data %>%
    group_by(!!sym(predictor_var)) %>%
    summarise(n_j = n(),
              y_j_mean = mean(get(response_var)))
  
  # Calcular I score
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Função para calcular o I score de todas as variáveis preditivas
calculate_all_I_scores <- function(data, response_var) {
  predictor_vars <- setdiff(names(data), response_var)
  I_scores <- sapply(predictor_vars, function(var) calculate_I_score_single(data, response_var, var))
  return(I_scores)
}


#Rank da Random Forest
ntree<-4000 #Numero de arvores dentro da floresta aleatoria.

set.seed(1) #Semente aleatoria para floresta aleatoria.

RF <- randomForest(fenotipo~., data=df,
                   ntree=ntree,
                   mtry=ncol(df)-1,
                   importance=TRUE)    

rank_RF<-sort(importance(RF)[,1],decreasing=TRUE)

#Rank do CAR Score

dados <- list() # Criando uma lista de dataframes
dados$x <- df[-length(names(df))] # Criando o primeiro elemento (x) da lista
dados$y <- df$fenotipo # Criando o primeiro elemento (y) da lista
xnames = names(df[-length(names(df))])

# Base de dados Completa com o pacote CAR
car = carscore(dados$x, dados$y, diagonal = FALSE) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_car <- xnames[order(car^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_car # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE) # Correlação Marginal.
rank_mcorr <- xnames[order(mcorr^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_mcorr # Rank dos SNPs ordenados de forma decrescente pelo Correlação Marginal.

# Rank da Correlação Parcial
# Base de dados Completa com o pacote CAR
pcor = pcor.shrink(cbind(dados$x,dados$y), lambda=0, verbose=FALSE)[-1,1]
rank_pcorr <- xnames[order(pcor^2, decreasing=TRUE)]
rank_pcorr

# Rank ISCORE
# Base de dados Completa com o pacote I SCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_I_scores_sorted<- sort(I_scores, decreasing = TRUE) # Ordena de forma decrescente o I SCORE
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores) # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.


# Construindo o Dataframe com todos os Ranks
rank_todos_4<-cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_4)<-c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_4)

# Exportando o dataframe com os primeiros 20 SNPs para todos Ranks
write.csv(rank_todos_4[1:20,], file = "Ranks_dados_4.csv", row.names = FALSE)

# Instalar e carregar o pacote xtable
if (!require("xtable")) {
  install.packages("xtable")
  library(xtable)
}

# Selecionando apenas as 20 primeiras linhas do dataframe
rank_todos_4_top20 <- head(rank_todos_4, 20)

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(rank_todos_4_top20, caption = "Rank dos SNPs por seis métodos de ornamento para os 20 primeiros SNPs na simulação 1.", label = "tab:ranks")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Ranks_dados_4.tex", include.rownames = FALSE)



################# Base de Dados Simulada 5 ###################
# Importando a base de dados simulada 5
nome_do_arquivo<-"dados_simul_5.csv"
dados_5 <- read.csv(nome_do_arquivo)
View(dados_5)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df<-dados_5


#Rank do Valor p

valor.p<-function(genotipo_fenotipo)
{
  valor_p_bruto<-vector();
  valor_p_ajustado<-vector();
  y<-genotipo_fenotipo[,ncol(genotipo_fenotipo)];
  m<-ncol(genotipo_fenotipo)-1;
  model_regression<-list();
  saida<-data.frame();
  
  for (i in 1:(ncol(genotipo_fenotipo)-1))
  {
    x<-genotipo_fenotipo[,i];
    model_regression[[i]]<-lm(y~x,data=genotipo_fenotipo);
    valor_p_bruto[i]<-ifelse(all(x==3)|all(x==2)|all(x==1),1,summary(model_regression[[i]])[[4]][2,4]);
    
  }
  valor_p_ajustado<- m*valor_p_bruto;  
  valor_p<-data.frame();
  valor_p<-cbind(valor_p_bruto,valor_p_ajustado);
  colnames(valor_p)<-c("Valor p bruto","Valor p ajustado")
  rownames(valor_p)<-names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo)-1)];
  return(valor_p)
}

valor_p<-as.data.frame(valor.p(df))
valor_p_order<-valor_p[order(valor_p[,1], decreasing = FALSE),]


# Rank ISCORE

# Função para calcular o I score de uma única variável preditiva
calculate_I_score_single <- function(data, response_var, predictor_var) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  # Criar partições baseadas na variável preditora
  partitions <- data %>%
    group_by(!!sym(predictor_var)) %>%
    summarise(n_j = n(),
              y_j_mean = mean(get(response_var)))
  
  # Calcular I score
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Função para calcular o I score de todas as variáveis preditivas
calculate_all_I_scores <- function(data, response_var) {
  predictor_vars <- setdiff(names(data), response_var)
  I_scores <- sapply(predictor_vars, function(var) calculate_I_score_single(data, response_var, var))
  return(I_scores)
}


#Rank da Random Forest
ntree<-4000 #Numero de arvores dentro da floresta aleatoria.

set.seed(1) #Semente aleatoria para floresta aleatoria.

RF <- randomForest(fenotipo~., data=df,
                   ntree=ntree,
                   mtry=ncol(df)-1,
                   importance=TRUE)    

rank_RF<-sort(importance(RF)[,1],decreasing=TRUE)

#Rank do CAR Score

dados <- list() # Criando uma lista de dataframes
dados$x <- df[-length(names(df))] # Criando o primeiro elemento (x) da lista
dados$y <- df$fenotipo # Criando o primeiro elemento (y) da lista
xnames = names(df[-length(names(df))])

# Base de dados Completa com o pacote CAR
car = carscore(dados$x, dados$y, diagonal = FALSE) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_car <- xnames[order(car^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_car # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE) # Correlação Marginal.
rank_mcorr <- xnames[order(mcorr^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_mcorr # Rank dos SNPs ordenados de forma decrescente pelo Correlação Marginal.

# Rank da Correlação Parcial
# Base de dados Completa com o pacote CAR
pcor = pcor.shrink(cbind(dados$x,dados$y), lambda=0, verbose=FALSE)[-1,1]
rank_pcorr <- xnames[order(pcor^2, decreasing=TRUE)]
rank_pcorr

# Rank ISCORE
# Base de dados Completa com o pacote I SCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_I_scores_sorted<- sort(I_scores, decreasing = TRUE) # Ordena de forma decrescente o I SCORE
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores) # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.


# Construindo o Dataframe com todos os Ranks
rank_todos_5<-cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_5)<-c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_5)

# Exportando o dataframe com os primeiros 20 SNPs para todos Ranks
write.csv(rank_todos_5[1:20,], file = "Ranks_dados_5.csv", row.names = FALSE)

# Instalar e carregar o pacote xtable
if (!require("xtable")) {
  install.packages("xtable")
  library(xtable)
}

# Selecionando apenas as 20 primeiras linhas do dataframe
rank_todos_5_top20 <- head(rank_todos_5, 20)

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(rank_todos_5_top20, caption = "Rank dos SNPs por seis métodos de ornamento para os 20 primeiros SNPs na simulação 5.", label = "tab:ranks")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Ranks_dados_5.tex", include.rownames = FALSE)



################# Base de Dados Simulada 6 ###################
# Importando a base de dados simulada 6
nome_do_arquivo<-"dados_simul_6.csv"
dados_6 <- read.csv(nome_do_arquivo)
View(dados_6)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df<-dados_6


#Rank do Valor p

valor.p<-function(genotipo_fenotipo)
{
  valor_p_bruto<-vector();
  valor_p_ajustado<-vector();
  y<-genotipo_fenotipo[,ncol(genotipo_fenotipo)];
  m<-ncol(genotipo_fenotipo)-1;
  model_regression<-list();
  saida<-data.frame();
  
  for (i in 1:(ncol(genotipo_fenotipo)-1))
  {
    x<-genotipo_fenotipo[,i];
    model_regression[[i]]<-lm(y~x,data=genotipo_fenotipo);
    valor_p_bruto[i]<-ifelse(all(x==3)|all(x==2)|all(x==1),1,summary(model_regression[[i]])[[4]][2,4]);
    
  }
  valor_p_ajustado<- m*valor_p_bruto;  
  valor_p<-data.frame();
  valor_p<-cbind(valor_p_bruto,valor_p_ajustado);
  colnames(valor_p)<-c("Valor p bruto","Valor p ajustado")
  rownames(valor_p)<-names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo)-1)];
  return(valor_p)
}

valor_p<-as.data.frame(valor.p(df))
valor_p_order<-valor_p[order(valor_p[,1], decreasing = FALSE),]


# Rank ISCORE

# Função para calcular o I score de uma única variável preditiva
calculate_I_score_single <- function(data, response_var, predictor_var) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  # Criar partições baseadas na variável preditora
  partitions <- data %>%
    group_by(!!sym(predictor_var)) %>%
    summarise(n_j = n(),
              y_j_mean = mean(get(response_var)))
  
  # Calcular I score
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Função para calcular o I score de todas as variáveis preditivas
calculate_all_I_scores <- function(data, response_var) {
  predictor_vars <- setdiff(names(data), response_var)
  I_scores <- sapply(predictor_vars, function(var) calculate_I_score_single(data, response_var, var))
  return(I_scores)
}


#Rank da Random Forest
ntree<-4000 #Numero de arvores dentro da floresta aleatoria.

set.seed(1) #Semente aleatoria para floresta aleatoria.

RF <- randomForest(fenotipo~., data=df,
                   ntree=ntree,
                   mtry=ncol(df)-1,
                   importance=TRUE)    

rank_RF<-sort(importance(RF)[,1],decreasing=TRUE)

#Rank do CAR Score

dados <- list() # Criando uma lista de dataframes
dados$x <- df[-length(names(df))] # Criando o primeiro elemento (x) da lista
dados$y <- df$fenotipo # Criando o primeiro elemento (y) da lista
xnames = names(df[-length(names(df))])

# Base de dados Completa com o pacote CAR
car = carscore(dados$x, dados$y, diagonal = FALSE) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_car <- xnames[order(car^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_car # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE) # Correlação Marginal.
rank_mcorr <- xnames[order(mcorr^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_mcorr # Rank dos SNPs ordenados de forma decrescente pelo Correlação Marginal.

# Rank da Correlação Parcial
# Base de dados Completa com o pacote CAR
pcor = pcor.shrink(cbind(dados$x,dados$y), lambda=0, verbose=FALSE)[-1,1]
rank_pcorr <- xnames[order(pcor^2, decreasing=TRUE)]
rank_pcorr

# Rank ISCORE
# Base de dados Completa com o pacote I SCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var) # Pensar depois em como criar a mesma estrutura da RF para o CARSCORE.
rank_I_scores_sorted<- sort(I_scores, decreasing = TRUE) # Ordena de forma decrescente o I SCORE
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores) # Rank dos SNPs ordenados de forma decrescente pelo CARSCORE.


# Construindo o Dataframe com todos os Ranks
rank_todos_6<-cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_6)<-c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_6)

# Exportando o dataframe com os primeiros 20 SNPs para todos Ranks
write.csv(rank_todos_6[1:20,], file = "Ranks_dados_6.csv", row.names = FALSE)

# Instalar e carregar o pacote xtable
if (!require("xtable")) {
  install.packages("xtable")
  library(xtable)
}

# Selecionando apenas as 20 primeiras linhas do dataframe
rank_todos_6_top20 <- head(rank_todos_6, 20)

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(rank_todos_6_top20, caption = "Rank dos SNPs por seis métodos de ornamento para os 20 primeiros SNPs na simulação 6.", label = "tab:ranks")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Ranks_dados_6.tex", include.rownames = FALSE)


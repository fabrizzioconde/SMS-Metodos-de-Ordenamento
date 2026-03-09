# Carregar pacotes necessários
library(dplyr)

# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Função para calcular o I score de múltiplas variáveis preditivas
calculate_I_score_multiple <- function(data, response_var, predictor_vars) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  # Criar partições baseadas nas variáveis preditoras
  partitions <- data %>%
    group_by(across(all_of(predictor_vars))) %>%
    summarise(n_j = n(),
              y_j_mean = mean(get(response_var)))
  
  # Calcular I score
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Função para calcular o I score de todas as combinações de variáveis preditivas
calculate_combinations_I_scores <- function(data, response_var, combination_size = 2) {
  predictor_vars <- setdiff(names(data), response_var)
  combinations <- combn(predictor_vars, combination_size, simplify = FALSE)
  I_scores <- sapply(combinations, function(vars) calculate_I_score_multiple(data, response_var, vars))
  names(I_scores) <- sapply(combinations, paste, collapse = "+")
  return(I_scores)
}

# I Score da Base de Dados Simulada 2
data <- read.csv(file = "dados_simul_2.csv")
#response_var <- "fenotipo"
#I_scores_oito <- calculate_combinations_I_scores(data, response_var, combination_size = 2)
#I_scores_sort_oito <- sort(I_scores_oito, decreasing = TRUE)
#View(I_scores_sort_oito)



####################################
# Carregar pacotes necessários
library(dplyr)

# Função para calcular o I score de múltiplas variáveis preditivas
calculate_I_score_multiple <- function(data, response_var, predictor_vars) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  # Criar partições baseadas nas variáveis preditivas
  partitions <- data %>%
    group_by(across(all_of(predictor_vars))) %>%
    summarise(n_j = n(),
              y_j_mean = mean(get(response_var)), .groups = 'drop')
  
  # Calcular I score
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Definir a variável resposta
response_var <- "fenotipo"

# Definir a combinação específica de oito variáveis preditivas
predictor_vars <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

# Calcular o I score para a combinação específica de oito variáveis preditivas
I_score_eight <- calculate_I_score_multiple(data, response_var, predictor_vars)

# Imprimir o I score
print(I_score_eight)

###############################

# Definir a combinação específica de oito variáveis preditivas
predictor_vars_2 <- c("SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8","SNP9")

# Calcular o I score para a combinação específica de oito variáveis preditivas
I_score_eight_2 <- calculate_I_score_multiple(data, response_var, predictor_vars_2)

# Imprimir o I score
print(I_score_eight_2)



###############################

# Definir a combinação específica de oito variáveis preditivas
predictor_vars_3 <- c("SNP1","SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8","SNP9")

# Calcular o I score para a combinação específica de oito variáveis preditivas
I_score_eight_3 <- calculate_I_score_multiple(data, response_var, predictor_vars_3)

# Imprimir o I score
print(I_score_eight_3)


###############################

# Definir a combinação específica de oito variáveis preditivas
predictor_vars_4 <- c("SNP1","SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8","SNP9","SNP20")

# Calcular o I score para a combinação específica de oito variáveis preditivas
I_score_eight_4 <- calculate_I_score_multiple(data, response_var, predictor_vars_4)

# Imprimir o I score
print(I_score_eight_4)


###############################

# Definir a combinação específica de oito variáveis preditivas
predictor_vars_5 <- c("SNP1","SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8","SNP9","SNP10","SNP11","SNP20")

# Calcular o I score para a combinação específica de oito variáveis preditivas
I_score_eight_5 <- calculate_I_score_multiple(data, response_var, predictor_vars_5)

# Imprimir o I score
print(I_score_eight_5)

###############################

# Definir a combinação específica de oito variáveis preditivas
predictor_vars_6 <- c("SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP9","SNP10","SNP13")

# Calcular o I score para a combinação específica de oito variáveis preditivas
I_score_eight_6 <- calculate_I_score_multiple(data, response_var, predictor_vars_6)

# Imprimir o I score
print(I_score_eight_6)



########################################################
# Novo código

# Carregar pacotes necessários
library(dplyr)

# Função para calcular o I score de múltiplas variáveis preditivas
calculate_I_score_multiple <- function(data, response_var, predictor_vars) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  # Criar partições baseadas nas variáveis preditivas
  partitions <- data %>%
    group_by(across(all_of(predictor_vars))) %>%
    summarise(n_j = n(),
              y_j_mean = mean(get(response_var)), .groups = 'drop')
  
  # Calcular I score
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Função para calcular o I score de todas as combinações de variáveis preditivas
calculate_combinations_I_scores <- function(data, response_var, combination_size = 8) {
  predictor_vars <- setdiff(names(data), response_var)
  combinations <- combn(predictor_vars, combination_size, simplify = FALSE)
  I_scores <- sapply(combinations, function(vars) calculate_I_score_multiple(data, response_var, vars))
  names(I_scores) <- sapply(combinations, paste, collapse = "+")
  return(I_scores)
}

# Exemplo de aplicação do método

dados <- read.csv(file = "dados_simul_2.csv")
var_sel<-c(1:20,101)
data <- dados[,var_sel]
View(data)

# Calcular o I score de todas as combinações de 8 variáveis preditivas
response_var <- "fenotipo"
I_scores_combinations_8 <- calculate_combinations_I_scores(data, response_var, combination_size = 8)
I_scores_combinations_8_sort <- sort(I_scores_oito, decreasing = TRUE)

# Imprimir os I scores
print(I_scores_combinations_8_sort)


###########################################
# Carregar pacotes necessários
library(dplyr)

# Função para calcular o I score de múltiplas variáveis preditivas
calculate_I_score_multiple <- function(data, response_var, predictor_vars) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  # Criar partições baseadas nas variáveis preditivas
  partitions <- data %>%
    group_by(across(all_of(predictor_vars))) %>%
    summarise(n_j = n(),
              y_j_mean = mean(.data[[response_var]]), .groups = 'drop')
  
  # Calcular I score
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Função para calcular o I score de todas as combinações de variáveis preditivas
calculate_combinations_I_scores <- function(data, response_var, combination_size = 8) {
  predictor_vars <- setdiff(names(data), response_var)
  combinations <- combn(predictor_vars, combination_size, simplify = FALSE)
  I_scores <- sapply(combinations, function(vars) calculate_I_score_multiple(data, response_var, vars))
  names(I_scores) <- sapply(combinations, paste, collapse = "+")
  return(I_scores)
}

# Exemplo de aplicação do método

dados <- read.csv(file = "dados_simul_2.csv")
var_sel<-c(1:20,101)
data <- dados[,var_sel]
View(data)

# Calcular o I score de todas as combinações de 8 variáveis preditivas
response_var <- "fenotipo"
I_scores_combinations_8 <- calculate_combinations_I_scores(data, response_var, combination_size = 8)
I_scores_combinations_8_sort <- sort(I_scores_oito, decreasing = TRUE)

# Imprimir os I scores
print(I_scores_combinations_8_sort)


# Carregar pacotes necessários
library(dplyr)

# Função para calcular o I score de múltiplas variáveis preditivas
calculate_I_score_multiple <- function(data, response_var, predictor_vars) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  # Criar partições baseadas nas variáveis preditivas
  partitions <- data %>%
    group_by(across(all_of(predictor_vars))) %>%
    summarise(n_j = n(),
              y_j_mean = mean(.data[[response_var]]), .groups = 'drop')
  
  # Calcular I score
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Função para calcular o I score de todas as combinações de variáveis preditivas
calculate_combinations_I_scores <- function(data, response_var, combination_size = 8) {
  predictor_vars <- setdiff(names(data), response_var)
  combinations <- combn(predictor_vars, combination_size, simplify = FALSE)
  I_scores <- sapply(combinations, function(vars) {
    I_score <- calculate_I_score_multiple(data, response_var, vars)
    combination_name <- paste(vars, collapse = "+")
    cat("Combinação:", combination_name, "- I score:", I_score, "\n")
    return(I_score)
  })
  names(I_scores) <- sapply(combinations, paste, collapse = "+")
  return(I_scores)
}

# Exemplo de aplicação do método

dados <- read.csv(file = "dados_simul_2.csv")
var_sel<-c(1:20,101)
data <- dados[,var_sel]
View(data)

# Calcular o I score de todas as combinações de 8 variáveis preditivas
response_var <- "fenotipo"
I_scores_combinations_8 <- calculate_combinations_I_scores(data, response_var, combination_size = 8)
I_scores_combinations_8_sort <- sort(I_scores_combinations_8, decreasing = TRUE)
View(I_scores_combinations_8_sort)
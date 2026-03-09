# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Carregar pacotes necessários
library(dplyr)

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

# Exemplo de aplicação do método
# Gerar dados fictícios
set.seed(123)
n <- 200
data <- data.frame(
  response = rnorm(n),
  predictor1 = rnorm(n),
  predictor2 = rnorm(n),
  predictor3 = rnorm(n),
  predictor4 = rnorm(n)
)

# Calcular o I score de todas as variáveis preditivas
response_var <- "response"
I_scores <- calculate_all_I_scores(data, response_var)

# Imprimir os I scores
print(I_scores)

# I Score da Base de Dados Simulada 1
data <- read.csv(file = "dados_simul_1.csv")
response_var <- "fenotipo"
I_scores_1 <- calculate_all_I_scores(data, response_var)
I_scores_sort_1 <- sort(I_scores_1, decreasing = TRUE)
View(I_scores_sort_1)

# I Score da Base de Dados Simulada 2
data <- read.csv(file = "dados_simul_2.csv")
response_var <- "fenotipo"
I_scores_2 <- calculate_all_I_scores(data, response_var)
I_scores_sort_2 <- sort(I_scores_2, decreasing = TRUE)
View(I_scores_sort_2)

# I Score da Base de Dados Simulada 3
data <- read.csv(file = "dados_simul_3.csv")
response_var <- "fenotipo"
I_scores_3 <- calculate_all_I_scores(data, response_var)
I_scores_sort_3 <- sort(I_scores_3, decreasing = TRUE)
View(I_scores_sort_3)

# I Score da Base de Dados Simulada 4
data <- read.csv(file = "dados_simul_4.csv")
response_var <- "fenotipo"
I_scores_4 <- calculate_all_I_scores(data, response_var)
I_scores_sort_4 <- sort(I_scores_4, decreasing = TRUE)
View(I_scores_sort_4)

# I Score da Base de Dados Simulada 5
data <- read.csv(file = "dados_simul_5.csv")
response_var <- "fenotipo"
I_scores_5 <- calculate_all_I_scores(data, response_var)
I_scores_sort_5 <- sort(I_scores_5, decreasing = TRUE)
View(I_scores_sort_5)

# I Score da Base de Dados Simulada 6
data <- read.csv(file = "dados_simul_6.csv")
response_var <- "fenotipo"
I_scores_6 <- calculate_all_I_scores(data, response_var)
I_scores_sort_6 <- sort(I_scores_6, decreasing = TRUE)
View(I_scores_sort_6)
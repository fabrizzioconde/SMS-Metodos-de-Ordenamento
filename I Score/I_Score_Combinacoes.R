# Carregar pacotes necessários
library(dplyr)

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

# Calcular o I score de todas as duplas de variáveis preditivas
response_var <- "response"
I_scores_duplas <- calculate_combinations_I_scores(data, response_var, combination_size = 2)

# Calcular o I score de todos os trios de variáveis preditivas
I_scores_trios <- calculate_combinations_I_scores(data, response_var, combination_size = 3)

# Imprimir os I scores
print(I_scores_duplas)
print(I_scores_trios)

# I Score da Base de Dados Simulada 2
data <- read.csv(file = "dados_simul_2.csv")
response_var <- "fenotipo"
I_scores_duplas_2 <- calculate_combinations_I_scores(data, response_var, combination_size = 2)
I_scores_sort_duplas_2 <- sort(I_scores_duplas_2, decreasing = TRUE)
View(I_scores_sort_duplas_2)
# Carregar pacotes necessários
library(caret)
library(dplyr)

# Função para calcular o I score em um problema de regressão
calculate_I_score <- function(data, response_var, predictor_vars) {
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

# Função de seleção de variáveis usando PR para regressão
partition_retention <- function(data, response_var, predictor_vars) {
  selected_vars <- predictor_vars
  best_I_score <- calculate_I_score(data, response_var, selected_vars)
  
  improving <- TRUE
  while(improving) {
    improving <- FALSE
    for (var in selected_vars) {
      current_vars <- setdiff(selected_vars, var)
      current_I_score <- calculate_I_score(data, response_var, current_vars)
      if (current_I_score > best_I_score) {
        best_I_score <- current_I_score
        selected_vars <- current_vars
        improving <- TRUE
      }
    }
  }
  return(selected_vars)
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

# Seleção de variáveis
predictor_vars <- colnames(data)[-1]  # Excluindo a variável resposta
selected_vars <- partition_retention(data, "response", predictor_vars)

# Imprimir variáveis selecionadas
print(selected_vars)

# Criar modelo de regressão linear com variáveis selecionadas
model <- lm(response ~ ., data = data[, c("response", selected_vars)])
summary(model)

data <- dados[[1]]
response_var <- 'fenotipo'
predictor_vars <- colnames(data)[-length(colnames(dados[[1]]))]

# Base de Dados Simulada 1
I_score <- calculate_I_score(data, response_var, predictor_vars)
I_score
# Nome dos pacotes
packages <- c("randomForest", "care", "dplyr", "xtable")

# Instalando os pacotes que ainda não foram instalados
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Carregando os pacotes
invisible(lapply(packages, library, character.only = TRUE))

# Importando a base de dados simulada 1
nome_do_arquivo <- "dados_simul_1.csv"
dados_1 <- read.csv(nome_do_arquivo)
View(dados_1)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df <- dados_1

# Definindo SNPs causais
snps_causais_1 <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

# Função para calcular o valor p
valor.p <- function(genotipo_fenotipo) {
  valor_p_bruto <- vector()
  valor_p_ajustado <- vector()
  y <- genotipo_fenotipo[, ncol(genotipo_fenotipo)]
  m <- ncol(genotipo_fenotipo) - 1
  model_regression <- list()
  saida <- data.frame()
  
  for (i in 1:(ncol(genotipo_fenotipo) - 1)) {
    x <- genotipo_fenotipo[, i]
    model_regression[[i]] <- lm(y ~ x, data = genotipo_fenotipo)
    valor_p_bruto[i] <- ifelse(all(x == 3) | all(x == 2) | all(x == 1), 1, summary(model_regression[[i]])[[4]][2, 4])
  }
  valor_p_ajustado <- m * valor_p_bruto  
  valor_p <- data.frame()
  valor_p <- cbind(valor_p_bruto, valor_p_ajustado)
  colnames(valor_p) <- c("Valor p bruto", "Valor p ajustado")
  rownames(valor_p) <- names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo) - 1)]
  return(valor_p)
}

valor_p <- as.data.frame(valor.p(df))
valor_p_order <- valor_p[order(valor_p[, 1], decreasing = FALSE), ]

# Função para calcular o I score de uma única variável preditiva
calculate_I_score_single <- function(data, response_var, predictor_var) {
  n <- nrow(data)
  y <- data[[response_var]]
  y_mean <- mean(y)
  s <- sd(y)
  
  partitions <- data %>%
    group_by(!!sym(predictor_var)) %>%
    summarise(n_j = n(), y_j_mean = mean(get(response_var)))
  
  I_score <- sum(partitions$n_j * (partitions$y_j_mean - y_mean)^2) / (n * s^2)
  return(I_score)
}

# Função para calcular o I score de todas as variáveis preditivas
calculate_all_I_scores <- function(data, response_var) {
  predictor_vars <- setdiff(names(data), response_var)
  I_scores <- sapply(predictor_vars, function(var) calculate_I_score_single(data, response_var, var))
  return(I_scores)
}

# Rank da Random Forest
ntree <- 4000

set.seed(1)

RF <- randomForest(fenotipo ~ ., data = df,
                   ntree = ntree,
                   mtry = ncol(df) - 1,
                   importance = TRUE)    

rank_RF <- sort(importance(RF)[, 1], decreasing = TRUE)

# Rank do CAR Score
dados <- list()
dados$x <- df[-length(names(df))]
dados$y <- df$fenotipo
xnames = names(df[-length(names(df))])

car = carscore(dados$x, dados$y, diagonal = FALSE)
rank_car <- xnames[order(car^2, decreasing = TRUE)]

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE)
rank_mcorr <- xnames[order(mcorr^2, decreasing = TRUE)]

# Rank da Correlação Parcial
pcor = pcor.shrink(cbind(dados$x, dados$y), lambda = 0, verbose = FALSE)[-1, 1]
rank_pcorr <- xnames[order(pcor^2, decreasing = TRUE)]

# Rank ISCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var)
rank_I_scores_sorted <- sort(I_scores, decreasing = TRUE)
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores)

# Construindo o Dataframe com todos os Ranks
rank_todos_1 <- cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_1) <- c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_1)

# Identificando as posições dos SNPs causais em cada ranking
positions_causais_1 <- data.frame(
  SNPs = snps_causais_1,
  Random_Forest = match(snps_causais_1, names(rank_RF)),
  Valor_p = match(snps_causais_1, rownames(valor_p_order)),
  CAR_Score = match(snps_causais_1, rank_car),
  Correlacao_Marginal = match(snps_causais_1, rank_mcorr),
  Correlacao_Parcial = match(snps_causais_1, rank_pcorr),
  I_Score = match(snps_causais_1, rank_I_scores)
)

# Garantindo a ordem correta dos SNPs causais
positions_causais_1$SNPs <- factor(positions_causais_1$SNPs, levels = snps_causais_1)
positions_causais_1 <- positions_causais_1[order(positions_causais_1$SNPs), ]

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(positions_causais_1, caption = "Posição dos SNPs causais por seis métodos de ordenamento na simulação 1.", label = "tab:positions_causais")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Positions_causais_dados_1.tex", include.rownames = FALSE)

########## Base de Dados Simulada 2 #########################

# Importando a base de dados simulada 2
nome_do_arquivo <- "dados_simul_2.csv"
dados_2 <- read.csv(nome_do_arquivo)
View(dados_2)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df <- dados_2

# Definindo SNPs causais
snps_causais_2 <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

# Função para calcular o valor p
valor_p <- as.data.frame(valor.p(df))
valor_p_order <- valor_p[order(valor_p[, 1], decreasing = FALSE), ]

# Rank da Random Forest
ntree <- 4000

set.seed(1)

RF <- randomForest(fenotipo ~ ., data = df,
                   ntree = ntree,
                   mtry = ncol(df) - 1,
                   importance = TRUE)    

rank_RF <- sort(importance(RF)[, 1], decreasing = TRUE)

# Rank do CAR Score
dados <- list()
dados$x <- df[-length(names(df))]
dados$y <- df$fenotipo
xnames = names(df[-length(names(df))])

car = carscore(dados$x, dados$y, diagonal = FALSE)
rank_car <- xnames[order(car^2, decreasing = TRUE)]

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE)
rank_mcorr <- xnames[order(mcorr^2, decreasing = TRUE)]

# Rank da Correlação Parcial
pcor = pcor.shrink(cbind(dados$x, dados$y), lambda = 0, verbose = FALSE)[-1, 1]
rank_pcorr <- xnames[order(pcor^2, decreasing = TRUE)]

# Rank ISCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var)
rank_I_scores_sorted <- sort(I_scores, decreasing = TRUE)
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores)

# Construindo o Dataframe com todos os Ranks
rank_todos_2 <- cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_2) <- c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_2)

# Identificando as posições dos SNPs causais em cada ranking
positions_causais_2 <- data.frame(
  SNPs = snps_causais_2,
  Random_Forest = match(snps_causais_2, names(rank_RF)),
  Valor_p = match(snps_causais_2, rownames(valor_p_order)),
  CAR_Score = match(snps_causais_2, rank_car),
  Correlacao_Marginal = match(snps_causais_2, rank_mcorr),
  Correlacao_Parcial = match(snps_causais_2, rank_pcorr),
  I_Score = match(snps_causais_2, rank_I_scores)
)

# Garantindo a ordem correta dos SNPs causais
positions_causais_2$SNPs <- factor(positions_causais_2$SNPs, levels = snps_causais_2)
positions_causais_2 <- positions_causais_2[order(positions_causais_2$SNPs), ]

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(positions_causais_2, caption = "Posição dos SNPs causais por seis métodos de ordenamento na simulação 2.", label = "tab:positions_causais")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Positions_causais_dados_2.tex", include.rownames = FALSE)


########## Base de Dados Simulada 3 #########################

# Importando a base de dados simulada 3
nome_do_arquivo <- "dados_simul_3.csv"
dados_3 <- read.csv(nome_do_arquivo)
View(dados_3)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df <- dados_3

# Definindo SNPs causais
snps_causais_3 <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8", "SNP9")

# Função para calcular o valor p
valor_p <- as.data.frame(valor.p(df))
valor_p_order <- valor_p[order(valor_p[, 1], decreasing = FALSE), ]

# Rank da Random Forest
ntree <- 4000

set.seed(1)

RF <- randomForest(fenotipo ~ ., data = df,
                   ntree = ntree,
                   mtry = ncol(df) - 1,
                   importance = TRUE)    

rank_RF <- sort(importance(RF)[, 1], decreasing = TRUE)

# Rank do CAR Score
dados <- list()
dados$x <- df[-length(names(df))]
dados$y <- df$fenotipo
xnames = names(df[-length(names(df))])

car = carscore(dados$x, dados$y, diagonal = FALSE)
rank_car <- xnames[order(car^2, decreasing = TRUE)]

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE)
rank_mcorr <- xnames[order(mcorr^2, decreasing = TRUE)]

# Rank da Correlação Parcial
pcor = pcor.shrink(cbind(dados$x, dados$y), lambda = 0, verbose = FALSE)[-1, 1]
rank_pcorr <- xnames[order(pcor^2, decreasing = TRUE)]

# Rank ISCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var)
rank_I_scores_sorted <- sort(I_scores, decreasing = TRUE)
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores)

# Construindo o Dataframe com todos os Ranks
rank_todos_3 <- cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_3) <- c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_3)

# Identificando as posições dos SNPs causais em cada ranking
positions_causais_3 <- data.frame(
  SNPs = snps_causais_3,
  Random_Forest = match(snps_causais_3, names(rank_RF)),
  Valor_p = match(snps_causais_3, rownames(valor_p_order)),
  CAR_Score = match(snps_causais_3, rank_car),
  Correlacao_Marginal = match(snps_causais_3, rank_mcorr),
  Correlacao_Parcial = match(snps_causais_3, rank_pcorr),
  I_Score = match(snps_causais_3, rank_I_scores)
)

# Garantindo a ordem correta dos SNPs causais
positions_causais_3$SNPs <- factor(positions_causais_3$SNPs, levels = snps_causais_3)
positions_causais_3 <- positions_causais_3[order(positions_causais_3$SNPs), ]

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(positions_causais_3, caption = "Posição dos SNPs causais por seis métodos de ordenamento na simulação 3.", label = "tab:positions_causais")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Positions_causais_dados_3.tex", include.rownames = FALSE)


########## Base de Dados Simulada 4 #########################

# Importando a base de dados simulada 4
nome_do_arquivo <- "dados_simul_4.csv"
dados_4 <- read.csv(nome_do_arquivo)
View(dados_4)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df <- dados_4

# Definindo SNPs causais
snps_causais_4 <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

# Função para calcular o valor p
valor_p <- as.data.frame(valor.p(df))
valor_p_order <- valor_p[order(valor_p[, 1], decreasing = FALSE), ]

# Rank da Random Forest
ntree <- 4000

set.seed(1)

RF <- randomForest(fenotipo ~ ., data = df,
                   ntree = ntree,
                   mtry = ncol(df) - 1,
                   importance = TRUE)    

rank_RF <- sort(importance(RF)[, 1], decreasing = TRUE)

# Rank do CAR Score
dados <- list()
dados$x <- df[-length(names(df))]
dados$y <- df$fenotipo
xnames = names(df[-length(names(df))])

car = carscore(dados$x, dados$y, diagonal = FALSE)
rank_car <- xnames[order(car^2, decreasing = TRUE)]

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE)
rank_mcorr <- xnames[order(mcorr^2, decreasing = TRUE)]

# Rank da Correlação Parcial
pcor = pcor.shrink(cbind(dados$x, dados$y), lambda = 0, verbose = FALSE)[-1, 1]
rank_pcorr <- xnames[order(pcor^2, decreasing = TRUE)]

# Rank ISCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var)
rank_I_scores_sorted <- sort(I_scores, decreasing = TRUE)
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores)

# Construindo o Dataframe com todos os Ranks
rank_todos_4 <- cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_4) <- c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_4)

# Identificando as posições dos SNPs causais em cada ranking
positions_causais_4 <- data.frame(
  SNPs = snps_causais_4,
  Random_Forest = match(snps_causais_4, names(rank_RF)),
  Valor_p = match(snps_causais_4, rownames(valor_p_order)),
  CAR_Score = match(snps_causais_4, rank_car),
  Correlacao_Marginal = match(snps_causais_4, rank_mcorr),
  Correlacao_Parcial = match(snps_causais_4, rank_pcorr),
  I_Score = match(snps_causais_4, rank_I_scores)
)

# Garantindo a ordem correta dos SNPs causais
positions_causais_4$SNPs <- factor(positions_causais_4$SNPs, levels = snps_causais_4)
positions_causais_4 <- positions_causais_4[order(positions_causais_4$SNPs), ]

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(positions_causais_4, caption = "Posição dos SNPs causais por seis métodos de ordenamento na simulação 4.", label = "tab:positions_causais")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Positions_causais_dados_4.tex", include.rownames = FALSE)


########## Base de Dados Simulada 5 #########################

# Importando a base de dados simulada 5
nome_do_arquivo <- "dados_simul_5.csv"
dados_5 <- read.csv(nome_do_arquivo)
View(dados_5)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df <- dados_5

# Definindo SNPs causais
snps_causais_5 <- c("SNP1", "SNP2", "SNP3", "SNP4")

# Função para calcular o valor p
valor_p <- as.data.frame(valor.p(df))
valor_p_order <- valor_p[order(valor_p[, 1], decreasing = FALSE), ]

# Rank da Random Forest
ntree <- 4000

set.seed(1)

RF <- randomForest(fenotipo ~ ., data = df,
                   ntree = ntree,
                   mtry = ncol(df) - 1,
                   importance = TRUE)    

rank_RF <- sort(importance(RF)[, 1], decreasing = TRUE)

# Rank do CAR Score
dados <- list()
dados$x <- df[-length(names(df))]
dados$y <- df$fenotipo
xnames = names(df[-length(names(df))])

car = carscore(dados$x, dados$y, diagonal = FALSE)
rank_car <- xnames[order(car^2, decreasing = TRUE)]

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE)
rank_mcorr <- xnames[order(mcorr^2, decreasing = TRUE)]

# Rank da Correlação Parcial
pcor = pcor.shrink(cbind(dados$x, dados$y), lambda = 0, verbose = FALSE)[-1, 1]
rank_pcorr <- xnames[order(pcor^2, decreasing = TRUE)]

# Rank ISCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var)
rank_I_scores_sorted <- sort(I_scores, decreasing = TRUE)
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores)

# Construindo o Dataframe com todos os Ranks
rank_todos_5 <- cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_5) <- c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_5)

# Identificando as posições dos SNPs causais em cada ranking
positions_causais_5 <- data.frame(
  SNPs = snps_causais_5,
  Random_Forest = match(snps_causais_5, names(rank_RF)),
  Valor_p = match(snps_causais_5, rownames(valor_p_order)),
  CAR_Score = match(snps_causais_5, rank_car),
  Correlacao_Marginal = match(snps_causais_5, rank_mcorr),
  Correlacao_Parcial = match(snps_causais_5, rank_pcorr),
  I_Score = match(snps_causais_5, rank_I_scores)
)

# Garantindo a ordem correta dos SNPs causais
positions_causais_5$SNPs <- factor(positions_causais_5$SNPs, levels = snps_causais_5)
positions_causais_5 <- positions_causais_5[order(positions_causais_5$SNPs), ]

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(positions_causais_5, caption = "Posição dos SNPs causais por seis métodos de ordenamento na simulação 5.", label = "tab:positions_causais")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Positions_causais_dados_5.tex", include.rownames = FALSE)


########## Base de Dados Simulada 6 #########################

# Importando a base de dados simulada 6
nome_do_arquivo <- "dados_simul_6.csv"
dados_6 <- read.csv(nome_do_arquivo)
View(dados_6)

# Salvando o dataframe de cada simulação em um dataframe padrão df
df <- dados_6

# Definindo SNPs causais
snps_causais_6 <- c("SNP1", "SNP10", "SNP20", "SNP30", "SNP40", "SNP50", "SNP60")

# Função para calcular o valor p
valor_p <- as.data.frame(valor.p(df))
valor_p_order <- valor_p[order(valor_p[, 1], decreasing = FALSE), ]

# Rank da Random Forest
ntree <- 4000

set.seed(1)

RF <- randomForest(fenotipo ~ ., data = df,
                   ntree = ntree,
                   mtry = ncol(df) - 1,
                   importance = TRUE)    

rank_RF <- sort(importance(RF)[, 1], decreasing = TRUE)

# Rank do CAR Score
dados <- list()
dados$x <- df[-length(names(df))]
dados$y <- df$fenotipo
xnames = names(df[-length(names(df))])

car = carscore(dados$x, dados$y, diagonal = FALSE)
rank_car <- xnames[order(car^2, decreasing = TRUE)]

# Rank da Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE)
rank_mcorr <- xnames[order(mcorr^2, decreasing = TRUE)]

# Rank da Correlação Parcial
pcor = pcor.shrink(cbind(dados$x, dados$y), lambda = 0, verbose = FALSE)[-1, 1]
rank_pcorr <- xnames[order(pcor^2, decreasing = TRUE)]

# Rank ISCORE
response_var <- "fenotipo"
I_scores <- calculate_all_I_scores(df, response_var)
rank_I_scores_sorted <- sort(I_scores, decreasing = TRUE)
rank_I_scores <- names(rank_I_scores_sorted)
View(rank_I_scores)

# Construindo o Dataframe com todos os Ranks
rank_todos_6 <- cbind(names(rank_RF), rownames(valor_p_order), rank_car, rank_mcorr, rank_pcorr, rank_I_scores)
colnames(rank_todos_6) <- c("Random Forest", "Valor p", "CAR Score", "Correlação Marginal", "Correlação Parcial", "I Score")
View(rank_todos_6)

# Identificando as posições dos SNPs causais em cada ranking
positions_causais_6 <- data.frame(
  SNPs = snps_causais_6,
  Random_Forest = match(snps_causais_6, names(rank_RF)),
  Valor_p = match(snps_causais_6, rownames(valor_p_order)),
  CAR_Score = match(snps_causais_6, rank_car),
  Correlacao_Marginal = match(snps_causais_6, rank_mcorr),
  Correlacao_Parcial = match(snps_causais_6, rank_pcorr),
  I_Score = match(snps_causais_6, rank_I_scores)
)

# Garantindo a ordem correta dos SNPs causais
positions_causais_6$SNPs <- factor(positions_causais_6$SNPs, levels = snps_causais_6)
positions_causais_6 <- positions_causais_6[order(positions_causais_6$SNPs), ]

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(positions_causais_6, caption = "Posição dos SNPs causais por seis métodos de ordenamento na simulação 6.", label = "tab:positions_causais")

# Salvar a tabela em um arquivo .tex
print(latex_table, file = "Positions_causais_dados_6.tex", include.rownames = FALSE)

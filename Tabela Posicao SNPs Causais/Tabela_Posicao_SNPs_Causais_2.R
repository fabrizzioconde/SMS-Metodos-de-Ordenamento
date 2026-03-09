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
colnames(rank_todos_1) <- c("RF", "Valor p", "CAR Score", "MAR Score", "PAR Score", "I Score")
View(rank_todos_1)

# Identificando as posições dos SNPs causais em cada ranking
positions_causais <- data.frame(
  SNPs = snps_causais_1,
  RF = match(snps_causais_1, names(rank_RF)),
  Valor_p = match(snps_causais_1, rownames(valor_p_order)),
  CAR_Score = match(snps_causais_1, rank_car),
  MAR_Score = match(snps_causais_1, rank_mcorr),
  PAR_Score = match(snps_causais_1, rank_pcorr),
  I_Score = match(snps_causais_1, rank_I_scores)
)

# Garantindo a ordem correta dos SNPs causais
positions_causais$SNPs <- factor(positions_causais$SNPs, levels = snps_causais_1)
positions_causais <- positions_causais[order(positions_causais$SNPs), ]

# Quebrar o cabeçalho em duas linhas com makecell
addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- c(
  "\\hline\n SNPs & \\makecell{RF} & \\makecell{Valor p} & \\makecell{CAR Score} & \\makecell{MAR Score} & \\makecell{PAR Score} & \\makecell{I Score} \\\\ \\hline\n"
)

# Converter o DataFrame para formato LaTeX usando xtable
latex_table <- xtable(positions_causais, caption = "Posição dos SNPs causais por seis métodos de ordenamento na simulação 1.", label = "tab:positions_causais")

# Adicionar pacote makecell e ajustar as margens ao preâmbulo do documento
cat("\\documentclass{article}\n\\usepackage{makecell}\n\\usepackage{geometry}\n\\geometry{a4paper, margin=1in}\n\\begin{document}\n", file = "Positions_causais_dados_1.tex")
print(latex_table, file = "Positions_causais_dados_1_6.tex", include.rownames = FALSE, add.to.row = addtorow, table.placement = "h", append = TRUE)
cat("\\end{document}\n", file = "Positions_causais_dados_1_6.tex", append = TRUE)

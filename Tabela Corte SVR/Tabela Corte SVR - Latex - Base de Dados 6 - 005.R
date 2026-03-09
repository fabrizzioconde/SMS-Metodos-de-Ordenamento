# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Instale e carregue os pacotes necessários
if (!require("xtable")) install.packages("xtable", dependencies = TRUE)
library(xtable)

# Criação do dataframe
df <- data.frame(
  SVR = c("SVR linear", "SVR radial $\\gamma = 0.001$", "SVR radial $\\gamma = 0.01$", "SVR radial $\\gamma = 0.1$", "SVR radial $\\gamma = 1$"),
  `Rank RF` = c(20, 40, 20, 20, 20),
  `CAR Score` = c(30, 70, 70, 20, 20),
  `MAR Score` = c(30, 70, 70, 20, 20),
  `PAR Score` = c(20, 20, 20, 60, 20),
  `I Score` = c(40, 70, 60, 20, 20)
)

# Adicionar título
title <- "Número de SNPs selecionados pela etapa de Corte para simulação 6 com percentual snps igual a 0,05."

# Converter dataframe para tabela LaTeX
latex_table <- xtable(df, caption = title)

# Definir a formatação dos números como inteiros
digits(latex_table) <- c(0, 0, 0, 0, 0, 0, 0)

# Salvar a tabela em um arquivo .tex
print.xtable(latex_table, file = "numero_de_snps_selecionados_simulacao6_005.tex", include.rownames = FALSE, sanitize.text.function = function(x){x})

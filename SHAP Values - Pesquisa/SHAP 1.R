# 1. Simulação dos Dados

# Instalação e carregamento dos pacotes necessários
# install.packages(c("randomForest", "fastshap", "ggplot2", "data.table", "reshape2", "ggforce", "viridis"))
library(randomForest)
library(fastshap)
library(ggplot2)
library(data.table)
library(reshape2)
library(ggforce)
library(viridis)

set.seed(42)  # Para reprodutibilidade

n_individuals <- 200
n_snps <- 100

# Gerando genótipos aleatórios (0, 1, 2)
X <- matrix(sample(0:2, n_individuals * n_snps, replace = TRUE), nrow = n_individuals, ncol = n_snps)
snp_columns <- paste0("SNP_", 1:n_snps)
df_snps <- as.data.frame(X)
colnames(df_snps) <- snp_columns

# Selecionando alguns SNPs para ter efeito no fenótipo
effect_snps <- c("SNP_5", "SNP_20", "SNP_50")
effects <- c(1.5, -2.0, 3.0)

# Calculando o fenótipo
noise <- rnorm(n_individuals, mean = 0, sd = 1)
phenotype <- as.matrix(df_snps[, effect_snps]) %*% effects + noise

df_snps$Phenotype <- phenotype

# 2. Construção do Modelo

X <- df_snps[, snp_columns]
y <- df_snps$Phenotype

# Divisão em conjunto de treinamento e teste
set.seed(42)
train_indices <- sample(1:n_individuals, size = 0.8 * n_individuals)
X_train <- X[train_indices, ]
X_test <- X[-train_indices, ]
y_train <- y[train_indices]
y_test <- y[-train_indices]

# Treinando o modelo
model <- randomForest(x = X_train, y = y_train, ntree = 100)

# 3. Cálculo dos Valores SHAP

# Função preditora para o pacote 'fastshap'
predict_function <- function(object, newdata) {
  predict(object, newdata = newdata)
}

# Calculando os valores SHAP
set.seed(42)
shap_values <- fastshap::explain(
  object = model,
  X = X_test,
  pred_wrapper = predict_function,
  nsim = 50
)

# 4. Análise dos Resultados

# Convertendo os valores SHAP para data frame
shap_df <- as.data.frame(shap_values)
shap_df$ID <- 1:nrow(shap_df)

# Cálculo da importância média absoluta
mean_abs_shap <- colMeans(abs(shap_values))

# Ordenando os SNPs por importância
importance_df <- data.frame(
  SNP = names(mean_abs_shap),
  Importance = mean_abs_shap
)
importance_df <- importance_df[order(-importance_df$Importance), ]

# Visualização dos resultados

# Gráfico de importância média absoluta
ggplot(importance_df, aes(x = reorder(SNP, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Importância dos SNPs (Valores SHAP)",
       x = "SNP",
       y = "Importância Média Absoluta") +
  theme_minimal()

# Gráfico de resumo dos valores SHAP para os principais SNPs
# Selecionando os SNPs mais importantes para visualização
top_snps <- importance_df$SNP[1:10]

# Criando um data frame para plotagem
plot_data <- melt(shap_df[, c(top_snps, "ID")], id.vars = "ID", variable.name = "SNP", value.name = "SHAP_value")
plot_data$Feature_value <- as.vector(as.matrix(X_test[, top_snps]))

ggplot(plot_data, aes(x = SHAP_value, y = SNP, color = Feature_value)) +
  geom_sina(alpha = 0.7) +
  scale_color_viridis(option = "D") +
  labs(title = "Valores SHAP para os Principais SNPs",
       x = "Valor SHAP",
       y = "SNP",
       color = "Valor do SNP") +
  theme_minimal()

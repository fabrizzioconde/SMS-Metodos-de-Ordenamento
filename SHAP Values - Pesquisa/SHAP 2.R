# Instalando e carregando pacotes
packages <- c("scrime","e1071","kernlab","randomForest", 
              "doParallel", "GA", "ggplot2","care", "fastshap", "data.table", "reshape2", "ggforce", "viridis")

# Instalando os pacotes que ainda nao foram instalados
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Carregando os pacotes
invisible(lapply(packages, library, character.only = TRUE))

# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Parâmetros para simulação

dados<-list();                   # Cria lista de dados.
genotipo<-list();                # Cria lista de genótipos.
num_individuos<-1000;            # Número de indivíduos na amostra
num_snp<-100;                    # Total de marcadores simulados
list.snp<-list(c(1,2,3),c(4,5,6),c(7,8,9))  # Indica quais são os marcadores causais
list.ia<-list(c(1,2,3),c(-1,2,3),c(1,-2,3))  # Construção das interações entre os marcadores causais

# Executa a simulação do genótipo e do fenótipo simultaneamente
simulacao<-simulateSNPglm(n.obs=num_individuos,
                          n.snp=num_snp,
                          list.ia=list.ia,
                          list.snp=list.snp,
                          beta0=640,
                          beta=c(3,3,3),
                          maf=c(0.1,0.4),
                          err.fun=rnorm,
                          rand=123)

# Definindo o genótipo
genotipo[[1]]<-as.data.frame(simulacao$x)

# Definindo o fenótipo
fenotipo<-as.data.frame(simulacao$y)
names(fenotipo)<-"fenotipo"
pdf(file="Histograma_fenotipo_3.pdf",height=5,width=9)
hist(fenotipo$fenotipo,xlab="Fenótipo simulado",ylab="Número de touros",col="gray",main="")
dev.off()
pdf(file="Boxplot_fenotipo_3.pdf",height=5,width=9)
boxplot(fenotipo$fenotipo, ylab="Fenótipo simulado")
dev.off()

# Testando a normalidade do Fenótipo
teste_normal<-shapiro.test(fenotipo$fenotipo) # Não é normal

# Definindo o dataframe genotipo com fenotipo
dados[[1]]<-as.data.frame(cbind(genotipo[[1]],fenotipo))
colnames(dados[[1]])[ncol(dados[[1]])]<-"fenotipo"

save(dados,file="Simulacao3.RData")

# Função para realizar para o SVR com k-fold.
validacao_cruzada<-function(data,folds,gamma,cost,epsilon,kernel) 
{
  datalength = nrow(data)
  index <- 1:datalength
  size = trunc(datalength/folds)
  set.seed(123)
  geral = matrix(sample(index), ncol=folds, byrow=TRUE)
  
  mse = vector('double', folds)
  mape = vector('double', folds)
  corr = vector('double', folds)
  r2 = vector('double', folds)
  r2_adj = vector('double', folds)
  
  for(i in 1:folds)
  {
    testv  = geral[,i];
    trainv = index[-testv]
    testset<-data.frame()
    trainset<-data.frame()
    
    testset  = as.data.frame(na.omit(data[testv,]))
    trainset = as.data.frame(na.omit(data[trainv,]))
    
    svmR.model <- svm(fenotipo~., data = trainset, kernel=kernel, gamma=gamma, cost=cost,epsilon=epsilon)
    seltestset = as.data.frame(testset[,-ncol(testset)]);
    names(seltestset) = names(testset)[-ncol(testset)];
    svmR.pred <- predict(svmR.model, seltestset) 
    
    mse[i] = sum((svmR.pred-testset$fenotipo)^2)/dim(testset)[1]
    mape[i] = (sum( abs( (testset$fenotipo - svmR.pred) / testset$fenotipo) ) / length(testset$fenotipo)) * 100;
    corr[i]<-cor(svmR.pred,testset$fenotipo,method="pearson")
    r2[i] <- cor(svmR.pred,testset$fenotipo,method="pearson")^2
    r2_adj[i] <- 1-(1-r2[i])*((nrow(testset)-1)/(nrow(testset)-ncol(testset)-1))
  }
  
  resultado<-c(mean(mse),sd(mse),mean(mape),sd(mape),mean(corr),sd(corr),mean(r2),sd(r2),mean(r2_adj),sd(r2_adj))
  
  return(resultado)
}

########

# Função para o cálculo do valor p bruto e ajustado pela correção de Bonferroni
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

################ Início do SMS##################

# Definição de listas para cada kernel do SVR
mean_svr_SHAP_list<-list() # Cria lista para a média do SVR sobre o ranking dos valores SHAP
GA<-list()
minimo<-list()
corte<-list()
snps_selec_corte<-list()
snps_selec_ref<-list()
percentual_snps<-0.95;

i=1 # Contador do kernel do SVR utilizado

# Random Forest
ntree<-4000           # Número de árvores dentro da floresta aleatória.

data_temp<-as.data.frame(dados[[1]]) # Transforma a base de dados em dataframe.

set.seed(1) # Semente aleatória para floresta aleatória.

RF<- randomForest(fenotipo~., data=data_temp,
                  ntree=ntree,
                  mtry=ncol(data_temp)-1,
                  importance=TRUE)    

# Cálculo dos valores SHAP
# Função preditora para o 'fastshap'
predict_function <- function(object, newdata) {
  predict(object, newdata = newdata)
}

# Calculando os valores SHAP
set.seed(42)  # Para reprodutibilidade
shap_values <- fastshap::explain(
  object = RF,
  X = data_temp[, -ncol(data_temp)],  # Excluindo a coluna 'fenotipo'
  pred_wrapper = predict_function,
  nsim = 50
)

# Calculando a importância média absoluta dos valores SHAP
mean_abs_shap <- colMeans(abs(shap_values))

# Ordenando os SNPs por importância
rank_SHAP <- sort(mean_abs_shap, decreasing = TRUE)
View(rank_SHAP)

# Gráfico de Importância dos Valores SHAP
library(ggplot2)

shap_importance_df <- data.frame(
  SNP = names(rank_SHAP),
  Importance = rank_SHAP
)

ggplot(shap_importance_df, aes(x = reorder(SNP, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Importância dos SNPs (Valores SHAP)",
       x = "SNP",
       y = "Importância Média Absoluta") +
  theme_minimal()



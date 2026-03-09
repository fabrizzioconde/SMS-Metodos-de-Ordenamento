# Base de dados Completa com o pacote CAR
dados$x = dados[[1]][-length(names(dados[[1]]))]
dados$y = dados[[1]]$fenotipo
xnames = names(dados[[1]][-length(names(dados[[1]]))])

# Rank CARSCORE (já foi feito)
car = carscore(dados$x, dados$y, diagonal = FALSE) # CARSCORE.
rank_car <- xnames[order(car^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_car
SMS_Completo_Corr_1_CAR

# Rank Correlação Marginal
mcorr = carscore(dados$x, dados$y, lambda = 0, diagonal = TRUE) # Correlação Marginal.
rank_mcorr <- xnames[order(car^2, decreasing=TRUE)] # Lista ordenada decrescente com os nomes dos SNPs
rank_mcorr # Rank dos SNPs ordenados de forma decrescente pelo Correlação Marginal.
# SMS_Completo_Corr_1_CORRMAR é o nome do arquivo de script do R

# Rank Correlação Parcial
pcor = pcor.shrink(cbind(dados$x,dados$y), lambda=0, verbose=FALSE)[-1,1]
rank_pcorr <- xnames[order(pcor^2, decreasing=TRUE)]
rank_pcorr
# SMS_Completo_Corr_1_CORRPAR é o nome do arquivo de script do R
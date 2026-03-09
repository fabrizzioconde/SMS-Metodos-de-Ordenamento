####ADAPTAÇÃO DO SMS PARA QUE NÃO OCORRA CORTES ACIMA DO NÚMERO DE SNPs do dataframe

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]] <- which.min(mean_svr_CAR_list[[i]])
corte[[i]] <- (minimo[[i]] + 1) * passo

# VERIFICAÇÃO CONDICIONAL PARA O CORTE
limite <- dim(dados[[1]][-length(dados[[1]])])[2]
corte[[i]] <- ifelse(corte[[i]] > limite, limite, corte[[i]])
